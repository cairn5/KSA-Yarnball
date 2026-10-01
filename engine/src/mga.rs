//! Flyby-sequence search: which planets to fly past, and when, from a departure planet to a target.
//!
//! Model (ported from the pykep prototype `mga_dp_search.py`): each leg is a Lambert arc between two
//! planet encounters, fixed by its two dates. Costs are the burn out of a 200 km parking orbit, a
//! powered-flyby burn wherever gravity alone can't turn the incoming v∞ into the outgoing one
//! (pykep's `fb_vel` model), an optional capture burn, and an optional cost per day of flight.
//!
//! Every cost depends only on neighbouring encounters, so for a fixed sequence dynamic programming
//! on a date grid finds the cheapest route exactly: per leg it keeps only the cheapest way found to
//! reach it. Sequences form a tree (each prefix's DP table extends to all its children), searched
//! depth first with leg tables shared between sequences. A branch is pruned once everything in it
//! already costs more than the best complete route plus a margin.

use std::collections::HashMap;

use wasm_bindgen::prelude::*;

use crate::bodies::{planets, Planet};
use crate::lambert::{conic_xy, lambert, norm, sub, transfer_angle, V3};
use crate::{AU_KM, DAY_S, MU_SUN_KM};

/// Parking orbit altitude for the departure burn, km.
pub const PARKING_ALT_KM: f64 = 200.0;
/// Eccentricity of the capture orbit (periapsis at the planet's lowest flyby altitude).
pub const CAPTURE_E: f64 = 0.9;
/// Heliocentric speed assumed for the quickest legs between different planets, km/s.
pub const FAST_TRANSFER_KMS: f64 = 30.0;
/// Routes launching closer together than this count as the same launch window, days.
pub const DISTINCT_LAUNCH_DAYS: f64 = 180.0;
/// For a fixed sequence, how many launch windows to report.
pub const FIXED_SEQUENCE_WINDOWS: usize = 5;
/// For a fixed sequence, report launch windows within this fraction of the best one. Wider than
/// the margin between sequences; it also drops "windows" that are just the edge of the date range.
pub const FIXED_SEQUENCE_MARGIN: f64 = 0.5;

/// Burn from a circular parking orbit to leave `p` with excess speed `vinf` (km/s).
pub fn departure_dv(p: &Planet, vinf: f64) -> f64 {
    let r = p.radius_km + PARKING_ALT_KM;
    (vinf * vinf + 2.0 * p.mu / r).sqrt() - (p.mu / r).sqrt()
}

/// Burn at periapsis to capture into an orbit of eccentricity `CAPTURE_E` around `p`, km/s.
pub fn capture_dv(p: &Planet, vinf: f64) -> f64 {
    let rp = p.safe_radius * p.radius_km;
    ((vinf * vinf + 2.0 * p.mu / rp).sqrt() - (p.mu * (1.0 + CAPTURE_E) / rp).sqrt()).abs()
}

/// Largest turn of v∞ (radians) gravity alone gives at the lowest allowed periapsis.
pub fn max_turn(p: &Planet, vinf: f64) -> f64 {
    let rp = p.safe_radius * p.radius_km;
    2.0 * (1.0 / (1.0 + rp * vinf * vinf / p.mu)).asin()
}

/// Powered-flyby burn at `p` turning incoming v∞ `vin` into outgoing `vout`, km/s (as pykep.fb_vel).
pub fn flyby_dv(vin: V3, vout: V3, p: &Planet) -> f64 {
    let (n1, n2) = (norm(vin), norm(vout));
    flyby_dv_with(vin, n1, max_turn(p, n1), vout, n2)
}

/// `flyby_dv` with the norms and maximum turn precomputed.
fn flyby_dv_with(vin: V3, n1: f64, alpha_max: f64, vout: V3, n2: f64) -> f64 {
    let c = ((vin[0] * vout[0] + vin[1] * vout[1] + vin[2] * vout[2]) / (n1 * n2)).clamp(-1.0, 1.0);
    let excess = (c.acos() - alpha_max).max(0.0);
    (n1 * n1 + n2 * n2 - 2.0 * n1 * n2 * excess.cos()).max(0.0).sqrt()
}

/// What to search for.
#[derive(Clone, Debug)]
pub struct Problem {
    /// Planet indices into `planets()`.
    pub from: usize,
    pub to: usize,
    /// Launch window, days after J2000.
    pub launch_start: f64,
    pub launch_end: f64,
    pub max_flybys: usize,
    /// Capture at the target (true) or fly past it (false).
    pub capture: bool,
    /// Cost of flight time, km/s per day.
    pub time_cost: f64,
    /// Longest leg, and the latest arrival after the launch window closes, days.
    pub max_tof: f64,
    /// Date grid step, days.
    pub step: f64,
    /// Keep routes within this fraction of the best score.
    pub margin: f64,
}

/// Cost breakdown of a route, km/s.
#[derive(Clone, Debug)]
pub struct Cost {
    pub vinf_launch: f64,
    pub departure: f64,
    pub flybys: Vec<f64>,
    pub arrival: f64,
    /// Flight-time penalty, in the same units.
    pub time: f64,
}

impl Cost {
    /// Spacecraft and launch Δv, without the time penalty.
    pub fn dv(&self) -> f64 {
        self.departure + self.flybys.iter().sum::<f64>() + self.arrival
    }

    /// What the search minimises.
    pub fn score(&self) -> f64 {
        self.dv() + self.time
    }
}

#[derive(Clone, Debug)]
pub struct Route {
    pub bodies: Vec<usize>,
    /// Encounter dates, days after J2000.
    pub dates: Vec<f64>,
    pub cost: Cost,
    /// Score as the grid search computed it (single precision), km/s.
    pub grid_score: f64,
}

impl Problem {
    /// Allowed duration of a leg from planet `a` to `b`, days, from their orbital periods. `None` if
    /// no duration fits.
    pub fn leg_bounds(&self, ps: &[Planet], a: usize, b: usize) -> Option<(f64, f64)> {
        let (pa, pb) = (&ps[a], &ps[b]);
        let (lo, hi) = if a == b {
            let p = pa.period_days();
            (0.6 * p, 2.2 * p)
        } else {
            // Shortest: a quarter of the Hohmann time, or crossing the gap between the orbits at
            // FAST_TRANSFER_KMS, whichever is shorter. Flyby-boosted outer legs are far quicker
            // than Hohmann: Voyager 2 took 3.6 years from Uranus to Neptune, against 61.
            let hohmann = 182.625 * (0.5 * (pa.a_au + pb.a_au)).powf(1.5);
            let crossing = (pa.a_au - pb.a_au).abs() * AU_KM / FAST_TRANSFER_KMS / DAY_S;
            ((0.25 * hohmann).min(crossing), 1.5 * hohmann)
        };
        let (lo, hi) = (lo.max(self.step), hi.min(self.max_tof));
        (hi > lo + self.step).then_some((lo, hi))
    }

    /// Cost of flying `bodies` with encounters on `dates`. `None` if a leg has no Lambert solution.
    pub fn route_cost(&self, ps: &[Planet], bodies: &[usize], dates: &[f64]) -> Option<Cost> {
        let n = bodies.len() - 1;
        let states: Vec<_> = bodies.iter().zip(dates).map(|(&b, &t)| ps[b].state_at(t)).collect();
        let mut dep = Vec::with_capacity(n);
        let mut arr = Vec::with_capacity(n);
        for k in 0..n {
            let tof = dates[k + 1] - dates[k];
            if tof <= 0.0 {
                return None;
            }
            let (v1, v2) = lambert(states[k].0, states[k + 1].0, tof * DAY_S, MU_SUN_KM)?;
            dep.push(sub(v1, states[k].1));
            arr.push(sub(v2, states[k + 1].1));
        }
        let vinf_launch = norm(dep[0]);
        let vinf_arrive = norm(arr[n - 1]);
        Some(Cost {
            vinf_launch,
            departure: departure_dv(&ps[bodies[0]], vinf_launch),
            flybys: (1..n).map(|k| flyby_dv(arr[k - 1], dep[k], &ps[bodies[k]])).collect(),
            arrival: if self.capture { capture_dv(&ps[bodies[n]], vinf_arrive) } else { 0.0 },
            time: self.time_cost * (dates[n] - dates[0]),
        })
    }

    /// Planets worth flying past: none further out than the departure or target, not the target
    /// itself, and Mercury only for Mercury missions.
    fn candidates(&self, ps: &[Planet]) -> Vec<usize> {
        let reach = ps[self.from].a_au.max(ps[self.to].a_au) * 1.05;
        let mercury = self.from == 0 || self.to == 0;
        (0..ps.len()).filter(|&k| k != self.to && ps[k].a_au <= reach && (k != 0 || mercury)).collect()
    }
}

/// Search statistics, reported as it runs.
#[derive(Clone, Debug, Default)]
pub struct Progress {
    /// Sequence prefix being searched.
    pub current: Vec<usize>,
    pub nodes: usize,
    pub pruned: usize,
    pub sequences: usize,
    pub lambert: usize,
    /// Best score found so far, km/s.
    pub best: f64,
}

/// One Lambert arc on the grid: v∞ at each end, NaN if there is no solution. Single precision is
/// plenty for the grid; `Problem::route_cost` recomputes reported routes in double.
#[derive(Clone, Copy)]
struct Arc {
    vd: [f32; 3],
    va: [f32; 3],
}

const NO_ARC: Arc = Arc { vd: [f32::NAN; 3], va: [f32::NAN; 3] };

/// All arcs from one planet to another: by departure grid index, then by duration index.
/// Rows are computed on first use.
struct Table {
    /// Shortest leg, grid steps.
    lo: usize,
    /// Number of durations.
    m: usize,
    rows: Vec<Option<Vec<Arc>>>,
}

/// DP state after a leg from `a` to `b`: `v[g * m + i]` is the cheapest score so far for a leg
/// departing on grid date `g` with duration index `i`; `back` holds the previous leg's duration index.
struct Layer {
    a: usize,
    b: usize,
    lo: usize,
    m: usize,
    v: Vec<f32>,
    back: Vec<u16>,
}

impl Layer {
    fn min(&self) -> f32 {
        self.v.iter().copied().fold(f32::INFINITY, f32::min)
    }
}

struct Search<'a> {
    prob: &'a Problem,
    ps: Vec<Planet>,
    n_grid: usize,
    n_launch: usize,
    /// Planet states on every grid date: `eph[planet][g]`, km and km/s.
    eph: Vec<Vec<(V3, V3)>>,
    tables: HashMap<(usize, usize), Option<Table>>,
    candidates: Vec<usize>,
    /// Restricts the search to one sequence (after the departure planet).
    only: Option<Vec<usize>>,
    /// Routes to keep per sequence, each from a different launch window.
    windows: usize,
    stack: Vec<Layer>,
    best: f64,
    routes: Vec<Route>,
    progress: Progress,
    report: &'a mut dyn FnMut(&Progress),
}

fn to_f32(v: V3) -> [f32; 3] {
    [v[0] as f32, v[1] as f32, v[2] as f32]
}
fn to_f64(v: [f32; 3]) -> V3 {
    [v[0] as f64, v[1] as f64, v[2] as f64]
}

impl<'a> Search<'a> {
    fn new(prob: &'a Problem, only: Option<Vec<usize>>, report: &'a mut dyn FnMut(&Progress)) -> Self {
        let ps = planets();
        let n_launch = ((prob.launch_end - prob.launch_start) / prob.step).floor() as usize + 1;
        let n_grid = ((prob.launch_end + prob.max_tof - prob.launch_start) / prob.step).floor() as usize + 1;
        let eph = ps
            .iter()
            .map(|p| (0..n_grid).map(|g| p.state_at(prob.launch_start + g as f64 * prob.step)).collect())
            .collect();
        let candidates = prob.candidates(&ps);
        Search {
            prob,
            ps,
            n_grid,
            n_launch,
            eph,
            tables: HashMap::new(),
            candidates,
            windows: if only.is_some() { FIXED_SEQUENCE_WINDOWS } else { 1 },
            only,
            stack: Vec::new(),
            best: f64::INFINITY,
            routes: Vec::new(),
            progress: Progress { best: f64::INFINITY, ..Default::default() },
            report,
        }
    }

    /// Scores at or above this can't lead to a route worth keeping.
    fn bound(&self) -> f32 {
        (self.best * (1.0 + self.prob.margin)) as f32
    }

    /// Makes sure the table for `key` exists (if the leg is possible at all) and has row `g`.
    fn ensure_row(&mut self, key: (usize, usize), g: usize) -> bool {
        let (prob, ps, n_grid) = (self.prob, &self.ps, self.n_grid);
        let table = self.tables.entry(key).or_insert_with(|| {
            prob.leg_bounds(ps, key.0, key.1).map(|(lo, hi)| {
                let lo_i = (lo / prob.step).ceil() as usize;
                let hi_i = (hi / prob.step).floor() as usize;
                Table { lo: lo_i, m: hi_i + 1 - lo_i, rows: vec![None; n_grid] }
            })
        });
        let Some(table) = table else { return false };
        if table.rows[g].is_none() {
            let (r1, vp1) = self.eph[key.0][g];
            let row = (0..table.m)
                .map(|i| {
                    let g2 = g + table.lo + i;
                    if g2 >= self.n_grid {
                        return NO_ARC;
                    }
                    let (r2, vp2) = self.eph[key.1][g2];
                    self.progress.lambert += 1;
                    match lambert(r1, r2, (table.lo + i) as f64 * prob.step * DAY_S, MU_SUN_KM) {
                        Some((v1, v2)) => Arc { vd: to_f32(sub(v1, vp1)), va: to_f32(sub(v2, vp2)) },
                        None => NO_ARC,
                    }
                })
                .collect();
            table.rows[g] = Some(row);
        }
        true
    }

    fn table(&self, key: (usize, usize)) -> &Table {
        self.tables[&key].as_ref().unwrap()
    }

    /// The first leg, from the departure planet to `b`, launching in the window.
    fn first_layer(&mut self, b: usize) -> Option<Layer> {
        let key = (self.prob.from, b);
        let bound = self.bound();
        let mut layer: Option<Layer> = None;
        for g in 0..self.n_launch {
            if !self.ensure_row(key, g) {
                return None;
            }
            let t = self.table(key);
            let l = layer.get_or_insert_with(|| Layer {
                a: key.0,
                b,
                lo: t.lo,
                m: t.m,
                v: vec![f32::INFINITY; self.n_grid * t.m],
                back: vec![0; self.n_grid * t.m],
            });
            let row = t.rows[g].as_ref().unwrap();
            for (i, arc) in row.iter().enumerate() {
                if arc.vd[0].is_nan() {
                    continue;
                }
                let days = (t.lo + i) as f64 * self.prob.step;
                let v = departure_dv(&self.ps[key.0], norm(to_f64(arc.vd))) + self.prob.time_cost * days;
                if (v as f32) < bound {
                    l.v[g * t.m + i] = v as f32;
                }
            }
        }
        layer
    }

    /// Extends `prev` (ending at a flyby of `prev.b`) with a leg to `c`.
    fn extend(&mut self, prev: &Layer, c: usize) -> Option<Layer> {
        let key_in = (prev.a, prev.b);
        let key_out = (prev.b, c);
        let flyby = self.ps[prev.b].clone();
        let bound = self.bound();
        let mut next: Option<Layer> = None;
        // Incoming arcs at one encounter date: (score, v∞, |v∞|, max turn, duration index).
        let mut incoming: Vec<(f64, V3, f64, f64, usize)> = Vec::new();

        for e in 0..self.n_grid {
            incoming.clear();
            let t_in = self.table(key_in);
            for i in 0..prev.m {
                if e < prev.lo + i {
                    break;
                }
                let g = e - prev.lo - i;
                let v = prev.v[g * prev.m + i];
                if v < bound {
                    let vin = to_f64(t_in.rows[g].as_ref().unwrap()[i].va);
                    let n1 = norm(vin);
                    incoming.push((v as f64, vin, n1, max_turn(&flyby, n1), i));
                }
            }
            if incoming.is_empty() {
                continue;
            }
            if !self.ensure_row(key_out, e) {
                return None;
            }
            let t_out = self.table(key_out);
            let l = next.get_or_insert_with(|| Layer {
                a: prev.b,
                b: c,
                lo: t_out.lo,
                m: t_out.m,
                v: vec![f32::INFINITY; self.n_grid * t_out.m],
                back: vec![0; self.n_grid * t_out.m],
            });
            let row = t_out.rows[e].as_ref().unwrap();
            for (j, arc) in row.iter().enumerate() {
                if arc.vd[0].is_nan() {
                    continue;
                }
                let vout = to_f64(arc.vd);
                let n2 = norm(vout);
                let (mut best, mut arg) = (f64::INFINITY, 0);
                for &(v, vin, n1, amax, i) in &incoming {
                    // The burn is at least the change in speed.
                    if v + (n1 - n2).abs() >= best {
                        continue;
                    }
                    let total = v + flyby_dv_with(vin, n1, amax, vout, n2);
                    if total < best {
                        best = total;
                        arg = i;
                    }
                }
                let total = best + self.prob.time_cost * (t_out.lo + j) as f64 * self.prob.step;
                if (total as f32) < bound {
                    l.v[e * t_out.m + j] = total as f32;
                    l.back[e * t_out.m + j] = arg as u16;
                }
            }
        }
        next
    }

    /// Records the best routes ending with the top layer, which arrives at the target: the best
    /// one, or with `windows` > 1 the best in each distinct launch window.
    fn arrive(&mut self) {
        let top = self.stack.last().unwrap();
        let t = self.table((top.a, top.b));
        let target = &self.ps[top.b];
        let bound = self.bound() as f64;
        let mut finals: Vec<(f64, usize, usize)> = Vec::new();
        for g in 0..self.n_grid {
            for i in 0..top.m {
                let v = top.v[g * top.m + i];
                if !v.is_finite() {
                    continue;
                }
                let arc = t.rows[g].as_ref().unwrap()[i];
                let arrival = if self.prob.capture { capture_dv(target, norm(to_f64(arc.va))) } else { 0.0 };
                let total = v as f64 + arrival;
                if total < bound {
                    finals.push((total, g, i));
                }
            }
        }
        self.progress.sequences += 1;
        finals.sort_by(|a, b| a.0.total_cmp(&b.0));

        let bodies: Vec<usize> = std::iter::once(self.prob.from).chain(self.stack.iter().map(|l| l.b)).collect();
        let mut launches: Vec<f64> = Vec::new();
        for (total, g, i) in finals {
            if launches.len() == self.windows {
                break;
            }
            let dates = self.backtrack(g, i);
            if launches.iter().any(|&l| (l - dates[0]).abs() < DISTINCT_LAUNCH_DAYS) {
                continue;
            }
            launches.push(dates[0]);
            if let Some(cost) = self.prob.route_cost(&self.ps, &bodies, &dates) {
                self.routes.push(Route { bodies: bodies.clone(), dates, cost, grid_score: total });
                self.best = self.best.min(total);
                self.progress.best = self.best;
            }
        }
    }

    /// Encounter dates of the route whose last leg departs on grid date `g` with duration index `i`,
    /// found by walking back through the layers.
    fn backtrack(&self, mut g: usize, mut i: usize) -> Vec<f64> {
        let n = self.stack.len();
        let mut idx = vec![0; n + 1];
        idx[n] = g + self.stack[n - 1].lo + i;
        for k in (0..n).rev() {
            idx[k] = g;
            if k > 0 {
                let mp = self.stack[k].back[g * self.stack[k].m + i] as usize;
                g -= self.stack[k - 1].lo + mp;
                i = mp;
            }
        }
        idx.iter().map(|&j| self.prob.launch_start + j as f64 * self.prob.step).collect()
    }

    /// Next planets to try after the top layer: the target first, so a bound comes early.
    fn children(&self) -> Vec<usize> {
        let depth = self.stack.len();
        if let Some(only) = &self.only {
            return only.get(depth).copied().into_iter().collect();
        }
        if depth > self.prob.max_flybys {
            return vec![];
        }
        let mut c = vec![self.prob.to];
        if depth < self.prob.max_flybys {
            c.extend(&self.candidates);
        }
        c
    }

    fn push_and_descend(&mut self, layer: Option<Layer>) {
        self.progress.nodes += 1;
        match layer {
            Some(l) if l.min() < self.bound() => {
                self.stack.push(l);
                self.progress.current = std::iter::once(self.prob.from).chain(self.stack.iter().map(|l| l.b)).collect();
                (self.report)(&self.progress);
                self.descend();
                self.stack.pop();
            }
            _ => self.progress.pruned += 1,
        }
    }

    fn descend(&mut self) {
        if self.stack.last().map(|l| l.b) == Some(self.prob.to) {
            self.arrive();
            return;
        }
        for c in self.children() {
            let layer = if self.stack.is_empty() {
                self.first_layer(c)
            } else {
                let top = self.stack.pop().unwrap();
                let next = self.extend(&top, c);
                self.stack.push(top);
                next
            };
            self.push_and_descend(layer);
        }
    }
}

/// Best route for each flyby sequence within `prob.margin` of the overall best, cheapest first.
pub fn search(prob: &Problem, report: &mut dyn FnMut(&Progress)) -> Vec<Route> {
    run(prob, None, report)
}

/// Best route in each of up to `FIXED_SEQUENCE_WINDOWS` distinct launch windows for one fixed
/// sequence (planet indices after the departure planet, ending at the target), cheapest first.
/// `prob.margin` and `prob.max_flybys` are ignored.
pub fn search_sequence(prob: &Problem, sequence: &[usize]) -> Vec<Route> {
    let prob = Problem { margin: FIXED_SEQUENCE_MARGIN, ..prob.clone() };
    run(&prob, Some(sequence.to_vec()), &mut |_| {})
}

fn run(prob: &Problem, only: Option<Vec<usize>>, report: &mut dyn FnMut(&Progress)) -> Vec<Route> {
    let mut s = Search::new(prob, only, report);
    s.descend();
    (s.report)(&s.progress);
    let cutoff = s.best * (1.0 + prob.margin);
    let mut routes: Vec<Route> = s.routes.into_iter().filter(|r| r.grid_score <= cutoff).collect();
    routes.sort_by(|a, b| a.grid_score.total_cmp(&b.grid_score));
    routes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unpowered_flyby_within_turn_limit_is_free() {
        let earth = &planets()[2];
        let vin = [5.0, 0.0, 0.0];
        let a = 0.5 * max_turn(earth, 5.0);
        let vout = [5.0 * a.cos(), 5.0 * a.sin(), 0.0];
        assert!(flyby_dv(vin, vout, earth) < 1e-9);
        // Same direction, different speed: the burn is the speed change.
        assert!((flyby_dv(vin, [7.0, 0.0, 0.0], earth) - 2.0).abs() < 1e-9);
    }

    fn mjd(y: i32, m: u32, d: u32) -> f64 {
        // Days from 2000-01-01 12:00 to y-m-d 00:00 (proleptic Gregorian).
        let days_from_civil = |y: i32, m: u32, d: u32| {
            let y = if m <= 2 { y - 1 } else { y };
            let era = y.div_euclid(400);
            let yoe = (y - era * 400) as i64;
            let mp = (m as i64 + 9) % 12;
            let doy = (153 * mp + 2) / 5 + d as i64 - 1;
            let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
            era as i64 * 146_097 + doe - 719_468
        };
        (days_from_civil(y, m, d) - days_from_civil(2000, 1, 1)) as f64 - 0.5
    }

    #[test]
    fn finds_cassinis_window() {
        // EVVEJS from an 18-month window around Cassini's launch (1997-10-15), as in the prototype.
        let prob = Problem {
            from: 2,
            to: 5,
            launch_start: mjd(1997, 1, 1),
            launch_end: mjd(1998, 7, 1),
            max_flybys: 4,
            capture: true,
            time_cost: 0.0,
            max_tof: 3650.0,
            step: 5.0,
            margin: 0.0,
        };
        let r = &search_sequence(&prob, &[1, 1, 2, 4, 5])[0];
        assert!(r.dates[0] > mjd(1997, 9, 1) && r.dates[0] < mjd(1997, 12, 31), "launch {}", r.dates[0]);
        assert!(r.dates[1] > mjd(1998, 3, 1) && r.dates[1] < mjd(1998, 7, 1), "first Venus {}", r.dates[1]);
        assert!(r.dates[3] > mjd(1999, 7, 1) && r.dates[3] < mjd(1999, 10, 1), "Earth {}", r.dates[3]);
        // Grid and continuous scoring agree.
        assert!((r.grid_score - r.cost.score()).abs() < 1e-3, "{} vs {}", r.grid_score, r.cost.score());
    }

    #[test]
    fn fixed_sequence_reports_each_launch_window() {
        // Direct Earth–Mars over five years: windows about 26 months apart.
        let prob = Problem {
            from: 2,
            to: 3,
            launch_start: mjd(2026, 1, 1),
            launch_end: mjd(2031, 1, 1),
            max_flybys: 0,
            capture: true,
            time_cost: 0.0,
            max_tof: 500.0,
            step: 5.0,
            margin: 0.0,
        };
        let routes = search_sequence(&prob, &[3]);
        assert!(routes.len() >= 2, "{} windows", routes.len());
        for (k, a) in routes.iter().enumerate() {
            for b in &routes[k + 1..] {
                assert!((a.dates[0] - b.dates[0]).abs() >= DISTINCT_LAUNCH_DAYS);
            }
        }
        assert!(routes.windows(2).all(|w| w[0].grid_score <= w[1].grid_score));
    }

    #[test]
    fn search_beats_or_matches_direct() {
        // Earth to Jupiter: flybys may help, but never make the best route worse than the direct one.
        let prob = Problem {
            from: 2,
            to: 4,
            launch_start: mjd(2030, 1, 1),
            launch_end: mjd(2032, 1, 1),
            max_flybys: 2,
            capture: true,
            time_cost: 0.0,
            max_tof: 2500.0,
            step: 10.0,
            margin: 0.1,
        };
        let direct = &search_sequence(&prob, &[4])[0];
        let routes = search(&prob, &mut |_| {});
        assert!(!routes.is_empty());
        assert!(routes[0].grid_score <= direct.grid_score + 1e-6);
        assert!(routes.iter().all(|r| r.grid_score <= routes[0].grid_score * 1.1 + 1e-6));
    }
}

/// A route found by `find_routes`, for the interface. Costs in m/s.
#[wasm_bindgen(getter_with_clone)]
pub struct FoundRoute {
    /// Planet indices into `planets()`, departure first.
    pub bodies: Vec<u32>,
    /// Encounter dates, days after J2000.
    pub dates: Vec<f64>,
    /// Launch energy, km²/s².
    pub c3: f64,
    pub departure: f64,
    pub flybys: Vec<f64>,
    pub arrival: f64,
    /// Departure + flybys + arrival.
    pub dv: f64,
    /// Δv plus the flight-time penalty: what the search ranks by.
    pub score: f64,
}

/// Searches flyby sequences from planet `from` to `to` (indices into `planets()`). Dates are days
/// after J2000; `time_cost` is m/s per year of flight. `flybys` fixes the planets flown past (empty
/// to search for the best sequence); with it fixed, the best route in each launch window comes
/// back instead. `progress` is called with a status string.
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn find_routes(
    from: usize,
    to: usize,
    launch_start: f64,
    launch_end: f64,
    max_flybys: usize,
    capture: bool,
    time_cost: f64,
    max_tof: f64,
    step: f64,
    margin: f64,
    flybys: Vec<u32>,
    progress: &js_sys::Function,
) -> Vec<FoundRoute> {
    let fixed = !flybys.is_empty();
    let prob = Problem {
        from,
        to,
        launch_start,
        launch_end,
        max_flybys,
        capture,
        time_cost: time_cost / 1000.0 / 365.25,
        max_tof,
        step,
        margin: if fixed { FIXED_SEQUENCE_MARGIN } else { margin },
    };
    let names: Vec<String> = planets().iter().map(|p| p.name.chars().next().unwrap().to_string()).collect();
    let mut report = |p: &Progress| {
        let seq: String = p.current.iter().map(|&k| names[k].as_str()).collect::<Vec<_>>().join("-");
        let best = if p.best.is_finite() { format!("{:.0} m/s", p.best * 1000.0) } else { "none yet".into() };
        let msg = format!(
            "{seq} · {} branches, {} pruned, {} sequences · {} Lambert arcs · best score {best}",
            p.nodes, p.pruned, p.sequences, p.lambert
        );
        let _ = progress.call1(&JsValue::NULL, &JsValue::from_str(&msg));
    };
    let only = fixed.then(|| flybys.iter().map(|&b| b as usize).chain([to]).collect());
    run(&prob, only, &mut report)
        .into_iter()
        .map(|r| FoundRoute {
            bodies: r.bodies.iter().map(|&b| b as u32).collect(),
            c3: r.cost.vinf_launch.powi(2),
            departure: r.cost.departure * 1000.0,
            flybys: r.cost.flybys.iter().map(|v| v * 1000.0).collect(),
            arrival: r.cost.arrival * 1000.0,
            dv: r.cost.dv() * 1000.0,
            score: r.cost.score() * 1000.0,
            dates: r.dates,
        })
        .collect()
}

/// Points along each leg of a route (`n + 1` per leg, legs one after another), interleaved x, y in AU.
#[wasm_bindgen]
pub fn route_xy(bodies: Vec<u32>, dates: Vec<f64>, n: usize) -> Vec<f64> {
    let ps = planets();
    let mut xy = Vec::new();
    for k in 0..bodies.len() - 1 {
        let (r1, _) = ps[bodies[k] as usize].state_at(dates[k]);
        let (r2, _) = ps[bodies[k + 1] as usize].state_at(dates[k + 1]);
        let tof = (dates[k + 1] - dates[k]) * DAY_S;
        match lambert(r1, r2, tof, MU_SUN_KM) {
            Some((v1, _)) => xy.extend(conic_xy(r1, v1, transfer_angle(r1, r2), MU_SUN_KM, n).into_iter().map(|c| c / AU_KM)),
            None => xy.extend(std::iter::repeat(f64::NAN).take(2 * (n + 1))),
        }
    }
    xy
}
