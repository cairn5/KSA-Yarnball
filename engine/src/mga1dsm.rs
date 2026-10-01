//! MGA-1DSM: multiple gravity assists with one deep-space burn per leg (Izzo 2010, "Global
//! optimization and space pruning for spacecraft trajectory design").
//!
//! Each leg starts with the craft's velocity after launch or after an unpowered flyby, coasts
//! for a fraction η of the leg, makes one burn, then follows a Lambert arc to the next planet.
//! The chromosome follows pykep's `mga_1dsm` with direct time encoding, so the two can exchange
//! chromosomes and be checked against each other:
//!
//! `[t0, u, v, V∞, η1, T1] + [β, rp/R, η, T]` for each later leg
//!
//! t0 in pykep's mjd2000 (days from 2000-01-01 00:00), V∞ in m/s, T in days, β in radians.
//! Launch direction: θ = 2πu, φ = acos(2v − 1) − π/2 in the ecliptic frame.
//!
//! Costs match the route finder: launch burn from a 200 km parking orbit, the deep-space
//! burns, an optional capture burn, and the cost of flight time.

use std::f64::consts::{PI, TAU};

use wasm_bindgen::prelude::*;

use crate::bodies::{planets, Planet};
use crate::kepler::propagate;
use crate::lambert::{conic_xy, cross, dot, lambert, norm, scale, sub, transfer_angle, V3};
use crate::mga::{capture_dv, departure_dv, Problem};
use crate::sade::{Rng, Sade};
use crate::{AU_KM, DAY_S, MU_SUN_KM};

/// pykep's mjd2000 starts at 2000-01-01 00:00, half a day before the engine's J2000 dates.
pub const PYKEP_EPOCH_OFFSET: f64 = 0.5;
/// Where in each leg its burn may happen, as a fraction of the leg (as in the pykep runs).
pub const ETA_BOUNDS: (f64, f64) = (0.01, 0.9);
/// Highest flyby periapsis, planet radii (as in the pykep runs).
pub const RP_MAX: f64 = 200.0;
/// Launch v∞ range, km/s (as in the pykep runs: C3 6–144); widened if a seed needs more.
pub const VINF_BOUNDS_KMS: (f64, f64) = (2.5, 12.0);
/// Fitness of a chromosome that can't be evaluated, km/s.
const PENALTY: f64 = 1e6;

/// Population size, members started from the seed (1 exact + nudged copies), and generations per
/// round, as in the pykep runs (pygmo sade with gen = 500).
pub const POP: usize = 30;
pub const SEEDED: usize = 6;
pub const GEN_PER_ROUND: usize = 500;

fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn angle(a: V3, b: V3) -> f64 {
    (dot(a, b) / (norm(a) * norm(b))).clamp(-1.0, 1.0).acos()
}

/// Velocity after an unpowered flyby: the incoming v∞ (relative to the planet) is turned by the
/// angle a hyperbola with periapsis `rp` gives, in the plane at angle `beta` about it. `beta` is
/// measured from the direction normal to both v∞ and the planet's velocity (pykep's convention).
pub fn fb_prop(v_in: V3, v_pla: V3, rp: f64, beta: f64, mu: f64) -> V3 {
    let rel = sub(v_in, v_pla);
    let vinf = norm(rel);
    let e = 1.0 + rp * vinf * vinf / mu;
    let delta = 2.0 * (1.0 / e).asin();
    let b1 = scale(rel, 1.0 / vinf);
    let c = cross(b1, v_pla);
    let b2 = scale(c, 1.0 / norm(c));
    let b3 = cross(b1, b2);
    let (sd, cd) = delta.sin_cos();
    let (sb, cb) = beta.sin_cos();
    let out = [0, 1, 2].map(|i| vinf * (cd * b1[i] + cb * sd * b2[i] + sb * sd * b3[i]));
    add(v_pla, out)
}

/// Cost breakdown of one MGA-1DSM trajectory, km/s.
#[derive(Clone, Debug)]
pub struct Breakdown {
    /// Encounter dates, days after J2000.
    pub dates: Vec<f64>,
    /// Burn dates, one per leg.
    pub dsm_dates: Vec<f64>,
    pub dsms: Vec<f64>,
    pub vinf_launch: f64,
    pub vinf_arrive: f64,
    pub departure: f64,
    pub arrival: f64,
    pub time: f64,
}

impl Breakdown {
    pub fn dv(&self) -> f64 {
        self.departure + self.dsms.iter().sum::<f64>() + self.arrival
    }

    pub fn score(&self) -> f64 {
        self.dv() + self.time
    }
}

/// One leg as flown: from the planet, coast, burn, Lambert arc to the next planet. km, km/s, s.
struct Leg {
    r0: V3,
    v0: V3,
    coast: f64,
    r_dsm: V3,
    v_arc: V3,
    r_end: V3,
}

pub struct Mga1Dsm {
    pub bodies: Vec<usize>,
    ps: Vec<Planet>,
    pub lb: Vec<f64>,
    pub ub: Vec<f64>,
    capture: bool,
    /// km/s per day.
    time_cost: f64,
}

impl Mga1Dsm {
    /// The model for `bodies` with the route finder's launch window, leg limits and costs.
    /// `vinf_max` (km/s) raises the launch v∞ bound if a seed needs more.
    pub fn new(prob: &Problem, bodies: &[usize], vinf_max: f64) -> Self {
        let ps = planets();
        let n = bodies.len() - 1;
        let tof = |k: usize| prob.leg_bounds(&ps, bodies[k], bodies[k + 1]).unwrap_or((prob.step, prob.max_tof));
        let (t0_lo, t0_hi) = (prob.launch_start + PYKEP_EPOCH_OFFSET, prob.launch_end + PYKEP_EPOCH_OFFSET);
        let vinf_hi = VINF_BOUNDS_KMS.1.max(1.2 * vinf_max) * 1000.0;
        let mut lb = vec![t0_lo, 0.0, 0.0, VINF_BOUNDS_KMS.0 * 1000.0, ETA_BOUNDS.0, tof(0).0];
        let mut ub = vec![t0_hi, 1.0, 1.0, vinf_hi, ETA_BOUNDS.1, tof(0).1];
        for k in 1..n {
            let (lo, hi) = tof(k);
            lb.extend([-TAU, ps[bodies[k]].safe_radius, ETA_BOUNDS.0, lo]);
            ub.extend([TAU, RP_MAX, ETA_BOUNDS.1, hi]);
        }
        Mga1Dsm { bodies: bodies.to_vec(), ps, lb, ub, capture: prob.capture, time_cost: prob.time_cost }
    }

    /// The model with no bounds or costs beyond the burns, for comparing with pykep.
    #[cfg(test)]
    fn bare(bodies: &[usize]) -> Self {
        Mga1Dsm { bodies: bodies.to_vec(), ps: planets(), lb: vec![], ub: vec![], capture: false, time_cost: 0.0 }
    }

    fn legs(&self, x: &[f64]) -> Option<(Breakdown, Vec<Leg>)> {
        let n = self.bodies.len() - 1;
        let t0 = x[0] - PYKEP_EPOCH_OFFSET;
        let tofs: Vec<f64> = (0..n).map(|k| x[5 + 4 * k]).collect();
        let mut dates = vec![t0];
        for &t in &tofs {
            dates.push(dates.last().unwrap() + t);
        }
        let states: Vec<(V3, V3)> = self.bodies.iter().zip(&dates).map(|(&b, &t)| self.ps[b].state_at(t)).collect();

        let theta = TAU * x[1];
        let phi = (2.0 * x[2] - 1.0).clamp(-1.0, 1.0).acos() - PI / 2.0;
        let vinf = x[3] / 1000.0;
        let vinf_vec = [vinf * phi.cos() * theta.cos(), vinf * phi.cos() * theta.sin(), vinf * phi.sin()];

        let mut legs = Vec::with_capacity(n);
        let mut dsms = Vec::with_capacity(n);
        let mut dsm_dates = Vec::with_capacity(n);
        let mut v_start = add(states[0].1, vinf_vec);
        let mut v_end = [0.0; 3];
        for k in 0..n {
            if k > 0 {
                let p = &self.ps[self.bodies[k]];
                let (rp, beta) = (x[7 + 4 * (k - 1)] * p.radius_km, x[6 + 4 * (k - 1)]);
                v_start = fb_prop(v_end, states[k].1, rp, beta, p.mu);
            }
            let eta = if k == 0 { x[4] } else { x[8 + 4 * (k - 1)] };
            let coast = eta * tofs[k] * DAY_S;
            let (r, v) = propagate(states[k].0, v_start, coast, MU_SUN_KM);
            let (v1, v2) = lambert(r, states[k + 1].0, (1.0 - eta) * tofs[k] * DAY_S, MU_SUN_KM)?;
            dsms.push(norm(sub(v1, v)));
            dsm_dates.push(dates[k] + eta * tofs[k]);
            legs.push(Leg { r0: states[k].0, v0: v_start, coast, r_dsm: r, v_arc: v1, r_end: states[k + 1].0 });
            v_end = v2;
        }

        let vinf_arrive = norm(sub(v_end, states[n].1));
        let target = &self.ps[self.bodies[n]];
        let b = Breakdown {
            departure: departure_dv(&self.ps[self.bodies[0]], vinf),
            arrival: if self.capture { capture_dv(target, vinf_arrive) } else { 0.0 },
            time: self.time_cost * (dates[n] - dates[0]),
            vinf_launch: vinf,
            vinf_arrive,
            dates,
            dsm_dates,
            dsms,
        };
        let finite = b.dsms.iter().all(|d| d.is_finite()) && b.arrival.is_finite();
        finite.then_some((b, legs))
    }

    pub fn evaluate(&self, x: &[f64]) -> Option<Breakdown> {
        self.legs(x).map(|(b, _)| b)
    }

    /// What SADE minimises: Δv plus the flight-time cost, km/s.
    pub fn fitness(&self, x: &[f64]) -> f64 {
        self.evaluate(x).map_or(PENALTY, |b| b.score())
    }

    /// Points along each leg (`2 (n + 1)` per leg: the coast, then the arc after the burn) and
    /// the burn positions, both interleaved x, y in AU.
    pub fn trace(&self, x: &[f64], n: usize) -> Option<(Vec<f64>, Vec<f64>)> {
        let (_, legs) = self.legs(x)?;
        let mut xy = Vec::new();
        let mut dsm = Vec::new();
        for leg in &legs {
            for k in 0..=n {
                let (r, _) = propagate(leg.r0, leg.v0, leg.coast * k as f64 / n as f64, MU_SUN_KM);
                xy.extend([r[0] / AU_KM, r[1] / AU_KM]);
            }
            let arc = conic_xy(leg.r_dsm, leg.v_arc, transfer_angle(leg.r_dsm, leg.r_end), MU_SUN_KM, n);
            xy.extend(arc.into_iter().map(|c| c / AU_KM));
            dsm.extend([leg.r_dsm[0] / AU_KM, leg.r_dsm[1] / AU_KM]);
        }
        Some((xy, dsm))
    }

    /// A chromosome that follows a flyby-only route with encounter `dates` (days after J2000),
    /// traced leg by leg as the model builds it, so each flyby is set from the actual incoming
    /// velocity. Where the route needed a powered flyby, the unpowered one can't match it exactly;
    /// the leg's burn then goes just after the flyby, where it does the same job.
    /// (Ported from `to_chromosome` in the pykep prototype `cassini_seeded_test.py`.)
    pub fn seed_from_dates(&self, dates: &[f64]) -> Vec<f64> {
        let n = self.bodies.len() - 1;
        let states: Vec<(V3, V3)> = self.bodies.iter().zip(dates).map(|(&b, &t)| self.ps[b].state_at(t)).collect();
        let tofs: Vec<f64> = dates.windows(2).map(|w| w[1] - w[0]).collect();
        // Where each flyby should send the craft: the flyby-only route's departure velocities.
        let v_dep: Vec<V3> = (0..n)
            .map(|k| {
                lambert(states[k].0, states[k + 1].0, tofs[k] * DAY_S, MU_SUN_KM)
                    .map_or(states[k].1, |(v1, _)| v1)
            })
            .collect();
        let eta_early = ETA_BOUNDS.0;

        // Launch: v∞ magnitude and direction in the (u, v) encoding.
        let vinf = sub(v_dep[0], states[0].1);
        let vmag = norm(vinf);
        let theta = vinf[1].atan2(vinf[0]).rem_euclid(TAU);
        let phi = (vinf[2] / vmag).asin();
        let eta0 = 0.5; // on the flyby-only arc the first burn is zero wherever it is
        let mut x = vec![dates[0] + PYKEP_EPOCH_OFFSET, theta / TAU, (1.0 - phi.sin()) / 2.0, vmag * 1000.0, eta0, tofs[0]];
        let leg_end = |r0: V3, v0: V3, eta: f64, k: usize| {
            let (r, v) = propagate(r0, v0, eta * tofs[k] * DAY_S, MU_SUN_KM);
            lambert(r, states[k + 1].0, (1.0 - eta) * tofs[k] * DAY_S, MU_SUN_KM).map_or(v, |(_, v2)| v2)
        };
        let mut v_end = leg_end(states[0].0, v_dep[0], eta0, 0);

        // Flybys: the turn needed sets the periapsis; then search the plane angle β that points
        // the outgoing v∞ along the route's next leg.
        for k in 1..n {
            let p = &self.ps[self.bodies[k]];
            let v_pla = states[k].1;
            let vin = sub(v_end, v_pla);
            let want = sub(v_dep[k], v_pla);
            let turn = angle(vin, want);
            let (rp_lo, rp_hi) = (p.safe_radius * p.radius_km, RP_MAX * p.radius_km);
            let rp = if turn < 1e-9 { rp_hi } else { (1.0 / (turn / 2.0).sin() - 1.0) * p.mu / dot(vin, vin) };
            let rp = rp.clamp(rp_lo, rp_hi);
            let miss = |beta: f64| angle(sub(fb_prop(v_end, v_pla, rp, beta, p.mu), v_pla), want);
            let best_of = |betas: &mut dyn Iterator<Item = f64>| {
                betas.map(|b| (miss(b), b)).min_by(|a, b| a.0.total_cmp(&b.0)).unwrap().1
            };
            let coarse = best_of(&mut (0..=720).map(|i| -PI + TAU * i as f64 / 720.0));
            let beta = best_of(&mut (0..=200).map(|i| coarse - 0.01 + 0.02 * i as f64 / 200.0));
            x.extend([beta, rp / p.radius_km, eta_early, tofs[k]]);

            let v_out = fb_prop(v_end, v_pla, rp, beta, p.mu);
            v_end = leg_end(states[k].0, v_out, eta_early, k);
        }
        self.clip(x)
    }

    fn clip(&self, x: Vec<f64>) -> Vec<f64> {
        x.into_iter().zip(self.lb.iter().zip(&self.ub)).map(|(v, (&l, &u))| v.clamp(l, u)).collect()
    }

    /// A seed with small random changes to every variable (scales as in the pykep prototype).
    pub fn nudge(&self, x: &[f64], rng: &mut Rng) -> Vec<f64> {
        let mut s = vec![3.0, 0.002, 0.002, 0.01 * x[3], 0.02, 3.0];
        for k in 1..self.bodies.len() - 1 {
            s.extend([0.05, 0.05 * x[7 + 4 * (k - 1)], 0.02, 3.0]);
        }
        self.clip(x.iter().zip(&s).map(|(v, s)| v + s * rng.normal()).collect())
    }
}

/// Seeded SADE on the MGA-1DSM model, started from a route-finder route. Run by the interface a
/// few generations at a time (one per Web Worker, as independent islands), so it can show
/// progress. Costs reported in m/s.
#[wasm_bindgen]
pub struct Refiner {
    model: Mga1Dsm,
    sade: Sade,
    seed: Vec<f64>,
    round_generation: usize,
    rounds: usize,
}

#[wasm_bindgen]
impl Refiner {
    /// `bodies` and `dates` (days after J2000) are a route-finder route; the other arguments are
    /// the search settings it was found with (`time_cost` in m/s per year). `seed` sets the
    /// island's random stream.
    #[wasm_bindgen(constructor)]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bodies: Vec<u32>,
        dates: Vec<f64>,
        launch_start: f64,
        launch_end: f64,
        capture: bool,
        time_cost: f64,
        max_tof: f64,
        step: f64,
        seed: u64,
    ) -> Refiner {
        let bodies: Vec<usize> = bodies.iter().map(|&b| b as usize).collect();
        let prob = Problem {
            from: bodies[0],
            to: *bodies.last().unwrap(),
            launch_start: launch_start.min(dates[0]),
            launch_end: launch_end.max(dates[0]),
            max_flybys: bodies.len() - 2,
            capture,
            time_cost: time_cost / 1000.0 / 365.25,
            max_tof,
            step,
            margin: 0.0,
        };
        // Find the seed's launch v∞ first, so the bounds can admit it.
        let probe = Mga1Dsm::new(&prob, &bodies, 0.0);
        let vinf = probe.seed_from_dates(&dates)[3] / 1000.0;
        let mut model = Mga1Dsm::new(&prob, &bodies, vinf);
        // Leg limits must admit the route's own legs.
        for (k, w) in dates.windows(2).enumerate() {
            let j = 5 + 4 * k;
            model.lb[j] = model.lb[j].min(w[1] - w[0]);
            model.ub[j] = model.ub[j].max(w[1] - w[0]);
        }
        let x0 = model.seed_from_dates(&dates);

        let mut rng = Rng::new(seed);
        let mut pop = vec![x0.clone()];
        while pop.len() < SEEDED {
            pop.push(model.nudge(&x0, &mut rng));
        }
        while pop.len() < POP {
            pop.push(model.lb.iter().zip(&model.ub).map(|(l, u)| l + rng.f64() * (u - l)).collect());
        }
        let sade = Sade::new(pop, model.lb.clone(), model.ub.clone(), rng, &mut |x| model.fitness(x));
        Refiner { model, sade, seed: x0, round_generation: 0, rounds: 0 }
    }

    /// Runs up to `generations` more generations. A round ends after 500 generations or when the
    /// population collapses; the next one redraws every F and CR, as pygmo does between evolves.
    pub fn evolve(&mut self, generations: usize) {
        let model = &self.model;
        let mut f = |x: &[f64]| model.fitness(x);
        for _ in 0..generations {
            self.sade.step(&mut f);
            self.round_generation += 1;
            if self.round_generation >= GEN_PER_ROUND || self.sade.collapsed() {
                self.rounds += 1;
                self.round_generation = 0;
                self.sade.new_round();
            }
        }
    }

    #[wasm_bindgen(getter)]
    pub fn generation(&self) -> usize {
        self.sade.generation
    }

    #[wasm_bindgen(getter)]
    pub fn evaluations(&self) -> usize {
        self.sade.evaluations
    }

    /// Rounds completed.
    #[wasm_bindgen(getter)]
    pub fn rounds(&self) -> usize {
        self.rounds
    }

    /// Best score so far (Δv plus flight-time cost), m/s.
    #[wasm_bindgen(getter)]
    pub fn best_score(&self) -> f64 {
        self.sade.best().1 * 1000.0
    }

    /// The seed converted from the route, as a chromosome.
    #[wasm_bindgen(getter)]
    pub fn seed(&self) -> Vec<f64> {
        self.seed.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn best(&self) -> Vec<f64> {
        self.sade.best().0.to_vec()
    }

    /// Breakdown and drawing of chromosome `x`, with `n + 1` points per coast and per arc.
    pub fn describe(&self, x: Vec<f64>, n: usize) -> Option<Trajectory> {
        let b = self.model.evaluate(&x)?;
        let (xy, dsm_xy) = self.model.trace(&x, n)?;
        Some(Trajectory {
            bodies: self.model.bodies.iter().map(|&b| b as u32).collect(),
            c3: b.vinf_launch.powi(2),
            departure: b.departure * 1000.0,
            dsms: b.dsms.iter().map(|d| d * 1000.0).collect(),
            arrival: b.arrival * 1000.0,
            dv: b.dv() * 1000.0,
            score: b.score() * 1000.0,
            dates: b.dates,
            dsm_dates: b.dsm_dates,
            xy,
            dsm_xy,
            chromosome: x,
        })
    }
}

/// An MGA-1DSM trajectory for the interface. Costs in m/s, dates in days after J2000, drawing in AU.
#[wasm_bindgen(getter_with_clone)]
pub struct Trajectory {
    pub bodies: Vec<u32>,
    pub dates: Vec<f64>,
    pub dsm_dates: Vec<f64>,
    pub dsms: Vec<f64>,
    pub c3: f64,
    pub departure: f64,
    pub arrival: f64,
    pub dv: f64,
    pub score: f64,
    /// Each leg's coast then arc, interleaved x, y.
    pub xy: Vec<f64>,
    /// One burn position per leg, interleaved x, y.
    pub dsm_xy: Vec<f64>,
    /// pykep-compatible chromosome.
    pub chromosome: Vec<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mga::search_sequence;

    #[test]
    fn seed_follows_a_flyby_only_route() {
        // Cassini's 1997 window from the grid search: converting it to a chromosome lands close
        // to the grid route's cost (the pykep runs saw 0.2–2%).
        let prob = Problem {
            from: 2,
            to: 5,
            launch_start: -1095.5,
            launch_end: -548.5,
            max_flybys: 4,
            capture: true,
            time_cost: 0.0,
            max_tof: 3650.0,
            step: 5.0,
            margin: 0.0,
        };
        let seq = [1, 1, 2, 4, 5];
        let route = &search_sequence(&prob, &seq)[0];
        let bodies: Vec<usize> = std::iter::once(2).chain(seq).collect();
        let model = Mga1Dsm::new(&prob, &bodies, 0.0);
        let x = model.seed_from_dates(&route.dates);
        let seeded = model.evaluate(&x).unwrap();
        let rel = (seeded.score() - route.cost.score()) / route.cost.score();
        assert!(rel.abs() < 0.05, "seed {} vs grid {}", seeded.score(), route.cost.score());
        // Encounter dates are kept.
        for (a, b) in seeded.dates.iter().zip(&route.dates) {
            assert!((a - b).abs() < 1e-9);
        }
    }

    #[test]
    fn refiner_improves_on_its_seed() {
        let mut r = Refiner::new(vec![2, 4, 5], vec![10_200.0, 10_700.0, 13_200.0], 9_800.0, 10_600.0, true, 0.0, 4_000.0, 10.0, 1);
        let start = r.model.fitness(&r.seed) * 1000.0;
        r.evolve(300);
        assert!(r.best_score() <= start + 1e-9);
        assert!(r.best_score() < start, "no improvement from {start}");
        let t = r.describe(r.best(), 20).unwrap();
        assert_eq!(t.xy.len(), 2 * 2 * 21 * 2);
        assert_eq!(t.dsm_xy.len(), 4);
    }

    /// The model's burns match pykep's mga_1dsm on the reference chromosomes.
    #[test]
    fn matches_pykep() {
        let data: serde_json::Value = serde_json::from_str(include_str!("../tests/data/pykep_reference.json")).unwrap();
        let names: Vec<String> = data["sequence"].as_array().unwrap().iter().map(|s| s.as_str().unwrap().to_string()).collect();
        let ps = planets();
        let bodies: Vec<usize> = names.iter().map(|n| ps.iter().position(|p| p.name.eq_ignore_ascii_case(n)).unwrap()).collect();
        let model = Mga1Dsm::bare(&bodies);
        let (mut checked, mut worst) = (0, 0.0f64);
        for case in data["mga_1dsm"].as_array().unwrap() {
            let x: Vec<f64> = case["x"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
            let want: Vec<f64> = case["dv"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
            // Near-180° Lambert arcs are rejected here but solved by pykep; skip those.
            let Some(b) = model.evaluate(&x) else { continue };
            let got = b.dsms.iter().chain([&b.vinf_arrive]);
            for (g, w) in got.zip(&want) {
                worst = worst.max((g - w).abs() / w.max(1.0));
            }
            checked += 1;
        }
        assert!(checked > 250, "only {checked} cases evaluated");
        assert!(worst < 1e-6, "worst relative burn error {worst:e}");
    }

    #[test]
    fn flyby_matches_pykep() {
        let data: serde_json::Value = serde_json::from_str(include_str!("../tests/data/pykep_reference.json")).unwrap();
        let venus = &planets()[1];
        let v3 = |v: &serde_json::Value| -> V3 { [0, 1, 2].map(|i| v[i].as_f64().unwrap()) };
        for case in data["flyby"].as_array().unwrap() {
            let out = fb_prop(v3(&case["v_in"]), v3(&case["v_pla"]), case["rp"].as_f64().unwrap(), case["beta"].as_f64().unwrap(), venus.mu);
            let want = v3(&case["v_out"]);
            assert!(norm(sub(out, want)) < 1e-9 * norm(want), "{out:?} vs {want:?}");
        }
    }
}
