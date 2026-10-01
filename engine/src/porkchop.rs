//! Porkchop plots: Lambert transfers between two planets over a grid of departure dates and flight times.

use wasm_bindgen::prelude::*;

use crate::bodies::Planet;
use crate::lambert::{conic_xy, lambert, transfer_angle};
use crate::{AU_KM, DAY_S, MU_SUN_KM};

/// Launch and arrival excess speeds over a grid. Row-major: `index = tof_index * n_dep + dep_index`.
/// Cells with no solution are NaN.
#[wasm_bindgen]
pub struct Porkchop {
    c3: Vec<f64>,
    vinf_arr: Vec<f64>,
}

#[wasm_bindgen]
impl Porkchop {
    /// Launch energy C3 = v∞², km²/s².
    pub fn c3(&self) -> Vec<f64> {
        self.c3.clone()
    }

    /// Arrival excess speed v∞, km/s.
    pub fn vinf_arr(&self) -> Vec<f64> {
        self.vinf_arr.clone()
    }
}

/// One Lambert transfer between two planets.
#[wasm_bindgen(getter_with_clone)]
pub struct Leg {
    /// Launch energy, km²/s².
    pub c3: f64,
    /// Arrival excess speed, km/s.
    pub vinf_arr: f64,
    /// Points along the arc, interleaved x, y in AU (ecliptic plan view).
    pub xy: Vec<f64>,
}

fn diff_norm(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Departure v∞ and arrival v∞ (km/s) for leaving `from` at `dep` and reaching `to` `tof` days later.
fn solve(from: &Planet, to: &Planet, dep: f64, tof: f64) -> Option<(f64, f64)> {
    let (r1, vp1) = from.state_at(dep);
    let (r2, vp2) = to.state_at(dep + tof);
    let (v1, v2) = lambert(r1, r2, tof * DAY_S, MU_SUN_KM)?;
    Some((diff_norm(v1, vp1), diff_norm(v2, vp2)))
}

/// Grid of transfers from `from` to `to`: `n_dep` departure dates across [dep_start, dep_end] and
/// `n_tof` flight times across [tof_min, tof_max]. Dates are days after J2000.
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn porkchop(
    from: &Planet,
    to: &Planet,
    dep_start: f64,
    dep_end: f64,
    n_dep: usize,
    tof_min: f64,
    tof_max: f64,
    n_tof: usize,
) -> Porkchop {
    let step = |lo: f64, hi: f64, n: usize, i: usize| if n > 1 { lo + (hi - lo) * i as f64 / (n - 1) as f64 } else { lo };
    let mut c3 = Vec::with_capacity(n_dep * n_tof);
    let mut vinf_arr = Vec::with_capacity(n_dep * n_tof);
    for j in 0..n_tof {
        let tof = step(tof_min, tof_max, n_tof, j);
        for i in 0..n_dep {
            let dep = step(dep_start, dep_end, n_dep, i);
            let (vd, va) = solve(from, to, dep, tof).unwrap_or((f64::NAN, f64::NAN));
            c3.push(vd * vd);
            vinf_arr.push(va);
        }
    }
    Porkchop { c3, vinf_arr }
}

/// The transfer leaving `from` at `dep` (days after J2000) and arriving at `to` `tof` days later,
/// with `n + 1` points along it.
#[wasm_bindgen]
pub fn lambert_leg(from: &Planet, to: &Planet, dep: f64, tof: f64, n: usize) -> Option<Leg> {
    let (r1, vp1) = from.state_at(dep);
    let (r2, vp2) = to.state_at(dep + tof);
    let (v1, v2) = lambert(r1, r2, tof * DAY_S, MU_SUN_KM)?;
    let xy = conic_xy(r1, v1, transfer_angle(r1, r2), MU_SUN_KM, n)
        .into_iter()
        .map(|c| c / AU_KM)
        .collect();
    Some(Leg { c3: diff_norm(v1, vp1).powi(2), vinf_arr: diff_norm(v2, vp2), xy })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::planets;

    #[test]
    fn earth_mars_2026_window() {
        // The late-2026 Earth–Mars window has a minimum C3 of roughly 9–12 km²/s².
        let p = planets();
        let (n_dep, n_tof) = (120, 80);
        let dep_start = 9_600.0; // 2026-04-14
        let grid = porkchop(&p[2], &p[3], dep_start, dep_start + 360.0, n_dep, 120.0, 400.0, n_tof);
        let best = grid.c3.iter().copied().filter(|v| v.is_finite()).fold(f64::MAX, f64::min);
        assert!(best > 7.0 && best < 14.0, "best C3 = {best}");
        let k = grid.c3.iter().position(|&v| v == best).unwrap();
        let dep = dep_start + 360.0 * (k % n_dep) as f64 / (n_dep - 1) as f64;
        // Best departure falls between September 2026 and January 2027.
        assert!(dep > 9_740.0 && dep < 9_860.0, "dep = {dep}");
    }
}
