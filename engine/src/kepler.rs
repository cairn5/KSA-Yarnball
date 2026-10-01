//! Two-body (Kepler) propagation in universal variables (Curtis, Algorithms 3.3 and 3.4), so the
//! same code covers elliptical, parabolic and hyperbolic orbits, forwards or backwards in time.

use crate::lambert::{dot, norm, stumpff_c, stumpff_s, V3};

/// Position and velocity after `dt` from `r0`, `v0` around a body with parameter `mu`.
/// Any consistent units, e.g. km, km/s, s and km³/s².
pub fn propagate(r0: V3, v0: V3, dt: f64, mu: f64) -> (V3, V3) {
    if dt == 0.0 {
        return (r0, v0);
    }
    let r0n = norm(r0);
    let vr0 = dot(r0, v0) / r0n;
    let alpha = 2.0 / r0n - dot(v0, v0) / mu; // 1 / semi-major axis
    let sqrt_mu = mu.sqrt();

    // Universal Kepler equation F(x) = 0. F increases with x (its slope is the radius), so a
    // bracket plus safeguarded Newton always converges. Far out on hyperbolic orbits the terms
    // overflow; there F is taken as infinite with the sign of x, which keeps the bracket valid.
    let f = |x: f64| {
        let z = alpha * x * x;
        let v = r0n * vr0 / sqrt_mu * x * x * stumpff_c(z) + (1.0 - alpha * r0n) * x.powi(3) * stumpff_s(z) + r0n * x
            - sqrt_mu * dt;
        if v.is_finite() { v } else { f64::INFINITY.copysign(x) }
    };
    let df = |x: f64| {
        let z = alpha * x * x;
        r0n * vr0 / sqrt_mu * x * (1.0 - z * stumpff_s(z)) + (1.0 - alpha * r0n) * x * x * stumpff_c(z) + r0n
    };

    let guess = (sqrt_mu * alpha.abs() * dt.abs()).max(1e-3 * r0n.sqrt());
    let (mut lo, mut hi) = if dt > 0.0 { (0.0, guess) } else { (-guess, 0.0) };
    for _ in 0..200 {
        if dt > 0.0 && f(hi) < 0.0 {
            lo = hi;
            hi *= 2.0;
        } else if dt < 0.0 && f(lo) > 0.0 {
            hi = lo;
            lo *= 2.0;
        } else {
            break;
        }
    }
    // Newton where it at least halves the bracket, bisection otherwise (Numerical Recipes' rtsafe).
    // On strongly hyperbolic orbits F grows exponentially, so plain Newton from a high guess
    // creeps down far too slowly.
    let mut x = 0.5 * (lo + hi);
    let mut step_old = hi - lo;
    let mut step = step_old;
    for _ in 0..200 {
        let fx = f(x);
        if fx == 0.0 {
            break;
        }
        if fx > 0.0 { hi = x } else { lo = x }
        let d = df(x);
        let next = x - fx / d;
        let newton = d.is_finite() && next > lo && next < hi && (2.0 * fx).abs() < (step_old * d).abs();
        step_old = step;
        if newton {
            step = x - next;
            x = next;
        } else {
            step = 0.5 * (hi - lo);
            x = lo + step;
        }
        if step.abs() <= 1e-13 * (1.0 + x.abs()) {
            break;
        }
    }

    let z = alpha * x * x;
    let lf = 1.0 - x * x / r0n * stumpff_c(z);
    let g = dt - x.powi(3) / sqrt_mu * stumpff_s(z);
    let r = [lf * r0[0] + g * v0[0], lf * r0[1] + g * v0[1], lf * r0[2] + g * v0[2]];
    let rn = norm(r);
    let fdot = sqrt_mu / (rn * r0n) * (z * x * stumpff_s(z) - x);
    let gdot = 1.0 - x * x / rn * stumpff_c(z);
    let v = [fdot * r0[0] + gdot * v0[0], fdot * r0[1] + gdot * v0[1], fdot * r0[2] + gdot * v0[2]];
    (r, v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curtis_example_3_7() {
        // Earth orbit, 60 minutes: r = (-3297.8, 7413.4) km, v = (-8.2977, -0.96309) km/s.
        let (r, v) = propagate([7000.0, -12124.0, 0.0], [2.6679, 4.6210, 0.0], 3600.0, 398_600.0);
        assert!((r[0] + 3297.8).abs() < 0.5 && (r[1] - 7413.4).abs() < 0.5, "r = {r:?}");
        assert!((v[0] + 8.2977).abs() < 1e-3 && (v[1] + 0.96309).abs() < 2e-3, "v = {v:?}");
    }

    #[test]
    fn matches_pykep() {
        let data: serde_json::Value = serde_json::from_str(include_str!("../tests/data/pykep_reference.json")).unwrap();
        let v3 = |v: &serde_json::Value| -> V3 { [0, 1, 2].map(|i| v[i].as_f64().unwrap()) };
        let mu = crate::MU_SUN_KM;
        for c in data["kepler"].as_array().unwrap() {
            let (r, v) = propagate(v3(&c["r"]), v3(&c["v"]), c["dt"].as_f64().unwrap(), mu);
            let (r1, v1) = (v3(&c["r1"]), v3(&c["v1"]));
            assert!(norm(crate::lambert::sub(r, r1)) < 1e-8 * norm(r1), "r {r:?} vs {r1:?}");
            assert!(norm(crate::lambert::sub(v, v1)) < 1e-8 * norm(v1), "v {v:?} vs {v1:?}");
        }
    }

    #[test]
    fn round_trip() {
        let mu = 1.327e11;
        let (r0, v0) = ([1.4e8, 3.0e7, 1.0e6], [-5.0, 35.0, 1.0]);
        for dt in [1e5, 3e7, -2e7, 2e8] {
            let (r1, v1) = propagate(r0, v0, dt, mu);
            let (r2, v2) = propagate(r1, v1, -dt, mu);
            for i in 0..3 {
                assert!((r2[i] - r0[i]).abs() < 1e-3, "dt {dt}: {r2:?}");
                assert!((v2[i] - v0[i]).abs() < 1e-9, "dt {dt}: {v2:?}");
            }
        }
    }
}
