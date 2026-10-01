//! Lambert's problem: the conic arc between two positions in a given time.
//!
//! Universal-variable formulation (Bate, Mueller & White; Curtis, Algorithm 5.2), prograde and
//! single-revolution only, solved by bisection on z. Multi-revolution arcs will need Izzo (2015).

use std::f64::consts::{PI, TAU};

pub(crate) type V3 = [f64; 3];

pub(crate) fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(crate) fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
pub(crate) fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
pub(crate) fn scale(a: V3, k: f64) -> V3 {
    [a[0] * k, a[1] * k, a[2] * k]
}
pub(crate) fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Stumpff function C(z).
pub(crate) fn stumpff_c(z: f64) -> f64 {
    if z.abs() < 1e-3 {
        0.5 - z / 24.0 + z * z / 720.0
    } else if z > 0.0 {
        (1.0 - z.sqrt().cos()) / z
    } else {
        ((-z).sqrt().cosh() - 1.0) / -z
    }
}

/// Stumpff function S(z).
pub(crate) fn stumpff_s(z: f64) -> f64 {
    if z.abs() < 1e-3 {
        1.0 / 6.0 - z / 120.0 + z * z / 5040.0
    } else if z > 0.0 {
        let s = z.sqrt();
        (s - s.sin()) / s.powi(3)
    } else {
        let s = (-z).sqrt();
        (s.sinh() - s) / s.powi(3)
    }
}

/// Prograde transfer angle from `r1` to `r2` about the ecliptic north pole, in [0, 2π).
pub fn transfer_angle(r1: V3, r2: V3) -> f64 {
    let c = (dot(r1, r2) / (norm(r1) * norm(r2))).clamp(-1.0, 1.0);
    let dnu = c.acos();
    if cross(r1, r2)[2] < 0.0 { TAU - dnu } else { dnu }
}

/// Velocities `(v1, v2)` at each end of the prograde, single-revolution arc from `r1` to `r2` taking `tof`.
/// Any consistent units, e.g. km, s and km^3/s^2. `None` near a 180° (or 0°) transfer, where the
/// transfer plane is undefined, or if no solution is found.
pub fn lambert(r1: V3, r2: V3, tof: f64, mu: f64) -> Option<(V3, V3)> {
    let (n1, n2) = (norm(r1), norm(r2));
    let dnu = transfer_angle(r1, r2);
    let a = dnu.sin() * (n1 * n2 / (1.0 - dnu.cos())).sqrt();
    if !(a.abs() > 1e-6 * (n1 + n2)) || tof <= 0.0 {
        return None;
    }

    let y = |z: f64| n1 + n2 + a * (z * stumpff_s(z) - 1.0) / stumpff_c(z).sqrt();
    let sqrt_mu_t = mu.sqrt() * tof;
    // Time-of-flight residual. It increases with z; where y < 0 (only when a > 0) time is too short.
    let f = |z: f64| {
        let yz = y(z);
        if yz < 0.0 {
            return -sqrt_mu_t;
        }
        (yz / stumpff_c(z)).powf(1.5) * stumpff_s(z) + a * yz.sqrt() - sqrt_mu_t
    };

    // Bracket z between strongly hyperbolic and the single-revolution limit 4π².
    let mut hi = 4.0 * PI * PI - 1e-6;
    let mut lo = -4.0 * PI * PI;
    while f(lo) > 0.0 {
        lo *= 2.0;
        if lo < -1e5 {
            return None;
        }
    }
    if f(hi) < 0.0 {
        return None;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 { hi = mid } else { lo = mid }
        if hi - lo < 1e-12 * (1.0 + mid.abs()) {
            break;
        }
    }

    let z = 0.5 * (lo + hi);
    let yz = y(z);
    let lf = 1.0 - yz / n1;
    let g = a * (yz / mu).sqrt();
    let gdot = 1.0 - yz / n2;
    let v1 = scale(sub(r2, scale(r1, lf)), 1.0 / g);
    let v2 = scale(sub(scale(r2, gdot), r1), 1.0 / g);
    Some((v1, v2))
}

/// `n + 1` points along the conic through `r1` with velocity `v1`, sweeping `dnu` radians of true
/// anomaly. Interleaved x, y in the units of `r1`.
pub fn conic_xy(r1: V3, v1: V3, dnu: f64, mu: f64, n: usize) -> Vec<f64> {
    let h = cross(r1, v1);
    let hn = norm(h);
    let r1n = norm(r1);
    let e_vec = sub(scale(cross(v1, h), 1.0 / mu), scale(r1, 1.0 / r1n));
    let e = norm(e_vec);
    let p = hn * hn / mu;
    let nu1 = if e < 1e-12 { 0.0 } else { (dot(h, cross(e_vec, r1)) / hn).atan2(dot(e_vec, r1)) };
    let ph = scale(r1, 1.0 / r1n);
    let qh = scale(cross(h, ph), 1.0 / hn);
    (0..=n)
        .flat_map(|k| {
            let d = dnu * k as f64 / n as f64;
            let r = p / (1.0 + e * (nu1 + d).cos());
            let (s, c) = d.sin_cos();
            [r * (c * ph[0] + s * qh[0]), r * (c * ph[1] + s * qh[1])]
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curtis_example_5_2() {
        let (v1, v2) =
            lambert([5000.0, 10000.0, 2100.0], [-14600.0, 2500.0, 7000.0], 3600.0, 398_600.0).unwrap();
        let close = |a: V3, b: V3| (0..3).all(|i| (a[i] - b[i]).abs() < 1e-3);
        assert!(close(v1, [-5.9925, 1.9254, 3.2456]), "v1 = {v1:?}");
        assert!(close(v2, [-3.3125, -4.1966, -0.38529]), "v2 = {v2:?}");
    }

    #[test]
    fn conserves_energy_and_momentum() {
        let mu = 1.327e11;
        let r1 = [1.5e8, 0.0, 0.0];
        let r2 = [-1.0e8, 1.8e8, 1.0e7];
        for days in [60.0, 200.0, 500.0] {
            let (v1, v2) = lambert(r1, r2, days * 86_400.0, mu).unwrap();
            let energy = |r: V3, v: V3| dot(v, v) / 2.0 - mu / norm(r);
            assert!((energy(r1, v1) - energy(r2, v2)).abs() < 1e-6 * energy(r1, v1).abs());
            let (h1, h2) = (cross(r1, v1), cross(r2, v2));
            assert!(norm(sub(h1, h2)) < 1e-6 * norm(h1));
        }
    }

    #[test]
    fn conic_ends_at_r2() {
        let mu = 1.327e11;
        let (r1, r2) = ([1.5e8, 0.0, 0.0], [-1.0e8, 1.8e8, 0.0]);
        let (v1, _) = lambert(r1, r2, 200.0 * 86_400.0, mu).unwrap();
        let xy = conic_xy(r1, v1, transfer_angle(r1, r2), mu, 100);
        let n = xy.len();
        assert!((xy[n - 2] - r2[0]).abs() < 1.0 && (xy[n - 1] - r2[1]).abs() < 1.0);
    }

    #[test]
    fn rejects_180_degree_transfer() {
        assert!(lambert([1.0, 0.0, 0.0], [-2.0, 0.0, 0.0], 5.0, 1.0).is_none());
    }
}
