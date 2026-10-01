//! Trajectory engine. Headless maths only; the TypeScript interface calls it via WebAssembly.

use wasm_bindgen::prelude::*;

pub mod bodies;
pub mod lambert;
pub mod mga;
pub mod porkchop;

/// Sun's gravitational parameter, m^3/s^2.
pub const MU_SUN: f64 = 1.327_124_400_18e20;
/// Astronomical unit, m.
pub const AU: f64 = 1.495_978_707e11;
/// Sun's gravitational parameter, km^3/s^2.
pub const MU_SUN_KM: f64 = MU_SUN * 1e-9;
/// Astronomical unit, km.
pub const AU_KM: f64 = AU * 1e-3;
/// Seconds per day.
pub const DAY_S: f64 = 86_400.0;

/// Result of a two-impulse Hohmann transfer between circular, coplanar orbits.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy)]
pub struct Transfer {
    /// Departure burn, m/s.
    pub dv1: f64,
    /// Arrival burn, m/s.
    pub dv2: f64,
    /// Time of flight, days.
    pub tof_days: f64,
}

#[wasm_bindgen]
impl Transfer {
    #[wasm_bindgen(getter)]
    pub fn total(&self) -> f64 {
        self.dv1 + self.dv2
    }
}

/// Hohmann transfer around the Sun between circular orbits of radius `r1_au` and `r2_au`.
#[wasm_bindgen]
pub fn hohmann(r1_au: f64, r2_au: f64) -> Transfer {
    hohmann_mu(r1_au * AU, r2_au * AU, MU_SUN)
}

/// Hohmann transfer for any central body. Radii in m, `mu` in m^3/s^2.
pub fn hohmann_mu(r1: f64, r2: f64, mu: f64) -> Transfer {
    let a = 0.5 * (r1 + r2);
    let v1 = (mu / r1).sqrt();
    let v2 = (mu / r2).sqrt();
    let vp = (mu * (2.0 / r1 - 1.0 / a)).sqrt();
    let va = (mu * (2.0 / r2 - 1.0 / a)).sqrt();
    let tof = std::f64::consts::PI * (a.powi(3) / mu).sqrt();
    Transfer {
        dv1: (vp - v1).abs(),
        dv2: (v2 - va).abs(),
        tof_days: tof / 86_400.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn earth_to_mars() {
        // Textbook values: ~2.94 km/s + ~2.65 km/s, ~259 days.
        let t = hohmann(1.0, 1.524);
        assert!((t.dv1 - 2945.0).abs() < 20.0, "dv1 = {}", t.dv1);
        assert!((t.dv2 - 2649.0).abs() < 20.0, "dv2 = {}", t.dv2);
        assert!((t.tof_days - 259.0).abs() < 2.0, "tof = {}", t.tof_days);
    }
}
