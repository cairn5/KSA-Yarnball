//! Planets and their orbital elements.

use std::f64::consts::TAU;

use wasm_bindgen::prelude::*;

use crate::{AU_KM, MU_SUN_KM};

/// A planet's Keplerian elements relative to the mean ecliptic and equinox of J2000.
/// Public fields are the J2000 values; `rates` move them linearly with time.
#[wasm_bindgen(getter_with_clone)]
#[derive(Debug, Clone)]
pub struct Planet {
    pub name: String,
    /// Semi-major axis, AU.
    pub a_au: f64,
    /// Eccentricity.
    pub e: f64,
    /// Inclination, degrees.
    pub i_deg: f64,
    /// Longitude of the ascending node, degrees.
    pub raan_deg: f64,
    /// Longitude of perihelion, degrees.
    pub lon_peri_deg: f64,
    /// Mean longitude, degrees.
    pub mean_lon_deg: f64,
    /// Gravitational parameter, km³/s².
    pub mu: f64,
    /// Equatorial radius, km.
    pub radius_km: f64,
    /// Lowest allowed flyby periapsis, in planet radii.
    pub safe_radius: f64,
    /// Rates of the six elements above, per Julian century.
    rates: [f64; 6],
}

#[wasm_bindgen]
impl Planet {
    /// Orbital period, days.
    #[wasm_bindgen(getter)]
    pub fn period_days(&self) -> f64 {
        365.25 * self.a_au.powf(1.5)
    }

    /// Heliocentric ecliptic position `[x, y, z]` in AU, `days` after J2000 (2000-01-01 12:00 TT).
    pub fn position_at(&self, days: f64) -> Vec<f64> {
        self.state_at(days).0.map(|c| c / AU_KM).to_vec()
    }

    /// `n + 1` points around the J2000 orbit (first repeated at the end), interleaved x, y in AU.
    pub fn orbit_xy(&self, n: usize) -> Vec<f64> {
        let el = self.elements_at(0.0);
        (0..=n)
            .flat_map(|k| {
                let ([x, y, _], _) = el.state(TAU * k as f64 / n as f64);
                [x / AU_KM, y / AU_KM]
            })
            .collect()
    }
}

impl Planet {
    /// Heliocentric ecliptic position (km) and velocity (km/s), `days` after J2000.
    pub fn state_at(&self, days: f64) -> ([f64; 3], [f64; 3]) {
        let el = self.elements_at(days);
        el.state(solve_kepler(el.m, el.e))
    }

    fn elements_at(&self, days: f64) -> Elements {
        let t = days / 36_525.0;
        let r = &self.rates;
        let raan = self.raan_deg + r[3] * t;
        let lon_peri = self.lon_peri_deg + r[4] * t;
        let mean_lon = self.mean_lon_deg + r[5] * t;
        Elements {
            a: (self.a_au + r[0] * t) * AU_KM,
            e: self.e + r[1] * t,
            i: (self.i_deg + r[2] * t).to_radians(),
            raan: raan.to_radians(),
            argp: (lon_peri - raan).to_radians(),
            m: (mean_lon - lon_peri).to_radians().rem_euclid(TAU),
        }
    }
}

/// Elliptical elements in km and radians.
struct Elements {
    a: f64,
    e: f64,
    i: f64,
    raan: f64,
    argp: f64,
    m: f64,
}

impl Elements {
    /// Position (km) and velocity (km/s) at eccentric anomaly `ea`.
    fn state(&self, ea: f64) -> ([f64; 3], [f64; 3]) {
        let (a, e) = (self.a, self.e);
        let b = a * (1.0 - e * e).sqrt();
        let (se, ce) = ea.sin_cos();
        let ea_dot = (MU_SUN_KM / a.powi(3)).sqrt() / (1.0 - e * ce);
        // Orbital plane, x towards perihelion.
        let (xp, yp) = (a * (ce - e), b * se);
        let (vxp, vyp) = (-a * se * ea_dot, b * ce * ea_dot);

        let (so, co) = self.raan.sin_cos();
        let (sw, cw) = self.argp.sin_cos();
        let (si, ci) = self.i.sin_cos();
        let p = [co * cw - so * sw * ci, so * cw + co * sw * ci, sw * si];
        let q = [-co * sw - so * cw * ci, -so * sw + co * cw * ci, cw * si];
        let rot = |x: f64, y: f64| [p[0] * x + q[0] * y, p[1] * x + q[1] * y, p[2] * x + q[2] * y];
        (rot(xp, yp), rot(vxp, vyp))
    }
}

/// Solves Kepler's equation `M = E - e sin E` for the eccentric anomaly E (elliptical orbits).
pub fn solve_kepler(m: f64, e: f64) -> f64 {
    let mut ea = if e < 0.8 { m } else { std::f64::consts::PI };
    for _ in 0..50 {
        let d = (ea - e * ea.sin() - m) / (1.0 - e * ea.cos());
        ea -= d;
        if d.abs() < 1e-14 {
            break;
        }
    }
    ea
}

/// JPL approximate planetary elements, valid 1800–2050 (Standish, Table 1). Earth is the Earth–Moon barycentre.
/// Order: a (AU), e, i, Ω, ϖ, L (degrees); then the same per Julian century.
const ELEMENTS: [(&str, [f64; 6], [f64; 6]); 8] = [
    ("Mercury",
        [0.38709927, 0.20563593, 7.00497902, 48.33076593, 77.45779628, 252.25032350],
        [0.00000037, 0.00001906, -0.00594749, -0.12534081, 0.16047689, 149472.67411175]),
    ("Venus",
        [0.72333566, 0.00677672, 3.39467605, 76.67984255, 131.60246718, 181.97909950],
        [0.00000390, -0.00004107, -0.00078890, -0.27769418, 0.00268329, 58517.81538729]),
    ("Earth",
        [1.00000261, 0.01671123, -0.00001531, 0.0, 102.93768193, 100.46457166],
        [0.00000562, -0.00004392, -0.01294668, 0.0, 0.32327364, 35999.37244981]),
    ("Mars",
        [1.52371034, 0.09339410, 1.84969142, 49.55953891, -23.94362959, -4.55343205],
        [0.00001847, 0.00007882, -0.00813131, -0.29257343, 0.44441088, 19140.30268499]),
    ("Jupiter",
        [5.20288700, 0.04838624, 1.30439695, 100.47390909, 14.72847983, 34.39644051],
        [-0.00011607, -0.00013253, -0.00183714, 0.20469106, 0.21252668, 3034.74612775]),
    ("Saturn",
        [9.53667594, 0.05386179, 2.48599187, 113.66242448, 92.59887831, 49.95424423],
        [-0.00125060, -0.00050991, 0.00193609, -0.28867794, -0.41897216, 1222.49362201]),
    ("Uranus",
        [19.18916464, 0.04725744, 0.77263783, 74.01692503, 170.95427630, 313.23810451],
        [-0.00196176, -0.00004397, -0.00242939, 0.04240589, 0.40805281, 428.48202785]),
    ("Neptune",
        [30.06992276, 0.00859048, 1.77004347, 131.78422574, 44.96476227, -55.12002969],
        [0.00026291, 0.00005105, 0.00035372, -0.00508664, -0.32241464, 218.45945325]),
];

/// Gravitational parameter (km³/s²), radius (km) and lowest flyby periapsis (planet radii), as in
/// pykep's jpl_lp planets. Flyby limits follow the Cassini and Voyager prototypes (rings at Saturn).
const PHYSICAL: [[f64; 3]; 8] = [
    [22_032.0, 2_440.0, 1.05],
    [324_859.0, 6_052.0, 1.05],
    [398_600.4418, 6_378.0, 1.15],
    [42_828.0, 3_397.0, 1.05],
    [126_686_534.0, 71_492.0, 1.7],
    [37_931_187.0, 60_330.0, 2.3],
    [5_793_939.0, 25_362.0, 2.1],
    [6_836_529.0, 24_622.0, 1.05],
];

/// The eight planets, Sun outwards.
#[wasm_bindgen]
pub fn planets() -> Vec<Planet> {
    ELEMENTS
        .iter()
        .zip(PHYSICAL)
        .map(|(&(name, [a_au, e, i_deg, raan_deg, lon_peri_deg, mean_lon_deg], rates), [mu, radius_km, safe_radius])| Planet {
            name: name.to_string(),
            a_au,
            e,
            i_deg,
            raan_deg,
            lon_peri_deg,
            mean_lon_deg,
            mu,
            radius_km,
            safe_radius,
            rates,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn norm(v: [f64; 3]) -> f64 {
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }

    #[test]
    fn periods() {
        let p = planets();
        assert!((p[2].period_days() - 365.25).abs() < 0.1);
        assert!((p[4].period_days() / 365.25 - 11.86).abs() < 0.02);
    }

    #[test]
    fn earth_at_j2000() {
        // Sun's apparent longitude on 2000-01-01 12:00 TT is ~280.4°, so Earth's is ~100.4°, at ~0.983 AU.
        let r = planets()[2].position_at(0.0);
        let lon = r[1].atan2(r[0]).to_degrees();
        assert!((lon - 100.4).abs() < 0.3, "lon = {lon}");
        assert!((r[0].hypot(r[1]) - 0.9833).abs() < 0.001);
    }

    #[test]
    fn earth_moves_one_orbit_per_year() {
        let earth = &planets()[2];
        let (r0, v0) = earth.state_at(0.0);
        let (r1, _) = earth.state_at(365.25636); // sidereal year
        let d = norm([r1[0] - r0[0], r1[1] - r0[1], r1[2] - r0[2]]);
        assert!(d < 1e5, "moved {d} km");
        // Perihelion speed is ~30.3 km/s, and Earth is near perihelion in early January.
        assert!((norm(v0) - 30.27).abs() < 0.05, "v = {}", norm(v0));
    }

    #[test]
    fn mercury_orbit_spans_perihelion_to_aphelion() {
        let p = &planets()[0];
        let r: Vec<f64> = p.orbit_xy(720).chunks(2).map(|c| c[0].hypot(c[1])).collect();
        let (min, max) = r.iter().fold((f64::MAX, 0.0f64), |(lo, hi), &v| (lo.min(v), hi.max(v)));
        // Plan-view radii are shrunk slightly by the 7° inclination.
        assert!((min - p.a_au * (1.0 - p.e)).abs() < 0.005, "min = {min}");
        assert!((max - p.a_au * (1.0 + p.e)).abs() < 0.005, "max = {max}");
    }
}
