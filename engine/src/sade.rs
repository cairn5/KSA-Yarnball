//! Self-adaptive differential evolution (SADE) with pygmo's defaults: mutation variant 2
//! (rand/1/exp), jDE adaptation (Brest et al. 2006), memory off, and a stop once the population
//! has collapsed (spread of fitness and of positions below 1e-6).
//!
//! It runs as a stepper, one generation at a time, so callers can report progress in between.

/// Small deterministic PRNG (SplitMix64), so runs repeat exactly on every platform.
#[derive(Clone)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1).
    pub fn f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in 0..n.
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    /// Standard normal (Box–Muller).
    pub fn normal(&mut self) -> f64 {
        let u = 1.0 - self.f64();
        let v = self.f64();
        (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()
    }
}

/// Collapse tolerances, as pygmo's sade defaults (ftol, xtol).
const TOL: f64 = 1e-6;

pub struct Sade {
    lb: Vec<f64>,
    ub: Vec<f64>,
    pop: Vec<Vec<f64>>,
    fit: Vec<f64>,
    f: Vec<f64>,
    cr: Vec<f64>,
    trial: Vec<f64>,
    pub rng: Rng,
    /// Generations run so far.
    pub generation: usize,
    /// Fitness evaluations so far.
    pub evaluations: usize,
}

impl Sade {
    /// Starts from `pop` (at least 4 members, each within `lb`..`ub`), evaluating each with `f`.
    pub fn new(pop: Vec<Vec<f64>>, lb: Vec<f64>, ub: Vec<f64>, rng: Rng, f: &mut impl FnMut(&[f64]) -> f64) -> Self {
        assert!(pop.len() >= 4, "SADE needs a population of at least 4");
        let fit: Vec<f64> = pop.iter().map(|x| f(x)).collect();
        let d = lb.len();
        let np = pop.len();
        let mut s = Sade {
            lb,
            ub,
            pop,
            fit,
            f: vec![0.0; np],
            cr: vec![0.0; np],
            trial: vec![0.0; d],
            rng,
            generation: 0,
            evaluations: np,
        };
        s.new_round();
        s
    }

    /// Starts a new round. With memory off, as in pygmo by default, every member's F and CR are
    /// redrawn: F uniformly in [0.1, 1], CR in [0, 1].
    pub fn new_round(&mut self) {
        for i in 0..self.pop.len() {
            self.f[i] = 0.1 + 0.9 * self.rng.f64();
            self.cr[i] = self.rng.f64();
        }
    }

    /// One generation: every member gets one trial, which replaces it if no worse.
    pub fn step(&mut self, f: &mut impl FnMut(&[f64]) -> f64) {
        let np = self.pop.len();
        let d = self.lb.len();
        for i in 0..np {
            // jDE: each member's parameters are redrawn with probability 0.1.
            let fi = if self.rng.f64() < 0.1 { 0.1 + 0.9 * self.rng.f64() } else { self.f[i] };
            let cri = if self.rng.f64() < 0.1 { self.rng.f64() } else { self.cr[i] };

            let pick = |not: &[usize], rng: &mut Rng| loop {
                let k = rng.below(np);
                if !not.contains(&k) {
                    return k;
                }
            };
            let a = pick(&[i], &mut self.rng);
            let b = pick(&[i, a], &mut self.rng);
            let c = pick(&[i, a, b], &mut self.rng);

            // rand/1 mutation with exponential crossover: a run of consecutive genes, at least one.
            self.trial.copy_from_slice(&self.pop[i]);
            let mut j = self.rng.below(d);
            let mut len = 0;
            loop {
                self.trial[j] = self.pop[a][j] + fi * (self.pop[b][j] - self.pop[c][j]);
                j = (j + 1) % d;
                len += 1;
                if !(self.rng.f64() < cri && len < d) {
                    break;
                }
            }
            // Genes outside the bounds are redrawn uniformly inside them.
            for j in 0..d {
                if !(self.trial[j] >= self.lb[j] && self.trial[j] <= self.ub[j]) {
                    self.trial[j] = self.lb[j] + self.rng.f64() * (self.ub[j] - self.lb[j]);
                }
            }

            let ft = f(&self.trial);
            self.evaluations += 1;
            if ft <= self.fit[i] {
                self.pop[i].copy_from_slice(&self.trial);
                self.fit[i] = ft;
                self.f[i] = fi;
                self.cr[i] = cri;
            }
        }
        self.generation += 1;
    }

    /// Index of the best member.
    fn best_index(&self) -> usize {
        (0..self.pop.len()).min_by(|&a, &b| self.fit[a].total_cmp(&self.fit[b])).unwrap()
    }

    /// The best member and its fitness.
    pub fn best(&self) -> (&[f64], f64) {
        let k = self.best_index();
        (&self.pop[k], self.fit[k])
    }

    /// Whether the population has collapsed: best and worst within the tolerances in both fitness
    /// and position, which is when pygmo's sade stops a round early.
    pub fn collapsed(&self) -> bool {
        let best = self.best_index();
        let worst = (0..self.pop.len()).max_by(|&a, &b| self.fit[a].total_cmp(&self.fit[b])).unwrap();
        let dx: f64 = self.pop[best].iter().zip(&self.pop[worst]).map(|(a, b)| (a - b).abs()).sum();
        let df = (self.fit[best] - self.fit[worst]).abs();
        dx < TOL && df < TOL
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_pop(lb: &[f64], ub: &[f64], np: usize, rng: &mut Rng) -> Vec<Vec<f64>> {
        (0..np).map(|_| lb.iter().zip(ub).map(|(l, u)| l + rng.f64() * (u - l)).collect()).collect()
    }

    #[test]
    fn minimises_rosenbrock() {
        let mut rosen = |x: &[f64]| (0..x.len() - 1).map(|i| 100.0 * (x[i + 1] - x[i] * x[i]).powi(2) + (1.0 - x[i]).powi(2)).sum();
        let (lb, ub) = (vec![-5.0; 4], vec![5.0; 4]);
        let mut rng = Rng::new(7);
        let pop = random_pop(&lb, &ub, 40, &mut rng);
        let mut s = Sade::new(pop, lb, ub, rng, &mut rosen);
        for round in 0..6 {
            if round > 0 {
                s.new_round();
            }
            for _ in 0..500 {
                s.step(&mut rosen);
                if s.collapsed() {
                    break;
                }
            }
        }
        assert!(s.best().1 < 1e-6, "f = {}", s.best().1);
    }

    #[test]
    fn never_loses_the_best() {
        let mut sphere = |x: &[f64]| x.iter().map(|v| (v - 3.0).powi(2)).sum();
        let (lb, ub) = (vec![-10.0; 6], vec![10.0; 6]);
        let mut rng = Rng::new(1);
        let mut pop = random_pop(&lb, &ub, 20, &mut rng);
        pop[0] = vec![3.0; 6];
        let mut s = Sade::new(pop, lb, ub, rng, &mut sphere);
        for _ in 0..50 {
            s.step(&mut sphere);
            assert_eq!(s.best().1, 0.0);
        }
    }
}
