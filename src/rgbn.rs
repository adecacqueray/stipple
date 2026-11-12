//! Recursive *adaptive* Gaussian blue noise.
//!
//! The energy is the one of GBN (Ahmed, Ren, Wonka 2022): every dot carries a Gaussian whose size follows
//! the local density (variance σ²/a, a = target density at the dot, mean 1), and the dots' kernels must
//! add up to the image. Its gradient is a pair repulsion with the combined kernel of the two dots, and an
//! attraction to the image blurred by the dot's own kernel (read from the pyramid in `field`).
//!
//! The speed is the one of rgbn: neighbours come from a grid sort (the sorted grid stretches with the
//! density, exactly like the kernels, so the same stencil serves dark and pale areas), and the dots are optimised
//! coarse to fine, each level cloned ×3 from the previous one.
use crate::field::{bilinear, Density, LEVELS, T0};
use crate::gridsort::GridSorter;
use crate::image::G;
use crate::rng::Rng;

const EXPANSION: f64 = 0.3;
const LR: f32 = 0.3;
/// Grid neighbourhood of the pair sum. A dot in a pale area next to a dark one has many more neighbours
/// within reach than a dot in a flat area: 7 keeps the density of pale areas within ~6 % (10: ~1 %, 3× slower).
const STENCIL: i32 = 7;
const ROOT_N_ITER: usize = 20;
const N_PER_STEP: usize = 10;
/// Kernels stop growing below this density (× the mean): keeps near-white areas finite.
const A_MIN: f32 = 1.0 / 32.0;
/// Longest move of a dot in one step, in local spacings.
const MAX_MOVE: f32 = 0.3;
const CUTOFF: f32 = 30.0;

pub trait Observer {
    fn progress(&mut self, _done: f64, _total: f64) {}
    fn level_done(&mut self, _pts: &[[f32; 2]]) {}
}
pub struct Silent;
impl Observer for Silent {}

fn half_ball(radius: i32) -> Vec<(i32, i32)> {
    let mut v = vec![];
    for a in -radius..=radius {
        for b in -radius..=radius {
            let d2 = a * a + b * b;
            if d2 > 0 && d2 <= radius * radius && (a > 0 || (a == 0 && b > 0)) { v.push((a, b)); }
        }
    }
    v
}

fn level_sizes(n: usize) -> Vec<usize> {
    let mut v = vec![n];
    let mut m = n;
    while m > 3000 { m = m / 3 + m % 3; v.push(m); }
    v.reverse();
    v
}

#[inline]
fn wrap(x: f32) -> f32 { x - x.floor() }

pub struct Rgbn {
    rng: Rng,
    density: Density,
    shifts: Vec<(i32, i32)>,
    /// (dot count, local density map a) of the last level optimised
    amap: Option<(usize, Vec<f32>)>,
}

impl Rgbn {
    pub fn new(rho: Vec<f64>, seed: u64) -> Self {
        Rgbn { rng: Rng::new(seed), density: Density::new(rho), shifts: half_ball(STENCIL), amap: None }
    }

    fn set_level(&mut self, n: usize) {
        if !matches!(&self.amap, Some((m, _)) if *m == n) {
            // density seen at the scale of the mean spacing of this level
            let mut a = self.density.smooth(0.25 / n as f64);
            for v in a.iter_mut() { *v = v.max(A_MIN); }
            self.amap = Some((n, a));
        }
    }

    #[inline]
    fn a_at(&self, p: [f32; 2]) -> f32 {
        let (i, w) = bilinear(p);
        let a = &self.amap.as_ref().unwrap().1;
        a[i[0]] * w[0] + a[i[1]] * w[1] + a[i[2]] * w[2] + a[i[3]] * w[3]
    }

    pub fn run(&mut self, n: usize, n_iter: usize, obs: &mut dyn Observer) -> Vec<[f32; 2]> {
        let sizes = level_sizes(n.max(2));
        let steps = |k: usize| if k == 0 { ROOT_N_ITER } else { n_iter } * N_PER_STEP;
        let total: f64 = sizes.iter().enumerate().map(|(k, &m)| (m * steps(k)) as f64).sum();
        let mut done = 0.0;
        let mut x: Vec<[f32; 2]> = vec![];
        for (k, &m) in sizes.iter().enumerate() {
            self.set_level(m);
            x = if k == 0 { self.draw(m) } else { self.clone(&x, m) };
            x = self.descend(x, steps(k), &mut || { done += (m * N_PER_STEP) as f64; obs.progress(done, total); });
            obs.level_done(&x);
        }
        x
    }

    pub fn refine(&mut self, x: Vec<[f32; 2]>, steps: usize, obs: &mut dyn Observer) -> Vec<[f32; 2]> {
        let n = x.len();
        if n == 0 { return x; }
        self.set_level(n);
        let mut done = 0;
        self.descend(x, steps, &mut || { done += N_PER_STEP; obs.progress(done.min(steps) as f64, steps as f64); })
    }

    /// Root level: dots drawn from the density (rejection sampling).
    fn draw(&mut self, m: usize) -> Vec<[f32; 2]> {
        let mut out = Vec::with_capacity(m);
        while out.len() < m {
            let p = [self.rng.f64() as f32, self.rng.f64() as f32];
            let c = ((p[0] * G as f32) as usize).min(G - 1) * G + ((p[1] * G as f32) as usize).min(G - 1);
            if self.density.rho[c] > self.rng.f64() as f32 * self.density.rho_max { out.push(p); }
        }
        out
    }

    /// Expand N/3 parents into N children on small randomly rotated triangles, sized by the local spacing.
    fn clone(&mut self, parents: &[[f32; 2]], n: usize) -> Vec<[f32; 2]> {
        let mut xv = parents.to_vec();
        self.rng.shuffle(&mut xv);
        let np = n / 3;
        let simplex = [(-0.75f64.sqrt(), -0.5), (0.75f64.sqrt(), -0.5), (0.0, 1.0)];
        let mut out: Vec<[f32; 2]> = xv[np..].to_vec();
        for p in &xv[..np] {
            let scale = EXPANSION / (n as f64 * self.a_at(*p) as f64).sqrt();
            let (s, c) = (self.rng.f64() * std::f64::consts::TAU).sin_cos();
            for &(a, b) in &simplex {
                out.push([wrap(p[0] + ((c * a - s * b) * scale) as f32), wrap(p[1] + ((s * a + c * b) * scale) as f32)]);
            }
        }
        out
    }

    /// Preconditioned gradient descent on the GBN energy. Dots are re-sorted onto the grid, and their
    /// kernel sizes re-read, every N_PER_STEP steps.
    fn descend(&mut self, x: Vec<[f32; 2]>, steps: usize, on_sort: &mut dyn FnMut()) -> Vec<[f32; 2]> {
        let n = x.len();
        let gi = (n as f64).sqrt().ceil() as usize;
        let cells = gi * gi;
        let sigma2 = 1.0 / n as f32;                 // pair kernel exp(-c r²/σ²), c = 2 a_i a_j / (a_i + a_j)
        let inv_s2 = 1.0 / sigma2;
        let pull = std::f32::consts::PI * sigma2 / 2.0;   // weight of the image term against the pair sum
        let sorter = GridSorter::new(gi);

        let mut grid: Vec<[f32; 2]> = vec![[f32::NAN; 2]; cells];
        grid[..n].copy_from_slice(&x);
        let mut ga = vec![1f32; cells];              // a per cell
        let mut lv = vec![(0usize, 0f32); cells];    // pyramid level and blend per cell
        let mut gx = vec![[0f32; 2]; cells];
        let mut tmp: Vec<[f32; 2]> = vec![[0.0; 2]; cells];
        let mut empty = vec![false; cells];
        let mut cb = vec![[0f32; 2]; gi];

        for step in 0..steps {
            if step % N_PER_STEP == 0 {
                for c in 0..cells {
                    let v = grid[c];
                    empty[c] = !(v[0].is_finite() && v[1].is_finite());
                    tmp[c] = if empty[c] { [self.rng.f64() as f32, self.rng.f64() as f32] } else { [wrap(v[0]), wrap(v[1])] };
                }
                let order = sorter.argsort(&tmp, &mut self.rng);
                for c in 0..cells {
                    let o = order[c] as usize;
                    grid[c] = if empty[o] { [f32::NAN; 2] } else { tmp[o] };
                    if empty[o] { continue; }
                    let a = self.a_at(grid[c]);
                    ga[c] = a;
                    // the dot's own kernel has variance σ²/(4a) per axis
                    let l = ((sigma2 as f64 / (4.0 * a as f64) / T0).log2() as f32).clamp(0.0, (LEVELS - 1) as f32 - 1e-3);
                    lv[c] = (l as usize, l.fract());
                }
                on_sort();
            }
            micro_grad(&grid, &ga, gi, &self.shifts, inv_s2, &mut gx, &mut cb);
            for c in 0..cells {
                let v = grid[c];
                if v[0].is_nan() { continue; }
                let a = ga[c];
                let (i, w) = bilinear(v);
                let (k, f) = lv[c];
                let (f0, f1) = (&self.density.grad[k], &self.density.grad[k + 1]);
                let mut at = [0f32; 2];
                for q in 0..4 {
                    for d in 0..2 { at[d] += w[q] * (f0[2 * i[q] + d] * (1.0 - f) + f1[2 * i[q] + d] * f); }
                }
                // descent direction, preconditioned by a⁻² so that every dot moves in units of its own spacing
                let s = LR / (a * a);
                let mut dx = s * (pull * at[0] - gx[c][0]);
                let mut dy = s * (pull * at[1] - gx[c][1]);
                let lim = MAX_MOVE * (sigma2 / a).sqrt();
                let m = (dx * dx + dy * dy).sqrt();
                if m > lim { dx *= lim / m; dy *= lim / m; }
                grid[c] = [wrap(v[0] + dx), wrap(v[1] + dy)];
            }
        }
        grid.into_iter().filter(|v| v[0].is_finite() && v[1].is_finite()).collect()
    }
}

#[inline(always)]
fn fast_exp(x: f32) -> f32 {
    let t = x * std::f32::consts::LOG2_E;
    let n = (t + 12_582_912.0) - 12_582_912.0;
    let f = (x - n * 0.693_145_75) - n * 1.428_606_8e-6;
    let p = 1.0 + f * (1.0 + f * (0.5 + f * (0.166_666_67 + f * (0.041_666_668 + f * (0.008_333_334 + f * 0.001_388_889)))));
    p * f32::from_bits(((n as i32 + 127) as u32) << 23)
}

#[inline(always)]
fn tdelta(d: f32) -> f32 {
    let d = if d > 0.5 { d - 1.0 } else { d };
    if d < -0.5 { d + 1.0 } else { d }
}

/// out = (x_q - x_p) · c² exp(-c r²/σ²): minus the pair force on p, plus the one on q.
#[inline(always)]
fn pair_segment(xp: &[[f32; 2]], xq: &[[f32; 2]], ap: &[f32], aq: &[f32], inv_s2: f32, out: &mut [[f32; 2]]) {
    for ((((a, b), &ai), &aj), o) in xp.iter().zip(xq).zip(ap).zip(aq).zip(out.iter_mut()) {
        let d0 = tdelta(b[0] - a[0]);
        let d1 = tdelta(b[1] - a[1]);
        let c = 2.0 * ai * aj / (ai + aj);
        let r = (d0 * d0 + d1 * d1) * inv_s2 * c;
        let w = c * c * fast_exp(-r.min(CUTOFF));
        *o = if r < CUTOFF { [d0 * w, d1 * w] } else { [0.0, 0.0] };   // NaN (empty cell) or far pair -> 0
    }
}

#[inline(never)]
fn micro_grad(grid: &[[f32; 2]], ga: &[f32], gi: usize, shifts: &[(i32, i32)], inv_s2: f32,
              gx: &mut [[f32; 2]], cb: &mut [[f32; 2]]) {
    for g in gx.iter_mut() { *g = [0.0; 2]; }
    for &(s0, s1) in shifts {
        let s = s1.rem_euclid(gi as i32) as usize;
        for i in 0..gi {
            let qi = (i + gi - s0 as usize) % gi;
            let (p, q) = (i * gi, qi * gi);
            pair_segment(&grid[p + s..p + gi], &grid[q..q + gi - s], &ga[p + s..p + gi], &ga[q..q + gi - s], inv_s2, &mut cb[s..]);
            pair_segment(&grid[p..p + s], &grid[q + gi - s..q + gi], &ga[p..p + s], &ga[q + gi - s..q + gi], inv_s2, &mut cb[..s]);
            for (g, c) in gx[p..p + gi].iter_mut().zip(cb.iter()) { g[0] += c[0]; g[1] += c[1]; }
            for (j, c) in cb.iter().enumerate() {
                let jq = if j >= s { j - s } else { j + gi - s };
                let g = &mut gx[q + jq]; g[0] -= c[0]; g[1] -= c[1];
            }
        }
    }
}
