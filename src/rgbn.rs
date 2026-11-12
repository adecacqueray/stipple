//! Port of blue_sampler.run.run_recursive._recursive_pipeline for D=2 with an image target
//! (what im2points -> sample_points(method="rgbn", targets=image) executes), plus `refine`, which
//! continues the finest level for more steps ("Optimise further" in the page).
use crate::field::{Density, Fft2};
use crate::gridsort::GridSorter;
use crate::image::G;
use crate::rng::Rng;

// 2D preset with target: S=0.5, expension_factor=0.3, LR_spatial=0.1, spatial_radius=7
const S: f64 = 0.5;
const EXPANSION: f64 = 0.3;
const LR_SPATIAL: f64 = 0.1;
const SPATIAL_RADIUS: i32 = 7;
const ROOT_N_ITER: usize = 20;
const N_PER_STEP: usize = 10;
/// exp(-d2/sigma2) is skipped when d2/sigma2 > CUTOFF (contribution < 1e-15, below f32 precision)
const CUTOFF: f32 = 34.0;

pub trait Observer {
    fn progress(&mut self, _done: f64, _total: f64) {}
    fn level_done(&mut self, _pts: &[[f32; 2]]) {}
}
pub struct Silent;
impl Observer for Silent {}

/// integers_in_half_ball(radius, 2)
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

/// The rgbn optimiser for one target density. Keeps its random stream and the target field of the
/// last level, so `refine` simply carries on where `run` stopped.
pub struct Rgbn {
    rng: Rng,
    fft: Fft2,
    density: Density,
    shifts: Vec<(i32, i32)>,
    /// (dot count, target field) of the last level optimised
    field: Option<(usize, Vec<f32>)>,
}

impl Rgbn {
    pub fn new(rho: Vec<f64>, seed: u64) -> Self {
        let mut fft = Fft2::new();
        let density = Density::new(rho, &mut fft);
        Rgbn { rng: Rng::new(seed), fft, density, shifts: half_ball(SPATIAL_RADIUS), field: None }
    }

    /// Full recursive pipeline: n dots, `n_iter` passes (of N_PER_STEP steps) per level.
    pub fn run(&mut self, n: usize, n_iter: usize, obs: &mut dyn Observer) -> Vec<[f32; 2]> {
        let sizes = level_sizes(n.max(2));
        let steps = |k: usize| if k == 0 { ROOT_N_ITER } else { n_iter } * N_PER_STEP;
        let total: f64 = sizes.iter().enumerate().map(|(k, &m)| (m * steps(k)) as f64).sum();
        let mut done = 0.0;
        let mut x: Vec<[f32; 2]> = vec![];
        for (k, &m) in sizes.iter().enumerate() {
            // initial points: random (root) or cloned from the coarser level
            x = if k == 0 {
                (0..m).map(|_| [self.rng.f64() as f32, self.rng.f64() as f32]).collect()
            } else {
                clone(&x, m, &mut self.rng)
            };
            let field = self.density.field(S * 2.0 / m as f64, &mut self.fft, &mut self.rng);
            self.field = Some((m, field));
            x = self.descend(x, steps(k), &mut || { done += (m * N_PER_STEP) as f64; obs.progress(done, total); });
            obs.level_done(&x);
        }
        x
    }

    /// `steps` more rgbn steps on existing dots, at the finest level (same field as the end of `run`).
    pub fn refine(&mut self, x: Vec<[f32; 2]>, steps: usize, obs: &mut dyn Observer) -> Vec<[f32; 2]> {
        let n = x.len();
        if n == 0 { return x; }
        if !matches!(&self.field, Some((m, _)) if *m == n) {
            let field = self.density.field(S * 2.0 / n as f64, &mut self.fft, &mut self.rng);
            self.field = Some((n, field));
        }
        let mut done = 0;
        self.descend(x, steps, &mut || { done += N_PER_STEP; obs.progress(done.min(steps) as f64, steps as f64); })
    }

    /// Gradient descent of the dots on the current field. The dots are re-sorted onto the grid
    /// every N_PER_STEP steps (`on_sort` is called then).
    fn descend(&mut self, x: Vec<[f32; 2]>, steps: usize, on_sort: &mut dyn FnMut()) -> Vec<[f32; 2]> {
        let n = x.len();
        let field = &self.field.as_ref().unwrap().1;
        let gi = (n as f64).sqrt().ceil() as usize;
        let gi = if gi * gi < n { gi + 1 } else { gi };
        let cells = gi * gi;
        let inv_s2 = (1.0 / (S * 2.0 / n as f64)) as f32;
        let lr = (LR_SPATIAL / S * 0.5) as f32; // lr_micro = lr_macro
        let sorter = GridSorter::new(gi);

        // prepare_points: first n cells real, others NaN (fictive)
        let mut grid: Vec<[f32; 2]> = vec![[f32::NAN; 2]; cells];
        grid[..n].copy_from_slice(&x);
        let mut gx = vec![[0f32; 2]; cells];
        let mut tmp: Vec<[f32; 2]> = vec![[0.0; 2]; cells];
        let mut empty = vec![false; cells];
        let mut cb = vec![[0f32; 2]; gi];

        for step in 0..steps {
            if step % N_PER_STEP == 0 {
                // gridify: empties get random positions for the sort, then back to NaN
                for c in 0..cells {
                    let v = grid[c];
                    empty[c] = !(v[0].is_finite() && v[1].is_finite());
                    tmp[c] = if empty[c] { [self.rng.f64() as f32, self.rng.f64() as f32] } else { [wrap(v[0]), wrap(v[1])] };
                }
                let order = sorter.argsort(&tmp, &mut self.rng);
                for c in 0..cells {
                    let o = order[c] as usize;
                    grid[c] = if empty[o] { [f32::NAN; 2] } else { tmp[o] };
                }
                on_sort();
            }
            // micro gradient: pair forces over the half-ball of grid shifts (action/reaction)
            micro_grad(&grid, gi, &self.shifts, inv_s2, &mut gx, &mut cb);
            // descent step with the target field (nearest cell lookup on the 512x512 grid)
            for c in 0..cells {
                let v = grid[c];
                if v[0].is_nan() { continue; }
                let ix = ((v[0] * G as f32).floor() as i32).rem_euclid(G as i32) as usize;
                let iy = ((v[1] * G as f32).floor() as i32).rem_euclid(G as i32) as usize;
                let f = 2 * (ix * G + iy);
                grid[c] = [
                    wrap(v[0] - lr * (gx[c][0] + field[f])),
                    wrap(v[1] - lr * (gx[c][1] + field[f + 1])),
                ];
            }
        }
        grid.into_iter().filter(|v| v[0].is_finite() && v[1].is_finite()).collect()
    }
}

/// Expand N//3 parents into N children with randomly rotated equilateral-triangle offsets.
fn clone(parents: &[[f32; 2]], n: usize, rng: &mut Rng) -> Vec<[f32; 2]> {
    let mut xv = parents.to_vec();
    rng.shuffle(&mut xv);
    let np = n / 3;
    let scale = EXPANSION / (n as f64).sqrt();
    // simplex(2): (-sqrt3/2, -1/2), (sqrt3/2, -1/2), (0, 1)
    let simplex = [(-0.75f64.sqrt(), -0.5), (0.75f64.sqrt(), -0.5), (0.0, 1.0)];
    let mut out: Vec<[f32; 2]> = xv[np..].to_vec();
    for p in &xv[..np] {
        let t = rng.f64() * std::f64::consts::TAU;
        let (s, c) = t.sin_cos();
        for &(a, b) in &simplex {
            let ox = (c * a - s * b) * scale;
            let oy = (s * a + c * b) * scale;
            let px = p[0] as f64 + ox; let py = p[1] as f64 + oy;
            out.push([(px - px.floor()) as f32, (py - py.floor()) as f32]);
        }
    }
    out
}

/// exp(x) for x in [-CUTOFF, 0], branch-free (auto-vectorises), rel. error ~1e-7.
#[inline(always)]
fn fast_exp(x: f32) -> f32 {
    let t = x * std::f32::consts::LOG2_E;
    let n = (t + 12_582_912.0) - 12_582_912.0; // round-to-nearest without SSE4.1
    let f = (x - n * 0.693_145_75) - n * 1.428_606_8e-6;
    let p = 1.0 + f * (1.0 + f * (0.5 + f * (0.166_666_67 + f * (0.041_666_668 + f * (0.008_333_334 + f * 0.001_388_889)))));
    p * f32::from_bits(((n as i32 + 127) as u32) << 23)
}

/// torus_delta for d in (-1, 1): d - round(d), as compare/select (vectorises everywhere)
#[inline(always)]
fn tdelta(d: f32) -> f32 {
    let d = if d > 0.5 { d - 1.0 } else { d };
    if d < -0.5 { d + 1.0 } else { d }
}

#[inline(always)]
fn pair_segment(xp: &[[f32; 2]], xq: &[[f32; 2]], inv_s2: f32, out: &mut [[f32; 2]]) {
    for ((a, b), o) in xp.iter().zip(xq).zip(out.iter_mut()) {
        let d0 = tdelta(b[0] - a[0]);
        let d1 = tdelta(b[1] - a[1]);
        let r = (d0 * d0 + d1 * d1) * inv_s2;
        let w = fast_exp(-r.min(CUTOFF));
        // NaN (fictive cell) or far pair -> 0, like clean_grad
        *o = if r < CUTOFF { [d0 * w, d1 * w] } else { [0.0, 0.0] };
    }
}

#[inline(never)]
fn micro_grad(grid: &[[f32; 2]], gi: usize, shifts: &[(i32, i32)], inv_s2: f32,
              gx: &mut [[f32; 2]], cb: &mut [[f32; 2]]) {
    for g in gx.iter_mut() { *g = [0.0; 2]; }
    for &(s0, s1) in shifts {
        let s = s1.rem_euclid(gi as i32) as usize;
        for i in 0..gi {
            let qi = (i + gi - s0 as usize) % gi;
            let (p, q) = (i * gi, qi * gi);
            // j in s..gi -> jq = j - s ; j in 0..s -> jq = j + gi - s
            pair_segment(&grid[p + s..p + gi], &grid[q..q + gi - s], inv_s2, &mut cb[s..]);
            pair_segment(&grid[p..p + s], &grid[q + gi - s..q + gi], inv_s2, &mut cb[..s]);
            for (g, c) in gx[p..p + gi].iter_mut().zip(cb.iter()) { g[0] += c[0]; g[1] += c[1]; }
            for (j, c) in cb.iter().enumerate() {
                let jq = if j >= s { j - s } else { j + gi - s };
                let g = &mut gx[q + jq]; g[0] -= c[0]; g[1] -= c[1];
            }
        }
    }
}
