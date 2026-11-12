//! Target side of adaptive Gaussian blue noise: the image seen through Gaussians of every size.
//! `grad[k]` is the gradient of (rho ⊛ g_t) for t = T0·2^k (g_t: unit-mass Gaussian of variance t per axis),
//! so a dot with its own kernel size reads its attraction between two levels of this pyramid.
use crate::image::G;
use rustfft::{num_complex::Complex, Fft, FftPlanner};
use std::sync::Arc;

type C = Complex<f64>;

pub struct Fft2 { fwd: Arc<dyn Fft<f64>>, inv: Arc<dyn Fft<f64>>, scratch: Vec<C>, t: Vec<C> }

impl Fft2 {
    pub fn new() -> Self {
        let mut p = FftPlanner::new();
        let fwd = p.plan_fft_forward(G);
        let inv = p.plan_fft_inverse(G);
        let n = fwd.get_inplace_scratch_len().max(inv.get_inplace_scratch_len());
        Fft2 { fwd, inv, scratch: vec![C::default(); n], t: vec![C::default(); G * G] }
    }
    fn transpose(&mut self, a: &mut [C]) {
        for i in 0..G { for j in 0..G { self.t[j * G + i] = a[i * G + j]; } }
        a.copy_from_slice(&self.t);
    }
    pub fn run(&mut self, a: &mut [C], inverse: bool) {
        let f = if inverse { self.inv.clone() } else { self.fwd.clone() };
        for _ in 0..2 {
            for row in a.chunks_exact_mut(G) { f.process_with_scratch(row, &mut self.scratch); }
            self.transpose(a);
        }
        if inverse { let s = 1.0 / (G * G) as f64; for v in a.iter_mut() { *v *= s; } }
    }
}

#[inline]
fn freq(k: usize) -> f64 { if k < G / 2 { k as f64 } else { k as f64 - G as f64 } }

/// Smallest kernel variance of the pyramid: half a pixel of standard deviation.
pub const T0: f64 = 0.25 / (G * G) as f64;
/// Pyramid levels (one octave of variance each): up to a standard deviation of ~45 pixels.
pub const LEVELS: usize = 14;

pub struct Density {
    fft: Fft2,
    /// FFT of the density, scaled to mean 1.
    rho_hat: Vec<C>,
    pub rho_max: f32,
    /// rho (mean 1) on the grid, for the initial draw.
    pub rho: Vec<f32>,
    /// grad[k][2·(ix·G + iy)..]: gradient of rho ⊛ g_(T0·2^k), interleaved x, y.
    pub grad: Vec<Vec<f32>>,
}

impl Density {
    pub fn new(rho: Vec<f64>) -> Self {
        let mut fft = Fft2::new();
        let mean = rho.iter().sum::<f64>() / (G * G) as f64;
        let rho: Vec<f32> = rho.iter().map(|&v| (v / mean) as f32).collect();
        let mut rho_hat: Vec<C> = rho.iter().map(|&v| C::new(v as f64, 0.0)).collect();
        fft.run(&mut rho_hat, false);
        let pi = std::f64::consts::PI;
        let mut grad = Vec::with_capacity(LEVELS);
        let mut buf = vec![C::default(); G * G];
        for k in 0..LEVELS {
            let t = T0 * (1u64 << k) as f64;
            // both components in one complex inverse FFT: real part d/dx, imaginary part d/dy
            for (i, b) in buf.iter_mut().enumerate() {
                let (a, c) = (freq(i / G), freq(i % G));
                let g = (-2.0 * pi * pi * t * (a * a + c * c)).exp();
                *b = rho_hat[i] * C::new(0.0, 2.0 * pi * g) * C::new(a, c);
            }
            fft.run(&mut buf, true);
            grad.push(buf.iter().flat_map(|z| [z.re as f32, z.im as f32]).collect());
        }
        let rho_max = rho.iter().cloned().fold(0.0, f32::max);
        Density { fft, rho_hat, rho_max, rho, grad }
    }

    /// rho ⊛ g_t on the grid (mean 1): the local density that sizes the kernels at one level.
    pub fn smooth(&mut self, t: f64) -> Vec<f32> {
        let pi = std::f64::consts::PI;
        let mut buf: Vec<C> = self.rho_hat.iter().enumerate().map(|(i, r)| {
            let (a, c) = (freq(i / G), freq(i % G));
            r * (-2.0 * pi * pi * t * (a * a + c * c)).exp()
        }).collect();
        self.fft.run(&mut buf, true);
        buf.iter().map(|z| z.re as f32).collect()
    }
}

/// Bilinear weights of a point of the unit torus on the grid (samples at cell centres).
#[inline]
pub fn bilinear(p: [f32; 2]) -> ([usize; 4], [f32; 4]) {
    let (u, v) = (p[0] * G as f32 - 0.5, p[1] * G as f32 - 0.5);
    let (fu, fv) = (u.floor(), v.floor());
    let (a, b) = (u - fu, v - fv);
    let m = G as i32;
    let (x0, y0) = ((fu as i32).rem_euclid(m) as usize, (fv as i32).rem_euclid(m) as usize);
    let (x1, y1) = ((x0 + 1) % G, (y0 + 1) % G);
    ([x0 * G + y0, x0 * G + y1, x1 * G + y0, x1 * G + y1], [(1.0 - a) * (1.0 - b), (1.0 - a) * b, a * (1.0 - b), a * b])
}
