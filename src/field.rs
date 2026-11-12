//! Port of blue_sampler.grad.fields (_field_from_density, _log_barier, spatial_grad lookup).
use crate::image::G;
use crate::rng::Rng;
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

pub struct Density { rho_max: f64, rho: Vec<f64>, rho_hat: Vec<C> }

impl Density {
    pub fn new(rho: Vec<f64>, fft: &mut Fft2) -> Self {
        let mut rho_hat: Vec<C> = rho.iter().map(|&v| C::new(v, 0.0)).collect();
        fft.run(&mut rho_hat, false);
        let rho_max = rho.iter().cloned().fold(f64::MIN, f64::max);
        Density { rho_max, rho, rho_hat }
    }

    /// Gradient field (G*G*2, f32, interleaved) for kernel width sigma2.
    /// field = grad( V * (rho.max - rho - log_barrier(V * rho)) ), V = exp(-r^2/sigma2).
    pub fn field(&self, sigma2: f64, fft: &mut Fft2, rng: &mut Rng) -> Vec<f32> {
        let pi = std::f64::consts::PI;
        let c = (G * G) as f64 * (pi * sigma2);
        let v_hat: Vec<f64> = (0..G * G)
            .map(|k| { let (a, b) = (freq(k / G), freq(k % G)); c * (-(pi * pi) * sigma2 * (a * a + b * b)).exp() })
            .collect();
        // rho_smooth
        let mut buf: Vec<C> = self.rho_hat.iter().zip(&v_hat).map(|(r, v)| r * v).collect();
        fft.run(&mut buf, true);
        let mean = buf.iter().map(|z| z.re).sum::<f64>() / (G * G) as f64;
        let alpha = mean * 1e-3;
        // anti_rho - barrier
        for (k, z) in buf.iter_mut().enumerate() {
            let barrier = alpha / ((z.re / mean).max(0.0) + 1e-4) * (1.0 + rng.f64());
            *z = C::new(self.rho_max - self.rho[k] - barrier, 0.0);
        }
        fft.run(&mut buf, false);
        // both gradient components in one complex inverse FFT (each is real up to the
        // Nyquist bin, where V_hat underflows to ~0 for every sigma2 used here)
        for k in 0..G * G {
            let (a, b) = (freq(k / G), freq(k % G));
            let g = C::new(0.0, 2.0 * pi * v_hat[k]);
            buf[k] = buf[k] * g * C::new(a, b);
        }
        fft.run(&mut buf, true);
        let mut out = vec![0f32; G * G * 2];
        for k in 0..G * G { out[2 * k] = buf[k].re as f32; out[2 * k + 1] = buf[k].im as f32; }
        out
    }
}
