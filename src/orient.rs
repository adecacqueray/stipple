//! Stroke directions: perpendicular to the photo's gradient, from a smoothed structure tensor.
//! Where the photo is flat the tensor is biased toward 45° hatching.
use crate::blur::gaussian3;
use std::f64::consts::PI;

pub struct Orient { pub m: usize, pub cx: Vec<f32>, pub cy: Vec<f32> }

/// `rgba`: the S×S padded square.
pub fn orient_field(rgba: &[u8], s: usize) -> Orient {
    let m = 512;
    let mut lum = vec![0f32; m * m];
    for y in 0..m {
        let sy = ((((y as f64 + 0.5) * s as f64 / m as f64) as usize).min(s - 1)) * s;
        for x in 0..m {
            let sx = (((x as f64 + 0.5) * s as f64 / m as f64) as usize).min(s - 1);
            let o = (sy + sx) * 4;
            lum[y * m + x] = ((0.2126 * rgba[o] as f64 + 0.7152 * rgba[o + 1] as f64 + 0.0722 * rgba[o + 2] as f64) / 255.0) as f32;
        }
    }
    let l = gaussian3(&lum, m, 1.2);
    let (mut jxx, mut jxy, mut jyy) = (vec![0f32; m * m], vec![0f32; m * m], vec![0f32; m * m]);
    let mut tr = 0f64;
    for y in 0..m {
        for x in 0..m {
            let i = y * m + x;
            let gx = (l[y * m + (x + 1).min(m - 1)] as f64 - l[y * m + x.saturating_sub(1)] as f64) / 2.0;
            let gy = (l[(y + 1).min(m - 1) * m + x] as f64 - l[y.saturating_sub(1) * m + x] as f64) / 2.0;
            jxx[i] = (gx * gx) as f32;
            jxy[i] = (gx * gy) as f32;
            jyy[i] = (gy * gy) as f32;
            tr += jxx[i] as f64 + jyy[i] as f64;
        }
    }
    let (jxx, jxy, jyy) = (gaussian3(&jxx, m, 5.0), gaussian3(&jxy, m, 5.0), gaussian3(&jyy, m, 5.0));
    let bias = 0.04 * tr / (m * m) as f64;     // normal of a 45° hatch is (1,1)/√2 → adds bias/2 to each term
    let (mut cx, mut cy) = (vec![0f32; m * m], vec![0f32; m * m]);
    for i in 0..m * m {
        let (a, b, c) = (jxx[i] as f64 + bias / 2.0, jxy[i] as f64 + bias / 2.0, jyy[i] as f64 + bias / 2.0);
        let phi = 0.5 * (2.0 * b).atan2(a - c) + PI / 2.0;   // gradient direction + 90°
        cx[i] = phi.cos() as f32;
        cy[i] = phi.sin() as f32;
    }
    Orient { m, cx, cy }
}

impl Orient {
    /// Unit direction at each dot (`xy` in square pixels of a square of side `s`).
    pub fn at(&self, xy: &[f64], s: f64) -> Vec<f32> {
        let n = xy.len() / 2;
        let k = self.m as f64 / s;
        let cl = |v: f64| ((v * k) as i64).clamp(0, self.m as i64 - 1) as usize;
        let mut out = vec![0f32; 2 * n];
        for i in 0..n {
            let j = cl(xy[2 * i + 1]) * self.m + cl(xy[2 * i]);
            out[2 * i] = self.cx[j];
            out[2 * i + 1] = self.cy[j];
        }
        out
    }
}
