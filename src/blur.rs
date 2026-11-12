//! Separable Gaussian approximation by three box passes (Kovesi), on an M×M grid with clamped edges.

/// Blur `src` (row-major M×M) with a Gaussian of standard deviation `s` grid pixels.
pub fn gaussian3(src: &[f32], m: usize, s: f64) -> Vec<f32> {
    let nb = 3.0;
    let w_ideal = (12.0 * s * s / nb + 1.0).sqrt();
    let mut wl = w_ideal.floor() as i64;
    if wl % 2 == 0 { wl -= 1; }
    let wlf = wl as f64;
    let mm = ((12.0 * s * s - nb * wlf * wlf - 4.0 * nb * wlf - 3.0 * nb) / (-4.0 * wlf - 4.0)).round() as i64;
    let mut a = src.to_vec();
    let mut b = vec![0f32; m * m];
    for p in 0..3 {
        let r = ((if p < mm { wl } else { wl + 2 }) - 1) / 2;
        box_h(&a, &mut b, m, r);
        box_v(&b, &mut a, m, r);
    }
    a
}

#[inline]
fn cl(i: i64, m: usize) -> usize { i.clamp(0, m as i64 - 1) as usize }

fn box_h(src: &[f32], dst: &mut [f32], m: usize, r: i64) {
    let inv = 1.0 / (2 * r + 1) as f64;
    for y in 0..m {
        let row = &src[y * m..(y + 1) * m];
        let mut acc = 0f64;
        for i in -r - 1..r { acc += row[cl(i, m)] as f64; }
        for x in 0..m as i64 {
            acc += row[cl(x + r, m)] as f64 - row[cl(x - r - 1, m)] as f64;
            dst[y * m + x as usize] = (acc * inv) as f32;
        }
    }
}

fn box_v(src: &[f32], dst: &mut [f32], m: usize, r: i64) {
    let inv = 1.0 / (2 * r + 1) as f64;
    for x in 0..m {
        let mut acc = 0f64;
        for i in -r - 1..r { acc += src[cl(i, m) * m + x] as f64; }
        for y in 0..m as i64 {
            acc += src[cl(y + r, m) * m + x] as f64 - src[cl(y - r - 1, m) * m + x] as f64;
            dst[y as usize * m + x] = (acc * inv) as f32;
        }
    }
}
