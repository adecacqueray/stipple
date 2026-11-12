//! Port of blue_sampler.grad.im2fields.im2field:
//! PIL .convert("L") + .resize((512,512), LANCZOS) + transpose/flip + invert + normalise.
pub const G: usize = 512;

/// PIL ImagingConvert rgb2l (ITU-R 601-2 luma, fixed point).
pub fn luma(r: u8, g: u8, b: u8) -> u8 {
    ((r as u32 * 19595 + g as u32 * 38470 + b as u32 * 7471 + 0x8000) >> 16) as u8
}

fn lanczos(x: f64) -> f64 {
    fn sinc(x: f64) -> f64 {
        if x == 0.0 { 1.0 } else { let a = x * std::f64::consts::PI; a.sin() / a }
    }
    if -3.0 <= x && x < 3.0 { sinc(x) * sinc(x / 3.0) } else { 0.0 }
}

const PRECISION_BITS: u32 = 32 - 8 - 2;

/// PIL precompute_coeffs + normalize_coeffs_8bpc. Returns (xmin, xsize) bounds and int coefficients.
fn coeffs(in_size: usize, out_size: usize) -> (Vec<(usize, usize)>, Vec<i32>, usize) {
    let scale = in_size as f64 / out_size as f64;
    let filterscale = scale.max(1.0);
    let support = 3.0 * filterscale;
    let ksize = support.ceil() as usize * 2 + 1;
    let mut bounds = Vec::with_capacity(out_size);
    let mut kk = vec![0i32; out_size * ksize];
    let ss = 1.0 / filterscale;
    let mut k = vec![0f64; ksize];
    for xx in 0..out_size {
        let center = (xx as f64 + 0.5) * scale;
        let xmin = ((center - support + 0.5) as i64).max(0) as usize;
        let xmax = ((center + support + 0.5) as i64).min(in_size as i64) as usize - xmin;
        let mut ww = 0.0;
        for x in 0..xmax {
            let w = lanczos((x as f64 + xmin as f64 - center + 0.5) * ss);
            k[x] = w;
            ww += w;
        }
        for x in 0..xmax {
            let v = if ww != 0.0 { k[x] / ww } else { k[x] };
            let f = v * (1u64 << PRECISION_BITS) as f64;
            kk[xx * ksize + x] = if v < 0.0 { (-0.5 + f) as i32 } else { (0.5 + f) as i32 };
        }
        bounds.push((xmin, xmax));
    }
    (bounds, kk, ksize)
}

#[inline]
fn clip8(ss: i64) -> u8 { (ss >> PRECISION_BITS).clamp(0, 255) as u8 }

/// Lanczos resize of an 8-bit L image (horizontal pass then vertical, like PIL).
fn resize_l(src: &[u8], w: usize, h: usize, ow: usize, oh: usize) -> Vec<u8> {
    let (hb, hk, hks) = coeffs(w, ow);
    let mut tmp = vec![0u8; ow * h];
    for y in 0..h {
        let row = &src[y * w..(y + 1) * w];
        for xx in 0..ow {
            let (xmin, xsz) = hb[xx];
            let mut ss: i64 = 1 << (PRECISION_BITS - 1);
            for x in 0..xsz { ss += row[xmin + x] as i64 * hk[xx * hks + x] as i64; }
            tmp[y * ow + xx] = clip8(ss);
        }
    }
    let (vb, vk, vks) = coeffs(h, oh);
    let mut out = vec![0u8; ow * oh];
    for yy in 0..oh {
        let (ymin, ysz) = vb[yy];
        for xx in 0..ow {
            let mut ss: i64 = 1 << (PRECISION_BITS - 1);
            for y in 0..ysz { ss += tmp[(ymin + y) * ow + xx] as i64 * vk[yy * vks + y] as i64; }
            out[yy * ow + xx] = clip8(ss);
        }
    }
    out
}

/// RGBA (w*h*4) -> density rho[ix*512 + iy], ix = image column, iy = upward row. Sums to 1.
pub fn im2field(rgba: &[u8], w: usize, h: usize, invert: bool) -> Vec<f64> {
    let l: Vec<u8> = rgba.chunks_exact(4).map(|p| luma(p[0], p[1], p[2])).collect();
    let small = resize_l(&l, w, h, G, G);
    let mut rho = vec![0f32; G * G];
    for ix in 0..G {
        for iy in 0..G {
            let v = small[(G - 1 - iy) * G + ix] as f32 / 255.0;
            rho[ix * G + iy] = if invert { 1.0 - v } else { v };
        }
    }
    let s: f32 = rho.iter().sum();
    rho.iter().map(|&v| (v / s) as f64).collect()
}
