//! Photo colour per dot: under the dot, or averaged over its Voronoi cell.

/// Mean colour of the 3×3 pixels under each dot.
pub fn under_dots(rgba: &[u8], s: usize, xy: &[f64]) -> Vec<u8> {
    let n = xy.len() / 2;
    let mut out = vec![0u8; 3 * n];
    for i in 0..n {
        let (cx, cy) = (xy[2 * i].round() as i64, xy[2 * i + 1].round() as i64);
        let (mut acc, mut k) = ([0u32; 3], 0u32);
        for dy in -1..=1 {
            for dx in -1..=1 {
                let (x, y) = (cx + dx, cy + dy);
                if x < 0 || y < 0 || x >= s as i64 || y >= s as i64 { continue; }
                let o = (y as usize * s + x as usize) * 4;
                for c in 0..3 { acc[c] += rgba[o + c] as u32; }
                k += 1;
            }
        }
        if k > 0 { for c in 0..3 { out[3 * i + c] = (acc[c] / k) as u8; } }
    }
    out
}

/// Mean colour over each Voronoi cell, sampled at the dot (counted twice) and halfway to every vertex.
pub fn over_cells(rgba: &[u8], s: usize, xy: &[f64], verts: &[f32], off: &[u32]) -> Vec<u8> {
    let n = xy.len() / 2;
    let mut out = vec![0u8; 3 * n];
    let at = |x: f64, y: f64, acc: &mut [f64; 4]| {
        let cl = |v: f64| (v as i64).clamp(0, s as i64 - 1) as usize;
        let o = (cl(y) * s + cl(x)) * 4;
        for c in 0..3 { acc[c] += rgba[o + c] as f64; }
        acc[3] += 1.0;
    };
    for i in 0..n {
        let (x, y) = (xy[2 * i], xy[2 * i + 1]);
        let mut acc = [0f64; 4];
        at(x, y, &mut acc);
        at(x, y, &mut acc);
        for v in off[i] as usize..off[i + 1] as usize {
            at((x + verts[2 * v] as f64) / 2.0, (y + verts[2 * v + 1] as f64) / 2.0, &mut acc);
        }
        for c in 0..3 { out[3 * i + c] = (acc[c] / acc[3]) as u8; }
    }
    out
}

/// Ink strength of each cell from its area. Dots were placed with density ∝ darkness^1.5, so
/// darkness ∝ (1/area)^(2/3) recovers the photo's tone; the 2nd-percentile area is full ink.
pub fn ink_from_area(area: &[f64]) -> Vec<f32> {
    let mut sorted: Vec<f64> = area.iter().copied().filter(|a| a.is_finite()).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let a_ref = sorted.get((sorted.len() as f64 * 0.02) as usize).copied().filter(|&a| a != 0.0).unwrap_or(1.0);
    area.iter().map(|&a| (a_ref / a).powf(2.0 / 3.0).min(1.0) as f32).collect()
}
