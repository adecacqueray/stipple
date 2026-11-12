//! native benchmark: bench <rgba.raw> <w> <h> <N> <seed> <out.f32> [reps]
use std::time::Instant;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let rgba = std::fs::read(&a[1]).unwrap();
    let (w, h, n, seed): (usize, usize, usize, u64) =
        (a[2].parse().unwrap(), a[3].parse().unwrap(), a[4].parse().unwrap(), a[5].parse().unwrap());
    let reps: usize = a.get(7).map(|s| s.parse().unwrap()).unwrap_or(3);
    let mut pts = vec![];
    let mut ts = vec![];
    for r in 0..reps {
        let t = Instant::now();
        let rho = stipple::image::im2field(&rgba, w, h, true);
        pts = stipple::rgbn::Rgbn::new(rho, seed + r as u64).run(n, 12, &mut stipple::rgbn::Silent);
        ts.push(t.elapsed().as_secs_f64());
    }
    let bytes: Vec<u8> = pts.iter().flat_map(|p| p.iter().flat_map(|v| v.to_le_bytes())).collect();
    std::fs::write(&a[6], bytes).unwrap();
    println!("{n} rust times {:?} npts {}", ts.iter().map(|t| format!("{t:.3}s")).collect::<Vec<_>>(), pts.len());
}
