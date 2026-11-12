//! Port of gridpoints.argsort(X, (I, I), level=2, init="kdtree", n_iter=40) for 2D.
//! Returns order such that X[order].reshape(I, I, 2) is monotonic along rows, columns
//! and both diagonals.
use crate::rng::Rng;

pub struct GridSorter {
    i: usize,
    /// fibers per direction: flat cell lists + offsets. dirs: (0,1),(1,0),(1,-1),(1,1)
    fibers: [(Vec<u32>, Vec<u32>); 4],
}

impl GridSorter {
    pub fn new(i: usize) -> Self {
        let mut f: [(Vec<u32>, Vec<u32>); 4] = Default::default();
        let push = |f: &mut (Vec<u32>, Vec<u32>), cells: &mut dyn Iterator<Item = (usize, usize)>| {
            f.1.push(f.0.len() as u32);
            for (r, c) in cells { f.0.push((r * i + c) as u32); }
        };
        for r in 0..i { push(&mut f[0], &mut (0..i).map(|c| (r, c))); }
        for c in 0..i { push(&mut f[1], &mut (0..i).map(|r| (r, c))); }
        // (1,-1): fibers r + c = s, increasing r
        for s in 0..2 * i - 1 {
            let r0 = s.saturating_sub(i - 1);
            let r1 = s.min(i - 1);
            push(&mut f[2], &mut (r0..=r1).map(|r| (r, s - r)));
        }
        // (1,1): fibers r - c = d, increasing r
        for d in -(i as isize - 1)..=(i as isize - 1) {
            let r0 = d.max(0) as usize;
            let r1 = (i as isize - 1).min(i as isize - 1 + d) as usize;
            push(&mut f[3], &mut (r0..=r1).map(|r| (r, (r as isize - d) as usize)));
        }
        for fi in f.iter_mut() { fi.1.push(fi.0.len() as u32); }
        GridSorter { i, fibers: f }
    }

    /// pts: (n = I*I) finite points. Returns order (cell -> point index).
    pub fn argsort(&self, pts: &[[f32; 2]], rng: &mut Rng) -> Vec<u32> {
        let n = pts.len();
        let eps = if n <= 1 << 16 { 1e-5 } else { 1e-10 };
        // nan_infs_to_num_and_tiebreaker
        let x: Vec<[f64; 2]> = pts.iter()
            .map(|p| [p[0] as f64 + eps * rng.f64(), p[1] as f64 + eps * rng.f64()]).collect();
        let mut order = vec![0u32; n];
        let mut idx: Vec<u32> = (0..n as u32).collect();
        kdtree(&x, &mut idx, 0, self.i, 0, self.i, 0, self.i, &mut order);

        let keys: [Vec<f64>; 4] = [
            x.iter().map(|p| p[1]).collect(),
            x.iter().map(|p| p[0]).collect(),
            x.iter().map(|p| p[0] - p[1]).collect(),
            x.iter().map(|p| p[0] + p[1]).collect(),
        ];
        for level_dirs in [&[0usize, 1, 2, 3][..], &[0, 1][..]] {
            for _ in 0..40 {
                let mut all = true;
                for &d in level_dirs { all &= self.sort_dir(&mut order, d, &keys[d]); }
                if all { return order; }
            }
        }
        order
    }

    fn sort_dir(&self, order: &mut [u32], d: usize, key: &[f64]) -> bool {
        let (cells, offs) = &self.fibers[d];
        let mut was_sorted = true;
        let mut buf: Vec<u32> = Vec::with_capacity(self.i);
        for w in offs.windows(2) {
            let fib = &cells[w[0] as usize..w[1] as usize];
            buf.clear();
            buf.extend(fib.iter().map(|&c| order[c as usize]));
            // insertion sort: fibers are nearly sorted after the first pass
            let mut changed = false;
            for a in 1..buf.len() {
                let v = buf[a];
                let kv = key[v as usize];
                let mut b = a;
                while b > 0 && key[buf[b - 1] as usize] > kv { buf[b] = buf[b - 1]; b -= 1; }
                if b != a { buf[b] = v; changed = true; }
            }
            if changed {
                was_sorted = false;
                for (k, &c) in fib.iter().enumerate() { order[c as usize] = buf[k]; }
            }
        }
        was_sorted
    }
}

fn next_axis(shape: [usize; 2], start: usize) -> Option<usize> {
    (0..2).map(|o| (start + o) % 2).find(|&a| shape[a] > 1)
}

/// kd-tree initial ordering (gridpoints.kdtree.kdtree_order), as a recursion over grid rectangles.
#[allow(clippy::too_many_arguments)]
fn kdtree(x: &[[f64; 2]], idx: &mut [u32], r0: usize, nr: usize, c0: usize, nc: usize,
          axis: usize, gi: usize, out: &mut [u32]) {
    let shape = [nr, nc];
    if idx.len() == 1 { out[r0 * gi + c0] = idx[0]; return; }
    let Some(axis) = next_axis(shape, axis) else { out[r0 * gi + c0] = idx[0]; return; };
    let n = idx.len();
    let size = shape[axis];
    let slab = n / size;
    let lsize = size / 2;
    let ln = lsize * slab;
    let cmp = |a: &u32, b: &u32| x[*a as usize][axis].partial_cmp(&x[*b as usize][axis]).unwrap();
    let mut child = shape; child[axis] = lsize;
    let nax = next_axis(child, (axis + 1) % 2).unwrap_or(0);
    let rect = |lo: usize, len: usize| if axis == 0 { (r0 + lo, len, c0, nc) } else { (r0, nr, c0 + lo, len) };
    if size % 2 == 0 {
        idx.select_nth_unstable_by(ln, cmp);
        let (l, r) = idx.split_at_mut(ln);
        let (a, b, c, d) = rect(0, lsize); kdtree(x, l, a, b, c, d, nax, gi, out);
        let (a, b, c, d) = rect(lsize, lsize); kdtree(x, r, a, b, c, d, nax, gi, out);
    } else {
        if ln + slab < n { idx.select_nth_unstable_by(ln + slab, cmp); }
        let (lm, r) = idx.split_at_mut(ln + slab);
        if ln > 0 {
            lm.select_nth_unstable_by(ln, cmp);
            let (l, _) = lm.split_at_mut(ln);
            let (a, b, c, d) = rect(0, lsize); kdtree(x, l, a, b, c, d, nax, gi, out);
            let (a, b, c, d) = rect(lsize + 1, lsize); kdtree(x, r, a, b, c, d, nax, gi, out);
        }
        let mut cm = shape; cm[axis] = 1;
        let max = next_axis(cm, (axis + 1) % 2).unwrap_or(0);
        let (a, b, c, d) = rect(lsize, 1);
        kdtree(x, &mut lm[ln..], a, b, c, d, max, gi, out);
    }
}
