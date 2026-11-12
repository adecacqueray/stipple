//! Voronoi cells of the dots, clipped to the square [0, S]², from a Delaunay triangulation.
//!
//! Four far sentinels around the square close every real cell, so each one is a convex polygon
//! that is then clipped to the square (Sutherland–Hodgman). They are 10 S away while every point
//! of the square is within √2 S of a dot, so the sentinels' own cells never reach the square and
//! the clipped diagram is exactly the one of the dots alone.
use delaunator::{next_halfedge, triangulate, Point, EMPTY};

pub struct Cells {
    /// Vertices of every cell, flattened (x, y); cell i spans `off[i]..off[i + 1]` (in vertices).
    pub verts: Vec<f32>,
    pub off: Vec<u32>,
    /// Cell areas (infinite when a dot has no cell, e.g. an exact duplicate).
    pub area: Vec<f64>,
    /// Every Voronoi edge once, clipped to the square: x0, y0, x1, y1.
    pub edges: Vec<f32>,
}

/// `xy`: dots in square pixels (x right, y down), `s`: side of the square.
pub fn cells(xy: &[f64], s: f64) -> Cells {
    let n = xy.len() / 2;
    let mut pts: Vec<Point> = (0..n).map(|i| Point { x: xy[2 * i], y: xy[2 * i + 1] }).collect();
    let (c, k) = (s / 2.0, 10.0 * s);
    for (sx, sy) in [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)] {
        pts.push(Point { x: c + sx * k, y: c + sy * k });
    }
    let tri = triangulate(&pts);
    let nt = tri.triangles.len() / 3;
    let cc: Vec<[f64; 2]> = (0..nt).map(|t| {
        let [a, b, d] = [&pts[tri.triangles[3 * t]], &pts[tri.triangles[3 * t + 1]], &pts[tri.triangles[3 * t + 2]]];
        circumcenter(a, b, d)
    }).collect();

    // an incoming half-edge for each point
    let mut inedge = vec![EMPTY; pts.len()];
    for e in 0..tri.triangles.len() {
        let p = tri.triangles[next_halfedge(e)];
        if tri.halfedges[e] == EMPTY || inedge[p] == EMPTY { inedge[p] = e; }
    }

    let mut out = Cells { verts: Vec::with_capacity(14 * n), off: Vec::with_capacity(n + 1), area: vec![f64::INFINITY; n], edges: vec![] };
    out.off.push(0);
    let mut ring: Vec<[f64; 2]> = Vec::with_capacity(16);
    let mut tmp: Vec<[f64; 2]> = Vec::with_capacity(16);
    for i in 0..n {
        ring.clear();
        let e0 = inedge[i];
        if e0 != EMPTY {
            let mut e = e0;
            let mut guard = 0;
            loop {
                ring.push(cc[e / 3]);
                e = next_halfedge(e);            // now outgoing from i
                e = tri.halfedges[e];            // incoming to i in the next triangle around it
                guard += 1;
                if e == e0 || e == EMPTY || guard > 64 { break; }
            }
            clip_rect(&mut ring, &mut tmp, s);
        }
        if ring.len() >= 3 {
            let mut a = 0.0;
            for k in 0..ring.len() {
                let (p, q) = (ring[k], ring[(k + 1) % ring.len()]);
                a += p[0] * q[1] - q[0] * p[1];
            }
            out.area[i] = a.abs() / 2.0;
            for v in &ring { out.verts.push(v[0] as f32); out.verts.push(v[1] as f32); }
        }
        out.off.push((out.verts.len() / 2) as u32);
    }

    // edges between two real dots: segment between the circumcentres of the two triangles sharing it
    for e in 0..tri.triangles.len() {
        let o = tri.halfedges[e];
        if o == EMPTY || o < e { continue; }
        let (p, q) = (tri.triangles[e], tri.triangles[next_halfedge(e)]);
        if p >= n || q >= n { continue; }
        if let Some([x0, y0, x1, y1]) = clip_segment(cc[e / 3], cc[o / 3], s) {
            out.edges.extend_from_slice(&[x0 as f32, y0 as f32, x1 as f32, y1 as f32]);
        }
    }
    out
}

fn circumcenter(a: &Point, b: &Point, c: &Point) -> [f64; 2] {
    let (dx, dy, ex, ey) = (b.x - a.x, b.y - a.y, c.x - a.x, c.y - a.y);
    let bl = dx * dx + dy * dy;
    let cl = ex * ex + ey * ey;
    let d = dx * ey - dy * ex;
    if d.abs() < 1e-12 * (bl + cl).max(1e-300) {
        return [(a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0];   // degenerate sliver
    }
    let d = 0.5 / d;
    [a.x + (ey * bl - dy * cl) * d, a.y + (dx * cl - ex * bl) * d]
}

/// Clip a convex polygon to [0, s]² in place.
fn clip_rect(poly: &mut Vec<[f64; 2]>, tmp: &mut Vec<[f64; 2]>, s: f64) {
    // each plane: inside if sign * coord <= bound
    for (axis, sign, bound) in [(0usize, -1.0, 0.0), (0, 1.0, s), (1, -1.0, 0.0), (1, 1.0, s)] {
        if poly.is_empty() { return; }
        tmp.clear();
        let inside = |p: &[f64; 2]| sign * p[axis] <= sign * bound;
        for k in 0..poly.len() {
            let (cur, prev) = (poly[k], poly[(k + poly.len() - 1) % poly.len()]);
            let (ci, pi) = (inside(&cur), inside(&prev));
            if ci != pi {
                let t = (bound - prev[axis]) / (cur[axis] - prev[axis]);
                let mut x = [prev[0] + t * (cur[0] - prev[0]), prev[1] + t * (cur[1] - prev[1])];
                x[axis] = bound;
                tmp.push(x);
            }
            if ci { tmp.push(cur); }
        }
        std::mem::swap(poly, tmp);
    }
}

/// Liang–Barsky clip of a segment to [0, s]².
fn clip_segment(a: [f64; 2], b: [f64; 2], s: f64) -> Option<[f64; 4]> {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let (mut t0, mut t1) = (0.0f64, 1.0f64);
    for (p, q) in [(-dx, a[0]), (dx, s - a[0]), (-dy, a[1]), (dy, s - a[1])] {
        if p == 0.0 {
            if q < 0.0 { return None; }
        } else {
            let r = q / p;
            if p < 0.0 { t0 = t0.max(r); } else { t1 = t1.min(r); }
        }
    }
    if t0 >= t1 { return None; }
    Some([a[0] + t0 * dx, a[1] + t0 * dy, a[0] + t1 * dx, a[1] + t1 * dy])
}
