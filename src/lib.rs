//! Stipple: recursive Gaussian blue-noise stippling (port of `blue_sampler.im2points`, method rgbn)
//! and all the geometry the web app draws from (Voronoi cells, stroke directions, colours).
//! The browser side only handles the interface and the drawing.
pub mod blur;
pub mod colour;
pub mod field;
pub mod gridsort;
pub mod image;
pub mod orient;
pub mod rgbn;
pub mod rng;
pub mod voronoi;

/// Everything the drawing needs about one set of dots, in square pixels (x right, y down).
#[derive(Default)]
pub struct Analysis {
    /// Photo colour under each dot (3×3 mean), RGB.
    pub rgb: Vec<u8>,
    /// Photo colour averaged over each Voronoi cell, RGB.
    pub cell_rgb: Vec<u8>,
    /// Mesh ink strength per cell, from its area (0..1).
    pub ink: Vec<f32>,
    /// Unit stroke direction per dot.
    pub dir: Vec<f32>,
    /// Voronoi cells: flattened vertices and per-dot offsets (n + 1 entries).
    pub verts: Vec<f32>,
    pub off: Vec<u32>,
    /// Voronoi edges, x0 y0 x1 y1.
    pub edges: Vec<f32>,
}

/// One image (padded to an S×S square) and its current dots, with the caches that depend on them.
#[derive(Default)]
pub struct Scene {
    pub px: Vec<u8>,
    pub s: usize,
    /// Dots in the unit square, y up.
    pub pts: Vec<[f32; 2]>,
    pub out: Analysis,
    /// The optimiser of the last stipple, kept so "optimise further" can carry on.
    rgbn: Option<rgbn::Rgbn>,
    orient: Option<orient::Orient>,
}

impl Scene {
    pub fn set_image(&mut self, rgba: &[u8], s: usize) {
        *self = Scene { px: rgba.to_vec(), s, ..Scene::default() };
    }

    pub fn stipple(&mut self, n: usize, n_iter: usize, seed: u64, invert: bool, obs: &mut dyn rgbn::Observer) {
        let mut opt = rgbn::Rgbn::new(image::im2field(&self.px, self.s, self.s, invert), seed);
        self.pts = opt.run(n, n_iter, obs);
        self.rgbn = Some(opt);
    }

    /// `steps` more rgbn steps on the current dots (no-op before the first stipple).
    pub fn refine(&mut self, steps: usize, obs: &mut dyn rgbn::Observer) {
        if let Some(opt) = &mut self.rgbn {
            self.pts = opt.refine(std::mem::take(&mut self.pts), steps, obs);
        }
    }

    /// Recompute `out` for the current dots.
    pub fn analyse(&mut self) {
        let s = self.s as f64;
        let xy: Vec<f64> = self.pts.iter().flat_map(|p| [p[0] as f64 * s, (1.0 - p[1] as f64) * s]).collect();
        let c = voronoi::cells(&xy, s);
        let orient = self.orient.get_or_insert_with(|| orient::orient_field(&self.px, self.s));
        self.out = Analysis {
            rgb: colour::under_dots(&self.px, self.s, &xy),
            cell_rgb: colour::over_cells(&self.px, self.s, &xy, &c.verts, &c.off),
            ink: colour::ink_from_area(&c.area),
            dir: orient.at(&xy, s),
            verts: c.verts,
            off: c.off,
            edges: c.edges,
        };
    }
}

// ── WebAssembly interface (raw exports, no wasm-bindgen) ───────────────────────
#[cfg(target_arch = "wasm32")]
mod wasm {
    use super::*;
    #[link(wasm_import_module = "env")]
    extern "C" {
        fn js_progress(done: f64, total: f64);
        fn js_level(ptr: *const f32, n: u32);
    }
    struct Js;
    impl rgbn::Observer for Js {
        fn progress(&mut self, d: f64, t: f64) { unsafe { js_progress(d, t) } }
        fn level_done(&mut self, p: &[[f32; 2]]) { unsafe { js_level(p.as_ptr() as *const f32, p.len() as u32) } }
    }

    static mut SCENE: Option<Scene> = None;
    #[allow(static_mut_refs)]
    fn scene() -> &'static mut Scene { unsafe { SCENE.get_or_insert_with(Scene::default) } }

    #[no_mangle]
    pub extern "C" fn alloc(len: usize) -> *mut u8 {
        let mut v = Vec::<u8>::with_capacity(len);
        let p = v.as_mut_ptr(); std::mem::forget(v); p
    }
    #[no_mangle]
    pub unsafe extern "C" fn dealloc(p: *mut u8, len: usize) { drop(Vec::from_raw_parts(p, 0, len)); }

    /// Copy the S×S RGBA square (already in wasm memory at `ptr`) into the scene.
    #[no_mangle]
    pub unsafe extern "C" fn set_image(ptr: *const u8, s: u32) {
        scene().set_image(std::slice::from_raw_parts(ptr, (s * s * 4) as usize), s as usize);
    }
    /// Stipple the current image; returns the number of dots.
    #[no_mangle]
    pub extern "C" fn stipple(n: u32, n_iter: u32, seed: u32, invert: u32) -> u32 {
        let sc = scene();
        sc.stipple(n as usize, n_iter as usize, seed as u64, invert != 0, &mut Js);
        sc.pts.len() as u32
    }
    /// `steps` more rgbn steps on the current dots; reports progress through js_progress.
    #[no_mangle]
    pub extern "C" fn refine(steps: u32) { scene().refine(steps as usize, &mut Js); }
    #[no_mangle]
    pub extern "C" fn analyse() { scene().analyse(); }

    /// Result buffers: 0 points (f32), 1 rgb (u8), 2 cell rgb (u8), 3 ink (f32), 4 dir (f32),
    /// 5 cell vertices (f32), 6 cell offsets (u32), 7 edges (f32). Lengths are in elements.
    #[no_mangle]
    pub extern "C" fn buf_ptr(id: u32) -> *const u8 {
        let sc = scene(); let o = &sc.out;
        match id {
            0 => sc.pts.as_ptr() as _, 1 => o.rgb.as_ptr(), 2 => o.cell_rgb.as_ptr(), 3 => o.ink.as_ptr() as _,
            4 => o.dir.as_ptr() as _, 5 => o.verts.as_ptr() as _, 6 => o.off.as_ptr() as _, 7 => o.edges.as_ptr() as _,
            _ => std::ptr::null(),
        }
    }
    #[no_mangle]
    pub extern "C" fn buf_len(id: u32) -> u32 {
        let sc = scene(); let o = &sc.out;
        (match id {
            0 => 2 * sc.pts.len(), 1 => o.rgb.len(), 2 => o.cell_rgb.len(), 3 => o.ink.len(),
            4 => o.dir.len(), 5 => o.verts.len(), 6 => o.off.len(), 7 => o.edges.len(), _ => 0,
        }) as u32
    }
}
