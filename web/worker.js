// Owns the Rust/WebAssembly module (src/): stippling, optimisation and the per-dot geometry,
// all off the main thread. The page only draws what this worker sends back.
let ex = null, imageKey = null;
const mem = () => ex.memory.buffer;

// result buffers exported by the wasm module (see buf_ptr in src/lib.rs)
const BUFS = [['pts', Float32Array], ['rgb', Uint8Array], ['cellRgb', Uint8Array], ['ink', Float32Array],
              ['dir', Float32Array], ['verts', Float32Array], ['off', Uint32Array], ['edges', Float32Array]];
const read = (id) => { const [, T] = BUFS[id]; return new T(mem(), ex.buf_ptr(id), ex.buf_len(id)).slice(); };

/** Analyse the current dots and post them with their geometry. */
function postResult(type, extra = {}) {
  ex.analyse();
  const out = {};
  BUFS.forEach(([name], id) => { out[name] = read(id); });
  const { pts, ...geo } = out;
  self.postMessage({ type, pts, geo, ...extra }, Object.values(out).map((a) => a.buffer));
}

function setImage(rgba, S, key) {
  if (key === imageKey) return;
  const p = ex.alloc(rgba.length);
  new Uint8Array(mem(), p, rgba.length).set(rgba);
  ex.set_image(p, S);
  ex.dealloc(p, rgba.length);
  imageKey = key;
}

let onProgress = () => {};
self.onmessage = async ({ data: m }) => {
  if (m.type === 'init') {
    const env = {
      js_progress: (done, total) => onProgress(done / total),
      js_level: (ptr, n) => {          // coarse levels of the stipple, shown as a preview
        const pts = new Float32Array(mem(), ptr, n * 2).slice();
        self.postMessage({ type: 'level', pts }, [pts.buffer]);
      },
    };
    ex = (await WebAssembly.instantiate(m.wasm, { env })).instance.exports;
    self.postMessage({ type: 'ready' });

  } else if (m.type === 'run') {
    const { rgba, S, key, n, iters, seed, invert } = m;
    setImage(rgba, S, key);
    onProgress = (v) => self.postMessage({ type: 'progress', value: v });
    const t0 = performance.now();
    ex.stipple(n, iters, seed >>> 0, invert ? 1 : 0);
    postResult('done', { ms: performance.now() - t0 });

  } else if (m.type === 'refine') {
    // more rgbn steps on the dots of the last stipple, in chunks so the page can redraw in between
    const { steps } = m, CHUNK = 10, t0 = performance.now();
    for (let done = 0; done < steps; done += CHUNK) {
      const k = Math.min(CHUNK, steps - done);
      onProgress = (v) => self.postMessage({ type: 'progress', value: (done + v * k) / steps });
      ex.refine(k);
      if (done + k < steps) postResult('refined');
    }
    postResult('refine-done', { ms: performance.now() - t0 });
  }
};
