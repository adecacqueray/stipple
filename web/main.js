// Interface and drawing only. All the maths (stippling, optimisation, Voronoi cells, colours,
// stroke directions) runs in Rust/WebAssembly inside worker.js; see src/ in the repository.
const ASSETS = { worker: 'web/worker.js', wasm: 'web/stipple.wasm', sample: 'web/zebra.jpg' };
const MAX_SIDE = 1600;              // working resolution sent to the stippler (it resamples to 512² anyway)
const ITERS = 12;                   // passes per level for the initial stipple (blue-sampler's default)
const OPT_STEPS = ITERS * 10;       // rgbn steps per "Optimise further" click (one more pass of the finest level)
const $ = (id) => document.getElementById(id);
const sheet = $('sheet'), ctx = sheet.getContext('2d');

const state = {
  img: null, name: 'zebra', imgId: 0, points: null, preview: null,
  geo: null,      // per-dot geometry from the worker: rgb, cellRgb, ink, dir, verts/off (cells), edges
  crop: null,     // {S, ox, oy, w, h}: image rect inside the padded square, in square pixels
  invert: true, worker: null, ready: null, busy: false, ms: 0, extra: 0, extraMs: 0,
};

// ── controls ────────────────────────────────────────────────────────────────
const sig2 = (x) => { const p = Math.pow(10, Math.floor(Math.log10(x)) - 1); return Math.round(x / p) * p; };
const nFromSlider = (v) => sig2(2000 * Math.pow(64, v / 100));   // 2,000 … 130,000, default 30,000
const fmt = (x) => x.toLocaleString('en-US');
const N = () => nFromSlider(+$('n').value);
const dotScale = () => +$('size').value / 100;
const tint = () => +$('tint').value / 100;
const gap = () => +$('gap').value / 100;
const gain = () => +$('gain').value / 100;
const styleNow = () => document.querySelector('input[name=style]:checked').value;
const SIZE_LABEL = { dots: 'Dot size', mesh: 'Dot size', strokes: 'Stroke width' };
function syncOutputs() {
  $('nOut').textContent = fmt(N());
  $('sizeOut').textContent = dotScale().toFixed(2);
  const t = tint();
  $('tintOut').textContent = t === 0 ? 'Black & white' : t === 1 ? 'Full colour' : t.toFixed(2);
  $('sizeLabel').textContent = SIZE_LABEL[styleNow()];
  $('gapOut').textContent = gap().toFixed(2);
  $('gainOut').textContent = '×' + gain().toFixed(2);
  const mesh = styleNow() === 'mesh'; $('sizeField').hidden = mesh; $('gapField').hidden = !mesh;
}
['n', 'size', 'tint', 'gap', 'gain'].forEach((id) => $(id).addEventListener('input', syncOutputs));
['size', 'tint', 'gap', 'gain'].forEach((id) => $(id).addEventListener('input', () => scheduleDraw()));
document.querySelectorAll('input[name=style]').forEach((r) => r.addEventListener('change', () => { syncOutputs(); scheduleDraw(); }));
document.querySelectorAll('input[name=ink]').forEach((r) => r.addEventListener('change', () => {
  state.invert = $('inkDark').checked; run();
}));
syncOutputs();

// ── worker ──────────────────────────────────────────────────────────────────
let wasmBytes = null;
function startWorker() {
  if (state.worker) state.worker.terminate();
  const w = new Worker(ASSETS.worker);
  state.worker = w;
  state.ready = (async () => {
    wasmBytes ??= await (await fetch(ASSETS.wasm)).arrayBuffer();
    await new Promise((res) => {
      w.onmessage = (e) => { if (e.data.type === 'ready') res(); };
      w.postMessage({ type: 'init', wasm: wasmBytes });
    });
  })();
  return state.ready;
}

// ── image loading ───────────────────────────────────────────────────────────
function loadImage(src, name) {
  return new Promise((res, rej) => {
    const im = new Image();
    im.onload = () => res(im);
    im.onerror = () => rej(new Error(`Couldn't read ${name}. Try a JPEG, PNG or WebP file.`));
    im.src = src;
  }).then((im) => { state.img = im; state.name = name; state.imgId++; $('imgname').textContent = name; });
}

/** Pad the image to a square with the "no ink" colour so dots stay isotropic, then grab RGBA. */
function squarePixels() {
  const im = state.img;
  const k = Math.min(1, MAX_SIDE / Math.max(im.naturalWidth, im.naturalHeight));
  const w = Math.round(im.naturalWidth * k), h = Math.round(im.naturalHeight * k);
  const S = Math.max(w, h), ox = (S - w) / 2, oy = (S - h) / 2;
  const c = new OffscreenCanvas(S, S), g = c.getContext('2d', { willReadFrequently: true });
  g.fillStyle = state.invert ? '#fff' : '#000';
  g.fillRect(0, 0, S, S);
  g.imageSmoothingQuality = 'high';
  g.drawImage(im, ox, oy, w, h);
  state.crop = { S, ox, oy, w, h };
  const rgba = g.getImageData(0, 0, S, S).data;
  return { rgba, S };
}

// ── dots in square-pixel coordinates (x right, y down) ──────────────────────
const cache = new WeakMap();
const geom = (pts) => {
  let xy = cache.get(pts);
  if (!xy) {
    const { S } = state.crop, n = pts.length / 2; xy = new Float64Array(2 * n);
    for (let i = 0; i < n; i++) { xy[2 * i] = pts[2 * i] * S; xy[2 * i + 1] = (1 - pts[2 * i + 1]) * S; }
    cache.set(pts, xy);
  }
  return xy;
};

// ── rendering ───────────────────────────────────────────────────────────────
function sizeSheet() {
  const { w, h } = state.crop;
  const cssW = sheet.getBoundingClientRect().width || 900;
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  sheet.width = Math.round(cssW * dpr);
  sheet.height = Math.round(cssW * dpr * h / w);
}
const inkRGB = () => state.invert ? [[0, 0, 0], [255, 255, 255]] : [[255, 255, 255], [0, 0, 0]];   // [ink, paper]
const css = (r, g, b) => `rgb(${r},${g},${b})`;

/** The photo's hue at full ink strength. Tone is already carried by the dot density (or the cell area in mesh),
 *  so the colour must not carry it a second time: stretch the colour away from the paper until its strongest
 *  channel is full ink. Neutral greys become pure ink (= black & white mode); a light tint and a deep one of the
 *  same hue give the same mark. Near the paper colour the hue is mostly noise, so it fades back to plain ink. */
function fullInk(r, g, b, ink, paper) {
  const dr = r - paper[0], dg = g - paper[1], db = b - paper[2], m = Math.max(Math.abs(dr), Math.abs(dg), Math.abs(db));
  if (m === 0) return ink;
  const s = 255 / m, w = Math.min(1, m / 40);
  return [ink[0] + (paper[0] + dr * s - ink[0]) * w, ink[1] + (paper[1] + dg * s - ink[1]) * w, ink[2] + (paper[2] + db * s - ink[2]) * w];
}

/** Final colour of each point's mark for the current style and colour slider, grouped into buckets of equal colour. */
function groups(pts, geo, style, t, k = gain()) {
  const n = pts.length / 2, [ink, paper] = inkRGB();
  const cv = style === 'mesh' ? geo.ink : null;
  const rgb = t > 0 ? (cv ? geo.cellRgb : geo.rgb) : null;
  const map = new Map();
  for (let i = 0; i < n; i++) {
    let [r, g, b] = ink;
    // hue only, at constant strength, blended in by the colour slider
    if (rgb) { const c = fullInk(rgb[3 * i], rgb[3 * i + 1], rgb[3 * i + 2], ink, paper); r += (c[0] - r) * t; g += (c[1] - g) * t; b += (c[2] - b) * t; }
    // mesh: the tone comes from the cell area, like the density does for dots
    if (cv) { const a = cv[i]; r = paper[0] + (r - paper[0]) * a; g = paper[1] + (g - paper[1]) * a; b = paper[2] + (b - paper[2]) * a; }
    const q = (v) => Math.min(255, Math.round(v * k / 4) * 4);   // contrast: multiply, clip at 255
    const key = css(q(r), q(g), q(b));
    let arr = map.get(key); if (!arr) map.set(key, arr = []); arr.push(i);
  }
  return map;
}
const spacing = (pts) => state.crop.S / Math.sqrt(pts.length / 2);

/** Draw points onto a 2D context of size W×H in the current style (a canvas, or svgContext() for the SVG export). */
function paint(g, W, H, pts, geo, { style, scale, t, preview }) {
  const { ox, oy, w, h } = state.crop, px = W / w, [ink, paper] = inkRGB();
  g.setTransform(1, 0, 0, 1, 0, 0);
  g.fillStyle = css(...paper); g.fillRect(0, 0, W, H);
  g.setTransform(px, 0, 0, px, -ox * px, -oy * px);      // square-pixel coordinates from here on
  const xy = geom(pts), n = pts.length / 2, sp = spacing(pts), minW = 0.35 / px;
  // skip marks lying in the padding around the photo (they would be clipped anyway)
  const inView = (i) => xy[2 * i] >= ox - sp && xy[2 * i] <= ox + w + sp && xy[2 * i + 1] >= oy - sp && xy[2 * i + 1] <= oy + h + sp;
  const all = () => Array.from({ length: n }, (_, i) => i);
  if (preview || style === 'dots') {
    const r = Math.max(minW, scale * 0.5 * sp);
    for (const [col, idx] of preview ? [[preview, all()]] : groups(pts, geo, 'dots', t)) {
      g.fillStyle = col; g.beginPath();
      for (const i of idx) if (inView(i)) { const X = xy[2 * i], Y = xy[2 * i + 1]; g.moveTo(X + r, Y); g.arc(X, Y, r, 0, Math.PI * 2); }
      g.fill();
    }
  } else if (style === 'mesh') {
    const { verts: v, off, edges } = geo;
    for (const [col, idx] of groups(pts, geo, 'mesh', t)) {
      g.fillStyle = col; g.strokeStyle = col; g.lineWidth = 1 / px; g.lineJoin = 'round'; g.beginPath();
      for (const i of idx) {
        const k0 = off[i], k1 = off[i + 1]; if (k1 - k0 < 3 || !inView(i)) continue;
        g.moveTo(v[2 * k0], v[2 * k0 + 1]); for (let k = k0 + 1; k < k1; k++) g.lineTo(v[2 * k], v[2 * k + 1]); g.closePath();
      }
      g.fill(); g.stroke();                                 // 1-device-px stroke in the cell's own colour hides antialiasing seams
    }
    g.strokeStyle = '#000'; g.lineWidth = gap() * sp;   // exact gap width, no minimum
    g.beginPath();
    for (let k = 0; k < edges.length; k += 4) { g.moveTo(edges[k], edges[k + 1]); g.lineTo(edges[k + 2], edges[k + 3]); }
    g.stroke();
  } else {
    const d = geo.dir, half = 0.75 * sp;
    g.lineWidth = Math.max(minW, scale * 0.6 * sp); g.lineCap = 'round';
    for (const [col, idx] of groups(pts, geo, 'strokes', t)) {
      g.strokeStyle = col; g.beginPath();
      for (const i of idx) if (inView(i)) { const X = xy[2 * i], Y = xy[2 * i + 1], dx = d[2 * i] * half, dy = d[2 * i + 1] * half; g.moveTo(X - dx, Y - dy); g.lineTo(X + dx, Y + dy); }
      g.stroke();
    }
  }
  g.setTransform(1, 0, 0, 1, 0, 0);
}

/** Stand-in for the bits of the canvas API that paint() uses; records the drawing as SVG instead. */
function svgContext(W, H) {
  const out = [], f = (v) => +v.toFixed(2);
  let d = '', group = false;
  return {
    fillStyle: '#000', strokeStyle: '#000', lineWidth: 1, lineJoin: 'miter', lineCap: 'butt',
    setTransform(...m) {
      if (group) out.push('</g>');
      group = m.join() !== '1,0,0,1,0,0';
      if (group) out.push(`<g transform="matrix(${m.map(f).join(' ')})">`);
    },
    fillRect(x, y, w, h) { out.push(`<rect x="${f(x)}" y="${f(y)}" width="${f(w)}" height="${f(h)}" fill="${this.fillStyle}"/>`); },
    beginPath() { d = ''; },
    moveTo(x, y) { d += `M${f(x)} ${f(y)}`; },
    lineTo(x, y) { d += `L${f(x)} ${f(y)}`; },
    closePath() { d += 'Z'; },
    arc(x, y, r) { r = f(r); d += `a${r} ${r} 0 1 0 ${-2 * r} 0a${r} ${r} 0 1 0 ${2 * r} 0`; },   // full circle; paint() moves to (x+r, y) first
    fill() { if (d) out.push(`<path fill="${this.fillStyle}" d="${d}"/>`); },
    stroke() { if (d) out.push(`<path fill="none" stroke="${this.strokeStyle}" stroke-width="${f(this.lineWidth)}" stroke-linejoin="${this.lineJoin}" stroke-linecap="${this.lineCap}" d="${d}"/>`); },
    toString() {
      return [`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}">`, ...out, group ? '</g>' : '', '</svg>'].join('\n');
    },
  };
}
const opts = () => ({ style: styleNow(), scale: dotScale(), t: tint() });

function draw() {
  if (!state.crop) return;
  if (state.points) return paint(ctx, sheet.width, sheet.height, state.points, state.geo, opts());
  if (state.preview) {
    const npb = getComputedStyle(document.documentElement).getPropertyValue('--npb').trim();
    paint(ctx, sheet.width, sheet.height, state.preview, null, { ...opts(), preview: npb });
  }
}
let drawQueued = false;
function scheduleDraw() { if (drawQueued) return; drawQueued = true; requestAnimationFrame(() => { drawQueued = false; draw(); }); }
new ResizeObserver(() => { if (state.crop) { sizeSheet(); draw(); } }).observe(sheet);

// ── run ─────────────────────────────────────────────────────────────────────
const setBusy = (b) => { state.busy = b; ['svg', 'png', 'opt'].forEach((id) => $(id).disabled = b || !state.points); };
const doneCaption = () => {
  const n = state.points.length / 2;
  let s = `${fmt(n)} dots, placed in ${(state.ms / 1000).toFixed(2)} s`;
  if (state.extra) s += ` · ${fmt(state.extra)} extra steps (${(state.extraMs / 1000).toFixed(1)} s)`;
  return s + '.';
};
function fail(e) {
  setBusy(false);
  $('caption').classList.add('error');
  $('caption').textContent = `The stippler stopped: ${e.message || 'unknown error'}. Try fewer dots.`;
}
// identifies the padded square sent to the worker (image + paper colour), so it is only copied in when it changes
const keyNow = () => `${state.imgId}|${state.invert}`;
async function run() {
  if (!state.img) return;
  if (state.busy) startWorker();            // a new request replaces the running one
  setBusy(true);
  $('caption').classList.remove('error');
  $('caption').textContent = 'Placing dots…';
  $('bar').style.width = '0%';
  state.points = null; state.geo = null; state.preview = null; state.extra = 0; state.extraMs = 0;
  const { rgba, S } = squarePixels();
  sizeSheet();
  ctx.fillStyle = css(...inkRGB()[1]); ctx.fillRect(0, 0, sheet.width, sheet.height);
  await state.ready;
  const n = N();
  const w = state.worker;
  w.onmessage = ({ data: m }) => {
    if (w !== state.worker) return;
    if (m.type === 'progress') $('bar').style.width = (100 * m.value).toFixed(1) + '%';
    else if (m.type === 'level' && m.pts.length / 2 < n) { state.preview = m.pts; scheduleDraw(); }
    else if (m.type === 'done') {
      state.points = m.pts; state.geo = m.geo; state.ms = m.ms; setBusy(false); draw();
      $('bar').style.width = '100%';
      $('caption').textContent = doneCaption();
    }
  };
  w.onerror = fail;
  w.postMessage({ type: 'run', rgba, S, key: keyNow(), n, iters: ITERS,
                  seed: crypto.getRandomValues(new Uint32Array(1))[0], invert: state.invert }, [rgba.buffer]);
}
$('go').addEventListener('click', run);

async function optimise() {
  if (!state.points || state.busy) return;
  setBusy(true);
  $('caption').classList.remove('error');
  $('caption').textContent = `Optimising… 0 / ${OPT_STEPS} steps`;
  $('bar').style.width = '0%';
  await state.ready;
  const w = state.worker;
  w.onmessage = ({ data: m }) => {
    if (w !== state.worker) return;
    if (m.type === 'progress') {
      $('bar').style.width = (100 * m.value).toFixed(1) + '%';
      $('caption').textContent = `Optimising… ${Math.round(m.value * OPT_STEPS)} / ${OPT_STEPS} steps`;
    } else if (m.type === 'refined') { state.points = m.pts; state.geo = m.geo; scheduleDraw(); }
    else if (m.type === 'refine-done') {
      state.points = m.pts; state.geo = m.geo; state.extra += OPT_STEPS; state.extraMs += m.ms; setBusy(false); draw();
      $('bar').style.width = '100%';
      $('caption').textContent = doneCaption();
    }
  };
  w.onerror = fail;
  w.postMessage({ type: 'refine', steps: OPT_STEPS });   // the worker still holds the image and the dots
}
$('opt').addEventListener('click', optimise);

// ── input ───────────────────────────────────────────────────────────────────
async function useFile(f) {
  if (!f) return;
  try { await loadImage(URL.createObjectURL(f), f.name); run(); }
  catch (e) { $('caption').classList.add('error'); $('caption').textContent = e.message; }
}
$('file').addEventListener('change', (e) => useFile(e.target.files[0]));
const plate = $('plate');
plate.addEventListener('dragover', (e) => { e.preventDefault(); plate.classList.add('dragging'); });
plate.addEventListener('dragleave', () => plate.classList.remove('dragging'));
plate.addEventListener('drop', (e) => { e.preventDefault(); plate.classList.remove('dragging'); useFile(e.dataTransfer.files[0]); });

// ── export ──────────────────────────────────────────────────────────────────
function download(blob, name) {
  const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
const base = () => state.name.replace(/\.[^.]+$/, '') + '-' + styleNow();
$('svg').addEventListener('click', () => {
  const { w, h } = state.crop, g = svgContext(w, h);
  paint(g, w, h, state.points, state.geo, opts());
  download(new Blob([String(g)], { type: 'image/svg+xml' }), base() + '.svg');
});
$('png').addEventListener('click', async () => {
  const { w, h } = state.crop, k = 3000 / Math.max(w, h);
  const c = new OffscreenCanvas(Math.round(w * k), Math.round(h * k));
  paint(c.getContext('2d'), c.width, c.height, state.points, state.geo, opts());
  download(await c.convertToBlob({ type: 'image/png' }), base() + '.png');
});

// ── start: stipple the sample right away ────────────────────────────────────
startWorker();
loadImage(ASSETS.sample, 'Sample: zebra')
  .then(() => { $('go').disabled = false; run(); })
  .catch((e) => { $('caption').classList.add('error'); $('caption').textContent = e.message; });
