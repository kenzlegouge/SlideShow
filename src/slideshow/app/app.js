'use strict';
// SlideShow viewer. Points are drawn with WebGL2 (millions of cells); D3 does the zoom, the overlay labels, the
// tables and the gene chart. Marker statistics are computed in the browser from per-(sample x cluster) sufficient
// statistics (stats.f32), never from the cell matrix. Bundle format: SPEC.md (format_version 1).

const BUNDLE = 'bundle/';
const FAST_N = 400000;          // points drawn while panning / zooming (cells are stored in random order)
const MIN_REP_CELLS = 20;       // a sample contributes a per-sample fold change only with >= 20 cells on each side
const GRID = 400, GRID_EXT = 1.2;
// 8 validated hues first, then extension hues. With ~30 classes colour alone cannot carry identity:
// the on-map labels, the hover tooltip and the cluster list do.
const PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948',
  '#17becf', '#8c564b', '#bcbd22', '#9467bd', '#393b79', '#ff9896', '#637939', '#c49c94',
  '#7b4173', '#98df8a', '#ad494a', '#5254a3', '#ffbb78', '#8ca252', '#a55194', '#9edae5',
  '#843c39', '#c5b0d5', '#bd9e39', '#aec7e8', '#8c6d31', '#f7b6d2', '#dbdb8d', '#6b6ecf'];
const CAT_LABEL = {};            // category name -> menu label, from meta.json
const FORMAT_VERSION = 1;

const $ = s => document.querySelector(s);
const fmtInt = d3.format(',');
const pct = v => (v * 100).toFixed(v < 0.0995 ? 1 : 0) + '%';
const esc = s => String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const cssVar = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const rgb = c => { const k = d3.color(c).rgb(); return [k.r, k.g, k.b]; };

const S = {
  meta: null, N: 0, genes: [], nG: 0, nCl: 0,
  cats: {}, catData: {}, contData: new Map(),
  emb: [], embIdx: 0, embFrom: 0, mix: 0,
  mode: { type: 'cat', name: 'cluster' }, hidden: new Set(),
  selected: new Set(), labels: {},
  keep: null, filters: new Map(), filterCounts: {},      // cell filter: keep = Uint8Array(N) of 0 / 1 or null; filters = Map(category -> Set of levels switched off)
  t: d3.zoomIdentity, interacting: false, dragging: false, W: 1, H: 1,      // dragging: a panel splitter is moved (layout.js)
  stats: null, markers: [], sort: { key: 'lfc', dir: -1 }, pickGrid: {},
};
let gl, prog, U = {}, buf = {}, tex = {}, zoom, svg, renderQueued = false;

// ---------------------------------------------------------------- loading
async function fetchBin(file, Type) {
  const r = await fetch(BUNDLE + file);
  if (!r.ok) throw new Error(`${file}: HTTP ${r.status}`);
  return new Type(await r.arrayBuffer());
}

async function init() {
  const load = $('#loading');
  try {
    const r = await fetch(BUNDLE + 'meta.json');
    if (!r.ok) throw new Error('bundle/meta.json not found. Build the bundle first (slideshow build).');
    const m = S.meta = await r.json();
    if (m.format !== 'slideshow-bundle' || m.format_version !== FORMAT_VERSION)
      throw new Error(`This viewer reads SlideShow bundles of format version ${FORMAT_VERSION}; the bundle is ${m.format} version ${m.format_version}. Rebuild it or update slideshow.`);
    for (const c of m.categories) CAT_LABEL[c.name] = c.label || c.name;
    S.rep = m.sample.label;                                         // what a replicate is called: 'sample', 'core', 'patient' ...
    S.N = m.n_cells; S.genes = m.genes; S.nG = m.genes.length;
    for (const c of m.categories) S.cats[c.name] = [].concat(c.levels);
    S.nCl = S.cats.cluster.length;
    S.nS = m.samples.length; S.nGr = m.groups.n.length;
    load.textContent = `Loading ${fmtInt(S.N)} cells…`;
    const embs = m.embeddings;                                  // any number of 2-D embeddings (UMAP, PaCMAP, t-SNE ...)
    const [cl, stats, ...ed] = await Promise.all([
      fetchBin('cat_cluster.u8', Uint8Array),
      fetchBin('stats.f32', Float32Array),
      ...embs.map(e => fetchBin(e.file, Float32Array))]);
    if (ed.some(e => e.length !== 2 * S.N) || cl.length !== S.N) throw new Error('Bundle files do not match meta.json (n_cells).');
    if (stats.length !== 2 * S.nGr * S.nG) throw new Error('stats.f32 does not match meta.json (groups x genes).');
    S.emb = embs.map((e, i) => ({ ...e, data: ed[i] }));
    S.catData.cluster = cl; S.stats = stats;
    $('#subtitle').textContent = `${m.title} · ${fmtInt(S.N)} cells · ${S.nG} genes · ${S.nCl} clusters (${m.cluster.key}) · ${S.nS} ${plural(S.rep)}`;
    document.title = `${m.title} · SlideShow`;
    await loadLabels();
    initGL(); initUI(); resize();
    if (window.SP) await SP.init();              // tissue view, when the bundle has spatial files
    await setMode({ type: 'cat', name: 'cluster' });
    load.remove();
    applyLayout();                               // the map pane no longer has to show the loading message (layout.js)
    // deep link: #c=3,7&gene=EPCAM opens with clusters 3 and 7 selected and the map coloured by EPCAM
    const h = new URLSearchParams(location.hash.slice(1));
    (h.get('c') || '').split(',').map(id => S.cats.cluster.indexOf(id)).filter(i => i >= 0).forEach(i => select(i, true));
    if (S.genes.includes(h.get('gene'))) await setMode({ type: 'gene', j: S.genes.indexOf(h.get('gene')) });
    if (S.cats[h.get('color')]) await setMode({ type: 'cat', name: h.get('color') });          // #color=leiden, #color=level1 ...
    if (h.get('color') === 'label') await setMode({ type: 'label' });                             // #color=label : my labels
    const hide = (h.get('hide') || '').split(',').filter(Boolean);                                // #hide=name1,name2 : legend entries hidden
    if (hide.length && !isCont() && !(S.mode.type === 'cat' && S.mode.name === 'cluster')) {
      const lab = i => (S.labels[S.cats.cluster[i]] || {}).label || 'unlabelled';
      if (S.mode.type === 'label') d3.range(S.nCl).filter(i => hide.includes(lab(i))).forEach(i => S.hidden.add(i));
      else hide.map(n => S.cats[S.mode.name].indexOf(n)).filter(i => i >= 0).forEach(i => S.hidden.add(i));
      updateLUT(); renderLegend(); requestRender();
    }
    const e0 = S.emb.findIndex(e => e.name === h.get('emb'));                                     // #emb=pacmap : embedding shown
    if (e0 > 0) setEmbedding(e0);
    if (window.SP && h.get('sample')) SP.openSample(h.get('sample'));                             // #sample=<id> opens that sample in the tissue view
  } catch (e) {
    load.textContent = 'Could not load the dataset.\n' + e.message;
    console.error(e);
  }
}

// ---------------------------------------------------------------- WebGL
const VS = `#version 300 es
in vec2 aP0; in vec2 aP1; in float aCl; in float aVal; in float aSel; in float aKeep;
uniform float uMix, uSize; uniform vec4 uXf, uDim; uniform int uPass, uCont, uUseSel;
uniform sampler2D uLUT, uMask;
out vec4 vCol;
void main() {
  vec2 p = mix(aP0, aP1, uMix);
  vec4 col = texelFetch(uLUT, ivec2(int(aVal + 0.5), 0), 0);
  float sel = texelFetch(uMask, ivec2(int(aCl + 0.5), 0), 0).r;
  if (uUseSel == 1) sel *= aSel;                 // cell-level selection (lasso in the tissue view)
  // foreground = selected cluster and, for a gene / covariate, a non-zero value. Drawn last, on top.
  bool fg = sel > 0.5 && (uCont == 0 || aVal > 0.5);
  bool drop = col.a < 0.01 || aKeep < 0.5 || (uPass == 0 && fg) || (uPass == 1 && !fg);       // aKeep: cell filter
  vCol = sel > 0.5 ? col : uDim;
  gl_Position = drop ? vec4(2.0, 2.0, 2.0, 1.0) : vec4(p.x * uXf.x + uXf.y, p.y * uXf.z + uXf.w, 0.0, 1.0);
  gl_PointSize = uSize;
}`;
const FS = `#version 300 es
precision mediump float;
in vec4 vCol; uniform float uAlpha; out vec4 o;
void main() {
  vec2 d = gl_PointCoord - 0.5;
  if (dot(d, d) > 0.25) discard;
  o = vec4(vCol.rgb, vCol.a * uAlpha);
}`;

function initGL() {
  const canvas = $('#gl');
  gl = canvas.getContext('webgl2', { alpha: false, antialias: false });
  if (!gl) throw new Error('WebGL2 is not available in this browser.');
  const sh = (type, src) => {
    const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
    return s;
  };
  prog = gl.createProgram();
  gl.attachShader(prog, sh(gl.VERTEX_SHADER, VS)); gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, FS));
  gl.linkProgram(prog);
  if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog));
  gl.useProgram(prog);
  for (const u of ['uMix', 'uSize', 'uXf', 'uDim', 'uPass', 'uCont', 'uLUT', 'uMask', 'uAlpha', 'uUseSel']) U[u] = gl.getUniformLocation(prog, u);
  gl.bindVertexArray(gl.createVertexArray());
  const attr = (name, data, size, type, usage) => {
    const b = gl.createBuffer(), loc = gl.getAttribLocation(prog, name);
    gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.bufferData(gl.ARRAY_BUFFER, data, usage || gl.STATIC_DRAW);
    gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, size, type, false, 0, 0);
    return b;
  };
  buf.p0 = attr('aP0', S.emb[0].data, 2, gl.FLOAT, gl.DYNAMIC_DRAW);      // embedding shown before a switch
  buf.p1 = attr('aP1', S.emb[0].data, 2, gl.FLOAT, gl.DYNAMIC_DRAW);      // embedding shown after it (uMix goes 0 -> 1)
  attr('aCl', S.catData.cluster, 1, gl.UNSIGNED_BYTE);
  buf.val = attr('aVal', S.catData.cluster, 1, gl.UNSIGNED_BYTE, gl.DYNAMIC_DRAW);
  buf.sel = attr('aSel', new Uint8Array(S.N).fill(1), 1, gl.UNSIGNED_BYTE, gl.DYNAMIC_DRAW);
  buf.keep = attr('aKeep', new Uint8Array(S.N).fill(1), 1, gl.UNSIGNED_BYTE, gl.DYNAMIC_DRAW);
  const mkTex = unit => {
    const t = gl.createTexture(); gl.activeTexture(gl.TEXTURE0 + unit); gl.bindTexture(gl.TEXTURE_2D, t);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    return t;
  };
  tex.lut = mkTex(0); tex.mask = mkTex(1);
  gl.uniform1i(U.uLUT, 0); gl.uniform1i(U.uMask, 1);
  gl.enable(gl.BLEND); gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
  gl.disable(gl.DEPTH_TEST);
  updateMask();
}

function setTex(unit, t, data) {
  gl.activeTexture(gl.TEXTURE0 + unit); gl.bindTexture(gl.TEXTURE_2D, t);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 256, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, data);
}

function updateMask() {
  const a = new Uint8Array(1024);
  for (let i = 0; i < 256; i++) a[i * 4] = (S.selected.size === 0 || isCont() || S.selected.has(i)) ? 255 : 0;   // genes and covariates are never dimmed
  setTex(1, tex.mask, a);
  if (window.SP) SP.onMask(a);
}

// cell-level selection shared by the two views: mask = Uint8Array(N) of 0 / 1, or null to clear
function setCellSel(mask) {
  S.cellSel = mask || null;
  if (mask) { gl.bindBuffer(gl.ARRAY_BUFFER, buf.sel); gl.bufferData(gl.ARRAY_BUFFER, mask, gl.DYNAMIC_DRAW); }
  requestRender();
  if (window.SP) SP.onCellSel();
}

// ---------------------------------------------------------------- cell filter (panel drawer)
// Cells outside the filter are not drawn in either view and cannot be picked or lassoed. The marker table, the composition,
// the gene chart and the cluster counts come from sums precomputed over ALL cells (stats.f32) and ignore the filter.
const catCodes = async name => S.catData[name] || (S.catData[name] = await fetchBin(`cat_${name}.u8`, Uint8Array));

async function applyFilters() {
  const tok = S.filterTok = (S.filterTok || 0) + 1;
  const active = [...S.filters].filter(([, off]) => off.size);
  await Promise.all(active.map(([name]) => catCodes(name)));
  if (tok !== S.filterTok) return;                                    // a newer click superseded this one
  let mask = null, shown = S.N;
  if (active.length) {
    mask = new Uint8Array(S.N).fill(1);
    for (const [name, off] of active) {
      const codes = S.catData[name], drop = new Uint8Array(256);
      off.forEach(i => { drop[i] = 1; });
      for (let i = 0; i < S.N; i++) if (drop[codes[i]]) mask[i] = 0;
    }
    shown = d3.sum(mask);
  }
  S.keep = mask;
  gl.bindBuffer(gl.ARRAY_BUFFER, buf.keep); gl.bufferData(gl.ARRAY_BUFFER, mask || new Uint8Array(S.N).fill(1), gl.DYNAMIC_DRAW);
  if (window.SP) SP.onKeep();
  $('#filtercount').textContent = active.length ? `${fmtInt(shown)} of ${fmtInt(S.N)} cells shown (${pct(shown / S.N)})` : '';
  document.body.classList.toggle('filtered', active.length > 0);
  requestRender();
}

function renderFilters() {
  const names = [...S.filters.keys()];
  d3.select('#filtervar').selectAll('option').data(['', ...S.meta.categories.map(c => c.name).filter(v => v !== 'cluster' && !S.filters.has(v)), 'cluster']).join('option')
    .attr('value', d => d).text(d => d ? (CAT_LABEL[d] || d) : 'add a filter on…');
  $('#filtervar').value = '';
  const blocks = d3.select('#filters').selectAll('div.fblock').data(names, d => d).join(en => {
    const b = en.append('div').attr('class', 'fblock'), h = b.append('div').attr('class', 'fhead');
    h.append('b'); h.append('span').attr('class', 'spacer');
    h.append('button').attr('class', 'fall').text('all'); h.append('button').attr('class', 'fnone').text('none');
    h.append('button').attr('class', 'fx').attr('title', 'Remove this filter').text('×');
    b.append('div').attr('class', 'flevels');
    return b;
  });
  blocks.select('b').text(d => CAT_LABEL[d] || d);
  blocks.select('.fall').on('click', (e, d) => { S.filters.get(d).clear(); renderFilters(); applyFilters(); });
  blocks.select('.fnone').on('click', (e, d) => { S.cats[d].forEach((lv, i) => S.filters.get(d).add(i)); renderFilters(); applyFilters(); });
  blocks.select('.fx').on('click', (e, d) => { S.filters.delete(d); renderFilters(); applyFilters(); });
  blocks.select('.flevels').each(function (name) {
    const off = S.filters.get(name), n = S.filterCounts[name];
    const rows = d3.select(this).selectAll('label').data(S.cats[name].map((lv, i) => ({ lv, i })), d => d.i).join(en => {
      const l = en.append('label'); l.append('input').attr('type', 'checkbox'); l.append('span').attr('class', 'nm'); l.append('span').attr('class', 'muted');
      return l;
    });
    rows.select('input').property('checked', d => !off.has(d.i)).on('change', function (e, d) { this.checked ? off.delete(d.i) : off.add(d.i); applyFilters(); });
    rows.select('.nm').text(d => d.lv);
    rows.select('.muted').text(d => pct(n[d.i] / S.N));
  });
}

async function addFilter(name) {
  const codes = await catCodes(name), n = new Float64Array(S.cats[name].length);
  for (let i = 0; i < S.N; i++) n[codes[i]]++;
  S.filterCounts[name] = n; S.filters.set(name, new Set());
  renderFilters();
}

const isCont = () => S.mode.type === 'cont' || S.mode.type === 'gene';
const ramp = () => d3.interpolateRgbBasis(matchMedia('(prefers-color-scheme: dark)').matches
  ? ['#0d366b', '#1c5cab', '#3987e5', '#86b6ef', '#cde2fb'] : ['#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b']);

// colour of each level of the current categorical colouring
function levelColors() {
  const m = S.mode;
  if (m.type === 'label') return S.cats.cluster.map((id, i) => labelColor(i));
  return S.cats[m.name].map((lv, i) => (lv === 'isolated' || lv === 'NA' || lv === 'Unknown') ? cssVar('--dim') : PAL[i % PAL.length]);
}
function labelColor(i) {
  const lab = (S.labels[S.cats.cluster[i]] || {}).label;
  if (!lab) return cssVar('--dim');
  return PAL[usedLabels().indexOf(lab) % PAL.length];
}
const usedLabels = () => [...new Set(S.cats.cluster.map(id => (S.labels[id] || {}).label).filter(Boolean))];
const clusterColor = i => S.cats.cluster[i] === 'isolated' ? cssVar('--dim') : PAL[i % PAL.length];

function updateLUT() {
  const a = new Uint8Array(1024);
  if (isCont()) {
    const z = rgb(cssVar('--zero')), r = ramp();
    a.set([z[0], z[1], z[2], 150], 0);
    for (let i = 1; i < 256; i++) a.set([...rgb(r((i - 1) / 254)), 255], i * 4);
  } else {
    levelColors().forEach((c, i) => a.set([...rgb(c), S.hidden.has(i) ? 0 : 255], i * 4));
  }
  setTex(0, tex.lut, a);
  if (window.SP) SP.onLUT(a);
}

function requestRender() { if (!renderQueued) { renderQueued = true; requestAnimationFrame(draw); } }

function draw() {
  renderQueued = false;
  if ($('#map').hidden) return;                  // switched off in the panel drawer; resize() redraws when it comes back
  const { W, H, t } = S, k0 = Math.min(W, H) / 2 * 0.95;
  const bg = rgb(cssVar('--surface')), dim = rgb(cssVar('--dim'));
  gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
  gl.clearColor(bg[0] / 255, bg[1] / 255, bg[2] / 255, 1); gl.clear(gl.COLOR_BUFFER_BIT);
  gl.uniform4f(U.uXf, 2 * k0 * t.k / W, 2 * (t.k * W / 2 + t.x) / W - 1, 2 * k0 * t.k / H, 1 - 2 * (t.k * H / 2 + t.y) / H);
  gl.uniform1f(U.uMix, S.mix);
  const dpr = window.devicePixelRatio || 1;
  gl.uniform1f(U.uSize, Math.max(1, Math.min(14, +$('#psize').value * Math.sqrt(t.k))) * dpr);
  gl.uniform4f(U.uDim, dim[0] / 255, dim[1] / 255, dim[2] / 255, 0.6);
  gl.uniform1i(U.uCont, isCont() ? 1 : 0); gl.uniform1i(U.uUseSel, S.cellSel ? 1 : 0);
  const alpha = Math.min(1, +$('#palpha').value * Math.sqrt(t.k));
  const n = S.interacting || S.dragging ? Math.min(S.N, FAST_N) : S.N;
  if (isCont() || S.selected.size || S.cellSel) { // background: unselected clusters, or zero-expression cells
    gl.uniform1i(U.uPass, 0); gl.uniform1f(U.uAlpha, alpha * 0.5); gl.drawArrays(gl.POINTS, 0, n);
  }
  gl.uniform1i(U.uPass, 1); gl.uniform1f(U.uAlpha, alpha); gl.drawArrays(gl.POINTS, 0, n);
}

function resize() {
  const el = $('#map'), dpr = window.devicePixelRatio || 1;
  if (!el.clientWidth || !el.clientHeight) return;            // hidden: keep the last size and view
  // keep the same embedding point in the middle, at the same zoom factor, when the pane changes size
  const centre = S.W > 1 && (el.clientWidth !== S.W || el.clientHeight !== S.H) ? toData(S.W / 2, S.H / 2) : null;
  S.W = el.clientWidth; S.H = el.clientHeight;
  gl.canvas.width = Math.round(S.W * dpr); gl.canvas.height = Math.round(S.H * dpr);
  if (centre && svg) {
    const [x, y] = toScreen0(centre[0], centre[1]), k = S.t.k;
    svg.call(zoom.transform, d3.zoomIdentity.translate(S.W / 2 - x * k, S.H / 2 - y * k).scale(k));
  }
  positionLabels(); requestRender();
}

// ---------------------------------------------------------------- coordinates, picking, zoom
const k0 = () => Math.min(S.W, S.H) / 2 * 0.95;
const toScreen0 = (x, y) => [x * k0() + S.W / 2, -y * k0() + S.H / 2];        // before the zoom transform
function toData(sx, sy) {
  const [x0, y0] = S.t.invert([sx, sy]);
  return [(x0 - S.W / 2) / k0(), -(y0 - S.H / 2) / k0()];
}
const embName = () => S.emb[S.embIdx].name;

// majority cluster per cell of a GRID x GRID raster of the embedding, built on first use
function grid() {
  const e = S.emb[S.embIdx];
  if (S.pickGrid[e.name]) return S.pickGrid[e.name];
  const cnt = new Uint16Array(GRID * GRID * S.nCl), cl = S.catData.cluster, d = e.data, s = GRID / (2 * GRID_EXT);
  for (let i = 0; i < S.N; i++) {
    const gx = Math.floor((d[2 * i] + GRID_EXT) * s), gy = Math.floor((d[2 * i + 1] + GRID_EXT) * s);
    if (gx < 0 || gy < 0 || gx >= GRID || gy >= GRID) continue;
    const j = (gy * GRID + gx) * S.nCl + cl[i];
    if (cnt[j] < 65535) cnt[j]++;
  }
  const top = new Uint8Array(GRID * GRID).fill(255);
  for (let b = 0; b < GRID * GRID; b++) {
    let best = 0, bi = 255;
    for (let c = 0; c < S.nCl; c++) { const v = cnt[b * S.nCl + c]; if (v > best) { best = v; bi = c; } }
    top[b] = bi;
  }
  return (S.pickGrid[e.name] = top);
}

function pick(sx, sy) {
  const g = grid(), [x, y] = toData(sx, sy), s = GRID / (2 * GRID_EXT);
  const gx = Math.floor((x + GRID_EXT) * s), gy = Math.floor((y + GRID_EXT) * s);
  const binPx = k0() * S.t.k / s, R = Math.max(0, Math.min(3, Math.ceil(6 / binPx)));
  let best = -1, bd = Infinity;
  for (let dy = -R; dy <= R; dy++) for (let dx = -R; dx <= R; dx++) {
    const X = gx + dx, Y = gy + dy;
    if (X < 0 || Y < 0 || X >= GRID || Y >= GRID) continue;
    const c = g[Y * GRID + X], dd = dx * dx + dy * dy;
    if (c !== 255 && dd < bd) { bd = dd; best = c; }
  }
  return best;
}

function zoomToSelection() {
  if (!S.selected.size || $('#map').hidden) return;
  const key = 'box_' + embName();
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (const i of S.selected) {
    const b = S.meta.clusters[i][key];
    x0 = Math.min(x0, b[0]); x1 = Math.max(x1, b[1]); y0 = Math.min(y0, b[2]); y1 = Math.max(y1, b[3]);
  }
  const [ax, ay] = toScreen0(x0, y1), [bx, by] = toScreen0(x1, y0);
  const k = Math.max(1, Math.min(40, 0.8 * Math.min(S.W / Math.max(1, bx - ax), S.H / Math.max(1, by - ay))));
  const tr = d3.zoomIdentity.translate(S.W / 2, S.H / 2).scale(k).translate(-(ax + bx) / 2, -(ay + by) / 2);
  svg.transition().duration(600).call(zoom.transform, tr);
}

function positionLabels() {
  if (!svg) return;
  const a = S.emb[S.embFrom].name, b = S.emb[S.embIdx].name, show = $('#showlabels').checked;
  svg.selectAll('text.cl')
    .attr('display', (d, i) => show && d.id !== 'isolated' && (!S.selected.size || S.selected.has(i)) ? null : 'none')
    .classed('sel', (d, i) => S.selected.has(i))
    .attr('transform', d => {
      const x = d['cx_' + a] * (1 - S.mix) + d['cx_' + b] * S.mix, y = d['cy_' + a] * (1 - S.mix) + d['cy_' + b] * S.mix;
      const p = S.t.apply(toScreen0(x, y));
      return `translate(${p[0]},${p[1]})`;
    })
    .text(d => { const l = (S.labels[d.id] || {}).label; return l ? `${d.id} · ${l}` : d.id; });
}

function setEmbedding(idx) {
  if (idx === S.embIdx || !S.emb[idx]) return;
  if (S.embTimer) S.embTimer.stop();                       // a switch during a switch starts from the previous target
  S.embFrom = S.embIdx; S.embIdx = idx; S.mix = 0;
  gl.bindBuffer(gl.ARRAY_BUFFER, buf.p0); gl.bufferData(gl.ARRAY_BUFFER, S.emb[S.embFrom].data, gl.DYNAMIC_DRAW);
  gl.bindBuffer(gl.ARRAY_BUFFER, buf.p1); gl.bufferData(gl.ARRAY_BUFFER, S.emb[idx].data, gl.DYNAMIC_DRAW);
  d3.selectAll('#embtoggle button').classed('on', (d, i) => i === idx);
  const ease = d3.easeCubicInOut;
  S.interacting = true;
  const timer = S.embTimer = d3.timer(el => {
    const u = Math.min(1, el / 900);
    S.mix = ease(u);
    positionLabels(); requestRender();
    if (u >= 1) { timer.stop(); S.interacting = false; requestRender(); if (S.selected.size) zoomToSelection(); }
  });
}

// ---------------------------------------------------------------- colouring
async function setMode(mode) {
  let data;
  if (mode.type === 'label') data = S.catData.cluster;
  else if (mode.type === 'cat') data = S.catData[mode.name] || (S.catData[mode.name] = await fetchBin(`cat_${mode.name}.u8`, Uint8Array));
  else {
    const file = mode.type === 'gene' ? `expr/${mode.j}.u8` : S.meta.continuous.find(c => c.name === mode.name).file;
    data = S.contData.get(file);
    if (!data) {
      data = await fetchBin(file, Uint8Array);
      S.contData.set(file, data);
      if (S.contData.size > 40) S.contData.delete(S.contData.keys().next().value);
    }
  }
  S.mode = mode; S.hidden.clear();
  gl.bindBuffer(gl.ARRAY_BUFFER, buf.val); gl.bufferData(gl.ARRAY_BUFFER, data, gl.DYNAMIC_DRAW);
  if (window.SP) SP.onValue(data);
  updateMask(); updateLUT(); renderLegend(); requestRender();
  $('#colorby').value = mode.type === 'gene' ? '' : mode.type === 'label' ? 'label' : `${mode.type}:${mode.name}`;
  if (mode.type !== 'gene') $('#gene').value = '';
  d3.selectAll('#markers tbody tr').classed('on', d => mode.type === 'gene' && d.j === mode.j);
  renderGeneChart();
}

function renderLegend() {
  const el = d3.select('#legend').html(''), m = S.mode;
  if (isCont()) {
    const r = ramp(), c = m.type === 'gene' ? { label: `${S.genes[m.j]} · log1p CP10K`, lo: 0, hi: S.meta.expr_hi[m.j] }
      : S.meta.continuous.find(x => x.name === m.name);
    el.append('div').text(c.label);
    el.append('div').attr('class', 'ramp').style('background', `linear-gradient(to right, ${d3.range(0, 1.01, 0.1).map(r).join(',')})`);
    el.append('div').attr('class', 'rampax').html(`<span>${d3.format('.3~g')(c.lo)}</span><span>${d3.format('.3~g')(c.hi)}${m.type === 'gene' ? ' (99.9th pct)' : ''}</span>`);
    if (m.type === 'gene') el.append('div').attr('class', 'item').html(`<span class="sw" style="background:${cssVar('--zero')}"></span>not detected`);
    return;
  }
  if (m.type === 'cat' && m.name === 'cluster') return;      // the cluster list is the legend
  const cols = levelColors();
  let items;
  if (m.type === 'label') {
    const lab = i => (S.labels[S.cats.cluster[i]] || {}).label || '';            // a label hides / shows all its clusters
    items = usedLabels().map(l => ({ name: l, color: PAL[usedLabels().indexOf(l) % PAL.length], ids: d3.range(S.nCl).filter(i => lab(i) === l) }));
    items.push({ name: 'unlabelled', color: cssVar('--dim'), ids: d3.range(S.nCl).filter(i => !lab(i)) });
    el.append('div').attr('class', 'muted').text('My labels · click to hide');
  } else {
    items = S.cats[m.name].map((name, i) => ({ name, color: cols[i], i }));
    el.append('div').attr('class', 'muted').text(`${CAT_LABEL[m.name]} · click to hide`);
  }
  const idsOf = d => d.ids || [d.i], isOff = d => idsOf(d).length > 0 && idsOf(d).every(i => S.hidden.has(i));
  el.selectAll('div.item').data(items).join('div').attr('class', d => 'item' + (isOff(d) ? ' off' : ''))
    .html(d => `<span class="sw" style="background:${d.color}"></span>${esc(d.name)}`)
    .on('click', (e, d) => {
      const off = isOff(d);
      idsOf(d).forEach(i => off ? S.hidden.delete(i) : S.hidden.add(i));
      updateLUT(); renderLegend(); requestRender();          // the colour table is shared: UMAP and tissue view hide together
    });
}

// ---------------------------------------------------------------- selection
function select(i, additive) {
  if (i < 0) S.selected.clear();
  else if (additive) S.selected.has(i) ? S.selected.delete(i) : S.selected.add(i);
  else { S.selected.clear(); S.selected.add(i); }
  updateMask(); positionLabels(); renderList(); renderDetail(); renderMarkers(true); renderGeneChart();
  requestRender();
  if (i >= 0) zoomToSelection();
}

// ---------------------------------------------------------------- marker statistics
// Sums over (sample x cluster) groups. stats.f32 = [sum CP10K | n expressing], each nGroups x nGenes.
function computeMarkers() {
  const { nG, nS, nGr, stats } = S, g = S.meta.groups;
  const refv = $('#ref').value, ref = refv === 'rest' ? null : +refv;
  const sIn = new Float64Array(nS * nG), sOut = new Float64Array(nS * nG);
  const zIn = new Float64Array(nG), zOut = new Float64Array(nG);
  const nIn = new Float64Array(nS), nOut = new Float64Array(nS);
  for (let k = 0; k < nGr; k++) {
    const c = g.cluster[k], s = g.sample[k], isIn = S.selected.has(c);
    if (!isIn && ref !== null && c !== ref) continue;
    const dst = isIn ? sIn : sOut, z = isIn ? zIn : zOut, o = k * nG, oz = (nGr + k) * nG, od = s * nG;
    (isIn ? nIn : nOut)[s] += g.n[k];
    for (let j = 0; j < nG; j++) { dst[od + j] += stats[o + j]; z[j] += stats[oz + j]; }
  }
  const NI = d3.sum(nIn), NO = d3.sum(nOut);
  const cores = d3.range(nS).filter(s => nIn[s] >= MIN_REP_CELLS && nOut[s] >= MIN_REP_CELLS);
  const rows = [];
  for (let j = 0; j < nG; j++) {
    let a = 0, b = 0;
    for (let s = 0; s < nS; s++) { a += sIn[s * nG + j]; b += sOut[s * nG + j]; }
    const per = cores.map(s => Math.log2(sIn[s * nG + j] / nIn[s] + 1) - Math.log2(sOut[s * nG + j] / nOut[s] + 1));
    rows.push({
      j, gene: S.genes[j], meanIn: a / NI, meanOut: NO ? b / NO : NaN,
      lfc: NO ? Math.log2(a / NI + 1) - Math.log2(b / NO + 1) : NaN,
      pctIn: zIn[j] / NI, pctOut: NO ? zOut[j] / NO : NaN,
      coreLfc: per.length ? d3.median(per) : NaN, coreUp: per.filter(v => v > 0).length, coreN: per.length,
    });
  }
  S.markers = rows; S.markerN = { NI, NO, cores: cores.length };
}

const plural = w => w.endsWith('s') ? w : w + 's';
const COLS = () => [
  { key: 'gene', label: 'Gene', tip: 'Click a row to colour the map by this gene' },
  { key: 'lfc', label: 'log2FC', tip: 'Pooled cells: log2(mean CP10K + 1) in the selection minus the same in the reference', f: d3.format('+.2f') },
  { key: 'pctIn', label: '% in', tip: 'Share of selected cells with at least one transcript', f: pct },
  { key: 'pctOut', label: '% ref', tip: 'Share of reference cells with at least one transcript. A high value means the gene is not specific (or leaks from neighbouring cells)', f: pct },
  { key: 'meanIn', label: 'mean in', tip: 'Mean CP10K in the selection', f: d3.format('.3~s') },
  { key: 'coreLfc', label: `${S.rep} log2FC`, tip: `Median of the per-${S.rep} log2FC (${plural(S.rep)} with >= ${MIN_REP_CELLS} cells on both sides). The ${S.rep}, not the cell, is the replicate`, f: d3.format('+.2f') },
  { key: 'coreUp', label: `${plural(S.rep)} up`, tip: `${plural(S.rep)} in which the per-${S.rep} log2FC is positive / ${plural(S.rep)} evaluated`, f: (v, d) => `${v}/${d.coreN}` },
];

function renderMarkers(recompute) {
  const sec = $('#markersec');
  sec.hidden = !S.selected.size;
  if (!S.selected.size) return;
  if (recompute) {
    const ref = d3.select('#ref'), cur = $('#ref').value;
    ref.selectAll('option').data(['rest', ...d3.range(S.nCl).filter(i => !S.selected.has(i))]).join('option')
      .attr('value', d => d).text(d => d === 'rest' ? 'all other cells' : clusterName(d));
    $('#ref').value = [...ref.node().options].some(o => o.value === cur) ? cur : 'rest';
    computeMarkers();
  }
  const { NI, NO, cores } = S.markerN;
  $('#markertitle').textContent = `${fmtInt(NI)} cells vs ${fmtInt(NO)} · ${cores} ${plural(S.rep)} evaluated`;
  const minp = (+$('#minpct').value || 0) / 100, q = $('#genefilter').value.trim().toUpperCase();
  const { key, dir } = S.sort;
  const rows = S.markers.filter(d => (q ? d.gene.toUpperCase().includes(q) : d.pctIn >= minp))
    .sort((a, b) => key === 'gene' ? dir * a.gene.localeCompare(b.gene) : dir * ((a[key] - b[key]) || 0) || a.j - b.j);
  const tbl = d3.select('#markers').selectAll('table').data([0]).join(en => {
    const t = en.append('table'); t.append('thead').append('tr'); t.append('tbody'); return t;
  });
  const cols = COLS();
  tbl.select('thead tr').selectAll('th').data(cols).join('th').attr('title', d => d.tip)
    .text(d => d.label + (d.key === key ? (dir < 0 ? ' ▾' : ' ▴') : ''))
    .on('click', (e, d) => { S.sort = { key: d.key, dir: S.sort.key === d.key ? -S.sort.dir : (d.key === 'gene' ? 1 : -1) }; renderMarkers(false); });
  tbl.select('tbody').selectAll('tr').data(rows, d => d.j).join('tr')
    .classed('on', d => S.mode.type === 'gene' && S.mode.j === d.j)
    .on('click', (e, d) => setMode({ type: 'gene', j: d.j }))
    .order()
    .selectAll('td').data(d => cols.map(c => ({ c, d }))).join('td')
    // fewer than half of the samples agree on the sign: the pooled fold change is carried by a few cores
    .classed('weak', ({ c, d }) => c.key === 'coreUp' && d.coreN > 0 && d.coreUp / d.coreN < 0.5)
    .text(({ c, d }) => c.f ? (Number.isNaN(d[c.key]) ? '–' : c.f(d[c.key], d)) : d[c.key]);
}

// mean CP10K and % expressing of one gene in every cluster
function geneByCluster(j) {
  const { nG, nGr, stats } = S, g = S.meta.groups;
  const sum = new Float64Array(S.nCl), nz = new Float64Array(S.nCl), n = new Float64Array(S.nCl);
  for (let k = 0; k < nGr; k++) { const c = g.cluster[k]; sum[c] += stats[k * nG + j]; nz[c] += stats[(nGr + k) * nG + j]; n[c] += g.n[k]; }
  return d3.range(S.nCl).map(i => ({ i, id: S.cats.cluster[i], mean: sum[i] / n[i], pct: nz[i] / n[i] }));
}

function renderGeneChart() {
  const sec = $('#genesec');
  sec.hidden = S.mode.type !== 'gene';
  if (sec.hidden) return;
  const j = S.mode.j, data = geneByCluster(j);
  $('#genetitle').innerHTML = `${esc(S.genes[j])} <span class="muted">mean CP10K per cluster (hover for % expressing)</span>`;
  const W = $('#genechart').clientWidth || 470, H = 150, m = { t: 6, r: 4, b: 30, l: 38 };
  const x = d3.scaleBand().domain(data.map(d => d.id)).range([m.l, W - m.r]).paddingInner(0.2);
  const y = d3.scaleLinear().domain([0, d3.max(data, d => d.mean) || 1]).nice().range([H - m.b, m.t]);
  const s = d3.select('#genechart').selectAll('svg').data([0]).join('svg').attr('viewBox', `0 0 ${W} ${H}`).attr('height', H);
  s.selectAll('*').remove();
  s.append('g').attr('class', 'axis').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(4, '~s').tickSize(-(W - m.l - m.r)))
    .call(g => g.select('.domain').remove());
  s.append('g').attr('class', 'axis').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickSize(0).tickPadding(5))
    .selectAll('text').attr('transform', 'rotate(-60)').style('text-anchor', 'end').attr('dx', '-0.3em').attr('dy', '0.2em')
    .text(d => d === 'isolated' ? 'iso' : d);
  s.append('g').selectAll('path').data(data).join('path')
    .attr('d', d => { const w = x.bandwidth(), h = y(0) - y(d.mean), r = Math.min(3, w / 2, h); const X = x(d.id), Y = y(0);
      return `M${X},${Y}v${-(h - r)}q0,${-r} ${r},${-r}h${w - 2 * r}q${r},0 ${r},${r}v${h - r}z`; })
    .attr('fill', d => S.selected.has(d.i) ? cssVar('--bar-strong') : cssVar('--bar'));
  s.append('g').selectAll('rect').data(data).join('rect')                 // full-height hit targets
    .attr('x', d => x(d.id)).attr('width', x.step()).attr('y', m.t).attr('height', H - m.b - m.t).attr('fill', 'transparent')
    .style('cursor', 'pointer')
    .on('mousemove', (e, d) => tip(e, `<b>${esc(clusterName(d.i))}</b><br>mean ${d3.format('.3~s')(d.mean)} CP10K · ${pct(d.pct)} of cells`))
    .on('mouseleave', () => tip(null))
    .on('click', (e, d) => select(d.i, e.shiftKey));
}

// ---------------------------------------------------------------- sidebar
function clusterName(i) {
  const id = S.cats.cluster[i], l = (S.labels[id] || {}).label;
  return `cluster ${id}` + (l ? ` · ${l}` : '');
}
function topOf(i, v) {
  const row = S.meta.crosstab[v][i], n = d3.sum(row), k = d3.maxIndex(row);
  return { name: S.cats[v][k], frac: n ? row[k] / n : 0 };
}

function renderList() {
  const rows = d3.select('#clusterlist').selectAll('div.crow').data(S.meta.clusters, d => d.id).join(en => {
    const r = en.append('div').attr('class', 'crow');
    r.append('span').attr('class', 'sw');
    r.append('b');
    r.append('input').attr('list', 'labelopts').attr('placeholder', 'label…')
      .on('click', e => e.stopPropagation())
      .on('change', function (e, d) { setLabel([d.id], { label: this.value.trim() }); });
    r.append('span').attr('class', 'num');
    r.append('span').attr('class', 'l1');
    r.on('click', (e, d) => select(S.cats.cluster.indexOf(d.id), e.shiftKey));
    return r;
  });
  rows.classed('sel', (d, i) => S.selected.has(i));
  rows.select('.sw').style('background', (d, i) => clusterColor(i));
  rows.select('b').text(d => d.id);
  rows.select('input').property('value', d => (S.labels[d.id] || {}).label || '');
  rows.select('.num').text(d => fmtInt(d.n));
  const sum0 = S.meta.summary[0];                    // e.g. an automatic annotation: its most frequent level per cluster
  rows.select('.l1').text((d, i) => { if (!sum0) return ''; const t = topOf(i, sum0); return `${t.name} ${pct(t.frac)}`; })
    .attr('title', sum0 ? `Most frequent ${CAT_LABEL[sum0]} in this cluster` : null);
  const nl = S.meta.clusters.filter(d => (S.labels[d.id] || {}).label).length;
  $('#nlabelled').textContent = `${nl}/${S.nCl} labelled`;
  d3.select('#labelopts').selectAll('option').data([...new Set([...usedLabels(), ...S.meta.summary.flatMap(v => S.cats[v])])]).join('option').attr('value', d => d);
}

function compBlock(v, sel, enrich) {
  const tab = S.meta.crosstab[v], lv = S.cats[v], tot = new Float64Array(lv.length), all = new Float64Array(lv.length);
  tab.forEach((row, i) => row.forEach((x, k) => { all[k] += x; if (sel.has(i)) tot[k] += x; }));
  const n = d3.sum(tot), N = d3.sum(all);
  const items = d3.range(lv.length).map(k => ({ name: lv[k], f: tot[k] / n, e: (tot[k] / n) / (all[k] / N) }))
    .filter(d => d.f > 0).sort((a, b) => b.f - a.f).slice(0, 6);
  return `<div><h3>${CAT_LABEL[v]}</h3>` + items.map(d =>
    `<div class="brow" title="${esc(d.name)}: ${pct(d.f)} of the selection; ${d.e.toFixed(2)}x its share of all cells">
       <span class="nm">${esc(d.name)}</span><span class="track"><span class="fill" style="display:block;width:${(d.f * 100).toFixed(1)}%"></span></span>
       <span class="v">${pct(d.f)}${enrich ? ` <span class="${d.e > 2 ? 'warn' : ''}">×${d.e.toFixed(1)}</span>` : ''}</span></div>`).join('') + '</div>';
}

function renderDetail() {
  const el = $('#detail'), sel = S.selected;
  if (!sel.size) { el.innerHTML = '<p class="muted">Click a cluster on the map or in the list. Shift-click adds a cluster to the selection; markers are then computed for the union. Esc clears.</p>'; $('#compsec').hidden = true; return; }
  const ids = [...sel].sort((a, b) => a - b), cls = ids.map(i => S.meta.clusters[i]);
  const n = d3.sum(cls, d => d.n), g = S.meta.groups, perCore = new Float64Array(S.nS);
  for (let k = 0; k < S.nGr; k++) if (sel.has(g.cluster[k])) perCore[g.sample[k]] += g.n[k];
  const nc = perCore.filter(v => v >= MIN_REP_CELLS).length, top = d3.maxIndex(perCore);
  const one = cls.length === 1 ? cls[0] : null;
  const labs = [...new Set(ids.map(i => (S.labels[S.cats.cluster[i]] || {}).label || ''))];
  const notes = [...new Set(ids.map(i => (S.labels[S.cats.cluster[i]] || {}).note || ''))];
  const topShare = perCore[top] / n;
  el.innerHTML = `
    <h2>${ids.length === 1 ? 'Cluster ' + esc(cls[0].id) : 'Clusters ' + cls.map(d => esc(d.id)).join(', ')}</h2>
    <div class="edit">
      <input id="lab" list="labelopts" placeholder="${labs.length > 1 ? 'mixed labels — typing here relabels all selected clusters' : 'cell type label'}" value="${labs.length === 1 ? esc(labs[0]) : ''}">
      <textarea id="note" placeholder="note (evidence, doubts, markers checked)">${notes.length === 1 ? esc(notes[0]) : ''}</textarea>
    </div>
    <div class="kv">
      <div><b>${fmtInt(n)}</b><span>cells (${pct(n / S.N)})</span></div>
      <div><b>${one && one.median_counts != null ? d3.format('.3~r')(one.median_counts) : '–'}</b><span>median counts</span></div>
      <div><b>${one ? d3.format('.3~r')(one.median_genes) : '–'}</b><span>median genes</span></div>
      <div><b>${nc}/${S.nS}</b><span>${plural(S.rep)} with ≥ ${MIN_REP_CELLS} cells</span></div>
    </div>
    <div class="small ${topShare > 0.25 ? 'warn' : 'muted'}">Largest ${esc(S.rep)}: ${esc(S.meta.samples[top].id)} holds ${pct(topShare)} of these cells.</div>
`;
  $('#compsec').hidden = false;
  $('#comp').innerHTML = `
    <div class="comp">
      ${S.meta.composition.map(v => compBlock(v, sel, true)).join('')}
    </div>
    <p class="muted small">×n = share in the selection relative to the share among all cells. A cluster dominated by one ${esc(S.rep)}, slide or batch may be technical rather than a cell type.</p>`;
  $('#lab').addEventListener('change', e => setLabel(ids.map(i => S.cats.cluster[i]), { label: e.target.value.trim() }));
  $('#note').addEventListener('change', e => setLabel(ids.map(i => S.cats.cluster[i]), { note: e.target.value.trim() }));
}

// ---------------------------------------------------------------- labels: state and persistence
const fingerprint = () => S.meta.clusters.map(d => `${d.id}:${d.n}`).join('|');
const LS_KEY = 'slideshow-labels';          // browser fallback when the server is unreachable; checked by fingerprint

function setLabel(ids, patch) {
  for (const id of ids) S.labels[id] = { label: '', note: '', ...S.labels[id], ...patch };
  renderList(); positionLabels();
  if (S.mode.type === 'label') { updateLUT(); renderLegend(); requestRender(); }
  if (S.selected.size) { $('#markertitle') && renderMarkers(true); if ($('#lab') && patch.label !== undefined) $('#lab').value = patch.label; }
  save();
}

function payload() {
  return { dataset: S.meta.title, source: S.meta.source, fingerprint: fingerprint(), saved: new Date().toISOString(),
    rows: S.meta.clusters.map(d => ({ cluster: d.id, label: (S.labels[d.id] || {}).label || '', note: (S.labels[d.id] || {}).note || '', n_cells: d.n })) };
}

let saveTimer;
function save() {
  try { localStorage.setItem(LS_KEY, JSON.stringify(payload())); } catch (e) { /* private window */ }
  clearTimeout(saveTimer);
  saveTimer = setTimeout(async () => {
    const st = $('#savestatus');
    try {
      const r = await fetch('api/annotations', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload()) });
      if (!r.ok) throw new Error(r.status);
      const j = await r.json();
      st.className = 'muted'; st.textContent = `saved ${new Date().toLocaleTimeString()} → ${j.csv}`;
    } catch (e) {
      st.className = 'warn'; st.textContent = 'NOT saved to disk (server unreachable) — use Export CSV';
    }
  }, 400);
}

async function loadLabels() {
  let p = null;
  try { const r = await fetch('api/annotations'); if (r.ok) p = await r.json(); } catch (e) { /* no server */ }
  if (!p || !p.rows) { try { p = JSON.parse(localStorage.getItem(LS_KEY)); } catch (e) { p = null; } }
  if (!p || !p.rows) return;
  if (p.fingerprint !== fingerprint() &&
      !confirm('Saved labels were made on a different clustering (cluster sizes differ).\nCluster numbers may not mean the same thing any more.\n\nApply them anyway?')) return;
  for (const r of p.rows) if (S.cats.cluster.includes(String(r.cluster)) && (r.label || r.note)) S.labels[r.cluster] = { label: r.label || '', note: r.note || '' };
}

function exportCSV() {
  const q = s => `"${String(s).replace(/"/g, '""')}"`;
  const csv = 'cluster,label,note,n_cells\n' + payload().rows.map(r => [q(r.cluster), q(r.label), q(r.note), r.n_cells].join(',')).join('\n') + '\n';
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([csv], { type: 'text/csv' })); a.download = 'cluster_labels.csv'; a.click();
  URL.revokeObjectURL(a.href);
}

// ---------------------------------------------------------------- UI wiring
function tip(e, html) {
  const t = $('#tooltip');
  if (!e) { t.hidden = true; return; }
  t.innerHTML = html; t.hidden = false;
  const w = t.offsetWidth, h = t.offsetHeight;
  t.style.left = Math.min(window.innerWidth - w - 8, e.clientX + 14) + 'px';
  t.style.top = Math.min(window.innerHeight - h - 8, e.clientY + 14) + 'px';
}

function initUI() {
  svg = d3.select('#overlay');
  zoom = d3.zoom().scaleExtent([0.5, 300]).clickDistance(4)
    .on('start', () => { S.interacting = true; })
    .on('zoom', e => { S.t = e.transform; positionLabels(); requestRender(); if (window.SP) SP.hideUmapMark(); })
    .on('end', () => { S.interacting = false; requestRender(); });
  svg.call(zoom).on('dblclick.zoom', null);
  svg.selectAll('text.cl').data(S.meta.clusters).join('text').attr('class', 'cl');
  svg.on('click', e => { if (e.defaultPrevented) return; const [x, y] = d3.pointer(e); select(pick(x, y), e.shiftKey); });
  let hoverRaf = 0;
  svg.on('mousemove', e => {
    if (hoverRaf || S.interacting) return;
    hoverRaf = requestAnimationFrame(() => {
      hoverRaf = 0;
      const [x, y] = d3.pointer(e), c = pick(x, y);
      if (window.SP) SP.markFromUmap(x, y);
      if (c < 0) return tip(null);
      const tops = S.meta.summary.map(v => { const t = topOf(c, v); return `${esc(CAT_LABEL[v])}: ${esc(t.name)} ${pct(t.frac)}`; });
      tip(e, `<b>${esc(clusterName(c))}</b><br>${fmtInt(S.meta.clusters[c].n)} cells` + (tops.length ? `<br><span class="muted">${tops.join('<br>')}</span>` : ''));
    });
  }).on('mouseleave', () => tip(null));

  const cb = d3.select('#colorby');
  cb.append('option').attr('value', '').attr('hidden', true).text('gene');
  cb.append('option').attr('value', 'cat:cluster').text(CAT_LABEL.cluster);
  cb.append('option').attr('value', 'label').text('My labels');
  for (const c of S.meta.categories) if (c.name !== 'cluster') cb.append('option').attr('value', 'cat:' + c.name).text(c.label);
  for (const c of S.meta.continuous) cb.append('option').attr('value', 'cont:' + c.name).text(c.label);
  cb.on('change', function () {
    const [type, name] = this.value.split(':');
    setMode(type === 'label' ? { type } : { type, name });
  });
  d3.select('#genelist').selectAll('option').data(S.genes).join('option').attr('value', d => d);
  $('#gene').addEventListener('change', e => {
    const q = e.target.value.trim().toUpperCase(), j = S.genes.findIndex(g => g.toUpperCase() === q);
    if (j >= 0) setMode({ type: 'gene', j }); else if (q) e.target.value = '';
  });
  d3.select('#embtoggle').selectAll('button').data(S.emb).join('button').text(d => d.short || d.name).attr('title', d => d.label)
    .classed('on', (d, i) => i === 0).on('click', (e, d) => setEmbedding(S.emb.indexOf(d)));
  if (S.emb.length < 2) $('#embtoggle').hidden = true;
  for (const id of ['#psize', '#palpha']) $(id).addEventListener('input', requestRender);
  $('#showlabels').addEventListener('change', positionLabels);
  $('#filtervar').addEventListener('change', e => { if (e.target.value) addFilter(e.target.value); });
  renderFilters();
  $('#reset').addEventListener('click', () => svg.transition().duration(500).call(zoom.transform, d3.zoomIdentity));
  $('#export').addEventListener('click', exportCSV);
  $('#ref').addEventListener('change', () => { computeMarkers(); renderMarkers(false); });
  $('#minpct').addEventListener('input', () => renderMarkers(false));
  $('#genefilter').addEventListener('input', () => renderMarkers(false));
  window.addEventListener('keydown', e => { if (e.key === 'Escape' && !/INPUT|TEXTAREA/.test(document.activeElement.tagName)) { if (S.cellSel) setCellSel(null); select(-1); } });
  new ResizeObserver(resize).observe($('#map'));
  new ResizeObserver(() => requestAnimationFrame(renderGeneChart)).observe($('#genechart'));    // side panel resized: redraw at its width
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { updateLUT(); renderLegend(); renderList(); renderGeneChart(); requestRender(); });
  renderList(); positionLabels();
}

init();
