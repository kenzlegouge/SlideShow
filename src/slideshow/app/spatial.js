'use strict';
// Tissue view: the cells of one slide in their coordinates (um, y down).
// Slide scale: one point per cell (centroid). Zoomed in: the segmentation polygons of the tiles in view, loaded on
// demand from spatial/tiles/<k>.bin and triangulated with earcut. Colours come from the same lookup textures as the
// UMAP (cluster, labels, gene, ...), so every colouring works here unchanged. Uses the globals of app.js (S, tip, ...).
// Polygons are whatever segmentation the bundle was built with (for Xenium: cell territories, not membranes).

window.SP = (() => {
  const POLY_PPU = 0.3;        // px per um above which polygons replace points (a 10 um cell is then 3 px wide)
  const OUTLINE_PPU = 0.8;     // px per um above which cell outlines are drawn
  const POLY_BUDGET = 60000, CACHE_BUDGET = 200000, FAST = 250000, BUCKET_UM = 16;   // polygons drawn at once / kept in memory
  const P = { ready: false, slide: null, t: d3.zoomIdentity, W: 1, H: 1, cores: new Map(), loading: new Set(),
    points: {}, grids: {}, val: null, interacting: false, lasso: false, umapGrid: {} };
  let gl, prog, U = {}, LOC = {}, tex = {}, svg, zoom, queued = false;

  const VS = `#version 300 es
in vec2 aPos; in float aCl; in float aVal; in float aSel; in float aKeep;
uniform vec4 uXf, uDim, uFlat; uniform float uSize; uniform int uPass, uCont, uUseSel;
uniform sampler2D uLUT, uMask;
out vec4 vCol;
void main() {
  gl_Position = vec4(aPos.x * uXf.x + uXf.y, aPos.y * uXf.z + uXf.w, 0.0, 1.0);
  gl_PointSize = uSize;
  if (aKeep < 0.5) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); return; }          // cell filter (app.js)
  if (uFlat.a > 0.0) { vCol = uFlat; return; }
  vec4 col = texelFetch(uLUT, ivec2(int(aVal + 0.5), 0), 0);
  float sel = texelFetch(uMask, ivec2(int(aCl + 0.5), 0), 0).r;
  if (uUseSel == 1) sel *= aSel;
  bool fg = sel > 0.5 && (uCont == 0 || aVal > 0.5);
  bool drop = col.a < 0.01 || (uPass == 0 && fg) || (uPass == 1 && !fg);
  vCol = sel > 0.5 ? col : uDim;
  if (drop) gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
}`;
  const FS = `#version 300 es
precision mediump float;
in vec4 vCol; uniform float uAlpha; uniform int uPoint; out vec4 o;
void main() {
  if (uPoint == 1) { vec2 d = gl_PointCoord - 0.5; if (dot(d, d) > 0.25) discard; }
  o = vec4(vCol.rgb, vCol.a * uAlpha);
}`;

  // ---------------------------------------------------------------- set-up
  async function init() {
    const pane = $('#tissue');
    if (!S.meta.spatial) { pane.remove(); applyLayout(); return; }   // bundle without a tissue view: the drawer says so (layout.js)
    const [xy, index] = await Promise.all([fetchBin('spatial_xy.f32', Float32Array),
      fetch(BUNDLE + S.meta.spatial.index).then(r => { if (!r.ok) throw new Error('spatial/index.json: HTTP ' + r.status); return r.json(); })]);
    if (xy.length !== 2 * S.N) throw new Error('spatial_xy.f32 does not match the bundle');
    P.slideCat = S.meta.spatial.slide_category; P.sampleCat = S.meta.sample.name;
    const [slideCode, sampleCode] = await Promise.all([catCodes(P.slideCat), catCodes(P.sampleCat)]);
    P.slides = S.cats[P.slideCat]; P.samples = S.cats[P.sampleCat];
    P.xy = xy; P.index = index; P.slideCode = slideCode; P.sampleCode = sampleCode;
    // bundle cells of each slide, in bundle (random) order
    const nSl = P.slides.length, cnt = new Uint32Array(nSl);
    for (let i = 0; i < S.N; i++) cnt[slideCode[i]]++;
    P.slideCells = P.slides.map((s, k) => new Uint32Array(cnt[k]));
    cnt.fill(0);
    for (let i = 0; i < S.N; i++) { const k = slideCode[i]; P.slideCells[k][cnt[k]++] = i; }
    P.samplesOf = {}; P.tilesOf = {};
    for (const c of index.samples) (P.samplesOf[c.slide] = P.samplesOf[c.slide] || []).push(c);
    index.tiles.forEach((t, k) => { t.id = k; (P.tilesOf[t.slide] = P.tilesOf[t.slide] || []).push(t); });
    P.hasPoly = index.tiles.length > 0;
    if (!P.hasPoly) { $('#spout').parentElement.hidden = true; $('#spnuc').parentElement.hidden = true; }

    LAYOUT.tissueReady = true; applyLayout();        // shown unless switched off in the panel drawer
    initGL();
    svg = d3.select('#spoverlay');
    zoom = d3.zoom().scaleExtent([0.5, 600]).clickDistance(4)
      .on('start', () => { P.interacting = true; })
      .on('zoom', e => { P.t = e.transform; overlay(); request(); })
      .on('end', () => { P.interacting = false; request(); });
    svg.call(zoom).on('dblclick.zoom', null);
    svg.append('g').attr('class', 'corelabels');
    svg.append('path').attr('class', 'lasso');
    svg.append('circle').attr('class', 'cellmark').attr('r', 9).attr('display', 'none');
    svg.append('g').attr('class', 'scalebar');
    d3.select('#overlay').append('circle').attr('class', 'cellmark').attr('r', 7).attr('display', 'none');
    svg.on('click', e => { if (e.defaultPrevented || P.lasso) return; const [x, y] = d3.pointer(e); const i = nearest(x, y); select(i < 0 ? -1 : S.catData.cluster[i], e.shiftKey); });
    let raf = 0;
    svg.on('mousemove', e => {
      if (raf || P.interacting || P.lasso) return;
      raf = requestAnimationFrame(() => {
        raf = 0;
        const [x, y] = d3.pointer(e), i = nearest(x, y);
        markUmap(i);
        if (i < 0) return tip(null);
        tip(e, `<b>${esc(clusterName(S.catData.cluster[i]))}</b><br>${esc(S.rep)} ${esc(P.samples[sampleCode[i]])}`);
      });
    }).on('mouseleave', () => { tip(null); markUmap(-1); });
    lassoEvents();

    const slides = P.slides.filter((s, k) => P.slideCells[k].length && index.slides[s]);
    if (slides.length <= 12) {                            // one button per slide, or a menu when there are many
      d3.select('#spslides').selectAll('button').data(slides).join('button').text(d => d).attr('title', d => d).on('click', (e, d) => setSlide(d));
    } else {
      d3.select('#spslides').append('select').attr('id', 'spslidesel').on('change', function () { setSlide(this.value); })
        .selectAll('option').data(slides).join('option').attr('value', d => d).text(d => d);
    }
    $('#spnuc').addEventListener('change', request);
    $('#spout').addEventListener('change', request);
    $('#spfit').addEventListener('click', () => svg.transition().duration(400).call(zoom.transform, d3.zoomIdentity));
    $('#splasso').addEventListener('click', () => setLasso(!P.lasso));
    $('#spclear').addEventListener('click', () => { setCellSel(null); info(); });
    new ResizeObserver(resize).observe(pane);          // also fires when the pane is shown, hidden or dragged wider (layout.js)
    P.ready = true;
    setSlide(slides[0]);
  }

  function initGL() {
    gl = $('#glsp').getContext('webgl2', { alpha: false, antialias: true });
    const sh = (type, src) => { const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s)); return s; };
    prog = gl.createProgram();
    gl.attachShader(prog, sh(gl.VERTEX_SHADER, VS)); gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, FS));
    ['aPos', 'aCl', 'aVal', 'aSel', 'aKeep'].forEach((a, i) => { gl.bindAttribLocation(prog, i, a); LOC[a] = i; });
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog));
    gl.useProgram(prog);
    for (const u of ['uXf', 'uDim', 'uFlat', 'uSize', 'uPass', 'uCont', 'uUseSel', 'uLUT', 'uMask', 'uAlpha', 'uPoint']) U[u] = gl.getUniformLocation(prog, u);
    const mk = unit => { const t = gl.createTexture(); gl.activeTexture(gl.TEXTURE0 + unit); gl.bindTexture(gl.TEXTURE_2D, t);
      for (const p of [gl.TEXTURE_MIN_FILTER, gl.TEXTURE_MAG_FILTER]) gl.texParameteri(gl.TEXTURE_2D, p, gl.NEAREST);
      for (const p of [gl.TEXTURE_WRAP_S, gl.TEXTURE_WRAP_T]) gl.texParameteri(gl.TEXTURE_2D, p, gl.CLAMP_TO_EDGE);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 256, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array(1024).fill(255)); return t; };
    tex.lut = mk(0); tex.mask = mk(1);
    gl.uniform1i(U.uLUT, 0); gl.uniform1i(U.uMask, 1);
    gl.enable(gl.BLEND); gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    gl.disable(gl.DEPTH_TEST);
  }
  const setTex = (unit, t, a) => { gl.activeTexture(gl.TEXTURE0 + unit); gl.bindTexture(gl.TEXTURE_2D, t);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 256, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, a); };

  // a vertex array with its buffers: specs = [attribute, data, components, GL type, usage]
  function geometry(specs, count, mode) {
    const vao = gl.createVertexArray(), bufs = {};
    gl.bindVertexArray(vao);
    for (const [name, data, size, type, usage] of specs) {
      const b = bufs[name] = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.bufferData(gl.ARRAY_BUFFER, data, usage || gl.STATIC_DRAW);
      gl.enableVertexAttribArray(LOC[name]); gl.vertexAttribPointer(LOC[name], size, type, false, 0, 0);
    }
    gl.bindVertexArray(null);
    return { vao, bufs, count, mode };
  }
  function setBuf(g, name, data) { gl.bindBuffer(gl.ARRAY_BUFFER, g.bufs[name]); gl.bufferData(gl.ARRAY_BUFFER, data, gl.DYNAMIC_DRAW); }
  function free(g) { if (!g) return; Object.values(g.bufs).forEach(b => gl.deleteBuffer(b)); gl.deleteVertexArray(g.vao); }
  function gather(idx, src) { const o = new Uint8Array(idx.length); if (src) for (let k = 0; k < idx.length; k++) o[k] = idx[k] >= 0 ? src[idx[k]] : 1; else o.fill(1); return o; }   // idx < 0: not a bundle cell, always drawn

  // ---------------------------------------------------------------- slide and cores
  function setSlide(name) {
    P.slide = name; P.k = P.slides.indexOf(name);
    d3.selectAll('#spslides button').classed('on', d => d === name);
    if ($('#spslidesel')) $('#spslidesel').value = name;
    P.bbox = P.index.slides[name].bbox;
    if (!P.points[name]) {
      const idx = P.slideCells[P.k], pos = new Float32Array(2 * idx.length);
      for (let k = 0; k < idx.length; k++) { pos[2 * k] = P.xy[2 * idx[k]]; pos[2 * k + 1] = P.xy[2 * idx[k] + 1]; }
      const g = geometry([['aPos', pos, 2, gl.FLOAT], ['aCl', gather(idx, S.catData.cluster), 1, gl.UNSIGNED_BYTE],
        ['aVal', gather(idx, P.val || S.catData.cluster), 1, gl.UNSIGNED_BYTE, gl.DYNAMIC_DRAW],
        ['aSel', gather(idx, S.cellSel), 1, gl.UNSIGNED_BYTE, gl.DYNAMIC_DRAW], ['aKeep', gather(idx, S.keep), 1, gl.UNSIGNED_BYTE, gl.DYNAMIC_DRAW]], idx.length, gl.POINTS);
      g.idx = idx; P.points[name] = g;
    }
    // sample names on the slide (only when a slide holds several samples, e.g. the cores of a TMA)
    const smp = (P.samplesOf[name] || []).length > 1 ? P.samplesOf[name] : [];
    const short = id => (id.startsWith(name) && /^[_\-. ]/.test(id.slice(name.length)) ? id.slice(name.length + 1) : id);
    const lab = svg.select('g.corelabels').selectAll('text').data(smp, d => d.id).join('text')
      .text(d => short(d.id)).on('click', (e, d) => { e.stopPropagation(); zoomToCore(d); });
    lab.selectAll('title').data(d => [d]).join('title').text(d => `${d.id}: ${fmtInt(d.n)} cells. Click to open.`);
    svg.call(zoom.transform, d3.zoomIdentity);
    resize();
  }

  const TOP = 40;                 // px kept free for the toolbar when the slide is fitted
  const midY = () => TOP + (P.H - TOP) / 2;
  const k0 = () => Math.min(P.W / (P.bbox[2] - P.bbox[0]), (P.H - TOP) / (P.bbox[3] - P.bbox[1])) * 0.94;
  const ppu = () => k0() * P.t.k;
  function toScreen(x, y) { const k = k0(); return P.t.apply([(x - (P.bbox[0] + P.bbox[2]) / 2) * k + P.W / 2, (y - (P.bbox[1] + P.bbox[3]) / 2) * k + midY()]); }
  function toUm(sx, sy) { const k = k0(), [a, b] = P.t.invert([sx, sy]); return [(a - P.W / 2) / k + (P.bbox[0] + P.bbox[2]) / 2, (b - midY()) / k + (P.bbox[1] + P.bbox[3]) / 2]; }
  function zoomToCore(c, instant) {
    const k = k0(), cx = (P.bbox[0] + P.bbox[2]) / 2, cy = (P.bbox[1] + P.bbox[3]) / 2;
    const ax = (c.bbox[0] - cx) * k + P.W / 2, bx = (c.bbox[2] - cx) * k + P.W / 2, ay = (c.bbox[1] - cy) * k + midY(), by = (c.bbox[3] - cy) * k + midY();
    const s = 0.92 * Math.min(P.W / (bx - ax), (P.H - TOP) / (by - ay));
    const tr = d3.zoomIdentity.translate(P.W / 2, midY()).scale(s).translate(-(ax + bx) / 2, -(ay + by) / 2);
    if (instant) svg.call(zoom.transform, tr); else svg.transition().duration(600).call(zoom.transform, tr);
  }
  function visibleTiles() {
    const [x0, y0] = toUm(0, 0), [x1, y1] = toUm(P.W, P.H);
    return (P.tilesOf[P.slide] || []).filter(c => c.bbox[0] < x1 && c.bbox[2] > x0 && c.bbox[1] < y1 && c.bbox[3] > y0);
  }
  const nPoly = tiles => d3.sum(tiles, t => t.n);

  async function loadTile(c) {
    if (P.cores.has(c.id) || P.loading.has(c.id)) return;
    P.loading.add(c.id); info();
    try {
      const r = await fetch(BUNDLE + 'spatial/' + c.file);
      if (!r.ok) throw new Error(r.status);
      const buf = await r.arrayBuffer(), dv = new DataView(buf);
      const n = dv.getInt32(0, true), x0 = dv.getFloat32(4, true), y0 = dv.getFloat32(8, true), q = dv.getFloat32(12, true);
      const own = new Int32Array(buf, 16, n), nvC = new Uint8Array(buf, 16 + 4 * n, n), nvN = new Uint8Array(buf, 16 + 5 * n, n);   // own: bundle index, -1 = context
      let off = 16 + 6 * n; off += off % 2;
      let sC = 0, sN = 0, nOwnV = 0, triMax = 0;
      for (let k = 0; k < n; k++) { sC += nvC[k]; sN += nvN[k]; if (own[k] >= 0) { nOwnV += nvC[k]; triMax += Math.max(0, nvC[k] - 2) * 3; } }
      const vc = new Uint16Array(buf, off, 2 * sC), vn = new Uint16Array(buf, off + 4 * sC, 2 * sN);
      const tri = new Float32Array(2 * triMax), tIdx = new Uint32Array(triMax);
      const edge = new Float32Array(4 * nOwnV), ctx = new Float32Array(4 * (sC - nOwnV)), nuc = new Float32Array(4 * sN);
      const edgeCell = new Int32Array(2 * nOwnV), nucCell = new Int32Array(2 * sN);       // bundle index of the cell each outline vertex belongs to (-1: none)
      let ne = 0, nn = 0;
      let t = 0, e = 0, x = 0, u = 0, pc = 0, pn = 0;
      const flat = [];
      for (let k = 0; k < n; k++) {
        const m = nvC[k], b = own[k];
        flat.length = 2 * m;
        for (let v = 0; v < 2 * m; v += 2) { flat[v] = x0 + q * vc[2 * pc + v]; flat[v + 1] = y0 + q * vc[2 * pc + v + 1]; }
        if (b >= 0) for (let v = 0; v < 2 * m; v++) edgeCell[ne++] = b;
        const dst = b >= 0 ? edge : ctx; let w = b >= 0 ? e : x;
        for (let v = 0; v < m; v++) { const a = 2 * v, z = 2 * ((v + 1) % m); dst[w++] = flat[a]; dst[w++] = flat[a + 1]; dst[w++] = flat[z]; dst[w++] = flat[z + 1]; }
        if (b >= 0) {
          e = w;
          const ix = earcut(flat);
          for (let j = 0; j < ix.length; j++) { tri[2 * t] = flat[2 * ix[j]]; tri[2 * t + 1] = flat[2 * ix[j] + 1]; tIdx[t++] = b; }
        } else x = w;
        pc += m;
        const mn = nvN[k];
        for (let v = 0; v < 2 * mn; v++) nucCell[nn++] = b;
        for (let v = 0; v < mn; v++) { const a = 2 * (pn + v), z = 2 * (pn + (v + 1) % mn);
          nuc[u++] = x0 + q * vn[a]; nuc[u++] = y0 + q * vn[a + 1]; nuc[u++] = x0 + q * vn[z]; nuc[u++] = y0 + q * vn[z + 1]; }
        pn += mn;
      }
      const idx = tIdx.subarray(0, t);
      const fill = geometry([['aPos', tri.subarray(0, 2 * t), 2, gl.FLOAT], ['aCl', gather(idx, S.catData.cluster), 1, gl.UNSIGNED_BYTE],
        ['aVal', gather(idx, P.val || S.catData.cluster), 1, gl.UNSIGNED_BYTE, gl.DYNAMIC_DRAW],
        ['aSel', gather(idx, S.cellSel), 1, gl.UNSIGNED_BYTE, gl.DYNAMIC_DRAW], ['aKeep', gather(idx, S.keep), 1, gl.UNSIGNED_BYTE, gl.DYNAMIC_DRAW]], t, gl.TRIANGLES);
      fill.idx = idx;
      const lineGeom = (v, n, cell) => { const g = geometry([['aPos', v.subarray(0, 2 * n), 2, gl.FLOAT], ['aKeep', gather(cell, S.keep), 1, gl.UNSIGNED_BYTE, gl.DYNAMIC_DRAW]], n, gl.LINES); g.idx = cell; return g; };
      P.cores.set(c.id, { fill, edge: lineGeom(edge, e / 2, edgeCell), ctx: lineGeom(ctx, x / 2, new Int32Array(x / 2).fill(-1)),
        nuc: lineGeom(nuc, u / 2, nucCell), n });
      while (nPoly([...P.cores.keys()].map(id => P.index.tiles[id])) > CACHE_BUDGET && P.cores.size > 1) {   // oldest tile out
        const [id, g] = P.cores.entries().next().value;
        if (id === c.id) break;
        ['fill', 'edge', 'ctx', 'nuc'].forEach(k => free(g[k])); P.cores.delete(id);
      }
    } catch (err) { console.error('tile', c.file, err); }
    P.loading.delete(c.id); info(); request();
  }

  // ---------------------------------------------------------------- drawing
  function request() { if (P.ready && !queued) { queued = true; requestAnimationFrame(draw); } }
  function resize() {
    const el = $('#tissue'); if (!P.ready || el.hidden) return;
    const dpr = window.devicePixelRatio || 1, W = Math.max(1, el.clientWidth), H = Math.max(1, el.clientHeight);
    // keep the same tissue point in the middle and the same magnification when the pane changes size
    const keep = P.slide && P.W > 1 && (W !== P.W || H !== P.H) ? { c: toUm(P.W / 2, midY()), k: ppu() } : null;
    P.W = W; P.H = H;
    gl.canvas.width = Math.round(P.W * dpr); gl.canvas.height = Math.round(P.H * dpr);
    if (keep) {
      const kb = k0(), tk = keep.k / kb, x = (keep.c[0] - (P.bbox[0] + P.bbox[2]) / 2) * kb + P.W / 2, y = (keep.c[1] - (P.bbox[1] + P.bbox[3]) / 2) * kb + midY();
      svg.call(zoom.transform, d3.zoomIdentity.translate(P.W / 2 - x * tk, midY() - y * tk).scale(tk));
    }
    overlay(); request();
  }
  function passes(g, alpha, point) {
    gl.bindVertexArray(g.vao);
    gl.uniform1i(U.uPoint, point ? 1 : 0); gl.uniform4f(U.uFlat, 0, 0, 0, 0);
    const n = point && P.interacting ? Math.min(g.count, FAST) : g.count;
    if (isCont() || S.selected.size || S.cellSel) { gl.uniform1i(U.uPass, 0); gl.uniform1f(U.uAlpha, alpha * 0.55); gl.drawArrays(g.mode, 0, n); }
    gl.uniform1i(U.uPass, 1); gl.uniform1f(U.uAlpha, alpha); gl.drawArrays(g.mode, 0, n);
  }
  function lines(g, c, a) {
    if (!g.count) return;
    gl.bindVertexArray(g.vao); gl.uniform1i(U.uPoint, 0); gl.uniform1f(U.uAlpha, 1);
    gl.uniform4f(U.uFlat, c[0] / 255, c[1] / 255, c[2] / 255, a); gl.drawArrays(gl.LINES, 0, g.count);
  }
  function draw() {
    queued = false;
    if ($('#tissue').hidden || !P.slide) return;
    const bg = rgb(cssVar('--surface')), dim = rgb(cssVar('--dim')), ink = rgb(cssVar('--text')), k = ppu(), dpr = window.devicePixelRatio || 1;
    gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
    gl.clearColor(bg[0] / 255, bg[1] / 255, bg[2] / 255, 1); gl.clear(gl.COLOR_BUFFER_BIT);
    const cx = (P.bbox[0] + P.bbox[2]) / 2, cy = (P.bbox[1] + P.bbox[3]) / 2, kb = k0(), t = P.t;
    // screen x = ((x - cx) kb + W/2) t.k + t.x ; clip = 2 sx / W - 1 ; y is not flipped (image convention)
    gl.uniform4f(U.uXf, 2 * kb * t.k / P.W, 2 * ((P.W / 2 - cx * kb) * t.k + t.x) / P.W - 1,
      -2 * kb * t.k / P.H, 1 - 2 * ((midY() - cy * kb) * t.k + t.y) / P.H);
    gl.uniform4f(U.uDim, dim[0] / 255, dim[1] / 255, dim[2] / 255, 0.6);
    gl.uniform1i(U.uCont, isCont() ? 1 : 0); gl.uniform1i(U.uUseSel, S.cellSel ? 1 : 0);
    const vis = k >= POLY_PPU ? visibleTiles() : [];
    const poly = vis.length > 0 && nPoly(vis) <= POLY_BUDGET;
    let complete = poly;
    if (poly) for (const c of vis) if (!P.cores.has(c.id)) { complete = false; loadTile(c); }
    if (!complete) {
      gl.uniform1f(U.uSize, Math.max(1.3, Math.min(26, k * 7)) * dpr);
      passes(P.points[P.slide], 0.85, true);
    }
    if (poly) for (const c of vis) {
      const g = P.cores.get(c.id); if (!g) continue;
      lines(g.ctx, dim, 0.9);                                         // cells of the core that are not in this bundle
      passes(g.fill, 0.92, false);
      if ($('#spout').checked && k >= OUTLINE_PPU) lines(g.edge, bg, 0.9);
      if ($('#spnuc').checked) lines(g.nuc, ink, 0.55);
    }
    gl.bindVertexArray(null);
    P.polyShown = poly && complete;
  }

  function overlay() {
    if (!P.slide) return;
    const k = ppu();
    svg.select('g.corelabels').attr('display', k < POLY_PPU * 2.5 ? null : 'none').selectAll('text')
      .attr('transform', d => { const p = toScreen((d.bbox[0] + d.bbox[2]) / 2, d.bbox[1]); return `translate(${p[0]},${p[1] - 4})`; });
    const um = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000].reduce((a, b) => (b * k <= 140 ? b : a), 1), w = um * k;
    const sb = svg.select('g.scalebar').attr('transform', `translate(${P.W - w - 16},${P.H - 18})`);
    sb.selectAll('line').data([0]).join('line').attr('x1', 0).attr('x2', w).attr('y1', 0).attr('y2', 0);
    sb.selectAll('text').data([0]).join('text').attr('x', w / 2).attr('y', -5).text(um >= 1000 ? `${um / 1000} mm` : `${um} µm`);
    svg.select('circle.cellmark').attr('display', 'none');
    info();
  }
  function info(extra) {
    if (!P.slide) return;
    const k = ppu(), vis = k >= POLY_PPU ? visibleTiles() : [];
    let s = `${P.slide} · ${fmtInt(P.slideCells[P.k].length)} cells`;
    if (!P.hasPoly) s += ' · no segmentation in this bundle (centroids only)';
    else if (P.loading.size) s += ' · loading segmentation…';
    else if (k < POLY_PPU) s += ' · zoom in to see the segmentation';
    else if (nPoly(vis) > POLY_BUDGET) s += ' · zoom in further for the segmentation';
    if (extra !== undefined) P.extra = extra;
    $('#spinfo').innerHTML = esc(s) + (P.extra ? `<br>${P.extra}` : '');
  }

  // ---------------------------------------------------------------- picking and cross-marking
  // bucket grid over the cells of the current slide (um) or over a UMAP embedding
  function grid(key, n, getX, getY, x0, y0, bucket, nx, ny) {
    if (P.grids[key]) return P.grids[key];
    const start = new Uint32Array(nx * ny + 1), cellB = new Uint32Array(n);
    for (let k = 0; k < n; k++) { const gx = Math.min(nx - 1, Math.max(0, Math.floor((getX(k) - x0) / bucket))), gy = Math.min(ny - 1, Math.max(0, Math.floor((getY(k) - y0) / bucket))); cellB[k] = gy * nx + gx; start[cellB[k] + 1]++; }
    for (let b = 0; b < nx * ny; b++) start[b + 1] += start[b];
    const fill = start.slice(0, nx * ny), order = new Uint32Array(n);
    for (let k = 0; k < n; k++) order[fill[cellB[k]]++] = k;
    return (P.grids[key] = { start, order, x0, y0, bucket, nx, ny });
  }
  function query(g, x, y, r, getX, getY) {
    const gx = Math.floor((x - g.x0) / g.bucket), gy = Math.floor((y - g.y0) / g.bucket), R = Math.ceil(r / g.bucket);
    let best = -1, bd = r * r;
    for (let j = gy - R; j <= gy + R; j++) for (let i = gx - R; i <= gx + R; i++) {
      if (i < 0 || j < 0 || i >= g.nx || j >= g.ny) continue;
      for (let p = g.start[j * g.nx + i]; p < g.start[j * g.nx + i + 1]; p++) {
        const k = g.order[p], dx = getX(k) - x, dy = getY(k) - y, d = dx * dx + dy * dy;
        if (d < bd) { bd = d; best = k; }
      }
    }
    return best;
  }
  function slideGrid() {
    const idx = P.slideCells[P.k], b = P.bbox;
    return grid('s:' + P.slide, idx.length, k => P.xy[2 * idx[k]], k => P.xy[2 * idx[k] + 1], b[0], b[1], BUCKET_UM,
      Math.ceil((b[2] - b[0]) / BUCKET_UM) + 1, Math.ceil((b[3] - b[1]) / BUCKET_UM) + 1);
  }
  function nearest(sx, sy) {                       // bundle index of the cell under the pointer in the tissue view
    const idx = P.slideCells[P.k], [x, y] = toUm(sx, sy);
    const k = query(slideGrid(), x, y, Math.max(8, 7 / ppu()), j => P.xy[2 * idx[j]], j => P.xy[2 * idx[j] + 1]);
    return k < 0 || (S.keep && !S.keep[idx[k]]) ? -1 : idx[k];            // a filtered-out cell cannot be picked
  }
  function markUmap(i) {                           // ring on the UMAP at the position of bundle cell i
    const m = d3.select('#overlay').select('circle.cellmark');
    if (i < 0) return m.attr('display', 'none');
    const d = S.emb[S.embIdx].data, p = S.t.apply(toScreen0(d[2 * i], d[2 * i + 1]));
    m.attr('display', null).attr('cx', p[0]).attr('cy', p[1]).raise();
  }
  function markFromUmap(sx, sy) {                  // called by the UMAP on hover: ring on the tissue if the cell is on this slide
    if (!P.ready || $('#tissue').hidden) return;
    const e = S.emb[S.embIdx], d = e.data, [x, y] = toData(sx, sy), res = 512, ext = 1.3, bucket = 2 * ext / res;
    const g = grid('u:' + e.name, S.N, k => d[2 * k], k => d[2 * k + 1], -ext, -ext, bucket, res, res);
    const i = query(g, x, y, 6 / (Math.min(S.W, S.H) / 2 * 0.95 * S.t.k), k => d[2 * k], k => d[2 * k + 1]);
    const m = svg.select('circle.cellmark');
    if (i < 0 || P.slideCode[i] !== P.k) { m.attr('display', 'none'); if (i >= 0) info(`hovered cell: slide ${esc(P.slides[P.slideCode[i]])}, ${esc(S.rep)} ${esc(P.samples[P.sampleCode[i]])}`); return; }
    const p = toScreen(P.xy[2 * i], P.xy[2 * i + 1]);
    m.attr('display', null).attr('cx', p[0]).attr('cy', p[1]);
    info(`hovered cell: ${esc(S.rep)} ${esc(P.samples[P.sampleCode[i]])}`);
  }

  // ---------------------------------------------------------------- lasso
  function setLasso(on) {
    P.lasso = on; $('#splasso').classList.toggle('on', on); svg.style('cursor', on ? 'crosshair' : null);
    if (on) svg.on('.zoom', null); else { svg.call(zoom).on('dblclick.zoom', null); svg.select('path.lasso').attr('d', null); }
  }
  function lassoEvents() {
    let pts = null;
    svg.on('pointerdown.lasso', e => { if (!P.lasso) return; pts = [d3.pointer(e)]; svg.node().setPointerCapture(e.pointerId); })
      .on('pointermove.lasso', e => { if (!pts) return; pts.push(d3.pointer(e)); svg.select('path.lasso').attr('d', 'M' + pts.map(p => p.join(',')).join('L') + 'Z'); })
      .on('pointerup.lasso', () => {
        if (!pts) return;
        const poly = pts.map(p => toUm(p[0], p[1])); pts = null;
        if (poly.length > 2) {
          const idx = P.slideCells[P.k], mask = new Uint8Array(S.N), byCl = new Map();
          const xs = poly.map(p => p[0]), ys = poly.map(p => p[1]), x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
          let n = 0;
          for (let k = 0; k < idx.length; k++) {
            const i = idx[k], x = P.xy[2 * i], y = P.xy[2 * i + 1];
            if (S.keep && !S.keep[i]) continue;                                   // filtered-out cells are not selectable
            if (x < x0 || x > x1 || y < y0 || y > y1 || !d3.polygonContains(poly, [x, y])) continue;
            mask[i] = 1; n++; const c = S.catData.cluster[i]; byCl.set(c, (byCl.get(c) || 0) + 1);
          }
          if (n) {
            setCellSel(mask);
            const top = [...byCl].sort((a, b) => b[1] - a[1]).slice(0, 6).map(([c, v]) => `${esc(S.cats.cluster[c])}: ${pct(v / n)}`).join(' · ');
            info(`lasso: <b>${fmtInt(n)}</b> cells · clusters ${top}`);
          }
        }
        setLasso(false);
      });
  }

  // ---------------------------------------------------------------- hooks called by app.js
  function eachGeom(f) { for (const g of Object.values(P.points)) f(g); for (const c of P.cores.values()) f(c.fill); }
  function eachKeepGeom(f) { eachGeom(f); for (const c of P.cores.values()) { f(c.edge); f(c.nuc); } }
  return {
    init, markFromUmap,
    openSample(id) { const c = P.ready && P.index.samples.find(d => d.id === id); if (c) { setSlide(c.slide); zoomToCore(c, true); } },
    onLUT(a) { if (P.ready) { setTex(0, tex.lut, a); request(); } else P.pendingLUT = a; },
    onMask(a) { if (P.ready) { setTex(1, tex.mask, a); request(); } },
    onValue(data) { P.val = data; if (!P.ready) return; eachGeom(g => setBuf(g, 'aVal', gather(g.idx, data))); request(); },
    onCellSel() { if (!P.ready) return; eachGeom(g => setBuf(g, 'aSel', gather(g.idx, S.cellSel))); if (!S.cellSel) info(''); request(); },
    onKeep() { if (!P.ready) return; eachKeepGeom(g => setBuf(g, 'aKeep', gather(g.idx, S.keep))); request(); },
    onRender() { request(); },
    hideUmapMark() { d3.select('#overlay').select('circle.cellmark').attr('display', 'none'); },
  };
})();
