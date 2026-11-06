'use strict';
// Panel layout of the dashboard: which of the three panels are shown (embedding map, tissue view, side panel) and how
// wide they are. The ☰ button opens a drawer on the left with one switch per panel; the bars between the panels are
// dragged to resize them (double-click a bar: default width). The layout is kept in this browser (localStorage). Each
// dashboard runs on its own port, i.e. its own origin, so each dashboard remembers its own layout.
// Loaded before spatial.js and app.js: at load it only touches the DOM. S and requestRender (app.js) are used in the
// drag handlers, which can only run once the page is loaded.

const LAYOUT_KEY = 'slideshow-layout';
const LAYOUT_DEFAULT = { umap: true, tissue: true, side: true, umapFrac: 0.5, sideW: 540 };   // umapFrac: share of the map + tissue width
const MIN_PANE_PX = 160;          // narrowest map or tissue pane when dragging
const MIN_SIDE_PX = 260;          // narrowest side panel when dragging
const LAYOUT = { ...LAYOUT_DEFAULT, tissueReady: false };     // tissueReady: set by spatial.js once the tissue view can draw
try { Object.assign(LAYOUT, JSON.parse(localStorage.getItem(LAYOUT_KEY))); } catch (e) { /* private window or nothing stored */ }
LAYOUT.tissueReady = false;

const byId = id => document.getElementById(id);

function saveLayout() {
  const keep = {};
  for (const k of Object.keys(LAYOUT_DEFAULT)) keep[k] = LAYOUT[k];
  try { localStorage.setItem(LAYOUT_KEY, JSON.stringify(keep)); } catch (e) { /* private window */ }
}

function applyLayout() {
  const map = byId('map'), tissue = byId('tissue'), side = byId('side');          // tissue is null when spatial.js removed it
  const tissueOn = LAYOUT.tissue && LAYOUT.tissueReady && tissue !== null;
  const tissueLoading = LAYOUT.tissue && !LAYOUT.tissueReady && tissue !== null;
  if (!LAYOUT.umap && !LAYOUT.side && !tissueOn && !tissueLoading) LAYOUT.umap = true;     // never an empty page
  const umapOn = LAYOUT.umap || byId('loading') !== null;   // the map pane shows the loading message (removed by app.js when done)
  const wrapOn = umapOn || tissueOn;

  map.hidden = !umapOn;
  if (tissue) tissue.hidden = !tissueOn;
  byId('mapwrap').hidden = !wrapOn;
  side.hidden = !LAYOUT.side;
  byId('split-views').hidden = !(umapOn && tissueOn);
  byId('split-side').hidden = !(wrapOn && LAYOUT.side);

  // map and tissue share the width of #views in the ratio umapFrac; the side panel has a fixed width unless it is alone
  map.style.flex = tissueOn ? `${LAYOUT.umapFrac} 1 0px` : '1 1 0px';
  if (tissue) tissue.style.flex = umapOn ? `${1 - LAYOUT.umapFrac} 1 0px` : '1 1 0px';
  side.style.flex = wrapOn ? `0 1 ${LAYOUT.sideW}px` : '1 1 0px';

  document.body.classList.toggle('no-umap', !umapOn);       // hides the toolbar controls that only act on the map
  // the colour legend (also the key of the tissue colours, and where categories are hidden) follows the visible pane
  const legend = byId('legend'), host = umapOn || !tissueOn ? map : tissue;
  if (legend.parentElement !== host) host.appendChild(legend);

  // drawer switches: the last visible panel cannot be switched off
  const on = { umap: LAYOUT.umap, tissue: tissueOn, side: LAYOUT.side };
  const nOn = on.umap + on.tissue + on.side;
  for (const box of document.querySelectorAll('#drawer input[data-panel]')) {
    const p = box.dataset.panel;
    box.checked = p === 'tissue' ? LAYOUT.tissue : on[p];
    box.disabled = (on[p] && nOn === 1) || (p === 'tissue' && !LAYOUT.tissueReady);
  }
  byId('tissuenote').textContent = tissue === null ? 'no spatial data in this bundle' : LAYOUT.tissueReady ? '' : 'loading…';
}

// ---------------------------------------------------------------- drawer
function setDrawer(open) {
  byId('drawer').classList.toggle('open', open);
  byId('menu').setAttribute('aria-expanded', open);
}
byId('menu').addEventListener('click', () => setDrawer(!byId('drawer').classList.contains('open')));
document.addEventListener('pointerdown', e => { if (!e.target.closest('#drawer, #menu')) setDrawer(false); });
// capture phase: Esc closes the drawer and does not also clear the cluster selection (keydown handler of app.js)
window.addEventListener('keydown', e => {
  if (e.key === 'Escape' && byId('drawer').classList.contains('open')) { setDrawer(false); e.stopPropagation(); }
}, true);
for (const box of document.querySelectorAll('#drawer input[data-panel]')) {
  box.addEventListener('change', () => { LAYOUT[box.dataset.panel] = box.checked; applyLayout(); saveLayout(); });
}
byId('layoutreset').addEventListener('click', () => { Object.assign(LAYOUT, LAYOUT_DEFAULT); applyLayout(); saveLayout(); });

// ---------------------------------------------------------------- splitters
// onMove(clientX) updates LAYOUT; the canvases follow through their ResizeObservers (app.js, spatial.js)
function dragSplitter(bar, onMove) {
  bar.addEventListener('pointerdown', e => {
    if (e.button !== 0) return;
    e.preventDefault();
    bar.setPointerCapture(e.pointerId);
    bar.classList.add('dragging'); document.body.classList.add('resizing');
    S.dragging = true;                                      // app.js draws the fast subset of the map meanwhile
  });
  bar.addEventListener('pointermove', e => { if (bar.hasPointerCapture(e.pointerId)) { onMove(e.clientX); applyLayout(); } });
  bar.addEventListener('lostpointercapture', () => {        // pointer up or cancelled
    bar.classList.remove('dragging'); document.body.classList.remove('resizing');
    S.dragging = false; requestRender(); saveLayout();
  });
}
dragSplitter(byId('split-views'), x => {
  const r = byId('views').getBoundingClientRect(), w = r.width - byId('split-views').offsetWidth;
  const lo = Math.min(0.5, MIN_PANE_PX / w);
  LAYOUT.umapFrac = Math.min(1 - lo, Math.max(lo, (x - r.left) / w));
});
dragSplitter(byId('split-side'), x => {
  const r = document.querySelector('main').getBoundingClientRect();
  LAYOUT.sideW = Math.round(Math.max(MIN_SIDE_PX, Math.min(r.width - MIN_PANE_PX, r.right - x)));
});
byId('split-views').addEventListener('dblclick', () => { LAYOUT.umapFrac = LAYOUT_DEFAULT.umapFrac; applyLayout(); saveLayout(); });
byId('split-side').addEventListener('dblclick', () => { LAYOUT.sideW = LAYOUT_DEFAULT.sideW; applyLayout(); saveLayout(); });

applyLayout();
