// ================= Reading aids shared by every mode viewer and band plot =================
// - how many cells a mode shape shows (M × N, user-chosen),
// - motion arrows (out-of-plane: up/down, in-plane: sideways),
// - acoustic / optical branch labels on dispersion curves.
const Aids = (() => {
const TAU = 2 * Math.PI, MAXC = 15;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const st = { M: 3, N: 3, arrows: true, branches: true };
try { const s = JSON.parse(localStorage.getItem('mmw-view') || 'null'); if (s) Object.assign(st, s); } catch (e) { /* storage unavailable */ }
const clampC = v => Math.max(1, Math.min(MAXC, Math.round(+v || 1)));
st.M = clampC(st.M); st.N = clampC(st.N);
function save() { try { localStorage.setItem('mmw-view', JSON.stringify(st)); } catch (e) { /* ignore */ } }

// cells shown in a mode shape, centred on the reference cell (0, 0)
function tiles() { const out = [], i0 = -Math.floor((st.M - 1) / 2), j0 = -Math.floor((st.N - 1) / 2); for (let j = j0; j < j0 + st.N; j++) for (let i = i0; i < i0 + st.M; i++) out.push([i, j]); return out; }
const label = () => `${st.M} × ${st.N}`;
const key = () => `${st.M}x${st.N}`;

// ---------- controls (any element with these classes, on any page) ----------
function sync() {
  const preset = st.M === st.N && [1, 2, 3, 4, 6, 9].includes(st.M) ? String(st.M) : 'custom';
  document.querySelectorAll('select.tilesel').forEach(s => { s.value = preset; });
  document.querySelectorAll('.tilecustom').forEach(el => { el.classList.toggle('hidden', preset !== 'custom'); const [m, n] = el.querySelectorAll('input'); if (m && document.activeElement !== m) m.value = st.M; if (n && document.activeElement !== n) n.value = st.N; });
  document.querySelectorAll('input.arrowchk').forEach(c => { c.checked = st.arrows; });
  document.querySelectorAll('input.branchchk').forEach(c => { c.checked = st.branches; });
}
function changed() { save(); sync(); window.dispatchEvent(new CustomEvent('mmview')); }
function bind() {
  document.querySelectorAll('select.tilesel').forEach(s => s.addEventListener('change', () => {
    if (s.value === 'custom') { document.querySelectorAll('.tilecustom').forEach(el => el.classList.remove('hidden')); return; }
    st.M = st.N = +s.value; changed(); }));
  document.querySelectorAll('.tilecustom input').forEach(inp => inp.addEventListener('change', () => {
    const box = inp.closest('.tilecustom'), [m, n] = box.querySelectorAll('input'); st.M = clampC(m.value); st.N = clampC(n.value); changed(); }));
  document.querySelectorAll('input.arrowchk').forEach(c => c.addEventListener('change', () => { st.arrows = c.checked; changed(); }));
  document.querySelectorAll('input.branchchk').forEach(c => c.addEventListener('change', () => { st.branches = c.checked; changed(); }));
  sync();
}
// ---------- colour theme: auto (system) / light / dark ----------
let theme = 'auto'; try { theme = localStorage.getItem('mmw-theme') || 'auto'; } catch (e) { /* ignore */ }
function applyTheme() { if (theme === 'auto') document.documentElement.removeAttribute('data-theme'); else document.documentElement.setAttribute('data-theme', theme);
  const b = document.getElementById('themeBtn'); if (b) b.textContent = theme === 'auto' ? 'Theme: auto' : theme === 'light' ? 'Theme: light' : 'Theme: dark'; }
applyTheme();
function bindTheme() { const b = document.getElementById('themeBtn'); if (!b) return; applyTheme();
  b.addEventListener('click', () => { theme = theme === 'auto' ? 'light' : theme === 'light' ? 'dark' : 'auto'; try { localStorage.setItem('mmw-theme', theme); } catch (e) { /* ignore */ } applyTheme(); window.dispatchEvent(new CustomEvent('mmview')); window.dispatchEvent(new Event('resize')); }); }
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', () => { bind(); bindTheme(); }); else { bind(); bindTheme(); }

// ---------- arrows ----------
function arrow(ctx, ax, ay, bx, by, col, lw = 1.7) {
  const dx = bx - ax, dy = by - ay, L = Math.hypot(dx, dy); if (L < 4) return;
  const hd = Math.min(9, 0.45 * L + 2), ux = dx / L, uy = dy / L;
  ctx.strokeStyle = col; ctx.fillStyle = col; ctx.lineWidth = lw; ctx.lineCap = 'round';
  ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(bx - ux * hd * 0.7, by - uy * hd * 0.7); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(bx, by); ctx.lineTo(bx - ux * hd - uy * hd * 0.45, by - uy * hd + ux * hd * 0.45); ctx.lineTo(bx - ux * hd + uy * hd * 0.45, by - uy * hd - ux * hd * 0.45); ctx.closePath(); ctx.fill();
}
const colOOP = () => css('--arr-oop') || '#CC3311', colIP = () => css('--arr-ip') || '#0077BB';
// a few representative points per cell: one per bin of a coarse grid (the one that moves most), plus the given ones (masses)
function pickNodes(pts, mags, h, must = []) {
  const best = new Map();
  pts.forEach((p, i) => { if (!p) return; const k = Math.floor(p[0] / h) + ',' + Math.floor(p[1] / h); const b = best.get(k); if (b == null || mags[i] > mags[b]) best.set(k, i); });
  const s = new Set([...best.values(), ...must]), mx = Math.max(1e-30, ...[...s].map(i => mags[i] || 0));
  return [...s].filter(i => (mags[i] || 0) > 0.12 * mx);
}
// small key in the top-right corner of a mode canvas
function arrowKey(ctx, w, kinds) {
  if (!st.arrows) return;
  ctx.save(); ctx.font = '12px system-ui, sans-serif'; ctx.textAlign = 'left';
  const rows = kinds.map(k => k === 'oop' ? ['Out-of-plane motion (up / down)', colOOP(), 'v'] : ['In-plane motion (sideways)', colIP(), 'h']);
  const tw = Math.max(...rows.map(r => ctx.measureText(r[0]).width)) + 38, x0 = w - tw - 10; let y = 12;
  ctx.fillStyle = css('--panel'); ctx.globalAlpha = 0.88; ctx.fillRect(x0 - 6, 4, tw + 10, rows.length * 18 + 8); ctx.globalAlpha = 1;
  for (const [t, c, d] of rows) { if (d === 'v') arrow(ctx, x0 + 12, y + 13, x0 + 12, y - 1, c, 1.6); else arrow(ctx, x0 + 2, y + 6, x0 + 24, y + 6, c, 1.6); ctx.fillStyle = css('--ink-2'); ctx.fillText(t, x0 + 32, y + 10); y += 18; }
  ctx.restore();
}

// ---------- acoustic / optical branches ----------
// Acoustic branches start at 0 Hz at Γ (the whole cell moves together at long wavelength);
// optical branches start above zero (parts of the cell move against each other).
// nAc: number of acoustic branches of the waves shown: 1 out-of-plane, 2 in-plane, 3 for all waves of a plate-like lattice.
function gammaIdx(ks) { if (!ks) return []; const s = Math.max(1e-30, ...ks.map(k => Math.hypot(k[0], k[1]))); return ks.map((k, i) => Math.hypot(k[0], k[1]) < 1e-6 * s ? i : -1).filter(i => i >= 0); }
function ranks(F, keep) { return F.map((row, i) => row.map((f, b) => ({ f, b })).filter(o => Number.isFinite(o.f) && (!keep || keep(i, o.b))).sort((a, c) => a.f - c.f).map(o => o.b)); }
function branchOf(F, ks, keep, nAc, i, b) {
  if (!gammaIdx(ks).length) return '';
  const r = ranks([F[i]], keep ? (_, bb) => keep(i, bb) : null)[0].indexOf(b); if (r < 0) return '';
  return r < nAc ? `${nAc > 1 ? `acoustic branch ${r + 1} of ${nAc}` : 'acoustic branch'}: the whole cell moves together` : `optical branch ${r - nAc + 1}: parts of the cell move against each other`;
}
// two-line label: bold name, plain-language explanation underneath; (x, y) is the bottom-left (or bottom-right) corner
function pill(ctx, txt, sub, x, y, col, align = 'left') {
  ctx.font = '600 12px system-ui, sans-serif'; const t1 = ctx.measureText(txt).width; ctx.font = '11.5px system-ui, sans-serif'; const t2 = sub ? ctx.measureText(sub).width : 0;
  const tw = Math.max(t1, t2), hh = sub ? 32 : 18, x0 = align === 'left' ? x : x - tw - 12, y0 = y - hh;
  ctx.globalAlpha = 0.92; ctx.fillStyle = css('--panel'); ctx.fillRect(x0, y0, tw + 12, hh); ctx.globalAlpha = 1;
  ctx.strokeStyle = col; ctx.lineWidth = 1; ctx.strokeRect(x0 + 0.5, y0 + 0.5, tw + 11, hh - 1); ctx.textAlign = 'left';
  ctx.fillStyle = col; ctx.font = '600 12px system-ui, sans-serif'; ctx.fillText(txt, x0 + 6, y0 + 13);
  if (sub) { ctx.fillStyle = css('--ink-2'); ctx.font = '11.5px system-ui, sans-serif'; ctx.fillText(sub, x0 + 6, y0 + 27); }
}
// draws a soft band behind the acoustic branch(es) and labels both kinds; call before the band dots
function drawBranches(ctx, X, Y, xs, F, ks, keep, nAc, opt = {}) {
  if (!st.branches) return; const G = gammaIdx(ks); if (!G.length) return;
  const R = ranks(F, keep), ac = opt.color || css('--ac') || '#44AA99', fmax = opt.fmax || Infinity;
  ctx.save(); ctx.globalAlpha = 0.28; ctx.strokeStyle = ac; ctx.lineWidth = 7; ctx.lineJoin = ctx.lineCap = 'round';
  for (let r = 0; r < nAc; r++) { ctx.beginPath(); let on = false; xs.forEach((x, i) => { const b = R[i][r]; const f = b == null ? NaN : F[i][b]; if (!Number.isFinite(f) || f > fmax) { on = false; return; } const px = X(x), py = Y(f); on ? ctx.lineTo(px, py) : ctx.moveTo(px, py); on = true; }); ctx.stroke(); }
  ctx.restore(); if (opt.noLabels) return;
  // labels: "Acoustic" near the middle of the first leg, "Optical" just off Γ at the lowest optical branch
  const g = G[0], dir = g < xs.length / 2 ? 1 : -1, n = xs.length, iA = Math.max(0, Math.min(n - 1, g + dir * Math.round(0.18 * n)));
  const bA = R[iA][nAc - 1]; if (bA != null && F[iA][bA] <= fmax) pill(ctx, opt.tag ? `Acoustic (${opt.tag})` : 'Acoustic', 'whole cell moves together', X(xs[iA]) + 6 * dir, Y(F[iA][bA]) - 7, ac, dir > 0 ? 'left' : 'right');
  const bO = R[g][nAc]; if (bO != null && F[g][bO] <= fmax) pill(ctx, opt.tag ? `Optical (${opt.tag})` : 'Optical', 'parts move against each other', X(xs[g]) + 18 * dir, Y(F[g][bO]) - 7, css('--ink-2'), dir > 0 ? 'left' : 'right');
}

// ---------- gap names: which resonator mode opens a gap, and the bands at its edges ----------
// share(g): mean resonator share of kinetic energy in the band just below the gap (0..1), or null when unknown.
// A gap whose lower band is mostly resonator motion is a local-resonance gap; the n-th such gap (in frequency) is
// credited to the n-th resonator mode. Otherwise the gap comes from the periodicity of the host (Bragg).
function tagGaps(gaps, share) {
  let n = 0; [...(gaps || [])].sort((a, c) => a.lo - c.lo).forEach(g => { const s = share ? share(g) : null;
    if (s == null || !Number.isFinite(s)) { g.kind = null; g.lrn = null; } else if (s > 0.5) { g.kind = 'lr'; g.lrn = ++n; } else { g.kind = 'bragg'; g.lrn = null; } g.share = s; });
  return gaps;
}
function gapName(g, full) {
  const bands = g.below ? ` · bands ${g.below}–${g.below + 1}` : '';
  if (full) return `Full gap${bands}`;
  if (g.kind === 'lr') return `Resonance gap · resonator mode ${g.lrn}${bands}`;
  if (g.kind === 'bragg') return `Bragg gap (host periodicity)${bands}`;
  return `Gap${bands}`;
}
function gapShort(g) { return g.kind === 'lr' ? `resonator mode ${g.lrn}` : g.kind === 'bragg' ? 'Bragg (host periodicity)' : ''; }
// labels drawn at the right-hand side of the plot, inside each gap (or just above a thin one)
function drawGapLabels(ctx, gaps, Y, xR, fmin, fmax, full) {
  if (!st.branches) return; let last = Infinity;
  ctx.save(); ctx.font = '600 11.5px system-ui, sans-serif';
  for (const g of [...(gaps || [])].sort((a, c) => a.lo - c.lo)) { if (g.hi < fmin || g.lo > fmax) continue;
    const yl = Y(Math.max(g.lo, fmin)), yh = Y(Math.min(g.hi, fmax)), thin = yl - yh < 18; let y = thin ? yh - 5 : (yl + yh) / 2 + 4;
    if (y > last - 19) y = last - 19; last = y;
    const t = gapName(g, full), tw = ctx.measureText(t).width, x0 = xR - tw - 18;
    ctx.globalAlpha = 0.92; ctx.fillStyle = css('--panel'); ctx.fillRect(x0, y - 12, tw + 12, 17); ctx.globalAlpha = 1;
    ctx.strokeStyle = css('--gap'); ctx.lineWidth = 1.2; ctx.strokeRect(x0 + 0.5, y - 11.5, tw + 11, 16);
    ctx.fillStyle = css('--ink'); ctx.textAlign = 'left'; ctx.fillText(t, x0 + 6, y); }
  ctx.restore();
}

// ---------- zoom on the frequency axis (band structure and complex band share it) ----------
// state per group as fractions of the full range; wheel zooms at the cursor, drag pans, double-click resets
const zs = {}, zDraw = {}, MT = 16, MB = 36;
function fz(g, lo, hi) { const z = zs[g]; return z ? [lo + z.a * (hi - lo), lo + z.b * (hi - lo)] : [lo, hi]; }
function zoomed(g) { return !!zs[g]; }
function zset(g, a, span) { span = Math.max(0.002, Math.min(1, span)); a = Math.max(0, Math.min(1 - span, a)); if (span >= 0.999) delete zs[g]; else zs[g] = { a, b: a + span }; (zDraw[g] || []).forEach(f => f()); syncZoom(g); }
function zoomBy(g, fac, t = 0.5) { const z = zs[g] || { a: 0, b: 1 }, span = z.b - z.a, c = z.a + t * span, ns = span * fac; zset(g, c - t * ns, ns); }
function zreset(g) { delete zs[g]; (zDraw[g] || []).forEach(f => f()); syncZoom(g); }
function syncZoom(g) { document.querySelectorAll(`.zoombar[data-zg="${g}"]`).forEach(el => el.classList.toggle('on', !!zs[g])); }
function zoomable(cv, g, redraw) {
  if (!cv || cv._zoom) return; cv._zoom = g; (zDraw[g] = zDraw[g] || []).push(redraw);
  const frac = e => { const r = cv.getBoundingClientRect(); return Math.max(0, Math.min(1, (r.height - MB - (e.clientY - r.top)) / (r.height - MT - MB))); };
  cv.addEventListener('wheel', e => { e.preventDefault(); zoomBy(g, e.deltaY < 0 ? 0.8 : 1.25, frac(e)); }, { passive: false });
  let d = null;
  cv.addEventListener('pointerdown', e => { if (e.button !== 0) return; d = { y: e.clientY, z: zs[g] ? { ...zs[g] } : { a: 0, b: 1 }, moved: false }; });
  cv.addEventListener('pointermove', e => { if (!d) return; const dy = e.clientY - d.y; if (!d.moved && Math.abs(dy) < 5) return; if (!d.moved) { d.moved = true; try { cv.setPointerCapture(e.pointerId); } catch (er) { /* ignore */ } }
    const r = cv.getBoundingClientRect(), span = d.z.b - d.z.a; if (span >= 0.999) return; zset(g, d.z.a + dy / (r.height - MT - MB) * span, span); cv.style.cursor = 'grabbing'; });
  const up = () => { if (d && d.moved) cv._dragged = true; d = null; cv.style.cursor = ''; };
  cv.addEventListener('pointerup', up); cv.addEventListener('pointercancel', up);
  cv.addEventListener('click', e => { if (cv._dragged) { cv._dragged = false; e.stopImmediatePropagation(); } }, true);
  cv.addEventListener('dblclick', e => { e.preventDefault(); zreset(g); });
}
function bindZoomBars() { document.querySelectorAll('.zoombar').forEach(el => { const g = el.dataset.zg; el.querySelectorAll('button').forEach(b => b.addEventListener('click', () => { const z = b.dataset.z; if (z === 'in') zoomBy(g, 0.6); else if (z === 'out') zoomBy(g, 1 / 0.6); else zreset(g); })); }); }
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', bindZoomBars); else bindZoomBars();
return { st, tiles, label, key, arrow, colOOP, colIP, pickNodes, arrowKey, branchOf, drawBranches, gammaIdx, sync, tagGaps, gapName, gapShort, drawGapLabels, fz, zoomed, zoomable, zreset };
})();
