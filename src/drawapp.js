// ================= "Draw cell" tab: sketch a unit cell (beams, shell walls, shell plates with holes, masses, resonators) =================
const DrawUI = (() => {
const B = window.MMW, $ = s => document.querySelector(s), TAU = 2 * Math.PI, css = B.css;
const num = id => parseFloat(document.getElementById(id).value);
const ST_KEY = 'mmw-draw-v1';
const D = { roleAuto: true, eview: { yaw: -0.5, elev: 0.85 }, mview: { yaw: -0.6, elev: 0.9 }, dr: null, tool: 'select', role: 'host', sel: null, pending: null, cur: null, view: null, hist: [], fut: [], run: null, md: null, mdKey: null, yaw: -0.6, elev: 0.9, show: 'flex', mesh: true };
const DEF = { beam: { b: 1, h: 5 }, wall: { t: 1, H: 5, z0: -2.5 }, plate: { t: 1, z: 0 }, mass: { m: 0.5, side: 4, z: 0 }, spring: { m: 0.5, f: 800, z: 0 } };
const TOOLS = [['select', 'Select / move'], ['shape', 'Shape'], ['beam', 'Beam'], ['wall', 'Shell wall'], ['plate', 'Free plate'], ['hole', 'Free hole'], ['mass', 'Point mass'], ['spring', 'Mass-spring'], ['template', 'Resonator template']];

// ---------- state ----------
const clone = o => JSON.parse(JSON.stringify(o));
function save() { try { localStorage.setItem(ST_KEY, JSON.stringify(D.dr)); } catch (e) {} }
function commit() { D.hist.push(JSON.stringify(D.dr)); if (D.hist.length > 80) D.hist.shift(); D.fut = []; D.pendingEdit = true; }
// new parts are resonator parts as soon as the cell has a frame (host parts); before that they are frame parts
function syncRole() { if (!D.roleAuto || !D.dr) return; const r = D.dr.parts.some(p => p.role !== 'res' && p.kind !== 'spring') ? 'res' : 'host'; if (r === D.role) return; D.role = r; document.querySelectorAll('[data-drole]').forEach(x => x.setAttribute('aria-checked', x.dataset.drole === r)); }
// which parts are the lattice frame: parts the Lattice box built, or else beams and walls lying along the boundary of the Wigner–Seitz cell
function isFramePart(p) {
  if (p.frame) return true; if (p.kind !== 'beam' && p.kind !== 'wall' && !(p.kind === 'shape' && (p.use === 'frame' || p.use === 'walls'))) return false;
  const a1 = D.dr.a1, a2 = D.dr.a2, L = Math.max(Math.hypot(...a1), Math.hypot(...a2)), R = []; for (let i = -2; i <= 2; i++) for (let j = -2; j <= 2; j++) if (i || j) R.push([i * a1[0] + j * a2[0], i * a1[1] + j * a2[1]]);
  const onWS = q => { const d0 = Math.hypot(q[0], q[1]); let dm = Infinity; for (const r of R) dm = Math.min(dm, Math.hypot(q[0] - r[0], q[1] - r[1])); return Math.abs(d0 - dm) < 0.02 * L; };
  return segsOf(p).every(([a, b]) => onWS([(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]));
}
function autoRoles() {
  const fr = D.dr.parts.map(isFramePart); if (!fr.some(Boolean)) { B.setStatus('No lattice frame found: mark the frame parts as Host and the rest as Resonator one by one (Select, then Role).'); return; }
  commit(); let nh = 0, nr = 0; D.dr.parts.forEach((p, i) => { if (p.kind === 'shape' && p.use === 'hole') return; if (fr[i]) { p.role = 'host'; nh++; } else if (p.kind !== 'spring') { p.role = 'res'; nr++; } }); changed();
  B.setStatus(`Roles set: ${nh} frame part${nh === 1 ? '' : 's'} as host (blue), ${nr} part${nr === 1 ? '' : 's'} as resonator (orange). Check the colours and change any part with Select → Role.`);
}
function changed(fromGen) { if (D.dr && D.dr.latFollow) fitLattice(true); latCheck(); syncRole(); if (D.pendingEdit && !fromGen && D.dr) D.dr._edited = true; D.pendingEdit = false; if (typeof genNote === 'function') genNote(); save(); D.md = null; renderProps(); drawEditor(); scheduleModel(); sweepParts(); if (B.drawKPath) B.drawKPath(); }
// ---------- lattice vectors from the frame outline ----------
// A shape tiles the plane by translation alone only if it is a parallelogram (square, rectangle, rhombus) or a
// centrally symmetric hexagon (regular hexagon included); then the lattice vectors map each side onto the opposite one.
function outlineOf(p) {
  if (!p) return null;
  if (p.kind === 'shape') return p.type === 'circle' ? { circle: true } : Drw.shapePoly(p);
  if (p.kind === 'plate') return p.poly.map(q => [q[0], q[1]]);
  return null;
}
function latticeFromOutline(P) {
  if (!P || P.circle) return { err: 'A circle cannot fill the plane without gaps. Draw a square or hexagon frame around it and fit to that.' };
  const n = P.length, sub = (a, b) => [a[0] - b[0], a[1] - b[1]], d = (u, v) => Math.hypot(u[0] - v[0], u[1] - v[1]);
  const size = Math.max(...P.map(q => d(q, P[0]))), tol = 1e-4 * size + 1e-7, r6 = v => +v.toFixed(6);
  let a1 = null, a2 = null, kind = '';
  if (n === 4 && d(sub(P[1], P[0]), sub(P[2], P[3])) < tol && d(sub(P[3], P[0]), sub(P[2], P[1])) < tol) { a1 = sub(P[1], P[0]); a2 = sub(P[3], P[0]); kind = Math.abs(a1[0] * a2[0] + a1[1] * a2[1]) < tol * size ? (Math.abs(Math.hypot(...a1) - Math.hypot(...a2)) < tol ? 'square' : 'rectangle') : 'parallelogram'; }
  else if (n === 6) { const c = [P.reduce((s, q) => s + q[0], 0) / 6, P.reduce((s, q) => s + q[1], 0) / 6];
    if ([0, 1, 2].every(k => d([P[k][0] + P[k + 3][0], P[k][1] + P[k + 3][1]], [2 * c[0], 2 * c[1]]) < tol)) { a1 = sub(P[0], P[4]); a2 = sub(P[1], P[5]); kind = 'hexagon'; } }
  if (!a1) {
    if (n === 3) return { err: 'A triangle cannot tile the plane by translation alone: neighbours must be turned upside down. Use Geometry → Build your own → Triangular lattice, whose cell holds two triangles.' };
    if (n === 4 || n === 6) return { err: `This ${n === 4 ? 'four' : 'six'}-sided outline is not ${n === 4 ? 'a parallelogram' : 'centrally symmetric'} (opposite sides must be equal and parallel), so copies of it cannot fill the plane.` };
    return { err: `A ${n}-sided outline cannot fill the plane by translation. Only squares, rectangles, parallelograms and (centrally symmetric) hexagons can; use one of those as the frame.` };
  }
  const pos = v => (v[0] < -tol || (Math.abs(v[0]) <= tol && v[1] < 0)) ? [-v[0], -v[1]] : v; a1 = pos(a1); a2 = pos(a2); // point to the right (or up)
  if (Math.abs(a2[1]) <= tol && Math.abs(a1[1]) > tol) [a1, a2] = [a2, a1]; // a horizontal vector becomes a₁
  if (a1[0] * a2[1] - a1[1] * a2[0] < 0) { if (Math.abs(a1[1]) <= tol) a2 = [-a2[0], -a2[1]]; else [a1, a2] = [a2, a1]; }
  return { a1: a1.map(r6), a2: a2.map(r6), kind };
}
// the part to fit: the selected shape or plate, else the frame part marked earlier, else the largest host outline
function frameOutlinePart() {
  const ps = D.dr.parts, ok = p => p && (p.kind === 'shape' || p.kind === 'plate');
  if (D.sel != null && ok(ps[D.sel])) return ps[D.sel];
  const marked = ps.find(p => p.latFrame && ok(p)); if (marked) return marked;
  let best = null, ba = 0; for (const p of ps) { if (!ok(p) || p.role === 'res') continue; const P = outlineOf(p); if (!P || P.circle) continue; const A = Math.abs(Drw.polyArea(P)); if (A > ba) { ba = A; best = p; } }
  return best;
}
function fitLattice(quiet, part) {
  const msg = $('#dLatMsg'), p = part || frameOutlinePart();
  if (!p) { if (!quiet && msg) msg.innerHTML = '<span class="err">No frame outline found.</span> Draw the frame as a shape (square, rectangle or hexagon, as plate, beam frame or wall frame) or as a plate outline, select it, then press the button.'; return false; }
  const L = latticeFromOutline(outlineOf(p));
  if (L.err) { if (msg) msg.innerHTML = `<span class="err">Cannot fit:</span> ${L.err}`; return false; }
  const same = D.dr.a1.every((v, i) => Math.abs(v - L.a1[i]) < 1e-9) && D.dr.a2.every((v, i) => Math.abs(v - L.a2[i]) < 1e-9);
  D.dr.parts.forEach(q => { delete q.latFrame; }); p.latFrame = true;
  if (!same) { if (!quiet) commit(); D.dr.a1 = L.a1; D.dr.a2 = L.a2; D.view = null; syncLat(); }
  if (msg) msg.innerHTML = `Fitted to the ${L.kind} frame: a₁ = (${L.a1.join(', ')}), a₂ = (${L.a2.join(', ')}) mm, |a₁| = ${Math.hypot(...L.a1).toFixed(3)} mm, |a₂| = ${Math.hypot(...L.a2).toFixed(3)} mm. Neighbouring cells share the frame's sides.`;
  return !same;
}
// the host frame outline whose copies should tile the plane (ignores the selection)
function hostFrame() {
  const ps = D.dr.parts, anyHost = ps.some(p => p.role !== 'res' && p.kind !== 'spring' && p.use !== 'hole'), ok = p => p && (p.kind === 'shape' || p.kind === 'plate') && (p.role !== 'res' || !anyHost) && p.use !== 'hole';
  const marked = ps.find(p => p.latFrame && ok(p)); if (marked) return marked;
  let best = null, ba = 0; for (const p of ps) { if (!ok(p)) continue; const P = outlineOf(p); if (!P || P.circle) continue; const A = Math.abs(Drw.polyArea(P)); if (A > ba) { ba = A; best = p; } }
  return best;
}
// same lattice? (bases may differ: each vector of one must be an integer combination of the other)
function sameLattice(a1, a2, b1, b2) {
  const det = a1[0] * a2[1] - a1[1] * a2[0]; if (Math.abs(det) < 1e-12) return false;
  const coef = v => [(v[0] * a2[1] - v[1] * a2[0]) / det, (a1[0] * v[1] - a1[1] * v[0]) / det], int = x => Math.abs(x - Math.round(x)) < 1e-4;
  const db = b1[0] * b2[1] - b1[1] * b2[0]; return Math.abs(Math.abs(db) - Math.abs(det)) < 1e-4 * Math.abs(det) && [b1, b2].every(v => coef(v).every(int));
}
// warning above the drawing when the frame and the lattice vectors do not fit together
function latCheck() {
  const el = document.getElementById('dLatWarn'); if (!el || !D.dr) return;
  const p = hostFrame(), L = p ? latticeFromOutline(outlineOf(p)) : null, bad = !!(L && !L.err && !sameLattice(D.dr.a1, D.dr.a2, L.a1, L.a2));
  el.classList.toggle('hidden', !bad); D.latPart = bad ? p : null;
  if (bad) el.querySelector('span').innerHTML = `<b>The lattice does not match the ${L.kind} frame:</b> copies of the cell overlap or leave gaps, so the walls do not join. For this frame a₁ = (${L.a1.join(', ')}) and a₂ = (${L.a2.join(', ')}) mm.`;
}
function undo() { if (!D.hist.length) return; D.fut.push(JSON.stringify(D.dr)); D.dr = JSON.parse(D.hist.pop()); D.sel = null; changed(); }
function redo() { if (!D.fut.length) return; D.hist.push(JSON.stringify(D.dr)); D.dr = JSON.parse(D.fut.pop()); D.sel = null; changed(); }
function importLeft() {
  const P = B.readParams(B.lastLat()), err = B.validate(P); if (err) { B.setStatus(err); return false; }
  try { const dr = Drw.fromSpec(B.buildSpec(P)); dr.mesh = D.dr ? D.dr.mesh : 2; dr.grid = D.dr ? D.dr.grid : 0.5;
    // the workbench's thin out-of-plane walls become beams; keep them
    dr._design = 'own'; dr._edited = true; dr._gen = ''; commit(); D.dr = dr; D.sel = null; D.view = null; syncLat(); changed(); B.setStatus('Left-panel cell imported into the drawing: every wall, arm and block is now an editable part.'); return true; } catch (e) { B.setStatus(e.message); return false; }
}
function syncLat() { const d = D.dr; $('#dA1x').value = d.a1[0]; $('#dA1y').value = d.a1[1]; $('#dA2x').value = d.a2[0]; $('#dA2y').value = d.a2[1]; $('#dMesh').value = d.mesh || 2; $('#dGrid').value = d.grid != null ? d.grid : 0.5; const fl = document.getElementById('dLatFollow'); if (fl) fl.checked = !!d.latFollow; }

// ---------- view transform ----------
function fitView(w, h) {
  const d = D.dr, c = [[0, 0], d.a1, [d.a1[0] + d.a2[0], d.a1[1] + d.a2[1]], d.a2]; for (const p of d.parts) try { for (const v of vertsOf(p)) c.push(v); } catch (e) {}
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity; for (const p of c) { x0 = Math.min(x0, p[0]); x1 = Math.max(x1, p[0]); y0 = Math.min(y0, p[1]); y1 = Math.max(y1, p[1]); }
  const pad = 0.55 * Math.max(x1 - x0, y1 - y0), sc = Math.min(w / (x1 - x0 + 2 * pad), h / (y1 - y0 + 2 * pad));
  D.view = { sc, cx: (x0 + x1) / 2, cy: (y0 + y1) / 2, w, h };
}
const toS = p => [D.view.w / 2 + (p[0] - D.view.cx) * D.view.sc, D.view.h / 2 - (p[1] - D.view.cy) * D.view.sc];
const toW = (x, y) => [D.view.cx + (x - D.view.w / 2) / D.view.sc, D.view.cy - (y - D.view.h / 2) / D.view.sc];
const imgs = () => { const d = D.dr, o = []; for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) o.push([i * d.a1[0] + j * d.a2[0], i * d.a1[1] + j * d.a2[1], i, j]); return o; };

// ---------- geometry of parts (2D) ----------
function vertsOf(p) { if (p.kind === 'shape') return Drw.shapePoly(p).concat([p.c]); if (p.kind === 'beam' || p.kind === 'wall') return p.p.map(q => [q[0], q[1]]); if (p.kind === 'plate') return p.poly.concat(...(p.holes || [])); return [[p.p[0], p.p[1]]]; }
function segsOf(p) { const s = []; if (p.kind === 'shape') { const r = Drw.shapePoly(p); for (let i = 0; i < r.length; i++) s.push([r[i], r[(i + 1) % r.length]]); return s; } if (p.kind === 'beam' || p.kind === 'wall') s.push([p.p[0], p.p[1]]); if (p.kind === 'plate') for (const r of [p.poly, ...(p.holes || [])]) for (let i = 0; i < r.length; i++) s.push([r[i], r[(i + 1) % r.length]]); return s; }
// handles a part shows when selected: shapes have a centre handle and a size/rotation handle; other parts their points
function handlesOf(p) { if (p.kind !== 'shape') return vertsOf(p); const poly = Drw.shapePoly(p); return [p.c, p.type === 'rect' ? poly[2] : poly[0]]; }
function setHandle(p, k, q) {
  if (p.kind !== 'shape') return setVertex(p, k, q);
  if (k === 0) { p.c = [q[0], q[1]]; return; }
  const dx = q[0] - p.c[0], dy = q[1] - p.c[1], r4 = v => +v.toFixed(4);
  if (p.type === 'rect') { const a = (p.ang || 0) * Math.PI / 180, lx = Math.cos(a) * dx + Math.sin(a) * dy, ly = -Math.sin(a) * dx + Math.cos(a) * dy; p.sx = r4(Math.max(1e-3, 2 * Math.abs(lx))); p.sy = r4(Math.max(1e-3, 2 * Math.abs(ly))); }
  else { p.R = r4(Math.max(1e-3, Math.hypot(dx, dy))); p.ang = +(Math.atan2(dy, dx) * 180 / Math.PI).toFixed(2); }
}
function centreOf(p) { if (p.kind === 'shape') return p.c.slice(); const v = p.kind === 'plate' ? p.poly : vertsOf(p); return [v.reduce((s, q) => s + q[0], 0) / v.length, v.reduce((s, q) => s + q[1], 0) / v.length]; }
function snap(w, exclude) {
  const tolW = 9 / D.view.sc; let best = null;
  for (const [sx, sy] of imgs()) D.dr.parts.forEach((p, pi) => { for (const v of vertsOf(p)) { if (exclude && exclude.pi === pi) continue; const q = [v[0] + sx, v[1] + sy], d = Math.hypot(q[0] - w[0], q[1] - w[1]); if (d < tolW && (!best || d < best.d)) best = { d, p: q, kind: 'vertex' }; } });
  if (D.pending && D.pending.pts.length) { const q = D.pending.pts[0], d = Math.hypot(q[0] - w[0], q[1] - w[1]); if (d < tolW && (!best || d < best.d)) best = { d, p: q.slice(), kind: 'vertex' }; }
  if (!best) for (const [sx, sy] of imgs()) D.dr.parts.forEach((p, pi) => { if (exclude && exclude.pi === pi) return; for (const [a, b] of segsOf(p)) { const A = [a[0] + sx, a[1] + sy], Bq = [b[0] + sx, b[1] + sy], r = Drw.segDist(w, A, Bq); if (r.d < 0.7 * tolW && (!best || r.d < best.d)) best = { d: r.d, p: [A[0] + (Bq[0] - A[0]) * r.t, A[1] + (Bq[1] - A[1]) * r.t], kind: 'edge' }; } });
  if (!best) { const g = D.dr.grid; best = { p: g > 0 ? [Math.round(w[0] / g) * g, Math.round(w[1] / g) * g] : w.slice(), kind: 'grid' }; }
  best.p = best.p.map(v => +v.toFixed(6)); return best;
}
function hit(w) {
  const tolW = 7 / D.view.sc; let best = null;
  D.dr.parts.forEach((p, pi) => {
    let d = Infinity;
    if (p.kind === 'beam' || p.kind === 'wall') d = Drw.segDist(w, p.p[0], p.p[1]).d - (p.kind === 'beam' ? p.b / 2 : p.t / 2);
    else if (p.kind === 'shape') { const poly = Drw.shapePoly(p), ed = Math.min(...segsOf(p).map(s => Drw.segDist(w, s[0], s[1]).d)); d = (p.use === 'plate' || p.use === 'hole') && Drw.pointInPoly(w, poly) ? Math.min(tolW * 0.85, ed) : ed - (p.use === 'frame' ? (p.b || 1) / 2 : p.use === 'walls' ? (p.t || 1) / 2 : 0); }
    else if (p.kind === 'plate') { d = Drw.pointInPoly(w, p.poly) && !(p.holes || []).some(hh => Drw.pointInPoly(w, hh)) ? tolW * 0.9 : Math.min(...segsOf(p).map(s => Drw.segDist(w, s[0], s[1]).d)); }
    else d = Math.hypot(w[0] - p.p[0], w[1] - p.p[1]) - (p.side || 2) / 2;
    const area = p.kind === 'plate' || (p.kind === 'shape' && (p.use === 'plate' || p.use === 'hole')); if (d < tolW && (!best || (area === best.area ? d < best.d : !area))) best = { d, pi, kind: p.kind, area };
  });
  return best;
}

// ---------- editor drawing ----------
function drawEditor() {
  const cv = $('#cvDraw'); if (!cv || $('#p-drw').classList.contains('hidden')) return;
  const { ctx, w, h } = B.prep(cv); if (!D.view || D.view.w !== w || D.view.h !== h) { const v = D.view; fitView(w, h); if (v && v.keep) Object.assign(D.view, { sc: v.sc, cx: v.cx, cy: v.cy, keep: true }); }
  const d = D.dr, host = css('--host'), res = css('--res'), ink2 = css('--ink-2'), panel = css('--panel');
  // grid
  const g = d.grid > 0 ? d.grid : 1; let step = g; while (step * D.view.sc < 8) step *= 2;
  const [wx0, wy1] = toW(0, 0), [wx1, wy0] = toW(w, h);
  ctx.strokeStyle = css('--grid'); ctx.lineWidth = 1; ctx.beginPath();
  for (let x = Math.ceil(wx0 / step) * step; x <= wx1; x += step) { const s = toS([x, 0])[0]; ctx.moveTo(s, 0); ctx.lineTo(s, h); }
  for (let y = Math.ceil(wy0 / step) * step; y <= wy1; y += step) { const s = toS([0, y])[1]; ctx.moveTo(0, s); ctx.lineTo(w, s); }
  ctx.stroke();
  // mesh preview (shell triangles) and parts, neighbours faded
  const order = imgs().sort((a, b) => (a[2] || a[3] ? 0 : 1) - (b[2] || b[3] ? 0 : 1));
  for (const [sx, sy, i, j] of order) {
    const main = !i && !j; ctx.globalAlpha = main ? 1 : 0.28;
    const S = p => toS([p[0] + sx, p[1] + sy]);
    d.parts.forEach((p, pi) => { const shp = p.kind === 'shape' && (p.use === 'plate' || p.use === 'hole'); if (p.kind !== 'plate' && !shp) return; const col = p.role === 'res' ? res : host;
      if (shp && p.use === 'hole') { ctx.beginPath(); Drw.shapePoly(p).forEach((q, k) => { const s = S(q); k ? ctx.lineTo(s[0], s[1]) : ctx.moveTo(s[0], s[1]); }); ctx.closePath(); ctx.fillStyle = panel; ctx.fill(); ctx.setLineDash([4, 3]); ctx.strokeStyle = main && D.sel === pi ? css('--ink') : ink2; ctx.lineWidth = main && D.sel === pi ? 2.5 : 1.2; ctx.stroke(); ctx.setLineDash([]); return; }
      ctx.beginPath(); for (const r of shp ? [Drw.shapePoly(p)] : [p.poly, ...(p.holes || [])]) r.forEach((q, k) => { const s = S(q); k ? ctx.lineTo(s[0], s[1]) : ctx.moveTo(s[0], s[1]); }), ctx.closePath();
      ctx.fillStyle = B.mix(col, panel, 0.62); ctx.fill('evenodd'); ctx.strokeStyle = col; ctx.lineWidth = main && D.sel === pi ? 2.5 : 1.2; ctx.stroke(); });
    if (main && D.mesh && D.md) { ctx.strokeStyle = B.mix(host, panel, 0.35); ctx.lineWidth = 0.5; ctx.beginPath(); for (const e of D.md.elems) if (e.type === 'shell' && Math.abs(e.Xd[0][2] - e.Xd[1][2]) + Math.abs(e.Xd[0][2] - e.Xd[2][2]) < 1e-9) { const P = e.Xd.map(x => toS([x[0] * 1000, x[1] * 1000])); ctx.moveTo(P[0][0], P[0][1]); ctx.lineTo(P[1][0], P[1][1]); ctx.lineTo(P[2][0], P[2][1]); ctx.closePath(); } ctx.stroke(); }
    d.parts.forEach((p, pi) => { const col = p.role === 'res' ? res : host, sel = main && D.sel === pi;
      if (p.kind === 'beam' || p.kind === 'wall') { const a = S(p.p[0]), b = S(p.p[1]); ctx.lineCap = 'butt'; ctx.strokeStyle = sel ? css('--ink') : col; ctx.lineWidth = Math.max(2, (p.kind === 'beam' ? p.b : p.t) * D.view.sc) + (sel ? 2 : 0); ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
        if (p.kind === 'wall') { ctx.strokeStyle = panel; ctx.lineWidth = 1; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke(); ctx.setLineDash([]); } }
      else if (p.kind === 'shape' && (p.use === 'frame' || p.use === 'walls')) { const poly = Drw.shapePoly(p).map(S); ctx.strokeStyle = sel ? css('--ink') : col; ctx.lineJoin = 'round'; ctx.lineWidth = Math.max(2, (p.use === 'frame' ? p.b : p.t) * D.view.sc) + (sel ? 2 : 0); ctx.beginPath(); poly.forEach((q, k) => k ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1])); ctx.closePath(); ctx.stroke();
        if (p.use === 'walls') { ctx.strokeStyle = panel; ctx.lineWidth = 1; ctx.setLineDash([3, 3]); ctx.stroke(); ctx.setLineDash([]); } }
      else if (p.kind === 'mass') { const s = S(p.p), r = Math.max(7, (p.side || 3) * D.view.sc); ctx.fillStyle = sel ? css('--ink') : col; ctx.fillRect(s[0] - r / 2, s[1] - r / 2, r, r); }
      else if (p.kind === 'spring') { const s = S(p.p); ctx.strokeStyle = sel ? css('--ink') : res; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(s[0], s[1], 8, 0, TAU); ctx.stroke(); ctx.fillStyle = res; ctx.beginPath(); ctx.arc(s[0], s[1], 3.5, 0, TAU); ctx.fill(); } });
  }
  ctx.globalAlpha = 1;
  // unit cell
  ctx.save(); ctx.setLineDash([6, 4]); ctx.strokeStyle = css('--ink'); ctx.lineWidth = 1.3; ctx.beginPath(); [[0, 0], d.a1, [d.a1[0] + d.a2[0], d.a1[1] + d.a2[1]], d.a2].forEach((p, k) => { const s = toS(p); k ? ctx.lineTo(s[0], s[1]) : ctx.moveTo(s[0], s[1]); }); ctx.closePath(); ctx.stroke(); ctx.restore();
  // selected part handles
  if (D.sel != null && d.parts[D.sel]) { ctx.fillStyle = panel; ctx.strokeStyle = css('--ink'); ctx.lineWidth = 1.5; const hs = handlesOf(d.parts[D.sel]); hs.forEach((v, k) => { const s = toS(v); if (d.parts[D.sel].kind === 'shape' && k === 1) { ctx.beginPath(); ctx.arc(s[0], s[1], 5, 0, TAU); ctx.fill(); ctx.stroke(); } else { ctx.fillRect(s[0] - 4, s[1] - 4, 8, 8); ctx.strokeRect(s[0] - 4, s[1] - 4, 8, 8); } }); }
  // pending polyline / rubber band
  if (D.pending) { const pts = D.pending.pts.concat(D.cur ? [D.cur.p] : []); ctx.strokeStyle = css('--gap'); ctx.lineWidth = 2; ctx.setLineDash([5, 3]); ctx.beginPath(); pts.forEach((q, k) => { const s = toS(q); k ? ctx.lineTo(s[0], s[1]) : ctx.moveTo(s[0], s[1]); }); ctx.stroke(); ctx.setLineDash([]); for (const q of D.pending.pts) { const s = toS(q); ctx.fillStyle = css('--gap'); ctx.beginPath(); ctx.arc(s[0], s[1], 3.5, 0, TAU); ctx.fill(); } }
  if (D.cur && D.tool !== 'select') { const s = toS(D.cur.p); ctx.strokeStyle = D.cur.kind === 'vertex' ? css('--danger') : D.cur.kind === 'edge' ? css('--gap') : ink2; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(s[0], s[1], 6, 0, TAU); ctx.stroke();
    ctx.fillStyle = ink2; ctx.font = '12px system-ui, sans-serif'; ctx.fillText(`${D.cur.p[0].toFixed(2)}, ${D.cur.p[1].toFixed(2)} mm${D.cur.kind === 'vertex' ? ' (point)' : D.cur.kind === 'edge' ? ' (on edge)' : ''}`, s[0] + 10, s[1] - 10); }
  ctx.fillStyle = ink2; ctx.font = '12.5px system-ui, sans-serif'; ctx.fillText(hintText(), 10, h - 10);
}
function hintText() {
  const t = D.tool;
  if (t === 'select') return 'Click a part to edit it. Drag a part to move it, a square handle to move a point, a round handle to resize or turn a shape. Drag the background to pan, wheel to zoom, Delete removes the part.';
  if (t === 'shape') return 'Click where the centre of the shape goes (or type the centre in the panel and press Add). Then drag its handles or edit its numbers.';
  if (t === 'beam' || t === 'wall') return D.pending ? 'Click the end point (keeps drawing from there). Esc or right-click to stop.' : 'Click the start point. Points snap to existing points (red), to edges (yellow) and to the grid.';
  if (t === 'plate' || t === 'hole') return D.pending ? `Click the next corner; click the first corner or press Enter to close the ${t === 'hole' ? 'hole' : 'plate'}.` : `Click the first corner of the ${t === 'hole' ? 'hole (inside a plate)' : 'plate outline'}.`;
  if (t === 'template') return 'Click the point where the resonator attaches (on a beam, wall, plate edge or point).';
  return 'Click to place it on a beam, wall or plate.';
}

// ---------- mouse ----------
let drag = null;
function bindEditor() {
  const cv = $('#cvDraw');
  const W = e => { const r = cv.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
  cv.addEventListener('contextmenu', e => { e.preventDefault(); finishChain(); });
  cv.addEventListener('wheel', e => { e.preventDefault(); const [x, y] = W(e), w0 = toW(x, y), f = Math.exp(-e.deltaY * 0.0015); D.view.sc *= f; const w1 = toW(x, y); D.view.cx += w0[0] - w1[0]; D.view.cy += w0[1] - w1[1]; D.view.keep = true; drawEditor(); }, { passive: false });
  cv.addEventListener('pointerdown', e => {
    if (e.button === 2) return; const [x, y] = W(e), w = toW(x, y); cv.setPointerCapture(e.pointerId);
    if (D.tool === 'select') {
      if (D.sel != null) { const p = D.dr.parts[D.sel]; const vs = handlesOf(p); for (let k = 0; k < vs.length; k++) { const s = toS(vs[k]); if (Math.hypot(s[0] - x, s[1] - y) < 8) { commit(); drag = { kind: 'vertex', k, pi: D.sel }; return; } } }
      const h = hit(w); if (h) { D.sel = h.pi; renderProps(); drawEditor(); drag = { kind: 'move', pi: h.pi, w0: w, orig: JSON.stringify(D.dr.parts[h.pi]), moved: false }; return; }
      D.sel = null; renderProps(); drag = { kind: 'pan', x, y, cx: D.view.cx, cy: D.view.cy }; drawEditor(); return;
    }
    const sp = snap(w); onClick(sp.p, e);
  });
  cv.addEventListener('pointermove', e => {
    const [x, y] = W(e), w = toW(x, y);
    if (drag && drag.kind === 'pan') { D.view.cx = drag.cx - (x - drag.x) / D.view.sc; D.view.cy = drag.cy + (y - drag.y) / D.view.sc; D.view.keep = true; drawEditor(); return; }
    if (drag && drag.kind === 'vertex') { const p = D.dr.parts[drag.pi], sp = snap(w, { pi: drag.pi }).p; setHandle(p, drag.k, sp); D.md = null; renderPropsLive(); drawEditor(); return; }
    if (drag && drag.kind === 'move') { const g = D.dr.grid > 0 ? D.dr.grid : 0; let dx = w[0] - drag.w0[0], dy = w[1] - drag.w0[1]; if (g) { dx = Math.round(dx / g) * g; dy = Math.round(dy / g) * g; }
      if (!drag.moved && Math.hypot(dx, dy) * D.view.sc < 3) return; if (!drag.moved) { commit(); drag.moved = true; }
      const q = JSON.parse(drag.orig); Drw.transformPart(q, ([a, b]) => [a + dx, b + dy]); D.dr.parts[drag.pi] = q; D.md = null; renderPropsLive(); drawEditor(); return; }
    D.cur = D.tool === 'select' ? null : snap(w); drawEditor();
  });
  cv.addEventListener('pointerup', () => { if (drag && (drag.kind === 'vertex' || (drag.kind === 'move' && drag.moved))) changed(); drag = null; });
  cv.addEventListener('pointerleave', () => { D.cur = null; drawEditor(); });
  cv.addEventListener('dblclick', () => { if (D.pending && (D.tool === 'plate' || D.tool === 'hole')) closePoly(); else finishChain(); });
  document.addEventListener('keydown', e => {
    if ($('#p-drw').classList.contains('hidden') || /INPUT|TEXTAREA|SELECT/.test(document.activeElement.tagName)) return;
    if (e.key === 'Escape') { finishChain(); D.pending = null; drawEditor(); }
    else if (e.key === 'Enter' && D.pending && (D.tool === 'plate' || D.tool === 'hole')) closePoly();
    else if ((e.key === 'Delete' || e.key === 'Backspace') && D.sel != null) { e.preventDefault(); delSel(); }
    else if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'z') { e.preventDefault(); e.shiftKey ? redo() : undo(); }
  });
}
function setVertex(p, k, q) {
  if (p.kind === 'beam') p.p[k] = [q[0], q[1], p.p[k][2] || 0];
  else if (p.kind === 'wall') p.p[k] = [q[0], q[1]];
  else if (p.kind === 'plate') { const rings = [p.poly, ...(p.holes || [])]; let c = k; for (const r of rings) { if (c < r.length) { r[c] = [q[0], q[1]]; return; } c -= r.length; } }
  else p.p = [q[0], q[1], p.p[2] || 0];
}
function onClick(p, e) {
  const t = D.tool, d = D.dr, role = D.role;
  if (t === 'beam' || t === 'wall') {
    if (!D.pending) { D.pending = { pts: [p] }; drawEditor(); return; }
    const a = D.pending.pts[0]; if (Math.hypot(p[0] - a[0], p[1] - a[1]) < 1e-6) return;
    commit();
    if (t === 'beam') { const z = num('dPz') || 0; d.parts.push({ kind: 'beam', p: [[a[0], a[1], z], [p[0], p[1], z]], b: num('dPb'), h: num('dPh'), role }); }
    else d.parts.push({ kind: 'wall', p: [a, p], t: num('dPt'), H: num('dPH'), z0: num('dPz0'), role });
    D.pending = { pts: [p] }; changed(); return;
  }
  if (t === 'plate' || t === 'hole') {
    if (D.pending && D.pending.pts.length >= 3) { const f = D.pending.pts[0]; if (Math.hypot(p[0] - f[0], p[1] - f[1]) < 1e-6) { closePoly(); return; } }
    if (!D.pending) D.pending = { pts: [] }; D.pending.pts.push(p); drawEditor(); return;
  }
  if (t === 'mass') { commit(); d.parts.push({ kind: 'mass', p: [p[0], p[1], num('dPmz') || 0], m: num('dPm'), side: num('dPside'), role }); changed(); return; }
  if (t === 'spring') { commit(); d.parts.push({ kind: 'spring', p: [p[0], p[1], num('dPsz') || 0], m: num('dPsm'), f: num('dPsf'), dir: [0, 0, 1], role: 'res' }); changed(); return; }
  if (t === 'template') { commit(); addTemplate(p); changed(); }
  if (t === 'shape') addShape(p);
}
function finishChain() { if (D.pending && (D.tool === 'beam' || D.tool === 'wall')) D.pending = null; drawEditor(); }
function closePoly() {
  const pts = D.pending.pts; D.pending = null;
  if (pts.length < 3) { B.setStatus('A plate or hole needs at least three corners.'); drawEditor(); return; }
  commit();
  if (D.tool === 'plate') D.dr.parts.push({ kind: 'plate', poly: pts, holes: [], z: num('dPpz'), t: num('dPpt'), role: D.role });
  else { const c = [pts.reduce((s, q) => s + q[0], 0) / pts.length, pts.reduce((s, q) => s + q[1], 0) / pts.length];
    const pl = D.dr.parts.findIndex(p => p.kind === 'plate' && Drw.pointInPoly(c, p.poly));
    if (pl < 0) { D.hist.pop(); B.setStatus('Draw the hole inside a plate.'); drawEditor(); return; }
    D.dr.parts[pl].holes = D.dr.parts[pl].holes || []; D.dr.parts[pl].holes.push(pts); }
  changed();
}
function delSel() { if (D.sel == null) return; commit(); D.dr.parts.splice(D.sel, 1); D.sel = null; changed(); }
function addTemplate(p) {
  const ty = $('#dTt').value, ang = num('dTa') * Math.PI / 180, z = num('dTz') || 0, d = [Math.cos(ang), Math.sin(ang)];
  const L1 = num('dTL'), b = num('dTb'), h = num('dTh'), bL = num('dTbL'), bW = num('dTbW'), bH = num('dTbH'), parts = D.dr.parts;
  const P = (q, s) => [+(q[0] + d[0] * s).toFixed(6), +(q[1] + d[1] * s).toFixed(6)];
  if (ty === 'spring') { parts.push({ kind: 'spring', p: [p[0], p[1], z], m: num('dTm'), f: num('dTf'), dir: [0, 0, 1], role: 'res' }); return; }
  let start = [p[0], p[1]];
  const stub = num('dTs') || 0;
  if (stub > 0) { const e = P(start, stub); parts.push({ kind: 'beam', p: [[start[0], start[1], z], [e[0], e[1], z]], b, h: num('dTsh') || bH, role: 'host' }); start = e; }
  if (ty === 'cant') { const L = Math.max(1e-3, L1 - bL / 2), e = P(start, L); parts.push({ kind: 'beam', p: [[start[0], start[1], z], [e[0], e[1], z]], b, h, role: 'res' }); const f = P(e, bL); parts.push({ kind: 'beam', p: [[e[0], e[1], z], [f[0], f[1], z]], b: bW, h: bH, role: 'res' }); }
  if (ty === 'L') { const e = P(start, L1), turn = num('dTturn') * Math.PI / 180, d2 = [Math.cos(ang + turn), Math.sin(ang + turn)], L2 = num('dTL2');
    const f = [+(e[0] + d2[0] * L2).toFixed(6), +(e[1] + d2[1] * L2).toFixed(6)], g = [+(f[0] + d2[0] * bL).toFixed(6), +(f[1] + d2[1] * bL).toFixed(6)];
    parts.push({ kind: 'beam', p: [[start[0], start[1], z], [e[0], e[1], z]], b, h, role: 'res' }, { kind: 'beam', p: [[e[0], e[1], z], [f[0], f[1], z]], b, h, role: 'res' }, { kind: 'beam', p: [[f[0], f[1], z], [g[0], g[1], z]], b: bW, h: bH, role: 'res' }); }
}

// ---------- properties panel: every number of the selected part, plus move / turn / scale / mirror / copy ----------
const MATF = [['E_GPa', 'Own E (blank = from material)', 'GPa'], ['nu', 'Own ν (blank = from material)', ''], ['rho', 'Own ρ (blank = from material)', 'kg/m³']];
function fieldsOf(p) {
  if (p.kind === 'beam') return [['x1', 'Start x₁', 'mm'], ['y1', 'Start y₁', 'mm'], ['x2', 'End x₂', 'mm'], ['y2', 'End y₂', 'mm'], ['L', 'Length (keeps start and direction)', 'mm'], ['b', 'Width b (in plane)', 'mm'], ['h', 'Height h (out of plane)', 'mm'], ['z', 'Height z (where it joins other parts)', 'mm'], ...MATF];
  if (p.kind === 'wall') return [['x1', 'Start x₁', 'mm'], ['y1', 'Start y₁', 'mm'], ['x2', 'End x₂', 'mm'], ['y2', 'End y₂', 'mm'], ['L', 'Length (keeps start and direction)', 'mm'], ['t', 'Thickness t', 'mm'], ['H', 'Height H', 'mm'], ['z0', 'Bottom edge z₀', 'mm'], ...MATF];
  if (p.kind === 'plate') return [['t', 'Thickness t', 'mm'], ['z', 'Height z (where it joins other parts)', 'mm'], ...MATF];
  if (p.kind === 'mass') return [['x', 'x', 'mm'], ['y', 'y', 'mm'], ['m', 'Mass', 'g'], ['side', 'Cube side (rotary inertia)', 'mm'], ['z', 'Height z', 'mm']];
  if (p.kind === 'spring') return [['x', 'x', 'mm'], ['y', 'y', 'mm'], ['m', 'Mass', 'g'], ['f', 'Frequency', 'Hz'], ['z', 'Height z', 'mm']];
  const f = [['cx', 'Centre x', 'mm'], ['cy', 'Centre y', 'mm']];
  if (p.type === 'ngon') f.push(['n', 'Number of sides', ''], ['R', 'Size: corner radius R', 'mm'], ['s', 'Size: side length', 'mm'], ['ri', 'Size: inner radius (to side)', 'mm']);
  if (p.type === 'circle') f.push(['R', 'Radius R', 'mm'], ['nseg', 'Segments on the circle', '']);
  if (p.type === 'rect') f.push(['sx', 'Width (before rotation)', 'mm'], ['sy', 'Height in plane (before rotation)', 'mm']);
  f.push(['ang', 'Rotation', 'deg']);
  if (p.use === 'plate') f.push(['t', 'Plate thickness t', 'mm'], ['z', 'Height z (where it joins other parts)', 'mm']);
  if (p.use === 'frame') f.push(['b', 'Beam width b', 'mm'], ['h', 'Beam height h', 'mm'], ['z', 'Height z (where it joins other parts)', 'mm']);
  if (p.use === 'walls') f.push(['t', 'Wall thickness t', 'mm'], ['H', 'Wall height H', 'mm'], ['z0', 'Bottom edge z₀', 'mm']);
  if (p.use !== 'hole') f.push(...MATF);
  return f;
}
const r6 = v => +(+v).toFixed(6);
function getF(p, k) {
  if (k === 'x1') return p.p[0][0]; if (k === 'y1') return p.p[0][1]; if (k === 'x2') return p.p[1][0]; if (k === 'y2') return p.p[1][1];
  if (k === 'L') return r6(Math.hypot(p.p[1][0] - p.p[0][0], p.p[1][1] - p.p[0][1]));
  if (k === 'x') return p.p[0]; if (k === 'y') return p.p[1];
  if (k === 'cx') return p.c[0]; if (k === 'cy') return p.c[1];
  if (k === 's') return r6(2 * p.R * Math.sin(Math.PI / p.n)); if (k === 'ri') return r6(p.R * Math.cos(Math.PI / p.n));
  if (k === 'ang') return p.ang || 0; if (k === 'nseg') return p.nseg || 24;
  if (k === 'z') return p.kind === 'beam' ? p.p[0][2] || 0 : p.kind === 'plate' || p.kind === 'shape' ? p.z || 0 : (p.p[2] || 0);
  return p[k];
}
function setF(p, k, v) {
  const ok = isFinite(v);
  if (['x1', 'y1', 'x2', 'y2'].includes(k)) { if (!ok) return; const i = k.endsWith('1') ? 0 : 1, c = k[0] === 'x' ? 0 : 1; p.p[i][c] = v; return; }
  if (k === 'L') { if (!(v > 0)) return; const dx = p.p[1][0] - p.p[0][0], dy = p.p[1][1] - p.p[0][1], L = Math.hypot(dx, dy) || 1; p.p[1][0] = r6(p.p[0][0] + dx / L * v); p.p[1][1] = r6(p.p[0][1] + dy / L * v); return; }
  if (k === 'x' || k === 'y') { if (ok) p.p[k === 'x' ? 0 : 1] = v; return; }
  if (k === 'cx' || k === 'cy') { if (ok) p.c[k === 'cx' ? 0 : 1] = v; return; }
  if (k === 's') { if (v > 0) p.R = r6(v / (2 * Math.sin(Math.PI / p.n))); return; }
  if (k === 'ri') { if (v > 0) p.R = r6(v / Math.cos(Math.PI / p.n)); return; }
  if (k === 'n') { if (v >= 3) p.n = Math.round(v); return; }
  if (k === 'z') { if (!ok) return; if (p.kind === 'beam') { p.p[0][2] = v; p.p[1][2] = v; } else if (p.kind === 'plate' || p.kind === 'shape') p.z = v; else p.p[2] = v; return; }
  if ((k === 'E_GPa' || k === 'rho') && !(v > 0)) { delete p[k]; return; }
  if (k === 'nu') { if (!(v > -1 && v < 0.5)) delete p[k]; else p[k] = v; return; }
  if (ok) p[k] = v;
}
const KIND_NAME = { beam: 'Beam', wall: 'Shell wall', plate: 'Free plate', mass: 'Point mass', spring: 'Mass-spring', shape: 'Shape' };
const SHAPE_NAME = { ngon: 'regular polygon', rect: 'rectangle', circle: 'circle' }, USE_NAME = { plate: 'shell plate', frame: 'frame of beams', walls: 'shell walls', hole: 'hole in a plate' };
function partName(p, i) { return `${KIND_NAME[p.kind]} ${i + 1}${p.kind === 'shape' ? ` (${p.type === 'ngon' ? p.n + '-sided polygon' : SHAPE_NAME[p.type]}, ${USE_NAME[p.use]})` : ''}`; }
const inp = (id, lab, val, u, attr = '') => `<div class="row"><label for="${id}">${lab}</label><div class="inp"><input id="${id}" type="number" step="any" value="${val ?? ''}" ${attr}><span>${u}</span></div></div>`;
function renderProps() {
  const el = $('#dProps'), p = D.sel != null ? D.dr.parts[D.sel] : null;
  if (!p) { el.innerHTML = '<p class="note">Nothing selected. Pick the Select tool and click a part to edit its numbers, or drag it to move it.</p>'; return; }
  const geo = p.kind === 'plate' ? `${p.poly.length} corners, ${(p.holes || []).length} hole${(p.holes || []).length === 1 ? '' : 's'}` : '';
  el.innerHTML = `<h3>${partName(p, D.sel)}</h3>${geo ? `<p class="note">${geo}</p>` : ''}` +
    (p.kind === 'shape' ? `<div class="row"><label for="dF_use">Use as</label><div class="inp"><select id="dF_use">${Object.entries(USE_NAME).map(([k, v]) => `<option value="${k}"${p.use === k ? ' selected' : ''}>${v[0].toUpperCase() + v.slice(1)}</option>`).join('')}</select></div></div>` : '') +
    fieldsOf(p).map(([k, lab, u]) => inp('dF_' + k, lab, getF(p, k), u, `data-f="${k}"`)).join('') +
    (p.kind === 'plate' ? `<div class="row"><label for="dF_poly">Corners (x, y per line)</label></div><textarea id="dF_poly" class="hjson" style="min-height:90px" spellcheck="false">${p.poly.map(q => q.join(', ')).join('\n')}</textarea>` : '') +
    ((p.kind === 'beam' || p.kind === 'plate' || (p.kind === 'shape' && (p.use === 'plate' || p.use === 'frame'))) ? `<div class="row"><label for="dF_zref">Thickness measured</label><div class="inp"><select id="dF_zref"><option value=""${!p.zref ? ' selected' : ''}>Centred on z</option><option value="bottom"${p.zref === 'bottom' ? ' selected' : ''}>Upward from z (sits on z)</option><option value="top"${p.zref === 'top' ? ' selected' : ''}>Downward from z (hangs from z)</option></select></div></div>` : '') +
    (['beam', 'wall', 'plate', 'shape'].includes(p.kind) && p.use !== 'hole' ? `<div class="row"><label for="dF_mat">Material</label><div class="inp"><select id="dF_mat"><option value=""${!p.mat ? ' selected' : ''}>From role (host or resonator)</option><option value="host"${p.mat === 'host' ? ' selected' : ''}>Host</option><option value="res"${p.mat === 'res' ? ' selected' : ''}>Resonator</option><option value="mem"${p.mat === 'mem' ? ' selected' : ''}>Membrane</option></select></div></div>` : '') +
    (p.kind !== 'spring' && p.use !== 'hole' ? `<div class="row"><label for="dF_role">Role</label><div class="inp"><select id="dF_role"><option value="host"${p.role !== 'res' ? ' selected' : ''}>Host (lattice)</option><option value="res"${p.role === 'res' ? ' selected' : ''}>Resonator</option></select></div></div>` : '') +
    `<h3 style="margin-top:12px">Move, turn, scale, copy</h3>
    <div class="row"><label>Move by dx, dy</label><div class="inp2"><input id="dX_dx" type="number" step="any" value="0"><input id="dX_dy" type="number" step="any" value="0"><span>mm</span></div></div>
    <div class="row"><label for="dX_rot">Turn by (about its centre)</label><div class="inp"><input id="dX_rot" type="number" step="any" value="60"><span>deg</span></div></div>
    <div class="row"><label for="dX_sc">Scale by (about its centre)</label><div class="inp"><input id="dX_sc" type="number" step="any" value="1.1"><span>×</span></div></div>
    <div class="row"><label>Copies around point x, y</label><div class="inp2"><input id="dX_px" type="number" step="any" value="0"><input id="dX_py" type="number" step="any" value="0"><span>mm</span></div></div>
    <div class="row"><label for="dX_n">Number of copies in the circle</label><div class="inp"><input id="dX_n" type="number" step="1" value="6"></div></div>
    <div class="xbtns"><button class="mini" type="button" data-x="move">Move</button><button class="mini" type="button" data-x="rot">Turn</button><button class="mini" type="button" data-x="sc">Scale</button><button class="mini" type="button" data-x="mx">Mirror ↔</button><button class="mini" type="button" data-x="my">Mirror ↕</button><button class="mini" type="button" data-x="dup">Duplicate</button><button class="mini" type="button" data-x="polar">Copy around point</button></div>` +
    `<div class="xbtns">${p.kind === 'plate' && (p.holes || []).length ? `<button class="mini" type="button" id="dDelHole">Remove last hole</button>` : ''}<button class="mini" type="button" id="dDel">Delete part</button></div>`;
  el.querySelectorAll('[data-f]').forEach(q => q.addEventListener('change', () => { commit(); setF(p, q.dataset.f, parseFloat(q.value)); changed(); }));
  const r = el.querySelector('#dF_role'); if (r) r.addEventListener('change', () => { commit(); p.role = r.value; changed(); });
  const mt = el.querySelector('#dF_mat'); if (mt) mt.addEventListener('change', () => { commit(); if (mt.value) p.mat = mt.value; else delete p.mat; changed(); });
  const zr = el.querySelector('#dF_zref'); if (zr) zr.addEventListener('change', () => { commit(); if (zr.value) p.zref = zr.value; else delete p.zref; changed(); });
  const us = el.querySelector('#dF_use'); if (us) us.addEventListener('change', () => { commit(); p.use = us.value; shapeDefaults(p); changed(); });
  const pl = el.querySelector('#dF_poly'); if (pl) pl.addEventListener('change', () => { const pts = pl.value.split('\n').map(l => l.split(/[\s,;]+/).filter(Boolean).map(Number)).filter(q => q.length >= 2 && q.every(isFinite)).map(q => [q[0], q[1]]);
    if (pts.length < 3) { B.setStatus('A plate needs at least three corners, one “x, y” per line.'); return; } commit(); p.poly = pts; changed(); });
  el.querySelectorAll('[data-x]').forEach(b => b.addEventListener('click', () => xform(b.dataset.x)));
  el.querySelector('#dDel').addEventListener('click', delSel);
  const dh = el.querySelector('#dDelHole'); if (dh) dh.addEventListener('click', () => { commit(); p.holes.pop(); changed(); });
}
// while dragging: refresh the numbers without rebuilding the panel
function renderPropsLive() { const p = D.sel != null ? D.dr.parts[D.sel] : null; if (!p) return; document.querySelectorAll('#dProps [data-f]').forEach(q => { if (document.activeElement !== q) { const v = getF(p, q.dataset.f); q.value = v == null ? '' : r6(v); } }); }
function xform(op) {
  const p = D.dr.parts[D.sel]; if (!p) return; const c = centreOf(p), deg = Math.PI / 180;
  const rot = (o, a) => ([x, y]) => [o[0] + Math.cos(a) * (x - o[0]) - Math.sin(a) * (y - o[1]), o[1] + Math.sin(a) * (x - o[0]) + Math.cos(a) * (y - o[1])];
  commit();
  if (op === 'move') { const dx = num('dX_dx') || 0, dy = num('dX_dy') || 0; Drw.transformPart(p, ([x, y]) => [x + dx, y + dy]); }
  if (op === 'rot') { const a = (num('dX_rot') || 0) * deg; Drw.transformPart(p, rot(c, a)); if (p.kind === 'shape') p.ang = +(((p.ang || 0) + num('dX_rot')) % 360).toFixed(4); }
  if (op === 'sc') { const f = num('dX_sc'); if (!(f > 0)) { D.hist.pop(); B.setStatus('Give a positive scale factor.'); return; }
    Drw.transformPart(p, ([x, y]) => [c[0] + f * (x - c[0]), c[1] + f * (y - c[1])]); if (p.kind === 'shape') { if (p.R) p.R = r6(p.R * f); if (p.sx) { p.sx = r6(p.sx * f); p.sy = r6(p.sy * f); } } }
  if (op === 'mx' || op === 'my') { Drw.transformPart(p, ([x, y]) => op === 'mx' ? [2 * c[0] - x, y] : [x, 2 * c[1] - y]);
    if (p.kind === 'shape') p.ang = +((op === 'mx' ? 180 - (p.ang || 0) : -(p.ang || 0)) % 360).toFixed(4); if (p.kind === 'plate') { p.poly.reverse(); (p.holes || []).forEach(h => h.reverse()); } }
  if (op === 'dup') { const q = clone(p), g = Math.max(D.dr.grid || 0, 0.5) * 2; Drw.transformPart(q, ([x, y]) => [x + g, y + g]); D.dr.parts.push(q); D.sel = D.dr.parts.length - 1; }
  if (op === 'polar') { const n = Math.round(num('dX_n')), o = [num('dX_px') || 0, num('dX_py') || 0]; if (!(n >= 2 && n <= 36)) { D.hist.pop(); B.setStatus('Number of copies: 2 to 36 (the part itself counts as one).'); return; }
    for (let k = 1; k < n; k++) { const q = clone(p), a = 360 / n * k; Drw.transformPart(q, rot(o, a * deg)); if (q.kind === 'shape') q.ang = +(((q.ang || 0) + a) % 360).toFixed(4); D.dr.parts.push(q); } }
  changed();
}
function shapeDefaults(p) {
  if (p.use === 'plate' && !(p.t > 0)) p.t = num('dSt') || 1;
  if (p.use === 'frame') { if (!(p.b > 0)) p.b = num('dSb') || 1; if (!(p.h > 0)) p.h = num('dSh') || 3; }
  if (p.use === 'walls') { if (!(p.t > 0)) p.t = num('dSt') || 1; if (!(p.H > 0)) p.H = num('dSH') || 3; if (p.z0 == null) p.z0 = num('dSz0') || 0; }
}
function newShape(c) {
  const ty = $('#dSty').value, mode = $('#dSm').value, v = num('dSv'), n = Math.round(num('dSn')) || 6, use = $('#dSu').value;
  const p = { kind: 'shape', type: ty, c: [r6(c[0]), r6(c[1])], ang: num('dSa') || 0, use, role: use === 'hole' ? 'host' : D.role };
  if (ty === 'ngon') { p.n = Math.max(3, n); p.R = mode === 's' ? r6(v / (2 * Math.sin(Math.PI / p.n))) : mode === 'ri' ? r6(v / Math.cos(Math.PI / p.n)) : v; }
  if (ty === 'circle') { p.R = mode === 'd' ? v / 2 : v; p.nseg = 24; }
  if (ty === 'rect') { p.sx = num('dSx'); p.sy = num('dSy'); }
  if (use === 'plate') { p.t = num('dSt'); p.z = num('dSz') || 0; }
  if (use === 'frame') { p.b = num('dSb'); p.h = num('dSh'); p.z = num('dSz') || 0; }
  if (use === 'walls') { p.t = num('dSt'); p.H = num('dSH'); p.z0 = num('dSz0') || 0; }
  if (use === 'plate' && $('#dSmem').checked) p.mat = 'mem';
  if ((use === 'plate' || use === 'frame') && $('#dSzr').value) p.zref = $('#dSzr').value;
  return p;
}
function addShape(c) {
  const p = newShape(c);
  try { Drw.expandParts([p], D.dr); } catch (e) { if (!/not inside any plate/.test(e.message)) { B.setStatus(e.message.replace('Shape 1', 'The new shape')); return; } }
  commit(); D.dr.parts.push(p); D.sel = D.dr.parts.length - 1;
  // the first frame shape of a cell sets the lattice: copies meet side to side, and keep doing so when it is resized
  const others = D.dr.parts.some(q => q !== p && q.role !== 'res' && q.kind !== 'spring' && q.use !== 'hole');
  if (!others && p.use !== 'hole' && !latticeFromOutline(outlineOf(p)).err) { p.role = 'host'; D.dr.latFollow = true; fitLattice(true, p); syncLat(); changed(); setTool('select'); B.setStatus(`Frame added (as host) and the lattice fitted to it (a₁ = (${D.dr.a1.join(', ')}), a₂ = (${D.dr.a2.join(', ')}) mm); it follows the frame when you resize it. Untick “Keep fitting” in the Lattice box to set a₁ and a₂ yourself.`); return; }
  changed(); setTool('select');
  B.setStatus('Shape added. Drag it to move it, drag the round handle to resize or turn it, or type its numbers on the right.');
}
function setTool(t) {
  D.tool = t; D.pending = null; document.querySelectorAll('[data-dtool]').forEach(b => b.setAttribute('aria-checked', b.dataset.dtool === t));
  document.querySelectorAll('[data-dnew]').forEach(b => b.classList.toggle('hidden', !b.dataset.dnew.split(' ').includes(t)));
  drawEditor();
}

// ---------- model (built lazily) and 3D view ----------
let mT = null;
function scheduleModel() { clearTimeout(mT); mT = setTimeout(() => { buildModel(); draw3D(); drawEditor(); }, 250); }
function readSolve() { const P = B.readParams(); return { P, mat: P.mat, beam: P.beam, kper: Math.round(num('kper')), nb: Math.max(8, Math.round(num('nb'))), fmax: num('fmax') || 0, fmin: Math.max(0, num('fmin') || 0), nfc: Math.min(150, Math.max(10, Math.round(num('nf')))), path: $('#path').value, pathErr: P.kpathErr }; }
function buildModel() {
  const s = readSolve(), key = JSON.stringify([D.dr, s.mat, s.beam]);
  if (D.md && D.mdKey === key) return D.md;
  try { D.md = Drw.build(D.dr, s.mat, { beam: s.beam }); D.mdKey = key; D.mdErr = null; } catch (e) { D.md = null; D.mdErr = e.message; }
  renderInfo(); return D.md;
}
function renderInfo() {
  const el = $('#dInfo'); if (!el) return;
  if (D.mdErr) { el.innerHTML = `<tr><td class="err" colspan="2">${D.mdErr}</td></tr>`; return; }
  const md = D.md; if (!md) { el.innerHTML = ''; return; }
  const nb = md.elems.filter(e => e.type === 'beam').length, ns = md.elems.length - nb, mh = md.elems.filter(e => e.role === 'host').reduce((s, e) => s + e.mass, 0) + md.pmass.filter(p => p.role === 'host').reduce((s, p) => s + p.m, 0);
  const mr = md.elems.filter(e => e.role === 'res').reduce((s, e) => s + e.mass, 0) + md.pmass.filter(p => p.role === 'res').reduce((s, p) => s + p.m, 0) + md.springs.reduce((s, p) => s + p.m, 0);
  const rows = [['Parts', D.dr.parts.length], ['Beam / shell elements', `${nb} / ${ns}`], ['Nodes, degrees of freedom', `${md.N}, ${6 * md.N + md.springs.length}`], ['Cell area', `${(md.S * 1e6).toFixed(1)} mm²`], ['Host mass per cell', B.fmtM(mh)], ['Resonator mass per cell', mr ? B.fmtM(mr) : '–'], ['Mass ratio', mr ? (mr / mh).toFixed(3) : '–']];
  if (!mr && D.dr.parts.length) rows.push(['<span class="warn">No resonator parts</span>', 'every part is host (blue): resonator share 0 %, host material everywhere. Use “Frame = host, rest = resonator” or set Role per part']);
  if (6 * md.N > 1500) rows.push(['<span class="warn">Large model</span>', 'raise the mesh size for speed']);
  if (md.merged && md.merged.length) rows.push(['Overlapping beams counted once', `part${md.merged.length > 1 ? 's' : ''} ${md.merged.map(i => i + 1).join(', ')}`]);
  if (md.floating && md.floating.length) rows.push(['<span class="err">Not connected</span>', `part${md.floating.length > 1 ? 's' : ''} ${md.floating.map(i => i + 1).join(', ')}`]);
  el.innerHTML = rows.map(r => `<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`).join('');
}
function draw3D() {
  const cv = $('#cvDraw3d'); if (!cv || $('#p-drw').classList.contains('hidden')) return;
  const { ctx, w, h } = B.prep(cv), md = D.md || buildModel();
  if (!md) { ctx.fillStyle = css('--danger'); ctx.font = '600 14px system-ui, sans-serif'; B.wrapText(ctx, D.mdErr || 'Draw a cell.', 14, 28, w - 28, 18); return; }
  const cy = Math.cos(D.yaw), sy = Math.sin(D.yaw), ce = Math.cos(D.elev), se = Math.sin(D.elev);
  const f = p => { const x = cy * p[0] - sy * p[1], y = sy * p[0] + cy * p[1]; return [x, y * ce + p[2] * se, -y * se + p[2] * ce]; };
  const items = [], tiles = [[0, 0], [1, 0], [0, 1], [1, 1]];
  for (const [i, j] of tiles) { const sh = [i * md.a1[0] + j * md.a2[0], i * md.a1[1] + j * md.a2[1], 0], main = !i && !j;
    for (const e of md.elems) { const X = e.Xd.map(x => [x[0] + sh[0], x[1] + sh[1], x[2] + (e.dz || 0)]); items.push({ e, X, main }); }
    for (const pm of md.pmass) { const q = D.dr.parts[pm.part], p = q && q.p ? [q.p[0] / 1000, q.p[1] / 1000, (q.p[2] || 0) / 1000] : md.nodes[pm.node]; items.push({ mass: true, X: [[p[0] + sh[0], p[1] + sh[1], p[2]]], side: pm.side, role: pm.role, main }); }
    for (const sp of md.springs) { const q = D.dr.parts[sp.part], p = q && q.p ? [q.p[0] / 1000, q.p[1] / 1000, (q.p[2] || 0) / 1000] : md.nodes[sp.node]; items.push({ mass: true, X: [[p[0] + sh[0], p[1] + sh[1], p[2]]], side: 0.003, role: 'res', main, spring: true }); } }
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity; for (const it of items) for (const p of it.X) { const q = f(p); x0 = Math.min(x0, q[0]); x1 = Math.max(x1, q[0]); y0 = Math.min(y0, q[1]); y1 = Math.max(y1, q[1]); }
  const sc = Math.min((w - 40) / Math.max(1e-9, x1 - x0), (h - 50) / Math.max(1e-9, y1 - y0)), ox = w / 2 - sc * (x0 + x1) / 2, oy = (h - 20) / 2 + sc * (y0 + y1) / 2;
  const Pj = p => { const q = f(p); return [ox + sc * q[0], oy - sc * q[1], q[2]]; };
  items.forEach(it => { it.z = it.X.reduce((s, p) => s + Pj(p)[2], 0) / it.X.length; }); items.sort((a, b) => a.z - b.z);
  const host = css('--host'), res = css('--res'), panel = css('--panel'), L = unitV([0.3, -0.5, 0.8]);
  for (const it of items) {
    const base = (it.e ? it.e.role : it.role) === 'res' ? res : host, col = it.main ? base : B.mix(base, panel, 0.55);
    if (it.mass) { const q = Pj(it.X[0]), s = Math.max(6, Math.min(30, (it.side || 0.003) * sc)); ctx.fillStyle = col; if (it.spring) { ctx.beginPath(); ctx.arc(q[0], q[1], s / 2, 0, TAU); ctx.fill(); } else ctx.fillRect(q[0] - s / 2, q[1] - s / 2, s, s); continue; }
    if (it.e.type === 'beam') { const a = Pj(it.X[0]), b = Pj(it.X[1]); ctx.strokeStyle = col; ctx.lineWidth = Math.max(1.5, Math.min(24, it.e.draw.b * sc)); ctx.lineCap = 'round'; ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke(); }
    else { const P = it.X.map(Pj), nn = unitV(crossV(subV(it.X[1], it.X[0]), subV(it.X[2], it.X[0]))), sh = 0.55 + 0.45 * Math.abs(nn[0] * L[0] + nn[1] * L[1] + nn[2] * L[2]);
      ctx.fillStyle = B.mix(panel, col, 0.35 + 0.5 * sh); ctx.strokeStyle = B.mix(panel, col, 0.25 + 0.4 * sh); ctx.lineWidth = 0.6; ctx.beginPath(); ctx.moveTo(P[0][0], P[0][1]); ctx.lineTo(P[1][0], P[1][1]); ctx.lineTo(P[2][0], P[2][1]); ctx.closePath(); ctx.fill(); ctx.stroke(); }
  }
  ctx.fillStyle = css('--ink-2'); ctx.font = '12.5px system-ui, sans-serif'; ctx.fillText('2 × 2 cells, mesh as solved; drag to rotate', 10, h - 10);
}
const subV = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], crossV = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]], unitV = a => { const L = Math.hypot(...a) || 1; return a.map(v => v / L); };
function bind3D() {
  const cv = $('#cvDraw3d'); let dg = null;
  cv.addEventListener('pointerdown', e => { dg = { x: e.clientX, y: e.clientY, yaw: D.yaw, elev: D.elev }; cv.setPointerCapture(e.pointerId); });
  cv.addEventListener('pointermove', e => { if (!dg) return; D.yaw = dg.yaw - (e.clientX - dg.x) * 0.01; D.elev = Math.max(0, Math.min(Math.PI / 2, dg.elev + (e.clientY - dg.y) * 0.01)); draw3D(); });
  cv.addEventListener('pointerup', () => { dg = null; });
}

// ---------- solve: FEM (full Bloch FE of the drawn cell) and/or PWE (homogenised plate) ----------
async function run() {
  if (B.isRunning()) return false;
  const s = readSolve();
  if (D.genErr && isGenDesign()) { B.setStatus('Cannot build this cell: ' + D.genErr); return false; }
  if (!(s.kper >= 3 && s.kper <= 60)) { B.setStatus('k-points per segment must be between 3 and 60.'); return false; }
  if (!(s.nb >= 4 && s.nb <= 60)) { B.setStatus('Bands must be between 4 and 60.'); return false; }
  if (s.path === 'custom' && s.pathErr) { B.setStatus(s.pathErr); return false; }
  B.startRun(); const t0 = performance.now(); let ok = false;
  try {
    D.md = null; const md = buildModel(); if (!md) throw new Error(D.mdErr);
    const Lat = latFor(md, s.path), path = Eng.kPath(Lat, s.kper), ell = Math.max(Math.hypot(...md.a1), Math.hypot(...md.a2));
    const R = await Res.solve({ meth: 'fem', Lat, path, plate: null, res: [], md, nfc: s.nfc, nbFE: s.nb, fmin: s.fmin, fmax: s.fmax, ell }, B.setStatus);
    D.run = R; D.runKey = runKey(); D.sel = null; D.mode = null;
    const nF = R.fe.gaps.length, nz = R.fe.flex.length, ni = R.fe.ip.length;
    B.setStatus(`Drawn cell solved in ${((performance.now() - t0) / 1000).toFixed(1)} s along ${path.labels.join('–')} (${Lat.path.mode === 'std' ? 'standard path' : Lat.path.mode === 'custom' ? 'your path' : Lat.path.mode + ' zone'}): ${nz} out-of-plane gap${nz === 1 ? '' : 's'}, ${ni} in-plane, ${nF} complete. Click a band point to see its mode.`);
    if (!md.elems.some(e => e.role === 'res') && !md.pmass.some(q => q.role === 'res') && !md.springs.length) B.setStatus($('#status').textContent + ' Note: no part is marked as resonator, so every band shows 0 % resonator share and the host material is used everywhere. Use “Frame = host, rest = resonator” in the Unit cell tab.');
    if (md.floating && md.floating.length) B.setStatus($('#status').textContent + ` Warning: part${md.floating.length > 1 ? 's' : ''} ${md.floating.map(i => i + 1).join(', ')} ${md.floating.length > 1 ? 'are' : 'is'} not connected to the lattice and float${md.floating.length > 1 ? '' : 's'} freely; make ${md.floating.length > 1 ? 'them' : 'it'} touch a host part (points snap when they meet).`);
    renderResults(); drawBands(); drawnAssessment(); ok = true;
  } catch (e) { B.setStatus('Solver stopped: ' + e.message); if (!/stopped by the user/.test(e.message)) console.error(e); }
  B.endRun(); return ok;
}
function drawBands() {
  if ($('#p-disp').classList.contains('hidden') || $('#dispDrawn').classList.contains('hidden')) return;
  D.geoB = Res.bands($('#cvDrawBand'), D.run, D.show, D.sel); Res.complex($('#cvDrawCplx'), D.run, D.show);
  const R = D.run; $('#dBandSub').textContent = R ? `${6 * R.inp.md.N + R.inp.md.springs.length} degrees of freedom, ${R.fe.F[0].length} bands` : '';
  if (!D.mode) { const { ctx } = B.prep($('#cvDrawMode')); ctx.fillStyle = css('--ink-2'); ctx.font = '14px system-ui, sans-serif'; ctx.fillText(R ? 'Click a point on the band structure.' : 'Run the analysis, then click a band point.', 16, 30); }
}
function renderResults() { if (D.run) Res.table($('#dGaps'), D.run); }
function pickMode(ev) {
  const R = D.run, G = D.geoB; if (!R || !G) return;
  const r = ev.currentTarget.getBoundingClientRect(), p = Viz.pickPoint(G.pickSets, G.X, G.Y, ev.clientX - r.left, ev.clientY - r.top); if (!p) return;
  try { D.sel = p; D.mode = Res.mode(R, p, D.lines); const [sg, u] = R.inp.path.seg[p.i], L = R.inp.path.labels;
    const brn = Res.branchName(R, p, D.show);
    $('#dModeSub').textContent = `Band ${p.b + 1}${brn ? ' (' + brn + ')' : ''}, ${B.fmtF(D.mode.f)}, ${u < 1e-9 ? L[sg] : u > 1 - 1e-9 ? L[sg + 1] : `${Math.round(100 * u)} % from ${L[sg]} to ${L[sg + 1]}`}`; drawBands(); }
  catch (e) { B.setStatus(e.message); }
}
window.addEventListener('mmview', () => { if (D.mode) D.mode._arw = null; drawBands(); });
function drawMode(ph) { if (!D.mode) return; const amp = num('dAmp') || 1; if (D.mode.kind === 'fe') Viz.modeFE($('#cvDrawMode'), D.mode, ph, amp, D.mview); else Viz.modePWE($('#cvDrawMode'), D.mode, ph, amp); }


// ---------- cell generator: ready-made cells B, C and "build your own" (lattice + resonator), all editable afterwards ----------
const DESIGN_FIX = { hexspokes: { lat: 'hex', res: 'spokes' }, hexbeam: { lat: 'hex', res: 'across' } };
function readGen() {
  const des = $('#design').value, fx = DESIGN_FIX[des] || {}, v = id => num(id), sel = id => $('#' + id).value, txt = id => $('#' + id).value.trim();
  return { des, lat: fx.lat || sel('gLat'), res: fx.res || sel('gRes'), L: v('gL'), w: v('gW'), H: v('gH'),
    cw: txt('gCw'), ct: v('gCt'), cL: v('gCL'), cb: v('gCb'), ch: v('gCh'), bl: v('gBl'), bw: v('gBw'), bh: v('gBh'),
    ms: sel('gMs'), mR: v('gMR'), mt: v('gMt'), mrho: v('gMrho'), sb: v('gSb'), sh: v('gSh'), ax: sel('gAx'), rr: v('gRr'), rw: v('gRw'), ks: v('gKs'), kt: v('gKt'),
    mm: v('gMm'), pm: v('gPm'), pf: v('gPf'), pw: txt('gPw'), zr: sel('gZr') };
}
// frame of the lattice: lattice vectors, frame beams, and the polygon (corners, centre) that holds the resonator
function latticeOf(G) {
  const L = G.L, s3 = Math.sqrt(3), pt = (x, y) => [r6(x), r6(y)];
  if (G.lat === 'square') { const h = L / 2, P = [pt(h, -h), pt(h, h), pt(-h, h), pt(-h, -h)]; return { a1: [L, 0], a2: [0, L], c: [0, 0], P, frame: [[P[3], P[0]], [P[2], P[3]]] }; }
  if (G.lat === 'tri') { const P = [pt(0, 0), pt(L, 0), pt(L / 2, s3 / 2 * L)]; return { a1: [L, 0], a2: [r6(L / 2), r6(s3 / 2 * L)], c: [r6(L / 2), r6(s3 / 6 * L)], P, frame: [[P[0], P[1]], [P[0], P[2]], [P[1], P[2]]] }; }
  if (G.lat === 'kagome') { const c = [1.5 * L, s3 / 2 * L], P = [0, 1, 2, 3, 4, 5].map(k => pt(c[0] + L * Math.cos(k * Math.PI / 3), c[1] + L * Math.sin(k * Math.PI / 3)));
    const n1 = pt(L, 0), n2 = pt(L / 2, s3 / 2 * L); return { a1: [2 * L, 0], a2: [L, r6(s3 * L)], c: c.map(r6), P, frame: [[[0, 0], n1], [[0, 0], n2], [n1, n2], [n1, pt(2 * L, 0)], [n2, pt(L, s3 * L)], [n2, pt(0, s3 * L)]] }; }
  if (G.lat === 'free') { const a1 = D.dr ? D.dr.a1.slice() : [14, 0], a2 = D.dr ? D.dr.a2.slice() : [7, 12.124356]; return { a1, a2, c: [r6((a1[0] + a2[0]) / 2), r6((a1[1] + a2[1]) / 2)], P: null, frame: [] }; }
  const P = [0, 1, 2, 3, 4, 5].map(k => pt(L * Math.cos(k * Math.PI / 3), L * Math.sin(k * Math.PI / 3)));
  return { a1: [r6(1.5 * L), r6(s3 / 2 * L)], a2: [0, r6(s3 * L)], c: [0, 0], P, frame: [[P[0], P[1]], [P[1], P[2]], [P[2], P[3]]] };
}
const wallList = (t, n) => { t = (t || '').toLowerCase(); if (!t || t === 'all') return [...Array(n).keys()]; return [...new Set(t.split(/[^0-9]+/).filter(Boolean).map(Number).filter(k => k >= 1 && k <= n).map(k => k - 1))]; };
function generate(G) {
  const pos = (v, name) => { if (!(v > 0)) throw new Error(`${name} must be a positive number.`); };
  if (G.lat !== 'free') { pos(G.L, 'The lattice size'); pos(G.w, 'The wall width'); pos(G.H, 'The wall height'); }
  const Lt = latticeOf(G), parts = [], c = Lt.c, P = Lt.P;
  for (const [a, b] of Lt.frame) parts.push({ kind: 'beam', p: [[a[0], a[1], 0], [b[0], b[1], 0]], b: G.w, h: G.H, role: 'host', frame: true });
  if (G.res !== 'none' && !P) throw new Error('With a free lattice, draw the frame and the resonator in the Unit cell tab (choose “None” as resonator here).');
  const n = P ? P.length : 0, dir = k => { const d = [P[k][0] - c[0], P[k][1] - c[1]], l = Math.hypot(...d); return [d[0] / l, d[1] / l]; }, ang0 = P ? Math.atan2(P[0][1] - c[1], P[0][0] - c[0]) * 180 / Math.PI : 0;
  const zr = G.zr === 'bottom' ? 'bottom' : undefined, rho = G.mrho > 0 ? G.mrho : undefined;
  const centreMass = () => { pos(G.mR, 'The centre mass radius'); pos(G.mt, 'The centre mass thickness');
    const sh = G.ms === 'circle' ? { type: 'circle', R: G.mR, nseg: 24 } : { type: 'ngon', n, R: G.mR };
    return Object.assign({ kind: 'shape', c: c.slice(), ang: r6(ang0), use: 'plate', t: G.mt, z: 0, role: 'res', rho, zref: zr }, sh); };
  const inR = P ? Math.min(...P.map((q, k) => { const r = P[(k + 1) % n], m = [(q[0] + r[0]) / 2 - c[0], (q[1] + r[1]) / 2 - c[1]]; return Math.hypot(...m); })) : 0;
  if (G.res === 'cant') { pos(G.cL, 'The arm length'); pos(G.cb, 'The arm width'); pos(G.ch, 'The arm thickness'); pos(G.bl, 'The tip block length'); pos(G.bw, 'The tip block width'); pos(G.bh, 'The tip block height');
    if (!(G.ct > 0 && G.ct < 1)) throw new Error('The position along the wall must lie between 0 and 1.');
    const ws = wallList(G.cw, n); if (!ws.length) throw new Error(`Walls: give numbers from 1 to ${n}, or “all”.`);
    for (const k of ws) { const A = P[k], Bq = P[(k + 1) % n], q = [A[0] + (Bq[0] - A[0]) * G.ct, A[1] + (Bq[1] - A[1]) * G.ct]; let nx = -(Bq[1] - A[1]), ny = Bq[0] - A[0]; const l = Math.hypot(nx, ny); nx /= l; ny /= l;
      if ((c[0] - q[0]) * nx + (c[1] - q[1]) * ny < 0) { nx = -nx; ny = -ny; }
      const at = d => [r6(q[0] + nx * d), r6(q[1] + ny * d), 0], s0 = G.w / 2, e = s0 + Math.max(1e-3, G.cL - G.bl / 2);
      const reach = (c[0] - q[0]) * nx + (c[1] - q[1]) * ny; if (s0 + G.cL + G.bl / 2 > 0.95 * reach) throw new Error(`The arm with its tip block reaches the centre of the cell (it may touch the other arms). Keep arm length + half the block length below ${(0.95 * reach - s0).toFixed(2)} mm.`);
      parts.push({ kind: 'beam', p: [at(0), at(s0)], b: G.cb, h: G.H, role: 'host' }, { kind: 'beam', p: [at(s0), at(e)], b: G.cb, h: G.ch, role: 'res' }, { kind: 'beam', p: [at(e), at(e + G.bl)], b: G.bw, h: G.bh, role: 'res', rho, zref: zr }); } }
  if (G.res === 'spokes') { pos(G.sb, 'The spoke width'); pos(G.sh, 'The spoke thickness'); const m = centreMass(); parts.push(m); if (G.mR >= 0.9 * inR) throw new Error('The centre mass is too large for the cell.');
    const ring = G.rr > G.mR && G.rr < inR;
    for (let k = 0; k < n; k++) { const d = dir(k), a = G.ms === 'circle' ? [c[0] + d[0] * G.mR, c[1] + d[1] * G.mR] : Drw.shapePoly(m)[k]; parts.push({ kind: 'beam', p: [[r6(a[0]), r6(a[1]), 0], [P[k][0], P[k][1], 0]], b: G.sb, h: G.sh, role: 'res' }); }
    if (ring) parts.push({ kind: 'shape', type: 'ngon', n, R: G.rr, c: c.slice(), ang: r6(ang0), use: 'frame', b: G.rw > 0 ? G.rw : G.sb, h: G.sh, z: 0, role: 'res' });
    if (G.ks > 0) { if (!(G.kt > 0)) throw new Error('Give the block thickness.');
      for (let k = 0; k < n; k++) { const d = dir(k), Rk = Math.hypot(P[k][0] - c[0], P[k][1] - c[1]), rb = ring ? G.rr : (G.mR + Rk) / 2;
        parts.push({ kind: 'shape', type: 'ngon', n: 6, R: G.ks, c: [r6(c[0] + d[0] * rb), r6(c[1] + d[1] * rb)], ang: r6(ang0), use: 'plate', t: G.kt, z: 0, role: 'res', rho, zref: zr }); } } }
  if (G.res === 'across') { pos(G.sb, 'The beam width'); pos(G.sh, 'The beam thickness'); parts.push(centreMass());
    let ends; if (G.ax === 'corner' || n % 2) { const k = 0, j = n % 2 ? null : n / 2; ends = j == null ? [P[0], [(P[1][0] + P[2][0]) / 2, (P[1][1] + P[2][1]) / 2]] : [P[k], P[j]]; }
    else ends = [0, n / 2].map(k => [(P[k][0] + P[k + 1][0]) / 2, (P[k][1] + P[k + 1][1]) / 2]);
    for (const e of ends) parts.push({ kind: 'beam', p: [[c[0], c[1], 0], [r6(e[0]), r6(e[1]), 0]], b: G.sb, h: G.sh, role: 'res' }); }
  if (G.res === 'membrane') { pos(G.mm, 'The membrane thickness'); const m = centreMass(); parts.push({ kind: 'plate', poly: P.map(q => q.slice()), holes: [Drw.shapePoly(m)], z: 0, t: G.mm, role: 'res', mat: 'mem' }, m); }
  if (G.res === 'spring') { pos(G.pm, 'The resonator mass'); pos(G.pf, 'The resonator frequency'); const ws = wallList(G.pw, n); if (!ws.length) throw new Error(`Walls: give numbers from 1 to ${n}, or “all”.`);
    for (const k of ws) { const A = P[k], Bq = P[(k + 1) % n]; parts.push({ kind: 'spring', p: [r6((A[0] + Bq[0]) / 2), r6((A[1] + Bq[1]) / 2), 0], m: G.pm, f: G.pf, dir: [0, 0, 1], role: 'res' }); } }
  const L = G.lat === 'free' ? Math.hypot(...Lt.a1) : G.L;
  return { a1: Lt.a1, a2: Lt.a2, mesh: G.lat === 'free' && D.dr ? D.dr.mesh || 2 : r6(Math.max(0.2, L / 5)), grid: D.dr && D.dr.grid != null ? D.dr.grid : 0.5, parts, _gen: JSON.stringify(G), _edited: false };
}
const isGenDesign = () => ['hexspokes', 'hexbeam', 'own'].includes($('#design').value);
function regenerate(force) {
  if (!isGenDesign()) return;
  const G = readGen(), key = JSON.stringify(G);
  if (!force && D.dr && D.dr._gen === key) { D.genErr = null; $('#gErr').classList.add('hidden'); genNote(); return; }
  if (!force && D.dr && D.dr._edited && D.dr._design === G.des) { genNote(); return; }
  let dr; try { dr = generate(G); } catch (e) { D.genErr = e.message; $('#gErr').textContent = 'Cannot build this cell: ' + e.message; $('#gErr').classList.remove('hidden'); B.setStatus('Geometry: ' + e.message); return; }
  D.genErr = null; $('#gErr').classList.add('hidden');
  dr._design = G.des; if (D.dr) commit(); D.dr = dr; D.sel = null; D.view = null; syncLat(); D.dr._edited = false; changed(true); genNote();
}
function genNote() { const el = $('#gEdited'); if (!el || !D.dr) return; el.classList.toggle('hidden', !(D.dr._edited && D.dr._gen !== JSON.stringify(readGen()))); }
// ---------- save / load ----------
function download() { const blob = new Blob([JSON.stringify(D.dr, null, 1)], { type: 'application/json' }), a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'unit-cell-drawing.json'; document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500); }
function load(file) { const r = new FileReader(); r.onload = () => { try { const o = JSON.parse(r.result); if (!Array.isArray(o.a1) || !Array.isArray(o.parts)) throw new Error('not a drawing file'); commit(); D.dr = o; D.sel = null; D.view = null; syncLat(); changed(); B.setStatus('Drawing loaded.'); } catch (e) { B.setStatus('Could not read that file: ' + e.message); } }; r.readAsText(file); }

// ---------- zone path of a drawn cell: symmetry found from the drawing itself ----------
function latFor(md, mode, dr) {
  const Lat = Hom.latticeFor({ a1: md.a1, a2: md.a2 }, 'auto');
  Lat.sym = Drw.symmetry(dr || D.dr);
  Lat.path = Eng.makePath(Lat, mode || 'std');
  return Lat;
}
function lattice() { // for the zone-path preview in the Solver panel
  const d = D.dr, a1 = [d.a1[0] / 1000, d.a1[1] / 1000], a2 = [d.a2[0] / 1000, d.a2[1] / 1000];
  if (!(Math.abs(a1[0] * a2[1] - a1[1] * a2[0]) > 1e-14)) throw new Error('The lattice vectors are parallel.');
  return latFor({ a1, a2 }, $('#path').value);
}
const runKey = () => JSON.stringify([D.dr, B.readParams()]);
const FAMD = { oop: 'flex', ip: 'ip', all: 'gaps' }, FAMD_LABEL = { oop: 'Out-of-plane', ip: 'In-plane', all: 'All' };
function famRows(F, Z, fam) { if (fam === 'all') return F; const want = fam === 'oop', rows = F.map((r, i) => r.filter((f, b) => (Z[i][b] >= 0.5) === want)), n = Math.min(...rows.map(r => r.length)); return n > 0 ? rows.map(r => r.slice(0, n)) : null; }

// ---------- sensitivity sweep of a drawn cell ----------
const SWEEP_FIELDS = { beam: [['b', 'width b'], ['h', 'height h']], wall: [['t', 'thickness t'], ['H', 'height H']], plate: [['t', 'thickness t']], mass: [['m', 'mass']], spring: [['m', 'mass'], ['f', 'frequency']],
  shape: p => [['size', 'size (scaled about its centre)'], ...(p.use === 'plate' ? [['t', 'thickness t']] : p.use === 'frame' ? [['b', 'beam width b'], ['h', 'beam height h']] : p.use === 'walls' ? [['t', 'wall thickness t'], ['H', 'wall height H']] : [])] };
function sweepParts() {
  const el = $('#spPart'); if (!el) return; const cur = el.value, opts = [];
  D.dr.parts.forEach((p, i) => { const f = SWEEP_FIELDS[p.kind]; for (const [k, lab] of (typeof f === 'function' ? f(p) : f || [])) opts.push(`<option value="${i}:${k}">${partName(p, i)}: ${lab}</option>`); });
  el.innerHTML = opts.join('') || '<option value="">(draw some parts first)</option>'; if ([...el.options].some(o => o.value === cur)) el.value = cur;
}
function scaledDrawing(dr, key, fac, partKey) {
  const d = clone(dr);
  if (key === 'a') { const S = ([x, y]) => [x * fac, y * fac]; d.a1 = S(d.a1).map(r6); d.a2 = S(d.a2).map(r6); d.mesh = r6((d.mesh || 2) * fac);
    for (const p of d.parts) { Drw.transformPart(p, S); if (p.kind === 'beam') p.b = r6(p.b * fac); if (p.kind === 'wall') p.t = r6(p.t * fac); if (p.kind === 'shape') { if (p.R) p.R = r6(p.R * fac); if (p.sx) { p.sx = r6(p.sx * fac); p.sy = r6(p.sy * fac); } if (p.use === 'frame') p.b = r6(p.b * fac); if (p.use === 'walls') p.t = r6(p.t * fac); } } }
  if (key === 'part') { const [i, k] = partKey.split(':'), p = d.parts[+i]; if (!p) throw new Error('Pick a part to vary.');
    if (k === 'size') { const c = p.c.slice(); Drw.transformPart(p, ([x, y]) => [c[0] + fac * (x - c[0]), c[1] + fac * (y - c[1])]); if (p.R) p.R = r6(p.R * fac); if (p.sx) { p.sx = r6(p.sx * fac); p.sy = r6(p.sy * fac); } }
    else p[k] = r6(p[k] * fac); }
  return d;
}
async function runSweep() {
  if (B.isRunning()) return;
  const s = readSolve(), key = $('#sp').value, rng = num('sr'), fam = $('#sFam').value; let n = Math.round(num('sn'));
  if (s.path === 'custom' && s.pathErr) { B.setStatus(s.pathErr); return; }
  if (!(rng > 1)) { B.setStatus('Range must be greater than 1.'); return; }
  if (!(n >= 3)) n = 3; if (n % 2 === 0) n += 1; n = Math.min(n, 21);
  const partKey = $('#spPart').value; if (key === 'part' && !partKey) { B.setStatus('Pick the part and the dimension to vary.'); return; }
  let md0; try { md0 = Drw.build(D.dr, s.mat, { beam: s.beam }); } catch (e) { B.setStatus(e.message); return; }
  const hasRes = md0.elems.some(e => e.role === 'res') || md0.pmass.some(p => p.role === 'res') || md0.springs.length;
  if ((key === 'resM' || key === 'resK') && !hasRes) { B.setStatus('This drawing has no resonator part, so there is no resonator mass or stiffness to vary. Mark parts as “resonator”.'); return; }
  B.startRun(); const t0 = performance.now(), samples = [], kper = Math.max(5, Math.round(s.kper / 2));
  try {
    const Lat = latFor(md0, s.path), path = Eng.kPath(Lat, kper), fref = s.fmax || (D.run ? D.run.fmax : 0) || 3000;
    for (let j = 0; j < n; j++) {
      const fac = rng ** (2 * j / (n - 1) - 1), sc = { hostE: 1, hostRho: 1, resE: 1, resRho: 1 };
      if (key === 'resM') sc.resRho = fac; if (key === 'resK') sc.resE = fac; if (key === 'hostRho') sc.hostRho = fac; if (key === 'hostE') sc.hostE = fac;
      const dr = (key === 'a' || key === 'part') ? scaledDrawing(D.dr, key, fac, partKey) : D.dr, md = Drw.build(dr, s.mat, { beam: s.beam, scale: sc });
      const L = key === 'a' ? latFor(md, s.path, dr) : Lat, pth = key === 'a' ? Eng.kPath(L, kper) : path, F = [], Z = [];
      for (let i = 0; i < pth.ks.length; i++) { const r = Drw.blochAt(md, pth.ks[i][0], pth.ks[i][1], s.nb, fref * 0.5); F.push(r.f); Z.push(r.z); if (i % 3 === 0) { B.setStatus(`Sweep: sample ${j + 1} of ${n}, k-point ${i + 1} of ${pth.ks.length} (${6 * md.N} DOF)`, (j + i / pth.ks.length) / n); await B.tick(); } }
      const rows = famRows(F, Z, fam); if (!rows) { samples.push({ fac, gaps: [], top: Math.min(...F.map(r => r[r.length - 1])), fr: null }); continue; }
      const top = Math.min(...rows.map(r => r[r.length - 1]));
      samples.push({ fac, gaps: Hom.gapsOf(rows).filter(g => g.lo < top), top, fr: null });
    }
    const nom = samples[(n - 1) / 2], cap = s.fmax > 0 ? s.fmax : Math.min(nom.top * 0.98, D.run && D.runKey === runKey() ? D.run.fmax : Infinity), fminS = Math.min(s.fmin, cap * 0.999);
    samples.forEach(q => { q.gaps = q.gaps.filter(g => g.lo < cap && g.hi > fminS); });
    const lab = key === 'part' ? $('#spPart').selectedOptions[0].textContent : $('#sp').selectedOptions[0].textContent;
    B.S.sweep = { key, label: `${lab} (${fam === 'all' ? 'all waves' : FAMD_LABEL[fam].toLowerCase() + ' waves'})`, samples, P: null, cap, fmin: fminS, drawn: true };
    B.setStatus(`Sweep done in ${((performance.now() - t0) / 1000).toFixed(1)} s (${n} designs, ${path.ks.length} k-points each).`);
    B.drawSweep(); B.renderSens();
  } catch (e) { B.setStatus('Sweep stopped: ' + e.message); if (!/stopped by the user/.test(e.message)) console.error(e); }
  B.endRun();
}

// ---------- edge waves of a drawn cell: strip of N cells, banded FE ----------
async function runEdge() {
  if (B.isRunning()) return;
  if (!D.run || D.runKey !== runKey()) { const ok = await run(); if (!ok) return; }
  const s = readSolve(), R = D.run, md = R.inp.md, Lat = R.inp.Lat;
  const N = Math.max(2, Math.min(16, Math.round(num('eN')))), dir = $('#eDir').value, nk = Math.max(5, Math.min(41, Math.round(num('eK')))), fam = $('#eFam').value;
  const w = Math.max(0.5, Math.min(N / 2, num('eW') || 2)), fmin = R.fmin, fmax = R.fmax, fref = R.fe.fref;
  B.startRun(); const t0 = performance.now();
  try {
    const mdS = Drw.buildStrip(D.dr, s.mat, N, dir, { beam: s.beam }), bS = Eng.recip(mdS.a1, mdS.a2).b1;
    if (6 * mdS.N > 30000) throw new Error(`The strip has ${6 * mdS.N} degrees of freedom; use fewer cells or a larger element size.`);
    const B1 = dir === 'a2' ? Lat.b2 : Lat.b1, B2 = dir === 'a2' ? Lat.b1 : Lat.b2, NU = 8, tol = 0.003 * fmax; let nbB = R.fe.F[0].length, fmaxE = fmax;
    const mus = [], proj = [], modes = []; let nR = 0;
    for (let i = 0; i < nk; i++) {
      const mu = Math.PI * i / (nk - 1); mus.push(mu);
      const per = []; let below = 0; // bulk bands of this family, enough of them to reach the top of the range
      for (let v = 0; v < NU; v++) { const nu = TAU * v / NU, k = [(mu * B1[0] + nu * B2[0]) / TAU, (mu * B1[1] + nu * B2[1]) / TAU];
        let r = Drw.blochAt(md, k[0], k[1], nbB, fref * 0.5), fr = famRows([r.f], [r.z], fam);
        while ((!fr || fr[0][fr[0].length - 1] < fmax) && nbB < Math.min(60, 6 * md.N)) { nbB = Math.min(60, 6 * md.N, Math.ceil(nbB * 1.6)); r = Drw.blochAt(md, k[0], k[1], nbB, fref * 0.5); fr = famRows([r.f], [r.z], fam); }
        if (fr) fmaxE = Math.min(fmaxE, fr[0][fr[0].length - 1]);
        below = Math.max(below, r.f.filter(f => f <= fmax * 1.05).length); per.push(fr); B.setStatus(`Edge analysis: bulk bands at k-point ${i + 1} of ${nk}`, (i + 0.4 * v / NU) / nk); await B.tick(); }
      const nbF = Math.min(...per.map(q => q ? q[0].length : 0)), lo = new Array(nbF).fill(Infinity), hi = new Array(nbF).fill(-Infinity);
      for (const q of per) for (let b = 0; b < nbF; b++) { lo[b] = Math.min(lo[b], q[0][b]); hi[b] = Math.max(hi[b], q[0][b]); }
      proj.push(lo.map((x, b) => [x, hi[b]]));
      const kr = [mu / TAU * bS[0], mu / TAU * bS[1]], nTot = 6 * mdS.N + mdS.springs.length;
      nR = Math.max(nR, Math.min(nTot, Math.ceil((below + 2) * N * 1.15)));
      let r = Drw.stripAt(mdS, kr[0], kr[1], nR, fref, w);
      while (r.f[r.f.length - 1] < fmax && nR < nTot) { nR = Math.min(nTot, Math.ceil(nR * 1.5)); B.setStatus(`Edge analysis: strip k-point ${i + 1} of ${nk}, ${nR} modes`, (i + 0.5) / nk); await B.tick(); r = Drw.stripAt(mdS, kr[0], kr[1], nR, fref, w); }
      const row = [];
      r.f.forEach((f, b) => { if (f > fmax * 1.02 || f < fmin * 0.98) return; if (fam === 'oop' && r.z[b] < 0.5) return; if (fam === 'ip' && r.z[b] >= 0.5) return;
        const inGap = !proj[i].some(([a, c]) => f >= a - tol && f <= c + tol); row.push({ f, bot: r.bot[b], top: r.top[b], L: r.bot[b] + r.top[b], inGap, mu, u: r.vec[b], k: kr }); });
      modes.push(row);
      B.setStatus(`Edge analysis: k-point ${i + 1} of ${nk} done (strip ${r.n} DOF, band ${r.bw})`, (i + 1) / nk); await B.tick();
    }
    const gaps = R.fe[FAMD[fam]] || [], fTop = Math.min(fmax, fmaxE * 0.995); modes.forEach((row, i) => { modes[i] = row.filter(q => q.f <= fTop * 1.02); });
    const uniS = Math.min(1, 2 * w / N);
    B.S.edge = { uni: uniS, Lthr: Math.max(0.6, uniS + 0.25), mus, proj, modes, N, dir, w, fmin, fmax: fTop, mdR: mdS, key: runKey(), family: fam, famLabel: FAMD_LABEL[fam], gaps, drawn: true };
    B.S.emode = null; B.S.esel = null;
    edgeAssessment();
    B.setStatus(`Edge analysis done in ${((performance.now() - t0) / 1000).toFixed(1)} s: strip of ${N} cells, ${6 * mdS.N} DOF.`);
    B.showTab('edge');
  } catch (e) { B.setStatus('Edge analysis stopped: ' + e.message); if (!/stopped by the user/.test(e.message)) console.error(e); }
  B.endRun();
}
function edgeAssessment() {
  const E = B.S.edge, el = $('#assess2'); if (!el) return;
  if (E && E.drawn && E.kind === 'finite') { el.innerHTML = '<h3>Edge waves of the drawn cell</h3>' + finiteSummary().html; return; }
  if (!E || !E.drawn) { if (B.S.zak && B.S.zak.drawn) { el.innerHTML = '<h3>Edge waves of the drawn cell</h3><p>Run the edge analysis to see the edge spectrum.</p>' + zakHtml(); return; } }
  if (!E || !E.drawn) { el.innerHTML = '<h3>Edge waves of the drawn cell</h3><p>Press “Run edge analysis”. The drawn cell is repeated N times across a strip that is periodic along the chosen edge and free on both long sides.</p>'; return; }
  const all = E.modes.flat().filter(q => q.f > 0.01 * E.fmax), edge = all.filter(q => q.inGap && q.L >= (E.Lthr || 0.6)), dn = E.dir === 'a2' ? 'a₂' : 'a₁', fmt = B.fmtF;
  const inG = edge.filter(q => E.gaps.some(g => q.f >= g.lo && q.f <= g.hi));
  const gapTxt = E.gaps.length ? `Bulk gaps for these waves: ${E.gaps.map(g => `${fmt(g.lo)} to ${fmt(g.hi)}`).join(', ')}.` : 'The bulk has no complete gap for these waves in the plotted range.';
  el.innerHTML = `<h3>Edge waves of the drawn cell</h3><p class="verdict">Edge waves (${E.famLabel.toLowerCase()} waves, free edges along ${dn}): ${edge.length ? 'yes' + (inG.length ? ', inside the band gap' : ', outside the band gap') : 'no'}.</p>` +
    (edge.length ? `<p>${edge.length} strip modes lie outside the projected bulk bands with at least ${Math.round(100 * (E.Lthr || 0.6))} % of their energy in the outer ${E.w} cells of the two edges, between ${fmt(Math.min(...edge.map(q => q.f)))} and ${fmt(Math.max(...edge.map(q => q.f)))}${inG.length ? `; ${inG.length} of them inside a complete bulk gap, where the bulk cannot carry them` : ''}. ${gapTxt}</p>` : `<p>No edge-localised strip mode falls outside the bulk bands. ${gapTxt}</p>`) +
    zakHtml() + `<p class="note">The strip is cut along the dashed cell outline exactly as drawn, so which walls end up on the edge depends on how the cell is drawn: move parts, or draw the cell differently, to try another termination. Click a point to see the strip mode.</p>`;
}
function drawEdgeMode(ph) {
  const E = B.S.edge, m = B.S.emode, cv = $('#cvEdgeMode'); if (!E || !E.drawn || !cv) return;
  if (!m) { const { ctx } = B.prep(cv); ctx.fillStyle = css('--ink-2'); ctx.font = '14px system-ui, sans-serif'; ctx.fillText('Click a point on the strip bands.', 16, 30); return; }
  if (!m.fe) m.fe = { kind: 'fe', md: E.mdR, vec: m.u, k: m.k, f: m.f, tiles: [[-1, 0], [0, 0], [1, 0]], caption: `${E.N} cells × 3 periods; grey: undeformed; drag to turn` };
  Viz.modeFE(cv, m.fe, ph == null ? B.S.phase || 0 : ph, 1, D.eview);
}

// ---------- written summary of a drawn-cell run ----------
function famRowsRs(R, fam) { const F = R.fe.F, Z = R.fe.Z, Rs = R.fe.Rs; if (fam === 'gaps') return { F, Rs };
  const want = fam === 'flex', idx = F.map((r, i) => r.map((f, b) => b).filter(b => (Z[i][b] >= 0.5) === want)), n = Math.min(...idx.map(r => r.length));
  return n > 0 ? { F: idx.map((r, i) => r.slice(0, n).map(b => F[i][b])), Rs: idx.map((r, i) => r.slice(0, n).map(b => Rs[i][b])) } : null; }
function drawnAssessment() {
  const el = $('#dAssess'); if (!el) return; const R = D.run;
  if (!R || !R.fe) { el.innerHTML = '<h3>What this geometry gives</h3><p>Run the analysis to get an assessment.</p>'; return; }
  const fmt = B.fmtF, flat = R.fe.flat, cp = R.fe.cplx || {};
  const fams = flat ? [['flex', 'Out-of-plane (flexural)', 'flex'], ['ip', 'In-plane', 'ip'], ['gaps', 'Complete (all waves)', null]] : [['gaps', 'Complete (all waves)', 'all']];
  const one = ([fk, lab, ck]) => { const gaps = R.fe[fk] || [], rr = famRowsRs(R, fk), rows = ck && cp[ck] ? cp[ck] : null;
    if (!gaps.length) return { yes: false, html: `<p><b>${lab} band gap: no.</b> No ${fk === 'gaps' ? 'complete gap' : 'gap for these waves'} between ${fmt(R.fmin)} and ${fmt(R.fmax)}.</p>` };
    const li = gaps.map(g => { const b = g.below - 1, share = rr ? rr.Rs.reduce((a, r) => a + (r[b] || 0), 0) / rr.Rs.length : 0;
      const inG = rows ? rows.filter(r => r.ev && r.ev.length && r.f >= g.lo && r.f <= g.hi) : [], att = inG.length ? Math.max(...inG.map(r => Math.min(...r.ev.map(q => q.att)))) : null;
      const noRes = !R.inp.md.elems.some(e => e.role === 'res') && !R.inp.md.pmass.some(q => q.role === 'res') && !R.inp.md.springs.length;
      const type = noRes ? 'type not known: no part is marked as resonator, so the energy split cannot be measured' : share > 0.5 ? `a local-resonance gap: the band below it carries ${Math.round(100 * share)} % of its kinetic energy in the resonators` : `a Bragg-type gap: the bands that bound it are mainly host waves (resonator share ${Math.round(100 * share)} %)`;
      return `<li><b>${fmt(g.lo)} to ${fmt(g.hi)}</b> (${(200 * (g.hi - g.lo) / (g.hi + g.lo)).toFixed(1)} % wide), ${type}${att != null && isFinite(att) ? `. Peak attenuation ${(8.686 * att).toFixed(1)} dB per cell along a₁` : ''}.</li>`; });
    return { yes: true, html: `<p><b>${lab} band gap: yes.</b> ${gaps.length === 1 ? 'One gap' : gaps.length + ' gaps'} between ${fmt(R.fmin)} and ${fmt(R.fmax)}:</p><ul>${li.join('')}</ul>` }; };
  const parts = fams.map(one), E = B.S.edge, edgeTxt = E && E.drawn && E.key === runKey() ? (E.kind === 'finite' ? finiteSummary() : null) : null;
  let edgeLine = 'not analysed yet (Edge waves tab)';
  if (E && E.drawn && E.key === runKey()) { if (E.kind === 'finite') edgeLine = edgeTxt.short; else { const edge = E.modes.flat().filter(q => q.f > 0.01 * E.fmax && q.inGap && q.L >= (E.Lthr || 0.6)); edgeLine = edge.length ? `yes, ${fmt(Math.min(...edge.map(q => q.f)))} to ${fmt(Math.max(...edge.map(q => q.f)))} (strip)` : 'none found (strip)'; } }
  const yn = v => v ? 'yes' : 'no';
  el.innerHTML = `<h3>What this geometry gives</h3><p class="verdict">${flat ? `Out-of-plane gap: ${yn(parts[0].yes)}. In-plane gap: ${yn(parts[1].yes)}. Complete gap (all waves): ${yn(parts[2].yes)}.` : `Complete gap (all waves): ${yn(parts[0].yes)}.`} Edge waves: ${edgeLine}.</p>${parts.map(p => p.html).join('')}` +
    (edgeTxt ? edgeTxt.html : (B.S.zak && B.S.zak.drawn && B.S.zak.key === runKey() ? zakHtml() : '')) + `<p class="note">${flat ? 'The cell is flat and symmetric about its mid-plane, so out-of-plane and in-plane waves are independent; the complete gap is the safe criterion for vibration of any kind.' : 'The cell is not flat (offset parts or walls), so out-of-plane and in-plane waves couple: only the complete gap is exact.'} The gap type compares the kinetic energy of the band below the gap in the resonators with that in the frame.</p>`;
}

// ---------- finite sample M × N with free edges ----------
async function runFinite(src) {
  if (B.isRunning()) return;
  let dr, gaps, fmaxR, fminR, mat, beam, drawn = !src, key;
  const fam = $('#eFam').value;
  if (drawn) { if (!D.run || D.runKey !== runKey()) { const ok = await run(); if (!ok) return; } const s = readSolve(); dr = D.dr; mat = s.mat; beam = s.beam; gaps = D.run.fe[FAMD[fam]] || []; fmaxR = D.run.fmax; fminR = D.run.fmin; key = runKey(); }
  else ({ dr, gaps, fmax: fmaxR, fmin: fminR, mat, beam, key } = src);
  const M = Math.max(1, Math.min(12, Math.round(num('eM')))), N = Math.max(1, Math.min(12, Math.round(num('eNn')))), nm = Math.max(6, Math.min(300, Math.round(num('eNm')) || 40)), shape = $('#eShape').value, drive = $('#eDrive').value;
  let fc = num('eFc'); if (!(fc > 0)) fc = gaps.length ? (gaps[0].lo + gaps[0].hi) / 2 : (fminR + fmaxR) / 2;
  B.startRun(); const t0 = performance.now();
  try {
    B.setStatus(`Finite sample ${M} × ${N}: building the model…`, 0.05); await B.tick();
    const md = Drw.buildFinite(dr, mat, M, N, { beam, shape });
    B.setStatus(`Finite sample ${M} × ${N} (${6 * md.N} degrees of freedom): factorising and finding ${nm} modes near ${B.fmtF(fc)}… this can take a while`, 0.3); await B.tick();
    const r = Drw.finiteModes(md, fc, nm);
    const uni = r.zoneShare[1] + r.zoneShare[2], mods = r.modes.filter(q => fam === 'all' || (fam === 'oop' ? q.z >= 0.5 : q.z < 0.5)).map(q => {
      const outer = q.edge + q.corner, type = q.corner >= Math.max(0.5, 2.5 * r.zoneShare[2]) ? 'corner' : outer >= Math.max(0.6, uni + 0.25) ? 'edge' : 'bulk';
      return Object.assign(q, { outer, type, inGap: gaps.some(g => q.f >= g.lo && q.f <= g.hi) }); });
    if (!mods.length) throw new Error('No modes of this polarisation near that frequency; try “All waves” or another centre frequency.');
    let drv = null; if (drive !== 'none') { B.setStatus('Finite sample: driven response…', 0.9); await B.tick(); drv = drivenResponse(md, r, drive, fam); }
    B.S.edge = { kind: 'finite', drawn, M, N, md, modes: mods, gaps, fc, fam, famLabel: FAMD_LABEL[fam], uni, zoneShare: r.zoneShare, key, n: r.n, drive: drv, shape: md.fin.shape };
    if (drv) { let bi = 0, bs = -1; drv.ratio.forEach((v, i) => { const sc = v > 0 ? v * (gaps.some(g => drv.freqs[i] >= g.lo && drv.freqs[i] <= g.hi) ? 2 : 1) : -1; if (sc > bs) { bs = sc; bi = i; } }); B.S.dsel = bi; }
    const pick = mods.findIndex(q => q.inGap && q.type !== 'bulk'); B.S.esel = pick >= 0 ? pick : mods.findIndex(q => q.type !== 'bulk'); if (B.S.esel < 0) B.S.esel = 0; B.S.emode = null;
    if (md.fin.shape !== shape) B.setStatus('A staggered rectangle needs a lattice whose a₁, a₂ and a₂ − a₁ have equal length (hexagonal or triangular); a parallelogram was used instead.');
    B.setStatus(`Finite sample done in ${((performance.now() - t0) / 1000).toFixed(1)} s: ${md.fin.shape === 'rect' ? `${N} columns of ${M} cells` : `${M} × ${N} cells`}, ${r.n} degrees of freedom, ${mods.length} modes from ${B.fmtF(mods[0].f)} to ${B.fmtF(mods[mods.length - 1].f)}.`);
    if (drawn) { edgeAssessment(); drawnAssessment(); } B.showTab('edge');
  } catch (e) { B.setStatus('Finite sample stopped: ' + e.message); if (!/stopped by the user/.test(e.message)) console.error(e); }
  B.endRun();
}
function finiteSummary() {
  const E = B.S.edge, fmt = B.fmtF; if (!E || E.kind !== 'finite') return { html: '', short: '', yes: null, where: '' };
  const inGap = E.modes.filter(q => q.inGap), ed = inGap.filter(q => q.type === 'edge'), co = inGap.filter(q => q.type === 'corner'), out = E.modes.filter(q => !q.inGap && q.type !== 'bulk');
  const list = a => a.slice(0, 12).map(q => fmt(q.f)).join(', ') + (a.length > 12 ? ` and ${a.length - 12} more` : '');
  const gapTxt = E.gaps.length ? E.gaps.map(g => `${fmt(g.lo)} to ${fmt(g.hi)}`).join(', ') : 'none in range';
  const yes = ed.length + co.length > 0;
  const html = `<p><b>Finite sample ${E.shape === 'rect' ? `${E.N} columns of ${E.M} cells (staggered rectangle)` : `${E.M} × ${E.N}`}, all edges free (${E.famLabel.toLowerCase()} waves):</b> ${E.modes.length} modes found near ${fmt(E.fc)}, from ${fmt(E.modes[0].f)} to ${fmt(E.modes[E.modes.length - 1].f)}. Bulk band gap${E.gaps.length === 1 ? '' : 's'}: ${gapTxt}.</p><ul>` +
    `<li><b>Edge modes inside the gap:</b> ${ed.length ? `${ed.length}, at ${list(ed)}` : 'none'}.</li>` +
    `<li><b>Corner modes inside the gap:</b> ${co.length ? `${co.length}, at ${list(co)}` : 'none'}.</li>` +
    `<li><b>Bulk-like modes inside the gap:</b> ${inGap.length - ed.length - co.length}${inGap.length - ed.length - co.length ? ' (usually resonators near the free edges that lost their neighbours, or a gap edge that the finite sample shifts)' : ''}.</li>` +
    (out.length ? `<li>Outside the gap, ${out.length} more modes sit mainly on edges or corners (${list(out)}); there they mix with bulk waves.</li>` : '') + `</ul>` +
    `<p class="note">A mode counts as an edge mode when at least ${Math.round(100 * Math.max(0.6, E.uni + 0.25))} % of its kinetic energy lies within half a cell of the sample’s boundary or in the four corner cells (together ${Math.round(100 * E.uni)} % of the sample), and as a corner mode when the corner cells alone carry at least ${Math.round(100 * Math.max(0.5, 2.5 * E.zoneShare[2]))} %. Small samples are mostly boundary; 6 × 6 or larger separates edge and bulk modes more clearly.</p>`;
  return { html: html + drivenSummary(E) + zakHtml(), yes, where: co.length ? ', in the corners' : ed.length ? ', on the edges, inside the gap' : '', short: yes ? `yes, ${ed.length} edge and ${co.length} corner modes in the gap (finite sample)` : 'none in the gap (finite sample)' };
}
let finGeom = null;
function edgeTitles(fin, driven) {
  const set = (id, t) => { const el = $(id); if (el) el.innerHTML = t; };
  set('#edgeT1', driven ? 'Does the energy follow the edge?' : fin ? 'Modes of the finite sample' : 'Strip bands'); set('#edgeT2', driven ? 'Where the energy goes' : fin ? 'Where the mode vibrates' : 'Strip mode shape');
  if (driven) { set('#edgeLede', 'A finite sample with every edge free is driven by a shaker on one cell, like a measurement. For each frequency the plot compares the vibration energy of edge cells with that of interior cells, both counted away from the shaker: above the dashed line the energy follows the edges. Click the curve to see where the energy goes at that frequency.'); set('#edgeLegend', '<span><span class="sw" style="background:var(--edge)"></span>Edge / interior energy ratio</span><span><span class="hatch"></span>Bulk band gap</span>'); return; }
  set('#edgeLede', fin ? 'A finite piece of M × N cells with every edge free is solved for its vibration modes near a chosen frequency (by default the centre of the bulk band gap). Each mode is placed by its frequency and by how much of its energy sits on the edges and corners: modes far to the right inside the hatched gap are edge or corner modes that the infinite lattice cannot carry. Click a mode to see where it vibrates.' : 'A strip N cells wide, repeating along the edge and free on both long sides, is solved for its own bands. Grey areas are the bulk bands projected onto the edge direction. A strip mode outside them that keeps most of its energy in the outer cells is an edge wave; one inside a complete gap is an edge wave the bulk cannot carry.');
  set('#edgeLegend', fin ? '<span><span class="sw" style="background:var(--evan)"></span>Bulk-like mode</span><span><span class="sw" style="background:var(--edge)"></span>Edge mode</span><span><span class="sw" style="background:var(--edge);box-shadow:0 0 0 2px var(--ink);border-radius:50%"></span>Corner mode</span><span><span class="hatch"></span>Bulk band gap</span>' : '<span><span class="sw" style="background:var(--grid);border-radius:2px;border:1px solid var(--rule)"></span>Bulk bands projected on the edge</span><span><span class="sw" style="background:var(--evan)"></span>Bulk-like strip mode</span><span><span class="sw" style="background:var(--edge)"></span>Edge-localised mode</span><span><span class="sw" style="background:var(--edge);box-shadow:0 0 0 2px var(--ink)"></span>Edge wave inside a bulk gap</span><span><span class="hatch"></span>Complete bulk band gap</span>');
}
function drawFinite() {
  const cv = $('#cvEdge'), E = B.S.edge; if (!cv || !E || E.kind !== 'finite') return;
  if ($('#eView').value === 'driven' && E.drive) { edgeTitles(true, true); drawDriven(); return; }
  const { ctx, w, h } = B.prep(cv), fmt = B.fmtF, mods = E.modes, f0 = mods[0].f, f1 = mods[mods.length - 1].f, pad = 0.06 * (f1 - f0 || f0), fmin = f0 - pad, fmax = f1 + pad;
  const m = { l: 58, r: 14, t: 16, b: 40 }, U = B.fUnit(fmax), X = x => m.l + x * (w - m.l - m.r), Y = f => h - m.b - (f - fmin) / (fmax - fmin) * (h - m.t - m.b);
  finGeom = { X, Y };
  ctx.save(); ctx.beginPath(); ctx.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b); ctx.clip();
  for (const g of E.gaps) { const y1 = Y(Math.min(g.hi, fmax)), y0 = Y(Math.max(g.lo, fmin)); if (y0 < m.t || y1 > h - m.b) continue; ctx.fillStyle = css('--gap-tint'); ctx.fillRect(m.l, y1, w - m.l - m.r, y0 - y1); ctx.globalAlpha = 0.45; ctx.fillStyle = B.hatch(ctx, css('--gap')); ctx.fillRect(m.l, y1, w - m.l - m.r, y0 - y1); ctx.globalAlpha = 1; }
  ctx.strokeStyle = css('--ink-2'); ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(X(E.uni), m.t); ctx.lineTo(X(E.uni), h - m.b); ctx.stroke(); ctx.setLineDash([]);
  const evan = css('--evan'), edge = css('--edge'), ink = css('--ink');
  mods.forEach((q, i) => { const x = X(q.outer), y = Y(q.f); ctx.strokeStyle = q.type === 'bulk' ? evan : edge; ctx.lineWidth = q.type === 'bulk' ? 1.5 : 2.5; ctx.beginPath(); ctx.moveTo(x - 9, y); ctx.lineTo(x + 9, y); ctx.stroke();
    if (q.type === 'corner') { ctx.fillStyle = edge; ctx.beginPath(); ctx.arc(x, y, 4.5, 0, TAU); ctx.fill(); ctx.strokeStyle = ink; ctx.lineWidth = 1.2; ctx.stroke(); }
    if (B.S.esel === i) { ctx.strokeStyle = ink; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(x, y, 9, 0, TAU); ctx.stroke(); } });
  ctx.restore();
  ctx.font = '13px system-ui, sans-serif'; ctx.fillStyle = css('--ink-2'); ctx.strokeStyle = css('--rule'); ctx.lineWidth = 1;
  for (const v of B.niceTicks(fmin / U.d, fmax / U.d, 6)) { const y = Y(v * U.d); if (y < m.t || y > h - m.b) continue; ctx.textAlign = 'right'; ctx.fillText(+v.toFixed(3), m.l - 7, y + 4); }
  ctx.textAlign = 'center'; for (const v of [0, 0.25, 0.5, 0.75, 1]) ctx.fillText(Math.round(100 * v) + ' %', X(v), h - m.b + 16);
  ctx.fillText('Share of the mode’s energy on the edges and corners of the sample', (m.l + w - m.r) / 2, h - 6);
  ctx.save(); ctx.translate(15, (m.t + h - m.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(`Frequency (${U.u})`, 0, 0); ctx.restore();
  ctx.textAlign = 'left'; ctx.fillStyle = css('--ink-2'); ctx.font = '12px system-ui, sans-serif'; ctx.fillText('dashed: evenly spread', X(E.uni) + 5, m.t + 12);
  ctx.strokeRect(m.l, m.t, w - m.l - m.r, h - m.t - m.b);
  $('#edgeSub').textContent = `${E.famLabel.toLowerCase()} waves, ${E.M} × ${E.N} cells, all edges free, ${E.n} DOF`; edgeTitles(true);
}
function drawFiniteMode() {
  const cv = $('#cvEdgeMode'), E = B.S.edge; if (!cv || !E || E.kind !== 'finite') return;
  if ($('#eView').value === 'driven' && E.drive) { drawDrivenMap(); return; }
  const { ctx, w, h } = B.prep(cv), q = E.modes[B.S.esel] || E.modes[0], md = E.md, [vr, vi] = q.vec;
  const amp = new Float64Array(md.N); let amax = 0; for (let i = 0; i < md.N; i++) { let a = 0; for (let c = 0; c < 3; c++) a += vr[6 * i + c] ** 2 + vi[6 * i + c] ** 2; amp[i] = a; amax = Math.max(amax, a); }
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity; for (const e of md.elems) for (const P of e.Xd) { x0 = Math.min(x0, P[0]); x1 = Math.max(x1, P[0]); y0 = Math.min(y0, P[1]); y1 = Math.max(y1, P[1]); }
  const sc = Math.min((w - 30) / (x1 - x0), (h - 60) / (y1 - y0)), ox = w / 2 - sc * (x0 + x1) / 2, oy = (h - 30) / 2 + sc * (y0 + y1) / 2, S = P => [ox + sc * P[0], oy - sc * P[1]];
  const hot = css('--edge'), base = css('--rule'), col = a => B.mix(base, hot, Math.min(1, Math.sqrt(a / (amax || 1))));
  const items = md.elems.map(e => ({ e, a: e.nd.reduce((s, [id]) => s + amp[id], 0) / e.nd.length })).sort((a, b) => a.a - b.a);
  for (const { e, a } of items) { const Pp = e.Xd.map(S); if (e.type === 'beam') { ctx.strokeStyle = col(a); ctx.lineWidth = Math.max(1.5, Math.min(8, (e.draw ? e.draw.b : 0.001) * sc)); ctx.lineCap = 'round'; ctx.beginPath(); ctx.moveTo(Pp[0][0], Pp[0][1]); ctx.lineTo(Pp[1][0], Pp[1][1]); ctx.stroke(); }
    else { ctx.fillStyle = col(a); ctx.beginPath(); ctx.moveTo(Pp[0][0], Pp[0][1]); ctx.lineTo(Pp[1][0], Pp[1][1]); ctx.lineTo(Pp[2][0], Pp[2][1]); ctx.closePath(); ctx.fill(); } }
  for (const pm of md.pmass) { const a = amp[pm.node], P = S(md.nodes[pm.node]); ctx.fillStyle = col(a); ctx.fillRect(P[0] - 3, P[1] - 3, 6, 6); }
  ctx.fillStyle = css('--ink-2'); ctx.font = '12.5px system-ui, sans-serif'; ctx.textAlign = 'left'; ctx.fillText('Top view; grey: still, purple: strongest vibration', 10, h - 8);
  $('#edgeModeSub').textContent = `${B.fmtF(q.f)}: ${q.type === 'corner' ? 'corner mode' : q.type === 'edge' ? 'edge mode' : 'bulk-like mode'}${q.inGap ? ' inside the bulk gap' : ''}; ${Math.round(100 * q.edge)} % of its energy on the edges, ${Math.round(100 * q.corner)} % in the corners, ${Math.round(100 * q.bulk)} % inside`;
}
function pickFinite(ev) {
  const E = B.S.edge, g = finGeom; if (!E || E.kind !== 'finite') return;
  if ($('#eView').value === 'driven' && E.drive && B.S.driveGeom) { const r = ev.currentTarget.getBoundingClientRect(), mx = ev.clientX - r.left, Gm = B.S.driveGeom; let bi = 0, bd = Infinity; E.drive.freqs.forEach((f, i) => { const d = Math.abs(Gm.X(f) - mx); if (d < bd) { bd = d; bi = i; } }); B.S.dsel = bi; drawDriven(); drawDrivenMap(); return; }
  if (!g) return;
  const r = ev.currentTarget.getBoundingClientRect(), mx = ev.clientX - r.left, my = ev.clientY - r.top; let best = -1, bd = 16 * 16;
  E.modes.forEach((q, i) => { const d = (g.X(q.outer) - mx) ** 2 + (g.Y(q.f) - my) ** 2; if (d < bd) { bd = d; best = i; } });
  if (best >= 0) { B.S.esel = best; drawFinite(); drawFiniteMode(); }
}

// ---------- driven finite sample: a shaker on one cell, modal frequency response, energy per cell ----------
function sampleFrame(md) { // coordinates with the column direction (a1) pointing up
  const u = md.fin.u, q = [u[1], -u[0]]; return x => [x[0] * q[0] + x[1] * q[1], x[0] * u[0] + x[1] * u[1]]; }
function shakerCell(md, where) {
  const T = sampleFrame(md), P = md.fin.cells.map(c => T(c.c)), xs = P.map(p => p[0]), ys = P.map(p => p[1]);
  const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys), tgt = { bottom: [(x0 + x1) / 2, y0], left: [x0, (y0 + y1) / 2], centre: [(x0 + x1) / 2, (y0 + y1) / 2], corner: [x0, y0] }[where] || [(x0 + x1) / 2, y0];
  let b = 0, d = Infinity; P.forEach((p, k) => { const dd = (p[0] - tgt[0]) ** 2 + (p[1] - tgt[1]) ** 2 * (where === 'bottom' ? 4 : 1); if (dd < d) { d = dd; b = k; } }); return b;
}
function drivenResponse(md, r, where, fam, zeta = 0.02) {
  const cells = md.fin.cells, sh = shakerCell(md, where), c0 = cells[sh].c, n = r.mdiag.length, NN = md.N;
  let fn = 0, dmin = Infinity; md.nodes.forEach((x, i) => { const d = Math.hypot(x[0] * 1000 - c0[0], x[1] * 1000 - c0[1]); if (d < dmin) { dmin = d; fn = i; } });
  const Fv = fam === 'ip' ? [[6 * fn, md.fin.u[0]], [6 * fn + 1, md.fin.u[1]]] : [[6 * fn + 2, 1]];
  const modes = r.modes, f0 = modes[0].f, f1 = modes[modes.length - 1].f, span = f1 - f0, lo = f0 + 0.08 * span, hi = f1 - 0.08 * span, NF = 120, freqs = [], maps = [];
  const pf = modes.map(q => { let a = 0, b = 0; for (const [d, w] of Fv) { a += q.vec[0][d] * w; b -= q.vec[1][d] * w; } return [a, b]; }); // conj(φ)·F
  const lmin = Math.min(...cells.map((c, k) => k === sh ? Infinity : Math.hypot(c.c[0] - c0[0], c.c[1] - c0[1])));
  const far = cells.map(c => Math.hypot(c.c[0] - c0[0], c.c[1] - c0[1]) >= 2.4 * lmin);
  const ratio = [], farShare = [];
  for (let s = 0; s < NF; s++) { const f = lo + (hi - lo) * s / (NF - 1), w = TAU * f, ur = new Float64Array(n), ui = new Float64Array(n);
    modes.forEach((q, m) => { const wn = TAU * q.f, dr = wn * wn - w * w, di = 2 * zeta * wn * w, den = dr * dr + di * di, [a, b] = pf[m], cr = (a * dr + b * di) / den, ci = (b * dr - a * di) / den, [vr, vi] = q.vec;
      for (let i = 0; i < 6 * NN; i++) { if (i % 6 > 2) continue; ur[i] += cr * vr[i] - ci * vi[i]; ui[i] += cr * vi[i] + ci * vr[i]; } });
    const E = new Float64Array(cells.length); for (let i = 0; i < 6 * NN; i++) { if (i % 6 > 2) continue; E[md.fin.cellOf[Math.floor(i / 6)]] += r.mdiag[i] * (ur[i] * ur[i] + ui[i] * ui[i]); }
    const tot = E.reduce((a, b) => a + b, 0) || 1; for (let k = 0; k < E.length; k++) E[k] /= tot;
    let se = 0, ne = 0, sb = 0, nbk = 0, sf = 0; cells.forEach((c, k) => { if (!far[k]) return; sf += E[k]; if (c.edge) { se += E[k]; ne++; } else { sb += E[k]; nbk++; } });
    freqs.push(f); maps.push(Float32Array.from(E)); ratio.push(ne && nbk && sb > 0 ? (se / ne) / (sb / nbk) : null); farShare.push(sf); }
  return { shaker: sh, freqs, maps, ratio, farShare, zeta, where, valid: [lo, hi], enoughInterior: far.some((v, k) => v && !cells[k].edge) };
}
function drawDriven() {
  const cv = $('#cvEdge'), E = B.S.edge, Dv = E && E.drive; if (!cv || !Dv) return;
  const { ctx, w, h } = B.prep(cv), m = { l: 64, r: 14, t: 16, b: 40 }, f0 = Dv.freqs[0], f1 = Dv.freqs[Dv.freqs.length - 1], U = B.fUnit(f1);
  const rs = Dv.ratio.filter(v => v > 0), rmax = Math.max(30, ...rs), rmin = Math.min(0.3, ...rs), X = f => m.l + (f - f0) / (f1 - f0) * (w - m.l - m.r), Y = v => h - m.b - (Math.log10(v) - Math.log10(rmin)) / (Math.log10(rmax) - Math.log10(rmin)) * (h - m.t - m.b);
  B.S.driveGeom = { X, f0, f1 };
  ctx.save(); ctx.beginPath(); ctx.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b); ctx.clip();
  for (const g of E.gaps) { const x0 = X(Math.max(g.lo, f0)), x1 = X(Math.min(g.hi, f1)); if (x1 <= x0) continue; ctx.fillStyle = css('--gap-tint'); ctx.fillRect(x0, m.t, x1 - x0, h - m.t - m.b); ctx.globalAlpha = 0.45; ctx.fillStyle = B.hatch(ctx, css('--gap')); ctx.fillRect(x0, m.t, x1 - x0, h - m.t - m.b); ctx.globalAlpha = 1; }
  ctx.strokeStyle = css('--ink-2'); ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(m.l, Y(1)); ctx.lineTo(w - m.r, Y(1)); ctx.stroke(); ctx.setLineDash([]);
  ctx.strokeStyle = css('--edge'); ctx.lineWidth = 2.2; ctx.beginPath(); let on = false; Dv.ratio.forEach((v, i) => { if (!(v > 0)) { on = false; return; } const x = X(Dv.freqs[i]), y = Y(v); on ? ctx.lineTo(x, y) : ctx.moveTo(x, y); on = true; }); ctx.stroke();
  const si = B.S.dsel ?? 0; ctx.strokeStyle = css('--ink'); ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(X(Dv.freqs[si]), m.t); ctx.lineTo(X(Dv.freqs[si]), h - m.b); ctx.stroke();
  ctx.restore();
  ctx.font = '13px system-ui, sans-serif'; ctx.fillStyle = css('--ink-2'); ctx.strokeStyle = css('--rule'); ctx.lineWidth = 1; ctx.textAlign = 'right';
  for (const v of [0.1, 0.3, 1, 3, 10, 30, 100, 300]) { if (v < rmin || v > rmax) continue; ctx.fillText(String(v), m.l - 7, Y(v) + 4); }
  ctx.textAlign = 'center'; for (const v of B.niceTicks(f0 / U.d, f1 / U.d, 6)) { const x = X(v * U.d); if (x < m.l || x > w - m.r) continue; ctx.fillText(+v.toFixed(3), x, h - m.b + 16); }
  ctx.fillText(`Frequency (${U.u})`, (m.l + w - m.r) / 2, h - 6);
  ctx.save(); ctx.translate(14, (m.t + h - m.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText('Energy per cell: edge / interior', 0, 0); ctx.restore();
  ctx.textAlign = 'left'; ctx.font = '12px system-ui, sans-serif'; ctx.fillText('dashed: as much on the edge as inside', m.l + 6, Y(1) - 5);
  ctx.strokeRect(m.l, m.t, w - m.l - m.r, h - m.t - m.b);
  $('#edgeSub').textContent = `shaker at the ${({ bottom: 'middle of the bottom edge', left: 'middle of the left edge', centre: 'centre', corner: 'bottom-left corner' })[Dv.where]}, ${Math.round(100 * Dv.zeta)} % damping; cells more than 2.4 cells from the shaker`;
}
function drawDrivenMap() {
  const cv = $('#cvEdgeMode'), E = B.S.edge, Dv = E && E.drive; if (!cv || !Dv) return;
  const { ctx, w, h } = B.prep(cv), md = E.md, T = sampleFrame(md), si = B.S.dsel ?? 0, Em = Dv.maps[si], emax = Math.max(...Em);
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity; for (const e of md.elems) for (const P of e.Xd) { const q = T(P); x0 = Math.min(x0, q[0]); x1 = Math.max(x1, q[0]); y0 = Math.min(y0, q[1]); y1 = Math.max(y1, q[1]); }
  const sc = Math.min((w - 30) / (x1 - x0), (h - 60) / (y1 - y0)), ox = w / 2 - sc * (x0 + x1) / 2, oy = (h - 30) / 2 + sc * (y0 + y1) / 2, S = P => { const q = T(P); return [ox + sc * q[0], oy - sc * q[1]]; };
  const lo = css('--rule'), hot = css('--edge'), col = e => { const t = Math.max(0, Math.min(1, (10 * Math.log10(e / emax) + 30) / 30)); return B.mix(lo, hot, t); };
  const items = md.elems.map(e => ({ e, v: Em[md.fin.cellOf[e.nd[0][0]]] })).sort((a, b) => a.v - b.v);
  for (const { e, v } of items) { const P = e.Xd.map(S), c = col(v); if (e.type === 'beam') { ctx.strokeStyle = c; ctx.lineWidth = Math.max(1.5, Math.min(8, (e.draw ? e.draw.b : 0.001) * sc)); ctx.lineCap = 'round'; ctx.beginPath(); ctx.moveTo(P[0][0], P[0][1]); ctx.lineTo(P[1][0], P[1][1]); ctx.stroke(); }
    else { ctx.fillStyle = c; ctx.beginPath(); ctx.moveTo(P[0][0], P[0][1]); ctx.lineTo(P[1][0], P[1][1]); ctx.lineTo(P[2][0], P[2][1]); ctx.closePath(); ctx.fill(); } }
  const sp = S(md.fin.cells[Dv.shaker].c.map(v => v / 1000)); ctx.fillStyle = '#fff'; ctx.strokeStyle = css('--ink'); ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(sp[0], sp[1], 7, 0, TAU); ctx.fill(); ctx.stroke();
  ctx.fillStyle = css('--ink-2'); ctx.font = '12.5px system-ui, sans-serif'; ctx.textAlign = 'left'; ctx.fillText('Purple: strongest; grey: 30 dB weaker; ○ shaker', 10, h - 8);
  const r = Dv.ratio[si]; $('#edgeModeSub').textContent = `${B.fmtF(Dv.freqs[si])}: edge cells carry ${r > 0 ? r.toFixed(1) + ' ×' : '–'} the energy of interior cells; ${Math.round(100 * Dv.farShare[si])} % of the energy is more than 2.4 cells from the shaker`;
}
function drivenSummary(E) {
  const Dv = E.drive; if (!Dv) return '';
  const fmt = B.fmtF, bands = []; let st = null;
  Dv.ratio.forEach((v, i) => { const on = v >= 3; if (on && st == null) st = i; if ((!on || i === Dv.ratio.length - 1) && st != null) { const e = on ? i : i - 1; bands.push([Dv.freqs[st], Dv.freqs[e], Math.max(...Dv.ratio.slice(st, e + 1))]); st = null; } });
  if (!Dv.enoughInterior) return '<p><b>Driven response:</b> the sample is too small to have interior cells far from the shaker; use at least 6 × 6 cells.</p>';
  return `<p><b>Driven response (shaker at the ${({ bottom: 'middle of the bottom edge', left: 'middle of the left edge', centre: 'centre', corner: 'corner' })[Dv.where]}):</b> ${bands.length ? `the energy follows the edges (edge cells at least 3 × the interior, away from the shaker) at ${bands.map(b => `${fmt(b[0])}–${fmt(b[1])} (up to ${b[2].toFixed(0)} ×)`).join(', ')}` : 'the energy does not stay on the edges anywhere in this range'}. The response is computed from the ${E.modes.length ? 'modes found' : 'modes'} and is reliable between ${fmt(Dv.valid[0])} and ${fmt(Dv.valid[1])}; find more modes to widen it.</p>`;
}

// ---------- topology check: Zak phase (Wilson loop) of the bands below each gap ----------
async function zakCheck() {
  if (B.isRunning()) return;
  const fam = $('#eFam').value, drawn = B.isDrawn ? B.isDrawn() : $('#lat').value === 'drawn';
  let md, dr, gaps, nbRun, key;
  if (drawn) { if (!D.run || D.runKey !== runKey()) { const ok = await run(); if (!ok) return; } md = D.run.inp.md; dr = D.dr; gaps = D.run.fe[FAMD[fam]] || []; nbRun = D.run.fe.F[0].length; key = runKey(); }
  else { const src = B.paramSource ? await B.paramSource(fam) : null; if (!src) return; dr = src.dr; gaps = src.gaps; key = src.key; nbRun = 16; md = Drw.build(dr, src.mat, { beam: src.beam }); }
  if (!gaps.length) { B.setStatus('No band gap of this wave type in the range: there is nothing to check. Run the analysis first, or choose another wave polarisation.'); return; }
  B.startRun(); const t0 = performance.now();
  try {
    const sym = Drw.symmetry(dr), n = 6 * md.N + md.springs.length, { Mr } = Drw.blochMats(md, 0, 0), mdiag = Float64Array.from({ length: n }, (_, i) => Mr[i][i]);
    const pos = Array.from({ length: n }, (_, i) => md.nodes[i < 6 * md.N ? Math.floor(i / 6) : md.springs[i - 6 * md.N].node]);
    const { b1, b2 } = Eng.recip(md.a1, md.a2), a1 = dr.a1, a2 = dr.a2;
    // reference point: the inversion centre nearest the middle of the drawn parts (cells are cut along the drawn outline)
    let ref = [0, 0]; if (sym.centre) { const v = []; for (const p of dr.parts) try { for (const q of vertsOf(p)) v.push(q); } catch (e) {}
      const cg = v.length ? [v.reduce((s, q) => s + q[0], 0) / v.length, v.reduce((s, q) => s + q[1], 0) / v.length] : [0, 0]; let best = Infinity;
      for (let i = -3; i <= 3; i++) for (let j = -3; j <= 3; j++) { const c = [sym.centre[0] + (i * a1[0] + j * a2[0]) / 2, sym.centre[1] + (i * a1[1] + j * a2[1]) / 2], d = Math.hypot(c[0] - cg[0], c[1] - cg[1]); if (d < best) { best = d; ref = c; } } }
    const refm = [ref[0] / 1000, ref[1] / 1000], fref = gaps[0].lo, Ns = 24;
    const dirs = [['edges parallel to a₂', b1], ['edges parallel to a₁', b2], ['edges parallel to a₂ − a₁', [b1[0] + b2[0], b1[1] + b2[1]]]];
    const want = q => fam === 'all' || (fam === 'oop' ? q >= 0.5 : q < 0.5), rows = [];
    let step = 0; const total = gaps.length * dirs.length * Ns;
    for (const g of gaps) { const fcut = (g.lo + g.hi) / 2, row = { g, dirs: [] };
      for (const [lab, bc] of dirs) {
        const V = [];
        for (let s = 0; s < Ns; s++) { const k = [bc[0] * s / Ns, bc[1] * s / Ns]; let nb = Math.max(nbRun, 8), r;
          for (;;) { r = Drw.blochAt(md, k[0], k[1], nb, fref, true); if (r.f[r.f.length - 1] > fcut || nb >= 6 * md.N) break; nb = Math.min(6 * md.N, nb + 8); }
          V.push(r.f.map((f, b) => [f, b]).filter(([f, b]) => f < fcut && want(r.z[b])).map(([, b]) => r.vec[b]));
          B.setStatus(`Topology check: gap ${B.fmtF(g.lo)}–${B.fmtF(g.hi)}, ${lab}`, ++step / total); if (s % 3 === 0) await B.tick(); }
        const nocc = V[0].length; if (!nocc || V.some(v => v.length !== nocc)) { row.dirs.push({ lab, ok: false }); continue; }
        V.push(V[0]); const dk = [bc[0] / Ns, bc[1] / Ns], wr = pos.map(x => Math.cos(-(dk[0] * x[0] + dk[1] * x[1]))), wi = pos.map(x => Math.sin(-(dk[0] * x[0] + dk[1] * x[1])));
        let ph = 0;
        for (let s = 0; s < Ns; s++) { const U = [];
          for (let i = 0; i < nocc; i++) { U.push([]); const [ar, ai] = V[s][i];
            for (let j = 0; j < nocc; j++) { const [br, bi] = V[s + 1][j]; let sr = 0, si = 0;
              for (let q = 0; q < n; q++) { const m = mdiag[q]; if (!m) continue; const cr = br[q] * wr[q] - bi[q] * wi[q], ci = br[q] * wi[q] + bi[q] * wr[q]; sr += m * (ar[q] * cr + ai[q] * ci); si += m * (ar[q] * ci - ai[q] * cr); }
              U[i].push([sr, si]); } }
          const [dr0, di0] = cdet(U); ph += Math.atan2(di0, dr0); }
        let z = (-ph - nocc * (bc[0] * refm[0] + bc[1] * refm[1])) % TAU; if (z < 0) z += TAU; if (z > TAU - 1e-9) z = 0;
        const q = sym.c2 ? (Math.abs(z - Math.PI) < Math.PI / 2 ? 1 : 0) : null, dev = sym.c2 ? Math.min(Math.abs(z - q * Math.PI), Math.abs(z - q * Math.PI - TAU)) : 0;
        row.dirs.push({ lab, ok: true, z, q, dev, nocc }); }
      rows.push(row); }
    B.S.zak = { key, drawn, fam, rows, c2: sym.c2, ref };
    B.setStatus(`Topology check done in ${((performance.now() - t0) / 1000).toFixed(0)} s.`);
    if (drawn) { edgeAssessment(); drawnAssessment(); } else if (B.renderAssessment) B.renderAssessment();
    B.showTab('edge');
  } catch (e) { B.setStatus('Topology check stopped: ' + e.message); if (!/stopped by the user/.test(e.message)) console.error(e); }
  B.endRun();
}
function cdet(A) { const N = A.length, a = A.map(r => r.map(z => z.slice())); let dr = 1, di = 0;
  for (let c = 0; c < N; c++) { let p = c; for (let r = c + 1; r < N; r++) if (Math.hypot(...a[r][c]) > Math.hypot(...a[p][c])) p = r; if (p !== c) { [a[p], a[c]] = [a[c], a[p]]; dr = -dr; di = -di; }
    const [pr, pi] = a[c][c], nr = dr * pr - di * pi; di = dr * pi + di * pr; dr = nr; const den = pr * pr + pi * pi; if (!den) return [0, 0];
    for (let r = c + 1; r < N; r++) { const [xr, xi] = a[r][c], fr = (xr * pr + xi * pi) / den, fi = (xi * pr - xr * pi) / den; for (let q = c; q < N; q++) { const [yr, yi] = a[c][q]; a[r][q][0] -= fr * yr - fi * yi; a[r][q][1] -= fr * yi + fi * yr; } } }
  return [dr, di]; }
function zakHtml() {
  const Z = B.S.zak; if (!Z) return '';
  const fmt = B.fmtF, any = Z.rows.some(r => r.dirs.some(d => d.ok && d.q === 1));
  const rows = Z.rows.map(r => `<li><b>Gap ${fmt(r.g.lo)}–${fmt(r.g.hi)}</b> (${r.dirs.find(d => d.ok) ? r.dirs.find(d => d.ok).nocc : '?'} band${(r.dirs.find(d => d.ok) || {}).nocc === 1 ? '' : 's'} below): ` +
    r.dirs.map(d => !d.ok ? `${d.lab}: not defined (the bands below the gap change along the loop)` : Z.c2 ? `${d.lab}: <b>${d.q ? 'π' : '0'}</b>${d.dev > 0.15 ? ` (off by ${(d.dev / Math.PI).toFixed(2)}π: refine the mesh)` : ''}` : `${d.lab}: ${(d.z / Math.PI).toFixed(2)}π`).join('; ') + '.</li>').join('');
  return `<h3 style="margin-top:14px">Topology check</h3><p><b>Chern number: 0</b> for every band, as in any passive structure (time reversal makes the Berry curvature odd in k). <b>Zak phase</b> of the bands below each gap, measured from the cell centre (${Z.c2 ? 'the inversion centre nearest the middle of the drawn cell' : 'the cell has no 180° symmetry, so the phase is not quantised'}):</p><ul>${rows}</ul>` +
    (Z.c2 ? `<p>${any ? '<b>π</b> means the vibration of those bands is centred on the cell boundary: a sample cut along the drawn cell outline must carry edge states in that gap on edges parallel to that direction (a polarisation edge state, protected by the 180° symmetry).' : 'All phases are <b>0</b>: the bands below the gaps vibrate around the cell centre, so no edge state is required by topology. Edge waves that the strip or the finite sample show are then ordinary (trivial) edge states of the cut: real, but movable by changing the edge.'}</p>` : '');
}

// ---------- init ----------
function init() {
  try { const s = localStorage.getItem(ST_KEY); if (s) D.dr = JSON.parse(s); } catch (e) {}
  if (!D.dr || !Array.isArray(D.dr.parts)) { D.dr = { a1: [14, 0], a2: [7, 12.124356], mesh: 2, grid: 0.5, parts: [] }; }
  $('#dTools').innerHTML = TOOLS.map(([t, lab]) => `<button type="button" role="radio" data-dtool="${t}" aria-checked="${t === 'select'}">${lab}</button>`).join('');
  document.querySelectorAll('[data-dtool]').forEach(b => b.addEventListener('click', () => setTool(b.dataset.dtool)));
  $('#dAutoRole').addEventListener('click', autoRoles);
  document.querySelectorAll('[data-drole]').forEach(b => b.addEventListener('click', () => { D.roleAuto = false; D.role = b.dataset.drole; document.querySelectorAll('[data-drole]').forEach(x => x.setAttribute('aria-checked', x === b)); }));
  document.querySelectorAll('[data-dshow]').forEach(b => b.addEventListener('click', () => { D.show = b.dataset.dshow; document.querySelectorAll('[data-dshow]').forEach(x => x.setAttribute('aria-checked', x === b)); drawBands(); }));
  ['dA1x', 'dA1y', 'dA2x', 'dA2y'].forEach(id => $('#' + id).addEventListener('change', () => { commit(); D.dr.latFollow = false; $('#dLatFollow').checked = false; D.dr.a1 = [num('dA1x'), num('dA1y')]; D.dr.a2 = [num('dA2x'), num('dA2y')]; D.view = null; changed(); }));
  $('#dMesh').addEventListener('change', () => { const v = num('dMesh'); if (v > 0) { D.dr.mesh = v; changed(); } });
  $('#dGrid').addEventListener('change', () => { D.dr.grid = Math.max(0, num('dGrid') || 0); save(); drawEditor(); });
  $('#dLatFit').addEventListener('click', () => { fitLattice(false); changed(); });
  $('#dLatWarnFit').addEventListener('click', () => { if (D.latPart) { fitLattice(false, D.latPart); changed(); } });
  $('#dLatFollow').addEventListener('change', e => { D.dr.latFollow = e.target.checked; if (e.target.checked) fitLattice(false); changed(); });
  $('#dLatSet').addEventListener('click', () => { const t = $('#dLatP').value, a = num('dLatA'), b = num('dLatB') || a; if (!(a > 0)) return; commit(); D.dr.latFollow = false; $('#dLatFollow').checked = false;
    if (t === 'square') { D.dr.a1 = [a, 0]; D.dr.a2 = [0, a]; } else if (t === 'hex') { D.dr.a1 = [a, 0]; D.dr.a2 = [a / 2, +(a * Math.sqrt(3) / 2).toFixed(6)]; } else { D.dr.a1 = [a, 0]; D.dr.a2 = [0, b]; }
    D.view = null; syncLat(); changed(); });
  $('#dClear').addEventListener('click', () => { commit(); D.dr.parts = []; D.sel = null; changed(); });
  $('#dUndo').addEventListener('click', undo); $('#dRedo').addEventListener('click', redo); $('#dFit').addEventListener('click', () => { D.view = null; drawEditor(); });
  $('#dMeshShow').addEventListener('change', () => { D.mesh = $('#dMeshShow').checked; drawEditor(); });
  $('#dSave').addEventListener('click', download); $('#dLoad').addEventListener('change', e => { if (e.target.files[0]) load(e.target.files[0]); e.target.value = ''; });
  $('#gRebuild').addEventListener('click', () => regenerate(true));
  Aids.zoomable($('#cvDrawBand'), 'drawn', drawBands); Aids.zoomable($('#cvDrawCplx'), 'drawn', drawBands); $('#cvDrawBand').addEventListener('click', pickMode); Viz.register('cvDrawMode', drawMode);
  { const cv = $('#cvDrawMode'); let dg = null; cv.addEventListener('pointerdown', e => { dg = { x: e.clientX, y: e.clientY, yaw: D.mview.yaw, elev: D.mview.elev }; cv.setPointerCapture(e.pointerId); }); cv.addEventListener('pointermove', e => { if (!dg) return; D.mview.yaw = dg.yaw - (e.clientX - dg.x) * 0.01; D.mview.elev = Math.max(0, Math.min(Math.PI / 2, dg.elev + (e.clientY - dg.y) * 0.01)); }); cv.addEventListener('pointerup', () => { dg = null; }); }
  $('#dTt').addEventListener('change', () => { const t = $('#dTt').value; document.querySelectorAll('[data-dtt]').forEach(el => el.classList.toggle('hidden', !el.dataset.dtt.split(' ').includes(t))); });
  $('#dTt').dispatchEvent(new Event('change'));
  const tgl = (sel, attr) => { const v = $(sel).value; document.querySelectorAll(`[${attr}]`).forEach(el => el.classList.toggle('hidden', !el.getAttribute(attr).split(' ').includes(v))); };
  $('#dSty').addEventListener('change', () => { tgl('#dSty', 'data-dst'); const m = $('#dSm'); [...m.options].forEach(o => { o.hidden = !o.dataset.for.split(' ').includes($('#dSty').value); }); const vis = [...m.options].find(o => !o.hidden); if (vis && m.selectedOptions[0].hidden) m.value = vis.value; });
  $('#dSu').addEventListener('change', () => tgl('#dSu', 'data-dsu'));
  $('#dSty').dispatchEvent(new Event('change')); $('#dSu').dispatchEvent(new Event('change'));
  $('#dSadd').addEventListener('click', () => { const x = num('dScx'), y = num('dScy'); if (!isFinite(x) || !isFinite(y)) { B.setStatus('Type the centre x and y of the shape.'); return; } addShape([x, y]); });
  { const cv = $('#cvEdgeMode'); let dg = null; cv.addEventListener('pointerdown', e => { if (!B.S.edge || !B.S.edge.drawn) return; dg = { x: e.clientX, y: e.clientY, yaw: D.eview.yaw, elev: D.eview.elev }; cv.setPointerCapture(e.pointerId); });
    cv.addEventListener('pointermove', e => { if (!dg) return; D.eview.yaw = dg.yaw - (e.clientX - dg.x) * 0.01; D.eview.elev = Math.max(0, Math.min(Math.PI / 2, dg.elev + (e.clientY - dg.y) * 0.01)); if (B.S.emode) B.S.emode.fe = B.S.emode.fe && Object.assign(B.S.emode.fe, { box: null }); drawEdgeMode(); }); cv.addEventListener('pointerup', () => { dg = null; }); }
  sweepParts();
  bindEditor(); bind3D(); syncLat(); setTool('select'); renderProps();
  setTimeout(redraw, 0);
}
function redraw() { drawEditor(); scheduleModel(); draw3D(); drawBands(); }
function onRail() { D.md = null; if ($('#lat').value === 'drawn') { regenerate(false); scheduleModel(); } }
init();
return { zakCheck, zakHtml, drivenSummary, select: i => { D.sel = i; renderProps(); drawEditor(); }, edgeTitles, runFinite, finiteSummary, drawFinite, drawFiniteMode, pickFinite, drawnAssessment, regenerate, generate, readGen, redraw, onRail, run, runSweep, runEdge, drawEdgeMode, edgeAssessment, lattice, importLeft, sweepParts, D };
})();
