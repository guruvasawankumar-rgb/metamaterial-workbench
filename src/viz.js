// ================= Shared plots for the homogenised and drawn-cell tabs: complex band, mode shapes =================
const Viz = (() => {
const B = window.MMW, TAU = 2 * Math.PI, css = B.css;
const ATT_MAX = 2.5, PROP = 0.015;

// complex band: left half attenuation Im(μ) (Np/cell), right half phase Re(μ); sets = [{rows:[{f, ev:[{att, ph}]}], color, shape:'dot'|'square', label}]
function complexPanel(cv, o) {
  const { ctx, w, h } = B.prep(cv), ink2 = css('--ink-2'), M = { t: 16, b: 36 }, fmin = o.fmin, fmax = o.fmax;
  if (!o.sets || !o.sets.length) { ctx.fillStyle = ink2; ctx.font = '15px system-ui, sans-serif'; ctx.fillText(o.empty || 'Runs together with the bands.', 20, 40); return null; }
  const Y = f => h - M.b - (f - fmin) / (fmax - fmin) * (h - M.t - M.b), l = 16, r = 12, mid = l + (w - l - r) * 0.44;
  const XA = a => mid - (a / ATT_MAX) * (mid - l), XP = p => mid + 6 + (p / Math.PI) * (w - r - mid - 6);
  for (const g of o.gaps || []) { const y1 = Y(Math.min(g.hi, fmax)), y0 = Y(Math.max(g.lo, fmin)); if (g.hi < fmin || g.lo > fmax || y0 < M.t) continue; ctx.fillStyle = css('--gap-tint'); ctx.fillRect(l, y1, w - l - r, y0 - y1); ctx.fillStyle = B.hatch(ctx, css('--gap')); ctx.globalAlpha = 0.4; ctx.fillRect(l, y1, w - l - r, y0 - y1); ctx.globalAlpha = 1; }
  ctx.strokeStyle = css('--grid'); ctx.lineWidth = 1; ctx.font = '13px system-ui, sans-serif'; ctx.fillStyle = ink2; ctx.textAlign = 'center';
  const U = B.fUnit(fmax);
  for (const v of B.niceTicks(fmin / U.d, fmax / U.d, 6)) { const y = Y(v * U.d); ctx.beginPath(); ctx.moveTo(l, y); ctx.lineTo(w - r, y); ctx.stroke(); }
  for (const a of [0.5, 1, 1.5, 2, 2.5]) { const x = XA(a); ctx.beginPath(); ctx.moveTo(x, M.t); ctx.lineTo(x, h - M.b); ctx.stroke(); ctx.fillText(a, x, h - M.b + 16); }
  [[0, '0'], [Math.PI / 2, 'π/2'], [Math.PI, 'π']].forEach(([p, s]) => { const x = XP(p); ctx.beginPath(); ctx.moveTo(x, M.t); ctx.lineTo(x, h - M.b); ctx.stroke(); ctx.fillText(s, x, h - M.b + 16); });
  ctx.strokeStyle = css('--ink'); ctx.lineWidth = 1.2; ctx.beginPath(); ctx.moveTo(mid + 3, M.t); ctx.lineTo(mid + 3, h - M.b); ctx.stroke();
  ctx.fillStyle = ink2; ctx.fillText('Attenuation Im(μ), Np/cell', (l + mid) / 2, h - 4); ctx.fillText('Phase Re(μ), rad/cell', (mid + w - r) / 2, h - 4);
  ctx.save(); ctx.beginPath(); ctx.rect(l, M.t - 4, w - l - r, h - M.t - M.b + 4); ctx.clip();
  for (const set of o.sets) for (const row of set.rows) { if (!row.ev) continue; const y = Y(row.f);
    for (const e of row.ev) { if (e.att > ATT_MAX * 1.05) continue; const prop = e.att < PROP; ctx.fillStyle = set.color; ctx.globalAlpha = prop ? 1 : 0.75;
      const dot = (x) => { if (set.shape === 'square') ctx.fillRect(x - 1.7, y - 1.7, 3.4, 3.4); else { ctx.beginPath(); ctx.arc(x, y, 1.7, 0, TAU); ctx.fill(); } };
      dot(XP(e.ph)); if (!prop) dot(XA(e.att)); } }
  ctx.globalAlpha = 1; ctx.restore();
  return { Y };
}
// least attenuation per frequency, for the assessment text
function minAtt(rows) { return rows.map(r => ({ f: r.f, a: r.ev && r.ev.length ? Math.min(...r.ev.map(e => e.att)) : NaN })); }

// ---------- mode shapes ----------
// diverging colour: -1 host colour, 0 panel, +1 resonator colour
function divColor(v) { const host = css('--host'), res = css('--res'), panel = css('--panel'); return v >= 0 ? B.mix(panel, res, Math.min(1, v)) : B.mix(panel, host, Math.min(1, -v)); }
function cellView(w, h, a1, a2, tiles, pad = 26) {
  const pts = []; for (const [i, j] of tiles) for (const [u, v] of [[0, 0], [1, 0], [1, 1], [0, 1]]) pts.push([(i + u) * a1[0] + (j + v) * a2[0], (i + u) * a1[1] + (j + v) * a2[1]]);
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity; for (const p of pts) { x0 = Math.min(x0, p[0]); x1 = Math.max(x1, p[0]); y0 = Math.min(y0, p[1]); y1 = Math.max(y1, p[1]); }
  const sc = Math.min((w - 2 * pad) / (x1 - x0), (h - 2 * pad - 16) / (y1 - y0)), ox = w / 2 - sc * (x0 + x1) / 2, oy = (h - 16) / 2 + sc * (y0 + y1) / 2;
  return { S: p => [ox + sc * p[0], oy - sc * p[1]], sc };
}
// mode of the homogenised plate. m: {fam, k:[kx,ky], G, W:[re,im] of length nw, a1, a2, lines:[[p0,p1]], res:[{pos, abs:[re,im] (z or [ux,uy]) }]}
function modePWE(cv, m, phase, amp) {
  const { ctx, w, h } = B.prep(cv), tiles = Aids.tiles();
  const V = cellView(w, h, m.a1, m.a2, tiles), Ng = m.G.length, cph = Math.cos(phase), sph = Math.sin(phase);
  const field = (x, y, off) => { let re = 0, im = 0; for (let g = 0; g < Ng; g++) { const t = -((m.k[0] + m.G[g][0]) * x + (m.k[1] + m.G[g][1]) * y), c = Math.cos(t), s = Math.sin(t), wr = m.W[0][off + g], wi = m.W[1][off + g]; re += wr * c - wi * s; im += wr * s + wi * c; } return [re, im]; };
  if (!m.grid || m.grid.key !== Aids.key()) { // sample once per tiling
    const N = Math.max(12, Math.round(64 / Math.sqrt(tiles.length / 9))), cells = []; let mx = 0;
    for (const [ti, tj] of tiles) for (let a = 0; a < N; a++) for (let b = 0; b < N; b++) { const u = ti + (a + 0.5) / N, v = tj + (b + 0.5) / N, x = u * m.a1[0] + v * m.a2[0], y = u * m.a1[1] + v * m.a2[1];
      if (m.fam === 'flex') { const f = field(x, y, 0); mx = Math.max(mx, Math.hypot(...f)); cells.push({ u, v, f }); }
      else if (a % 6 === 3 && b % 6 === 3) { const fu = field(x, y, 0), fv = field(x, y, Ng); mx = Math.max(mx, Math.hypot(...fu, ...fv)); cells.push({ u, v, fu, fv, x, y }); } }
    m.grid = { N, cells, mx: mx || 1, key: Aids.key() };
  }
  const G = m.grid, du = 1 / G.N;
  if (m.fam === 'flex') {
    for (const c of G.cells) { const val = (c.f[0] * cph - c.f[1] * sph) / G.mx; ctx.fillStyle = divColor(val * amp);
      const p = [[c.u - du / 2, c.v - du / 2], [c.u + du / 2, c.v - du / 2], [c.u + du / 2, c.v + du / 2], [c.u - du / 2, c.v + du / 2]].map(([u, v]) => V.S([u * m.a1[0] + v * m.a2[0], u * m.a1[1] + v * m.a2[1]]));
      ctx.beginPath(); p.forEach((q, i) => i ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1])); ctx.closePath(); ctx.fill(); }
  } else {
    const L = 0.11 * Math.hypot(...m.a1) * V.sc * amp, col = Aids.st.arrows ? Aids.colIP() : css('--ink');
    for (const c of G.cells) { const ux = (c.fu[0] * cph - c.fu[1] * sph) / G.mx, uy = (c.fv[0] * cph - c.fv[1] * sph) / G.mx, p = V.S([c.x, c.y]); Aids.arrow(ctx, p[0], p[1], p[0] + L * ux, p[1] - L * uy, col, 1.4); }
  }
  // host geometry and resonators
  ctx.strokeStyle = css('--ink-2'); ctx.globalAlpha = 0.55; ctx.lineWidth = 1;
  for (const [ti, tj] of tiles) { const R = [ti * m.a1[0] + tj * m.a2[0], ti * m.a1[1] + tj * m.a2[1]]; ctx.beginPath(); for (const [p0, p1] of m.lines) { const a = V.S([p0[0] + R[0], p0[1] + R[1]]), b = V.S([p1[0] + R[0], p1[1] + R[1]]); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); } ctx.stroke(); }
  ctx.globalAlpha = 1;
  let rmx = 0; for (const r of m.res) rmx = Math.max(rmx, Math.hypot(...r.abs));
  for (const [ti, tj] of tiles) { const R = [ti * m.a1[0] + tj * m.a2[0], ti * m.a1[1] + tj * m.a2[1]], t = -(m.k[0] * R[0] + m.k[1] * R[1]);
    for (const r of m.res) { const p = V.S([r.pos[0] + R[0], r.pos[1] + R[1]]);
      if (m.fam === 'flex') { const zr = r.abs[0] * Math.cos(t + phase) - r.abs[1] * Math.sin(t + phase), v = zr / Math.max(rmx, G.mx); ctx.fillStyle = divColor(v * amp); ctx.strokeStyle = css('--ink'); ctx.lineWidth = 1.2; const s = 6 + 5 * Math.min(1, Math.abs(v)); ctx.fillRect(p[0] - s, p[1] - s, 2 * s, 2 * s); ctx.strokeRect(p[0] - s, p[1] - s, 2 * s, 2 * s); }
      else { const c = Math.cos(t + phase), s = Math.sin(t + phase), ux = (r.abs[0] * c - r.abs[1] * s) / Math.max(rmx, G.mx), uy = (r.abs[2] * c - r.abs[3] * s) / Math.max(rmx, G.mx), L = 0.11 * Math.hypot(...m.a1) * V.sc * amp; ctx.fillStyle = css('--res'); ctx.fillRect(p[0] + L * ux - 5, p[1] - L * uy - 5, 10, 10); } } }
  ctx.fillStyle = css('--ink-2'); ctx.font = '12.5px system-ui, sans-serif'; ctx.fillText(m.fam === 'flex' ? `Top view: colour = deflection (warm up, cool down), squares = resonators; ${Aids.label()} cells` : `In-plane displacement (arrows), resonators (squares); ${Aids.label()} cells`, 10, h - 8);
}
// FE mode of a drawn / 3D model. m: {md, k:[kx,ky], vec:[re,im]}
function modeFE(cv, m, phase, amp, view) {
  const { ctx, w, h } = B.prep(cv), md = m.md, [vr, vi] = m.vec, flat = md.nodes.every(p => Math.abs(p[2] - md.nodes[0][2]) < 1e-9);
  const tiles = m.tiles || Aids.tiles(), tkey = m.tiles ? 'own' : Aids.key();
  const yaw = view.yaw, elev = view.elev, cy = Math.cos(yaw), sy = Math.sin(yaw), ce = Math.cos(elev), se = Math.sin(elev);
  const f = p => { const x = cy * p[0] - sy * p[1], y = sy * p[0] + cy * p[1]; return [x, y * ce + p[2] * se, -y * se + p[2] * ce]; };
  let umax = 0; for (let i = 0; i < md.N; i++) umax = Math.max(umax, Math.hypot(vr[6 * i], vi[6 * i], vr[6 * i + 1], vi[6 * i + 1], vr[6 * i + 2], vi[6 * i + 2]));
  md.springs.forEach((sp, j) => { const J = 6 * md.N + j; umax = Math.max(umax, Math.hypot(vr[J], vi[J])); });
  const scale = amp * 0.12 * md.Lref / (umax || 1);
  const disp = (id, R) => { const t = m.k[0] * R[0] + m.k[1] * R[1] + phase, c = Math.cos(t), s = Math.sin(t); return [0, 1, 2].map(a => scale * (vr[6 * id + a] * c - vi[6 * id + a] * s)); };
  // draw every element where it was drawn (Xd), not at the folded node position, so each cell looks exactly like the drawing;
  // the lattice shift between the drawn point and its folded node sets the Bloch phase of that point
  const shiftOf = (id, P) => [P[0] - md.nodes[id][0], P[1] - md.nodes[id][1]];
  if (!md.drawPos) { md.drawPos = new Array(md.N); for (const e of md.elems) { const Xd = e.Xd || e.X; e.nd.forEach(([id], a) => { if (!md.drawPos[id]) md.drawPos[id] = Xd[a]; }); } }
  const items = [];
  for (const [ti, tj] of tiles) { const T = [ti * md.a1[0] + tj * md.a2[0], ti * md.a1[1] + tj * md.a2[1]];
    for (const e of md.elems) { const Xd = e.Xd || e.X;
      items.push({ ghost: true, X: Xd.map(P => [P[0] + T[0], P[1] + T[1], P[2] + (e.dz || 0)]), e });
      const X = e.nd.map(([id], a) => { const P = Xd[a], sh = shiftOf(id, P), u = disp(id, [T[0] + sh[0], T[1] + sh[1]]); return [P[0] + T[0] + u[0], P[1] + T[1] + u[1], P[2] + (e.dz || 0) + u[2]]; }); items.push({ e, X, main: !ti && !tj }); }
    md.springs.forEach((sp, j) => { const P = md.drawPos[sp.node] || md.nodes[sp.node], sh = shiftOf(sp.node, P), R = [T[0] + sh[0], T[1] + sh[1]], J = 6 * md.N + j, t = m.k[0] * R[0] + m.k[1] * R[1] + phase, c = Math.cos(t), s = Math.sin(t), q = scale * (vr[J] * c - vi[J] * s), u = disp(sp.node, R);
      items.push({ mass: true, X: [[P[0] + T[0] + u[0] + sp.n[0] * q, P[1] + T[1] + u[1] + sp.n[1] * q, P[2] + u[2] + sp.n[2] * q]], main: !ti && !tj }); }); }
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  if (!m.box) { for (const it of items) for (const p of it.X) { const q = f(p); x0 = Math.min(x0, q[0]); x1 = Math.max(x1, q[0]); y0 = Math.min(y0, q[1]); y1 = Math.max(y1, q[1]); } m.box = { x0, x1, y0, y1, yaw, elev, tkey }; }
  ({ x0, x1, y0, y1 } = m.box); if (m.box.yaw !== yaw || m.box.elev !== elev || m.box.tkey !== tkey) { m.box = null; return modeFE(cv, m, phase, amp, view); }
  const pad = 0.08 * Math.max(x1 - x0, y1 - y0), sc = Math.min((w - 30) / (x1 - x0 + 2 * pad), (h - 40) / (y1 - y0 + 2 * pad)), ox = w / 2 - sc * (x0 + x1) / 2, oy = (h - 20) / 2 + sc * (y0 + y1) / 2;
  const P = p => { const q = f(p); return [ox + sc * q[0], oy - sc * q[1], q[2]]; };
  items.forEach(it => { it.z = it.X.reduce((s, p) => s + P(p)[2], 0) / it.X.length; }); items.sort((a, b) => (a.ghost ? -1e9 : a.z) - (b.ghost ? -1e9 : b.z));
  const host = css('--host'), res = css('--res'), panel = css('--panel');
  for (const it of items) {
    if (it.ghost) { ctx.strokeStyle = css('--grid'); ctx.lineWidth = 1; ctx.beginPath(); it.X.forEach((p, i) => { const q = P(p); i ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1]); }); if (it.X.length > 2) ctx.closePath(); ctx.stroke(); continue; }
    if (it.mass) { const q = P(it.X[0]); ctx.fillStyle = it.main ? res : B.mix(res, panel, 0.5); ctx.fillRect(q[0] - 5, q[1] - 5, 10, 10); continue; }
    const col0 = it.e.role === 'res' ? res : host, col = it.main ? col0 : B.mix(col0, panel, 0.5);
    if (it.e.type === 'beam') { const a = P(it.X[0]), b = P(it.X[1]); ctx.strokeStyle = col; ctx.lineWidth = Math.max(1.5, Math.min(16, (it.e.draw ? it.e.draw.b : 0.001) * sc)); ctx.lineCap = 'round'; ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke(); }
    else { const Q = it.X.map(P); ctx.fillStyle = B.mix(panel, col0, it.main ? 0.6 : 0.3); ctx.strokeStyle = B.mix(panel, col0, it.main ? 0.9 : 0.45); ctx.lineWidth = 0.5; ctx.beginPath(); ctx.moveTo(Q[0][0], Q[0][1]); ctx.lineTo(Q[1][0], Q[1][1]); ctx.lineTo(Q[2][0], Q[2][1]); ctx.closePath(); ctx.fill(); ctx.stroke(); }
  }
  // motion arrows from the rest position: red up / down (out-of-plane), blue sideways (in-plane)
  if (Aids.st.arrows) {
    const cellL = Math.hypot(md.a1[0], md.a1[1]), many = tiles.length > 25, ak = (many ? 'm' : 'f') + md.N;
    if (!m._arw || m._arw.k !== ak) {
      const mag = i => Math.hypot(vr[6 * i], vi[6 * i], vr[6 * i + 1], vi[6 * i + 1], vr[6 * i + 2], vi[6 * i + 2]), mags = Array.from({ length: md.N }, (_, i) => mag(i));
      const ids = Aids.pickNodes(Array.from({ length: md.N }, (_, i) => md.drawPos[i] || null), mags, cellL / (many ? 2 : 4), md.springs.map(sp => sp.node));
      const ref = Math.max(1e-30, ...ids.map(i => mags[i])), cz = i => Math.hypot(vr[6 * i + 2], vi[6 * i + 2]), cxy = i => Math.hypot(vr[6 * i], vi[6 * i], vr[6 * i + 1], vi[6 * i + 1]);
      m._arw = { k: ak, ids, z: ids.map(i => cz(i) > 0.15 * ref), xy: ids.map(i => cxy(i) > 0.15 * ref) }; }
    const A = m._arw, cO = Aids.colOOP(), cI = Aids.colIP(), lw = many ? 1.2 : 1.7; let anyZ = false, anyXY = false;
    for (const [ti, tj] of tiles) { const T = [ti * md.a1[0] + tj * md.a2[0], ti * md.a1[1] + tj * md.a2[1]];
      A.ids.forEach((id, n) => { const P0 = md.drawPos[id]; if (!P0) return; const sh = shiftOf(id, P0), u = disp(id, [T[0] + sh[0], T[1] + sh[1]]), r = [P0[0] + T[0], P0[1] + T[1], P0[2]], a = P(r);
        if (A.z[n]) { const b = P([r[0], r[1], r[2] + u[2]]); Aids.arrow(ctx, a[0], a[1], b[0], b[1], cO, lw); anyZ = true; }
        if (A.xy[n]) { const b = P([r[0] + 2.5 * u[0], r[1] + 2.5 * u[1], r[2]]); Aids.arrow(ctx, a[0], a[1], b[0], b[1], cI, lw); anyXY = true; } }); }
    const kinds = []; if (anyZ) kinds.push('oop'); if (anyXY) kinds.push('ip'); if (kinds.length) Aids.arrowKey(ctx, w, kinds);
  }
  ctx.fillStyle = css('--ink-2'); ctx.font = '12.5px system-ui, sans-serif'; ctx.textAlign = 'left'; ctx.fillText(m.caption || `Full FE mode, ${Aids.label()} cells, rest shape in grey; drag to rotate`, 10, h - 8);
}
// animation registry
const live = new Map();
function register(id, fn) { live.set(id, fn); }
let phase = 0;
function loop() { phase = (phase + 0.07) % TAU; for (const [id, fn] of live) { const el = document.getElementById(id); if (el && el.offsetParent !== null) fn(phase); } requestAnimationFrame(loop); }
requestAnimationFrame(loop);
// nearest band point to a click: sets = [{F, xs, tag}], X, Y mapping
function pickPoint(sets, X, Y, mx, my) { let best = null; for (const s of sets) s.F.forEach((row, i) => row.forEach((f, b) => { const d = Math.hypot(X(s.xs[i]) - mx, Y(f) - my); if (d < 12 && (!best || d < best.d)) best = { d, tag: s.tag, i, b, f }; })); return best; }
return { complexPanel, minAtt, modePWE, modeFE, register, pickPoint, ATT_MAX };
})();

// ================= Shared solver + results for the homogenised (PWE) and full (FEM) models =================
const Res = (() => {
const B = window.MMW, TAU = 2 * Math.PI, css = B.css;
const FAMN = { flex: 'Flexural (out-of-plane)', ip: 'In-plane', all: 'All waves' };
function polarBands(F, Z, want) { const rows = F.map((r, i) => r.filter((f, b) => (Z[i][b] >= 0.5) === want)), n = Math.min(...rows.map(r => r.length)); return n > 0 ? rows.map(r => r.slice(0, n)) : null; }
// inp: {meth:'pwe'|'fem'|'both', Lat, path, plate, res, md (FE model in drawn format), NF, NFc, nfc, nbFE, fmin, fmax, r0, extraMass, ell, label}
async function solve(inp, status) {
  const R = { inp, pwe: null, fe: null }, path = inp.path, doP = inp.meth !== 'fem' && inp.plate, doF = inp.meth !== 'pwe' && inp.md;
  if (inp.meth !== 'fem' && !inp.plate) throw new Error('The plane-wave method needs a homogenised plate.');
  if (inp.meth !== 'pwe' && !inp.md) throw new Error('No FE model is available for this geometry.');
  const fl = (inp.res || []).flatMap(r => (r.use || []).map(q => q.f)), fr = fl.length ? Math.min(...fl) : 0;
  let steps = (doP ? 2 : 0) + (doF ? 1 : 0) + 1, step = 0; const tick = async (txt, u) => { status(txt, (step + u) / steps); await B.tick(); };
  if (doP) { R.pwe = {};
    for (const fam of ['flex', 'ip']) {
      const P = Hom.pweSetup(inp.Lat, inp.plate, inp.res, { fam, NF: inp.NF, r0: inp.r0, extraMass: inp.extraMass || 0 }), nb = Math.min(P.modes.length + 10, P.nw + P.modes.length), F = [], Fr = [];
      for (let i = 0; i < path.ks.length; i++) { const r = Hom.pweAt(P, path.ks[i][0], path.ks[i][1], nb, true); F.push(r.f); Fr.push(r.frac); if (i % 2 === 0) await tick(`Plane waves, ${fam === 'flex' ? 'flexural' : 'in-plane'} (${P.G.length} plane waves + ${P.modes.length} resonator modes): k-point ${i + 1} of ${path.ks.length}`, i / path.ks.length); }
      R.pwe[fam] = { F, Fr, setup: P, nb, top: Math.min(...F.map(r => r[r.length - 1])) };
      step++; }
  }
  if (doF) {
    const md = inp.md, fref = inp.fmax || (fr ? 3 * fr : 3000), F = [], Z = [], Rs = [];
    for (let i = 0; i < path.ks.length; i++) { const r = Drw.blochAt(md, path.ks[i][0], path.ks[i][1], inp.nbFE, fref * 0.5); F.push(r.f); Z.push(r.z); Rs.push(r.res); await tick(`Full Bloch FE (${6 * md.N + md.springs.length} DOF): k-point ${i + 1} of ${path.ks.length}`, (i + 1) / path.ks.length); }
    const zs = new Set(md.nodes.map(p => Math.round(p[2] * 1e7))), flat = zs.size === 1 && !md.hasOffset;
    R.fe = { F, Z, Rs, top: Math.min(...F.map(r => r[r.length - 1])), flat, fref };
    step++;
  }
  // frequency range
  const tops = []; if (R.pwe) tops.push(R.pwe.flex.top, R.pwe.ip.top); if (R.fe) tops.push(R.fe.top);
  const topMin = Math.min(...tops); R.fmax = inp.fmax || Math.min(topMin * 0.98, fr > 0 ? 3 * fr : topMin * 0.98); R.fmin = Math.min(inp.fmin || 0, R.fmax * 0.99);
  const ell = inp.ell, cut = g => g.lo < R.fmax && g.hi > R.fmin;
  if (R.pwe) for (const fam of ['flex', 'ip']) { const x = R.pwe[fam]; x.gaps = Hom.gapsOf(x.F).filter(g => g.lo < x.top).filter(cut).map(g => Object.assign(g, { ratio: Hom.plateWavelength(inp.plate, 0.5 * (g.lo + g.hi), fam) / ell })); }
  if (R.pwe) { R.pwe.full = []; for (const a of R.pwe.flex.gaps) for (const b of R.pwe.ip.gaps) { const lo = Math.max(a.lo, b.lo), hi = Math.min(a.hi, b.hi); if (hi > lo * 1.001) R.pwe.full.push({ lo, hi, ratio: Math.min(a.ratio, b.ratio) }); } }
  if (R.fe) { const x = R.fe; x.gaps = Hom.gapsOf(x.F).filter(g => g.lo < x.top).filter(cut); for (const [fam, want] of [['flex', true], ['ip', false]]) { const Fb = polarBands(x.F, x.Z, want); x[fam] = Fb ? Hom.gapsOf(Fb).filter(g => g.lo < x.top).filter(cut) : []; } }
  // complex bands
  const nf = inp.nfc, fs = []; for (let i = 1; i <= nf; i++) fs.push(R.fmin + (R.fmax - R.fmin) * i / nf);
  if (R.pwe) for (const fam of ['flex', 'ip']) { const P = Hom.pweSetup(inp.Lat, inp.plate, inp.res, { fam, NF: inp.NFc, r0: inp.r0, extraMass: inp.extraMass || 0 }), rows = [];
    for (const [i, f] of fs.entries()) { rows.push({ f, ev: Hom.pweComplexAt(P, inp.Lat, f) }); if (i % 4 === 0) await tick(`Complex band, plane waves (${fam === 'flex' ? 'flexural' : 'in-plane'}): ${i + 1} of ${nf} frequencies`, (fam === 'ip' ? 0.5 : 0) + 0.25 * i / nf); }
    R.pwe[fam].cplx = rows; }
  if (R.fe) { const md = R.fe && inp.md, sets = R.fe.flat ? { flex: i => i >= 6 * md.N ? Math.abs(md.springs[i - 6 * md.N].n[2]) > 0.5 : [2, 3, 4].includes(i % 6), ip: i => i >= 6 * md.N ? Math.abs(md.springs[i - 6 * md.N].n[2]) <= 0.5 : [0, 1, 5].includes(i % 6) } : { all: null };
    R.fe.cplx = {};
    for (const key of Object.keys(sets)) { const W = Drw.wfeModel(md, sets[key]), rows = [];
      for (const [i, f] of fs.entries()) { rows.push({ f, ev: Drw.wfeAt(W, f) }); if (i % 4 === 0) await tick(`Complex band, FE (${key === 'all' ? 'all waves' : key === 'flex' ? 'flexural' : 'in-plane'}): ${i + 1} of ${nf} frequencies`, 0.5 + 0.25 * i / nf); }
      R.fe.cplx[key] = rows; } }
  return R;
}
const M = { l: 58, r: 14, t: 16, b: 36 };
function bands(cv, R, view, sel) {
  const { ctx, w, h } = B.prep(cv), ink2 = css('--ink-2');
  if (!R) { ctx.fillStyle = ink2; ctx.font = '15px system-ui, sans-serif'; ctx.fillText('Run the model to compute the bands.', 20, 40); return null; }
  const path = R.inp.path, xs = path.xs, xmax = xs[xs.length - 1], [fmin, fmax] = Aids.fz(cv._zoom || cv.id, R.fmin, R.fmax), U = B.fUnit(fmax);
  const X = u => M.l + u / xmax * (w - M.l - M.r), Y = f => h - M.b - (f - fmin) / (fmax - fmin) * (h - M.t - M.b);
  const gaps = gapList(R, view);
  for (const g of gaps.fe.concat(gaps.fe.length ? [] : gaps.pwe)) { if (g.hi < fmin || g.lo > fmax) continue; const y1 = Y(Math.min(g.hi, fmax)), y0 = Y(Math.max(g.lo, fmin)); ctx.fillStyle = css('--gap-tint'); ctx.fillRect(M.l, y1, w - M.l - M.r, y0 - y1); ctx.fillStyle = B.hatch(ctx, css('--gap')); ctx.globalAlpha = 0.45; ctx.fillRect(M.l, y1, w - M.l - M.r, y0 - y1); ctx.globalAlpha = 1; }
  ctx.font = '13px system-ui, sans-serif'; ctx.fillStyle = ink2; ctx.strokeStyle = css('--grid'); ctx.lineWidth = 1;
  for (const v of B.niceTicks(fmin / U.d, fmax / U.d, 6)) { const y = Y(v * U.d); ctx.beginPath(); ctx.moveTo(M.l, y); ctx.lineTo(w - M.r, y); ctx.stroke(); ctx.textAlign = 'right'; ctx.fillText(+v.toFixed(3), M.l - 7, y + 4); }
  ctx.textAlign = 'center'; path.ticks.forEach((t, i) => { const x = X(t); ctx.strokeStyle = css('--rule'); ctx.beginPath(); ctx.moveTo(x, M.t); ctx.lineTo(x, h - M.b); ctx.stroke(); ctx.fillStyle = css('--ink'); ctx.font = '600 14px system-ui, sans-serif'; ctx.fillText(path.labels[i], x, h - M.b + 18); });
  ctx.save(); ctx.translate(15, (M.t + h - M.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillStyle = ink2; ctx.font = '13px system-ui, sans-serif'; ctx.fillText(`Frequency (${U.u})`, 0, 0); ctx.restore();
  ctx.fillStyle = ink2; ctx.font = '13px system-ui, sans-serif'; ctx.fillText('Wave vector along the zone path', (M.l + w - M.r) / 2, h - 4);
  ctx.save(); ctx.beginPath(); ctx.rect(M.l, M.t - 4, w - M.l - M.r, h - M.t - M.b + 4); ctx.clip();
  const pickSets = [], fam = view === 'all' ? null : view, host = css('--host'), res = css('--res'), edge = css('--edge');
  if (R.pwe) for (const f of fam ? [fam] : ['flex', 'ip']) { const x = R.pwe[f], col = f === 'flex' ? res : edge;
    ctx.strokeStyle = col; ctx.globalAlpha = 0.55; ctx.lineWidth = 1.2; ctx.setLineDash([5, 4]);
    for (let b = 0; b < x.nb; b++) { ctx.beginPath(); xs.forEach((u, i) => { const y = Y(x.F[i][b]); i ? ctx.lineTo(X(u), y) : ctx.moveTo(X(u), y); }); ctx.stroke(); }
    ctx.setLineDash([]); ctx.globalAlpha = 1;
    Aids.drawBranches(ctx, X, Y, xs, x.F, path.ks, null, f === 'flex' ? 1 : 2, { fmax, tag: fam || R.fe ? null : (f === 'flex' ? 'out-of-plane' : 'in-plane'), noLabels: !!R.fe }); ctx.fillStyle = col;
    for (let b = 0; b < x.nb; b++) xs.forEach((u, i) => { const fq = x.F[i][b]; if (fq > fmax * 1.02 || fq < fmin * 0.98) return; ctx.fillRect(X(u) - 2.2, Y(fq) - 2.2, 4.4, 4.4); });
    pickSets.push({ F: x.F, xs, tag: 'pwe-' + f }); }
  if (R.fe) { const x = R.fe, keep = (i, b) => !fam || (fam === 'flex' ? x.Z[i][b] >= 0.5 : x.Z[i][b] < 0.5), Fk = x.F.map((r, i) => r.map((f, b) => keep(i, b) ? f : NaN));
    Aids.drawBranches(ctx, X, Y, xs, Fk, path.ks, null, !fam ? 3 : fam === 'flex' ? 1 : 2, { fmax });
    for (let b = 0; b < x.F[0].length; b++) xs.forEach((u, i) => { const fq = Fk[i][b]; if (!(fq <= fmax * 1.02) || fq < fmin * 0.98) return; ctx.fillStyle = B.mix(host, res, (x.Rs[i][b] || 0) * 0.85); ctx.beginPath(); ctx.arc(X(u), Y(fq), 2.7, 0, TAU); ctx.fill(); });
    pickSets.push({ F: Fk, xs, tag: 'fe' }); }
  tagRes(R); { const gl = gapList(R, view); Aids.drawGapLabels(ctx, gl.fe.length ? gl.fe : gl.pwe, Y, w - M.r, fmin, fmax, view === 'all' && !!(R.fe && R.fe.flat)); }
  if (sel) { const set = sel.tag === 'fe' ? R.fe.F : R.pwe[sel.tag.slice(4)].F, fq = set[sel.i][sel.b]; ctx.strokeStyle = css('--ink'); ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(X(xs[sel.i]), Y(fq), 7, 0, TAU); ctx.stroke(); }
  ctx.restore();
  return { X, Y, pickSets };
}
// which resonator mode opens each FE gap (energy share of the band below it), see Aids.tagGaps
function tagRes(R) {
  if (!R || !R.fe || !R.fe.Rs) return; const x = R.fe, md = R.inp.md, hasRes = md.elems.some(e => e.role === 'res') || (md.pmass || []).some(q => q.role === 'res') || md.springs.length > 0;
  const mean = (Rs, b) => Rs.reduce((a, r) => a + (r[b] || 0), 0) / Rs.length;
  Aids.tagGaps(x.gaps, g => hasRes && g.below ? mean(x.Rs, g.below - 1) : null);
  for (const [fam, want] of [['flex', true], ['ip', false]]) { if (!x[fam]) continue; const Rp = x.Rs.map((r, i) => r.filter((_, b) => (x.Z[i][b] >= 0.5) === want)); Aids.tagGaps(x[fam], g => hasRes && g.below ? mean(Rp, g.below - 1) : null); }
}
function gapList(R, view) {
  const pick = (o, k) => (o ? (k === 'all' ? (o.gaps || o.full || []) : o[k] || []) : []);
  if (!R) return { fe: [], pwe: [] };
  const fe = R.fe ? (view === 'all' ? R.fe.gaps : R.fe[view]) : [], pwe = R.pwe ? (view === 'all' ? R.pwe.full : R.pwe[view].gaps) : [];
  return { fe: fe || [], pwe: pwe || [] };
}
function complex(cv, R, view) {
  if (!R) return Viz.complexPanel(cv, { sets: [], empty: 'Runs together with the bands.' });
  const sets = [];
  if (R.fe) { const c = R.fe.cplx; const rows = c.all || (view === 'all' ? null : c[view]); if (rows) sets.push({ rows, color: css('--host'), shape: 'dot' }); else if (view === 'all') { sets.push({ rows: c.flex, color: css('--host'), shape: 'dot' }, { rows: c.ip, color: css('--host'), shape: 'dot' }); } }
  if (R.pwe) for (const f of view === 'all' ? ['flex', 'ip'] : [view]) sets.push({ rows: R.pwe[f].cplx, color: f === 'flex' ? css('--res') : css('--edge'), shape: 'square' });
  const g = gapList(R, view); const [zl, zh] = Aids.fz(cv._zoom || cv.id, R.fmin, R.fmax); return Viz.complexPanel(cv, { sets, fmin: zl, fmax: zh, gaps: g.fe.length ? g.fe : g.pwe });
}
function table(el, R) {
  if (!R) return; tagRes(R);
  const ov = (a, b) => Math.min(a.hi, b.hi) - Math.max(a.lo, b.lo), match = (G, g) => { let best = null; for (const q of G || []) if (ov(q, g) > 0 && (!best || ov(q, g) > ov(best, g))) best = q; return best; };
  const verdict = r => r >= 4 ? '<span class="ok">long-wave, reliable</span>' : r >= 2 ? '<span class="warn">use with care</span>' : '<span class="err">Bragg range: use FEM</span>';
  const rng = g => `${B.fmtF(g.lo)} – ${B.fmtF(g.hi)}`, wid = g => `${(200 * (g.hi - g.lo) / (g.hi + g.lo)).toFixed(1)} %`;
  const rows = [], both = R.fe && R.pwe;
  for (const fam of ['flex', 'ip', 'all']) {
    const fe = R.fe ? (fam === 'all' ? R.fe.gaps : R.fe[fam]) : [], pw = R.pwe ? (fam === 'all' ? R.pwe.full : R.pwe[fam].gaps) : [], used = new Set(), lab = fam === 'all' ? '<b>Complete (all waves)</b>' : FAMN[fam];
    for (const g of fe) { const q = match(pw, g), op = fam === 'all' && R.fe.flat ? '' : Aids.gapShort(g); if (q) used.add(q); rows.push(`<tr><td>${lab}</td><td>${rng(g)}${op ? `<br><span class="sub">${op}</span>` : ''}</td><td>${wid(g)}</td>${R.pwe ? `<td>${q ? rng(q) : 'no gap'}</td><td>${q ? verdict(q.ratio) : ''}</td>` : ''}</tr>`); }
    for (const q of pw) if (!used.has(q)) rows.push(`<tr><td>${lab}</td>${R.fe ? '<td>no gap</td><td></td>' : `<td>${rng(q)}</td><td>${wid(q)}</td>`}${R.fe ? `<td>${rng(q)}</td><td>${verdict(q.ratio)}</td>` : `<td>${verdict(q.ratio)}</td>`}</tr>`);
  }
  const head = both ? '<th>Waves</th><th>FEM (full Bloch FE)</th><th>Width</th><th>PWE (homogenised plate)</th><th>PWE validity</th>' : R.fe ? '<th>Waves</th><th>FEM (full Bloch FE)</th><th>Width</th>' : '<th>Waves</th><th>PWE (homogenised plate)</th><th>Width</th><th>Validity</th>';
  el.innerHTML = rows.length ? `<table><thead><tr>${head}</tr></thead><tbody>${rows.join('')}</tbody></table>` : `<p class="empty">No band gap between ${B.fmtF(R.fmin)} and ${B.fmtF(R.fmax)}.</p>`;
  const notes = [];
  if (R.fe) notes.push(R.fe.flat ? 'The FE cell lies in one plane, so flexural and in-plane waves are exactly separate, in the bands and in the complex band.' : 'The FE cell is not flat: flexural and in-plane waves couple. The FE rows per polarisation sort modes by their dominant motion and are approximate; the complete-gap row is exact, and the FE complex band shows all waves together.');
  if (R.pwe) notes.push('PWE validity compares the plate wavelength at the gap centre with the cell size: long-wave local-resonance gaps are reliable, gaps near the Bragg range need FEM.');
  notes.push('Click a point on the bands to animate its mode.');
  el.innerHTML += `<p class="note" style="padding:0 12px 10px">${notes.join(' ')}</p>`;
}
// mode of a picked point. lines: host geometry (top view, metres) for the PWE plot
function mode(R, p, lines) {
  const k = R.inp.path.ks[p.i];
  if (p.tag === 'fe') { const r = Drw.blochAt(R.inp.md, k[0], k[1], Math.max(p.b + 1, R.fe.F[0].length), R.fe.fref * 0.5, true); return { kind: 'fe', md: R.inp.md, k, vec: r.vec[p.b], f: r.f[p.b], z: r.z[p.b] }; }
  const fam = p.tag.slice(4), x = R.pwe[fam], P = x.setup, r = Hom.pweAt(P, k[0], k[1], Math.max(p.b + 1, x.nb), false, true), [yr, yi] = r.vec[p.b];
  const Ng = P.G.length, nw = P.nw, res = [];
  P.res.forEach((Rs, ri) => {
    const base = [0, 1].slice(0, fam === 'flex' ? 1 : 2).map(o => { let re = 0, im = 0; for (const a of Rs.att) for (let g = 0; g < Ng; g++) { const t = -((k[0] + P.G[g][0]) * a[0] + (k[1] + P.G[g][1]) * a[1]), c = Math.cos(t), s = Math.sin(t); re += yr[o * Ng + g] * c - yi[o * Ng + g] * s; im += yr[o * Ng + g] * s + yi[o * Ng + g] * c; } return [re / Rs.att.length, im / Rs.att.length]; });
    P.modes.forEach((q, j) => { if (q.ri !== ri) return; base.forEach((bv, o) => { bv[0] += q.c[o] * yr[nw + j]; bv[1] += q.c[o] * yi[nw + j]; }); });
    const pos = Rs.pos ? [Rs.pos[0], Rs.pos[1]] : Rs.att[0];
    res.push({ pos, abs: fam === 'flex' ? base[0] : [base[0][0], base[0][1], base[1][0], base[1][1]] });
  });
  return { kind: 'pwe', fam, k, G: P.G, W: [yr, yi], a1: R.inp.Lat.a1, a2: R.inp.Lat.a2, lines, res, f: r.f[p.b] };
}
function hostLines(md) { // top-view segments of a drawn / 3D FE model: beams and free shell edges
  const segs = [], cnt = new Map(), key = (a, b) => { const k1 = a.map(v => Math.round(v * 1e7)).join(','), k2 = b.map(v => Math.round(v * 1e7)).join(','); return k1 < k2 ? k1 + '|' + k2 : k2 + '|' + k1; };
  for (const e of md.elems) { const X = e.Xd || e.X; if (e.type === 'beam') segs.push([X[0], X[1]]); else for (const [a, b] of [[0, 1], [1, 2], [2, 0]]) { const k = key(X[a], X[b]); cnt.set(k, (cnt.get(k) || { n: 0, s: [X[a], X[b]] })); cnt.get(k).n++; } }
  for (const v of cnt.values()) if (v.n === 1) segs.push(v.s);
  return segs.map(([a, b]) => [[a[0], a[1]], [b[0], b[1]]]);
}
// acoustic / optical name of a picked band point (view: 'flex' | 'ip' | 'all')
function branchName(R, p, view) {
  const ks = R.inp.path.ks;
  if (p.tag === 'fe') { const x = R.fe, fam = view === 'all' ? null : view, keep = (i, b) => !fam || (fam === 'flex' ? x.Z[i][b] >= 0.5 : x.Z[i][b] < 0.5); return Aids.branchOf(x.F, ks, keep, !fam ? 3 : fam === 'flex' ? 1 : 2, p.i, p.b); }
  const f = p.tag.slice(4); return Aids.branchOf(R.pwe[f].F, ks, null, f === 'flex' ? 1 : 2, p.i, p.b);
}
return { solve, bands, complex, table, mode, hostLines, gapList, branchName };
})();
