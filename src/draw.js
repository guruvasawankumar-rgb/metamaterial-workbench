// ================= Drawn unit cells: beams, shell walls, shell plates with holes, masses, springs =================
// Geometry is stored as coordinates (mm in the editor, metres here). Building a model meshes every part, welds
// coincident points, and identifies periodic images, so a part may cross the cell boundary freely.
// Solvers: full 3D Bloch FE (beam + flat shell elements, shift-invert block Lanczos) and the homogenised plate
// (homog.js) with resonator modes from the drawn resonator parts.
const Drw = (() => {
const Eg = (typeof Eng !== 'undefined') ? Eng : require('./engine.js');
const Hm = (typeof Hom !== 'undefined') ? Hom : require('./homog.js');
const TAU = 2 * Math.PI;
const zeros = (n, m) => { const A = new Array(n); for (let i = 0; i < n; i++) A[i] = new Float64Array(m); return A; };
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const nrm = a => Math.hypot(a[0], a[1], a[2]), unit = a => { const L = nrm(a); return [a[0] / L, a[1] / L, a[2] / L]; };

// ---------- flat shell triangle: CST membrane + drilling penalty, DKT bending; 6 DOF per node ----------
function shellTri(P, t, E, nu, rho) {
  const e1 = unit(sub(P[1], P[0])), nz = unit(cross(sub(P[1], P[0]), sub(P[2], P[0]))), e2 = cross(nz, e1);
  const loc = P.map(p => { const d = sub(p, P[0]); return [dot(d, e1), dot(d, e2)]; });
  const [x1, y1] = loc[0], [x2, y2] = loc[1], [x3, y3] = loc[2];
  const A2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1), A = A2 / 2;
  if (!(A > 1e-18)) throw new Error('A shell element has zero area: check the drawn plates and walls.');
  const Kl = zeros(18, 18), Ml = zeros(18, 18);
  // membrane CST, local DOFs u = 6i, v = 6i+1
  const b = [y2 - y3, y3 - y1, y1 - y2], c = [x3 - x2, x1 - x3, x2 - x1];
  const Cm = E * t / (1 - nu * nu), Dm = [[Cm, Cm * nu, 0], [Cm * nu, Cm, 0], [0, 0, Cm * (1 - nu) / 2]];
  const Bm = [new Float64Array(18), new Float64Array(18), new Float64Array(18)];
  for (let i = 0; i < 3; i++) { Bm[0][6 * i] = b[i] / A2; Bm[1][6 * i + 1] = c[i] / A2; Bm[2][6 * i] = c[i] / A2; Bm[2][6 * i + 1] = b[i] / A2; }
  const addBDB = (Bq, D, w) => { for (let r = 0; r < 3; r++) for (let s = 0; s < 3; s++) { const d = D[r][s] * w; if (!d) continue; const Br = Bq[r], Bs = Bq[s]; for (let i = 0; i < 18; i++) { if (!Br[i]) continue; const v = Br[i] * d; for (let j = 0; j < 18; j++) if (Bs[j]) Kl[i][j] += v * Bs[j]; } } };
  addBDB(Bm, Dm, A);
  // drilling: kd Σ (θz_i − ω)², ω = (v,x − u,y)/2
  const G = E / (2 * (1 + nu)), kd = G * t * A; // Hughes–Brezzi-type drilling penalty, γ = G
  for (let i = 0; i < 3; i++) { const g = new Float64Array(18); g[6 * i + 5] = 1; for (let j = 0; j < 3; j++) { g[6 * j + 1] -= 0.5 * b[j] / A2; g[6 * j] += 0.5 * c[j] / A2; }
    for (let p = 0; p < 18; p++) if (g[p]) for (let q = 0; q < 18; q++) if (g[q]) Kl[p][q] += kd / 3 * g[p] * g[q]; }
  // DKT bending (Batoz 1980), local DOFs w = 6i+2, θx = 6i+3 (= w,y), θy = 6i+4 (= −w,x)
  const Db = E * t ** 3 / (12 * (1 - nu * nu)), Dbm = [[Db, Db * nu, 0], [Db * nu, Db, 0], [0, 0, Db * (1 - nu) / 2]];
  const xij = (i, j) => loc[i][0] - loc[j][0], yij = (i, j) => loc[i][1] - loc[j][1];
  const ed = [[1, 2], [2, 0], [0, 1]], P_ = [], q_ = [], r_ = [], t_ = [];
  for (const [i, j] of ed) { const xx = xij(i, j), yy = yij(i, j), l2 = xx * xx + yy * yy; P_.push(-6 * xx / l2); q_.push(3 * xx * yy / l2); t_.push(-6 * yy / l2); r_.push(3 * yy * yy / l2); }
  // derivatives of Hx, Hy with respect to ξ, η (Batoz, Bathe & Ho); P4,P5,P6 for edges 23, 31, 12
  const [P4, P5, P6] = P_, [q4, q5, q6] = q_, [r4, r5, r6] = r_, [t4, t5, t6] = t_;
  const x31 = xij(2, 0), x12 = xij(0, 1), y31 = yij(2, 0), y12 = yij(0, 1);
  const gp = [[1 / 6, 1 / 6], [2 / 3, 1 / 6], [1 / 6, 2 / 3]];
  for (const [xi, eta] of gp) {
    const Hxx = [P6 * (1 - 2 * xi) + (P5 - P6) * eta, q6 * (1 - 2 * xi) - (q5 + q6) * eta, -4 + 6 * (xi + eta) + r6 * (1 - 2 * xi) - eta * (r5 + r6),
      -P6 * (1 - 2 * xi) + eta * (P4 + P6), q6 * (1 - 2 * xi) - eta * (q6 - q4), -2 + 6 * xi + r6 * (1 - 2 * xi) + eta * (r4 - r6),
      -eta * (P5 + P4), eta * (q4 - q5), -eta * (r5 - r4)];
    const Hyx = [t6 * (1 - 2 * xi) + eta * (t5 - t6), 1 + r6 * (1 - 2 * xi) - eta * (r5 + r6), -q6 * (1 - 2 * xi) + eta * (q5 + q6),
      -t6 * (1 - 2 * xi) + eta * (t4 + t6), -1 + r6 * (1 - 2 * xi) + eta * (r4 - r6), -q6 * (1 - 2 * xi) - eta * (q4 - q6),
      -eta * (t4 + t5), eta * (r4 - r5), -eta * (q4 - q5)];
    const Hxe = [-P5 * (1 - 2 * eta) - xi * (P6 - P5), q5 * (1 - 2 * eta) - xi * (q5 + q6), -4 + 6 * (xi + eta) + r5 * (1 - 2 * eta) - xi * (r5 + r6),
      xi * (P4 + P6), xi * (q4 - q6), -xi * (r6 - r4),
      P5 * (1 - 2 * eta) - xi * (P4 + P5), q5 * (1 - 2 * eta) + xi * (q4 - q5), -2 + 6 * eta + r5 * (1 - 2 * eta) + xi * (r4 - r5)];
    const Hye = [-t5 * (1 - 2 * eta) - xi * (t6 - t5), 1 + r5 * (1 - 2 * eta) - xi * (r5 + r6), -q5 * (1 - 2 * eta) + xi * (q5 + q6),
      xi * (t4 + t6), xi * (r4 - r6), -xi * (q4 - q6),
      t5 * (1 - 2 * eta) - xi * (t4 + t5), -1 + r5 * (1 - 2 * eta) + xi * (r4 - r5), -q5 * (1 - 2 * eta) - xi * (q4 - q5)];
    const Bb = [new Float64Array(18), new Float64Array(18), new Float64Array(18)];
    for (let k = 0; k < 9; k++) {
      const node = Math.floor(k / 3), dof = 6 * node + 2 + (k % 3);
      Bb[0][dof] = (y31 * Hxx[k] + y12 * Hxe[k]) / A2;
      Bb[1][dof] = (-x31 * Hyx[k] - x12 * Hye[k]) / A2;
      Bb[2][dof] = (-x31 * Hxx[k] - x12 * Hxe[k] + y31 * Hyx[k] + y12 * Hye[k]) / A2;
    }
    addBDB(Bb, Dbm, A / 3);
  }
  // lumped mass
  const m3 = rho * t * A / 3, J3 = rho * t ** 3 / 12 * A / 3, Jd = 1e-6 * m3 * A;
  for (let i = 0; i < 3; i++) { Ml[6 * i][6 * i] = Ml[6 * i + 1][6 * i + 1] = Ml[6 * i + 2][6 * i + 2] = m3; Ml[6 * i + 3][6 * i + 3] = Ml[6 * i + 4][6 * i + 4] = J3; Ml[6 * i + 5][6 * i + 5] = Jd; }
  // local -> global; DKT rotations: θx about e1, θy about e2 (positive by the right-hand rule) -> rotation vector θx e1 + θy e2
  const R = [[e1[0], e2[0], nz[0]], [e1[1], e2[1], nz[1]], [e1[2], e2[2], nz[2]]];
  const rot = Al => { const B = zeros(18, 18), C = zeros(18, 18);
    for (let i = 0; i < 18; i++) for (let jb = 0; jb < 6; jb++) for (let a = 0; a < 3; a++) { let v = 0; for (let q = 0; q < 3; q++) v += Al[i][3 * jb + q] * R[a][q]; B[i][3 * jb + a] = v; }
    for (let ib = 0; ib < 6; ib++) for (let a = 0; a < 3; a++) for (let j = 0; j < 18; j++) { let v = 0; for (let q = 0; q < 3; q++) v += R[a][q] * B[3 * ib + q][j]; C[3 * ib + a][j] = v; }
    for (let i = 0; i < 18; i++) for (let j = 0; j < i; j++) { const v = 0.5 * (C[i][j] + C[j][i]); C[i][j] = C[j][i] = v; }
    return C; };
  return { K: rot(Kl), M: rot(Ml), area: A, mass: rho * t * A };
}

// ---------- 2D geometry helpers ----------
function pointInPoly(p, poly) { let inside = false; for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) { const a = poly[i], b = poly[j]; if (((a[1] > p[1]) !== (b[1] > p[1])) && (p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0])) inside = !inside; } return inside; }
function segDist(p, a, b) { const dx = b[0] - a[0], dy = b[1] - a[1], L2 = dx * dx + dy * dy; let t = L2 > 0 ? ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2 : 0; t = Math.max(0, Math.min(1, t)); return { d: Math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy), t }; }
function polyArea(p) { let s = 0; for (let i = 0; i < p.length; i++) { const a = p[i], b = p[(i + 1) % p.length]; s += a[0] * b[1] - b[0] * a[1]; } return s / 2; }
// Bowyer–Watson Delaunay triangulation of 2D points
function delaunay(pts) {
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity; for (const p of pts) { x0 = Math.min(x0, p[0]); y0 = Math.min(y0, p[1]); x1 = Math.max(x1, p[0]); y1 = Math.max(y1, p[1]); }
  const d = Math.max(x1 - x0, y1 - y0) * 20 + 1e-9, cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
  const P = pts.concat([[cx - d, cy - d], [cx + d, cy - d], [cx, cy + d]]), n = pts.length;
  const circ = (a, b, c) => { const A = P[a], B = P[b], C = P[c], D = 2 * (A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]));
    const ux = ((A[0] ** 2 + A[1] ** 2) * (B[1] - C[1]) + (B[0] ** 2 + B[1] ** 2) * (C[1] - A[1]) + (C[0] ** 2 + C[1] ** 2) * (A[1] - B[1])) / D;
    const uy = ((A[0] ** 2 + A[1] ** 2) * (C[0] - B[0]) + (B[0] ** 2 + B[1] ** 2) * (A[0] - C[0]) + (C[0] ** 2 + C[1] ** 2) * (B[0] - A[0])) / D;
    return [ux, uy, (A[0] - ux) ** 2 + (A[1] - uy) ** 2]; };
  let tris = [[n, n + 1, n + 2, circ(n, n + 1, n + 2)]];
  for (let i = 0; i < n; i++) {
    const p = P[i], bad = [], keep = [];
    for (const t of tris) ((p[0] - t[3][0]) ** 2 + (p[1] - t[3][1]) ** 2 < t[3][2] * (1 + 1e-12) ? bad : keep).push(t);
    const edges = new Map();
    for (const t of bad) for (const [a, b] of [[t[0], t[1]], [t[1], t[2]], [t[2], t[0]]]) { const k = a < b ? a + ',' + b : b + ',' + a; edges.set(k, edges.has(k) ? null : [a, b]); }
    tris = keep;
    for (const e of edges.values()) if (e) { const c = circ(e[0], e[1], i); if (isFinite(c[2])) tris.push([e[0], e[1], i, c]); }
  }
  return tris.filter(t => t[0] < n && t[1] < n && t[2] < n).map(t => [t[0], t[1], t[2]]);
}
// mesh a polygon with holes at spacing h; boundary segments are subdivided uniformly (so shared edges conform);
// extra interior points (hard points) are kept as mesh nodes
function meshPlate(poly, holes, h, hard) {
  const rings = [poly, ...holes].map(r => polyArea(r) < 0 ? r.slice().reverse() : r.slice());
  let segs = []; // boundary segments with their sample points
  const bpts = [];
  rings.forEach(r => { for (let i = 0; i < r.length; i++) { const a = r[i], b = r[(i + 1) % r.length], L = Math.hypot(b[0] - a[0], b[1] - a[1]), n = Math.max(1, Math.ceil(L / h - 1e-9));
    for (let k = 0; k < n; k++) bpts.push([a[0] + (b[0] - a[0]) * k / n, a[1] + (b[1] - a[1]) * k / n]); segs.push([a, b]); } });
  const inside = p => pointInPoly(p, poly) && !holes.some(hh => pointInPoly(p, hh));
  const ipts = [];
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity; for (const p of poly) { x0 = Math.min(x0, p[0]); y0 = Math.min(y0, p[1]); x1 = Math.max(x1, p[0]); y1 = Math.max(y1, p[1]); }
  const dy = h * Math.sqrt(3) / 2;
  for (let j = 0, y = y0 + dy / 2; y < y1; j++, y += dy) for (let x = x0 + (j % 2 ? h / 2 : 0) + h / 4; x < x1; x += h) { const p = [x, y]; if (!inside(p)) continue; if (segs.some(s => segDist(p, s[0], s[1]).d < 0.55 * h)) continue; ipts.push(p); }
  for (const p of hard || []) { if (!inside(p)) continue; if (bpts.some(q => Math.hypot(q[0] - p[0], q[1] - p[1]) < 1e-9 * h)) continue; for (let k = ipts.length - 1; k >= 0; k--) if (Math.hypot(ipts[k][0] - p[0], ipts[k][1] - p[1]) < 0.5 * h) ipts.splice(k, 1); ipts.push(p); }
  // conforming Delaunay: split boundary segments that are missing from the triangulation
  let pts, tris;
  for (let it = 0; it < 8; it++) {
    pts = bpts.concat(ipts); tris = delaunay(pts).filter(t => { const c = [(pts[t[0]][0] + pts[t[1]][0] + pts[t[2]][0]) / 3, (pts[t[0]][1] + pts[t[1]][1] + pts[t[2]][1]) / 3]; return inside(c); });
    const have = new Set(); for (const t of tris) for (const [a, b] of [[t[0], t[1]], [t[1], t[2]], [t[2], t[0]]]) have.add(a < b ? a + ',' + b : b + ',' + a);
    // boundary neighbours in bpts order per ring
    const missing = []; let off = 0;
    rings.forEach(r => { const cnt = []; let tot = 0; for (let i = 0; i < r.length; i++) { const a = r[i], b = r[(i + 1) % r.length], n = Math.max(1, Math.ceil(Math.hypot(b[0] - a[0], b[1] - a[1]) / h - 1e-9)); cnt.push(n); tot += n; }
      for (let k = 0; k < tot; k++) { const a = off + k, b = off + (k + 1) % tot; if (!have.has(a < b ? a + ',' + b : b + ',' + a)) missing.push([a, b]); } off += tot; });
    if (!missing.length) break;
    // remove interior points too close to a missing edge
    for (const [a, b] of missing) { const A = pts[a], B = pts[b], L = Math.hypot(B[0] - A[0], B[1] - A[1]); for (let k = ipts.length - 1; k >= 0; k--) if (segDist(ipts[k], A, B).d < 0.6 * L) ipts.splice(k, 1); }
  }
  return { pts, tris };
}

// ---------- build a model from the drawing ----------
// ---------- shapes: standard figures (regular polygon, rectangle, circle) stored by their numbers ----------
// {kind:'shape', type:'ngon'|'rect'|'circle', c:[x,y], ang (deg), n, R (circumradius), sx, sy (rectangle), nseg (circle),
//  use:'plate'|'frame'|'walls'|'hole', t, z (plate), b, h (frame beams), t, H, z0 (walls), role, mat, E_GPa, nu, rho}
function shapePoly(s) {
  const r6 = v => +v.toFixed(6), a0 = (s.ang || 0) * Math.PI / 180, c = s.c || [0, 0], out = [];
  if (s.type === 'rect') { const hx = s.sx / 2, hy = s.sy / 2, ca = Math.cos(a0), sa = Math.sin(a0);
    for (const [x, y] of [[-hx, -hy], [hx, -hy], [hx, hy], [-hx, hy]]) out.push([r6(c[0] + ca * x - sa * y), r6(c[1] + sa * x + ca * y)]); return out; }
  const n = s.type === 'circle' ? Math.max(8, Math.round(s.nseg || 24)) : Math.max(3, Math.round(s.n || 6));
  for (let k = 0; k < n; k++) { const t = a0 + TAU * k / n; out.push([r6(c[0] + s.R * Math.cos(t)), r6(c[1] + s.R * Math.sin(t))]); }
  return out;
}
function shapeCheck(s, i) {
  const nm = `Shape ${i + 1}`;
  if (s.type === 'rect' ? !(s.sx > 0 && s.sy > 0) : !(s.R > 0)) throw new Error(`${nm} needs a positive size.`);
  if (s.type === 'ngon' && !(s.n >= 3)) throw new Error(`${nm} needs at least 3 sides.`);
  if (s.use === 'plate' && !(s.t > 0)) throw new Error(`${nm} (plate) needs a thickness.`);
  if (s.use === 'frame' && !(s.b > 0 && s.h > 0)) throw new Error(`${nm} (frame) needs beam width b and height h.`);
  if (s.use === 'walls' && !(s.t > 0 && s.H > 0)) throw new Error(`${nm} (shell walls) needs thickness t and height H.`);
}
// shapes -> primitive parts; holes are cut from the plate (or plate image) that contains their centre. src = index of the drawn part
function expandParts(parts, dr) {
  const out = [], holes = [], keep = ['role', 'mat', 'E_GPa', 'nu', 'rho', 'zref'];
  const base = (p, i) => { const q = { src: i }; for (const k of keep) if (p[k] != null && p[k] !== '') q[k] = p[k]; return q; };
  parts.forEach((p, i) => {
    if (p.kind !== 'shape') { out.push(Object.assign({}, p, { src: i })); return; }
    shapeCheck(p, i); const poly = shapePoly(p), z = p.z || 0;
    if (p.use === 'hole') { holes.push({ poly, c: p.c || [0, 0], i }); return; }
    if (p.use === 'plate') { out.push(Object.assign(base(p, i), { kind: 'plate', poly, holes: [], z, t: p.t })); return; }
    for (let k = 0; k < poly.length; k++) { const a = poly[k], b = poly[(k + 1) % poly.length];
      if (p.use === 'walls') out.push(Object.assign(base(p, i), { kind: 'wall', p: [a, b], t: p.t, H: p.H, z0: p.z0 || 0 }));
      else out.push(Object.assign(base(p, i), { kind: 'beam', p: [[a[0], a[1], z], [b[0], b[1], z]], b: p.b, h: p.h })); }
  });
  if (holes.length) {
    const a1 = dr ? dr.a1 : [0, 0], a2 = dr ? dr.a2 : [0, 0];
    const copied = new Set();
    for (const hl of holes) { let done = false;
      for (const pl of out) { if (pl.kind !== 'plate' || done) continue;
        for (let i = -1; i <= 1 && !done; i++) for (let j = -1; j <= 1 && !done; j++) { const sx = i * a1[0] + j * a2[0], sy = i * a1[1] + j * a2[1];
          if (pointInPoly([hl.c[0] - sx, hl.c[1] - sy], pl.poly)) { if (!copied.has(pl)) { pl.holes = (pl.holes || []).slice(); copied.add(pl); } pl.holes.push(hl.poly.map(q => [+(q[0] - sx).toFixed(6), +(q[1] - sy).toFixed(6)])); done = true; } } }
      if (!done) throw new Error(`Shape ${hl.i + 1} is a hole, but its centre is not inside any plate.`); }
  }
  return out;
}
// move / scale / rotate / mirror a drawn part in place (mm, degrees)
function transformPart(p, T) { // T: function [x,y] -> [x,y]
  const f = q => { const r = T([q[0], q[1]]); return q.length > 2 ? [+r[0].toFixed(6), +r[1].toFixed(6), q[2]] : [+r[0].toFixed(6), +r[1].toFixed(6)]; };
  if (p.kind === 'beam' || p.kind === 'wall') p.p = p.p.map(f);
  else if (p.kind === 'plate') { p.poly = p.poly.map(f); p.holes = (p.holes || []).map(r => r.map(f)); }
  else if (p.kind === 'shape') p.c = f(p.c || [0, 0]);
  else p.p = f(p.p);
}

// ---------- symmetry of a drawing: which point operations (about any centre) map the periodic pattern onto itself ----------
function symmetry(dr) {
  const parts = expandParts(dr.parts || [], dr), a1 = dr.a1, a2 = dr.a2, det = a1[0] * a2[1] - a1[1] * a2[0];
  const L = Math.max(Math.hypot(...a1), Math.hypot(...a2)), tol = 2e-5 * L, segs = [], pts = [];
  const tg = (p, extra) => [p.kind, p.role === 'res' ? 'r' : 'h', p.mat || '', p.E_GPa || '', p.rho || '', extra].join(':');
  for (const p of parts) {
    if (p.kind === 'beam') segs.push({ A: p.p[0], B: p.p[1], tag: tg(p, [p.b, p.h, p.p[0][2] || 0, p.p[1][2] || 0].join(',')) });
    else if (p.kind === 'wall') segs.push({ A: p.p[0], B: p.p[1], tag: tg(p, [p.t, p.H, p.z0 || 0].join(',')) });
    else if (p.kind === 'plate') for (const r of [p.poly, ...(p.holes || [])]) for (let i = 0; i < r.length; i++) segs.push({ A: r[i], B: r[(i + 1) % r.length], tag: tg(p, [p.t, p.z || 0, p.zref || ''].join(',')) });
    else if (p.kind === 'mass') pts.push({ P: p.p, tag: tg(p, [p.m, p.side, p.p[2] || 0].join(',')) });
    else if (p.kind === 'spring') pts.push({ P: p.p, tag: tg(p, [p.m, p.f, p.k, p.p[2] || 0].join(',')) });
  }
  // split every segment at the end points of the others (and their periodic images) so T-joints do not matter
  const ends = []; for (const s of segs) ends.push(s.A, s.B); for (const q of pts) ends.push(q.P);
  const imgs = []; for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) imgs.push([i * a1[0] + j * a2[0], i * a1[1] + j * a2[1]]);
  const pieces = [];
  for (const s of segs) { const dx = s.B[0] - s.A[0], dy = s.B[1] - s.A[1], L2 = dx * dx + dy * dy, ts = [0, 1]; if (!(L2 > 0)) continue;
    for (const e of ends) for (const sh of imgs) { const x = e[0] + sh[0] - s.A[0], y = e[1] + sh[1] - s.A[1], t = (x * dx + y * dy) / L2; if (t > 1e-6 && t < 1 - 1e-6 && Math.abs(x * dy - y * dx) / Math.sqrt(L2) < tol) ts.push(t); }
    ts.sort((p, q) => p - q); for (let k = 0; k < ts.length - 1; k++) if (ts[k + 1] - ts[k] > 1e-6) pieces.push({ A: [s.A[0] + dx * ts[k], s.A[1] + dy * ts[k]], B: [s.A[0] + dx * ts[k + 1], s.A[1] + dy * ts[k + 1]], tag: s.tag }); }
  const red = (x, y) => { let u = (x * a2[1] - y * a2[0]) / det, v = (a1[0] * y - a1[1] * x) / det; u -= Math.floor(u + 1e-7); v -= Math.floor(v + 1e-7); return [u * a1[0] + v * a2[0], u * a1[1] + v * a2[1]]; };
  const q = v => Math.round(v / tol), set = new Set();
  const segKey = (A, B) => { const m = red((A[0] + B[0]) / 2, (A[1] + B[1]) / 2); let dx = (B[0] - A[0]) / 2, dy = (B[1] - A[1]) / 2; if (dx < -tol || (Math.abs(dx) <= tol && dy < 0)) { dx = -dx; dy = -dy; } return [m, [dx, dy]]; };
  const nb3 = [-1, 0, 1], addSeg = (tag, m, d) => { for (const i of nb3) for (const j of nb3) for (const k of nb3) for (const l of nb3) set.add(`${tag}|${q(m[0]) + i},${q(m[1]) + j}|${q(d[0]) + k},${q(d[1]) + l}`); };
  for (const s of pieces) { const [m, d] = segKey(s.A, s.B); addSeg(s.tag, m, d); if (Math.abs(d[0]) <= 2 * tol) addSeg(s.tag, m, [-d[0], -d[1]]); }
  for (const p of pts) { const m = red(p.P[0], p.P[1]); for (const i of nb3) for (const j of nb3) set.add(`${p.tag}|${q(m[0]) + i},${q(m[1]) + j}`); }
  const has = (A, B, tag) => { const [m, d] = segKey(A, B); return set.has(`${tag}|${q(m[0])},${q(m[1])}|${q(d[0])},${q(d[1])}`); };
  const hasPt = (P, tag) => { const m = red(P[0], P[1]); return set.has(`${tag}|${q(m[0])},${q(m[1])}`); };
  const test = M => { // linear map M (2x2) about the origin, then any translation t
    const T = (p, t) => [M[0] * p[0] + M[1] * p[1] + t[0], M[2] * p[0] + M[3] * p[1] + t[1]];
    const ref = pieces.length ? pieces[0] : null, cands = [];
    if (ref) { const r0 = T([(ref.A[0] + ref.B[0]) / 2, (ref.A[1] + ref.B[1]) / 2], [0, 0]);
      for (const s of pieces) if (s.tag === ref.tag) cands.push([(s.A[0] + s.B[0]) / 2 - r0[0], (s.A[1] + s.B[1]) / 2 - r0[1]]); }
    else if (pts.length) { const r0 = T(pts[0].P, [0, 0]); for (const p of pts) if (p.tag === pts[0].tag) cands.push([p.P[0] - r0[0], p.P[1] - r0[1]]); }
    return cands.find(t => pieces.every(s => has(T(s.A, t), T(s.B, t), s.tag)) && pts.every(p => hasPt(T(p.P, t), p.tag))) || null;
  };
  const c6 = Math.cos(Math.PI / 3), s6 = Math.sin(Math.PI / 3);
  const c2 = test([-1, 0, 0, -1]);
  return { mx: !!test([-1, 0, 0, 1]), my: !!test([1, 0, 0, -1]), r90: !!test([0, -1, 1, 0]), r60: !!test([c6, -s6, s6, c6]), c2: !!c2, centre: c2 ? [c2[0] / 2, c2[1] / 2] : null };
}

// ---------- strip (supercell) for edge waves: N copies across, periodic along the edge, free on both long sides ----------
function stripDrawing(dr, N, dir) {
  const e = dir === 'a2' ? dr.a2 : dr.a1, c = dir === 'a2' ? dr.a1 : dr.a2, prim = expandParts(dr.parts || [], dr), parts = [];
  for (let i = 0; i < N; i++) for (const p of prim) { const q = JSON.parse(JSON.stringify(p)); q.src = p.src; q.cell = i; transformPart(q, ([x, y]) => [x + (i + 2) * c[0], y + (i + 2) * c[1]]); parts.push(q); } // offset by two cells so no node wraps across the strip
  const span = N + 4;
  return { a1: e.slice(), a2: [c[0] * span, c[1] * span], mesh: dr.mesh, parts, strip: { N, span, e, c } };
}
function buildStrip(dr, mat, N, dir, opt = {}) {
  const sd = stripDrawing(dr, N, dir), md = build(sd, mat, opt), a1 = md.a1, a2 = md.a2, det = a1[0] * a2[1] - a1[1] * a2[0], span = sd.strip.span;
  const depth = new Float64Array(md.N);
  md.nodes.forEach((X, i) => { let v = (a1[0] * X[1] - a1[1] * X[0]) / det * span; if (v > span - 2) v -= span; depth[i] = v; });
  md.Lref = Math.max(Math.hypot(...a1), Math.hypot(sd.strip.c[0], sd.strip.c[1]) / 1000);
  md.strip = { N, dir, depth, dmin: Math.min(...depth), dmax: Math.max(...depth) };
  return md;
}

// rigid offset of the element's mid-surface (or axis) by dz along z from its nodes: u_mid = u + θ × (0, 0, dz); A ← Tᵀ A T
function offsetMat(A, nn, dz) {
  for (let a = 0; a < nn; a++) { const o = 6 * a; for (const r of A) { r[o + 3] -= dz * r[o + 1]; r[o + 4] += dz * r[o]; } }
  for (let a = 0; a < nn; a++) { const o = 6 * a, R0 = A[o], R1 = A[o + 1], R3 = A[o + 3], R4 = A[o + 4]; for (let j = 0; j < R0.length; j++) { R3[j] -= dz * R1[j]; R4[j] += dz * R0[j]; } }
}
// drawing (mm, g, Hz): { a1, a2, mesh, parts:[ {kind:'beam', p:[[x,y,z],[x,y,z]], b, h, role}, {kind:'wall', p:[[x,y],[x,y]], z0, H, t, role},
//   {kind:'plate', poly:[[x,y]..], holes:[[[x,y]..]], z, t, role}, {kind:'mass', p:[x,y,z], m, side, role}, {kind:'spring', p:[x,y,z], m, f, dir} ], E?, rho? per part }
function build(dr, mat, opt = {}) {
  const mm = v => v / 1000, a1 = [mm(dr.a1[0]), mm(dr.a1[1])], a2 = [mm(dr.a2[0]), mm(dr.a2[1])];
  const det = a1[0] * a2[1] - a1[1] * a2[0]; if (!(Math.abs(det) > 1e-14)) throw new Error('The lattice vectors are parallel.');
  const h = mm(dr.mesh || 2), Lref = Math.max(Math.hypot(...a1), Math.hypot(...a2)), tol = 1e-6 * Lref, tim = opt.beam !== 'eb';
  const parts = expandParts(dr.parts || [], dr); if (!parts.some(p => p.role !== 'res' && ['beam', 'wall', 'plate'].includes(p.kind))) throw new Error('Draw at least one host beam, wall or plate.');
  // segments for T-junction splitting (beams, wall bases) and hard points
  const P3 = q => [mm(q[0]), mm(q[1]), mm(q[2] || 0)];
  const hardPts = [];
  for (const p of parts) { if (p.kind === 'beam') hardPts.push(P3(p.p[0]), P3(p.p[1])); if (p.kind === 'wall') { const z0 = mm(p.z0 || 0), z1 = z0 + mm(p.H); for (const q of p.p) hardPts.push([mm(q[0]), mm(q[1]), z0], [mm(q[0]), mm(q[1]), z1]); }
    if (p.kind === 'mass' || p.kind === 'spring') hardPts.push(P3(p.p)); if (p.kind === 'plate') for (const q of p.poly) hardPts.push([mm(q[0]), mm(q[1]), mm(p.z || 0)]); }
  const imgs = []; for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) imgs.push([i * a1[0] + j * a2[0], i * a1[1] + j * a2[1]]);
  // beams that cross other beams or plate outlines at the same height are joined there
  { const segs = [];
    for (const p of parts) { if (p.kind === 'beam' && Math.abs((p.p[0][2] || 0) - (p.p[1][2] || 0)) < 1e-9) segs.push([mm(p.p[0][0]), mm(p.p[0][1]), mm(p.p[1][0]), mm(p.p[1][1]), mm(p.p[0][2] || 0), 1]);
      if (p.kind === 'plate') for (const r of [p.poly, ...(p.holes || [])]) for (let i = 0; i < r.length; i++) { const a = r[i], b = r[(i + 1) % r.length]; segs.push([mm(a[0]), mm(a[1]), mm(b[0]), mm(b[1]), mm(p.z || 0), 0]); } }
    for (let i = 0; i < segs.length; i++) for (let j = i + 1; j < segs.length; j++) { const A = segs[i], Bs = segs[j]; if (!(A[5] || Bs[5]) || Math.abs(A[4] - Bs[4]) > tol) continue;
      for (const s of imgs) { const x1 = A[0], y1 = A[1], dx = A[2] - x1, dy = A[3] - y1, x3 = Bs[0] + s[0], y3 = Bs[1] + s[1], ex = Bs[2] + s[0] - x3, ey = Bs[3] + s[1] - y3, den = dx * ey - dy * ex;
        if (Math.abs(den) < 1e-12 * Lref * Lref) continue; const t = ((x3 - x1) * ey - (y3 - y1) * ex) / den, u = ((x3 - x1) * dy - (y3 - y1) * dx) / den;
        if (t > 1e-6 && t < 1 - 1e-6 && u > 1e-6 && u < 1 - 1e-6) hardPts.push([x1 + t * dx, y1 + t * dy, A[4]]); } } }
  // split a 3D segment at hard points lying on it (any periodic image)
  const splitSeg = (A, B) => { const ts = [0, 1], d = sub(B, A), L2 = dot(d, d);
    for (const q of hardPts) for (const s of imgs) { const Q = [q[0] + s[0], q[1] + s[1], q[2]], t = dot(sub(Q, A), d) / L2; if (t <= 1e-6 || t >= 1 - 1e-6) continue; const e = nrm(sub(Q, [A[0] + t * d[0], A[1] + t * d[1], A[2] + t * d[2]])); if (e < tol * 10) ts.push(t); }
    // merge split points closer than the joining tolerance, so no piece collapses to a point when nodes are welded
    const L = Math.sqrt(L2), dt = 50 * tol / L, srt = [...new Set(ts.map(t => +t.toFixed(10)))].sort((x, y) => x - y), out = [0];
    for (const t of srt) { if (t <= 0) continue; if (t >= 1) break; if (t - out[out.length - 1] > dt && 1 - t > dt) out.push(t); }
    out.push(1); return out; };
  const nodes = [], key = new Map(), raw = []; // raw element list with 3D coordinates
  const elems = [];
  // material of a part: its slot (host, resonator, membrane; by default from its role), then any value typed on the part
  // thickness reference: 'mid' (centred on z, default), 'bottom' (part sits on z), 'top' (part hangs below z); the part joins others at z
  const zoff = (p, th) => p.zref === 'bottom' ? th / 2 : p.zref === 'top' ? -th / 2 : 0;
  const slotOf = p => { const k = p.mat || (p.role === 'res' ? 'res' : 'host'); return k === 'mem' ? (mat.mem || mat.res || mat) : k === 'res' ? (mat.res || mat) : mat; };
  const SC = Object.assign({ hostE: 1, hostRho: 1, resE: 1, resRho: 1 }, opt.scale || {}), fE = p => p.role === 'res' ? SC.resE : SC.hostE, fR = p => p.role === 'res' ? SC.resRho : SC.hostRho;
  const props = p => { const m = slotOf(p); return { E: (p.E_GPa > 0 ? p.E_GPa * 1e9 : m.E) * fE(p), rho: (p.rho > 0 ? p.rho : m.rho) * fR(p), nu: p.nu > -1 && p.nu < 0.5 && p.nu !== '' && p.nu != null ? p.nu : m.nu }; };
  parts.forEach((p, pi) => {
    const role = p.role === 'res' ? 'res' : 'host', pr = props(p);
    if (p.kind === 'beam') {
      if (!(p.b > 0 && p.h > 0)) throw new Error(`Beam ${(p.src ?? pi) + 1} needs width b and height h.`);
      const A = P3(p.p[0]), B = P3(p.p[1]); if (nrm(sub(B, A)) < tol) throw new Error(`Beam ${(p.src ?? pi) + 1} has zero length: both of its ends are at the same point. Select it (Select tool) and delete it, or drag one end away.`);
      const ts = splitSeg(A, B), sec = { shape: 'rect', b: mm(p.b), h: mm(p.h), E: pr.E, rho: pr.rho, nu: pr.nu }, bp = Hm.beamProps(sec, mat);
      for (let k = 0; k < ts.length - 1; k++) { const p0 = A.map((v, c) => v + (B[c] - v) * ts[k]), p1 = A.map((v, c) => v + (B[c] - v) * ts[k + 1]), L = nrm(sub(p1, p0)), n = Math.max(opt.neMin || 2, Math.ceil(L / h - 1e-9));
        for (let s = 0; s < n; s++) raw.push({ type: 'beam', X: [p0.map((v, c) => v + (p1[c] - v) * s / n), p0.map((v, c) => v + (p1[c] - v) * (s + 1) / n)], bp, role, part: p.src ?? pi, draw: { b: sec.b, h: sec.h }, dz: zoff(p, mm(p.h)) }); }
    } else if (p.kind === 'wall') {
      if (!(p.t > 0 && p.H > 0)) throw new Error(`Wall ${(p.src ?? pi) + 1} needs thickness t and height H.`);
      const z0 = mm(p.z0 || 0), H = mm(p.H), A = [mm(p.p[0][0]), mm(p.p[0][1]), z0], B = [mm(p.p[1][0]), mm(p.p[1][1]), z0], ts = splitSeg(A, B), nH = Math.max(2, Math.ceil(H / h - 1e-9));
      for (let k = 0; k < ts.length - 1; k++) { const p0 = A.map((v, c) => v + (B[c] - v) * ts[k]), p1 = A.map((v, c) => v + (B[c] - v) * ts[k + 1]), nL = Math.max(1, Math.ceil(nrm(sub(p1, p0)) / h - 1e-9));
        const at = (i, j) => [p0[0] + (p1[0] - p0[0]) * i / nL, p0[1] + (p1[1] - p0[1]) * i / nL, z0 + H * j / nH];
        for (let i = 0; i < nL; i++) for (let j = 0; j < nH; j++) { const q = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)];
          const pair = (i + j) % 2 ? [[q[0], q[1], q[2]], [q[0], q[2], q[3]]] : [[q[0], q[1], q[3]], [q[1], q[2], q[3]]];
          for (const X of pair) raw.push({ type: 'shell', X, t: mm(p.t), pr, role, part: p.src ?? pi }); } }
    } else if (p.kind === 'plate') {
      if (!(p.t > 0)) throw new Error(`Plate ${(p.src ?? pi) + 1} needs a thickness.`);
      if (!(p.poly && p.poly.length >= 3)) throw new Error(`Plate ${(p.src ?? pi) + 1} needs at least three corners.`);
      const z = mm(p.z || 0), poly = p.poly.map(q => [mm(q[0]), mm(q[1])]), holes = (p.holes || []).filter(hh => hh.length >= 3).map(hh => hh.map(q => [mm(q[0]), mm(q[1])]));
      // insert hard points that lie on the outline (T-junctions) into the polygon rings
      const ins = ring => { const out = []; for (let i = 0; i < ring.length; i++) { const a = [ring[i][0], ring[i][1], z], b = [ring[(i + 1) % ring.length][0], ring[(i + 1) % ring.length][1], z], ts = splitSeg(a, b); for (const t of ts.slice(0, -1)) out.push([a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]); } return out; };
      const interior = []; for (const q of hardPts) for (const s of imgs) if (Math.abs(q[2] - z) < tol) interior.push([q[0] + s[0], q[1] + s[1]]);
      const m = meshPlate(ins(poly), holes.map(ins), h, interior);
      const dz = zoff(p, mm(p.t)); for (const t of m.tris) raw.push({ type: 'shell', X: t.map(i => [m.pts[i][0], m.pts[i][1], z]), t: mm(p.t), pr, role, part: p.src ?? pi, dz });
    }
  });
  // weld: canonical master position inside the cell, integer shift
  const frac = X => { const u = (X[0] * a2[1] - X[1] * a2[0]) / det, v = (a1[0] * X[1] - a1[1] * X[0]) / det; return [u, v]; };
  const nodeRef = X => { const [u, v] = frac(X); const s1 = Math.floor(u + 1e-7), s2 = Math.floor(v + 1e-7); const x = X[0] - s1 * a1[0] - s2 * a2[0], y = X[1] - s1 * a1[1] - s2 * a2[1];
    const k = `${Math.round(x / tol)},${Math.round(y / tol)},${Math.round(X[2] / tol)}`; let id = key.get(k);
    if (id == null) { // tolerant lookup
      const kx = Math.round(x / tol), ky = Math.round(y / tol), kz = Math.round(X[2] / tol);
      outer: for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) { const kk = `${kx + dx},${ky + dy},${kz + dz}`; if (key.has(kk)) { id = key.get(kk); break outer; } }
      if (id == null) { id = nodes.length; nodes.push([x, y, X[2]]); } key.set(k, id); }
    return [id, [s1, s2]]; };
  const skipped = new Set(), merged = new Set(), beamSeen = new Set(); let hasOffset = false;
  for (const r of raw) {
    const refs = r.X.map(nodeRef), base = refs[0][1], nd = refs.map(([id, s]) => [id, [s[0] - base[0], s[1] - base[1]]]);
    const X = nd.map(([id, s]) => [nodes[id][0] + s[0] * a1[0] + s[1] * a2[0], nodes[id][1] + s[0] * a1[1] + s[1] * a2[1], nodes[id][2]]);
    // skip a piece whose corners were welded into one point (two parts drawn almost on top of each other)
    const same = (i, j) => nd[i][0] === nd[j][0] && nd[i][1][0] === nd[j][1][0] && nd[i][1][1] === nd[j][1][1];
    if (r.type === 'beam' ? same(0, 1) : (same(0, 1) || same(1, 2) || same(0, 2))) { skipped.add(r.part); continue; }
    if (r.type === 'shell') { const e1 = sub(X[1], X[0]), e2 = sub(X[2], X[0]); if (nrm(cross(e1, e2)) < 1e-6 * nrm(e1) * nrm(e2)) { skipped.add(r.part); continue; } }
    if (r.type === 'beam') {
      // the same beam piece drawn twice (e.g. a full hexagon frame whose walls are shared with the neighbouring cells): keep it once
      const kk2 = (() => { const [i0, s0] = nd[0], [i1, s1] = nd[1], m0 = `${i1}:0,0|${i0}:${s0[0] - s1[0]},${s0[1] - s1[1]}`, m1 = `${i0}:0,0|${i1}:${s1[0] - s0[0]},${s1[1] - s0[1]}`; return m0 < m1 ? m0 : m1; })();
      const kk3 = kk2 + `|${r.bp.A}|${r.bp.Iy}|${r.bp.Iz}|${r.bp.E}|${r.bp.rho}|${r.dz || 0}`;
      if (beamSeen.has(kk3)) { merged.add(r.part); continue; } beamSeen.add(kk3); }
    let fe; try { if (r.type === 'beam') { fe = Hm.frame3D(X[0], X[1], r.bp, tim); fe.mass = r.bp.rho * r.bp.A * fe.L; } else fe = shellTri(X, r.t, r.pr.E, r.pr.nu, r.pr.rho); }
    catch (e) { throw new Error(`Part ${r.part + 1}: ${e.message} Two points of the drawing are almost on top of each other; move one of them, or snap it onto the other.`); }
    if (r.dz) { offsetMat(fe.K, r.X.length, r.dz); offsetMat(fe.M, r.X.length, r.dz); hasOffset = true; }
    elems.push({ type: r.type, nd, K: fe.K, M: fe.M, mass: fe.mass, role: r.role, part: r.part, X, Xd: r.X, draw: r.draw, t: r.t, dz: r.dz || 0 });
  }
  const pmass = [], springs = [];
  parts.forEach((p, pi) => {
    if (p.kind === 'mass') { if (!(p.m > 0)) throw new Error(`Mass ${(p.src ?? pi) + 1} needs a mass.`); const [id] = nodeRef(P3(p.p)); const m = p.m / 1000 * fR(p), s = mm(p.side || Math.cbrt(p.m / 1000 / slotOf(p).rho) * 1000); pmass.push({ node: id, m, J: m * s * s / 6, side: s, role: p.role === 'res' ? 'res' : 'host', part: p.src ?? pi }); }
    if (p.kind === 'spring') { if (!(p.m > 0 && (p.f > 0 || p.k > 0))) throw new Error(`Mass-spring ${(p.src ?? pi) + 1} needs a mass and f (or k).`); const [id] = nodeRef(P3(p.p)); const m0 = p.m / 1000, m = m0 * SC.resRho, k = (p.k > 0 ? p.k : (TAU * p.f) ** 2 * m0) * SC.resE; springs.push({ node: id, m, k, f: Math.sqrt(k / m) / TAU, n: unit(p.dir || [0, 0, 1]), part: p.src ?? pi }); }
  });
  let N = nodes.length; const used = new Uint8Array(N); for (const e of elems) for (const [id] of e.nd) used[id] = 1;
  for (const pm of pmass) if (!used[pm.node]) throw new Error('A point mass is not on any drawn part: place it on a beam, wall or plate node.');
  for (const sp of springs) if (!used[sp.node]) throw new Error('A mass-spring is not on any drawn part: place it on a beam, wall or plate node.');
  if (used.some(u => !u)) { // nodes left only by skipped pieces: renumber
    const map = new Int32Array(N).fill(-1), keep = []; for (let i = 0; i < N; i++) if (used[i]) { map[i] = keep.length; keep.push(nodes[i]); }
    for (const e of elems) e.nd = e.nd.map(([id, sh]) => [map[id], sh]); for (const q of pmass) q.node = map[q.node]; for (const q of springs) q.node = map[q.node];
    nodes.length = 0; nodes.push(...keep); N = nodes.length; }
  const zs = nodes.map(p => p[2]), Henv = Math.max(1e-6, Math.max(...zs) - Math.min(...zs) + Math.max(0, ...parts.filter(p => p.kind === 'beam').map(p => mm(p.h)), ...parts.filter(p => p.kind === 'plate').map(p => mm(p.t))));
  // parts that are not connected to any host part float freely (zero-frequency modes): report them
  const par = Int32Array.from({ length: N }, (_, i) => i), fnd = i => { while (par[i] !== i) { par[i] = par[par[i]]; i = par[i]; } return i; };
  for (const e of elems) for (let k = 1; k < e.nd.length; k++) { const a = fnd(e.nd[0][0]), b = fnd(e.nd[k][0]); if (a !== b) par[a] = b; }
  const hostC = new Set(elems.filter(e => e.role === 'host').map(e => fnd(e.nd[0][0])));
  const floating = [...new Set(elems.filter(e => !hostC.has(fnd(e.nd[0][0]))).map(e => e.part))].sort((a, b) => a - b);
  return { a1, a2, S: Math.abs(det), nodes, elems, pmass, springs, N, Henv, h, Lref, floating, merged: [...merged].sort((a, b) => a - b), hasOffset };
}

// ---------- complex Hermitian helpers (separate real / imaginary arrays) ----------
function hermChol(Ar, Ai, n) { // in place: lower L with real diagonal, A = L L^H
  for (let j = 0; j < n; j++) {
    const Rj = Ar[j], Ij = Ai[j]; let d = Rj[j]; for (let k = 0; k < j; k++) d -= Rj[k] * Rj[k] + Ij[k] * Ij[k];
    if (!(d > 0)) return false; const ljj = Math.sqrt(d); Rj[j] = ljj; Ij[j] = 0;
    for (let i = j + 1; i < n; i++) { const Ri = Ar[i], Ii = Ai[i]; let sr = Ri[j], si = Ii[j];
      for (let k = 0; k < j; k++) { const ar = Ri[k], ai = Ii[k], br = Rj[k], bi = -Ij[k]; sr -= ar * br - ai * bi; si -= ar * bi + ai * br; }
      Ri[j] = sr / ljj; Ii[j] = si / ljj; }
  }
  return true;
}
function hermSolve(Lr, Li, n, br, bi) { // L L^H x = b, in place
  for (let i = 0; i < n; i++) { const Ri = Lr[i], Ii = Li[i]; let sr = br[i], si = bi[i]; for (let k = 0; k < i; k++) { sr -= Ri[k] * br[k] - Ii[k] * bi[k]; si -= Ri[k] * bi[k] + Ii[k] * br[k]; } br[i] = sr / Ri[i]; bi[i] = si / Ri[i]; }
  for (let i = n - 1; i >= 0; i--) { let sr = br[i], si = bi[i]; for (let k = i + 1; k < n; k++) { const lr = Lr[k][i], li = -Li[k][i]; sr -= lr * br[k] - li * bi[k]; si -= lr * bi[k] + li * br[k]; } br[i] = sr / Lr[i][i]; bi[i] = si / Lr[i][i]; }
}
// Bloch matrices of a model
function blochMats(md, kx, ky) {
  const n = 6 * md.N + md.springs.length, Kr = zeros(n, n), Ki = zeros(n, n), Mr = zeros(n, n), Mi = zeros(n, n);
  for (const e of md.elems) {
    const nn = e.nd.length, ph = e.nd.map(([, s]) => { const R0 = s[0] * md.a1[0] + s[1] * md.a2[0], R1 = s[0] * md.a1[1] + s[1] * md.a2[1], a = kx * R0 + ky * R1; return [Math.cos(a), Math.sin(a)]; });
    for (let A = 0; A < nn; A++) for (let Bq = 0; Bq < nn; Bq++) {
      // u_B in cell s_B = u_B e^{i k·R_B}; entry (A,B) gets e^{-i k R_A} e^{i k R_B}
      const cr = ph[A][0] * ph[Bq][0] + ph[A][1] * ph[Bq][1], ci = ph[A][0] * ph[Bq][1] - ph[A][1] * ph[Bq][0];
      const oA = 6 * e.nd[A][0], oB = 6 * e.nd[Bq][0];
      for (let a = 0; a < 6; a++) { const KA = e.K[6 * A + a], MA = e.M[6 * A + a], kr = Kr[oA + a], ki = Ki[oA + a], mr = Mr[oA + a], mi = Mi[oA + a];
        for (let b = 0; b < 6; b++) { const k = KA[6 * Bq + b], m = MA[6 * Bq + b]; if (k) { kr[oB + b] += k * cr; ki[oB + b] += k * ci; } if (m) { mr[oB + b] += m * cr; mi[oB + b] += m * ci; } } }
    }
  }
  for (const pm of md.pmass) { for (let a = 0; a < 3; a++) Mr[6 * pm.node + a][6 * pm.node + a] += pm.m; for (let a = 3; a < 6; a++) Mr[6 * pm.node + a][6 * pm.node + a] += pm.J; }
  md.springs.forEach((sp, j) => { const J = 6 * md.N + j; Kr[J][J] += sp.k; Mr[J][J] += sp.m; for (let a = 0; a < 3; a++) { Mr[6 * sp.node + a][6 * sp.node + a] += sp.mr != null ? sp.mr : sp.m; const v = sp.m * sp.n[a]; Mr[6 * sp.node + a][J] += v; Mr[J][6 * sp.node + a] += v; } });
  return { Kr, Ki, Mr, Mi, n };
}
// lowest nb eigenpairs of K(k) u = λ M(k) u by shift-invert block Krylov (M-orthonormal basis, full re-orthogonalisation)
function lowestModes(Kr, Ki, Mr, Mi, n, nb, sigma, opt = {}) {
  const Ar = Kr.map((r, i) => { const o = Float64Array.from(r); for (let j = 0; j < n; j++) o[j] += sigma * Mr[i][j]; return o; }), Ai = Ki.map((r, i) => { const o = Float64Array.from(r); for (let j = 0; j < n; j++) o[j] += sigma * Mi[i][j]; return o; });
  if (!hermChol(Ar, Ai, n)) throw new Error('The shifted stiffness is not positive definite: part of the drawing may be a mechanism (free to move without strain).');
  const Mx = (xr, xi) => { const yr = new Float64Array(n), yi = new Float64Array(n); for (let i = 0; i < n; i++) { const mr = Mr[i], mi = Mi[i]; let sr = 0, si = 0; for (let j = 0; j < n; j++) { const a = mr[j], b = mi[j]; if (!a && !b) continue; sr += a * xr[j] - b * xi[j]; si += a * xi[j] + b * xr[j]; } yr[i] = sr; yi[i] = si; } return [yr, yi]; };
  return lowestModesOp(n, nb, sigma, (br, bi) => hermSolve(Ar, Ai, n, br, bi), Mx, opt);
}
// shift-invert block Krylov core: solve(b) overwrites b with (K + σM)⁻¹ b; Mx(x) returns M x
function lowestModesOp(n, nb, sigma, solve, Mx, opt = {}) {
  const ip = (ar, ai, br, bi) => { let sr = 0, si = 0; for (let i = 0; i < n; i++) { sr += ar[i] * br[i] + ai[i] * bi[i]; si += ar[i] * bi[i] - ai[i] * br[i]; } return [sr, si]; }; // a^H b
  const bs = Math.min(opt.block || 6, n), m = Math.min(n, opt.m || Math.max(2 * nb + 24, 36));
  const Q = [], MQ = []; // basis vectors and M times them
  let seed = 12345; const rnd = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648 - 0.5; };
  const addVecRel = (vr, vi) => {
    for (let pass = 0; pass < 2; pass++) for (let q = 0; q < Q.length; q++) { const [cr, ci] = ip(MQ[q][0], MQ[q][1], vr, vi); for (let i = 0; i < n; i++) { const qr = Q[q][0][i], qi = Q[q][1][i]; vr[i] -= cr * qr - ci * qi; vi[i] -= cr * qi + ci * qr; } }
    const [mr, mi] = Mx(vr, vi), [nn] = ip(vr, vi, mr, mi); if (!(nn > 0) || !isFinite(nn)) return false;
    const s = 1 / Math.sqrt(nn); for (let i = 0; i < n; i++) { vr[i] *= s; vi[i] *= s; mr[i] *= s; mi[i] *= s; }
    Q.push([vr, vi]); MQ.push([mr, mi]); return true; };
  let block = [];
  for (let b = 0; b < bs; b++) { const vr = new Float64Array(n), vi = new Float64Array(n); for (let i = 0; i < n; i++) { vr[i] = rnd(); vi[i] = rnd(); } if (addVecRel(vr, vi)) block.push(Q.length - 1); }
  const W = []; // Op applied to Q[q]: W[q] = A^{-1} M Q[q]
  while (Q.length < m && block.length) {
    const next = [];
    for (const q of block) { if (Q.length >= m) break; const wr = Float64Array.from(MQ[q][0]), wi = Float64Array.from(MQ[q][1]); solve(wr, wi); W[q] = [Float64Array.from(wr), Float64Array.from(wi)];
      const [m0r, m0i] = Mx(wr, wi), [n0] = ip(wr, wi, m0r, m0i);
      for (let pass = 0; pass < 2; pass++) for (let p = 0; p < Q.length; p++) { const [cr, ci] = ip(MQ[p][0], MQ[p][1], wr, wi); for (let i = 0; i < n; i++) { const qr = Q[p][0][i], qi = Q[p][1][i]; wr[i] -= cr * qr - ci * qi; wi[i] -= cr * qi + ci * qr; } }
      const [mr, mi] = Mx(wr, wi), [nn] = ip(wr, wi, mr, mi);
      if (!(nn > 1e-20 * n0)) continue;
      const s = 1 / Math.sqrt(nn); for (let i = 0; i < n; i++) { wr[i] *= s; wi[i] *= s; mr[i] *= s; mi[i] *= s; }
      Q.push([wr, wi]); MQ.push([mr, mi]); next.push(Q.length - 1); }
    block = next;
  }
  const mq = Q.length;
  for (let q = 0; q < mq; q++) if (!W[q]) { const wr = Float64Array.from(MQ[q][0]), wi = Float64Array.from(MQ[q][1]); solve(wr, wi); W[q] = [wr, wi]; }
  // H = Q^H M W  (Hermitian, eigenvalues θ = 1/(λ + σ))
  const Hr = zeros(mq, mq), Hi = zeros(mq, mq);
  for (let p = 0; p < mq; p++) for (let q = 0; q < mq; q++) { const [a, b] = ip(MQ[p][0], MQ[p][1], W[q][0], W[q][1]); Hr[p][q] = a; Hi[p][q] = b; }
  for (let p = 0; p < mq; p++) for (let q = 0; q < p; q++) { const a = 0.5 * (Hr[p][q] + Hr[q][p]), b = 0.5 * (Hi[p][q] - Hi[q][p]); Hr[p][q] = Hr[q][p] = a; Hi[p][q] = b; Hi[q][p] = -b; }
  const N2 = 2 * mq, C = zeros(N2, N2);
  for (let i = 0; i < mq; i++) for (let j = 0; j < mq; j++) { C[i][j] = C[i + mq][j + mq] = Hr[i][j]; C[i][j + mq] = -Hi[i][j]; C[i + mq][j] = Hi[i][j]; }
  const d = new Float64Array(N2), e = new Float64Array(N2);
  Eg.tred2(C, N2, d, e, true); Eg.tqli(d, e, N2, C);
  const ord = Array.from({ length: N2 }, (_, i) => i).sort((p, q) => opt.nearest ? Math.abs(d[q]) - Math.abs(d[p]) : d[q] - d[p]); // largest θ first (nearest the shift)
  const out = { lam: [], vec: [] };
  for (let k = 0; k < Math.min(nb, mq); k++) {
    const j = ord[2 * k], th = d[j]; out.lam.push(th > 0 || (opt.nearest && th) ? 1 / th - sigma : Infinity);
    if (opt.vectors) { const zr = new Float64Array(n), zi = new Float64Array(n); for (let q = 0; q < mq; q++) { const cr = C[q][j], ci = C[q + mq][j]; if (!cr && !ci) continue; for (let i = 0; i < n; i++) { zr[i] += cr * Q[q][0][i] - ci * Q[q][1][i]; zi[i] += cr * Q[q][1][i] + ci * Q[q][0][i]; } } out.vec.push([zr, zi]); }
  }
  return out;
}
// Bloch FE at one k: frequencies (Hz) and, per band, kinetic-energy shares {z, res}
function blochAt(md, kx, ky, nb, fref, wantVec) {
  const { Kr, Ki, Mr, Mi, n } = blochMats(md, kx, ky), sigma = (TAU * Math.max(1, fref || 100)) ** 2 * 0.05;
  const r = lowestModes(Kr, Ki, Mr, Mi, n, Math.min(nb, n), sigma, { vectors: true });
  if (!md.resDof) { md.resDof = new Uint8Array(n); for (const e of md.elems) if (e.role === 'res') for (const [id] of e.nd) for (let a = 0; a < 6; a++) md.resDof[6 * id + a] = 1; for (let j = 0; j < md.springs.length; j++) md.resDof[6 * md.N + j] = 1; }
  const f = [], z = [], res = [];
  r.lam.forEach((lam, b) => { f.push(Math.sqrt(Math.max(0, lam)) / TAU); const [vr, vi] = r.vec[b]; let tot = 0, ez = 0, er = 0;
    for (let i = 0; i < n; i++) { const w = Mr[i][i] * (vr[i] * vr[i] + vi[i] * vi[i]); tot += w; if (i < 6 * md.N && i % 6 === 2) ez += w; if (md.resDof[i]) er += w; }
    // in-plane share = translations x, y; z share relative to all translations
    let tr = 0; for (let i = 0; i < 6 * md.N; i++) if (i % 6 < 3) tr += Mr[i][i] * (vr[i] * vr[i] + vi[i] * vi[i]);
    z.push(tr > 0 ? ez / tr : 0); res.push(tot > 0 ? er / tot : 0); });
  return wantVec ? { f, z, res, vec: r.vec } : { f, z, res };
}

// ---------- banded solver for long strips: nodes ordered by depth, banded complex Cholesky, sparse M ----------
function stripAt(md, kx, ky, nb, fref, w) {
  const NN = md.N, ns = md.springs.length, n = 6 * NN + ns, order = Array.from({ length: NN }, (_, i) => i).sort((a, b) => md.strip.depth[a] - md.strip.depth[b]);
  const sprAt = Array.from({ length: NN }, () => []); md.springs.forEach((sp, j) => sprAt[sp.node].push(j));
  const pos = new Int32Array(n); let g = 0; // original dof -> banded index
  for (const id of order) { for (let a = 0; a < 6; a++) pos[6 * id + a] = g++; for (const j of sprAt[id]) pos[6 * NN + j] = g++; }
  const rows = Array.from({ length: n }, () => new Map()); // banded index -> Map(col -> [kr, ki, mr, mi])
  const add = (i, j, kr, ki, mr, mi) => { const r = rows[i]; let v = r.get(j); if (!v) { v = [0, 0, 0, 0]; r.set(j, v); } v[0] += kr; v[1] += ki; v[2] += mr; v[3] += mi; };
  for (const e of md.elems) {
    const nn = e.nd.length, ph = e.nd.map(([, s]) => { const R0 = s[0] * md.a1[0] + s[1] * md.a2[0], R1 = s[0] * md.a1[1] + s[1] * md.a2[1], a = kx * R0 + ky * R1; return [Math.cos(a), Math.sin(a)]; });
    for (let A = 0; A < nn; A++) for (let Bq = 0; Bq < nn; Bq++) {
      const cr = ph[A][0] * ph[Bq][0] + ph[A][1] * ph[Bq][1], ci = ph[A][0] * ph[Bq][1] - ph[A][1] * ph[Bq][0], oA = 6 * e.nd[A][0], oB = 6 * e.nd[Bq][0];
      for (let a = 0; a < 6; a++) { const KA = e.K[6 * A + a], MA = e.M[6 * A + a], I = pos[oA + a];
        for (let b = 0; b < 6; b++) { const k = KA[6 * Bq + b], m = MA[6 * Bq + b]; if (k || m) add(I, pos[oB + b], k * cr, k * ci, m * cr, m * ci); } }
    }
  }
  for (const pm of md.pmass) for (let a = 0; a < 6; a++) { const I = pos[6 * pm.node + a]; add(I, I, 0, 0, a < 3 ? pm.m : pm.J, 0); }
  md.springs.forEach((sp, j) => { const J = pos[6 * NN + j]; add(J, J, sp.k, 0, sp.m, 0); for (let a = 0; a < 3; a++) { const I = pos[6 * sp.node + a], v = sp.m * sp.n[a]; add(I, I, 0, 0, sp.mr != null ? sp.mr : sp.m, 0); if (v) { add(I, J, 0, 0, v, 0); add(J, I, 0, 0, v, 0); } } });
  const sigma = (TAU * Math.max(1, fref || 100)) ** 2 * 0.05;
  let bw = 0; for (let i = 0; i < n; i++) for (const j of rows[i].keys()) if (i - j > bw) bw = i - j;
  const W = bw + 1, Lr = new Float64Array(n * W), Li = new Float64Array(n * W); // L[i][j] at i*W + (j - i + bw)
  const cr = [], cc = [], vr = [], vi = []; // sparse M (CSR)
  for (let i = 0; i < n; i++) { const cols = [], re = [], im = [];
    for (const [j, v] of rows[i]) { if (j <= i) { Lr[i * W + j - i + bw] += v[0] + sigma * v[2]; Li[i * W + j - i + bw] += v[1] + sigma * v[3]; } if (v[2] || v[3]) { cols.push(j); re.push(v[2]); im.push(v[3]); } }
    cc.push(Int32Array.from(cols)); vr.push(Float64Array.from(re)); vi.push(Float64Array.from(im)); rows[i] = null; }
  // banded Hermitian Cholesky
  for (let j = 0; j < n; j++) {
    const oj = j * W - j + bw, k0 = Math.max(0, j - bw); let d = Lr[oj + j];
    for (let k = k0; k < j; k++) d -= Lr[oj + k] * Lr[oj + k] + Li[oj + k] * Li[oj + k];
    if (!(d > 0)) throw new Error('The shifted stiffness of the strip is not positive definite: part of the drawing may be a mechanism.');
    const l = Math.sqrt(d); Lr[oj + j] = l; Li[oj + j] = 0;
    const iEnd = Math.min(n - 1, j + bw);
    for (let i = j + 1; i <= iEnd; i++) { const oi = i * W - i + bw; let sr = Lr[oi + j], si = Li[oi + j];
      for (let k = Math.max(k0, i - bw); k < j; k++) { const ar = Lr[oi + k], ai = Li[oi + k], br = Lr[oj + k], bi = -Li[oj + k]; sr -= ar * br - ai * bi; si -= ar * bi + ai * br; }
      Lr[oi + j] = sr / l; Li[oi + j] = si / l; }
  }
  const solve = (br, bi) => {
    for (let i = 0; i < n; i++) { const oi = i * W - i + bw; let sr = br[i], si = bi[i]; for (let k = Math.max(0, i - bw); k < i; k++) { const lr = Lr[oi + k], li = Li[oi + k]; sr -= lr * br[k] - li * bi[k]; si -= lr * bi[k] + li * br[k]; } br[i] = sr / Lr[oi + i]; bi[i] = si / Lr[oi + i]; }
    for (let i = n - 1; i >= 0; i--) { let sr = br[i], si = bi[i]; const kEnd = Math.min(n - 1, i + bw); for (let k = i + 1; k <= kEnd; k++) { const o = k * W - k + bw + i, lr = Lr[o], li = -Li[o]; sr -= lr * br[k] - li * bi[k]; si -= lr * bi[k] + li * br[k]; } const d = Lr[i * W + bw]; br[i] = sr / d; bi[i] = si / d; }
  };
  const Mx = (xr, xi) => { const yr = new Float64Array(n), yi = new Float64Array(n); for (let i = 0; i < n; i++) { const c = cc[i], a = vr[i], b = vi[i]; let sr = 0, si = 0; for (let q = 0; q < c.length; q++) { const j = c[q]; sr += a[q] * xr[j] - b[q] * xi[j]; si += a[q] * xi[j] + b[q] * xr[j]; } yr[i] = sr; yi[i] = si; } return [yr, yi]; };
  const r = lowestModesOp(n, Math.min(nb, n), sigma, solve, Mx, { vectors: true });
  // back to the standard dof numbering; energy shares
  const diag = new Float64Array(n); for (let i = 0; i < n; i++) { const c = cc[i]; for (let q = 0; q < c.length; q++) if (c[q] === i) diag[i] = vr[i][q]; }
  const dep = md.strip.depth, lo = md.strip.dmin + w, hi = md.strip.dmax - w, dofDepth = i => i < 6 * NN ? dep[Math.floor(i / 6)] : dep[md.springs[i - 6 * NN].node];
  const out = { f: [], z: [], bot: [], top: [], vec: [] };
  r.lam.forEach((lam, b) => { const [ur, ui] = r.vec[b], xr = new Float64Array(n), xi = new Float64Array(n); for (let i = 0; i < n; i++) { xr[i] = ur[pos[i]]; xi[i] = ui[pos[i]]; }
    let tot = 0, eb = 0, et = 0, ez = 0, tr = 0;
    for (let i = 0; i < n; i++) { const wgt = diag[pos[i]] * (xr[i] * xr[i] + xi[i] * xi[i]); tot += wgt; const d = dofDepth(i); if (d < lo) eb += wgt; if (d > hi) et += wgt; if (i < 6 * NN && i % 6 < 3) { tr += wgt; if (i % 6 === 2) ez += wgt; } }
    out.f.push(Math.sqrt(Math.max(0, lam)) / TAU); out.z.push(tr > 0 ? ez / tr : 0); out.bot.push(tot > 0 ? eb / tot : 0); out.top.push(tot > 0 ? et / tot : 0); out.vec.push([xr, xi]); });
  out.n = n; out.bw = bw; return out;
}

// ---------- finite sample: M × N copies of the cell, every edge free; modes near a chosen frequency ----------
function buildFinite(dr, mat, M, N, opt = {}) {
  // shape 'para': M × N cells along a1, a2; 'rect': N columns along a1, each M cells long, columns staggered so the sample is rectangular
  const prim = expandParts(dr.parts || [], dr), parts = [], a1 = dr.a1, a2 = dr.a2, L1 = Math.hypot(...a1), L2 = Math.hypot(...a2);
  const d12 = Math.hypot(a2[0] - a1[0], a2[1] - a1[1]), rect = opt.shape === 'rect' && Math.abs(d12 - L1) < 0.02 * L1 && Math.abs(L2 - L1) < 0.02 * L1;
  const ij = []; if (rect) { for (let p = 0; p < N; p++) for (let r = 0; r < M; r++) ij.push([r - Math.floor(p / 2), p]); } else for (let i = 0; i < M; i++) for (let j = 0; j < N; j++) ij.push([i, j]);
  const o1 = 2 + Math.ceil(N / 2), o2 = 2, cells = [];
  for (const [i, j] of ij) { const sx = (i + o1) * a1[0] + (j + o2) * a2[0], sy = (i + o1) * a1[1] + (j + o2) * a2[1]; cells.push({ c: [sx, sy], ij: [i, j] });
    for (const p of prim) { const q = JSON.parse(JSON.stringify(p)); q.src = p.src; transformPart(q, ([x, y]) => [x + sx, y + sy]); parts.push(q); } }
  const sM = M + N + 8, sN = N + 6, md = build({ a1: [a1[0] * sM, a1[1] * sM], a2: [a2[0] * sN, a2[1] * sN], mesh: dr.mesh, parts }, mat, opt);
  for (const e of md.elems) for (const [, s] of e.nd) if (s[0] || s[1]) throw new Error('Internal: the finite sample wrapped around.');
  // neighbour directions of the lattice (the shortest of ±a1, ±a2, ±(a1−a2), ±(a1+a2)) and the cells that miss some of them
  const cand = [a1, a2, [a1[0] - a2[0], a1[1] - a2[1]], [a1[0] + a2[0], a1[1] + a2[1]]], lmin = Math.min(...cand.map(v => Math.hypot(...v)));
  const nb = []; for (const v of cand) if (Math.hypot(...v) < 1.02 * lmin) nb.push(v, [-v[0], -v[1]]);
  const tol = 0.05 * lmin, at = new Map(cells.map((c, k) => [Math.round(c.c[0] / tol) + ',' + Math.round(c.c[1] / tol), k]));
  const find = (x, y) => { for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) { const k = at.get((Math.round(x / tol) + dx) + ',' + (Math.round(y / tol) + dy)); if (k != null) return k; } return -1; };
  for (const c of cells) { c.miss = nb.filter(v => find(c.c[0] + v[0], c.c[1] + v[1]) < 0); c.corner = nb.length === 4 ? c.miss.length >= 2 && !(c.miss.length === 2 && Math.abs(c.miss[0][0] + c.miss[1][0]) + Math.abs(c.miss[0][1] + c.miss[1][1]) < 1e-6 * lmin) : c.miss.length > nb.length / 2; c.edge = c.miss.length > 0; }
  const X = md.nodes.map(x => [x[0] * 1000, x[1] * 1000]), cellOf = new Int32Array(md.N), zone = new Uint8Array(md.N);
  X.forEach((x, i) => { let b = 0, d = Infinity; cells.forEach((c, k) => { const dd = (x[0] - c.c[0]) ** 2 + (x[1] - c.c[1]) ** 2; if (dd < d) { d = dd; b = k; } }); cellOf[i] = b;
    const c = cells[b]; zone[i] = c.corner ? 2 : c.miss.some(v => (x[0] - c.c[0]) * v[0] + (x[1] - c.c[1]) * v[1] > 0) ? 1 : 0; });
  // order the nodes along the longer extent of the sample, so the band covers the shorter one
  const u = [a1[0] / L1, a1[1] / L1], pr = x => x[0] * u[0] + x[1] * u[1], qr = x => -x[0] * u[1] + x[1] * u[0], ext = f => Math.max(...cells.map(c => f(c.c))) - Math.min(...cells.map(c => f(c.c)));
  const along = ext(pr) >= ext(qr) ? pr : qr, key = Float64Array.from(X, x => along(x));
  md.Lref = Math.max(L1, L2) / 1000;
  md.fin = { M, N, shape: rect ? 'rect' : 'para', cells, cellOf, zone, key, u };
  return md;
}
// real banded LDLᵀ of K − s M (indefinite), nodes ordered across the shorter side; Krylov for the modes nearest s
function finiteModes(md, fc, nm, status) {
  const NN = md.N, ns = md.springs.length, n = 6 * NN + ns, key = md.fin.key;
  const order = Array.from({ length: NN }, (_, i) => i).sort((a, b) => key[a] - key[b]);
  const sprAt = Array.from({ length: NN }, () => []); md.springs.forEach((sp, j) => sprAt[sp.node].push(j));
  const pos = new Int32Array(n); let g = 0; for (const id of order) { for (let a = 0; a < 6; a++) pos[6 * id + a] = g++; for (const j of sprAt[id]) pos[6 * NN + j] = g++; }
  const rows = Array.from({ length: n }, () => new Map()), add = (i, j, k, m) => { const r = rows[i]; let v = r.get(j); if (!v) { v = [0, 0]; r.set(j, v); } v[0] += k; v[1] += m; };
  for (const e of md.elems) { const nn = e.nd.length; for (let A = 0; A < nn; A++) for (let Bq = 0; Bq < nn; Bq++) { const oA = 6 * e.nd[A][0], oB = 6 * e.nd[Bq][0];
    for (let a = 0; a < 6; a++) { const KA = e.K[6 * A + a], MA = e.M[6 * A + a], I = pos[oA + a]; for (let b = 0; b < 6; b++) { const k = KA[6 * Bq + b], m = MA[6 * Bq + b]; if (k || m) add(I, pos[oB + b], k, m); } } } }
  for (const pm of md.pmass) for (let a = 0; a < 6; a++) { const I = pos[6 * pm.node + a]; add(I, I, 0, a < 3 ? pm.m : pm.J); }
  md.springs.forEach((sp, j) => { const J = pos[6 * NN + j]; add(J, J, sp.k, sp.m); for (let a = 0; a < 3; a++) { const I = pos[6 * sp.node + a], v = sp.m * sp.n[a]; add(I, I, 0, sp.mr != null ? sp.mr : sp.m); if (v) { add(I, J, 0, v); add(J, I, 0, v); } } });
  let bw = 0; for (let i = 0; i < n; i++) for (const j of rows[i].keys()) if (i - j > bw) bw = i - j;
  if (n * (bw + 1) > 4e7) throw new Error(`The sample is too large for the browser (${n} degrees of freedom, band ${bw}). Use fewer cells or a larger element size.`);
  const W = bw + 1, Mcol = [], Mval = [], diag = new Float64Array(n);
  for (let i = 0; i < n; i++) { const c = [], v = []; for (const [j, q] of rows[i]) if (q[1]) { c.push(j); v.push(q[1]); if (j === i) diag[i] = q[1]; } Mcol.push(Int32Array.from(c)); Mval.push(Float64Array.from(v)); }
  let s = (TAU * fc) ** 2, L = null, D = null;
  for (let attempt = 0; attempt < 4 && !L; attempt++) { // shifted LDLᵀ; nudge the shift if it lands on an eigenvalue
    L = new Float64Array(n * W); D = new Float64Array(n);
    for (let i = 0; i < n; i++) for (const [j, q] of rows[i]) if (j <= i) L[i * W + j - i + bw] += q[0] - s * q[1];
    let ok = true;
    for (let j = 0; j < n && ok; j++) { const oj = j * W - j + bw, k0 = Math.max(0, j - bw); let d = L[oj + j]; for (let k = k0; k < j; k++) d -= L[oj + k] * L[oj + k] * D[k];
      if (!(Math.abs(d) > 1e-13 * Math.abs(L[oj + j] || 1))) { ok = false; break; } D[j] = d; L[oj + j] = 1;
      const iEnd = Math.min(n - 1, j + bw); for (let i = j + 1; i <= iEnd; i++) { const oi = i * W - i + bw; let t = L[oi + j]; for (let k = Math.max(k0, i - bw); k < j; k++) t -= L[oi + k] * L[oj + k] * D[k]; L[oi + j] = t / d; } }
    if (!ok) { L = null; s *= 1 + 1e-3 * (attempt + 1); } }
  if (!L) throw new Error('Could not factorise the sample at this frequency; try a slightly different centre frequency.');
  for (let i = 0; i < n; i++) rows[i] = null;
  const solveR = b => { for (let i = 0; i < n; i++) { const oi = i * W - i + bw; let t = b[i]; for (let k = Math.max(0, i - bw); k < i; k++) t -= L[oi + k] * b[k]; b[i] = t; }
    for (let i = 0; i < n; i++) b[i] /= D[i];
    for (let i = n - 1; i >= 0; i--) { let t = b[i]; const kEnd = Math.min(n - 1, i + bw); for (let k = i + 1; k <= kEnd; k++) t -= L[k * W - k + bw + i] * b[k]; b[i] = t; } };
  const solve = (br, bi) => { solveR(br); solveR(bi); };
  const Mx = (xr, xi) => { const yr = new Float64Array(n), yi = new Float64Array(n); for (let i = 0; i < n; i++) { const c = Mcol[i], v = Mval[i]; let a = 0, b = 0; for (let q = 0; q < c.length; q++) { a += v[q] * xr[c[q]]; b += v[q] * xi[c[q]]; } yr[i] = a; yi[i] = b; } return [yr, yi]; };
  const r = lowestModesOp(n, Math.min(nm, n), -s, solve, Mx, { vectors: true, nearest: true, block: 4 });
  // corner and edge zones: within one cell of the sample's boundary
  // zones: corner cells (2); the outer half of every other boundary cell (1); the rest (0)
  const zoneOf = id => md.fin.zone[id];
  const zone = new Uint8Array(n); for (let i = 0; i < n; i++) zone[i] = zoneOf(i < 6 * NN ? Math.floor(i / 6) : md.springs[i - 6 * NN].node);
  const out = [];
  r.lam.forEach((lam, b) => { if (!isFinite(lam) || lam < (0.02 * TAU * fc) ** 2) return; const [ur, ui] = r.vec[b], xr = new Float64Array(n), xi = new Float64Array(n); for (let i = 0; i < n; i++) { xr[i] = ur[pos[i]]; xi[i] = ui[pos[i]]; }
    let tot = 0, ez = [0, 0, 0], zt = 0, tr = 0;
    for (let i = 0; i < n; i++) { const w = diag[pos[i]] * (xr[i] * xr[i] + xi[i] * xi[i]); tot += w; ez[zone[i]] += w; if (i < 6 * NN && i % 6 < 3) { tr += w; if (i % 6 === 2) zt += w; } }
    out.push({ f: Math.sqrt(lam) / TAU, bulk: ez[0] / tot, edge: ez[1] / tot, corner: ez[2] / tot, z: tr > 0 ? zt / tr : 0, vec: [xr, xi] }); });
  out.sort((a, b) => a.f - b.f);
  // share of the sample area in each zone, for reference
  const area = [0, 0, 0]; for (let i = 0; i < NN; i++) area[zoneOf(i)]++;
  const mdiag = new Float64Array(n); for (let i = 0; i < n; i++) mdiag[i] = diag[pos[i]];
  return { modes: out, n, bw, zoneShare: area.map(a => a / NN), mdiag };
}

// ---------- homogenisation and resonator reduction of a drawn model ----------
function homogenize(md) {
  const host = md.elems.filter(e => e.role === 'host'), used = new Map();
  for (const e of host) for (const [id] of e.nd) if (!used.has(id)) used.set(id, used.size);
  const nm = used.size, nd = 6 * nm, S = md.S;
  const mHost = host.reduce((s, e) => s + e.mass, 0) + md.pmass.filter(p => p.role === 'host').reduce((s, p) => s + p.m, 0);
  const zref = host.reduce((s, e) => s + e.mass * e.X.reduce((a, x) => a + x[2], 0) / e.X.length, 0) / host.reduce((s, e) => s + e.mass, 0);
  let x0 = 0, y0 = 0; for (const id of used.keys()) { x0 += md.nodes[id][0]; y0 += md.nodes[id][1]; } x0 /= nm; y0 /= nm;
  const Fat = p => { const x = p[0] - x0, y = p[1] - y0, z = p[2] - zref; return [[x, 0, y / 2, -z * x, 0, -z * y / 2], [0, y, x / 2, 0, -z * y, -z * x / 2], [0, 0, 0, x * x / 2, y * y / 2, x * y / 2], [0, 0, 0, 0, y, x / 2], [0, 0, 0, -x, 0, -y / 2], [0, 0, 0, 0, 0, 0]]; };
  const Kpp = zeros(nd, nd), Kpf = zeros(nd, 6), Kff = zeros(6, 6);
  for (const e of host) {
    const ne = 6 * e.nd.length, F = e.X.flatMap(Fat), idx = []; for (const [id] of e.nd) for (let a = 0; a < 6; a++) idx.push(6 * used.get(id) + a);
    const KF = zeros(ne, 6); for (let a = 0; a < ne; a++) { const Ka = e.K[a]; for (let c = 0; c < 6; c++) { let v = 0; for (let b = 0; b < ne; b++) v += Ka[b] * F[b][c]; KF[a][c] = v; } }
    for (let a = 0; a < ne; a++) { const ra = Kpp[idx[a]], Ka = e.K[a]; for (let b = 0; b < ne; b++) ra[idx[b]] += Ka[b]; for (let c = 0; c < 6; c++) Kpf[idx[a]][c] += KF[a][c]; }
    for (let r = 0; r < 6; r++) for (let c = 0; c < 6; c++) { let v = 0; for (let a = 0; a < ne; a++) v += F[a][r] * KF[a][c]; Kff[r][c] += v; }
  }
  const n = nd - 3, Kr = zeros(n, n); for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) Kr[a][b] = 0.5 * (Kpp[a + 3][b + 3] + Kpp[b + 3][a + 3]);
  // Cholesky with a check
  let dmax = 0; for (let i = 0; i < n; i++) dmax = Math.max(dmax, Kr[i][i]);
  for (let j = 0; j < n; j++) { const Rj = Kr[j]; let s = Rj[j]; for (let k = 0; k < j; k++) s -= Rj[k] * Rj[k]; if (!(s > 1e-14 * dmax)) throw new Error('The host frame has a mechanism or is disconnected (a part not joined to the rest, or joined only at a point that can rotate freely).'); const dj = Math.sqrt(s); Rj[j] = dj; for (let i = j + 1; i < n; i++) { const Ri = Kr[i]; let t = Ri[j]; for (let k = 0; k < j; k++) t -= Ri[k] * Rj[k]; Ri[j] = t / dj; } }
  const X = []; for (let c = 0; c < 6; c++) { const b = new Float64Array(n); for (let a = 0; a < n; a++) b[a] = Kpf[a + 3][c];
    for (let i = 0; i < n; i++) { let t = b[i]; for (let k = 0; k < i; k++) t -= Kr[i][k] * b[k]; b[i] = t / Kr[i][i]; } for (let i = n - 1; i >= 0; i--) { let t = b[i]; for (let k = i + 1; k < n; k++) t -= Kr[k][i] * b[k]; b[i] = t / Kr[i][i]; } X.push(b); }
  const Q = zeros(6, 6); for (let r = 0; r < 6; r++) for (let c = 0; c < 6; c++) { let v = Kff[r][c]; for (let a = 0; a < n; a++) v -= Kpf[a + 3][r] * X[c][a]; Q[r][c] = v / S; }
  for (let r = 0; r < 6; r++) for (let c = 0; c < r; c++) { const v = 0.5 * (Q[r][c] + Q[c][r]); Q[r][c] = Q[c][r] = v; }
  const A = [0, 1, 2].map(r => [0, 1, 2].map(c => Q[r][c])), B = [0, 1, 2].map(r => [0, 1, 2].map(c => Q[r][c + 3])), D = [0, 1, 2].map(r => [0, 1, 2].map(c => Q[r + 3][c + 3]));
  const inv3 = M => { const [a, b, c] = M[0], [d, e, f] = M[1], [g, h, i] = M[2], dt = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g); return [[(e * i - f * h) / dt, (c * h - b * i) / dt, (b * f - c * e) / dt], [(f * g - d * i) / dt, (a * i - c * g) / dt, (c * d - a * f) / dt], [(d * h - e * g) / dt, (b * g - a * h) / dt, (a * e - b * d) / dt]]; };
  const Ai = inv3(A), AiB = [0, 1, 2].map(r => [0, 1, 2].map(c => Ai[r][0] * B[0][c] + Ai[r][1] * B[1][c] + Ai[r][2] * B[2][c]));
  const Dfree = [0, 1, 2].map(r => [0, 1, 2].map(c => D[r][c] - (B[0][r] * AiB[0][c] + B[1][r] * AiB[1][c] + B[2][r] * AiB[2][c])));
  return { A, B, D, Dfree, Q, mCell: mHost, mArea: mHost / S, S, zref, Henv: md.Henv, rhoStar: mHost / (S * md.Henv), cond: 0, nElem: host.length, nNodes: nm };
}
// resonator components (connected 'res' elements and masses), clamped where they touch the host; modal reduction
function resonators(md, nModes, baseR = 0) {
  const hostNode = new Set(); for (const e of md.elems) if (e.role === 'host') for (const [id] of e.nd) hostNode.add(id);
  const resEl = md.elems.filter(e => e.role === 'res'), adj = new Map();
  resEl.forEach((e, k) => { for (const [id] of e.nd) if (!hostNode.has(id)) { if (!adj.has(id)) adj.set(id, []); adj.get(id).push(k); } });
  const seen = new Set(), out = [];
  for (const start of adj.keys()) {
    if (seen.has(start)) continue;
    const off = new Map([[start, [0, 0]]]), stack = [start], els = new Set(); seen.add(start);
    while (stack.length) { const i = stack.pop(), o = off.get(i); for (const k of adj.get(i)) { els.add(k); const e = resEl[k], me = e.nd.find(r => r[0] === i)[1];
      for (const [id, s] of e.nd) { const so = [o[0] + s[0] - me[0], o[1] + s[1] - me[1]]; if (!off.has(id)) off.set(id, so); if (!hostNode.has(id) && !seen.has(id)) { seen.add(id); stack.push(id); } } } }
    // local mesh
    const lid = new Map(), pos = [], root = [], rootIds = [];
    const L = (id, s) => { const k = `${id}|${s[0]}|${s[1]}`; if (!lid.has(k)) { lid.set(k, pos.length); pos.push([md.nodes[id][0] + s[0] * md.a1[0] + s[1] * md.a2[0], md.nodes[id][1] + s[0] * md.a1[1] + s[1] * md.a2[1], md.nodes[id][2]]); if (hostNode.has(id)) { root.push(pos.length - 1); rootIds.push(id); } } return lid.get(k); };
    const elsL = [];
    for (const k of els) { const e = resEl[k], me = e.nd.find(([id]) => !hostNode.has(id)) || e.nd[0]; const base = off.get(me[0]), b0 = [base[0] - me[1][0], base[1] - me[1][1]]; elsL.push({ e, ids: e.nd.map(([id, s]) => L(id, [b0[0] + s[0], b0[1] + s[1]])) }); }
    if (!root.length) throw new Error('A resonator part is not connected to the host: it must share a point with a host beam, wall or plate.');
    let clampIdx = root.slice(), hostIn = 0;
    // the host compliance around the joint matters for a resonator hung from the middle of a member; one joined at frame
    // junctions or along many points is clamped rigidly (Bloch FE comparisons of both kinds of cell)
    const hdeg = id => md.elems.reduce((c, e) => c + (e.role === 'host' && e.nd.some(r => r[0] === id) ? 1 : 0), 0);
    const rigid = root.length > 4 || rootIds.every(id => hdeg(id) >= 3);
    if (baseR > 0 && !rigid) { // carry the host around the attachment (within baseR) into the resonator model, clamp where it meets the rest
      const rootPos = root.map(r => pos[r]), inc = [], near = X => Math.min(...rootPos.map(q => nrm(sub(X, q))));
      md.elems.forEach((e, ei) => { if (e.role !== 'host') return; let best = null;
        for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) { const sh = [i * md.a1[0] + j * md.a2[0], i * md.a1[1] + j * md.a2[1], 0], Xs = e.X.map(x => [x[0] + sh[0], x[1] + sh[1], x[2]]), dmax = Math.max(...Xs.map(near)); if (!best || dmax < best.d) best = { d: dmax, Xs, i, j }; }
        if (best.d <= baseR * (1 + 1e-9)) inc.push({ e, ei, Xs: best.Xs, i: best.i, j: best.j }); });
      const incSet = new Set(inc.map(q => q.ei)), useOut = new Set(); md.elems.forEach((e, ei) => { if (e.role === 'host' && !incSet.has(ei)) for (const [id] of e.nd) useOut.add(id); });
      const keyP = X => X.map(v => Math.round(v * 1e9)).join(','), pidx = new Map(pos.map((p, i) => [keyP(p), i]));
      const getP = (X, id) => { const k = keyP(X); if (!pidx.has(k)) { pidx.set(k, pos.length); pos.push(X); } const i = pidx.get(k); if (useOut.has(id)) boundary.add(i); return i; };
      const boundary = new Set();
      for (const q of inc) { elsL.push({ e: q.e, ids: q.e.nd.map(([id], a) => getP(q.Xs[a], id)) }); hostIn += q.e.mass; }
      if (inc.length) clampIdx = [...boundary];
      if (!clampIdx.length) throw new Error('The resonator base region covers the whole host: reduce the base radius.');
    }
    const masses = md.pmass.filter(p => p.role === 'res' && off.has(p.node) && !hostNode.has(p.node)).map(p => ({ i: L(p.node, off.get(p.node)), m: p.m, J: p.J }));
    const N = pos.length, nd = 6 * N, K = zeros(nd, nd), M = zeros(nd, nd);
    for (const { e, ids } of elsL) for (let a = 0; a < ids.length * 6; a++) for (let b = 0; b < ids.length * 6; b++) { const I = 6 * ids[Math.floor(a / 6)] + a % 6, J = 6 * ids[Math.floor(b / 6)] + b % 6; K[I][J] += e.K[a][b]; M[I][J] += e.M[a][b]; }
    for (const pm of masses) { for (let a = 0; a < 3; a++) M[6 * pm.i + a][6 * pm.i + a] += pm.m; for (let a = 3; a < 6; a++) M[6 * pm.i + a][6 * pm.i + a] += pm.J; }
    const cl = new Set(); for (const r of clampIdx) for (let a = 0; a < 6; a++) cl.add(6 * r + a);
    const free = []; for (let i = 0; i < nd; i++) if (!cl.has(i)) free.push(i);
    const nf = free.length, Kf = zeros(nf, nf), Mf = zeros(nf, nf); for (let a = 0; a < nf; a++) for (let b = 0; b < nf; b++) { Kf[a][b] = K[free[a]][free[b]]; Mf[a][b] = M[free[a]][free[b]]; }
    { let mx = 0; for (let a = 0; a < nf; a++) mx = Math.max(mx, Mf[a][a]); for (let a = 0; a < nf; a++) if (!(Mf[a][a] > 1e-7 * mx)) Mf[a][a] = 1e-7 * mx; }
    const MR = [0, 1, 2].map(c => { const v = new Float64Array(nd); for (let i = 0; i < nd; i++) { let s = 0; for (let q = 0; q < N; q++) s += M[i][6 * q + c]; v[i] = s; } return v; });
    let mTot = 0; for (let q = 0; q < N; q++) mTot += MR[2][6 * q + 2];
    const eg = Hm.symGenEig(Kf, Mf, nf), modes = [];
    for (let j = 0; j < Math.min(nModes, nf); j++) { const phi = eg.vec[j], G = [0, 1, 2].map(c => { let s = 0; for (let a = 0; a < nf; a++) s += phi[a] * MR[c][free[a]]; return s; }), m = G[0] ** 2 + G[1] ** 2 + G[2] ** 2, w2 = Math.max(0, eg.lam[j]);
      modes.push({ f: Math.sqrt(w2) / TAU, m, k: w2 * m, n: m > 0 ? G.map(v => v / Math.sqrt(m)) : [0, 0, 1], share: m / mTot }); }
    out.push({ att: clampIdx.map(r => [pos[r][0], pos[r][1]]), pos: pos[root[0]], dir: [0, 0, 1], mTot, modes, hostIn, label: `Resonator ${out.length + 1}`, kind: 'drawn', parts: [...new Set([...els].map(k => resEl[k].part))] });
  }
  for (const sp of md.springs) { const p = md.nodes[sp.node]; out.push({ att: [[p[0], p[1]]], pos: p, dir: sp.n, mTot: sp.m, modes: [{ f: sp.f, m: sp.m, k: sp.k, n: sp.n, share: 1 }], label: 'Mass on spring', kind: 'spring' }); }
  return out;
}

// ---------- complex band of a drawn / 3D model: wave finite elements along a1 (phase along a2 = 0) ----------
// D(λ) = A0 + λ A+ + λ⁻¹ A− with λ = e^{iμ}, μ = k·a1; interior DOFs (not touched by cross-boundary terms) are condensed
// at each frequency, then the quadratic in λ is solved through λ = (1 + s)/(1 − s) as in the 2D workbench.
function wfeModel(md, keep) {
  const n = 6 * md.N + md.springs.length, K0 = zeros(n, n), M0 = zeros(n, n), Kp = new Map(), Mp = new Map(), C = new Set();
  const addP = (map, i, j, v) => { const k = i * n + j; map.set(k, (map.get(k) || 0) + v); };
  for (const e of md.elems) { const nn = e.nd.length;
    for (let A = 0; A < nn; A++) for (let Bq = 0; Bq < nn; Bq++) { const dd = e.nd[Bq][1][0] - e.nd[A][1][0]; if (Math.abs(dd) > 1) throw new Error('An element spans more than one cell along a₁: refine the mesh.');
      const oA = 6 * e.nd[A][0], oB = 6 * e.nd[Bq][0];
      for (let a = 0; a < 6; a++) for (let b = 0; b < 6; b++) { const k = e.K[6 * A + a][6 * Bq + b], m = e.M[6 * A + a][6 * Bq + b]; if (!k && !m) continue;
        if (dd === 0) { K0[oA + a][oB + b] += k; M0[oA + a][oB + b] += m; } else if (dd === 1) { addP(Kp, oA + a, oB + b, k); addP(Mp, oA + a, oB + b, m); C.add(oA + a); C.add(oB + b); } } } }
  for (const pm of md.pmass) { for (let a = 0; a < 3; a++) M0[6 * pm.node + a][6 * pm.node + a] += pm.m; for (let a = 3; a < 6; a++) M0[6 * pm.node + a][6 * pm.node + a] += pm.J; }
  md.springs.forEach((sp, j) => { const J = 6 * md.N + j; K0[J][J] += sp.k; M0[J][J] += sp.m; for (let a = 0; a < 3; a++) { M0[6 * sp.node + a][6 * sp.node + a] += sp.mr != null ? sp.mr : sp.m; const v = sp.m * sp.n[a]; M0[6 * sp.node + a][J] += v; M0[J][6 * sp.node + a] += v; } });
  // DOF scaling: rotations x Lref
  const Lr = md.Lref || 0.01, Sv = new Float64Array(n).fill(1); for (let i = 0; i < 6 * md.N; i++) if (i % 6 >= 3) Sv[i] = Lr;
  const ok = i => !keep || keep(i), cl = [...C].filter(ok).sort((a, b) => a - b), cset = new Set(cl), il = []; for (let i = 0; i < n; i++) if (!cset.has(i) && ok(i)) il.push(i);
  return { n, K0, M0, Kp, Mp, cl, il, Sv };
}
function wfeAt(W, f) {
  const w2 = (TAU * f) ** 2, { n, K0, M0, cl, il, Sv } = W, nc = cl.length, ni = il.length;
  if (!nc) return [];
  const D = (i, j) => (K0[i][j] - w2 * M0[i][j]) / (Sv[i] * Sv[j]);
  const DII = zeros(ni, ni), X = zeros(ni, nc);
  for (let a = 0; a < ni; a++) { for (let b = 0; b < ni; b++) DII[a][b] = D(il[a], il[b]); for (let b = 0; b < nc; b++) X[a][b] = D(il[a], cl[b]); }
  if (ni && !Eg.luSolve(DII, X, ni, nc)) return null;
  const R0 = zeros(nc, nc), Rp = zeros(nc, nc), Rm = zeros(nc, nc), pos = new Map(cl.map((c, i) => [c, i]));
  for (let a = 0; a < nc; a++) for (let b = 0; b < nc; b++) { let v = D(cl[a], cl[b]); for (let k = 0; k < ni; k++) v -= D(cl[a], il[k]) * X[k][b]; R0[a][b] = v; }
  for (const [key, k] of W.Kp) { const i = Math.floor(key / n), j = key % n; if (!pos.has(i) || !pos.has(j)) continue; const m = W.Mp.get(key) || 0, v = (k - w2 * m) / (Sv[i] * Sv[j]); Rp[pos.get(i)][pos.get(j)] += v; Rm[pos.get(j)][pos.get(i)] += v; }
  // (Rm + λ R0 + λ² Rp) u = 0, Möbius λ = (1 + s)/(1 − s)
  const Q2 = zeros(nc, nc), Q1 = zeros(nc, nc), Q0 = zeros(nc, nc);
  for (let i = 0; i < nc; i++) for (let j = 0; j < nc; j++) { const A2 = Rp[i][j], A1 = R0[i][j], A0 = Rm[i][j]; Q2[i][j] = A2 - A1 + A0; Q1[i][j] = 2 * (A2 - A0); Q0[i][j] = A2 + A1 + A0; }
  const Rh = zeros(nc, 2 * nc); for (let i = 0; i < nc; i++) for (let j = 0; j < nc; j++) { Rh[i][j] = Q0[i][j]; Rh[i][nc + j] = Q1[i][j]; }
  if (!Eg.luSolve(Q2, Rh, nc, 2 * nc)) return null;
  const T = zeros(2 * nc, 2 * nc); for (let i = 0; i < nc; i++) { T[i][nc + i] = 1; for (let j = 0; j < nc; j++) { T[nc + i][j] = -Rh[i][j]; T[nc + i][nc + j] = -Rh[i][nc + j]; } }
  for (let i = 0; i < 2 * nc; i++) for (let j = 0; j < 2 * nc; j++) if (!isFinite(T[i][j])) return null;
  const ev = Eg.eigNonsym(T, 2 * nc); if (!ev) return null;
  return ev.map(([sr, si]) => { const num = Math.hypot(1 + sr, si), den = Math.hypot(1 - sr, si); if (den < 1e-14 || num < 1e-14) return { att: 99, ph: 0 };
    let ph = Math.atan2(si, 1 + sr) - Math.atan2(-si, 1 - sr); ph = Math.abs(Math.atan2(Math.sin(ph), Math.cos(ph))); return { att: Math.abs(Math.log(num / den)), ph }; });
}
// a 3D frame model of homog.js (bloch3DModel) in the drawn-model format, so the same FE tools apply
function fromBloch3D(m3, cell) {
  const a1 = cell.a1, a2 = cell.a2, sh = n => [n[0] * a1[0] + n[1] * a2[0], n[0] * a1[1] + n[1] * a2[1], 0];
  const elems = m3.els.map(e => { const Xq = e.X[1]; return { type: 'beam', nd: [[e.p, [0, 0]], [e.q, e.n]], K: e.K, M: e.M, mass: e.mass, role: 'host', X: e.X, Xd: e.X, draw: e.draw }; });
  const seen = new Set(), springs = m3.springs.map(sp => { const first = !seen.has(sp.ri); seen.add(sp.ri); return { node: sp.node, m: sp.m, k: sp.k, n: sp.n, f: Math.sqrt(sp.k / sp.m) / TAU, mr: first ? m3.res[sp.ri].mTot : 0, ri: sp.ri }; });
  const S = Math.abs(a1[0] * a2[1] - a1[1] * a2[0]);
  return { a1, a2, S, nodes: m3.nodes, N: m3.N, elems, pmass: [], springs, Henv: cell.Henv, h: 0, Lref: Math.max(Math.hypot(...a1), Math.hypot(...a2)) };
}
// workbench spec (SI, flat) -> drawing (mm)
function fromSpec(spec) {
  const mm = v => +(v * 1000).toPrecision(9), parts = [];
  for (const b of spec.beams) { const sec = spec.sections[b.sec], R = Eg.shiftVec(spec, b.n), p0 = spec.nodes[b.i], q = spec.nodes[b.j];
    parts.push({ kind: 'beam', p: [[mm(p0.x), mm(p0.y), 0], [mm(q.x + R[0]), mm(q.y + R[1]), 0]], b: mm(sec.shape === 'circ' ? sec.h : sec.b), h: mm(sec.h), role: b.role === 'res' ? 'res' : 'host', rho: sec.rho, E_GPa: sec.E ? sec.E / 1e9 : undefined }); }
  spec.nodes.forEach(p => { if (p.m > 0) parts.push({ kind: 'mass', p: [mm(p.x), mm(p.y), 0], m: p.m * 1000, side: mm(p.side || 0), role: p.role === 'res' ? 'res' : 'host' }); });
  for (const sp of spec.springs || []) { const p = spec.nodes[sp.node]; parts.push({ kind: 'spring', p: [mm(p.x), mm(p.y), 0], m: sp.m * 1000, f: +(Math.sqrt(sp.k / sp.m) / TAU).toPrecision(6), dir: [0, 0, 1], role: 'res' }); }
  return { a1: [mm(spec.a1[0]), mm(spec.a1[1])], a2: [mm(spec.a2[0]), mm(spec.a2[1])], mesh: 2, parts };
}

return { polyArea, buildFinite, finiteModes, offsetMat, shapePoly, expandParts, transformPart, symmetry, stripDrawing, buildStrip, stripAt, lowestModesOp, wfeModel, wfeAt, fromBloch3D, fromSpec, shellTri, meshPlate, delaunay, build, blochMats, lowestModes, blochAt, homogenize, resonators, hermChol, hermSolve, pointInPoly, segDist, TAU };
})();
if (typeof module !== 'undefined') module.exports = Drw;
