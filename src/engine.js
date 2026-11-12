// ================= Metamaterial engine: out-of-plane beam lattices =================
const Eng = (() => {
const TAU = 2 * Math.PI;
function mat(n, m) { const A = new Array(n); for (let i = 0; i < n; i++) A[i] = new Float64Array(m); return A; }

function sectionProps(s) {
  if (s.shape === 'circ') { const d = s.h; const A = Math.PI * d * d / 4, I = Math.PI * d ** 4 / 64, J = Math.PI * d ** 4 / 32; return { A, I, Iz: I, J, Ip: J }; }
  const b = s.b, h = s.h, A = b * h, I = b * h ** 3 / 12, Iz = h * b ** 3 / 12;
  const aa = Math.max(b, h), bb = Math.min(b, h);
  const J = aa * bb ** 3 * (1 / 3 - 0.21 * (bb / aa) * (1 - bb ** 4 / (12 * aa ** 4)));
  return { A, I, Iz, J, Ip: b * h * (b * b + h * h) / 12 };
}

// Bending block for one plane, DOFs [v1, th1, v2, th2] (th positive like dv/dx).
// Timoshenko: interdependent interpolation (exact statics, no shear locking), consistent mass with rotary inertia
// (Friedman & Kosmatka 1993), Phi = 12 EI / (kappa G A L^2). Euler-Bernoulli: Phi = 0, no rotary inertia.
function bendBlock(L, EI, kGA, rA, rI, tim) {
  const L2 = L * L, P = tim ? 12 * EI / (kGA * L2) : 0, d = 1 + P, d2 = d * d;
  const kf = EI / (d * L2 * L);
  const K = [[12, 6 * L, -12, 6 * L], [6 * L, (4 + P) * L2, -6 * L, (2 - P) * L2], [-12, -6 * L, 12, -6 * L], [6 * L, (2 - P) * L2, -6 * L, (4 + P) * L2]].map(r => r.map(v => v * kf));
  const P2 = P * P;
  const a = 13 / 35 + 7 * P / 10 + P2 / 3, b = (11 / 210 + 11 * P / 120 + P2 / 24) * L, c = 9 / 70 + 3 * P / 10 + P2 / 6,
        e = (13 / 420 + 3 * P / 40 + P2 / 24) * L, f = (1 / 105 + P / 60 + P2 / 120) * L2, g = (1 / 140 + P / 60 + P2 / 120) * L2;
  const mt = rA * L / d2;
  const M = [[a, b, c, -e], [b, f, e, -g], [c, e, a, -b], [-e, -g, -b, f]].map(r => r.map(v => v * mt));
  if (tim && rI > 0) {
    const r1 = 6 / 5, r2 = (1 / 10 - P / 2) * L, r3 = (2 / 15 + P / 6 + P2 / 3) * L2, r4 = (-1 / 30 - P / 6 + P2 / 6) * L2, mr = rI / (d2 * L);
    const Mr = [[r1, r2, -r1, r2], [r2, r3, -r2, r4], [-r1, -r2, r1, -r2], [r2, r4, -r2, r3]];
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) M[i][j] += mr * Mr[i][j];
  }
  return { K, M };
}
function shearCoef(shape, nu) { return shape === 'circ' ? 6 * (1 + nu) / (7 + 6 * nu) : 10 * (1 + nu) / (12 + 11 * nu); }

// Out-of-plane (grillage) beam, Timoshenko or Euler-Bernoulli bending, with St-Venant torsion.
// Nodal DOFs: [w, dw/dx, dw/dy]. Local: [w, psi (slope along axis), chi (slope across = twist)]
function beamKM(L, c, s, E, G, rho, sec, tim) {
  const GJ = G * sec.J, rIp = rho * sec.Ip;
  const Kl = mat(6, 6), Ml = mat(6, 6), bi = [0, 1, 3, 4];
  const bb = bendBlock(L, E * sec.I, sec.kappa * G * sec.A, rho * sec.A, rho * sec.I, tim);
  for (let r = 0; r < 4; r++) for (let q = 0; q < 4; q++) { Kl[bi[r]][bi[q]] = bb.K[r][q]; Ml[bi[r]][bi[q]] = bb.M[r][q]; }
  const kt = GJ / L, mt = rIp * L / 6;
  Kl[2][2] = kt; Kl[5][5] = kt; Kl[2][5] = -kt; Kl[5][2] = -kt;
  Ml[2][2] = 2 * mt; Ml[5][5] = 2 * mt; Ml[2][5] = mt; Ml[5][2] = mt;
  const T = mat(6, 6);
  for (const o of [0, 3]) { T[o][o] = 1; T[o + 1][o + 1] = c; T[o + 1][o + 2] = s; T[o + 2][o + 1] = -s; T[o + 2][o + 2] = c; }
  return { K: tAT(T, Kl), M: tAT(T, Ml) };
}
// In-plane frame element: axial bar + Euler-Bernoulli bending in the plane.
// Nodal DOFs: [u, v, theta_z]. Local: [u_axial, v_transverse, theta_z]. Rotary inertia rho*Iz lumped.
function beamKMip(L, c, s, E, G, rho, sec, tim) {
  const EA = E * sec.A, rA = rho * sec.A;
  const Kl = mat(6, 6), Ml = mat(6, 6), bi = [1, 2, 4, 5];
  const bb = bendBlock(L, E * sec.Iz, sec.kappa * G * sec.A, rA, rho * sec.Iz, tim);
  for (let r = 0; r < 4; r++) for (let q = 0; q < 4; q++) { Kl[bi[r]][bi[q]] = bb.K[r][q]; Ml[bi[r]][bi[q]] = bb.M[r][q]; }
  const ka = EA / L, ma = rA * L / 6;
  Kl[0][0] = Kl[3][3] = ka; Kl[0][3] = Kl[3][0] = -ka;
  Ml[0][0] = Ml[3][3] = 2 * ma; Ml[0][3] = Ml[3][0] = ma;
  if (!tim) { const jr = rho * sec.Iz * L / 2; Ml[2][2] += jr; Ml[5][5] += jr; }   // Timoshenko has consistent rotary inertia
  const T = mat(6, 6);
  for (const o of [0, 3]) { T[o][o] = c; T[o][o + 1] = s; T[o + 1][o] = -s; T[o + 1][o + 1] = c; T[o + 2][o + 2] = 1; }
  return { K: tAT(T, Kl), M: tAT(T, Ml) };
}
function tAT(T, A) {
  const B = mat(6, 6), C = mat(6, 6);
  for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) { let v = 0; for (let k = 0; k < 6; k++) v += A[i][k] * T[k][j]; B[i][j] = v; }
  for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) { let v = 0; for (let k = 0; k < 6; k++) v += T[k][i] * B[k][j]; C[i][j] = v; }
  return C;
}

// ---------- reciprocal lattice, Brillouin zone and k-paths ----------
function recip(a1, a2) {
  const det = a1[0] * a2[1] - a1[1] * a2[0];
  return { b1: [TAU * a2[1] / det, -TAU * a2[0] / det], b2: [-TAU * a1[1] / det, TAU * a1[0] / det], area: Math.abs(det) };
}
function bzPolygon(b1, b2) {
  const R = 4 * (Math.hypot(b1[0], b1[1]) + Math.hypot(b2[0], b2[1]));
  let poly = [[-R, -R], [R, -R], [R, R], [-R, R]];
  for (let m = -2; m <= 2; m++) for (let n = -2; n <= 2; n++) {
    if (!m && !n) continue;
    const G = [m * b1[0] + n * b2[0], m * b1[1] + n * b2[1]], c = (G[0] * G[0] + G[1] * G[1]) / 2, out = [];
    for (let i = 0; i < poly.length; i++) {
      const P = poly[i], Q = poly[(i + 1) % poly.length], dp = P[0] * G[0] + P[1] * G[1] - c, dq = Q[0] * G[0] + Q[1] * G[1] - c;
      if (dp <= 0) out.push(P);
      if ((dp < 0 && dq > 0) || (dp > 0 && dq < 0)) { const t = dp / (dp - dq); out.push([P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])]); }
    }
    poly = out;
  }
  const clean = []; const tol = 1e-9 * Math.hypot(b1[0], b1[1]);
  for (const p of poly) if (!clean.some(q => Math.hypot(p[0] - q[0], p[1] - q[1]) < tol)) clean.push(p);
  return clean;
}
const SUB = '₀₁₂₃₄₅₆₇₈₉';
// Which point operations about the cell centre map the infinite pattern onto itself?
function symmetryInfo(spec) {
  const tol = 1e-5 * Math.hypot(spec.a1[0], spec.a1[1]), key = (x, y) => `${Math.round(x / tol)},${Math.round(y / tol)}`;
  const near = (x, y) => { const i = Math.round(x / tol), j = Math.round(y / tol), out = []; for (let di = -1; di <= 1; di++) for (let dj = -1; dj <= 1; dj++) out.push(`${i + di},${j + dj}`); return out; };
  const segs = new Set(), pts = new Set(), base = [];
  const C = spec.centre;
  for (let i = -3; i <= 3; i++) for (let j = -3; j <= 3; j++) {
    const R = shiftVec(spec, [i, j]);
    spec.beams.forEach(b => {
      const Rj = shiftVec(spec, b.n), p = [spec.nodes[b.i].x + R[0] - C[0], spec.nodes[b.i].y + R[1] - C[1]], q = [spec.nodes[b.j].x + Rj[0] + R[0] - C[0], spec.nodes[b.j].y + Rj[1] + R[1] - C[1]];
      const sc = (spec.sections || {})[b.sec] || {}, tag = `${sc.shape}:${(+sc.b).toPrecision(6)}:${(+sc.h).toPrecision(6)}|${b.role}`, k1 = key(...p), k2 = key(...q);
      segs.add(tag + '|' + (k1 < k2 ? k1 + '|' + k2 : k2 + '|' + k1));
      if (Math.abs(i) <= 1 && Math.abs(j) <= 1) base.push({ p, q, tag });
    });
    spec.nodes.forEach(nd => { if (nd.m > 0) pts.add('m' + (+nd.m).toPrecision(6) + '|' + key(nd.x + R[0] - C[0], nd.y + R[1] - C[1])); });
    spec.springs.forEach(sp => { const nd = spec.nodes[sp.node]; pts.add('s' + (+sp.k).toPrecision(6) + '|' + key(nd.x + R[0] - C[0], nd.y + R[1] - C[1])); });
  }
  const basePts = [];
  for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) { const R = shiftVec(spec, [i, j]);
    spec.nodes.forEach(nd => { if (nd.m > 0) basePts.push({ t: 'm' + (+nd.m).toPrecision(6), x: nd.x + R[0] - C[0], y: nd.y + R[1] - C[1] }); });
    spec.springs.forEach(sp => { const nd = spec.nodes[sp.node]; basePts.push({ t: 's' + (+sp.k).toPrecision(6), x: nd.x + R[0] - C[0], y: nd.y + R[1] - C[1] }); }); }
  const test = T => base.every(s => { const A = near(...T(s.p)), Bn = near(...T(s.q)); return A.some(k1 => Bn.some(k2 => segs.has(s.tag + '|' + (k1 < k2 ? k1 + '|' + k2 : k2 + '|' + k1)))); })
                    && basePts.every(p => near(...T([p.x, p.y])).some(k => pts.has(p.t + '|' + k)));
  const rot = a => p => [p[0] * Math.cos(a) - p[1] * Math.sin(a), p[0] * Math.sin(a) + p[1] * Math.cos(a)];
  return { mx: test(p => [-p[0], p[1]]), my: test(p => [p[0], -p[1]]), r90: test(rot(Math.PI / 2)), r60: test(rot(Math.PI / 3)) };
}
// user-defined zone path: {units:'frac'|'cart', pts:[{label, x, y}]}; frac = multiples of b1, b2; cart = (kx, ky) in units of π/a, a = |a1|
let customPath = null;
function setCustomPath(C) { customPath = C; }
function parseKPath(text) {
  const pts = [];
  for (const raw of String(text || '').split(/\n|;/)) { const line = raw.trim(); if (!line || line.startsWith('#')) continue;
    const m = line.match(/^(\S+)[\s,]+([-+]?[0-9.]+(?:e[-+]?\d+)?(?:\s*\/\s*[0-9.]+)?)[\s,]+([-+]?[0-9.]+(?:e[-+]?\d+)?(?:\s*\/\s*[0-9.]+)?)$/i);
    if (!m) throw new Error(`Custom k-path: cannot read “${line}”. Write one point per line: label x y (e.g. “X 0.5 0” or “K 2/3 1/3”).`);
    const v = t => { const [a, b] = t.split('/').map(Number); return b ? a / b : a; };
    pts.push({ label: m[1], x: v(m[2]), y: v(m[3]) }); }
  if (pts.length < 2) throw new Error('Custom k-path: give at least two points.');
  for (let i = 1; i < pts.length; i++) if (pts[i].x === pts[i - 1].x && pts[i].y === pts[i - 1].y) throw new Error(`Custom k-path: points ${i} and ${i + 1} are the same.`);
  return pts;
}
function makePath(Lat, mode) {
  const { b1, b2, bz } = Lat, comb = (p, q) => [p * b1[0] + q * b2[0], p * b1[1] + q * b2[1]];
  if (mode === 'custom' && customPath && customPath.pts && customPath.pts.length >= 2) {
    const a = Math.hypot(Lat.a1[0], Lat.a1[1]), toK = q => customPath.units === 'cart' ? [q.x * Math.PI / a, q.y * Math.PI / a] : comb(q.x, q.y);
    return { labels: customPath.pts.map(q => q.label), pts: customPath.pts.map(toK), mode: 'custom' };
  }
  if (mode === 'custom') mode = 'auto';
  const sy = Lat.sym || { mx: true, my: true, r90: true, r60: true };
  const forced = mode === 'stdforce'; if (forced) mode = 'std';
  if (!forced && mode === 'std' && !((Lat.std === 'hex' && sy.r60 && (sy.mx || sy.my)) || (Lat.std === 'square' && sy.r90 && (sy.mx || sy.my)))) mode = 'auto'; // the standard path only where the cell has the symmetry it assumes
  if (mode === 'auto') {
    if (Lat.std === 'square' && sy.r90 && sy.mx) mode = 'std';
    else if (Lat.std === 'hex' && sy.r60 && (sy.mx || sy.my)) mode = 'std';
    else mode = (sy.mx && sy.my) ? 'quarter' : 'half';
  }
  if (mode === 'std' && Lat.std === 'square') return { labels: ['Γ', 'X', 'M', 'Γ'], pts: [[0, 0], comb(0.5, 0), comb(0.5, 0.5), [0, 0]], mode };
  if (mode === 'std' && Lat.std === 'hex') return { labels: ['Γ', 'X', 'K', 'Γ'], pts: [[0, 0], comb(0.5, 0), comb(2 / 3, 1 / 3), [0, 0]], mode };
  if (mode === 'std') mode = 'quarter';
  const eps = 1e-7 * Math.hypot(b1[0], b1[1]);
  const hit = d => { let best = Infinity; for (let i = 0; i < bz.length; i++) { const P = bz[i], Q = bz[(i + 1) % bz.length], ex = Q[0] - P[0], ey = Q[1] - P[1], den = d[0] * ey - d[1] * ex; if (Math.abs(den) < 1e-30) continue; const s = (P[0] * ey - P[1] * ex) / den, u = (P[0] * d[1] - P[1] * d[0]) / den; if (s > 0 && u >= -1e-9 && u <= 1 + 1e-9) best = Math.min(best, s); } return [d[0] * best, d[1] * best]; };
  const ang = p => Math.atan2(p[1], p[0]);
  const vl = n => (k, i) => n > 1 ? 'V' + SUB[i + 1] : 'V';
  if (mode === 'half') {
    const X = hit([1, 0]), Xm = hit([-1, 0]), V = bz.filter(p => p[1] > eps).sort((p, q) => ang(p) - ang(q));
    return { labels: ['Γ', 'X', ...V.map(vl(V.length)), 'X′', 'Γ'], pts: [[0, 0], X, ...V, Xm, [0, 0]], mode };
  }
  const X = hit([1, 0]), Y = hit([0, 1]), V = bz.filter(p => p[0] > eps && p[1] > eps).sort((p, q) => ang(p) - ang(q));
  return { labels: ['Γ', 'X', ...V.map(vl(V.length)), 'Y', 'Γ'], pts: [[0, 0], X, ...V, Y, [0, 0]], mode: 'quarter' };
}

// ---------- parametric unit cells (SI units) ----------
// spec: {a1,a2,centre,std,nodes:[{x,y,role,m,J,side}],beams:[{i,j,n,sec,role,wall}],springs:[{node,k,m}],sections:{name:{shape,b,h}},walls:[{p0,p1,wall}]}
function baseCell(G) {
  const s3 = Math.sqrt(3), N = (x, y) => ({ x, y, role: 'host', m: 0, J: 0, side: 0 });
  let a1, a2, nodes, beams, centre, std = null;
  const B = (i, j, n, wall, sec = 'host') => ({ i, j, n, sec, role: 'host', wall });
  if (G.lattice === 'square') { const a = G.a; a1 = [a, 0]; a2 = [0, a]; nodes = [N(0, 0)]; beams = [B(0, 0, [1, 0], 1), B(0, 0, [0, 1], 2)]; centre = [a / 2, a / 2]; std = 'square'; }
  else if (G.lattice === 'triangular') { const a = G.a; a1 = [a, 0]; a2 = [a / 2, a * s3 / 2]; nodes = [N(0, 0)]; beams = [B(0, 0, [1, 0], 1), B(0, 0, [0, 1], 2), B(0, 0, [-1, 1], 3)]; centre = [(a1[0] + a2[0]) / 3, (a1[1] + a2[1]) / 3]; std = 'hex'; }
  else if (G.lattice === 'kagome') { const L = G.a / 2; a1 = [G.a, 0]; a2 = [G.a / 2, G.a * s3 / 2]; nodes = [N(0, 0), N(L, 0), N(L / 2, s3 * L / 2)];
    beams = [B(0, 1, [0, 0], 1), B(0, 2, [0, 0], 2), B(1, 2, [0, 0], 3), B(1, 0, [1, 0], 4), B(2, 0, [0, 1], 5), B(2, 1, [-1, 1], 6)]; centre = [(a1[0] + a2[0]) / 2, (a1[1] + a2[1]) / 2]; std = 'hex'; }
  else { // hexagon (honeycomb) with vertical walls h, inclined walls l at angle theta (deg); theta < 0 is re-entrant
    const l = G.l, h = G.h, th = G.theta * Math.PI / 180, c = Math.cos(th), s = Math.sin(th);
    a1 = [2 * l * c, 0]; a2 = [l * c, h + l * s]; nodes = [N(0, 0), N(0, h)];
    beams = [B(0, 1, [0, 0], 1, 'hostV'), B(1, 0, [0, 1], 2), B(1, 0, [-1, 1], 3)]; centre = [l * c, h / 2];
    if (Math.abs(h - l) < 1e-9 * l && Math.abs(G.theta - 30) < 1e-6) std = 'hex';
  }
  const spec = { a1, a2, centre, std, nodes, beams, springs: [] };
  spec.walls = beams.map(b => { const p0 = [nodes[b.i].x, nodes[b.i].y], R = shiftVec(spec, b.n); return { p0, p1: [nodes[b.j].x + R[0], nodes[b.j].y + R[1]], wall: b.wall }; });
  return spec;
}
function shiftVec(spec, n) { return [n[0] * spec.a1[0] + n[1] * spec.a2[0], n[0] * spec.a1[1] + n[1] * spec.a2[1]]; }
function checkBase(G) {
  if (G.lattice === 'hexagon') {
    if (!(Math.abs(G.theta) < 85)) return 'Cell angle must lie between −85° and 85°.';
    const s = Math.sin(G.theta * Math.PI / 180);
    if (!(G.h + 2 * G.l * s > 1e-9 * G.l)) return 'These wall lengths and angle make the re-entrant cell close on itself: increase h or reduce |θ|.';
  }
  return null;
}
// Resonators. R: {type:'none'|'spring'|'arm'|'central', attach:'junction'|'mid'|'t', t, walls:null|[numbers], orient:'in'|'out'|'both'|'angle', angle,
//  shape:'straight'|'L', L1, L2, turn, mt, side, k, m, M, Mside, ligTo:'corners'|'mids'}
function addResonators(base, R) {
  const spec = JSON.parse(JSON.stringify(base));
  if (!R || R.type === 'none') return spec;
  const nodes = spec.nodes, splits = new Map();
  const wallBeams = spec.beams.map((b, k) => k);
  const wallGeom = k => { const b = spec.beams[k], p0 = [nodes[b.i].x, nodes[b.i].y], R0 = shiftVec(spec, b.n), p1 = [nodes[b.j].x + R0[0], nodes[b.j].y + R0[1]]; return { p0, p1 }; };
  const split = (k, t) => {
    const b = spec.beams[k];
    if (t <= 1e-6) return b.i;
    if (!splits.has(k)) splits.set(k, []);
    const list = splits.get(k), ex = list.find(q => Math.abs(q.t - t) < 1e-6); if (ex) return ex.node;
    const g = wallGeom(k); nodes.push({ x: g.p0[0] + t * (g.p1[0] - g.p0[0]), y: g.p0[1] + t * (g.p1[1] - g.p0[1]), role: 'host', m: 0, J: 0, side: 0 });
    list.push({ t, node: nodes.length - 1 }); return nodes.length - 1;
  };
  const C = spec.centre, scale = Math.hypot(spec.a1[0], spec.a1[1]);
  const nearestCentre = p => { let best = null, bd = Infinity; for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) { const Rr = shiftVec(spec, [i, j]), c = [C[0] + Rr[0], C[1] + Rr[1]]; const d = Math.hypot(p[0] - c[0], p[1] - c[1]) - (i === 0 && j === 0 ? 1e-7 * scale : 0); if (d < bd - 1e-9 * scale) { bd = d; best = c; } } return best; };
  const unit = v => { const L = Math.hypot(v[0], v[1]); return [v[0] / L, v[1] / L]; };
  const rot = (v, deg) => { const a = deg * Math.PI / 180; return [v[0] * Math.cos(a) - v[1] * Math.sin(a), v[0] * Math.sin(a) + v[1] * Math.cos(a)]; };
  const wantWall = k => !R.walls || R.walls.includes(spec.beams[k].wall);
  const sites = []; // {node, pos, dirs:[...]}
  if (R.type === 'spring' || R.type === 'arm') {
    if (R.attach === 'junction') {
      base.nodes.forEach((p, i) => { const pos = [p.x, p.y], c = nearestCentre(pos); let din = unit([c[0] - pos[0], c[1] - pos[1]]); sites.push({ node: i, pos, din }); });
    } else {
      const t = R.attach === 'mid' ? 0.5 : R.t;
      for (const k of wallBeams) if (wantWall(k)) {
        const g = wallGeom(k), d = unit([g.p1[0] - g.p0[0], g.p1[1] - g.p0[1]]), pos = [g.p0[0] + t * (g.p1[0] - g.p0[0]), g.p0[1] + t * (g.p1[1] - g.p0[1])];
        let nrm = [-d[1], d[0]]; const c = nearestCentre(pos); if (nrm[0] * (c[0] - pos[0]) + nrm[1] * (c[1] - pos[1]) < 0) nrm = [-nrm[0], -nrm[1]];
        sites.push({ node: split(k, t), pos, din: nrm, wallDir: d, offset: R.wallHalf ? R.wallHalf(spec.beams[k].sec) : 0 });
      }
    }
  }
  if (R.type === 'spring') for (const s of sites) spec.springs.push({ node: s.node, k: R.k, m: R.m });
  if (R.type === 'arm') {
    for (const s of sites) {
      let dirs;
      if (R.orient === 'angle') dirs = [[Math.cos(R.angle * Math.PI / 180), Math.sin(R.angle * Math.PI / 180)]];
      else if (R.orient === 'out') dirs = [[-s.din[0], -s.din[1]]];
      else if (R.orient === 'both') dirs = [s.din, [-s.din[0], -s.din[1]]];
      else dirs = [s.din];
      dirs.forEach((d0, side) => {
        const segs = [{ len: R.L1, d: d0 }];
        if (R.shape === 'L' && R.L2 > 0) segs.push({ len: R.L2, d: rot(d0, side ? -R.turn : R.turn) });
        const block = R.tipForm === 'block';
        if (block && R.blockPos === 'centred') { segs[segs.length - 1].len -= R.bL / 2; if (!(segs[segs.length - 1].len > 1e-6)) throw new Error('The tip block is longer than the arm can carry when centred on its tip: shorten the block or lengthen the arm.'); }
        let prev = s.node, pos = s.pos.slice();
        if (s.offset > 0) { // the part of the arm inside the wall is wall material: a stiff stub of wall height up to the wall face
          pos = [pos[0] + s.offset * d0[0], pos[1] + s.offset * d0[1]];
          nodes.push({ x: pos[0], y: pos[1], role: 'host', m: 0, J: 0, side: 0 });
          spec.beams.push({ i: prev, j: nodes.length - 1, n: [0, 0], sec: 'stub', role: 'host' }); prev = nodes.length - 1;
        }
        for (const sg of segs) {
          pos = [pos[0] + sg.len * sg.d[0], pos[1] + sg.len * sg.d[1]];
          nodes.push({ x: pos[0], y: pos[1], role: 'res', m: 0, J: 0, side: 0 });
          spec.beams.push({ i: prev, j: nodes.length - 1, n: [0, 0], sec: 'arm', role: 'res' }); prev = nodes.length - 1;
        }
        if (block) {
          const d = segs[segs.length - 1].d; pos = [pos[0] + R.bL * d[0], pos[1] + R.bL * d[1]];
          nodes.push({ x: pos[0], y: pos[1], role: 'res', m: 0, J: 0, side: 0 });
          spec.beams.push({ i: prev, j: nodes.length - 1, n: [0, 0], sec: 'tip', role: 'res' });
        } else Object.assign(nodes[prev], { m: R.mt, J: R.mt * R.side * R.side / 6, side: R.side });
      });
    }
  }
  if (R.type === 'central') {
    nodes.push({ x: C[0], y: C[1], role: 'res', m: R.M, J: R.M * R.Mside * R.Mside / 6, side: R.Mside }); const cn = nodes.length - 1;
    const cand = [];
    if (R.ligTo === 'mids') { for (const k of wallBeams) { const g = wallGeom(k), mid = [(g.p0[0] + g.p1[0]) / 2, (g.p0[1] + g.p1[1]) / 2]; for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) { const Rr = shiftVec(spec, [i, j]); cand.push({ k, n: [i, j], d: Math.hypot(mid[0] + Rr[0] - C[0], mid[1] + Rr[1] - C[1]) }); } } }
    else base.nodes.forEach((p, idx) => { for (let i = -2; i <= 2; i++) for (let j = -2; j <= 2; j++) { const Rr = shiftVec(spec, [i, j]); cand.push({ node: idx, n: [i, j], d: Math.hypot(p.x + Rr[0] - C[0], p.y + Rr[1] - C[1]) }); } });
    const dmin = Math.min(...cand.map(c => c.d));
    for (const c of cand) if (c.d < dmin * (1 + 1e-6)) {
      const j = c.node != null ? c.node : split(c.k, 0.5);
      spec.beams.push({ i: cn, j, n: c.n, sec: 'lig', role: 'res' });
    }
  }
  // rebuild split host walls
  const out = [];
  spec.beams.forEach((b, k) => {
    if (!splits.has(k)) { out.push(b); return; }
    const pts = [{ node: b.i, n: [0, 0] }, ...splits.get(k).sort((p, q) => p.t - q.t).map(q => ({ node: q.node, n: [0, 0] })), { node: b.j, n: b.n }];
    for (let s = 0; s < pts.length - 1; s++) out.push({ i: pts[s].node, j: pts[s + 1].node, n: [pts[s + 1].n[0] - pts[s].n[0], pts[s + 1].n[1] - pts[s].n[1]], sec: b.sec, role: b.role, wall: b.wall });
  });
  spec.beams = out;
  return spec;
}
function scaleSpec(spec, f) {
  const s = JSON.parse(JSON.stringify(spec)), m = v => [v[0] * f, v[1] * f];
  s.a1 = m(s.a1); s.a2 = m(s.a2); s.centre = m(s.centre); s.nodes.forEach(p => { p.x *= f; p.y *= f; });
  if (s.walls) s.walls.forEach(w => { w.p0 = m(w.p0); w.p1 = m(w.p1); });
  return s;
}
// JSON exchange in mm, g, N/m
function specToJSON(spec) {
  const r = v => +(+v).toPrecision(7), mm = v => r(v * 1000), pair = v => [mm(v[0]), mm(v[1])];
  const o = { units: 'mm, g, g·mm², N/m', a1: pair(spec.a1), a2: pair(spec.a2), centre: pair(spec.centre), symmetry: spec.std || 'none', sections: {} };
  for (const k in spec.sections) { const s = spec.sections[k], q = s.shape === 'circ' ? { shape: 'circ', d: mm(s.h) } : { shape: 'rect', b: mm(s.b), h: mm(s.h) }; if (s.rho) q.rho = r(s.rho); if (s.E) q.E_GPa = r(s.E / 1e9); o.sections[k] = q; }
  o.nodes = spec.nodes.map(p => { const q = { x: mm(p.x), y: mm(p.y), role: p.role }; if (p.m) { q.mass = r(p.m * 1000); q.J = r(p.J * 1e9); q.side = mm(p.side); } return q; });
  o.beams = spec.beams.map(b => { const q = { i: b.i, j: b.j, shift: b.n, section: b.sec, role: b.role }; if (b.wall) q.wall = b.wall; return q; });
  o.springs = spec.springs.map(s => ({ node: s.node, k: r(s.k), mass: r(s.m * 1000) }));
  const line = x => JSON.stringify(x);
  return '{\n' + Object.keys(o).map(k => {
    const v = o[k];
    if (Array.isArray(v) && v.length && typeof v[0] === 'object') return ` "${k}": [\n${v.map(e => '  ' + line(e)).join(',\n')}\n ]`;
    if (k === 'sections') return ` "sections": {\n${Object.keys(v).map(s => `  "${s}": ${line(v[s])}`).join(',\n')}\n }`;
    return ` "${k}": ${line(v)}`;
  }).join(',\n') + '\n}';
}
function specFromJSON(txt) {
  let o; try { o = JSON.parse(txt); } catch (e) { throw new Error('The cell definition is not valid JSON: ' + e.message); }
  const vec = (v, name) => { if (!Array.isArray(v) || v.length !== 2 || !v.every(Number.isFinite)) throw new Error(`"${name}" must be a pair of numbers in mm.`); return [v[0] / 1000, v[1] / 1000]; };
  const a1 = vec(o.a1, 'a1'), a2 = vec(o.a2, 'a2');
  if (Math.abs(a1[0] * a2[1] - a1[1] * a2[0]) < 1e-12) throw new Error('a1 and a2 must not be parallel.');
  if (!o.sections || typeof o.sections !== 'object') throw new Error('Add a "sections" object, for example {"host": {"shape": "rect", "b": 4, "h": 2}}.');
  const sections = {};
  for (const k in o.sections) { const s = o.sections[k]; if (s.shape === 'circ') { if (!(s.d > 0)) throw new Error(`Section "${k}" needs a diameter d > 0.`); sections[k] = { shape: 'circ', b: s.d / 1000, h: s.d / 1000 }; } else { if (!(s.b > 0 && s.h > 0)) throw new Error(`Section "${k}" needs b > 0 and h > 0.`); sections[k] = { shape: 'rect', b: s.b / 1000, h: s.h / 1000 }; }
    if (s.rho != null) { if (!(s.rho > 0)) throw new Error(`Section "${k}": rho must be positive (kg/m³).`); sections[k].rho = s.rho; }
    if (s.E_GPa != null) { if (!(s.E_GPa > 0)) throw new Error(`Section "${k}": E_GPa must be positive.`); sections[k].E = s.E_GPa * 1e9; } }
  if (!Array.isArray(o.nodes) || !o.nodes.length) throw new Error('Add at least one node.');
  const nodes = o.nodes.map((p, i) => { if (!Number.isFinite(p.x) || !Number.isFinite(p.y)) throw new Error(`Node ${i} needs numeric x and y.`); const m = (p.mass || 0) / 1000; return { x: p.x / 1000, y: p.y / 1000, role: p.role === 'res' ? 'res' : 'host', m, J: p.J != null ? p.J * 1e-9 : m * ((p.side || 0) / 1000) ** 2 / 6, side: (p.side || 0) / 1000 }; });
  if (!Array.isArray(o.beams) || !o.beams.length) throw new Error('Add at least one beam.');
  const beams = o.beams.map((b, k) => {
    const n = b.shift || [0, 0];
    if (!(Number.isInteger(b.i) && Number.isInteger(b.j) && b.i >= 0 && b.j >= 0 && b.i < nodes.length && b.j < nodes.length)) throw new Error(`Beam ${k} refers to a node that does not exist.`);
    if (!(Array.isArray(n) && n.length === 2 && n.every(Number.isInteger))) throw new Error(`Beam ${k}: "shift" must be two integers, the cell offset of node j.`);
    if (!sections[b.section]) throw new Error(`Beam ${k} uses section "${b.section}", which is not defined.`);
    return { i: b.i, j: b.j, n, sec: b.section, role: b.role === 'res' ? 'res' : 'host', wall: b.wall };
  });
  const springs = (o.springs || []).map((s, k) => { if (!(s.node >= 0 && s.node < nodes.length && s.k > 0 && s.mass > 0)) throw new Error(`Spring ${k} needs a valid node, k > 0 and mass > 0.`); return { node: s.node, k: s.k, m: s.mass / 1000 }; });
  const centre = o.centre ? vec(o.centre, 'centre') : [nodes.reduce((s, p) => s + p.x, 0) / nodes.length, nodes.reduce((s, p) => s + p.y, 0) / nodes.length];
  const spec = { a1, a2, centre, std: (o.symmetry === 'square' || o.symmetry === 'hex') ? o.symmetry : null, nodes, beams, springs, sections };
  spec.walls = beams.filter(b => b.role === 'host' && b.wall).map(b => { const R = shiftVec(spec, b.n); return { p0: [nodes[b.i].x, nodes[b.i].y], p1: [nodes[b.j].x + R[0], nodes[b.j].y + R[1]], wall: b.wall }; });
  return spec;
}

// ---------- FE model from a spec ----------
// opt: {ne, scale:{hostE,hostRho,resK,resM}, path}
function buildModel(spec, matP, opt) {
  const sc = Object.assign({ hostE: 1, hostRho: 1, resK: 1, resM: 1 }, opt.scale || {});
  const family = opt.family === 'ip' ? 'ip' : 'oop', ip = family === 'ip', tim = opt.beam !== 'eb';
  const rc = recip(spec.a1, spec.a2);
  const Lat = { a1: spec.a1, a2: spec.a2, b1: rc.b1, b2: rc.b2, area: rc.area, std: spec.std, centre: spec.centre, hex: spec.std === 'hex' };
  Lat.bz = bzPolygon(rc.b1, rc.b2); Lat.sym = symmetryInfo(spec); Lat.path = makePath(Lat, opt.path || 'auto');
  const secs = {}; for (const k in spec.sections) { secs[k] = sectionProps(spec.sections[k]); secs[k].kappa = shearCoef(spec.sections[k].shape, matP.nu); }
  const nodes = spec.nodes.map(p => ({ x: p.x, y: p.y, kind: p.role === 'res' ? 'res' : 'junction', role: p.role }));
  const geo = spec.beams.map(b => { const R = shiftVec(spec, b.n), p0 = [spec.nodes[b.i].x, spec.nodes[b.i].y], p1 = [spec.nodes[b.j].x + R[0], spec.nodes[b.j].y + R[1]]; return { p0, p1, len: Math.hypot(p1[0] - p0[0], p1[1] - p0[1]) }; });
  const hl = geo.filter((g, k) => spec.beams[k].role === 'host').map(g => g.len), Lref = (hl.length ? hl : geo.map(g => g.len)).reduce((s, v) => s + v, 0) / Math.max(1, (hl.length || geo.length));
  const elems = [], chains = []; let hostMass = 0, resMass = 0;
  spec.beams.forEach((b, k) => {
    const g = geo[k]; if (!(g.len > 1e-9)) throw new Error(`Beam ${k} has zero length.`);
    const sec = secs[b.sec]; if (!sec) throw new Error(`Section "${b.sec}" is not defined.`);
    const sraw = spec.sections[b.sec] || {}, res = b.role === 'res', E = (sraw.E || matP.E) * (res ? sc.resK : sc.hostE), rho = (sraw.rho || matP.rho) * (res ? sc.resM : sc.hostRho), G = E / (2 * (1 + (sraw.nu || matP.nu)));
    const ne = Math.max(2, Math.min(16, Math.round(opt.ne * g.len / Lref)));
    const dx = g.p1[0] - g.p0[0], dy = g.p1[1] - g.p0[1], c = dx / g.len, s = dy / g.len;
    const chain = [b.i];
    for (let t = 1; t < ne; t++) { nodes.push({ x: g.p0[0] + dx * t / ne, y: g.p0[1] + dy * t / ne, kind: res ? 'res' : 'beam', role: b.role }); chain.push(nodes.length - 1); }
    chain.push(b.j);
    const km = ip ? beamKMip(g.len / ne, c, s, E, G, rho, sec, tim) : beamKM(g.len / ne, c, s, E, G, rho, sec, tim), eidx = [];
    for (let t = 0; t < ne; t++) { elems.push({ p: chain[t], q: chain[t + 1], n: t === ne - 1 ? b.n : [0, 0], K: km.K, M: km.M, res }); eidx.push(elems.length - 1); }
    chains.push({ beam: b, chain, eidx, len: g.len, c, s, p0: g.p0, p1: g.p1, role: b.role, sec: spec.sections[b.sec] });
    if (res) resMass += rho * sec.A * g.len; else hostMass += rho * sec.A * g.len;
  });
  const pmass = [], massNodes = [];
  spec.nodes.forEach((p, i) => {
    if (!(p.m > 0)) return;
    const f = p.role === 'res' ? sc.resM : sc.hostRho, m = p.m * f, J = (p.J || 0) * f;
    if (ip) { pmass.push({ dof: 3 * i, m }, { dof: 3 * i + 1, m }); if (J > 0) pmass.push({ dof: 3 * i + 2, m: J }); }
    else { pmass.push({ dof: 3 * i, m }); if (J > 0) pmass.push({ dof: 3 * i + 1, m: J }, { dof: 3 * i + 2, m: J }); }
    massNodes.push({ node: i, side: p.side, m }); if (p.role === 'res') resMass += m; else hostMass += m;
  });
  // a mass on a spring acts along z for out-of-plane waves; in-plane it is taken as isotropic (same k along x and y)
  const springs = [];
  for (const s of spec.springs) for (const c of (ip ? [0, 1] : [0])) springs.push({ node: s.node, k: s.k * sc.resK, m: s.m * sc.resM, c });
  for (const sp of springs) if (sp.c === 0) resMass += sp.m;
  const N = nodes.length, n = 3 * N + springs.length;
  const deg = new Int32Array(N); for (const e of elems) { deg[e.p]++; deg[e.q]++; }
  for (let i = 0; i < N; i++) if (!deg[i]) throw new Error(`Node ${i} is not connected to any beam.`);
  springs.forEach((sp, i) => sp.dof = 3 * N + i);
  const resDofs = [];
  nodes.forEach((p, i) => { if (p.role === 'res') resDofs.push(3 * i, 3 * i + 1, 3 * i + 2); });
  for (const sp of springs) resDofs.push(sp.dof);
  const K0 = mat(n, n), M0 = mat(n, n), shifted = [];
  for (const e of elems) { if (e.n[0] || e.n[1]) { shifted.push(e); continue; } addElem(K0, M0, e, 1, 0); }
  for (const sp of springs) { const w = 3 * sp.node + sp.c; K0[w][w] += sp.k; K0[sp.dof][sp.dof] += sp.k; K0[w][sp.dof] -= sp.k; K0[sp.dof][w] -= sp.k; M0[sp.dof][sp.dof] += sp.m; }
  for (const pm of pmass) M0[pm.dof][pm.dof] += pm.m;
  return { spec, family, beam: tim ? 'timoshenko' : 'eb', P: { a: Lref }, Lat, nodes, elems, chains, springs, pmass, massNodes, resDofs, N, n, K0, M0, shifted, hostMass, resMass };
}
function addElem(K, M, e, cr, ci, Ki, Mi) {
  const dp = 3 * e.p, dq = 3 * e.q;
  for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) {
    K[dp + a][dp + b] += e.K[a][b]; M[dp + a][dp + b] += e.M[a][b];
    K[dq + a][dq + b] += e.K[3 + a][3 + b]; M[dq + a][dq + b] += e.M[3 + a][3 + b];
    K[dp + a][dq + b] += e.K[a][3 + b] * cr; M[dp + a][dq + b] += e.M[a][3 + b] * cr;
    K[dq + a][dp + b] += e.K[3 + a][b] * cr; M[dq + a][dp + b] += e.M[3 + a][b] * cr;
    if (Ki) {
      Ki[dp + a][dq + b] += e.K[a][3 + b] * ci; Mi[dp + a][dq + b] += e.M[a][3 + b] * ci;
      Ki[dq + a][dp + b] -= e.K[3 + a][b] * ci; Mi[dq + a][dp + b] -= e.M[3 + a][b] * ci;
    }
  }
}
function assembleBloch(md, kx, ky) {
  const n = md.n, Kr = mat(n, n), Mr = mat(n, n), Ki = mat(n, n), Mi = mat(n, n);
  for (let i = 0; i < n; i++) { Kr[i].set(md.K0[i]); Mr[i].set(md.M0[i]); }
  const { a1, a2 } = md.Lat;
  for (const e of md.shifted) {
    const Rx = e.n[0] * a1[0] + e.n[1] * a2[0], Ry = e.n[0] * a1[1] + e.n[1] * a2[1], ph = kx * Rx + ky * Ry;
    addElem(Kr, Mr, e, Math.cos(ph), Math.sin(ph), Ki, Mi);
  }
  return { Kr, Ki, Mr, Mi };
}

// ---------- symmetric eigen (tred2 / tqli, Numerical-Recipes style, 0-based) ----------
function tred2(a, n, d, e, vec, hs) {
  for (let i = n - 1; i > 0; i--) {
    const l = i - 1; let h = 0, scale = 0;
    if (l > 0) {
      for (let k = 0; k <= l; k++) scale += Math.abs(a[i][k]);
      if (scale === 0) e[i] = a[i][l];
      else {
        for (let k = 0; k <= l; k++) { a[i][k] /= scale; h += a[i][k] * a[i][k]; }
        let f = a[i][l]; let g = f >= 0 ? -Math.sqrt(h) : Math.sqrt(h);
        e[i] = scale * g; h -= f * g; a[i][l] = f - g; f = 0;
        for (let j = 0; j <= l; j++) {
          if (vec) a[j][i] = a[i][j] / h;
          g = 0;
          for (let k = 0; k <= j; k++) g += a[j][k] * a[i][k];
          for (let k = j + 1; k <= l; k++) g += a[k][j] * a[i][k];
          e[j] = g / h; f += e[j] * a[i][j];
        }
        const hh = f / (h + h);
        for (let j = 0; j <= l; j++) {
          f = a[i][j]; e[j] = g = e[j] - hh * f;
          for (let k = 0; k <= j; k++) a[j][k] -= (f * e[k] + g * a[i][k]);
        }
      }
    } else e[i] = a[i][l];
    d[i] = h; if (hs) hs[i] = h;
  }
  d[0] = 0; e[0] = 0;
  if (vec) {
    for (let i = 0; i < n; i++) {
      const l = i - 1;
      if (d[i]) for (let j = 0; j <= l; j++) {
        let g = 0;
        for (let k = 0; k <= l; k++) g += a[i][k] * a[k][j];
        for (let k = 0; k <= l; k++) a[k][j] -= g * a[k][i];
      }
      d[i] = a[i][i]; a[i][i] = 1;
      for (let j = 0; j <= l; j++) a[j][i] = a[i][j] = 0;
    }
  } else for (let i = 0; i < n; i++) d[i] = a[i][i];
}
function tqli(d, e, n, z) {
  for (let i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0;
  for (let l = 0; l < n; l++) {
    let iter = 0, m;
    do {
      for (m = l; m < n - 1; m++) { const dd = Math.abs(d[m]) + Math.abs(d[m + 1]); if (Math.abs(e[m]) <= 1e-15 * dd) break; }
      if (m !== l) {
        if (iter++ === 60) throw new Error('tqli: no convergence');
        let g = (d[l + 1] - d[l]) / (2 * e[l]); let r = Math.hypot(g, 1);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? r : -r));
        let s = 1, c = 1, p = 0, i;
        for (i = m - 1; i >= l; i--) {
          let f = s * e[i]; const b = c * e[i];
          e[i + 1] = (r = Math.hypot(f, g));
          if (r === 0) { d[i + 1] -= p; e[m] = 0; break; }
          s = f / r; c = g / r; g = d[i + 1] - p; r = (d[i] - g) * s + 2 * c * b;
          d[i + 1] = g + (p = s * r); g = c * r - b;
          if (z) for (let k = 0; k < n; k++) { f = z[k][i + 1]; z[k][i + 1] = s * z[k][i] + c * f; z[k][i] = c * z[k][i] - s * f; }
        }
        if (r === 0 && i >= l) continue;
        d[l] -= p; e[l] = g; e[m] = 0;
      }
    } while (m !== l);
  }
}

// inverse iteration on the symmetric tridiagonal (diag d, sub-diagonal e[i] between i-1 and i)
function tridiagInv(d, e, n, lam, span) {
  const sig = lam + 1e-11 * (span || 1);
  const dd = new Float64Array(n), du = new Float64Array(n), du2 = new Float64Array(n), dl = new Float64Array(n), mul = new Float64Array(n), piv = new Uint8Array(n);
  for (let i = 0; i < n; i++) { dd[i] = d[i] - sig; du[i] = i < n - 1 ? e[i + 1] : 0; dl[i] = i > 0 ? e[i] : 0; }
  for (let k = 0; k < n - 1; k++) {
    const sub = dl[k + 1];
    if (Math.abs(sub) > Math.abs(dd[k])) {
      const r0 = dd[k], r1 = du[k], r2 = du2[k];
      dd[k] = sub; du[k] = dd[k + 1]; du2[k] = du[k + 1];
      const m = r0 / dd[k]; dd[k + 1] = r1 - m * du[k]; du[k + 1] = r2 - m * du2[k]; piv[k] = 1; mul[k] = m;
    } else { const m = dd[k] === 0 ? 0 : sub / dd[k]; dd[k + 1] -= m * du[k]; piv[k] = 0; mul[k] = m; }
  }
  const tiny = 1e-300 + 1e-15 * (span || 1);
  for (let i = 0; i < n; i++) if (Math.abs(dd[i]) < tiny) dd[i] = tiny;
  let x = new Float64Array(n); for (let i = 0; i < n; i++) x[i] = 1 + 0.37 * Math.sin(1.7 * i + 0.3);
  for (let it = 0; it < 3; it++) {
    const bv = x;
    for (let k = 0; k < n - 1; k++) { if (piv[k]) { const t = bv[k]; bv[k] = bv[k + 1]; bv[k + 1] = t; } bv[k + 1] -= mul[k] * bv[k]; }
    for (let i = n - 1; i >= 0; i--) { let t = bv[i]; if (i + 1 < n) t -= du[i] * bv[i + 1]; if (i + 2 < n) t -= du2[i] * bv[i + 2]; bv[i] = t / dd[i]; }
    let nr = 0; for (let i = 0; i < n; i++) nr += bv[i] * bv[i]; nr = Math.sqrt(nr) || 1; for (let i = 0; i < n; i++) bv[i] /= nr;
  }
  return x;
}
// Hermitian generalized problem K u = w^2 M u via real 2n embedding
function solveK(md, kx, ky, nb, wantVec, wantFrac) {
  const { Kr, Ki, Mr, Mi } = assembleBloch(md, kx, ky);
  const n = md.n, N = 2 * n, A = mat(N, N), B = mat(N, N);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    A[i][j] = A[i + n][j + n] = Kr[i][j]; A[i][j + n] = -Ki[i][j]; A[i + n][j] = Ki[i][j];
    B[i][j] = B[i + n][j + n] = Mr[i][j]; B[i][j + n] = -Mi[i][j]; B[i + n][j] = Mi[i][j];
  }
  // Cholesky (lower)
  for (let j = 0; j < N; j++) {
    const Bj = B[j]; let s = Bj[j]; for (let k = 0; k < j; k++) s -= Bj[k] * Bj[k];
    if (!(s > 0)) throw new Error('Mass matrix is not positive definite: every DOF needs mass.');
    const dj = Math.sqrt(s); Bj[j] = dj;
    for (let i = j + 1; i < N; i++) { const Bi = B[i]; let t = Bi[j]; for (let k = 0; k < j; k++) t -= Bi[k] * Bj[k]; Bi[j] = t / dj; }
  }
  // Y = L^-1 A ; C = L^-1 Y^T
  const Y = mat(N, N);
  for (let i = 0; i < N; i++) { const Li = B[i], Yi = Y[i]; Yi.set(A[i]); for (let k = 0; k < i; k++) { const l = Li[k]; if (l) { const Yk = Y[k]; for (let c = 0; c < N; c++) Yi[c] -= l * Yk[c]; } } const inv = 1 / Li[i]; for (let c = 0; c < N; c++) Yi[c] *= inv; }
  const C = mat(N, N);
  for (let i = 0; i < N; i++) { const Li = B[i], Ci = C[i]; for (let c = 0; c < N; c++) Ci[c] = Y[c][i]; for (let k = 0; k < i; k++) { const l = Li[k]; if (l) { const Ck = C[k]; for (let c = 0; c < N; c++) Ci[c] -= l * Ck[c]; } } const inv = 1 / Li[i]; for (let c = 0; c < N; c++) Ci[c] *= inv; }
  for (let i = 0; i < N; i++) for (let j = 0; j < i; j++) { const v = 0.5 * (C[i][j] + C[j][i]); C[i][j] = C[j][i] = v; }
  const vec = wantVec || wantFrac;
  const d = new Float64Array(N), e = new Float64Array(N), hs = new Float64Array(N);
  tred2(C, N, d, e, false, hs);
  const d0 = d.slice(), e0 = e.slice();
  tqli(d, e, N, null);
  const ev = Array.from(d).sort((p, q) => p - q);
  const out = { f: [], frac: [], vecs: [] };
  const take = Math.min(nb, n);
  for (let b = 0; b < take; b++) {
    const lam = ev[2 * b]; out.f.push(Math.sqrt(Math.max(0, lam)) / TAU);
    if (!vec) continue;
    const z = tridiagInv(d0, e0, N, lam, Math.abs(ev[N - 1]) + Math.abs(ev[0]));
    // apply stored Householder reflections: z <- P_i z for ascending i
    for (let i = 2; i < N; i++) { const h = hs[i]; if (!h) continue; const u = C[i]; let s = 0; for (let k = 0; k < i; k++) s += u[k] * z[k]; s /= h; for (let k = 0; k < i; k++) z[k] -= s * u[k]; }
    const y = z;
    for (let i = N - 1; i >= 0; i--) { let t = y[i]; for (let k = i + 1; k < N; k++) t -= B[k][i] * y[k]; y[i] = t / B[i][i]; }
    if (wantFrac) {
      let part = 0; const R = md.resDofs;
      for (const i of R) for (const j of R) { const m = Mr[i][j]; if (m) part += m * (y[i] * y[j] + y[i + n] * y[j + n]); }
      out.frac.push(Math.min(1, Math.max(0, part)));
    }
    if (wantVec) out.vecs.push({ re: y.slice(0, n), im: y.slice(n) });
  }
  return out;
}

function kPath(Lat, nper) {
  const pts = Lat.path.pts, ks = [], xs = [], ticks = [0], seg = [], lens = [];
  for (let s = 0; s < pts.length - 1; s++) lens.push(Math.hypot(pts[s + 1][0] - pts[s][0], pts[s + 1][1] - pts[s][1]));
  const avg = lens.reduce((a, b) => a + b, 0) / lens.length; let acc = 0;
  for (let s = 0; s < pts.length - 1; s++) {
    const p = pts[s], q = pts[s + 1], m = Math.max(3, Math.round(nper * lens[s] / avg));
    for (let t = (s === 0 ? 0 : 1); t <= m; t++) { const u = t / m; ks.push([p[0] + (q[0] - p[0]) * u, p[1] + (q[1] - p[1]) * u]); xs.push(acc + lens[s] * u); seg.push([s, u]); }
    acc += lens[s]; ticks.push(acc);
  }
  return { ks, xs, ticks, labels: Lat.path.labels, seg };
}

function findGaps(F, nb) { // F[k][b]
  const gaps = [];
  for (let b = 0; b < nb - 1; b++) {
    let hi = -Infinity, lo = Infinity;
    for (const row of F) { if (row[b] > hi) hi = row[b]; if (row[b + 1] < lo) lo = row[b + 1]; }
    if (lo > hi * 1.01 && lo - hi > 1e-9) gaps.push({ lo: hi, hi: lo, below: b + 1 });
  }
  return gaps;
}

// ---------- WFE transfer matrix along a1 (phase along a2 = 0) ----------
function buildWFE(md) {
  const N = md.N, cross = md.chains.filter(ch => ch.beam.n[0] !== 0), nc = cross.length, nb = 3 * nc;
  const cutA = new Map(), cutB = new Map(), elemSide = new Map();
  cross.forEach((ch, r) => {
    const ci = Math.floor(ch.eidx.length / 2), node = ch.chain[ci], Lb = 3 * r, Rb = nb + 3 * r;
    if (ch.beam.n[0] === 1) { cutA.set(node, Rb); cutB.set(node, Lb); } else { cutA.set(node, Lb); cutB.set(node, Rb); }
    ch.eidx.forEach((ei, t) => elemSide.set(ei, t < ci ? 'A' : 'B'));
  });
  const map = new Int32Array(N).fill(-1); let next = 2 * nb;
  for (let i = 0; i < N; i++) if (!cutA.has(i)) { map[i] = next; next += 3; }
  const sBase = next, nw = next + md.springs.length;
  const dof = (node, side) => cutA.has(node) ? (side === 'B' ? cutB.get(node) : cutA.get(node)) : map[node];
  const K = mat(nw, nw), M = mat(nw, nw);
  md.elems.forEach((e, ei) => {
    const side = elemSide.get(ei) || 'A', dp = dof(e.p, side), dq = dof(e.q, side);
    const ix = [dp, dp + 1, dp + 2, dq, dq + 1, dq + 2];
    for (let a = 0; a < 6; a++) for (let b = 0; b < 6; b++) { K[ix[a]][ix[b]] += e.K[a][b]; M[ix[a]][ix[b]] += e.M[a][b]; }
  });
  md.springs.forEach((sp, i) => { const w = dof(sp.node, 'A') + (sp.c || 0), r = sBase + i; K[w][w] += sp.k; K[r][r] += sp.k; K[w][r] -= sp.k; K[r][w] -= sp.k; M[r][r] += sp.m; });
  for (const pm of md.pmass) { const node = Math.floor(pm.dof / 3), c = pm.dof % 3, d0 = dof(node, 'A'); M[d0 + c][d0 + c] += pm.m; }
  // DOF scaling: slopes and rotations * Lref so all DOFs are lengths
  const Lref = md.P.a, S = new Float64Array(nw).fill(1), rotC = md.family === 'ip' ? [2] : [1, 2];
  for (let i = 0; i < N; i++) {
    const bases = cutA.has(i) ? [cutA.get(i), cutB.get(i)] : [map[i]];
    for (const b of bases) for (const c of rotC) S[b + c] = Lref;
  }
  for (let i = 0; i < nw; i++) for (let j = 0; j < nw; j++) { K[i][j] /= S[i] * S[j]; M[i][j] /= S[i] * S[j]; }
  return { K, M, nw, nb, nc };
}
function luSolve(A, Bm, n, m) { // A n×n, Bm n×m, in place; returns X in Bm or null
  for (let k = 0; k < n; k++) {
    let p = k, mx = Math.abs(A[k][k]);
    for (let i = k + 1; i < n; i++) { const v = Math.abs(A[i][k]); if (v > mx) { mx = v; p = i; } }
    if (mx === 0 || !isFinite(mx)) return null;
    if (p !== k) { let t = A[p]; A[p] = A[k]; A[k] = t; t = Bm[p]; Bm[p] = Bm[k]; Bm[k] = t; }
    const Ak = A[k], piv = Ak[k], Bk = Bm[k];
    for (let i = k + 1; i < n; i++) {
      const Ai = A[i], f = Ai[k] / piv; if (!f) continue;
      for (let j = k; j < n; j++) Ai[j] -= f * Ak[j];
      const Bi = Bm[i]; for (let j = 0; j < m; j++) Bi[j] -= f * Bk[j];
    }
  }
  for (let i = n - 1; i >= 0; i--) { const Bi = Bm[i], Ai = A[i]; for (let k = i + 1; k < n; k++) { const a = Ai[k]; if (a) { const Bk = Bm[k]; for (let j = 0; j < m; j++) Bi[j] -= a * Bk[j]; } } for (let j = 0; j < m; j++) Bi[j] /= Ai[i]; }
  return Bm;
}
// --- nonsymmetric eigenvalues: balance + Hessenberg + shifted QR (NR, 1-based internally)
function eigNonsym(A0, n) {
  const a = mat(n + 1, n + 1);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) a[i + 1][j + 1] = A0[i][j];
  const RADIX = 2, sq = 4; let last = 0;
  while (!last) {
    last = 1;
    for (let i = 1; i <= n; i++) {
      let r = 0, c = 0;
      for (let j = 1; j <= n; j++) if (j !== i) { c += Math.abs(a[j][i]); r += Math.abs(a[i][j]); }
      if (c && r) {
        let g = r / RADIX, f = 1; const s = c + r;
        while (c < g) { f *= RADIX; c *= sq; }
        g = r * RADIX;
        while (c > g) { f /= RADIX; c /= sq; }
        if ((c + r) / f < 0.95 * s) { last = 0; g = 1 / f; for (let j = 1; j <= n; j++) a[i][j] *= g; for (let j = 1; j <= n; j++) a[j][i] *= f; }
      }
    }
  }
  for (let m = 2; m < n; m++) {
    let x = 0, i = m;
    for (let j = m; j <= n; j++) if (Math.abs(a[j][m - 1]) > Math.abs(x)) { x = a[j][m - 1]; i = j; }
    if (i !== m) { for (let j = m - 1; j <= n; j++) { const t = a[i][j]; a[i][j] = a[m][j]; a[m][j] = t; } for (let j = 1; j <= n; j++) { const t = a[j][i]; a[j][i] = a[j][m]; a[j][m] = t; } }
    if (x) for (i = m + 1; i <= n; i++) { let y = a[i][m - 1]; if (y !== 0) { y /= x; a[i][m - 1] = y; for (let j = m; j <= n; j++) a[i][j] -= y * a[m][j]; for (let j = 1; j <= n; j++) a[j][m] += y * a[j][i]; } }
  }
  for (let i = 1; i <= n; i++) for (let j = 1; j < i - 1; j++) a[i][j] = 0;
  const wr = new Float64Array(n + 1), wi = new Float64Array(n + 1);
  let anorm = 0;
  for (let i = 1; i <= n; i++) for (let j = Math.max(i - 1, 1); j <= n; j++) anorm += Math.abs(a[i][j]);
  let nn = n, t = 0;
  const SIGN = (x, y) => y >= 0 ? Math.abs(x) : -Math.abs(x);
  while (nn >= 1) {
    let its = 0, l;
    do {
      for (l = nn; l >= 2; l--) { let s = Math.abs(a[l - 1][l - 1]) + Math.abs(a[l][l]); if (s === 0) s = anorm; if (Math.abs(a[l][l - 1]) + s === s) { a[l][l - 1] = 0; break; } }
      let x = a[nn][nn];
      if (l === nn) { wr[nn] = x + t; wi[nn--] = 0; }
      else {
        let y = a[nn - 1][nn - 1], w = a[nn][nn - 1] * a[nn - 1][nn];
        if (l === nn - 1) {
          const p = 0.5 * (y - x), q = p * p + w; let z = Math.sqrt(Math.abs(q)); x += t;
          if (q >= 0) { z = p + SIGN(z, p); wr[nn - 1] = wr[nn] = x + z; if (z) wr[nn] = x - w / z; wi[nn - 1] = wi[nn] = 0; }
          else { wr[nn - 1] = wr[nn] = x + p; wi[nn - 1] = -(wi[nn] = z); }
          nn -= 2;
        } else {
          if (its === 60) return null;
          if (its === 10 || its === 20) { t += x; for (let i = 1; i <= nn; i++) a[i][i] -= x; const s = Math.abs(a[nn][nn - 1]) + Math.abs(a[nn - 1][nn - 2]); y = x = 0.75 * s; w = -0.4375 * s * s; }
          ++its;
          let m, p, q, r, z, s;
          for (m = nn - 2; m >= l; m--) {
            z = a[m][m]; r = x - z; s = y - z;
            p = (r * s - w) / a[m + 1][m] + a[m][m + 1]; q = a[m + 1][m + 1] - z - r - s; r = a[m + 2][m + 1];
            s = Math.abs(p) + Math.abs(q) + Math.abs(r); p /= s; q /= s; r /= s;
            if (m === l) break;
            const u = Math.abs(a[m][m - 1]) * (Math.abs(q) + Math.abs(r));
            const v = Math.abs(p) * (Math.abs(a[m - 1][m - 1]) + Math.abs(z) + Math.abs(a[m + 1][m + 1]));
            if (u + v === v) break;
          }
          for (let i = m + 2; i <= nn; i++) { a[i][i - 2] = 0; if (i !== m + 2) a[i][i - 3] = 0; }
          for (let k = m; k <= nn - 1; k++) {
            if (k !== m) {
              p = a[k][k - 1]; q = a[k + 1][k - 1]; r = 0; if (k !== nn - 1) r = a[k + 2][k - 1];
              if ((x = Math.abs(p) + Math.abs(q) + Math.abs(r)) !== 0) { p /= x; q /= x; r /= x; }
            }
            if ((s = SIGN(Math.sqrt(p * p + q * q + r * r), p)) !== 0) {
              if (k === m) { if (l !== m) a[k][k - 1] = -a[k][k - 1]; } else a[k][k - 1] = -s * x;
              p += s; x = p / s; y = q / s; z = r / s; q /= p; r /= p;
              for (let j = k; j <= nn; j++) { p = a[k][j] + q * a[k + 1][j]; if (k !== nn - 1) { p += r * a[k + 2][j]; a[k + 2][j] -= p * z; } a[k + 1][j] -= p * y; a[k][j] -= p * x; }
              const mmin = nn < k + 3 ? nn : k + 3;
              for (let i = l; i <= mmin; i++) { p = x * a[i][k] + y * a[i][k + 1]; if (k !== nn - 1) { p += z * a[i][k + 2]; a[i][k + 2] -= p * r; } a[i][k + 1] -= p * q; a[i][k] -= p; }
            }
          }
        }
      }
    } while (l < nn - 1);
  }
  const out = []; for (let i = 1; i <= n; i++) out.push([wr[i], wi[i]]); return out;
}
function wfeAt(W, f) {
  const w2 = (TAU * f) ** 2, nw = W.nw, nb = W.nb, b2 = 2 * nb, nI = nw - b2;
  const D = (i, j) => W.K[i][j] - w2 * W.M[i][j];
  const DII = mat(nI, nI), X = mat(nI, b2);
  for (let i = 0; i < nI; i++) { for (let j = 0; j < nI; j++) DII[i][j] = D(b2 + i, b2 + j); for (let j = 0; j < b2; j++) X[i][j] = D(b2 + i, j); }
  if (nI && !luSolve(DII, X, nI, b2)) return null;
  const Dc = mat(b2, b2);
  for (let i = 0; i < b2; i++) for (let j = 0; j < b2; j++) { let v = D(i, j); for (let k = 0; k < nI; k++) v -= D(i, b2 + k) * X[k][j]; Dc[i][j] = v; }
  const blk = (r0, c0) => { const B = mat(nb, nb); for (let i = 0; i < nb; i++) for (let j = 0; j < nb; j++) B[i][j] = Dc[r0 + i][c0 + j]; return B; };
  const DLL = blk(0, 0), DLR = blk(0, nb), DRL = blk(nb, 0), DRR = blk(nb, nb);
  // Quadratic eigenproblem  (D_RL + λ(D_LL+D_RR) + λ² D_LR) u_L = 0, solved through the
  // Möbius map λ = (1+s)/(1-s), which keeps it regular even when D_LR is rank-deficient.
  const Q2 = mat(nb, nb), Q1 = mat(nb, nb), Q0 = mat(nb, nb);
  for (let i = 0; i < nb; i++) for (let j = 0; j < nb; j++) {
    const A2 = DLR[i][j], A1 = DLL[i][j] + DRR[i][j], A0 = DRL[i][j];
    Q2[i][j] = A2 - A1 + A0; Q1[i][j] = 2 * (A2 - A0); Q0[i][j] = A2 + A1 + A0;
  }
  const R = mat(nb, 2 * nb); for (let i = 0; i < nb; i++) for (let j = 0; j < nb; j++) { R[i][j] = Q0[i][j]; R[i][nb + j] = Q1[i][j]; }
  if (!luSolve(Q2, R, nb, 2 * nb)) return null;
  const T = mat(b2, b2);
  for (let i = 0; i < nb; i++) { T[i][nb + i] = 1; for (let j = 0; j < nb; j++) { T[nb + i][j] = -R[i][j]; T[nb + i][nb + j] = -R[i][nb + j]; } }
  for (let i = 0; i < b2; i++) for (let j = 0; j < b2; j++) if (!isFinite(T[i][j])) return null;
  const ev = eigNonsym(T, b2); if (!ev) return null;
  return ev.map(([sr, si]) => {
    const num = Math.hypot(1 + sr, si), den = Math.hypot(1 - sr, si);
    if (den < 1e-14 || num < 1e-14) return { att: 99, ph: 0 };
    let ph = Math.atan2(si, 1 + sr) - Math.atan2(-si, 1 - sr);
    ph = Math.abs(Math.atan2(Math.sin(ph), Math.cos(ph)));
    return { att: Math.abs(Math.log(num / den)), ph };
  });
}
// ---------- strip (ribbon) model for edge waves ----------
function ribbonSpec(spec, N, dir) {
  const swap = dir === 'a2', A1 = swap ? spec.a2 : spec.a1, A2 = swap ? spec.a1 : spec.a2, nn = spec.nodes.length;
  const nodes = [], beams = [], springs = [];
  for (let m = 0; m < N; m++) spec.nodes.forEach(p => nodes.push(Object.assign({}, p, { x: p.x + m * A2[0], y: p.y + m * A2[1] })));
  for (let m = 0; m < N; m++) {
    spec.beams.forEach(b => { const n1 = swap ? b.n[1] : b.n[0], n2 = swap ? b.n[0] : b.n[1], mj = m + n2; if (mj < 0 || mj >= N) return; beams.push(Object.assign({}, b, { i: m * nn + b.i, j: mj * nn + b.j, n: [n1, 0] })); });
    spec.springs.forEach(sp => springs.push(Object.assign({}, sp, { node: m * nn + sp.node })));
  }
  const used = new Uint8Array(nodes.length); beams.forEach(b => { used[b.i] = 1; used[b.j] = 1; });
  const map = new Int32Array(nodes.length).fill(-1), keep = []; nodes.forEach((p, i) => { if (used[i]) { map[i] = keep.length; keep.push(p); } });
  return { a1: A1, a2: [A2[0] * (N + 3), A2[1] * (N + 3)], centre: [spec.centre[0] + (N - 1) / 2 * A2[0], spec.centre[1] + (N - 1) / 2 * A2[1]], std: null,
    nodes: keep, beams: beams.map(b => Object.assign({}, b, { i: map[b.i], j: map[b.j] })), springs: springs.filter(sp => used[sp.node]).map(sp => Object.assign({}, sp, { node: map[sp.node] })),
    sections: spec.sections, walls: [], ribbon: { N, A1, A2 } };
}
function nodeDepth(md, A1, A2) { const det = A1[0] * A2[1] - A1[1] * A2[0]; return md.nodes.map(p => (A1[0] * p.y - A1[1] * p.x) / det); }
function massDiag(md) { const d = new Float64Array(md.n); for (let i = 0; i < md.n; i++) d[i] = md.M0[i][i]; for (const e of md.shifted) for (let a = 0; a < 3; a++) { d[3 * e.p + a] += e.M[a][a]; d[3 * e.q + a] += e.M[3 + a][3 + a]; } return d; }
function edgeShares(md, u, depth, mdiag, w) {
  let smin = Infinity, smax = -Infinity; for (const s of depth) { if (s < smin) smin = s; if (s > smax) smax = s; }
  let tot = 0, bot = 0, top = 0;
  const add = (dof, s) => { const e = mdiag[dof] * (u.re[dof] ** 2 + u.im[dof] ** 2); tot += e; if (s < smin + w) bot += e; if (s > smax - w) top += e; };
  for (let i = 0; i < md.N; i++) for (let a = 0; a < 3; a++) add(3 * i + a, depth[i]);
  for (const sp of md.springs) add(sp.dof, depth[sp.node]);
  return { bot: bot / tot, top: top / tot };
}
return { setCustomPath, parseKPath, luSolve, tred2, tqli, tridiagInv, recip, bzPolygon, shiftVec, bendBlock, shearCoef, ribbonSpec, nodeDepth, massDiag, edgeShares, symmetryInfo, baseCell, checkBase, addResonators, scaleSpec, specToJSON, specFromJSON, buildModel, makePath, solveK, kPath, findGaps, buildWFE, wfeAt, eigNonsym, sectionProps, TAU };
})();
if (typeof module !== 'undefined') module.exports = Eng;
