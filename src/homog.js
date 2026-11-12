// ============ Homogenisation: equivalent Kirchhoff plate + plane-wave expansion with modal resonators ============
// Static energetic homogenisation of a periodic 3D beam frame (6 DOF per node) under the six Kirchhoff-Love
// macro strains, for oblique lattices, per-beam sections
// and Timoshenko beams. Resonators of any geometry are reduced to their clamped-base modes (modal mass, stiffness and
// direction of the effective-mass vector) and coupled to the plate in a plane-wave expansion.
const Hom = (() => {
const Eg = (typeof Eng !== 'undefined') ? Eng : require('./engine.js');
const TAU = 2 * Math.PI;
const zeros = (n, m) => { const A = new Array(n); for (let i = 0; i < n; i++) A[i] = new Float64Array(m); return A; };
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm3 = a => Math.hypot(a[0], a[1], a[2]);
const unit = a => { const L = norm3(a); if (!(L > 0)) throw new Error('A direction vector has zero length.'); return [a[0] / L, a[1] / L, a[2] / L]; };

// ---------- beam properties ----------
function beamProps(sec, mat) {
  const sp = Eg.sectionProps(sec), E = sec.E || mat.E, rho = sec.rho || mat.rho, nu = sec.nu || mat.nu;
  return { E, G: E / (2 * (1 + nu)), rho, A: sp.A, Iy: sp.I, Iz: sp.Iz, J: sec.J || sp.J, Ip: sp.Ip, kappa: Eg.shearCoef(sec.shape, nu) };
}

// ---------- 3D frame element, DOF order [ux uy uz rx ry rz] at each end ----------
// local x along the beam, local y = vref × x, local z = x × y (vref defaults to global z, as stiffness_global.m).
// Bending of w (local z) uses Iy = b h³/12, bending of v (local y) uses Iz = h b³/12.
function frame3D(p0, p1, pr, tim, vref) {
  const d = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]], L = norm3(d);
  if (!(L > 1e-12)) throw new Error('A beam has zero length.');
  const ex = [d[0] / L, d[1] / L, d[2] / L];
  let up = vref ? unit(vref) : [0, 0, 1];
  if (Math.abs(dot(ex, up)) > 0.99) { up = [1, 0, 0]; if (Math.abs(dot(ex, up)) > 0.99) up = [0, 1, 0]; }
  const ey = unit(cross(up, ex)), ez = cross(ex, ey);
  const Kl = zeros(12, 12), Ml = zeros(12, 12);
  const put = (A, idx, B, s) => { for (let a = 0; a < idx.length; a++) for (let b = 0; b < idx.length; b++) A[idx[a]][idx[b]] += (s ? s[a] * s[b] : 1) * B[a][b]; };
  const ka = pr.E * pr.A / L, ma = pr.rho * pr.A * L / 6;
  put(Kl, [0, 6], [[ka, -ka], [-ka, ka]]); put(Ml, [0, 6], [[2 * ma, ma], [ma, 2 * ma]]);
  const kt = pr.G * pr.J / L, mt = pr.rho * pr.Ip * L / 6;
  put(Kl, [3, 9], [[kt, -kt], [-kt, kt]]); put(Ml, [3, 9], [[2 * mt, mt], [mt, 2 * mt]]);
  const kGA = pr.kappa * pr.G * pr.A;
  const bv = Eg.bendBlock(L, pr.E * pr.Iz, kGA, pr.rho * pr.A, pr.rho * pr.Iz, tim);
  put(Kl, [1, 5, 7, 11], bv.K); put(Ml, [1, 5, 7, 11], bv.M);
  const bw = Eg.bendBlock(L, pr.E * pr.Iy, kGA, pr.rho * pr.A, pr.rho * pr.Iy, tim), sg = [1, -1, 1, -1];
  put(Kl, [2, 4, 8, 10], bw.K, sg); put(Ml, [2, 4, 8, 10], bw.M, sg);
  const R = [[ex[0], ey[0], ez[0]], [ex[1], ey[1], ez[1]], [ex[2], ey[2], ez[2]]];
  const rot = A => { // T A T^T, T = blkdiag(R, R, R, R)
    const B = zeros(12, 12), C = zeros(12, 12);
    for (let i = 0; i < 12; i++) for (let jb = 0; jb < 4; jb++) for (let a = 0; a < 3; a++) { let v = 0; for (let b = 0; b < 3; b++) v += A[i][3 * jb + b] * R[a][b]; B[i][3 * jb + a] = v; }
    for (let ib = 0; ib < 4; ib++) for (let a = 0; a < 3; a++) for (let j = 0; j < 12; j++) { let v = 0; for (let b = 0; b < 3; b++) v += R[a][b] * B[3 * ib + b][j]; C[3 * ib + a][j] = v; }
    for (let i = 0; i < 12; i++) for (let j = 0; j < i; j++) { const k = 0.5 * (C[i][j] + C[j][i]); C[i][j] = C[j][i] = k; }
    return C;
  };
  return { K: rot(Kl), M: rot(Ml), L, ex, ey, ez };
}

// ---------- small dense linear algebra ----------
function cholesky(A, n) { // in place, lower; false if not positive definite
  let dmax = 0; for (let i = 0; i < n; i++) dmax = Math.max(dmax, Math.abs(A[i][i]));
  for (let j = 0; j < n; j++) {
    const Aj = A[j]; let s = Aj[j]; for (let k = 0; k < j; k++) s -= Aj[k] * Aj[k];
    if (!(s > 1e-13 * dmax)) return false;
    const dj = Math.sqrt(s); Aj[j] = dj;
    for (let i = j + 1; i < n; i++) { const Ai = A[i]; let t = Ai[j]; for (let k = 0; k < j; k++) t -= Ai[k] * Aj[k]; Ai[j] = t / dj; }
  }
  return true;
}
function cholSolve(L, n, b) { // L L^T x = b, b overwritten
  for (let i = 0; i < n; i++) { let t = b[i]; for (let k = 0; k < i; k++) t -= L[i][k] * b[k]; b[i] = t / L[i][i]; }
  for (let i = n - 1; i >= 0; i--) { let t = b[i]; for (let k = i + 1; k < n; k++) t -= L[k][i] * b[k]; b[i] = t / L[i][i]; }
  return b;
}
function inv3(A) {
  const [a, b, c] = A[0], [d, e, f] = A[1], [g, h, i] = A[2];
  const det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
  if (!(Math.abs(det) > 0)) return null;
  return [[(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det], [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det], [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det]];
}
// dense real generalized symmetric eigenproblem K x = lam M x, all pairs, ascending, M-normalised vectors
function symGenEig(K, M, n) {
  const L = M.map(r => Float64Array.from(r));
  if (!cholesky(L, n)) throw new Error('A resonator mass matrix is not positive definite.');
  const Y = zeros(n, n);
  for (let i = 0; i < n; i++) { const Yi = Y[i]; Yi.set(K[i]); for (let k = 0; k < i; k++) { const l = L[i][k]; if (l) { const Yk = Y[k]; for (let c = 0; c < n; c++) Yi[c] -= l * Yk[c]; } } const inv = 1 / L[i][i]; for (let c = 0; c < n; c++) Yi[c] *= inv; }
  const C = zeros(n, n);
  for (let i = 0; i < n; i++) { const Ci = C[i]; for (let c = 0; c < n; c++) Ci[c] = Y[c][i]; for (let k = 0; k < i; k++) { const l = L[i][k]; if (l) { const Ck = C[k]; for (let c = 0; c < n; c++) Ci[c] -= l * Ck[c]; } } const inv = 1 / L[i][i]; for (let c = 0; c < n; c++) Ci[c] *= inv; }
  for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) { const v = 0.5 * (C[i][j] + C[j][i]); C[i][j] = C[j][i] = v; }
  const d = new Float64Array(n), e = new Float64Array(n);
  Eg.tred2(C, n, d, e, true); Eg.tqli(d, e, n, C);
  const ord = Array.from({ length: n }, (_, i) => i).sort((p, q) => d[p] - d[q]);
  const lam = ord.map(i => d[i]), vec = ord.map(j => { const y = new Float64Array(n); for (let i = 0; i < n; i++) y[i] = C[i][j]; for (let i = n - 1; i >= 0; i--) { let t = y[i]; for (let k = i + 1; k < n; k++) t -= L[k][i] * y[k]; y[i] = t / L[i][i]; } return y; });
  return { lam, vec };
}

// ---------- cell description ----------
// cell: { kind, a1:[x,y], a2:[x,y], nodes:[[x,y,z]], links:[{i,j,n:[n1,n2],sec}], sections:{name:{shape,b,h,E?,rho?,J?}},
//         resonators:[...], Henv, label }  (SI units)
function linkGeom(cell, l) {
  const a1 = cell.a1, a2 = cell.a2, q = cell.nodes[l.j], p0 = cell.nodes[l.i];
  return { p0, p1: [q[0] + l.n[0] * a1[0] + l.n[1] * a2[0], q[1] + l.n[0] * a1[1] + l.n[1] * a2[1], q[2]] };
}
function cellArea(cell) { return Math.abs(cell.a1[0] * cell.a2[1] - cell.a1[1] * cell.a2[0]); }

// ---------- static homogenisation ----------
function homogenize(cell, mat, opt = {}) {
  const tim = opt.beam === 'tim', nm = cell.nodes.length, nd = 6 * nm, S = cellArea(cell);
  if (!(S > 0)) throw new Error('The lattice vectors a1 and a2 are parallel.');
  if (!cell.links.length) throw new Error('The cell has no beams.');
  const els = cell.links.map((l, k) => {
    const sec = cell.sections[l.sec]; if (!sec) throw new Error(`Beam ${k} uses section "${l.sec}", which is not defined.`);
    if (l.i < 0 || l.j < 0 || l.i >= nm || l.j >= nm) throw new Error(`Beam ${k} refers to a node that does not exist.`);
    const g = linkGeom(cell, l), pr = beamProps(sec, mat), L = norm3([g.p1[0] - g.p0[0], g.p1[1] - g.p0[1], g.p1[2] - g.p0[2]]);
    return { l, p0: g.p0, p1: g.p1, pr, L, m: pr.rho * pr.A * L };
  });
  const deg = new Int32Array(nm); for (const e of els) { deg[e.l.i]++; deg[e.l.j]++; }
  for (let i = 0; i < nm; i++) if (!deg[i]) throw new Error(`Node ${i} is not connected to any beam.`);
  const mCell = els.reduce((s, e) => s + e.m, 0), zref = els.reduce((s, e) => s + e.m * (e.p0[2] + e.p1[2]) / 2, 0) / mCell;
  const x0 = cell.nodes.reduce((s, p) => s + p[0], 0) / nm, y0 = cell.nodes.reduce((s, p) => s + p[1], 0) / nm;
  const Fat = p => { const x = p[0] - x0, y = p[1] - y0, z = p[2] - zref;
    return [[x, 0, y / 2, -z * x, 0, -z * y / 2], [0, y, x / 2, 0, -z * y, -z * x / 2], [0, 0, 0, x * x / 2, y * y / 2, x * y / 2], [0, 0, 0, 0, y, x / 2], [0, 0, 0, -x, 0, -y / 2], [0, 0, 0, 0, 0, 0]]; };
  const Kpp = zeros(nd, nd), Kpf = zeros(nd, 6), Kff = zeros(6, 6);
  for (const e of els) {
    const fe = frame3D(e.p0, e.p1, e.pr, tim), F0 = Fat(e.p0), F1 = Fat(e.p1), F = [...F0, ...F1];
    const idx = []; for (let a = 0; a < 6; a++) idx.push(6 * e.l.i + a); for (let a = 0; a < 6; a++) idx.push(6 * e.l.j + a);
    const KF = zeros(12, 6);
    for (let a = 0; a < 12; a++) for (let c = 0; c < 6; c++) { let v = 0; for (let b = 0; b < 12; b++) v += fe.K[a][b] * F[b][c]; KF[a][c] = v; }
    for (let a = 0; a < 12; a++) { const ra = Kpp[idx[a]]; for (let b = 0; b < 12; b++) ra[idx[b]] += fe.K[a][b]; for (let c = 0; c < 6; c++) Kpf[idx[a]][c] += KF[a][c]; }
    for (let r = 0; r < 6; r++) for (let c = 0; c < 6; c++) { let v = 0; for (let a = 0; a < 12; a++) v += F[a][r] * KF[a][c]; Kff[r][c] += v; }
  }
  // gauge: the three fluctuating translations of node 0 are fixed; its rotations stay free
  const free = []; for (let i = 3; i < nd; i++) free.push(i);
  const n = free.length, Kr = zeros(n, n);
  for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) Kr[a][b] = 0.5 * (Kpp[free[a]][free[b]] + Kpp[free[b]][free[a]]);
  if (!cholesky(Kr, n)) throw new Error('The periodic frame has a mechanism (zero-energy fluctuation) or is disconnected: its micro-stiffness is not positive definite.');
  let dmin = Infinity, dmax = 0; for (let i = 0; i < n; i++) { const v = Kr[i][i] ** 2; dmin = Math.min(dmin, v); dmax = Math.max(dmax, v); }
  const X = []; for (let c = 0; c < 6; c++) { const b = new Float64Array(n); for (let a = 0; a < n; a++) b[a] = Kpf[free[a]][c]; X.push(cholSolve(Kr, n, b)); }
  const Q = zeros(6, 6);
  for (let r = 0; r < 6; r++) for (let c = 0; c < 6; c++) { let v = Kff[r][c]; for (let a = 0; a < n; a++) v -= Kpf[free[a]][r] * X[c][a]; Q[r][c] = v / S; }
  for (let r = 0; r < 6; r++) for (let c = 0; c < r; c++) { const v = 0.5 * (Q[r][c] + Q[c][r]); Q[r][c] = Q[c][r] = v; }
  const A = [0, 1, 2].map(r => [0, 1, 2].map(c => Q[r][c])), B = [0, 1, 2].map(r => [0, 1, 2].map(c => Q[r][c + 3])), D = [0, 1, 2].map(r => [0, 1, 2].map(c => Q[r + 3][c + 3]));
  const Ai = inv3(A); if (!Ai) throw new Error('The membrane stiffness A of the equivalent plate is singular: the frame can stretch freely in its plane.');
  const AiB = [0, 1, 2].map(r => [0, 1, 2].map(c => Ai[r][0] * B[0][c] + Ai[r][1] * B[1][c] + Ai[r][2] * B[2][c]));
  const Dfree = [0, 1, 2].map(r => [0, 1, 2].map(c => D[r][c] - (B[0][r] * AiB[0][c] + B[1][r] * AiB[1][c] + B[2][r] * AiB[2][c])));
  for (let r = 0; r < 3; r++) for (let c = 0; c < r; c++) { const v = 0.5 * (Dfree[r][c] + Dfree[c][r]); Dfree[r][c] = Dfree[c][r] = v; }
  const Henv = cell.Henv || 1;
  return { A, B, D, Dfree, Q, mCell, mArea: mCell / S, S, zref, Henv, rhoStar: mCell / (S * Henv), cond: dmax / dmin, nElem: els.length, nNodes: nm };
}

// ---------- resonators: clamped-base modes of any beam/mass geometry ----------
// G: { nodes:[[x,y,z]], beams:[{p,q,pr,vref}], root:[node indices, clamped], masses:[{node,m,J}] }
function frameModes(G, opt = {}) {
  const tim = opt.beam === 'tim', ne = opt.ne || 8, nodes = G.nodes.map(p => p.slice()), els = [];
  for (const b of G.beams) {
    const p0 = nodes[b.p], p1 = nodes[b.q]; let prev = b.p;
    for (let t = 1; t <= ne; t++) {
      let cur; if (t === ne) cur = b.q; else { nodes.push([p0[0] + (p1[0] - p0[0]) * t / ne, p0[1] + (p1[1] - p0[1]) * t / ne, p0[2] + (p1[2] - p0[2]) * t / ne]); cur = nodes.length - 1; }
      els.push({ p: prev, q: cur, pr: b.pr, vref: b.vref }); prev = cur;
    }
  }
  const N = nodes.length, nd = 6 * N, K = zeros(nd, nd), M = zeros(nd, nd);
  for (const e of els) {
    const fe = frame3D(nodes[e.p], nodes[e.q], e.pr, tim, e.vref), idx = [];
    for (let a = 0; a < 6; a++) idx.push(6 * e.p + a); for (let a = 0; a < 6; a++) idx.push(6 * e.q + a);
    for (let a = 0; a < 12; a++) for (let b = 0; b < 12; b++) { K[idx[a]][idx[b]] += fe.K[a][b]; M[idx[a]][idx[b]] += fe.M[a][b]; }
  }
  for (const pm of G.masses || []) { for (let a = 0; a < 3; a++) M[6 * pm.node + a][6 * pm.node + a] += pm.m; if (pm.J > 0) for (let a = 3; a < 6; a++) M[6 * pm.node + a][6 * pm.node + a] += pm.J; }
  const clamped = new Set(); for (const r of G.root) for (let a = 0; a < 6; a++) clamped.add(6 * r + a);
  const free = []; for (let i = 0; i < nd; i++) if (!clamped.has(i)) free.push(i);
  const nf = free.length; if (!nf) throw new Error('A resonator has no free part.');
  const Kf = zeros(nf, nf), Mf = zeros(nf, nf);
  for (let a = 0; a < nf; a++) for (let b = 0; b < nf; b++) { Kf[a][b] = K[free[a]][free[b]]; Mf[a][b] = M[free[a]][free[b]]; }
  // rigid translation of the whole resonator (root included) along x, y, z
  const MR = [0, 1, 2].map(c => { const v = new Float64Array(nd); for (let i = 0; i < nd; i++) { let s = 0; for (let n = 0; n < N; n++) s += M[i][6 * n + c]; v[i] = s; } return v; });
  let mTot = 0; for (let n = 0; n < N; n++) mTot += MR[2][6 * n + 2];
  const { lam, vec } = symGenEig(Kf, Mf, nf), nmo = Math.min(opt.nModes || 6, nf), modes = [];
  for (let j = 0; j < nmo; j++) {
    const phi = vec[j], Gm = [0, 1, 2].map(c => { let s = 0; for (let a = 0; a < nf; a++) s += phi[a] * MR[c][free[a]]; return s; });
    const m = Gm[0] ** 2 + Gm[1] ** 2 + Gm[2] ** 2, w2 = Math.max(0, lam[j]);
    modes.push({ f: Math.sqrt(w2) / TAU, m, k: w2 * m, n: m > 0 ? Gm.map(v => v / Math.sqrt(m)) : [0, 0, 1], share: m / mTot });
  }
  return { modes, mTot, nodes, els };
}
// Resonator record used by the PWE: { att:[[x,y]], pos:[x,y,z], mTot, modes:[{f,m,k,n}], draw:{beams:[{p0,p1,b,h,vref}], masses:[{p,side}]}, dir }
function segmentResonator(r, mat, opt) {
  const at = r.at, dir = unit(r.dir || [0, 0, 1]);
  let ax = r.axis ? unit(r.axis) : null;
  if (!ax || Math.abs(dot(ax, dir)) > 0.99) { const t = Math.abs(dir[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0]; ax = unit(cross(dir, cross(t, dir))); if (Math.abs(dir[2]) > 0.9) ax = unit(cross(dir, [0, 1, 0])); }
  ax = unit(cross(cross(dir, ax), dir)); // make the axis normal to the vibration direction
  const nodes = [at.slice()], beams = [], draw = { beams: [], masses: [] };
  let p = at.slice();
  r.segments.forEach((s, k) => {
    if (!(s.L > 0 && s.b > 0 && s.h > 0)) throw new Error(`Resonator segment ${k + 1} needs L, b and h greater than zero.`);
    const q = [p[0] + ax[0] * s.L, p[1] + ax[1] * s.L, p[2] + ax[2] * s.L]; nodes.push(q);
    const sec = { shape: 'rect', b: s.b, h: s.h, E: s.E, rho: s.rho };
    beams.push({ p: nodes.length - 2, q: nodes.length - 1, pr: beamProps(sec, mat), vref: dir });
    draw.beams.push({ p0: p, p1: q, b: s.b, h: s.h, vref: dir }); p = q;
  });
  const masses = [];
  if (r.tipMass > 0) { masses.push({ node: nodes.length - 1, m: r.tipMass, J: r.tipJ || 0 }); draw.masses.push({ p, side: r.tipSide || Math.cbrt(r.tipMass / mat.rho) }); }
  const fm = frameModes({ nodes, beams, root: [0], masses }, opt);
  return { att: r.att || [[at[0], at[1]]], pos: at, dir, axis: ax, mTot: fm.mTot, modes: fm.modes, draw, label: r.label, foot: r.foot, kind: 'segments' };
}
function makeResonator(r, mat, opt) {
  if (!Array.isArray(r.at) || r.at.length < 2) throw new Error('Each resonator needs "at": [x, y, z].');
  if (r.at.length === 2) r.at = [r.at[0], r.at[1], 0];
  if (r.segments && r.segments.length) return segmentResonator(r, mat, opt);
  const dir = unit(r.dir || [0, 0, 1]), draw = { beams: [], masses: [{ p: r.at, side: r.side || 0.004, spring: dir }] };
  if (r.modes && r.modes.length) {
    const modes = r.modes.map((q, k) => { if (!(q.f > 0 && q.m > 0)) throw new Error(`Resonator mode ${k + 1} needs f > 0 and m > 0.`); const n = unit(q.n || dir); return { f: q.f, m: q.m, k: (TAU * q.f) ** 2 * q.m, n, share: q.m / (r.mTot || q.m) }; });
    const mTot = r.mTot || modes.reduce((s, q) => Math.max(s, q.m), 0);
    return { att: r.att || [[r.at[0], r.at[1]]], pos: r.at, dir, mTot, modes, draw, label: r.label, foot: r.foot, kind: 'modes' };
  }
  if (!(r.m > 0)) throw new Error('A mass-spring resonator needs a mass greater than zero.');
  const f = r.f > 0 ? r.f : (r.k > 0 ? Math.sqrt(r.k / r.m) / TAU : NaN);
  if (!(f > 0)) throw new Error('A mass-spring resonator needs f_Hz or k_N_m.');
  return { att: r.att || [[r.at[0], r.at[1]]], pos: r.at, dir, mTot: r.m, modes: [{ f, m: r.m, k: (TAU * f) ** 2 * r.m, n: dir, share: 1 }], draw, label: r.label, foot: r.foot, kind: 'spring' };
}
// Reduce a resonator to the modes the PWE uses.
// 'modal': the lowest modes, each with its own effective mass (|Γ|²) and direction; the whole mass moves rigidly with the host.
// 'first': the first mode along the resonator direction with the TOTAL resonator mass, k = ω1² m (lumped convention).
function pickModes(R, conv, nModes, minShare = 0.002) {
  if (conv === 'first') {
    const d = R.dir || [0, 0, 1]; const q = R.modes.find(md => (dot(md.n, d)) ** 2 > 0.5 && md.m > 0) || R.modes[0];
    return [{ f: q.f, m: R.mTot, k: (TAU * q.f) ** 2 * R.mTot, n: d, share: 1, base: q }];
  }
  return R.modes.filter(q => q.share >= minShare).slice(0, nModes);
}

// host frame + resonators from a workbench spec (flat cell, z = 0)
function cellFromSpec(spec, mat, opt = {}) {
  const hostIdx = new Map(), nodes = [];
  const useHost = i => { if (!hostIdx.has(i)) { hostIdx.set(i, nodes.length); nodes.push([spec.nodes[i].x, spec.nodes[i].y, 0]); } return hostIdx.get(i); };
  const links = [], hostMassNodes = [];
  spec.beams.forEach(b => { if (b.role === 'res') return; links.push({ i: useHost(b.i), j: useHost(b.j), n: b.n.slice(), sec: b.sec }); });
  spec.nodes.forEach((p, i) => { if (p.role !== 'res' && p.m > 0) hostMassNodes.push({ node: useHost(i), m: p.m }); });
  const sections = {}; for (const k in spec.sections) sections[k] = Object.assign({}, spec.sections[k]);
  const Henv = Math.max(...links.map(l => sections[l.sec].h));
  const cell = { kind: 'workbench', a1: spec.a1.slice(), a2: spec.a2.slice(), nodes, links, sections, Henv, std: spec.std, spec, extraMass: hostMassNodes.reduce((s, q) => s + q.m, 0), label: 'Current workbench cell' };
  // resonator components (connected 'res' beams), with absolute positions through cell shifts
  const resBeams = spec.beams.filter(b => b.role === 'res'), comps = [], seen = new Set();
  const adj = new Map(); const addAdj = (i, e) => { if (!adj.has(i)) adj.set(i, []); adj.get(i).push(e); };
  resBeams.forEach((b, k) => { addAdj(b.i, k); addAdj(b.j, k); });
  const isRes = i => spec.nodes[i].role === 'res';
  for (const start of spec.nodes.keys()) {
    if (!isRes(start) || seen.has(start) || !adj.has(start)) continue;
    const off = new Map([[start, [0, 0]]]), stack = [start], bset = new Set(); seen.add(start);
    while (stack.length) { const i = stack.pop(); for (const k of adj.get(i) || []) { bset.add(k); const b = resBeams[k], o = off.get(i); const other = b.i === i ? b.j : b.i, s = b.i === i ? [o[0] + b.n[0], o[1] + b.n[1]] : [o[0] - b.n[0], o[1] - b.n[1]]; if (isRes(other) && !off.has(other)) { off.set(other, s); seen.add(other); stack.push(other); } } }
    comps.push({ off, beams: [...bset] });
  }
  const res = [];
  comps.forEach((c, ci) => {
    const nodesL = [], key = new Map(), root = [], masses = [];
    const pos = (i, o) => { const R = Eg.shiftVec(spec, o); return [spec.nodes[i].x + R[0], spec.nodes[i].y + R[1], 0]; };
    const get = (i, o) => { const k = `${i}|${o[0]}|${o[1]}`; if (!key.has(k)) { key.set(k, nodesL.length); nodesL.push(pos(i, o)); if (!isRes(i)) root.push(nodesL.length - 1); else if (spec.nodes[i].m > 0) masses.push({ node: nodesL.length - 1, m: spec.nodes[i].m, J: spec.nodes[i].J || 0, side: spec.nodes[i].side }); } return key.get(k); };
    const beams = [], draw = { beams: [], masses: [] };
    for (const k of c.beams) {
      const b = resBeams[k]; let oi, oj;
      if (isRes(b.i)) { oi = c.off.get(b.i); oj = [oi[0] + b.n[0], oi[1] + b.n[1]]; } else { oj = c.off.get(b.j); oi = [oj[0] - b.n[0], oj[1] - b.n[1]]; }
      const p = get(b.i, oi), q = get(b.j, oj), sec = spec.sections[b.sec];
      beams.push({ p, q, pr: beamProps(sec, mat) }); draw.beams.push({ p0: nodesL[p], p1: nodesL[q], b: sec.b, h: sec.h });
    }
    masses.forEach(m => draw.masses.push({ p: nodesL[m.node], side: m.side || 0.003 }));
    if (!root.length) throw new Error('A resonator is not attached to the host frame.');
    let clamp = root, hostIn = 0;
    if (opt.base === 'wall') { // carry the host beams around the root up to the nearest junctions (degree >= 3) into the resonator model
      const hdeg = new Map(); spec.beams.forEach(b => { if (b.role === 'res') return; hdeg.set(b.i, (hdeg.get(b.i) || 0) + 1); hdeg.set(b.j, (hdeg.get(b.j) || 0) + 1); });
      const keyOf = new Map(); key.forEach((v, k) => keyOf.set(v, k));
      const stack = root.map(r => { const [i, o1, o2] = keyOf.get(r).split('|').map(Number); return [i, [o1, o2]]; }), done = new Set(), used = new Set(); clamp = [];
      while (stack.length) {
        const [i, o] = stack.pop(), kk = `${i}|${o[0]}|${o[1]}`; if (done.has(kk)) continue; done.add(kk);
        const carries = spec.beams.some(b => (b.i === i || b.j === i) && (b.role === 'res' || b.sec === 'stub'));
        if (((hdeg.get(i) || 0) >= 3 && !carries) || spec.nodes[i].role === 'res') { clamp.push(get(i, o)); continue; }
        spec.beams.forEach((b, bi) => { if (b.role === 'res') return; let oj, other;
          if (b.i === i) { other = b.j; oj = [o[0] + b.n[0], o[1] + b.n[1]]; } else if (b.j === i) { other = b.i; oj = [o[0] - b.n[0], o[1] - b.n[1]]; } else return;
          const ek = [bi, ...(b.i === i ? o : oj)].join('|'); if (used.has(ek)) return; used.add(ek);
          const p = get(i, o), q = get(other, oj), sec = spec.sections[b.sec], pr = beamProps(sec, mat);
          beams.push({ p, q, pr }); draw.beams.push({ p0: nodesL[p], p1: nodesL[q], b: sec.b, h: sec.h, host: true }); hostIn += pr.rho * pr.A * Math.hypot(nodesL[q][0] - nodesL[p][0], nodesL[q][1] - nodesL[p][1]);
          stack.push([other, oj]); });
      }
      clamp = [...new Set(clamp)];
    }
    const fm = frameModes({ nodes: nodesL, beams, root: clamp, masses }, opt);
    const att = clamp.map(i => [nodesL[i][0], nodesL[i][1]]);
    res.push({ att, pos: nodesL[root[0]], dir: [0, 0, 1], mTot: fm.mTot, modes: fm.modes, draw, label: `Resonator ${ci + 1}`, kind: 'frame', hostIn });
  });
  for (const sp of spec.springs) {
    const p = spec.nodes[sp.node], f = Math.sqrt(sp.k / sp.m) / TAU, at = [p.x, p.y, 0];
    res.push({ att: [[p.x, p.y]], pos: at, dir: [0, 0, 1], mTot: sp.m, modes: [[0, 0, 1], [1, 0, 0], [0, 1, 0]].map(n => ({ f, m: sp.m, k: sp.k, n, share: 1 })), draw: { beams: [], masses: [{ p: at, side: 0.003, spring: [0, 0, 1] }] }, label: 'Mass on spring', kind: 'spring' });
  }
  cell.res = res; cell.extraMass -= res.reduce((s, r) => s + (r.hostIn || 0), 0);
  return cell;
}
// generic JSON cell (mm, g, GPa)
function cellFromJSON(txt, mat, opt) {
  let o; try { o = JSON.parse(txt); } catch (e) { throw new Error('The generic cell is not valid JSON: ' + e.message); }
  const mm = v => v / 1000, v2 = (v, name) => { if (!Array.isArray(v) || v.length !== 2 || !v.every(Number.isFinite)) throw new Error(`"${name}" must be [x, y] in mm.`); return [mm(v[0]), mm(v[1])]; };
  const a1 = v2(o.a1, 'a1'), a2 = v2(o.a2, 'a2');
  if (!o.sections) throw new Error('Add "sections", for example {"wall": {"shape": "rect", "b": 1, "h": 5}}.');
  const sections = {};
  for (const k in o.sections) { const s = o.sections[k];
    if (s.shape === 'circ') { if (!(s.d > 0)) throw new Error(`Section "${k}" needs d > 0.`); sections[k] = { shape: 'circ', b: mm(s.d), h: mm(s.d) }; }
    else { if (!(s.b > 0 && s.h > 0)) throw new Error(`Section "${k}" needs b > 0 and h > 0 (mm).`); sections[k] = { shape: 'rect', b: mm(s.b), h: mm(s.h) }; }
    if (s.E_GPa > 0) sections[k].E = s.E_GPa * 1e9; if (s.rho > 0) sections[k].rho = s.rho; if (s.J_mm4 > 0) sections[k].J = s.J_mm4 * 1e-12; }
  if (!Array.isArray(o.nodes) || !o.nodes.length) throw new Error('Add "nodes": [[x, y, z], ...] in mm.');
  const nodes = o.nodes.map((p, i) => { if (!Array.isArray(p) || p.length < 2 || !p.every(Number.isFinite)) throw new Error(`Node ${i} must be [x, y, z] in mm.`); return [mm(p[0]), mm(p[1]), mm(p[2] || 0)]; });
  if (!Array.isArray(o.beams) || !o.beams.length) throw new Error('Add "beams": [{"i": 0, "j": 1, "shift": [0, 0], "section": "wall"}, ...].');
  const links = o.beams.map((b, k) => { const n = b.shift || [0, 0];
    if (!(Number.isInteger(b.i) && Number.isInteger(b.j) && b.i >= 0 && b.j >= 0 && b.i < nodes.length && b.j < nodes.length)) throw new Error(`Beam ${k} refers to a node that does not exist.`);
    if (!(Array.isArray(n) && n.length === 2 && n.every(Number.isInteger))) throw new Error(`Beam ${k}: "shift" must be two integers.`);
    const sec = b.section || Object.keys(sections)[0]; if (!sections[sec]) throw new Error(`Beam ${k} uses section "${sec}", which is not defined.`);
    if (b.i === b.j && !n[0] && !n[1]) throw new Error(`Beam ${k} joins node ${b.i} to itself.`);
    return { i: b.i, j: b.j, n, sec }; });
  const zs = nodes.map(p => p[2]), Henv = o.thickness_mm > 0 ? mm(o.thickness_mm) : (Math.max(...zs) - Math.min(...zs)) + Math.max(...Object.values(sections).map(s => s.h));
  const cell = { kind: 'generic', a1, a2, nodes, links, sections, Henv, label: o.name || 'Generic cell' };
  cell.resonators = (o.resonators || []).map((r, k) => {
    const q = { at: Array.isArray(r.at) ? r.at.map(mm) : null, dir: r.dir, axis: r.axis, label: r.label || `Resonator ${k + 1}`, foot: r.footprint_mm > 0 ? mm(r.footprint_mm) : undefined };
    if (!q.at) throw new Error(`Resonator ${k + 1} needs "at": [x, y, z] in mm.`);
    if (r.attach) q.att = r.attach.map(p => [mm(p[0]), mm(p[1])]);
    if (r.segments) q.segments = r.segments.map(s => ({ L: mm(s.L), b: mm(s.b), h: mm(s.h), E: s.E_GPa > 0 ? s.E_GPa * 1e9 : undefined, rho: s.rho > 0 ? s.rho : undefined }));
    if (r.tip_mass_g > 0) { q.tipMass = r.tip_mass_g / 1000; q.tipSide = r.tip_side_mm > 0 ? mm(r.tip_side_mm) : undefined; }
    if (r.mass_g > 0) q.m = r.mass_g / 1000; if (r.f_Hz > 0) q.f = r.f_Hz; if (r.k_N_m > 0) q.k = r.k_N_m;
    if (r.modes) { q.modes = r.modes.map(md => ({ f: md.f_Hz, m: md.m_g / 1000, n: md.dir })); q.mTot = r.mass_g > 0 ? r.mass_g / 1000 : undefined; }
    return q;
  });
  return cell;
}
function cellToJSON(cell, res) {
  const r6 = v => +(+v).toPrecision(6), mm = v => r6(v * 1000);
  const o = { name: cell.label, units: 'mm, g, Hz, GPa, kg/m3', a1: cell.a1.map(mm), a2: cell.a2.map(mm), thickness_mm: mm(cell.Henv), sections: {} };
  for (const k in cell.sections) { const s = cell.sections[k]; const q = s.shape === 'circ' ? { shape: 'circ', d: mm(s.h) } : { shape: 'rect', b: mm(s.b), h: mm(s.h) }; if (s.E) q.E_GPa = r6(s.E / 1e9); if (s.rho) q.rho = r6(s.rho); if (s.J) q.J_mm4 = r6(s.J * 1e12); o.sections[k] = q; }
  o.nodes = cell.nodes.map(p => p.map(mm));
  o.beams = cell.links.map(l => ({ i: l.i, j: l.j, shift: l.n, section: l.sec }));
  o.resonators = (res || []).map(R => {
    const q = { label: R.label, at: R.pos.map(mm) };
    if (R.att.length > 1 || Math.hypot(R.att[0][0] - R.pos[0], R.att[0][1] - R.pos[1]) > 1e-9) q.attach = R.att.map(p => p.map(mm));
    if (R.src && R.src.segments) { q.dir = R.dir.map(r6); q.axis = R.axis.map(r6); q.segments = R.src.segments.map(s => { const z = { L: mm(s.L), b: mm(s.b), h: mm(s.h) }; if (s.E) z.E_GPa = r6(s.E / 1e9); if (s.rho) z.rho = r6(s.rho); return z; }); if (R.src.tipMass) q.tip_mass_g = r6(R.src.tipMass * 1000); }
    else { q.mass_g = r6(R.mTot * 1000); q.modes = R.modes.map(md => ({ f_Hz: r6(md.f), m_g: r6(md.m * 1000), dir: md.n.map(r6) })); }
    if (R.foot) q.footprint_mm = mm(R.foot);
    return q;
  });
  const line = x => JSON.stringify(x);
  return '{\n' + Object.keys(o).map(k => {
    const v = o[k];
    if (k === 'sections') return ` "sections": {\n${Object.keys(v).map(s => `  "${s}": ${line(v[s])}`).join(',\n')}\n }`;
    if (Array.isArray(v) && v.length && typeof v[0] === 'object') return ` "${k}": [\n${v.map(e => '  ' + line(e)).join(',\n')}\n ]`;
    return ` "${k}": ${line(v)}`;
  }).join(',\n') + '\n}';
}
function buildResonators(cell, mat, opt) {
  if (cell.res) return cell.res;
  return (cell.resonators || []).map(r => { const R = makeResonator(Object.assign({}, r), mat, opt); R.src = r; return R; });
}

// ---------- zone path ----------
function latticeFor(cell, mode) {
  const rc = Eg.recip(cell.a1, cell.a2), L1 = Math.hypot(...cell.a1), L2 = Math.hypot(...cell.a2);
  const cosA = (cell.a1[0] * cell.a2[0] + cell.a1[1] * cell.a2[1]) / (L1 * L2), eqL = Math.abs(L1 - L2) < 1e-6 * L1;
  let std = cell.std || null;
  if (!std && eqL && Math.abs(cosA) < 1e-6) std = 'square';
  if (!std && eqL && Math.abs(Math.abs(cosA) - 0.5) < 1e-6) std = 'hex';
  const Lat = { a1: cell.a1, a2: cell.a2, b1: rc.b1, b2: rc.b2, area: rc.area, std };
  Lat.bz = Eg.bzPolygon(rc.b1, rc.b2);
  const comb = (p, q) => [p * rc.b1[0] + q * rc.b2[0], p * rc.b1[1] + q * rc.b2[1]];
  if (!std && Math.abs(cosA) < 1e-6) { // rectangular
    Lat.path = mode === 'gxmg' ? { labels: ['Γ', 'X', 'M', 'Γ'], pts: [[0, 0], comb(0.5, 0), comb(0.5, 0.5), [0, 0]], mode } : { labels: ['Γ', 'X', 'M', 'Y', 'Γ'], pts: [[0, 0], comb(0.5, 0), comb(0.5, 0.5), comb(0, 0.5), [0, 0]], mode: 'full' };
  } else {
    Lat.sym = cell.spec ? Eg.symmetryInfo(cell.spec) : { mx: false, my: false, r90: false, r60: false };
    if (std && !cell.spec) Lat.sym = { mx: true, my: true, r90: std === 'square', r60: std === 'hex' };
    Lat.path = Eg.makePath(Lat, std ? 'std' : 'half');
  }
  return Lat;
}

// ---------- Bessel J1 (for the finite attachment footprint) ----------
function besselJ1(x) {
  const ax = Math.abs(x);
  if (ax < 8) { const y = x * x;
    const a = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
    const b = 144725228442.0 + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y)))); return a / b; }
  const z = 8 / ax, y = z * z, xx = ax - 2.356194491;
  const a = 1 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
  const b = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
  const v = Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * a - z * Math.sin(xx) * b); return x < 0 ? -v : v;
}
const form = x => (x < 1e-8 ? 1 : 2 * besselJ1(x) / x);

// ---------- plane-wave expansion ----------
// fam 'flex': Kirchhoff plate w with D_free, resonators coupled through n_z.
// fam 'ip'  : membrane (u, v) with A, resonators coupled through (n_x, n_y).
function pweSetup(Lat, plate, res, opt) {
  const NF = opt.NF, G = [];
  for (let p = -NF; p <= NF; p++) for (let q = -NF; q <= NF; q++) G.push([p * Lat.b1[0] + q * Lat.b2[0], p * Lat.b1[1] + q * Lat.b2[1]]);
  const modes = [];
  res.forEach((R, ri) => { for (const q of R.use) { const c = opt.fam === 'flex' ? [q.n[2]] : [q.n[0], q.n[1]]; if (Math.hypot(...c) < 1e-3) continue; modes.push({ ri, m: q.m, k: q.k, f: q.f, c }); } });
  return { G, modes, res, plate, fam: opt.fam, nw: (opt.fam === 'flex' ? 1 : 2) * G.length, r0: opt.r0 || 0, S: plate.S, mCell: plate.mCell + (opt.extraMass || 0) };
}
function pweMats(P, kx, ky, formAt) {
  const { G, modes, res, plate, fam, nw, S } = P, Ng = G.length, nm = modes.length, n = nw + nm, flex = fam === 'flex';
  const Kr = zeros(n, n), Mr = zeros(n, n), Mi = zeros(n, n);
  const Kx = new Float64Array(Ng), Ky = new Float64Array(Ng);
  for (let g = 0; g < Ng; g++) { Kx[g] = kx + G[g][0]; Ky[g] = ky + G[g][1]; }
  const D = plate.Dfree, A = plate.A;
  if (flex) for (let g = 0; g < Ng; g++) { const x = Kx[g], y = Ky[g];
    Kr[g][g] = S * (D[0][0] * x ** 4 + D[1][1] * y ** 4 + 2 * (D[0][1] + 2 * D[2][2]) * x * x * y * y + 4 * D[0][2] * x ** 3 * y + 4 * D[1][2] * x * y ** 3); }
  else for (let g = 0; g < Ng; g++) { const x = Kx[g], y = Ky[g], u = g, v = Ng + g;
    Kr[u][u] = S * (A[0][0] * x * x + 2 * A[0][2] * x * y + A[2][2] * y * y);
    Kr[v][v] = S * (A[2][2] * x * x + 2 * A[1][2] * x * y + A[1][1] * y * y);
    Kr[u][v] = Kr[v][u] = S * (A[0][2] * x * x + (A[0][1] + A[2][2]) * x * y + A[1][2] * y * y); }
  for (let i = 0; i < nw; i++) Mr[i][i] = P.mCell;
  // c_r(G) = footprint form factor × mean over the attachment points of exp(-i (k+G)·x_a)
  const cr = res.map(R => { const re = new Float64Array(Ng), im = new Float64Array(Ng), r0 = R.foot != null ? R.foot : P.r0;
    for (let g = 0; g < Ng; g++) { let a = 0, b = 0; for (const p of R.att) { const ph = Kx[g] * p[0] + Ky[g] * p[1]; a += Math.cos(ph); b -= Math.sin(ph); } const f = form((formAt ? Math.hypot(G[g][0], G[g][1]) : Math.hypot(Kx[g], Ky[g])) * r0) / R.att.length; re[g] = a * f; im[g] = b * f; }
    return { re, im }; });
  const blocks = flex ? 1 : 2;
  res.forEach((R, ri) => { const { re, im } = cr[ri], m = R.mTot; // host-block rigid mass: m conj(c_g) c_g'
    for (let bk = 0; bk < blocks; bk++) { const o = bk * Ng;
      for (let g = 0; g < Ng; g++) for (let h = 0; h < Ng; h++) { Mr[o + g][o + h] += m * (re[g] * re[h] + im[g] * im[h]); Mi[o + g][o + h] += m * (re[g] * im[h] - im[g] * re[h]); } } });
  modes.forEach((q, j) => { const J = nw + j, { re, im } = cr[q.ri];
    Kr[J][J] = q.k; Mr[J][J] = q.m;
    for (let bk = 0; bk < blocks; bk++) { const o = bk * Ng, s = q.m * q.c[bk]; if (!s) continue;
      for (let g = 0; g < Ng; g++) { Mr[o + g][J] += s * re[g]; Mi[o + g][J] -= s * im[g]; Mr[J][o + g] += s * re[g]; Mi[J][o + g] += s * im[g]; } } });
  return { Kr, Mr, Mi, n, nw, cr };
}
function pweAt(P, kx, ky, nb, wantFrac, wantVec) {
  const m = pweMats(P, kx, ky);
  return hermGenEig(m.Kr, null, m.Mr, m.Mi, m.n, nb, (wantFrac || wantVec) ? { from: m.nw, vec: wantVec } : null);
}
// Hermitian generalized eigenproblem via the real 2n embedding (as Eng.solveK)
function hermGenEig(Kr, Ki, Mr, Mi, n, nb, frac) {
  const N = 2 * n, A = zeros(N, N), B = zeros(N, N);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const k = Kr[i][j], ki = Ki ? Ki[i][j] : 0, m = Mr[i][j], mi = Mi ? Mi[i][j] : 0;
    A[i][j] = A[i + n][j + n] = k; A[i][j + n] = -ki; A[i + n][j] = ki;
    B[i][j] = B[i + n][j + n] = m; B[i][j + n] = -mi; B[i + n][j] = mi;
  }
  if (!cholesky(B, N)) throw new Error('The plane-wave mass matrix is not positive definite.');
  const Y = zeros(N, N);
  for (let i = 0; i < N; i++) { const Li = B[i], Yi = Y[i]; Yi.set(A[i]); for (let k = 0; k < i; k++) { const l = Li[k]; if (l) { const Yk = Y[k]; for (let c = 0; c < N; c++) Yi[c] -= l * Yk[c]; } } const inv = 1 / Li[i]; for (let c = 0; c < N; c++) Yi[c] *= inv; }
  const C = zeros(N, N);
  for (let i = 0; i < N; i++) { const Li = B[i], Ci = C[i]; for (let c = 0; c < N; c++) Ci[c] = Y[c][i]; for (let k = 0; k < i; k++) { const l = Li[k]; if (l) { const Ck = C[k]; for (let c = 0; c < N; c++) Ci[c] -= l * Ck[c]; } } const inv = 1 / Li[i]; for (let c = 0; c < N; c++) Ci[c] *= inv; }
  for (let i = 0; i < N; i++) for (let j = 0; j < i; j++) { const v = 0.5 * (C[i][j] + C[j][i]); C[i][j] = C[j][i] = v; }
  const d = new Float64Array(N), e = new Float64Array(N), hs = new Float64Array(N);
  Eg.tred2(C, N, d, e, false, hs);
  const d0 = d.slice(), e0 = e.slice();
  Eg.tqli(d, e, N, null);
  const ev = Array.from(d).sort((p, q) => p - q), take = Math.min(nb, n), out = { f: [], frac: [] };
  for (let b = 0; b < take; b++) {
    const lam = ev[2 * b]; out.f.push(Math.sqrt(Math.max(0, lam)) / TAU);
    if (!frac) continue;
    const z = Eg.tridiagInv(d0, e0, N, lam, Math.abs(ev[N - 1]) + Math.abs(ev[0]));
    for (let i = 2; i < N; i++) { const h = hs[i]; if (!h) continue; const u = C[i]; let s = 0; for (let k = 0; k < i; k++) s += u[k] * z[k]; s /= h; for (let k = 0; k < i; k++) z[k] -= s * u[k]; }
    const y = z; for (let i = N - 1; i >= 0; i--) { let t = y[i]; for (let k = i + 1; k < N; k++) t -= B[k][i] * y[k]; y[i] = t / B[i][i]; }
    let part = 0, tot = 0; for (let i = 0; i < n; i++) { const e = Mr[i][i] * (y[i] * y[i] + y[i + n] * y[i + n]); tot += e; if (frac.idx ? frac.idx.has(i) : i >= frac.from) part += e; }
    out.frac.push(tot > 0 ? Math.min(1, Math.max(0, part / tot)) : 0);
    if (frac.vec) { if (!out.vec) out.vec = []; out.vec.push([y.slice(0, n), y.slice(n, 2 * n)]); }
  }
  return out;
}
// ---------- full 3D Bloch FE of the frame (reference for the homogenised model) ----------
// Host beams meshed with ne elements, resonators attached at the nearest point of the frame as the same modal
// mass-spring sets used by the PWE. Returns a model for bloch3DAt.
function bloch3DModel(cell, res, mat, opt = {}) {
  const tim = opt.beam === 'tim', ne = Math.max(1, opt.ne || 3), a1 = cell.a1, a2 = cell.a2;
  const nodes = cell.nodes.map(p => p.slice()), els = [];
  const geo = cell.links.map(l => linkGeom(cell, l));
  // attachment sites: nearest point on any link (over neighbouring images)
  const cuts = cell.links.map(() => []), site = [];
  res.forEach((R, ri) => {
    const a = [R.att[0][0], R.att[0][1], R.pos[2] || 0]; let best = { d: Infinity };
    geo.forEach((g, li) => { for (let s1 = -1; s1 <= 1; s1++) for (let s2 = -1; s2 <= 1; s2++) {
      const sh = [s1 * a1[0] + s2 * a2[0], s1 * a1[1] + s2 * a2[1], 0], q0 = g.p0.map((v, c) => v + sh[c]), q1 = g.p1.map((v, c) => v + sh[c]);
      const d = q1.map((v, c) => v - q0[c]), L2 = dot(d, d); let t = L2 > 0 ? dot(a.map((v, c) => v - q0[c]), d) / L2 : 0; t = Math.min(1, Math.max(0, t));
      const dist = norm3(a.map((v, c) => v - q0[c] - t * d[c])); if (dist < best.d) best = { d: dist, li, t }; } });
    cuts[best.li].push({ t: best.t, ri }); site.push(best);
  });
  const resNode = new Array(res.length);
  cell.links.forEach((l, li) => {
    const g = geo[li], pr = beamProps(cell.sections[l.sec], mat);
    const ts = new Set(); for (let t = 1; t < ne; t++) ts.add(+(t / ne).toFixed(9));
    for (const c of cuts[li]) ts.add(+c.t.toFixed(9));
    const list = [...ts].filter(t => t > 1e-9 && t < 1 - 1e-9).sort((x, y) => x - y);
    const chain = [{ node: l.i, t: 0 }];
    for (const t of list) { nodes.push(g.p0.map((v, c) => v + (g.p1[c] - v) * t)); chain.push({ node: nodes.length - 1, t }); }
    chain.push({ node: l.j, t: 1, shift: l.n });
    for (const c of cuts[li]) { const hit = chain.reduce((b, x) => (Math.abs(x.t - c.t) < Math.abs(b.t - c.t) ? x : b), chain[0]); resNode[c.ri] = { node: hit.node, shift: hit.shift || [0, 0] }; }
    for (let s = 0; s < chain.length - 1; s++) {
      const P0 = s === 0 ? g.p0 : nodes[chain[s].node], P1 = s === chain.length - 2 ? g.p1 : nodes[chain[s + 1].node];
      const fe = frame3D(P0, P1, pr, tim);
      els.push({ p: chain[s].node, q: chain[s + 1].node, n: s === chain.length - 2 ? l.n : [0, 0], K: fe.K, M: fe.M, X: [P0, P1], draw: { b: cell.sections[l.sec].b, h: cell.sections[l.sec].h }, mass: pr.rho * pr.A * fe.L });
    }
  });
  const N = nodes.length, springs = [];
  res.forEach((R, ri) => { const rn = resNode[ri]; for (const q of R.use) springs.push({ node: rn.node, m: q.m, k: q.k, n: q.n, ri }); });
  const n = 6 * N + springs.length, zIdx = new Set(); for (let i = 0; i < N; i++) zIdx.add(6 * i + 2);
  return { N, n, nodes, els, springs, res, resNode, zIdx, site, ne };
}
function bloch3DAt(md, cell, kx, ky, nb, want) {
  const { n, N, els, springs, res, resNode } = md, Kr = zeros(n, n), Ki = zeros(n, n), Mr = zeros(n, n), Mi = zeros(n, n);
  for (const e of els) {
    const R = [e.n[0] * cell.a1[0] + e.n[1] * cell.a2[0], e.n[0] * cell.a1[1] + e.n[1] * cell.a2[1]], ph = kx * R[0] + ky * R[1], c = Math.cos(ph), s = Math.sin(ph);
    for (let a = 0; a < 12; a++) for (let b = 0; b < 12; b++) {
      const I = a < 6 ? 6 * e.p + a : 6 * e.q + a - 6, J = b < 6 ? 6 * e.p + b : 6 * e.q + b - 6, k = e.K[a][b], m = e.M[a][b];
      let pr = 1, pi = 0; if (a < 6 && b >= 6) { pr = c; pi = s; } else if (a >= 6 && b < 6) { pr = c; pi = -s; }
      Kr[I][J] += k * pr; Ki[I][J] += k * pi; Mr[I][J] += m * pr; Mi[I][J] += m * pi;
    }
  }
  res.forEach((Rz, ri) => { const nd = resNode[ri].node; for (let a = 0; a < 3; a++) Mr[6 * nd + a][6 * nd + a] += Rz.mTot; });
  springs.forEach((sp, j) => { const J = 6 * N + j; Kr[J][J] += sp.k; Mr[J][J] += sp.m;
    for (let a = 0; a < 3; a++) { const v = sp.m * sp.n[a]; if (!v) continue; Mr[6 * sp.node + a][J] += v; Mr[J][6 * sp.node + a] += v; } });
  const frac = want === 'z' ? { idx: md.zIdx } : (want === 'res' ? { from: 6 * N } : null);
  return hermGenEig(Kr, Ki, Mr, Mi, n, nb, frac);
}

// ---------- complex eigenvalues: Hessenberg by stabilised elimination, then single-shift QR with Givens rotations ----------
function cEig(Ar0, Ai0, n) {
  const ar = Ar0.map(r => Float64Array.from(r)), ai = Ai0.map(r => Float64Array.from(r));
  const ab = (i, j) => Math.abs(ar[i][j]) + Math.abs(ai[i][j]);
  for (let m = 1; m < n - 1; m++) {
    let piv = m, mx = ab(m, m - 1); for (let i = m + 1; i < n; i++) if (ab(i, m - 1) > mx) { mx = ab(i, m - 1); piv = i; }
    if (piv !== m) { [ar[piv], ar[m]] = [ar[m], ar[piv]]; [ai[piv], ai[m]] = [ai[m], ai[piv]]; for (let j = 0; j < n; j++) { let t = ar[j][piv]; ar[j][piv] = ar[j][m]; ar[j][m] = t; t = ai[j][piv]; ai[j][piv] = ai[j][m]; ai[j][m] = t; } }
    if (mx === 0) continue;
    const pr = ar[m][m - 1], pi = ai[m][m - 1], pd = pr * pr + pi * pi;
    for (let i = m + 1; i < n; i++) {
      const xr = ar[i][m - 1], xi = ai[i][m - 1]; if (!xr && !xi) continue;
      const yr = (xr * pr + xi * pi) / pd, yi = (xi * pr - xr * pi) / pd;
      for (let j = m - 1; j < n; j++) { const br = ar[m][j], bi = ai[m][j]; ar[i][j] -= yr * br - yi * bi; ai[i][j] -= yr * bi + yi * br; }
      for (let j = 0; j < n; j++) { const br = ar[j][i], bi = ai[j][i]; ar[j][m] += yr * br - yi * bi; ai[j][m] += yr * bi + yi * br; }
    }
  }
  for (let i = 2; i < n; i++) for (let j = 0; j < i - 1; j++) { ar[i][j] = 0; ai[i][j] = 0; }
  const out = []; let hi = n - 1, its = 0, tot = 0, nrm = 0;
  for (let i = 0; i < n; i++) for (let j = Math.max(0, i - 1); j < n; j++) nrm += ab(i, j);
  const cs = new Float64Array(n), snr = new Float64Array(n), sni = new Float64Array(n);
  while (hi >= 0) {
    if (hi === 0) { out.push([ar[0][0], ai[0][0]]); break; }
    let lo = hi; while (lo > 0) { const s0 = ab(lo, lo) + ab(lo - 1, lo - 1) || nrm; if (ab(lo, lo - 1) <= 1e-15 * s0) { ar[lo][lo - 1] = ai[lo][lo - 1] = 0; break; } lo--; }
    if (lo === hi) { out.push([ar[hi][hi], ai[hi][hi]]); hi--; its = 0; continue; }
    if (++tot > 60 * n) return null;
    its++;
    // Wilkinson shift from the trailing 2x2
    const a = [ar[hi - 1][hi - 1], ai[hi - 1][hi - 1]], b = [ar[hi - 1][hi], ai[hi - 1][hi]], c = [ar[hi][hi - 1], ai[hi][hi - 1]], d = [ar[hi][hi], ai[hi][hi]];
    let mu;
    if (its % 11 === 10) mu = [d[0] + 0.75 * Math.abs(ar[hi][hi - 1]), d[1]];
    else { const tr = [(a[0] + d[0]) / 2, (a[1] + d[1]) / 2], df = [(a[0] - d[0]) / 2, (a[1] - d[1]) / 2], bc = [b[0] * c[0] - b[1] * c[1], b[0] * c[1] + b[1] * c[0]];
      const z = [df[0] * df[0] - df[1] * df[1] + bc[0], 2 * df[0] * df[1] + bc[1]], r = Math.hypot(z[0], z[1]), sq = [Math.sqrt((r + z[0]) / 2), Math.sign(z[1] || 1) * Math.sqrt(Math.max(0, (r - z[0]) / 2))];
      const m1 = [tr[0] + sq[0], tr[1] + sq[1]], m2 = [tr[0] - sq[0], tr[1] - sq[1]];
      mu = Math.hypot(m1[0] - d[0], m1[1] - d[1]) < Math.hypot(m2[0] - d[0], m2[1] - d[1]) ? m1 : m2; }
    for (let i = lo; i <= hi; i++) { ar[i][i] -= mu[0]; ai[i][i] -= mu[1]; }
    for (let k = lo; k < hi; k++) {
      const xr = ar[k][k], xi = ai[k][k], yr = ar[k + 1][k], yi = ai[k + 1][k], ax = Math.hypot(xr, xi), r = Math.hypot(ax, Math.hypot(yr, yi));
      let c0, sr, si; if (r === 0) { c0 = 1; sr = 0; si = 0; } else if (ax === 0) { c0 = 0; sr = 1; si = 0; } else { c0 = ax / r; const ur = xr / ax, ui = xi / ax; sr = (ur * yr + ui * yi) / r; si = (ui * yr - ur * yi) / r; }
      cs[k] = c0; snr[k] = sr; sni[k] = si;
      for (let j = k; j <= hi; j++) { const pr = ar[k][j], pi = ai[k][j], qr = ar[k + 1][j], qi = ai[k + 1][j];
        ar[k][j] = c0 * pr + sr * qr - si * qi; ai[k][j] = c0 * pi + sr * qi + si * qr;
        ar[k + 1][j] = -(sr * pr + si * pi) + c0 * qr; ai[k + 1][j] = -(sr * pi - si * pr) + c0 * qi; }
    }
    for (let k = lo; k < hi; k++) { const c0 = cs[k], sr = snr[k], si = sni[k];
      for (let i = lo; i <= Math.min(k + 2, hi); i++) { const pr = ar[i][k], pi = ai[i][k], qr = ar[i][k + 1], qi = ai[i][k + 1];
        ar[i][k] = c0 * pr + sr * qr + si * qi; ai[i][k] = c0 * pi + sr * qi - si * qr;
        ar[i][k + 1] = -(sr * pr - si * pi) + c0 * qr; ai[i][k + 1] = -(sr * pi + si * pr) + c0 * qi; } }
    for (let i = lo; i <= hi; i++) { ar[i][i] += mu[0]; ai[i][i] += mu[1]; }
  }
  return out;
}
// complex band of the homogenised plate at frequency f along a1 (phase along a2 held at 0, the Γ–X leg):
// (sum_j s^j C_j) u = 0 with s the complex wavenumber along b1; the resonator phases are factored out exactly for a
// single attachment point (several points use their mean position); the footprint factor is taken at |G|.
function pweComplexAt(P, Lat, f) {
  const { G, fam, S } = P, Ng = G.length, flex = fam === 'flex', w2 = (TAU * f) ** 2;
  const bl = Math.hypot(...Lat.b1), e = [Lat.b1[0] / bl, Lat.b1[1] / bl], alpha = bl / 2, d = flex ? 4 : 2;
  const { modes, res } = P, nm = modes.length, nw = (flex ? 1 : 2) * Ng, n = nw + nm, blocks = flex ? 1 : 2;
  const Cr = [], Ci = []; for (let j = 0; j <= d; j++) { Cr.push(zeros(n, n)); Ci.push(zeros(n, n)); }
  // mass terms, expanded to second order in s for resonators attached at several points (exact for one point)
  const addM = (j, i, k, vr, vi) => { Cr[j][i][k] -= w2 * vr; Ci[j][i][k] -= w2 * vi; };
  for (let i = 0; i < nw; i++) addM(0, i, i, P.mCell, 0);
  modes.forEach((q, j) => { const J = nw + j; Cr[0][J][J] += q.k; addM(0, J, J, q.m, 0); });
  const fm = G.map(g => Math.hypot(g[0], g[1]));
  res.forEach((R, ri) => {
    const r0 = R.foot != null ? R.foot : P.r0, na = R.att.length, xb = [R.att.reduce((t, p) => t + p[0], 0) / na, R.att.reduce((t, p) => t + p[1], 0) / na];
    const dl = R.att.map(p => (p[0] - xb[0]) * e[0] + (p[1] - xb[1]) * e[1]), F = fm.map(x => form(x * r0));
    const ph = (g, p) => G[g][0] * p[0] + G[g][1] * p[1];
    // e^{i theta} ≈ 1 + i theta − theta²/2, theta = s·d (d = e·δ)
    const ser = dd => [[1, 0], [0, dd], [-dd * dd / 2, 0]];
    for (let bk = 0; bk < blocks; bk++) { const o = bk * Ng;
      for (let a = 0; a < na; a++) for (let b = 0; b < na; b++) { const sr = ser(alpha * (dl[a] - dl[b])), mm = R.mTot / (na * na);
        for (let g = 0; g < Ng; g++) for (let h = 0; h < Ng; h++) { const t = ph(g, R.att[a]) - ph(h, R.att[b]), cr = Math.cos(t) * F[g] * F[h] * mm, ci = Math.sin(t) * F[g] * F[h] * mm;
          for (let j = 0; j < 3; j++) { if (na === 1 && j) break; const [er, ei] = sr[j]; addM(j, o + g, o + h, cr * er - ci * ei, cr * ei + ci * er); } } } }
    modes.forEach((q, jm) => { if (q.ri !== ri) return; const J = nw + jm;
      for (let bk = 0; bk < blocks; bk++) { const o = bk * Ng, sc = q.m * q.c[bk] / na; if (!sc) continue;
        for (let a = 0; a < na; a++) { const sp = ser(alpha * dl[a]), sm = ser(-alpha * dl[a]);
          for (let g = 0; g < Ng; g++) { const t = ph(g, R.att[a]), cr = Math.cos(t) * F[g] * sc, ci = Math.sin(t) * F[g] * sc;
            for (let j = 0; j < 3; j++) { if (na === 1 && j) break; let [er, ei] = sp[j]; addM(j, o + g, J, cr * er - ci * ei, cr * ei + ci * er); [er, ei] = sm[j]; addM(j, J, o + g, cr * er + ci * ei, -ci * er + cr * ei); } } } } });
  });
  const pm = (a, b) => { const o = new Array(a.length + b.length - 1).fill(0); a.forEach((x, i) => b.forEach((y, j) => { o[i + j] += x * y; })); return o; };
  const add = (a, b, s = 1) => { const o = a.slice(); b.forEach((y, i) => { o[i] = (o[i] || 0) + s * y; }); return o; };
  const D = P.plate.Dfree, A = P.plate.A;
  for (let g = 0; g < Ng; g++) {
    const kx = [G[g][0], e[0] * alpha], ky = [G[g][1], e[1] * alpha];
    if (flex) { const x2 = pm(kx, kx), y2 = pm(ky, ky); let q = pm(pm(x2, x2), [D[0][0]]); q = add(q, pm(pm(y2, y2), [D[1][1]])); q = add(q, pm(pm(x2, y2), [2 * (D[0][1] + 2 * D[2][2])])); q = add(q, pm(pm(x2, pm(kx, ky)), [4 * D[0][2]])); q = add(q, pm(pm(y2, pm(kx, ky)), [4 * D[1][2]]));
      q.forEach((v, j) => { Cr[j][g][g] += S * v; }); }
    else { const xx = pm(kx, kx), yy = pm(ky, ky), xy = pm(kx, ky), u = g, v = Ng + g;
      const kuu = add(add(pm(xx, [A[0][0]]), pm(xy, [2 * A[0][2]])), pm(yy, [A[2][2]])), kvv = add(add(pm(xx, [A[2][2]]), pm(xy, [2 * A[1][2]])), pm(yy, [A[1][1]])), kuv = add(add(pm(xx, [A[0][2]]), pm(xy, [A[0][1] + A[2][2]])), pm(yy, [A[1][2]]));
      kuu.forEach((x, j) => { Cr[j][u][u] += S * x; }); kvv.forEach((x, j) => { Cr[j][v][v] += S * x; }); kuv.forEach((x, j) => { Cr[j][u][v] += S * x; Cr[j][v][u] += S * x; }); }
  }
  // Möbius map s = alpha (1 + t) / (1 - t): Q(t) = sum_j C_j (1 + t)^j (1 - t)^(d - j)
  const Qr = [], Qi = []; for (let i = 0; i <= d; i++) { Qr.push(zeros(n, n)); Qi.push(zeros(n, n)); }
  for (let j = 0; j <= d; j++) { let b = [1]; for (let q = 0; q < j; q++) b = pm(b, [1, 1]); for (let q = 0; q < d - j; q++) b = pm(b, [1, -1]);
    b.forEach((cf, i) => { if (!cf) return; for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) { Qr[i][r][c] += cf * Cr[j][r][c]; Qi[i][r][c] += cf * Ci[j][r][c]; } }); }
  // monic: Q_d^{-1} Q_i via a complex LU in the real 2n embedding
  const N2 = 2 * n, Lm = zeros(N2, N2), Rh = zeros(N2, 2 * n * d);
  for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) { Lm[r][c] = Lm[r + n][c + n] = Qr[d][r][c]; Lm[r][c + n] = -Qi[d][r][c]; Lm[r + n][c] = Qi[d][r][c]; }
  for (let i = 0; i < d; i++) for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) { Rh[r][i * n + c] = Qr[i][r][c]; Rh[r + n][i * n + c] = Qi[i][r][c]; Rh[r][d * n + i * n + c] = -Qi[i][r][c]; Rh[r + n][d * n + i * n + c] = Qr[i][r][c]; }
  if (!Eg.luSolve(Lm, Rh, N2, 2 * n * d)) return null;
  // companion matrix (complex), size d n
  const m = d * n, Tr = zeros(m, m), Ti = zeros(m, m);
  for (let i = 0; i < (d - 1) * n; i++) Tr[i][i + n] = 1;
  for (let i = 0; i < d; i++) for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) { Tr[(d - 1) * n + r][i * n + c] = -Rh[r][i * n + c]; Ti[(d - 1) * n + r][i * n + c] = -Rh[r + n][i * n + c]; }
  const ev = cEig(Tr, Ti, m); if (!ev) return null;
  const a1e = e[0] * Lat.a1[0] + e[1] * Lat.a1[1], out = [];
  for (const [tr, ti] of ev) { const den = (1 - tr) ** 2 + ti * ti; if (den < 1e-14) continue; if (Math.abs((1 - tr * tr - ti * ti) / den) > 1 + 1e-6) continue; // keep Re(s) in the first zone
    const sr = alpha * (1 - tr * tr - ti * ti) / den, si = alpha * 2 * ti / den, mr = sr * a1e, mi = si * a1e; // μ = s (e·a1)
    const ph = Math.abs(Math.atan2(Math.sin(mr), Math.cos(mr))); out.push({ att: Math.abs(mi), ph }); }
  return out;
}
// gaps between consecutive bands over the whole path; minRel = minimum relative width
function gapsOf(F, minRel = 0.002) {
  const nb = F[0].length, gaps = [];
  for (let b = 0; b < nb - 1; b++) {
    let hi = -Infinity, lo = Infinity; for (const r of F) { hi = Math.max(hi, r[b]); lo = Math.min(lo, r[b + 1]); }
    if (lo > hi && (lo - hi) > minRel * 0.5 * (lo + hi)) gaps.push({ lo: hi, hi: lo, below: b + 1 });
  }
  return gaps;
}
// plate wavelength at frequency f (bare homogenised plate), worst direction, for the validity check
function plateWavelength(plate, f, fam) {
  const w = TAU * f, m = plate.mArea; let kmax = 0;
  for (let t = 0; t < 36; t++) { const a = Math.PI * t / 36, c = Math.cos(a), s = Math.sin(a);
    if (fam === 'flex') { const D = plate.Dfree, Dk = D[0][0] * c ** 4 + D[1][1] * s ** 4 + 2 * (D[0][1] + 2 * D[2][2]) * c * c * s * s + 4 * D[0][2] * c ** 3 * s + 4 * D[1][2] * c * s ** 3; if (Dk > 0) kmax = Math.max(kmax, Math.pow(m * w * w / Dk, 0.25)); }
    else { const A = plate.A, k11 = A[0][0] * c * c + 2 * A[0][2] * c * s + A[2][2] * s * s, k22 = A[2][2] * c * c + 2 * A[1][2] * c * s + A[1][1] * s * s, k12 = A[0][2] * c * c + (A[0][1] + A[2][2]) * c * s + A[1][2] * s * s;
      const lmin = 0.5 * (k11 + k22) - Math.sqrt(0.25 * (k11 - k22) ** 2 + k12 * k12); if (lmin > 0) kmax = Math.max(kmax, w * Math.sqrt(m / lmin)); } }
  return kmax > 0 ? TAU / kmax : Infinity;
}

return { cEig, pweMats, pweComplexAt, bloch3DModel, bloch3DAt, frame3D, homogenize, frameModes, makeResonator, segmentResonator, pickModes, cellFromSpec, cellFromJSON, cellToJSON, buildResonators, latticeFor, pweSetup, pweAt, gapsOf, plateWavelength, besselJ1, beamProps, symGenEig, cellArea, linkGeom, TAU };
})();
if (typeof module !== 'undefined') module.exports = Hom;
