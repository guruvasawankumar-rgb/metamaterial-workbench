// Validation of the homogenisation + plane-wave module (src/homog.js).  Run:  node tests/validate_homog.js
// References: closed-form grillage and plate results and the workbench Bloch FE.
const Eng = require('../src/engine.js'), Hom = require('../src/homog.js');
const sc = v => (Array.isArray(v) ? (Array.isArray(v[0]) ? v[0][0] : v[0]) : v);
const TAU = 2 * Math.PI, deg = Math.PI / 180;
let fails = 0, n = 0;
function check(name, got, want, tol) {
  const err = Math.abs(got - want) / Math.max(Math.abs(want), 1e-30), ok = err <= tol; n++;
  if (!ok) fails++;
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${name}: ${(+got).toPrecision(7)} vs ${(+want).toPrecision(7)}  (rel. err ${err.toExponential(1)}, tol ${tol})`);
}
const t0 = Date.now();

console.log('— 1. Homogenisation, closed form (square grillage of beams, spacing a)');
const a = 0.02, al = { E: 70e9, nu: 0.33, rho: 2700 }, sq = { shape: 'rect', b: 0.002, h: 0.004 };
const grill = { a1: [a, 0], a2: [0, a], nodes: [[0, 0, 0]], links: [{ i: 0, j: 0, n: [1, 0], sec: 'w' }, { i: 0, j: 0, n: [0, 1], sec: 'w' }], sections: { w: sq }, Henv: 0.004 };
{
  const pl = Hom.homogenize(grill, al, { beam: 'eb' }), pr = Hom.beamProps(sq, al);
  check('D11 = EI/a', pl.Dfree[0][0], al.E * pr.Iy / a, 1e-9);
  check('D66 = GJ/(2a)  (twist κxy = 2 w,xy)', pl.Dfree[2][2], pr.G * pr.J / (2 * a), 1e-9);
  check('A11 = EA/a', pl.A[0][0], al.E * pr.A / a, 1e-9);
  check('mass per area = 2ρA/a', pl.mArea, 2 * al.rho * pr.A / a, 1e-12);
  const Lat = Hom.latticeFor(grill, 'std'), kx = 0.3 * Math.PI / a;
  check('PWE bare plate f(k) = √(D k⁴/m)/2π', Hom.pweAt(Hom.pweSetup(Lat, pl, [], { fam: 'flex', NF: 2 }), kx, 0, 1).f[0], Math.sqrt(pl.Dfree[0][0] * kx ** 4 / pl.mArea) / TAU, 1e-9);
  check('PWE bare membrane, shear wave √(A66/m) k/2π', Hom.pweAt(Hom.pweSetup(Lat, pl, [], { fam: 'ip', NF: 2 }), kx, 0, 1).f[0], kx * Math.sqrt(pl.A[2][2] / pl.mArea) / TAU, 1e-9);
  const m = 2e-3, f0 = 500, res = [{ att: [[a / 2, a / 2]], mTot: m, use: [{ m, k: (TAU * f0) ** 2 * m, f: f0, n: [0, 0, 1] }] }];
  check('Γ, rigid plate + resonator: f0·√(1 + m/M)', Hom.pweAt(Hom.pweSetup(Lat, pl, res, { fam: 'flex', NF: 0 }), 0, 0, 2).f[1], f0 * Math.sqrt(1 + m / pl.mCell), 1e-9);
}

console.log('— 2. Homogenisation vs full Bloch FE, long-wave flexural branch');
{ // your 14 mm hexagon, bare (2D frame, workbench FE)
  const G = { lattice: 'hexagon', l: 0.014 / Math.sqrt(3), h: 0.014 / Math.sqrt(3), theta: 30 };
  const base = Eng.baseCell(G); base.sections = { host: { shape: 'rect', b: 0.001, h: 0.005 }, hostV: { shape: 'rect', b: 0.001, h: 0.005 } };
  const md = Eng.buildModel(base, al, { ne: 4, path: 'std', family: 'oop', beam: 'tim' });
  const pl = Hom.homogenize(Hom.cellFromSpec(base, al, {}), al, { beam: 'tim' }), k = 0.005 * Math.PI / 0.014;
  check('hexagon 14 mm: plate vs workbench Bloch FE', Math.sqrt(pl.Dfree[0][0] * k ** 4 / pl.mArea) / TAU, Eng.solveK(md, k, 0, 1).f[0], 1e-3);
  check('hexagon is isotropic: D11 = D22', pl.Dfree[1][1], pl.Dfree[0][0], 1e-9);
}

console.log('— 3. Homogenised plate + PWE vs full Bloch FE, local-resonance gaps');
{ // your hexagon with the 6 x 6 x 5 mm block arm: workbench Timoshenko Bloch FE gives 2091.8–2961.9 Hz
  const G = { lattice: 'hexagon', l: 0.014 / Math.sqrt(3), h: 0.014 / Math.sqrt(3), theta: 30 };
  const sections = { host: { shape: 'rect', b: 0.001, h: 0.005 }, hostV: { shape: 'rect', b: 0.001, h: 0.005 }, arm: { shape: 'rect', b: 0.002, h: 0.001 }, stub: { shape: 'rect', b: 0.002, h: 0.005 }, tip: { shape: 'rect', b: 0.006, h: 0.005 }, lig: { shape: 'rect', b: 0.0015, h: 0.001 } };
  const R = { type: 'arm', attach: 'mid', t: 0.3, walls: [1], orient: 'in', angle: 0, tipForm: 'block', bL: 0.006, blockPos: 'centred', wallHalf: () => 0.0005, shape: 'straight', L1: 0.0065 };
  const spec = Eng.addResonators(Eng.baseCell(G), R); spec.sections = sections;
  const al2 = { E: 68.9e9, nu: 0.33, rho: 2700 }, cell = Hom.cellFromSpec(spec, al2, { beam: 'tim', ne: 8, nModes: 6, base: 'wall' });
  const pl = Hom.homogenize(cell, al2, { beam: 'tim' }); cell.res.forEach(r => r.use = Hom.pickModes(r, 'modal', 4));
  const Lat = Hom.latticeFor(cell, 'std'), P = Hom.pweSetup(Lat, pl, cell.res, { fam: 'flex', NF: 5, r0: 0.0005, extraMass: cell.extraMass });
  const F = Eng.kPath(Lat, 12).ks.map(k => Hom.pweAt(P, k[0], k[1], 6).f), gp = Hom.gapsOf(F)[0];
  check('hexagon: gap lower edge vs Bloch FE 2091.8 Hz', gp.lo, 2091.8, 0.015);
  check('hexagon: gap upper edge vs Bloch FE 2961.9 Hz', gp.hi, 2961.9, 0.015);
}

console.log('— 4. Drawn cells (draw.js): shell element and 3D Bloch FE');
{ const Drw = require('../src/draw.js'), al3 = { E: 70e9, nu: 0.3, rho: 2700 }, t = 1;
  const md = Drw.build({ a1: [10, 0], a2: [0, 10], mesh: 2.5, parts: [{ kind: 'plate', poly: [[0, 0], [10, 0], [10, 10], [0, 10]], z: 0, t, role: 'host' }] }, al3, {});
  const pl = Drw.homogenize(md), Db = al3.E * 1e-9 / 12 / (1 - 0.09);
  check('shell plate patch test: D11 = Et³/12(1−ν²)', pl.Dfree[0][0], Db, 1e-9); check('shell plate patch test: D66 = Gt³/12', pl.Dfree[2][2], Db * 0.35, 1e-9); check('shell plate: A11 = Et/(1−ν²)', pl.A[0][0], al3.E * 1e-3 / 0.91, 1e-9);
  const md2 = Drw.build({ a1: [10, 0], a2: [0, 10], mesh: 1.25, parts: [{ kind: 'plate', poly: [[0, 0], [10, 0], [10, 10], [0, 10]], z: 0, t, role: 'host' }] }, al3, {}), k = 0.5 * Math.PI / 0.01;
  check('shell plate Bloch FE flexural wave, k = π/2a', Drw.blochAt(md2, k, 0, 2, 3000).f[0], Math.sqrt(Db / (2700e-3)) * k * k / TAU, 3e-3);
  // your hexagon drawn as beams (imported from the left panel): must reproduce the workbench Bloch FE gaps
  const G = { lattice: 'hexagon', l: 0.014 / Math.sqrt(3), h: 0.014 / Math.sqrt(3), theta: 30 };
  const sections = { host: { shape: 'rect', b: 0.001, h: 0.005 }, hostV: { shape: 'rect', b: 0.001, h: 0.005 }, arm: { shape: 'rect', b: 0.002, h: 0.001 }, stub: { shape: 'rect', b: 0.002, h: 0.005 }, tip: { shape: 'rect', b: 0.006, h: 0.005 }, lig: { shape: 'rect', b: 0.0015, h: 0.001 } };
  const R = { type: 'arm', attach: 'mid', t: 0.3, walls: [1], orient: 'in', angle: 0, tipForm: 'block', bL: 0.006, blockPos: 'centred', wallHalf: () => 0.0005, shape: 'straight', L1: 0.0065 };
  const spec = Eng.addResonators(Eng.baseCell(G), R); spec.sections = sections;
  const mdh = Drw.build(Drw.fromSpec(spec), { E: 68.9e9, nu: 0.33, rho: 2700 }, {}), Lat = Hom.latticeFor({ a1: mdh.a1, a2: mdh.a2 }, 'std');
  const rs = Eng.kPath(Lat, 8).ks.map(q => Drw.blochAt(mdh, q[0], q[1], 8, 3000)), Fz = rs.map(r => r.f.filter((f, b) => r.z[b] >= 0.5)), nz = Math.min(...Fz.map(r => r.length)), gz = Hom.gapsOf(Fz.map(r => r.slice(0, nz)))[0];
  check('drawn hexagon (3D FE): gap lower edge vs workbench 2091.8 Hz', gz.lo, 2091.8, 2e-3); check('drawn hexagon (3D FE): gap upper edge vs workbench 2961.9 Hz', gz.hi, 2961.9, 2e-3);
}
console.log('— 5. Complex bands (phase and attenuation along Γ–X)');
{ const Drw = require('../src/draw.js'), G = { lattice: 'hexagon', l: 0.014 / Math.sqrt(3), h: 0.014 / Math.sqrt(3), theta: 30 };
  const sections = { host: { shape: 'rect', b: 0.001, h: 0.005 }, hostV: { shape: 'rect', b: 0.001, h: 0.005 }, arm: { shape: 'rect', b: 0.002, h: 0.001 }, stub: { shape: 'rect', b: 0.002, h: 0.005 }, tip: { shape: 'rect', b: 0.006, h: 0.005 }, lig: { shape: 'rect', b: 0.0015, h: 0.001 } };
  const R = { type: 'arm', attach: 'mid', t: 0.3, walls: [1], orient: 'in', angle: 0, tipForm: 'block', bL: 0.006, blockPos: 'centred', wallHalf: () => 0.0005, shape: 'straight', L1: 0.0065 };
  const spec = Eng.addResonators(Eng.baseCell(G), R); spec.sections = sections; const al4 = { E: 68.9e9, nu: 0.33, rho: 2700 };
  // FE: drawn-cell wave finite elements against the workbench WFE (out-of-plane DOFs of the flat cell)
  const md = Drw.build(Drw.fromSpec(spec), al4, {}), W = Drw.wfeModel(md, i => i >= 6 * md.N || [2, 3, 4].includes(i % 6)), Wb = Eng.buildWFE(Eng.buildModel(spec, al4, { ne: 4, path: 'std', family: 'oop', beam: 'tim' }));
  const minA = l => Math.min(...l.map(e => e.att));
  check('FE complex band, attenuation at 2.5 kHz (in the gap): drawn cell vs workbench', minA(Drw.wfeAt(W, 2500)), minA(Eng.wfeAt(Wb, 2500)), 0.01);
  // PWE: at a frequency on the real band, the complex solution must contain that real wave number
  const cell = Hom.cellFromSpec(spec, al4, { beam: 'tim', ne: 8, nModes: 6, base: 'wall' }), pl = Hom.homogenize(cell, al4, { beam: 'tim' }); cell.res.forEach(r => r.use = Hom.pickModes(r, 'modal', 4));
  const Lat = Hom.latticeFor(cell, 'std'), P = Hom.pweSetup(Lat, pl, cell.res, { fam: 'flex', NF: 2, r0: 0, extraMass: cell.extraMass }), X = [0.5 * Lat.b1[0], 0.5 * Lat.b1[1]], f = Hom.pweAt(P, 0.3 * X[0], 0.3 * X[1], 1).f[0];
  const ph = Hom.pweComplexAt(P, Lat, f).filter(e => e.att < 1e-6).map(e => e.ph);
  check('PWE complex band reproduces the real band (phase 0.3π at its frequency)', ph.reduce((b, v) => (Math.abs(v - 0.3 * Math.PI) < Math.abs(b - 0.3 * Math.PI) ? v : b), 0), 0.3 * Math.PI, 1e-4);
  const aP = minA(Hom.pweComplexAt(Hom.pweSetup(Lat, pl, cell.res, { fam: 'flex', NF: 2, r0: 0.0005, extraMass: cell.extraMass }), Lat, 2500));
  check('in-gap attenuation at 2.5 kHz: PWE vs FE (Np per cell)', aP, minA(Eng.wfeAt(Wb, 2500)), 0.15);
}
console.log('— 6. Hexagonal frame templates: central mass on six spokes, on a membrane (PWE vs FE)');
{ const Drw = require('../src/draw.js'), L = 10, r6 = v => +v.toFixed(6), V = (rad, k) => [r6(rad * Math.cos(k * Math.PI / 3)), r6(rad * Math.sin(k * Math.PI / 3)), 0], al5 = { E: 68.9e9, nu: 0.33, rho: 2700 };
  for (const ty of ['spokes', 'membrane']) {
    const parts = []; for (const k of [0, 1, 2]) parts.push({ kind: 'beam', p: [V(L, k), V(L, k + 1)], b: 1, h: 3, role: 'host' });
    const hub = [0, 1, 2, 3, 4, 5].map(k => V(2.5, k).slice(0, 2));
    if (ty === 'spokes') for (let k = 0; k < 6; k++) parts.push({ kind: 'beam', p: [V(2.5, k), V(L, k)], b: 0.6, h: 0.6, role: 'res' });
    else parts.push({ kind: 'plate', poly: [0, 1, 2, 3, 4, 5].map(k => V(L, k).slice(0, 2)), holes: [hub], z: 0, t: 0.1, role: 'res' });
    parts.push({ kind: 'plate', poly: hub, holes: [], z: 0, t: 2, role: 'res' });
    const md = Drw.build({ a1: [15, r6(Math.sqrt(3) / 2 * L)], a2: [0, r6(Math.sqrt(3) * L)], mesh: 2, parts }, al5, {}), Lat = Hom.latticeFor({ a1: md.a1, a2: md.a2 }, 'auto'), ks = Eng.kPath(Lat, 5).ks;
    const FE = ks.map(k => Drw.blochAt(md, k[0], k[1], 10, 3000)), Fz = FE.map(r => r.f.filter((f, b) => r.z[b] >= 0.5)), nz = Math.min(...Fz.map(r => r.length)), gF = Hom.gapsOf(Fz.map(r => r.slice(0, nz)))[0];
    const pl = Drw.homogenize(md), res = Drw.resonators(md, 12, 0.005); res.forEach(r => r.use = Hom.pickModes(r, 'modal', 4));
    const P = Hom.pweSetup(Lat, pl, res, { fam: 'flex', NF: 4, r0: 0.0005, extraMass: -res.reduce((s, r) => s + (r.hostIn || 0), 0) }), gP = Hom.gapsOf(ks.map(k => Hom.pweAt(P, k[0], k[1], 8).f))[0];
    check(`template ${ty}: gap lower edge, PWE vs FE`, gP.lo, gF.lo, 0.03); check(`template ${ty}: gap upper edge, PWE vs FE`, gP.hi, gF.hi, 0.05);
  }
}
console.log(`\n${fails ? fails + ' of ' + n + ' checks FAILED' : 'All ' + n + ' checks passed'} in ${((Date.now() - t0) / 1000).toFixed(1)} s`);
process.exitCode = fails ? 1 : 0;
