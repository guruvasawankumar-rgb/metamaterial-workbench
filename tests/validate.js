// Validation of the computational engine against analytical results.
// Usage: node tests/validate.js      (exit code 0 when every check passes)
const path = require('path');
const Eng = require(path.join(__dirname, '..', 'src', 'engine.js'));
const TAU = 2 * Math.PI, mat = { E: 68.9e9, nu: 0.33, rho: 2700 }, G = mat.E / (2 * (1 + mat.nu));
let allOk = true;
const check = (name, err, tol, detail) => { const ok = err <= tol; allOk = allOk && ok; console.log(`${ok ? 'PASS' : 'FAIL'}  ${name}: ${detail} (error ${(100 * err).toFixed(3)} %, limit ${(100 * tol).toFixed(2)} %)`); };
const nearest = (arr, x) => arr.reduce((p, q) => Math.abs(q - x) < Math.abs(p - x) ? q : p);

// 1. Timoshenko element: flexural wave on an infinite straight beam vs the exact Timoshenko dispersion relation
function beamSpec(a, b, h) {
  return { a1: [a, 0], a2: [0, a], centre: [a / 2, 0], std: null, nodes: [{ x: 0, y: 0, role: 'host', m: 0, J: 0, side: 0 }],
    beams: [{ i: 0, j: 0, n: [1, 0], sec: 'host', role: 'host' }], springs: [], sections: { host: { shape: 'rect', b, h } }, walls: [] };
}
function exactTimoshenko(k, I, A, kap) {
  const a2 = mat.rho ** 2 * I / (kap * G), a1 = mat.rho * A + k * k * mat.rho * I * (1 + mat.E / (kap * G)), a0 = mat.E * I * k ** 4;
  return Math.sqrt((a1 - Math.sqrt(a1 * a1 - 4 * a2 * a0)) / (2 * a2)) / TAU;
}
{
  const a = 0.008, b = 0.001, h = 0.005, I = b * h ** 3 / 12, A = b * h, kap = Eng.shearCoef('rect', mat.nu);
  for (const fam of ['oop', 'ip']) {
    const sp = fam === 'oop' ? beamSpec(a, b, h) : beamSpec(a, h, b);
    const md = Eng.buildModel(sp, mat, { ne: 8, family: fam, beam: 'tim' });
    let err = 0;
    for (const fr of [0.1, 0.4, 0.7, 1.0]) { const k = fr * Math.PI / a, fe = exactTimoshenko(k, I, A, kap); err = Math.max(err, Math.abs(nearest(Eng.solveK(md, k, 0, 3).f, fe) / fe - 1)); }
    check(`Timoshenko beam dispersion, ${fam}`, err, 0.005, `1 x 5 mm beam, k from 0.1 to 1.0 pi/a`);
  }
  // Euler-Bernoulli option: exact EB dispersion w = k^2 sqrt(EI / rho A)
  const md = Eng.buildModel(beamSpec(a, b, h), mat, { ne: 8, family: 'oop', beam: 'eb' }); let err = 0;
  for (const fr of [0.1, 0.5, 1.0]) { const k = fr * Math.PI / a, fe = k * k * Math.sqrt(mat.E * I / (mat.rho * A)) / TAU; err = Math.max(err, Math.abs(nearest(Eng.solveK(md, k, 0, 3).f, fe) / fe - 1)); }
  check('Euler-Bernoulli beam dispersion', err, 0.005, 'same beam, Timoshenko switched off');
}
// 2. In-plane: square lattice, long-wavelength longitudinal wave c = sqrt(E / (2 rho))
{
  const a = 0.01, s = Eng.baseCell({ lattice: 'square', a }); s.sections = { host: { shape: 'rect', b: 0.0005, h: 0.002 } };
  const md = Eng.buildModel(s, mat, { ne: 8, family: 'ip', beam: 'tim' }), k = 2, fL = Math.sqrt(mat.E / (2 * mat.rho)) * k / TAU;
  check('In-plane longitudinal wave, square lattice', Math.abs(nearest(Eng.solveK(md, k, 0, 4).f, fL) / fL - 1), 0.002, `c = sqrt(E/2rho)`);
}
// 3. In-plane: regular honeycomb shear wave vs Gibson-Ashby, G* = 0.5774 E (t/l)^3, rho* = 1.1547 rho t/l
{
  const l = 0.008083, tl = 0.05, s = Eng.baseCell({ lattice: 'hexagon', l, h: l, theta: 30 });
  s.sections = { host: { shape: 'rect', b: tl * l, h: 0.005 }, hostV: { shape: 'rect', b: tl * l, h: 0.005 } };
  const md = Eng.buildModel(s, mat, { ne: 8, family: 'ip', beam: 'eb' }), k = 2;
  const fT = Math.sqrt(0.5774 * mat.E * tl ** 3 / (1.1547 * mat.rho * tl)) * k / TAU;
  check('In-plane shear wave, honeycomb vs Gibson-Ashby', Math.abs(Eng.solveK(md, k, 0, 3).f[0] / fT - 1), 0.01, 't/l = 0.05');
}
// 4. Out-of-plane: bare square grillage at long wavelength vs equivalent plate, D = E I / a per unit width in each direction
{
  const a = 0.05, b = 0.004, h = 0.002, s = Eng.baseCell({ lattice: 'square', a }); s.sections = { host: { shape: 'rect', b, h } };
  const md = Eng.buildModel(s, mat, { ne: 4, family: 'oop', beam: 'eb' }), k = 5, EI = mat.E * b * h ** 3 / 12, rA = mat.rho * b * h;
  const fp = k * k * Math.sqrt(EI / (2 * rA)) / TAU;
  check('Out-of-plane grillage vs equivalent plate', Math.abs(Eng.solveK(md, k, 0, 2).f[0] / fp - 1), 0.002, 'square lattice, k = 5 rad/m');
}
// 5. Complex band (wave finite element) vs Bloch bands: every real Bloch solution appears as a propagating wave
const l = 0.014 / Math.sqrt(3);
const secs = { host: { shape: 'rect', b: 0.001, h: 0.005 }, hostV: { shape: 'rect', b: 0.001, h: 0.005 }, arm: { shape: 'rect', b: 0.002, h: 0.001 }, stub: { shape: 'rect', b: 0.002, h: 0.005 }, tip: { shape: 'rect', b: 0.006, h: 0.005 }, lig: { shape: 'rect', b: 0.0015, h: 0.001 } };
const userR = { type: 'arm', attach: 'mid', walls: [1], orient: 'in', shape: 'straight', L1: 0.0065, tipForm: 'block', bL: 0.006, blockPos: 'centred', wallHalf: () => 0.0005 };
const designs = [
  ['your hexagon with block arm', Eng.addResonators(Eng.baseCell({ lattice: 'hexagon', l, h: l, theta: 30 }), userR)],
  ['square with mass-on-spring', Eng.addResonators(Eng.baseCell({ lattice: 'square', a: 0.014 }), { type: 'spring', attach: 'junction', k: 3e4, m: 0.001 })],
  ['re-entrant with central mass', Eng.addResonators(Eng.baseCell({ lattice: 'hexagon', l: 0.008, h: 0.012, theta: -20 }), { type: 'central', M: 0.002, Mside: 0.004, ligTo: 'corners' })]
];
for (const [name, spec] of designs) for (const fam of ['oop', 'ip']) {
  spec.sections = secs;
  const md = Eng.buildModel(spec, mat, { ne: 4, family: fam, beam: 'tim' }), W = Eng.buildWFE(md), b1 = md.Lat.b1; let me = 0, standing = [];
  const cut = md.chains.filter(ch => ch.beam.n[0] !== 0).map(ch => ch.chain[Math.floor(ch.eidx.length / 2)]);
  for (const mu of [0.4, 1.3, 2.6]) { const r = Eng.solveK(md, b1[0] * mu / TAU, b1[1] * mu / TAU, 8, true);
    r.f.forEach((f, i) => { if (f < 1) return;
      // standing modes whose cut faces do not move carry nothing across the cell boundary: invisible to the transfer method (see Method text)
      const v = r.vecs[i], amp = q => Math.hypot(v.re[q], v.im[q]); let face = 0, all = 0;
      for (let q = 0; q < 3 * md.N; q++) all = Math.max(all, amp(q)); for (const nd of cut) for (let c = 0; c < 3; c++) face = Math.max(face, amp(3 * nd + c));
      if (face < 1e-8 * all) { standing.push(f); return; }
      const ev = Eng.wfeAt(W, f); if (!ev) return; me = Math.max(me, Math.min(...ev.map(e => Math.abs(e.ph - mu) + e.att))); }); }
  const sf = [...new Set(standing.map(f => f.toFixed(0)))];
  check(`Complex band vs Bloch, ${name}, ${fam}`, me, 1e-5, 'phase per cell' + (sf.length ? `; standing mode with motionless cut faces at ${sf.join(', ')} Hz excluded` : ''));
}
// 6. Reference result for your sample (Timoshenko, Gamma-X-K-Gamma)
{
  const spec = designs[0][1]; spec.sections = secs;
  const md = Eng.buildModel(spec, mat, { ne: 4, path: 'stdforce', family: 'oop', beam: 'tim' }), p = Eng.kPath(md.Lat, 16);
  const g = Eng.findGaps(p.ks.map(k => Eng.solveK(md, k[0], k[1], 10).f), 10)[0];
  check('Reference: first out-of-plane gap of your sample', Math.max(Math.abs(g.lo / 2091.8 - 1), Math.abs(g.hi / 2961.9 - 1)), 0.001, `${g.lo.toFixed(1)}-${g.hi.toFixed(1)} Hz (expected 2091.8-2961.9 Hz)`);
  // the one-sided arm breaks the six-fold symmetry: 'std' now falls back to the half zone, whose longer path finds the band maximum
  const md2 = Eng.buildModel(spec, mat, { ne: 4, path: 'std', family: 'oop', beam: 'tim' }), p2 = Eng.kPath(md2.Lat, 16);
  check('Standard path falls back to the half zone for this asymmetric cell', md2.Lat.path.mode === 'half' ? 0 : 1, 0, md2.Lat.path.labels.join('-'));
  const g2 = Eng.findGaps(p2.ks.map(k => Eng.solveK(md2, k[0], k[1], 10).f), 10)[0];
  check('Gap on the half zone lies inside the standard-path gap', (g2.lo >= g.lo - 1e-6 && g2.hi <= g.hi + 1e-6) ? 0 : 1, 0, `${g2.lo.toFixed(1)}-${g2.hi.toFixed(1)} Hz`);
}
console.log(allOk ? '\nALL CHECKS PASSED' : '\nSOME CHECKS FAILED'); process.exit(allOk ? 0 : 1);
