// Command-line runner: the same computation as the workbench page, without a browser.
// Usage:  node examples/run_analysis.js [design.json] [output folder]
// Writes, for out-of-plane (oop) and in-plane (ip) waves:
//   bands_<fam>.csv     Bloch frequencies along the zone path, and the resonator energy share of each point
//   complex_<fam>.csv   complex band along Gamma-X: phase Re(mu) and attenuation Im(mu) of every wave at each frequency
//   gaps.csv            complete gaps of each polarisation and the full gap (both), with their type
const fs = require('fs'), path = require('path');
const Eng = require(path.join(__dirname, '..', 'src', 'engine.js'));
const TAU = 2 * Math.PI;
const cfgFile = process.argv[2] || path.join(__dirname, 'design.json');
const outDir = process.argv[3] || path.join(process.cwd(), 'results');
const C = JSON.parse(fs.readFileSync(cfgFile, 'utf8'));
fs.mkdirSync(outDir, { recursive: true });
const mm = v => v / 1000;

// ---- design -> cell specification (identical to buildSpec/readParams in src/app.js)
function buildSpec(C) {
  if (C.custom_cell_file) return Eng.specFromJSON(fs.readFileSync(path.resolve(path.dirname(cfgFile), C.custom_cell_file), 'utf8'));
  const G = { lattice: C.lattice, a: mm(C.lattice_constant_mm || 0), l: mm(C.hexagon ? C.hexagon.l_mm : 0), h: mm(C.hexagon ? C.hexagon.h_mm : 0), theta: C.hexagon ? C.hexagon.theta_deg : 30 };
  const w = C.walls, vf = w.vertical_wall_width_factor || 1, r = C.resonator || { type: 'none' };
  const sections = {
    host: { shape: w.shape, b: mm(w.width_mm), h: mm(w.height_mm) },
    hostV: { shape: w.shape, b: w.shape === 'rect' ? mm(w.width_mm) * vf : mm(w.width_mm), h: w.shape === 'circ' ? mm(w.height_mm) * vf : mm(w.height_mm) },
    arm: { shape: 'rect', b: mm(r.arm_width_mm || 1), h: mm(r.arm_thickness_mm || 1) },
    stub: { shape: 'rect', b: mm(r.arm_width_mm || 1), h: mm(w.height_mm) },
    tip: { shape: 'rect', b: mm(r.block_width_mm || 1), h: mm(r.block_height_mm || 1) },
    lig: { shape: 'rect', b: mm(r.ligament_width_mm || 1), h: mm(r.ligament_thickness_mm || 1) }
  };
  if (r.block_density_kg_m3) sections.tip.rho = r.block_density_kg_m3;
  const R = { type: r.type, attach: r.attach || 'mid', t: r.position_along_wall || 0.5, walls: r.walls || null, orient: r.orient || 'in', angle: r.angle_deg || 0,
    tipForm: r.tip || 'block', bL: mm(r.block_length_mm || 0), blockPos: r.block_position || 'centred',
    wallHalf: sec => { const q = sections[sec] || sections.host; return (q.shape === 'circ' ? q.h : q.b) / 2; },
    shape: r.shape || 'straight', L1: mm(r.arm_length_mm || 0), L2: mm(r.second_leg_mm || 0), turn: r.turn_deg || 90,
    mt: (r.tip_mass_g || 0) / 1000, side: mm(r.cube_side_mm || 0), m: (r.mass_g || 0) / 1000, k: r.k_N_per_m || 0,
    M: (r.central_mass_g || 0) / 1000, Mside: mm(r.central_cube_side_mm || 0), ligTo: r.ligaments_to || 'corners' };
  const spec = Eng.addResonators(Eng.baseCell(G), R); spec.sections = sections; return spec;
}

const spec = buildSpec(C), mat = { E: C.material.E_GPa * 1e9, nu: C.material.nu, rho: C.material.rho_kg_m3 }, sv = C.solver;
const results = {};
for (const fam of ['oop', 'ip']) {
  const t0 = Date.now();
  // 1. finite-element model of the unit cell (Timoshenko or Euler-Bernoulli beams)
  const md = Eng.buildModel(spec, mat, { ne: sv.elements_per_wall, path: sv.path, family: fam, beam: sv.beam });
  // 2. Bloch eigenproblem [K(k) - w^2 M(k)] u = 0 at every k-point on the zone path
  const p = Eng.kPath(md.Lat, sv.kpoints_per_segment), nb = Math.min(sv.bands, md.n - 1);
  const F = [], Rf = [];
  for (const k of p.ks) { const r = Eng.solveK(md, k[0], k[1], nb, false, true); F.push(r.f); Rf.push(r.frac); }
  const fmax = Math.min(sv.fmax_Hz[fam], Math.min(...F.map(r => r[r.length - 1])) * 0.98);
  // 3. complete gaps: max of band n below min of band n+1 everywhere on the path
  const gaps = Eng.findGaps(F, nb).filter(g => g.lo < fmax).map(g => {
    const share = Rf.reduce((a, r) => a + r[g.below - 1], 0) / Rf.length;
    return Object.assign(g, { type: share > 0.5 ? 'local resonance' : 'Bragg', resonator_share: share });
  });
  // 4. complex band (wave finite element along a1): phase and attenuation per cell
  const W = Eng.buildWFE(md), rows = [];
  for (let i = 1; i <= sv.complex_band_samples; i++) { const f = fmax * i / sv.complex_band_samples; rows.push({ f, ev: Eng.wfeAt(W, f) }); }
  results[fam] = { md, p, F, Rf, gaps, fmax, rows };
  // --- write CSV
  let csv = 'k_index,path_coordinate_1_per_m,point,' + Array.from({ length: nb }, (_, b) => `f${b + 1}_Hz`).join(',') + ',' + Array.from({ length: nb }, (_, b) => `res_share${b + 1}`).join(',') + '\n';
  p.ks.forEach((k, i) => { const [sg, u] = p.seg[i]; const lab = u < 1e-9 ? p.labels[sg] : u > 1 - 1e-9 ? p.labels[sg + 1] : '';
    csv += `${i},${p.xs[i].toFixed(4)},${lab},${F[i].map(v => v.toFixed(3)).join(',')},${Rf[i].map(v => v.toFixed(4)).join(',')}\n`; });
  fs.writeFileSync(path.join(outDir, `bands_${fam}.csv`), csv);
  let cc = 'f_Hz,phase_rad_per_cell,attenuation_Np_per_cell,attenuation_dB_per_cell\n';
  for (const r of rows) if (r.ev) for (const e of r.ev) if (e.att < 5) cc += `${r.f.toFixed(3)},${e.ph.toFixed(5)},${e.att.toFixed(5)},${(8.686 * e.att).toFixed(4)}\n`;
  fs.writeFileSync(path.join(outDir, `complex_${fam}.csv`), cc);
  console.log(`${fam === 'oop' ? 'Out-of-plane' : 'In-plane    '}: ${md.n} complex DOFs, path ${p.labels.join('-')}, ${nb} bands up to ${Math.round(fmax)} Hz in ${((Date.now() - t0) / 1000).toFixed(1)} s`);
  gaps.forEach(g => console.log(`    gap ${g.lo.toFixed(1)}-${g.hi.toFixed(1)} Hz (${(200 * (g.hi - g.lo) / (g.hi + g.lo)).toFixed(1)} %), ${g.type}`));
  if (!gaps.length) console.log('    no complete gap');
}
// 5. full band gap: intersection of the two polarisations' gaps (over the range both cover)
const top = Math.min(results.oop.fmax, results.ip.fmax), full = [];
for (const a of results.oop.gaps) for (const b of results.ip.gaps) { const lo = Math.max(a.lo, b.lo), hi = Math.min(a.hi, b.hi, top); if (hi > lo) full.push({ lo, hi }); }
console.log(full.length ? `Full gap (both polarisations): ${full.map(g => `${g.lo.toFixed(1)}-${g.hi.toFixed(1)} Hz`).join(', ')}` : `Full gap (both polarisations): none below ${Math.round(top)} Hz`);
let gc = 'polarisation,between_bands,f_low_Hz,f_high_Hz,relative_width_percent,type,resonator_energy_share\n';
for (const fam of ['oop', 'ip']) for (const g of results[fam].gaps) gc += `${fam},${g.below}-${g.below + 1},${g.lo.toFixed(2)},${g.hi.toFixed(2)},${(200 * (g.hi - g.lo) / (g.hi + g.lo)).toFixed(2)},${g.type},${g.resonator_share.toFixed(3)}\n`;
for (const g of full) gc += `both,,${g.lo.toFixed(2)},${g.hi.toFixed(2)},${(200 * (g.hi - g.lo) / (g.hi + g.lo)).toFixed(2)},full gap,\n`;
fs.writeFileSync(path.join(outDir, 'gaps.csv'), gc);
console.log('CSV files written to', outDir);
