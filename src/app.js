(() => {
const $ = s => document.querySelector(s);
const TAU = Eng.TAU;
const S = { P: null, md: null, band: null, cplx: null, gaps: [], sel: null, mode: null, hoverF: null, running: false, sweep: null, playing: true, phase: 0, view: 'oop', fam: null, complete: [], shared: null };
const FAMS = ['oop', 'ip'], FAM_LABEL = { oop: 'Out-of-plane', ip: 'In-plane' }, FAM_WAVE = { oop: 'flexural (out-of-plane)', ip: 'in-plane' };
const MATS = { al: [68.9, 0.33, 2700], ti: [113.8, 0.342, 4430], ss: [193, 0.29, 8000], cf: [60, 0.3, 1550], pa: [1.7, 0.4, 1010], pla: [3.5, 0.36, 1240], resin: [2.6, 0.35, 1250], esp: [0.05, 0.3, 400] };
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const num = id => parseFloat(document.getElementById(id).value);

function parseWalls(t) { t = (t || '').trim().toLowerCase(); if (!t || t === 'all') return null; const w = t.split(/[^0-9]+/).filter(Boolean).map(Number); return w.length ? w : null; }
function readParams(latOverride) {
  const lat = latOverride || $('#lat').value, shape = $('#hshape').value, b = num('hb') / 1000, hh = num('hh') / 1000, vf = num('vf');
  const P = {
    lattice: lat, ne: parseInt($('#ne').value, 10), path: $('#path').value, beam: $('#beam').value,
    mat: { E: num('E') * 1e9, nu: num('nu'), rho: num('rho'), res: $('#rmatp').value === 'same' ? null : { E: num('rE') * 1e9, nu: num('rnu'), rho: num('rrho') }, mem: $('#mmatp').value === 'same' ? null : { E: num('mE') * 1e9, nu: num('mnu'), rho: num('mrho') } },
    solver: { nb: Math.round(num('nb')), kper: Math.round(num('kper')), nf: Math.round(num('nf')), fmax: num('fmax'), fmin: Math.max(0, num('fmin') || 0) }
  };
  if (P.path === 'custom') { P.kpath = [$('#kpU').value, $('#kpT').value]; try { Eng.setCustomPath({ units: $('#kpU').value, pts: Eng.parseKPath($('#kpT').value) }); } catch (e) { Eng.setCustomPath(null); P.kpathErr = e.message; } }
  if (lat === 'custom') { P.custom = $('#json').value; return P; }
  if (lat === 'drawn') return P;
  P.G = { lattice: lat, a: num('a') / 1000, l: num('hl') / 1000, h: num('hv') / 1000, theta: num('hth') };
  P.sections = { host: { shape, b, h: hh }, hostV: { shape, b: shape === 'rect' ? b * vf : b, h: shape === 'circ' ? hh * vf : hh },
    arm: { shape: 'rect', b: num('rb') / 1000, h: num('rh') / 1000 }, stub: { shape: 'rect', b: num('rb') / 1000, h: hh }, tip: { shape: 'rect', b: num('bW') / 1000, h: num('bH') / 1000 }, lig: { shape: 'rect', b: num('cb') / 1000, h: num('ch') / 1000 } };
  P.vf = vf;
  if (P.mat.res) for (const k of ['arm', 'tip', 'lig']) Object.assign(P.sections[k], { E: P.mat.res.E, rho: P.mat.res.rho, nu: P.mat.res.nu });
  if (num('bRho') > 0) P.sections.tip.rho = num('bRho');
  P.R = { type: $('#rtype').value, attach: $('#ratt').value, t: num('rt'), walls: parseWalls($('#rwalls').value), orient: $('#rori').value, angle: num('rang'),
    tipForm: $('#tipf').value, bL: num('bL') / 1000, blockPos: $('#bPos').value,
    wallHalf: sec => { const q = P.sections[sec] || P.sections.host; return (q.shape === 'circ' ? q.h : q.b) / 2; },
    shape: $('#rshape').value, L1: num('rL') / 1000, L2: num('rL2') / 1000, turn: num('rturn'), mt: num('rmt') / 1000, side: num('rside') / 1000,
    m: num('rm') / 1000, k: num('rk'), M: num('cM') / 1000, Mside: num('cS') / 1000, ligTo: $('#cTo').value };
  return P;
}
function validate(P) {
  const pos = (v, name) => (!(v > 0) ? `${name} must be a positive number.` : null);
  const chk = [pos(P.mat.E, 'Young’s modulus'), pos(P.mat.rho, 'Density'), (P.mat.nu > -1 && P.mat.nu < 0.5) ? null : 'Poisson’s ratio must lie between −1 and 0.5.'];
  if (P.mat.res) chk.push(pos(P.mat.res.E, 'Resonator Young’s modulus'), pos(P.mat.res.rho, 'Resonator density'), (P.mat.res.nu > -1 && P.mat.res.nu < 0.5) ? null : 'Resonator Poisson’s ratio must lie between −1 and 0.5.');
  if (P.mat.mem) chk.push(pos(P.mat.mem.E, 'Membrane Young’s modulus'), pos(P.mat.mem.rho, 'Membrane density'), (P.mat.mem.nu > -1 && P.mat.mem.nu < 0.5) ? null : 'Membrane Poisson’s ratio must lie between −1 and 0.5.');
  if (!(P.solver.nb >= 2 && P.solver.nb <= 40)) chk.push('Bands to compute must be between 2 and 40.');
  if (!(P.solver.kper >= 4 && P.solver.kper <= 80)) chk.push('k-points per segment must be between 4 and 80.');
  if (!(P.solver.nf >= 40 && P.solver.nf <= 2000)) chk.push('Frequency samples must be between 40 and 2000.');
  if (!(num('fmin') >= 0) && $('#fmin').value !== '') chk.push('Min frequency cannot be negative.');
  if (P.solver.fmax > 0 && P.solver.fmin >= P.solver.fmax) chk.push('Min frequency must be below max frequency.');
  if (P.kpathErr) chk.push(P.kpathErr);
  if (P.lattice === 'custom' || P.lattice === 'drawn') return chk.find(Boolean) || null;
  const G = P.G, S = P.sections, R = P.R;
  if (G.lattice === 'hexagon') chk.push(pos(G.l, 'Inclined wall length'), pos(G.h, 'Vertical wall length'), pos(P.vf, 'Vertical wall width factor'), Eng.checkBase(G)); else chk.push(pos(G.a, 'Lattice constant'));
  chk.push(pos(S.host.h, P.sections.host.shape === 'circ' ? 'Rod diameter' : 'Wall thickness'), S.host.shape === 'rect' ? pos(S.host.b, 'Wall width') : null);
  if (R.type === 'spring') chk.push(pos(R.m, 'Resonator mass'), pos(R.k, 'Spring stiffness'));
  if (R.type === 'arm') { chk.push(pos(R.L1, 'Arm length'), pos(S.arm.b, 'Arm width'), pos(S.arm.h, 'Arm thickness'));
    if (R.tipForm === 'block') chk.push(pos(R.bL, 'Block length'), pos(S.tip.b, 'Block width'), pos(S.tip.h, 'Block height')); else chk.push(R.mt >= 0 ? null : 'Tip mass cannot be negative.', R.side >= 0 ? null : 'Cube side cannot be negative.');
    if (R.shape === 'L') chk.push(pos(R.L2, 'Second leg length')); if (R.orient === 'angle' && !isFinite(R.angle)) chk.push('Enter an arm angle.'); }
  if ((R.type === 'arm' || R.type === 'spring') && R.attach === 't' && !(R.t > 0 && R.t < 1)) chk.push('Position along the wall must lie strictly between 0 and 1.');
  if (R.type === 'central') chk.push(pos(R.M, 'Central mass'), pos(S.lig.b, 'Ligament width'), pos(S.lig.h, 'Ligament thickness'), R.Mside >= 0 ? null : 'Cube side cannot be negative.');
  return chk.find(Boolean) || null;
}
function buildSpec(P) {
  if (P.lattice === 'custom') return Eng.specFromJSON(P.custom);
  const base = Eng.baseCell(P.G);
  if (P.R.walls) { const have = base.beams.map(b => b.wall); const bad = P.R.walls.filter(w => !have.includes(w)); if (bad.length) throw new Error(`This cell has walls 1 to ${have.length}; wall ${bad[0]} does not exist.`); }
  const spec = Eng.addResonators(base, P.R); spec.sections = P.sections; return spec;
}
function makeModel(P, scale, specOverride, family) { const spec = specOverride || buildSpec(P); return Eng.buildModel(spec, P.mat, { ne: P.ne, path: P.path, scale, family, beam: P.beam }); }
function resInfo(P, md, fam = 'oop') {
  const Ib = sec => fam === 'ip' ? sec.Iz : sec.I;
  const mu = md && md.hostMass ? md.resMass / md.hostMass : 0; let fr = null, approx = false;
  if (P.lattice !== 'custom') {
    const R = P.R, E = (P.mat.res || P.mat).E, rho = (P.mat.res || P.mat).rho;
    if (R.type === 'spring') fr = Math.sqrt(R.k / R.m) / TAU;
    else if (R.type === 'arm') {
      const sec = Eng.sectionProps(P.sections.arm); let L = R.L1 + (R.shape === 'L' ? R.L2 : 0), mt = R.mt; approx = R.shape === 'L';
      if (R.tipForm === 'block') { const tp = P.sections.tip, bm = (tp.rho || rho) * tp.b * tp.h * R.bL; mt = bm; L = R.blockPos === 'centred' ? L : L + R.bL / 2; approx = true; }
      fr = Math.sqrt(3 * E * Ib(sec) / L ** 3 / (mt + 0.24 * rho * sec.A * L)) / TAU;
    }
    else if (R.type === 'central' && md) { const lig = md.chains.filter(c => c.beam.sec === 'lig'); if (lig.length) { const sec = Eng.sectionProps(P.sections.lig), L = lig[0].len, k = lig.length * 12 * E * Ib(sec) / L ** 3; fr = Math.sqrt(k / (R.M + 0.37 * lig.length * rho * sec.A * L)) / TAU; approx = true; } }
  }
  return { fr, mu, approx };
}
function freqCap(P, md, topMin, fam = 'oop') {
  if (P.solver.fmax > 0) return P.solver.fmax;
  const { fr, mu } = resInfo(P, md, fam);
  return fr ? Math.min(topMin * 0.98, 3 * fr * Math.sqrt(1 + mu)) : topMin * 0.98;
}
function fmtF(f) { if (!isFinite(f)) return '–'; if (f === 0) return '0 Hz'; return f >= 1000 ? `${(f / 1000).toFixed(f >= 1e4 ? 1 : 2)} kHz` : `${f.toFixed(f >= 100 ? 0 : 1)} Hz`; }
function fmtM(kg) { return kg >= 1 ? `${kg.toFixed(3)} kg` : `${(kg * 1000).toFixed(kg * 1000 >= 10 ? 1 : 2)} g`; }

// ---------- UI wiring ----------
const DESIGN_NOTE = {
  hexcant: 'Honeycomb of beams with a cantilever arm and tip block, set by the numbers below. To change it freely (move parts, add shapes or plates), press “Edit this cell freely”.',
  hexspokes: 'Hexagonal frame; a centre mass held by six spokes to the corners, an inner ring and six small blocks. Change the numbers below, or edit any part in the Unit cell tab.',
  hexbeam: 'Hexagonal frame; one beam across the cell carries a centre mass. Change the numbers below, or edit any part in the Unit cell tab.',
  own: 'Choose a lattice and a resonator, then set their sizes. The cell appears in the Unit cell tab, where every part can be moved, resized or added to.',
  custom: 'A planar beam cell typed as nodes and beams (JSON).' };
const LAT_LABEL = { hex: 'Hexagon side L', square: 'Lattice constant a', tri: 'Triangle side a', kagome: 'Triangle side L (hexagon side L)', free: 'Size' };
function syncDesign() { const d = $('#design').value; $('#lat').value = d === 'hexcant' ? 'hexagon' : d === 'custom' ? 'custom' : 'drawn'; }
function syncVisibility() {
  syncDesign();
  const des = $('#design').value, gen = ['hexspokes', 'hexbeam', 'own'].includes(des);
  $('#designNote').textContent = DESIGN_NOTE[des] || '';
  $('#genBox').classList.toggle('hidden', !gen); $('#gLatRow').classList.toggle('hidden', des !== 'own'); $('#gResRow').classList.toggle('hidden', des !== 'own');
  { const gl = des === 'own' ? $('#gLat').value : 'hex', gr = des === 'hexspokes' ? 'spokes' : des === 'hexbeam' ? 'across' : $('#gRes').value;
    document.querySelectorAll('[data-glat]').forEach(el => el.classList.toggle('hidden', !el.dataset.glat.split(' ').includes(gl)));
    document.querySelectorAll('[data-gres]').forEach(el => el.classList.toggle('hidden', !el.dataset.gres.split(' ').includes(gr)));
    $('#gLlab').textContent = LAT_LABEL[gl] || 'Size'; }
  document.querySelectorAll('[data-ek]').forEach(el => el.classList.toggle('hidden', el.dataset.ek !== $('#eKind').value));
  const lat = $('#lat').value, t = $('#rtype').value, custom = lat === 'custom', drawn = lat === 'drawn';
  if (!drawn) S.lastLat = lat;
  $('#drawnBox').classList.toggle('hidden', true); $('#latNote').classList.toggle('hidden', drawn);
  $('#rmatBox').classList.toggle('hidden', $('#rmatp').value === 'same');
  $('#mmatBox').classList.toggle('hidden', $('#mmatp').value === 'same');
  $('#aBox').classList.toggle('hidden', lat === 'hexagon' || custom || drawn);
  $('#hexBox').classList.toggle('hidden', lat !== 'hexagon');
  $('#customBox').classList.toggle('hidden', !custom);
  $('#toCustom').classList.toggle('hidden', custom || drawn); $('#editAsDraw').classList.toggle('hidden', drawn);
  $('#kpathBox').classList.toggle('hidden', $('#path').value !== 'custom');
  document.querySelectorAll('option[data-drawn]').forEach(o => { o.hidden = !drawn; o.disabled = !drawn; if (!drawn && o.selected) o.parentElement.selectedIndex = 0; });
  $('#spPartBox').classList.toggle('hidden', $('#sp').value !== 'part');
  $('#hostSet').classList.toggle('hidden', custom || drawn);
  $('#resSet').classList.toggle('hidden', custom || drawn);
  ['#ne'].forEach(q => $(q).closest('.row').classList.toggle('hidden', drawn));
  $('#vfRow').classList.toggle('hidden', lat !== 'hexagon');
  $('#attachBox').classList.toggle('hidden', !(t === 'arm' || t === 'spring'));
  $('#tRow').classList.toggle('hidden', $('#ratt').value !== 't');
  $('#wallRow').classList.toggle('hidden', $('#ratt').value === 'junction');
  $('#armBox').classList.toggle('hidden', t !== 'arm');
  $('#angRow').classList.toggle('hidden', $('#rori').value !== 'angle');
  $('#lBox').classList.toggle('hidden', $('#rshape').value !== 'L');
  const fromFace = $('#ratt').value !== 'junction';
  $('#rLlab').textContent = ($('#rshape').value === 'L' ? 'First leg length' : 'Arm length') + (fromFace ? ' from wall face' : '');
  $('#blockBox').classList.toggle('hidden', $('#tipf').value !== 'block');
  $('#pointBox').classList.toggle('hidden', $('#tipf').value !== 'point');
  $('#springBox').classList.toggle('hidden', t !== 'spring');
  $('#centralBox').classList.toggle('hidden', t !== 'central');
  const circ = $('#hshape').value === 'circ';
  $('#hbRow').classList.toggle('hidden', circ);
  $('#hhLabel').textContent = circ ? 'Rod diameter d' : 'Thickness h (out of plane)';
}
function editAsDrawing() {
  if (typeof DrawUI === 'undefined') return; if (isDrawn()) { showTab('geo'); return; }
  if (!DrawUI.importLeft()) return; $('#design').value = 'own'; $('#gLat').value = 'free'; $('#gRes').value = 'none'; onEdit(); showTab('geo');
  setStatus('The cell is now editable: every wall, arm and mass is a part. Move or resize parts with the mouse or type their numbers; add shapes, plates and shell walls. Your parametric settings are kept: choose the cell type again to return to them.');
}
$('#editAsDraw').addEventListener('click', editAsDrawing); $('#editAsDraw2').addEventListener('click', editAsDrawing);
$('#sp').addEventListener('change', syncVisibility); $('#eKind').addEventListener('change', () => { syncVisibility(); drawEdge(); });
$('#eView').addEventListener('change', () => { drawEdge(); drawEdgeMode(); });
$('#eZak').addEventListener('click', () => DrawUI.zakCheck());
$('#kpFill').addEventListener('click', () => {
  try { const P = readParams(), md = isDrawn() ? null : makeModel(Object.assign({}, P, { path: 'auto' }), undefined, undefined, 'oop'); let L;
    if (md) L = md.Lat; else { const keep = $('#path').value; $('#path').value = 'auto'; try { L = DrawUI.lattice(); } finally { $('#path').value = keep; } }
    const fr = v => { for (const d of [1, 2, 3, 4, 6]) { const n = Math.round(v * d); if (Math.abs(n / d - v) < 1e-6) return d === 1 ? String(n) : `${n}/${d}`; } return (+v.toFixed(4)).toString(); };
    $('#kpU').value = 'frac';
    $('#kpT').value = L.path.pts.map((k, i) => `${L.path.labels[i]} ${fr((k[0] * L.a1[0] + k[1] * L.a1[1]) / TAU)} ${fr((k[0] * L.a2[0] + k[1] * L.a2[1]) / TAU)}`).join('\n');
    onEdit(); } catch (e) { setStatus(e.message); }
});
$('#setLC').addEventListener('click', () => { const lc = num('hLC'); if (lc > 0) { const l = (lc / Math.sqrt(3)).toFixed(4); $('#hl').value = l; $('#hv').value = l; $('#hth').value = 30; onEdit(); } });
$('#toCustom').addEventListener('click', () => {
  const P = readParams(), err = validate(P); if (err) { setStatus(err); return; }
  try { $('#json').value = Eng.specToJSON(buildSpec(P)); $('#design').value = 'custom'; onEdit(); setStatus('Design copied into the custom editor. Edit nodes, beams and sections there.'); } catch (e) { setStatus(e.message); }
});
$('#matp').addEventListener('change', () => { const m = MATS[$('#matp').value]; if (m) { $('#E').value = m[0]; $('#nu').value = m[1]; $('#rho').value = m[2]; } onEdit(); });
['E', 'nu', 'rho'].forEach(id => document.getElementById(id).addEventListener('input', () => { $('#matp').value = 'custom'; }));
$('#rmatp').addEventListener('change', () => { const m = MATS[$('#rmatp').value]; if (m) { $('#rE').value = m[0]; $('#rnu').value = m[1]; $('#rrho').value = m[2]; } onEdit(); });
['rE', 'rnu', 'rrho'].forEach(id => document.getElementById(id).addEventListener('input', () => { $('#rmatp').value = 'custom'; }));
$('#mmatp').addEventListener('change', () => { const m = MATS[$('#mmatp').value]; if (m) { $('#mE').value = m[0]; $('#mnu').value = m[1]; $('#mrho').value = m[2]; } onEdit(); });
['mE', 'mnu', 'mrho'].forEach(id => document.getElementById(id).addEventListener('input', () => { $('#mmatp').value = 'custom'; }));
$('#tune').addEventListener('click', () => { const m = num('rm') / 1000, f = num('ftarget'); if (m > 0 && f > 0) { $('#rk').value = (((TAU * f) ** 2) * m).toPrecision(5); onEdit(); } });
document.querySelectorAll('.rail input, .rail select, .rail textarea').forEach(el => el.addEventListener('input', onEdit));
function onEdit() { const wasDrawn = S.wasDrawn, dr = isDrawn(); S.wasDrawn = dr; syncVisibility(); drawKPath(); if (wasDrawn !== dr) { const cur = tabs.find(t => $('#t-' + t).getAttribute('aria-selected') === 'true'); showTab(cur || 'geo'); } drawGeometry(); if (typeof DrawUI !== 'undefined') DrawUI.onRail(); }

const tabs = ['geo', 'disp', 'sens', 'edge', 'meth'];
tabs.forEach(t => $('#t-' + t).addEventListener('click', () => showTab(t)));
document.querySelector('.tabs').addEventListener('keydown', e => {
  const i = tabs.findIndex(t => $('#t-' + t).getAttribute('aria-selected') === 'true');
  if (e.key === 'ArrowRight' || e.key === 'ArrowLeft') { const j = (i + (e.key === 'ArrowRight' ? 1 : tabs.length - 1)) % tabs.length; showTab(tabs[j]); $('#t-' + tabs[j]).focus(); }
});
function isDrawn() { return $('#lat').value === 'drawn'; }
function showTab(t) {
  tabs.forEach(x => { $('#t-' + x).setAttribute('aria-selected', x === t); $('#p-' + x).classList.toggle('hidden', x !== t); });
  const dr = isDrawn();
  if (t === 'geo' && dr) { $('#p-geo').classList.add('hidden'); $('#p-drw').classList.remove('hidden'); } else $('#p-drw').classList.add('hidden');
  $('#dispParam').classList.toggle('hidden', dr); $('#dispDrawn').classList.toggle('hidden', !dr);
  if (t === 'disp') { $('#dispTitle').textContent = dr ? 'Dispersion of the drawn cell' : VIEW_TEXT[S.view][0]; $('#dispLede').textContent = dr ? 'Full Bloch finite elements of the drawn cell (beam and shell elements, six degrees of freedom per node). Left: bands coloured by resonator energy share; right: complex band along the first path leg. Click a band point to animate its mode.' : VIEW_TEXT[S.view][1]; }
  if (t === 'edge') { if (dr && typeof DrawUI !== 'undefined') DrawUI.edgeAssessment(); else renderAssessment(); }
  try { history.replaceState(null, '', '#' + t); } catch (e) {}
  redrawAll();
}
function setStatus(txt, frac) { $('#status').textContent = txt; $('#prog').style.width = frac == null ? '0' : `${Math.round(frac * 100)}%`; }
const tick = () => new Promise((r, j) => setTimeout(() => (S.stop ? j(new Error('stopped by the user')) : r()), 0));
function startRun() { S.running = true; S.stop = false; ['#run', '#sweep', '#edgeRun'].forEach(q => { const el = $(q); if (el) el.disabled = true; }); $('#stop').disabled = false; }
function endRun() { S.running = false; ['#run', '#sweep', '#edgeRun'].forEach(q => { const el = $(q); if (el) el.disabled = false; }); $('#stop').disabled = true; }
$('#stop').addEventListener('click', () => { if (S.running) { S.stop = true; setStatus('Stopping…'); } });

// ---------- canvas helpers ----------
function prep(cv) {
  const dpr = window.devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight;
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
  const ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}
function hex2rgb(c) { c = c.replace('#', ''); if (c.length === 3) c = c.split('').map(x => x + x).join(''); return [0, 2, 4].map(i => parseInt(c.substr(i, 2), 16)); }
function mix(a, b, t) { const A = hex2rgb(a), B = hex2rgb(b); return `rgb(${A.map((v, i) => Math.round(v + (B[i] - v) * t)).join(',')})`; }
let hatchCache = {};
function hatch(ctx, color) {
  if (!hatchCache[color]) { const p = document.createElement('canvas'); p.width = p.height = 7; const c = p.getContext('2d'); c.strokeStyle = color; c.lineWidth = 1.2; c.beginPath(); for (const o of [-7, 0, 7]) { c.moveTo(o, 7); c.lineTo(o + 7, 0); } c.stroke(); hatchCache[color] = p; }
  return ctx.createPattern(hatchCache[color], 'repeat');
}
function niceTicks(lo, hi, n) {
  const span = hi - lo; if (!(span > 0)) return [lo];
  const raw = span / n, mag = 10 ** Math.floor(Math.log10(raw)), st = [1, 2, 2.5, 5, 10].map(s => s * mag).find(s => s >= raw);
  const out = []; for (let v = Math.ceil(lo / st) * st; v <= hi + 1e-9 * span; v += st) out.push(+v.toPrecision(12)); return out;
}
function fUnit(fmax) { return fmax >= 2000 ? { d: 1000, u: 'kHz' } : { d: 1, u: 'Hz' }; }

// ---------- geometry ----------
function drawGeometry() {
  if (isDrawn()) return;
  const cv = $('#cvGeo'); if ($('#p-geo').classList.contains('hidden')) { infoOnly(); return; }
  const { ctx, w, h } = prep(cv);
  const P = readParams(), err = validate(P);
  if (P.lattice === 'custom') $('#jsonNote').innerHTML = 'Lengths in mm, masses in g, J in g·mm². “shift” is the cell offset [n₁, n₂] of node j; role “res” marks resonator parts.';
  const fail = msg => { if (P.lattice === 'custom') $('#jsonNote').innerHTML = `<span class="err">${msg}</span>`; ctx.fillStyle = css('--danger'); ctx.font = '600 15px system-ui, sans-serif'; wrapText(ctx, msg, 20, 34, w - 40, 20); updateKV(null, msg); };
  if (err) return fail(err);
  let md; try { md = makeModel(P); } catch (e) { return fail(e.message); }
  const L = md.Lat, spec = md.spec, C = spec.centre;
  const corner = (u, v) => [C[0] + u * L.a1[0] + v * L.a2[0], C[1] + u * L.a1[1] + v * L.a2[1]];
  let reg;
  if (Math.abs(L.a1[1]) < 1e-12) { const W = 1.5 * Math.abs(L.a1[0]), H = 1.5 * Math.abs(L.a2[1]); reg = [[C[0] - W, C[1] - H], [C[0] + W, C[1] - H], [C[0] + W, C[1] + H], [C[0] - W, C[1] + H]]; }
  else reg = [[-1.5, -1.5], [1.5, -1.5], [1.5, 1.5], [-1.5, 1.5]].map(c => corner(c[0], c[1]));
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity; for (const p of reg) { x0 = Math.min(x0, p[0]); x1 = Math.max(x1, p[0]); y0 = Math.min(y0, p[1]); y1 = Math.max(y1, p[1]); }
  const pad = 22, inset = Math.min(170, w * 0.3), sc = Math.min((w - 2 * pad - inset) / (x1 - x0), (h - 2 * pad) / (y1 - y0));
  const ox = pad + ((w - inset - 2 * pad) - sc * (x1 - x0)) / 2 - sc * x0, oy = h - pad - ((h - 2 * pad) - sc * (y1 - y0)) / 2 + sc * y0;
  const X = x => ox + sc * x, Y = y => oy - sc * y;
  const host = css('--host'), res = css('--res'), ink2 = css('--ink-2');
  const wid = sec => Math.max(1.4, (sec.shape === 'circ' ? sec.h : sec.b) * sc);
  ctx.save(); ctx.beginPath(); reg.forEach((p, i) => i ? ctx.lineTo(X(p[0]), Y(p[1])) : ctx.moveTo(X(p[0]), Y(p[1]))); ctx.closePath(); ctx.clip();
  const tiles = []; for (let i = -4; i <= 4; i++) for (let j = -3; j <= 3; j++) tiles.push([i, j]);
  const shift = (i, j) => [i * L.a1[0] + j * L.a2[0], i * L.a1[1] + j * L.a2[1]];
  // true-scale outlines: every beam is a rectangle of its real width; wall joints are filled to the walls' outer edges
  const panel = css('--panel'), hostL = mix(host, panel, 0.45), resL = mix(res, panel, 0.45);
  const halfW = ch => (ch.sec.shape === 'circ' ? ch.sec.h : ch.sec.b) / 2;
  const rectOf = (ch, Rx, Ry) => { const w2 = halfW(ch), nx = -ch.s * w2, ny = ch.c * w2, a = [ch.p0[0] + Rx, ch.p0[1] + Ry], b = [ch.p1[0] + Rx, ch.p1[1] + Ry];
    return [[a[0] + nx, a[1] + ny], [b[0] + nx, b[1] + ny], [b[0] - nx, b[1] - ny], [a[0] - nx, a[1] - ny]]; };
  const poly = (pts, col) => { ctx.fillStyle = col; ctx.strokeStyle = col; ctx.lineWidth = 0.7; ctx.beginPath(); pts.forEach((p, k) => k ? ctx.lineTo(X(p[0]), Y(p[1])) : ctx.moveTo(X(p[0]), Y(p[1]))); ctx.closePath(); ctx.fill(); ctx.stroke(); };
  const hull = Q => { const P2 = Q.slice().sort((a, b) => a[0] - b[0] || a[1] - b[1]), cr = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]), lo = [], up = [];
    for (const p of P2) { while (lo.length >= 2 && cr(lo[lo.length - 2], lo[lo.length - 1], p) <= 0) lo.pop(); lo.push(p); }
    for (const p of P2.slice().reverse()) { while (up.length >= 2 && cr(up[up.length - 2], up[up.length - 1], p) <= 0) up.pop(); up.push(p); }
    return lo.slice(0, -1).concat(up.slice(0, -1)); };
  const walls = md.chains.filter(ch => ch.role === 'host' && ch.beam.sec !== 'stub'), joints = [];
  spec.nodes.forEach((nd, idx) => { if (nd.role !== 'host') return; const pts = [];
    for (const ch of walls) { const r4 = rectOf(ch, 0, 0);
      if (ch.beam.i === idx) pts.push(r4[0], r4[3]);
      if (ch.beam.j === idx) { const Rn = shift(ch.beam.n[0], ch.beam.n[1]); pts.push([r4[1][0] - Rn[0], r4[1][1] - Rn[1]], [r4[2][0] - Rn[0], r4[2][1] - Rn[1]]); } }
    if (pts.length >= 6) joints.push(hull(pts)); });
  const drawOrder = tiles.filter(([i, j]) => i || j).concat([[0, 0]]);
  for (const [i, j] of drawOrder) for (const pass of ['res', 'host']) {
    const [Rx, Ry] = shift(i, j), main = i === 0 && j === 0;
    for (const ch of md.chains) { if (ch.role !== pass || ch.beam.sec === 'stub') continue; poly(rectOf(ch, Rx, Ry), pass === 'res' ? (main ? res : resL) : (main ? host : hostL)); }
    if (pass === 'host') for (const jt of joints) poly(jt.map(p => [p[0] + Rx, p[1] + Ry]), main ? host : hostL);
  }
  for (const [i, j] of tiles) {
    const [Rx, Ry] = shift(i, j); ctx.globalAlpha = (i === 0 && j === 0) ? 1 : 0.6;
    ctx.fillStyle = host;
    ctx.fillStyle = res; ctx.strokeStyle = res;
    for (const mn of md.massNodes) { const p = spec.nodes[mn.node], s2 = Math.max(6, mn.side * sc); ctx.fillRect(X(p.x + Rx) - s2 / 2, Y(p.y + Ry) - s2 / 2, s2, s2); }
    for (const sp of md.springs) { const nd = spec.nodes[sp.node], px = X(nd.x + Rx), py = Y(nd.y + Ry), r = 6 + 2.5 * Math.cbrt(sp.m * 1000);
      ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(px, py, r, 0, TAU); ctx.stroke(); ctx.beginPath(); ctx.arc(px, py, r * 0.42, 0, TAU); ctx.fill(); }
  }
  ctx.restore(); ctx.globalAlpha = 1;
  ctx.save(); ctx.setLineDash([5, 4]); ctx.strokeStyle = css('--ink'); ctx.lineWidth = 1.2; ctx.beginPath();
  [[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]].forEach((c, i) => { const p = corner(c[0], c[1]); i ? ctx.lineTo(X(p[0]), Y(p[1])) : ctx.moveTo(X(p[0]), Y(p[1])); });
  ctx.closePath(); ctx.stroke(); ctx.restore();
  // wall numbers
  ctx.font = '600 12px system-ui, sans-serif'; ctx.textAlign = 'center';
  for (const wl of (spec.walls || [])) { const mx = X((wl.p0[0] + wl.p1[0]) / 2), my = Y((wl.p0[1] + wl.p1[1]) / 2), dx = wl.p1[0] - wl.p0[0], dy = wl.p1[1] - wl.p0[1], ln = Math.hypot(dx, dy) || 1, ofx = -dy / ln * 14, ofy = -dx / ln * 14;
    ctx.fillStyle = css('--ink'); ctx.beginPath(); ctx.arc(mx + ofx, my + ofy, 9, 0, TAU); ctx.fill(); ctx.fillStyle = css('--panel'); ctx.fillText(wl.wall, mx + ofx, my + ofy + 4); }
  ctx.textAlign = 'left';
  // lattice vectors
  const o = corner(-0.5, -0.5); ctx.strokeStyle = ink2; ctx.fillStyle = ink2; ctx.lineWidth = 1.3; ctx.font = 'italic 14px Georgia, serif';
  [[L.a1, 'a₁'], [L.a2, 'a₂']].forEach(([v, lab]) => { const xa = X(o[0]), ya = Y(o[1]), xb = X(o[0] + v[0] * 0.92), yb = Y(o[1] + v[1] * 0.92); ctx.beginPath(); ctx.moveTo(xa, ya); ctx.lineTo(xb, yb); ctx.stroke(); const an = Math.atan2(yb - ya, xb - xa); ctx.beginPath(); ctx.moveTo(xb, yb); ctx.lineTo(xb - 8 * Math.cos(an - 0.4), yb - 8 * Math.sin(an - 0.4)); ctx.lineTo(xb - 8 * Math.cos(an + 0.4), yb - 8 * Math.sin(an + 0.4)); ctx.fill(); ctx.fillStyle = css('--panel'); const tx = xb + 8 * Math.cos(an - Math.PI / 2) + 4, ty = yb + 8 * Math.sin(an - Math.PI / 2) - 2; ctx.fillRect(tx - 2, ty - 13, 20, 17); ctx.fillStyle = ink2; ctx.fillText(lab, tx, ty); });
  // Brillouin zone inset
  const bx = w - inset / 2 - 8, by = 26 + inset / 2, br = inset * 0.4, rmax = Math.max(...L.bz.map(p => Math.hypot(p[0], p[1])));
  const K = p => [bx + br * p[0] / rmax, by - br * p[1] / rmax];
  ctx.fillStyle = css('--panel-2'); ctx.strokeStyle = css('--rule'); ctx.lineWidth = 1; ctx.beginPath(); L.bz.forEach((p, i) => { const q = K(p); i ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1]); }); ctx.closePath(); ctx.fill(); ctx.stroke();
  ctx.strokeStyle = css('--gap'); ctx.lineWidth = 2.5; ctx.beginPath(); L.path.pts.forEach((p, i) => { const q = K(p); i ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1]); }); ctx.stroke();
  ctx.fillStyle = css('--ink'); ctx.font = '600 12px system-ui, sans-serif';
  L.path.pts.slice(0, -1).forEach((p, i) => { const q = K(p); ctx.beginPath(); ctx.arc(q[0], q[1], 2.5, 0, TAU); ctx.fill(); const ux = p[0] / (Math.hypot(p[0], p[1]) || 1), uy = p[1] / (Math.hypot(p[0], p[1]) || 1); ctx.fillText(L.path.labels[i], q[0] + (i === 0 ? -13 : 6 * ux + 2), q[1] + (i === 0 ? 14 : -6 * uy + 4)); });
  ctx.fillStyle = ink2; ctx.font = '13px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.fillText('Brillouin zone', bx, by + br + 26); ctx.textAlign = 'left';
  updateKV(md, null, P);
}
function wrapText(ctx, txt, x, y, maxW, lh) { const words = txt.split(' '); let line = ''; for (const wd of words) { const t = line ? line + ' ' + wd : wd; if (ctx.measureText(t).width > maxW && line) { ctx.fillText(line, x, y); y += lh; line = wd; } else line = t; } ctx.fillText(line, x, y); }
function infoOnly() { if (isDrawn()) return; try { const P = readParams(), err = validate(P); if (err) return updateKV(null, err); updateKV(makeModel(P), null, P); } catch (e) { updateKV(null, e.message); } }
function updateKV(md, err, P) {
  const kv = $('#kv');
  if (!md) { kv.innerHTML = `<tr><td class="err" colspan="2">${err || 'Check the parameters.'}</td></tr>`; $('#latNote').textContent = ''; $('#resNote').textContent = ''; return; }
  if (!P) P = readParams();
  const { fr, mu, approx } = resInfo(P, md), L = md.Lat, spec = md.spec;
  const wl = [...new Set((spec.walls || []).map(w => Math.hypot(w.p1[0] - w.p0[0], w.p1[1] - w.p0[1])).map(v => (v * 1000).toFixed(2)))];
  const rows = [['Wall lengths', wl.length ? wl.join(', ') + ' mm' : '–'], ['Walls per cell', (spec.walls || []).length || '–'], ['Cell area', `${(L.area * 1e6).toFixed(0)} mm²`],
    ['Host mass per cell', fmtM(md.hostMass)], ['Resonator mass per cell', md.resMass ? fmtM(md.resMass) : '–'], ['Mass ratio μ', md.resMass ? mu.toFixed(3) : '–'],
    ['Areal density', `${((md.hostMass + md.resMass) / L.area).toFixed(2)} kg/m²`]];
  if (P.lattice === 'hexagon') {
    const G = P.G, regular = Math.abs(G.h - G.l) < 1e-9 * G.l && Math.abs(G.theta - 30) < 1e-6;
    rows.unshift(['Cell type', G.theta < 0 ? 'Re-entrant (auxetic)' : regular ? 'Regular hexagon' : 'Stretched hexagon'], ['Lattice constant |a₁|', `${(Math.hypot(L.a1[0], L.a1[1]) * 1000).toFixed(2)} mm`]);
    const bw = P.sections.host.shape === 'rect' ? P.sections.host.b : P.sections.host.h;
    rows.splice(3, 0, ['In-plane solidity of walls', `${(100 * bw * (G.h * P.vf + 2 * G.l) / L.area).toFixed(1)} %`]);
  }
  if (fr) { const fi = resInfo(P, md, 'ip').fr;
    rows.push([`Resonator fᵣ, out-of-plane${approx ? ' (approx.)' : ''}`, fmtF(fr)], ['Lumped gap estimate, out-of-plane', `${fmtF(fr)} to ${fmtF(fr * Math.sqrt(1 + mu))}`]);
    if (fi) rows.push([`Resonator fᵣ, in-plane${approx ? ' (approx.)' : ''}`, fmtF(fi)], ['Lumped gap estimate, in-plane', `${fmtF(fi)} to ${fmtF(fi * Math.sqrt(1 + mu))}`]); }
  const sy = L.sym, symTxt = [sy.r60 ? '6-fold' : sy.r90 ? '4-fold' : null, (sy.mx || sy.my) ? `mirror ${[sy.mx ? 'x' : '', sy.my ? 'y' : ''].filter(Boolean).join(' and ')}` : null].filter(Boolean).join(', ') || 'none found';
  rows.push(['Cell symmetry', symTxt], ['Brillouin-zone path', L.path.labels.join('–')], ['Degrees of freedom', `${md.n} (complex)`]);
  kv.innerHTML = rows.map(r => `<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`).join('');
  $('#latNote').textContent = P.lattice === 'custom' ? `${spec.nodes.length} nodes, ${spec.beams.length} beams.` : `${(spec.walls || []).length} walls per cell; numbers on the drawing.`;
  const fiN = fr ? resInfo(P, md, 'ip').fr : null;
  $('#resNote').innerHTML = fr ? `Isolated resonance ${approx ? '≈ ' : ''}<b>${fmtF(fr)}</b> out-of-plane${fiN ? `, <b>${fmtF(fiN)}</b> in-plane` : ''}; mass ratio <b>${mu.toFixed(2)}</b>.` : (md.resMass ? `Mass ratio <b>${mu.toFixed(2)}</b>.` : '');
  const jn = $('#jsonNote'); if (P.lattice === 'custom') jn.innerHTML = `Valid: ${spec.nodes.length} nodes, ${spec.beams.length} beams. Lengths in mm, masses in g, J in g·mm². “shift” is the cell offset [n₁, n₂] of node j; role “res” marks resonator parts.`;
}

// ---------- analysis ----------
async function bandsFor(md, path, nb0, target, label, p0, p1) {
  const nbMax = Math.min(60, md.n - 1); let nb = Math.min(nb0, nbMax), F, R, topMin;
  for (;;) {
    F = []; R = [];
    for (let i = 0; i < path.ks.length; i++) {
      const r = Eng.solveK(md, path.ks[i][0], path.ks[i][1], nb + 1, false, true); F.push(r.f); R.push(r.frac);
      if (i % 3 === 0) { setStatus(`${label} Bloch bands (${nb}): k-point ${i + 1} of ${path.ks.length}`, p0 + (p1 - p0) * i / path.ks.length); await tick(); }
    }
    topMin = Math.min(...F.map(r => r[r.length - 1]));
    if (!(target > topMin) || nb >= nbMax) break;
    nb = Math.min(nbMax, Math.ceil(nb * 1.6) + 2);
  }
  return { F, R, nb, topMin };
}
function intersectGaps(A, B) {
  const out = [];
  for (const a of A) for (const b of B) { const lo = Math.max(a.lo, b.lo), hi = Math.min(a.hi, b.hi); if (hi > lo * 1.001) out.push({ lo, hi }); }
  return out.sort((x, y) => x.lo - y.lo);
}
function bandRanges(x) { const nb = x.F[0].length, out = []; for (let b = 0; b < nb; b++) { let lo = Infinity, hi = -Infinity; for (const r of x.F) { lo = Math.min(lo, r[b]); hi = Math.max(hi, r[b]); } out.push({ b: b + 1, lo, hi }); } return out; }
async function runAnalysis(switchTab = true) {
  if (S.running) return;
  if (isDrawn()) { if (typeof DrawUI !== 'undefined') { const ok = await DrawUI.run(); if (ok && switchTab) showTab('disp'); } return; }
  const P = readParams(), err = validate(P);
  if (err) { setStatus(err); return; }
  startRun();
  const t0 = performance.now();
  try {
    const fam = {}, user = P.solver.fmax > 0 ? P.solver.fmax : 0, fmin = P.solver.fmin;
    // pass 1: each polarisation with the requested number of bands, to find its natural frequency range
    for (const [j, f] of FAMS.entries()) {
      const md = makeModel(P, undefined, undefined, f), path = Eng.kPath(md.Lat, P.solver.kper);
      const bb = await bandsFor(md, path, P.solver.nb, user, FAM_LABEL[f], 0.3 * j, 0.3 * j + 0.3);
      fam[f] = Object.assign({ md, path, family: f }, bb, { cap: freqCap(P, md, bb.topMin, f) });
    }
    // pass 2: both polarisations must reach the shared range used by the combined view
    const shared = user || Math.max(fam.oop.cap, fam.ip.cap);
    for (const f of FAMS) { const x = fam[f]; if (x.topMin < shared && x.nb < Math.min(60, x.md.n - 1)) Object.assign(x, await bandsFor(x.md, x.path, Math.ceil(x.nb * 1.6) + 2, shared, FAM_LABEL[f], 0.6, 0.7)); }
    if (fmin >= Math.min(fam.oop.cap, fam.ip.cap, shared)) throw new Error(`Min frequency must be below the max frequency (${fmtF(Math.min(fam.oop.cap, fam.ip.cap))}).`);
    for (const f of FAMS) {
      const x = fam[f], fmax = user || x.cap;
      x.fmax = fmax;
      x.gapsAll = Eng.findGaps(x.F, x.F[0].length).filter(g => g.lo < shared && g.hi > fmin);
      x.gaps = x.gapsAll.filter(g => g.lo < fmax);
      x.band = { F: x.F.map(r => r.slice(0, x.nb)), R: x.R.map(r => r.slice(0, x.nb)), path: x.path, fmax, fmin, nb: x.nb, incomplete: fmax > x.topMin * 1.001 ? x.topMin : null, family: f };
    }
    // complex band of each polarisation over its own range
    for (const [j, f] of FAMS.entries()) {
      const x = fam[f], W = Eng.buildWFE(x.md), rows = [];
      for (let i = 1; i <= P.solver.nf; i++) {
        const fr = fmin + (x.fmax - fmin) * i / P.solver.nf; rows.push({ f: fr, ev: Eng.wfeAt(W, fr) });
        if (i % 20 === 0) { setStatus(`${FAM_LABEL[f]} complex band: ${i} of ${P.solver.nf} frequencies`, 0.7 + 0.15 * j + 0.15 * i / P.solver.nf); await tick(); }
      }
      x.cplx = { rows, dir: x.md.Lat.std ? 'Γ–X' : 'b₁ (normal to a₂)' };
    }
    S.fam = fam; S.shared = { fmin, fmax: shared };
    S.complete = intersectGaps(fam.oop.gapsAll, fam.ip.gapsAll).filter(g => g.lo < shared && g.hi > fmin);
    S.P = P; S.sel = null; S.mode = null; S.hoverF = null;
    applyView();
    $('#modeSub').textContent = 'Click a point on the band structure'; $('#readout').textContent = 'Hover over either plot to read frequency, phase and attenuation.';
    const inc = FAMS.filter(f => fam[f].band.incomplete);
    setStatus(`Done in ${((performance.now() - t0) / 1000).toFixed(1)} s: out-of-plane ${fam.oop.gaps.length} gap${fam.oop.gaps.length === 1 ? '' : 's'}, in-plane ${fam.ip.gaps.length}, full gap (both) ${S.complete.length}${inc.length ? `; ${inc.map(f => FAM_LABEL[f].toLowerCase()).join(' and ')} bands near the top may be missing, raise “Bands to compute”` : ''}`);
    renderGapList(); S.bandKey = JSON.stringify(P); renderAssessment(); if (switchTab) showTab('disp'); else redrawAll();
  } catch (e) { setStatus('Solver stopped: ' + e.message); if (!/stopped by the user/.test(e.message)) console.error(e); }
  endRun();
}
function applyView() {
  const F = S.fam && S.fam[S.view === 'ip' ? 'ip' : 'oop']; if (!F) return;
  S.md = F.md; S.band = F.band; S.cplx = F.cplx; S.gaps = S.view === 'both' ? S.complete : F.gaps;
}
const VIEW_TEXT = {
  oop: ['Flexural (out-of-plane) dispersion', 'Left: real Bloch bands along the zone boundary path, coloured by how much kinetic energy sits in the resonators. Right: the complex band along the first path leg, phase change per cell and attenuation per cell against the same frequency axis. Hover to read a frequency; click a band point to animate its mode.'],
  ip: ['In-plane dispersion', 'In-plane waves: every node moves in x and y and rotates about z, so the beams stretch and bend in the plane of the lattice. Left: Bloch bands coloured by resonator energy share. Right: complex band along the first path leg. Click a band point to animate the in-plane mode, seen from above.'],
  both: ['Both polarisations: full band gap', 'Out-of-plane and in-plane bands on one frequency axis. A full band gap, hatched, is a range where neither polarisation propagates in any direction on the path; only there does the lattice block vibration of every kind. Right: where each polarisation propagates and where it is stopped.']
};
const LEGEND = {
  single: `<span>Host<span class="ramp"></span>Resonator energy share</span><span><span class="hatch"></span>Complete band gap (all directions on path)</span><span><span class="sw" style="background:var(--host)"></span>Propagating wave</span><span><span class="sw" style="background:var(--evan)"></span>Evanescent wave</span>`,
  both: `<span><span class="sw" style="background:var(--host)"></span>Out-of-plane bands</span><span><span class="sw" style="background:var(--edge);border-radius:1px"></span>In-plane bands</span><span><span class="hatch"></span>Full band gap (both polarisations)</span>`
};
function setView(v) {
  S.view = v; S.sel = null; S.mode = null; S.hoverF = null;
  document.querySelectorAll('.seg [data-view]').forEach(b => b.setAttribute('aria-checked', b.dataset.view === v));
  $('#dispTitle').textContent = VIEW_TEXT[v][0]; $('#dispLede').textContent = VIEW_TEXT[v][1];
  $('#dispLegend').innerHTML = v === 'both' ? LEGEND.both : LEGEND.single;
  $('#cplxTitle').textContent = v === 'both' ? 'Where each polarisation is stopped' : 'Complex band';
  $('#modeSub').textContent = 'Click a point on the band structure';
  $('#readout').textContent = v === 'both' ? 'Hover over either plot to see which polarisations propagate at that frequency.' : 'Hover over either plot to read frequency, phase and attenuation.';
  applyView(); renderGapList(); redrawAll();
}
document.querySelectorAll('.seg [data-view]').forEach(b => b.addEventListener('click', () => setView(b.dataset.view)));
$('#run').addEventListener('click', () => runAnalysis(true));

function renderGapList() {
  const el = $('#gaps'); if (!S.fam) return;
  const sh = S.shared, rows = [];
  for (const f of FAMS) { const R = S.fam[f].band.R; Aids.tagGaps(S.fam[f].gaps, g => g.below ? R.reduce((a, r) => a + (r[g.below - 1] || 0), 0) / R.length : null); }
  const add = (lab, g, full) => rows.push(`<tr><td>${lab}</td><td>${(!full && Aids.gapShort(g)) || '–'}</td><td>${g.below ? `${g.below} and ${g.below + 1}` : '–'}</td><td>${fmtF(g.lo)}</td><td>${fmtF(g.hi)}</td><td>${(200 * (g.hi - g.lo) / (g.hi + g.lo)).toFixed(1)} %</td></tr>`);
  const showOop = S.view !== 'ip', showIp = S.view !== 'oop';
  if (showOop) S.fam.oop.gaps.forEach(g => add('Out-of-plane', g));
  if (showIp) S.fam.ip.gaps.forEach(g => add('In-plane', g));
  if (S.view === 'both') S.complete.forEach(g => add('<b>Full (both)</b>', g, true));
  const fams = S.view === 'both' ? FAMS : [S.view];
  $('#gapSub').textContent = fams.map(f => `${FAM_LABEL[f].toLowerCase()} ${S.fam[f].nb} bands`).join(', ') + `, ${S.fam.oop.path.ks.length} k-points`;
  const { fr } = resInfo(S.P, S.fam.oop.md, 'oop'), fi = resInfo(S.P, S.fam.oop.md, 'ip').fr;
  const range = S.view === 'both' ? sh : S.fam[S.view];
  el.innerHTML = (rows.length ? `<table><thead><tr><th>Polarisation</th><th>Opened by</th><th>Between bands</th><th>From</th><th>To</th><th>Relative width</th></tr></thead><tbody>${rows.join('')}</tbody></table>`
    : `<p class="empty">No complete gap between ${fmtF(range.fmin)} and ${fmtF(range.fmax)} along the ${S.fam.oop.path.labels.join('–')} path. Try a larger mass ratio, or check the complex band for a directional gap.</p>`) +
    (fr ? `<p class="note" style="padding:0 12px 10px">Isolated resonator frequency ${fmtF(fr)} out-of-plane${fi ? `, ${fmtF(fi)} in-plane` : ''}. A gap that opens just above it, with flat resonator-coloured bands at its lower edge, is a local-resonance gap; gaps near the zone-boundary frequency of the bare beams are Bragg gaps.</p>` : '');
}

// ---------- band plot ----------
const M = { l: 58, r: 14, t: 16, b: 36 };
function drawBand() {
  if (S.view === 'both' && S.fam) return drawBandBoth();
  const cv = $('#cvBand'); const { ctx, w, h } = prep(cv); const B = S.band;
  const ink2 = css('--ink-2'), grid = css('--grid');
  if (!B) { ctx.fillStyle = ink2; ctx.font = '15px system-ui, sans-serif'; ctx.fillText('Press “Run analysis” to compute the bands.', 20, 40); return; }
  const xs = B.path.xs, xmax = xs[xs.length - 1], [fmin, fmax] = Aids.fz('main', B.fmin || 0, B.fmax), U = fUnit(fmax);
  const X = x => M.l + x / xmax * (w - M.l - M.r), Y = f => h - M.b - (f - fmin) / (fmax - fmin) * (h - M.t - M.b);
  S.bandGeom = { X, Y, w, h };
  for (const g of S.gaps) { if (g.hi < fmin || g.lo > fmax) continue; const y1 = Y(Math.min(g.hi, fmax)), y0 = Y(Math.max(g.lo, fmin)); ctx.fillStyle = css('--gap-tint'); ctx.fillRect(M.l, y1, w - M.l - M.r, y0 - y1); ctx.fillStyle = hatch(ctx, css('--gap')); ctx.globalAlpha = 0.45; ctx.fillRect(M.l, y1, w - M.l - M.r, y0 - y1); ctx.globalAlpha = 1; }
  ctx.font = '13px system-ui, sans-serif'; ctx.fillStyle = ink2; ctx.strokeStyle = grid; ctx.lineWidth = 1;
  for (const v of niceTicks(fmin / U.d, fmax / U.d, 6)) { const y = Y(v * U.d); ctx.beginPath(); ctx.moveTo(M.l, y); ctx.lineTo(w - M.r, y); ctx.stroke(); ctx.textAlign = 'right'; ctx.fillText(+v.toFixed(3), M.l - 7, y + 4); }
  ctx.textAlign = 'center';
  B.path.ticks.forEach((t, i) => { const x = X(t); ctx.strokeStyle = css('--rule'); ctx.beginPath(); ctx.moveTo(x, M.t); ctx.lineTo(x, h - M.b); ctx.stroke(); ctx.fillStyle = css('--ink'); ctx.font = '600 14px system-ui, sans-serif'; ctx.fillText(B.path.labels[i], x, h - M.b + 18); });
  ctx.save(); ctx.translate(15, (M.t + h - M.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillStyle = ink2; ctx.font = '13px system-ui, sans-serif'; ctx.fillText(`Frequency (${U.u})`, 0, 0); ctx.restore();
  ctx.fillStyle = ink2; ctx.font = '13px system-ui, sans-serif'; ctx.fillText('Wave vector along the zone path', (M.l + w - M.r) / 2, h - 4);
  const host = css('--host'), res = css('--res');
  ctx.save(); ctx.beginPath(); ctx.rect(M.l, M.t - 4, w - M.l - M.r, h - M.t - M.b + 4); ctx.clip();
  ctx.globalAlpha = 0.28; ctx.strokeStyle = ink2; ctx.lineWidth = 1;
  for (let b = 0; b < B.nb; b++) { ctx.beginPath(); xs.forEach((x, i) => { const y = Y(B.F[i][b]); i ? ctx.lineTo(X(x), y) : ctx.moveTo(X(x), y); }); ctx.stroke(); }
  ctx.globalAlpha = 1;
  Aids.drawBranches(ctx, X, Y, xs, B.F, B.path.ks, null, (B.family || 'oop') === 'ip' ? 2 : 1, { fmax });
  for (let b = 0; b < B.nb; b++) xs.forEach((x, i) => { const f = B.F[i][b]; if (f > fmax * 1.02) return; ctx.fillStyle = mix(host, res, B.R[i][b] || 0); ctx.beginPath(); ctx.arc(X(x), Y(f), 2.6, 0, TAU); ctx.fill(); });
  Aids.tagGaps(S.gaps, g => g.below ? B.R.reduce((a, r) => a + (r[g.below - 1] || 0), 0) / B.R.length : null); Aids.drawGapLabels(ctx, S.gaps, Y, w - M.r, fmin, fmax);
  if (S.sel) { const f = B.F[S.sel.ki][S.sel.b]; ctx.strokeStyle = css('--ink'); ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(X(xs[S.sel.ki]), Y(f), 7, 0, TAU); ctx.stroke(); }
  if (S.hoverF != null) hoverLine(ctx, M.l, w - M.r, Y(S.hoverF));
  ctx.restore();
  $('#bandSub').textContent = `${S.P.lattice[0].toUpperCase() + S.P.lattice.slice(1)} cell, ${FAM_LABEL[B.family || 'oop'].toLowerCase()}, ${B.nb} bands`;
}
function drawBandBoth() {
  const cv = $('#cvBand'); const { ctx, w, h } = prep(cv);
  const ink2 = css('--ink-2'), grid = css('--grid'), O = S.fam.oop.band, path = O.path, xs = path.xs, xmax = xs[xs.length - 1];
  const [fmin, fmax] = Aids.fz('main', S.shared.fmin, S.shared.fmax), U = fUnit(fmax);
  const X = x => M.l + x / xmax * (w - M.l - M.r), Y = f => h - M.b - (f - fmin) / (fmax - fmin) * (h - M.t - M.b);
  S.bandGeom = { X, Y, w, h, both: true };
  for (const g of S.complete) { if (g.hi < fmin || g.lo > fmax) continue; const y1 = Y(Math.min(g.hi, fmax)), y0 = Y(Math.max(g.lo, fmin)); ctx.fillStyle = css('--gap-tint'); ctx.fillRect(M.l, y1, w - M.l - M.r, y0 - y1); ctx.fillStyle = hatch(ctx, css('--gap')); ctx.globalAlpha = 0.45; ctx.fillRect(M.l, y1, w - M.l - M.r, y0 - y1); ctx.globalAlpha = 1; }
  ctx.font = '13px system-ui, sans-serif'; ctx.fillStyle = ink2; ctx.strokeStyle = grid; ctx.lineWidth = 1;
  for (const v of niceTicks(fmin / U.d, fmax / U.d, 6)) { const y = Y(v * U.d); ctx.beginPath(); ctx.moveTo(M.l, y); ctx.lineTo(w - M.r, y); ctx.stroke(); ctx.textAlign = 'right'; ctx.fillText(+v.toFixed(3), M.l - 7, y + 4); }
  ctx.textAlign = 'center';
  path.ticks.forEach((t, i) => { const x = X(t); ctx.strokeStyle = css('--rule'); ctx.beginPath(); ctx.moveTo(x, M.t); ctx.lineTo(x, h - M.b); ctx.stroke(); ctx.fillStyle = css('--ink'); ctx.font = '600 14px system-ui, sans-serif'; ctx.fillText(path.labels[i], x, h - M.b + 18); });
  ctx.save(); ctx.translate(15, (M.t + h - M.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillStyle = ink2; ctx.font = '13px system-ui, sans-serif'; ctx.fillText(`Frequency (${U.u})`, 0, 0); ctx.restore();
  ctx.fillStyle = ink2; ctx.font = '13px system-ui, sans-serif'; ctx.fillText('Wave vector along the zone path', (M.l + w - M.r) / 2, h - 4);
  ctx.save(); ctx.beginPath(); ctx.rect(M.l, M.t - 4, w - M.l - M.r, h - M.t - M.b + 4); ctx.clip();
  for (const [f, col] of [['ip', css('--edge')], ['oop', css('--host')]]) {
    const F = S.fam[f].F, nb = F[0].length;
    ctx.globalAlpha = 0.3; ctx.strokeStyle = col; ctx.lineWidth = 1;
    for (let b = 0; b < nb; b++) { ctx.beginPath(); xs.forEach((x, i) => { const y = Y(F[i][b]); i ? ctx.lineTo(X(x), y) : ctx.moveTo(X(x), y); }); ctx.stroke(); }
    ctx.globalAlpha = 1; Aids.drawBranches(ctx, X, Y, xs, F, path.ks, null, f === 'ip' ? 2 : 1, { fmax, tag: f === 'ip' ? 'in-plane' : 'out-of-plane' }); ctx.fillStyle = col;
    for (let b = 0; b < nb; b++) xs.forEach((x, i) => { const fq = F[i][b]; if (fq > fmax * 1.02) return; if (f === 'oop') { ctx.beginPath(); ctx.arc(X(x), Y(fq), 2.5, 0, TAU); ctx.fill(); } else ctx.fillRect(X(x) - 2.2, Y(fq) - 2.2, 4.4, 4.4); });
  }
  Aids.drawGapLabels(ctx, S.complete, Y, w - M.r, fmin, fmax, true);
  if (S.sel && S.sel.fam) { const fq = S.fam[S.sel.fam].F[S.sel.ki][S.sel.b]; ctx.strokeStyle = css('--ink'); ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(X(xs[S.sel.ki]), Y(fq), 7, 0, TAU); ctx.stroke(); }
  if (S.hoverF != null) hoverLine(ctx, M.l, w - M.r, Y(S.hoverF));
  ctx.restore();
  $('#bandSub').textContent = `${S.P.lattice[0].toUpperCase() + S.P.lattice.slice(1)} cell, out-of-plane ● and in-plane ■`;
}
function drawGapMap() {
  const cv = $('#cvCplx'); const { ctx, w, h } = prep(cv);
  const ink2 = css('--ink-2'), [fmin, fmax] = Aids.fz('main', S.shared.fmin, S.shared.fmax), U = fUnit(fmax);
  const Y = f => h - M.b - (f - fmin) / (fmax - fmin) * (h - M.t - M.b);
  S.cplxGeom = { Y, h, both: true };
  const l = 16, r = 12, cw = (w - l - r) / 3, top = Y(fmax), bot = Y(fmin);
  ctx.strokeStyle = css('--grid'); ctx.lineWidth = 1;
  for (const v of niceTicks(fmin / U.d, fmax / U.d, 6)) { const y = Y(v * U.d); ctx.beginPath(); ctx.moveTo(l, y); ctx.lineTo(w - r, y); ctx.stroke(); }
  const cols = [['Out-of-plane', ['oop'], css('--host'), S.fam.oop.gapsAll], ['In-plane', ['ip'], css('--edge'), S.fam.ip.gapsAll], ['Full gap', ['oop', 'ip'], css('--evan'), S.complete]];
  ctx.save(); ctx.beginPath(); ctx.rect(l, M.t - 2, w - l - r, h - M.t - M.b + 4); ctx.clip();
  cols.forEach(([lab, fams, col, gaps], i) => {
    const x0 = l + i * cw + 10, x1 = x0 + cw - 20;
    ctx.fillStyle = css('--panel-2'); ctx.fillRect(x0, top, x1 - x0, bot - top);
    ctx.globalAlpha = i < 2 ? 0.55 : 0.4;
    for (const f of fams) { ctx.fillStyle = i < 2 ? col : (f === 'oop' ? css('--host') : css('--edge')); for (const bnd of bandRanges(S.fam[f])) { if (bnd.lo > fmax) continue; const ya = Y(Math.min(bnd.hi, fmax)), yb = Y(Math.max(bnd.lo, fmin)); if (yb > ya) ctx.fillRect(x0, ya, x1 - x0, yb - ya); } }
    ctx.globalAlpha = 1;
    for (const g of gaps) { const ya = Y(Math.min(g.hi, fmax)), yb = Y(Math.max(g.lo, fmin)); if (yb <= ya) continue;
      ctx.fillStyle = css('--gap-tint'); ctx.fillRect(x0, ya, x1 - x0, yb - ya); ctx.fillStyle = hatch(ctx, css('--gap')); ctx.globalAlpha = 0.6; ctx.fillRect(x0, ya, x1 - x0, yb - ya); ctx.globalAlpha = 1;
      ctx.strokeStyle = css('--gap'); ctx.lineWidth = 1.2; ctx.strokeRect(x0 + 0.5, ya + 0.5, x1 - x0 - 1, yb - ya - 1);
      if (yb - ya > 26) { ctx.fillStyle = css('--ink'); ctx.font = '600 12px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.fillText(fmtF(g.lo), (x0 + x1) / 2, yb - 6); ctx.fillText(fmtF(g.hi), (x0 + x1) / 2, ya + 14); } }
    for (const f of fams) if (S.fam[f].topMin < fmax) { const ya = Y(fmax), yb = Y(S.fam[f].topMin); ctx.fillStyle = css('--paper'); ctx.fillRect(x0, ya, x1 - x0, yb - ya); }
  });
  if (S.hoverF != null) hoverLine(ctx, l, w - r, Y(S.hoverF));
  ctx.restore();
  ctx.fillStyle = ink2; ctx.font = '13px system-ui, sans-serif'; ctx.textAlign = 'center';
  cols.forEach(([lab], i) => ctx.fillText(lab, l + i * cw + cw / 2, h - M.b + 16));
  ctx.fillText('Shaded: propagating bands. Hatched: stopped.', (l + w - r) / 2, h - 4);
  $('#cplxSub').textContent = 'all directions on the zone path';
}
function hoverLine(ctx, x0, x1, y) { ctx.save(); ctx.setLineDash([4, 3]); ctx.strokeStyle = css('--ink'); ctx.globalAlpha = 0.6; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x1, y); ctx.stroke(); ctx.restore(); }

// ---------- complex band plot ----------
const ATT_MAX = 2.5, PROP = 0.015;
function drawCplx() {
  if (S.view === 'both' && S.fam) return drawGapMap();
  const cv = $('#cvCplx'); const { ctx, w, h } = prep(cv); const C = S.cplx, B = S.band;
  const ink2 = css('--ink-2'), grid = css('--grid');
  if (!C) { ctx.fillStyle = ink2; ctx.font = '15px system-ui, sans-serif'; ctx.fillText('Runs together with the bands.', 20, 40); return; }
  const [fmin, fmax] = Aids.fz('main', B.fmin || 0, B.fmax), Y = f => h - M.b - (f - fmin) / (fmax - fmin) * (h - M.t - M.b);
  const l = 16, r = 12, mid = l + (w - l - r) * 0.44;
  const XA = a => mid - (a / ATT_MAX) * (mid - l), XP = p => mid + 6 + (p / Math.PI) * (w - r - mid - 6);
  S.cplxGeom = { Y, h };
  // directional gaps (all waves evanescent)
  const evan = C.rows.map(row => row.ev ? Math.min(...row.ev.map(e => e.att)) > PROP : false);
  ctx.fillStyle = css('--gap-tint');
  let start = null; const df = (fmax - fmin) / C.rows.length;
  C.rows.forEach((row, i) => { if (evan[i] && start == null) start = row.f - df / 2; if ((!evan[i] || i === C.rows.length - 1) && start != null) { const end = evan[i] ? row.f + df / 2 : row.f - df / 2; ctx.fillRect(l, Y(end), w - l - r, Y(start) - Y(end)); start = null; } });
  for (const g of S.gaps) { if (g.hi < fmin || g.lo > fmax) continue; const y1 = Y(Math.min(g.hi, fmax)), y0 = Y(Math.max(g.lo, fmin)); ctx.fillStyle = hatch(ctx, css('--gap')); ctx.globalAlpha = 0.45; ctx.fillRect(l, y1, w - l - r, y0 - y1); ctx.globalAlpha = 1; }
  ctx.strokeStyle = grid; ctx.lineWidth = 1; ctx.font = '13px system-ui, sans-serif'; ctx.fillStyle = ink2; ctx.textAlign = 'center';
  const U = fUnit(fmax);
  for (const v of niceTicks(fmin / U.d, fmax / U.d, 6)) { const y = Y(v * U.d); ctx.beginPath(); ctx.moveTo(l, y); ctx.lineTo(w - r, y); ctx.stroke(); }
  for (const a of [0.5, 1, 1.5, 2, 2.5]) { const x = XA(a); ctx.beginPath(); ctx.moveTo(x, M.t); ctx.lineTo(x, h - M.b); ctx.stroke(); ctx.fillText(a, x, h - M.b + 16); }
  [[0, '0'], [Math.PI / 2, 'π/2'], [Math.PI, 'π']].forEach(([p, s]) => { const x = XP(p); ctx.beginPath(); ctx.moveTo(x, M.t); ctx.lineTo(x, h - M.b); ctx.stroke(); ctx.fillText(s, x, h - M.b + 16); });
  ctx.strokeStyle = css('--ink'); ctx.lineWidth = 1.2; ctx.beginPath(); ctx.moveTo(mid + 3, M.t); ctx.lineTo(mid + 3, h - M.b); ctx.stroke();
  ctx.fillStyle = ink2; ctx.fillText('Attenuation Im(μ), Np/cell', (l + mid) / 2, h - 4); ctx.fillText('Phase Re(μ), rad/cell', (mid + w - r) / 2, h - 4);
  ctx.save(); ctx.beginPath(); ctx.rect(l, M.t - 4, w - l - r, h - M.t - M.b + 4); ctx.clip();
  const host = css('--host'), ev = css('--evan');
  for (const row of C.rows) { if (!row.ev) continue; const y = Y(row.f);
    for (const e of row.ev) {
      if (e.att > ATT_MAX * 1.05) continue;
      const prop = e.att < PROP;
      ctx.fillStyle = prop ? host : ev; ctx.globalAlpha = prop ? 1 : 0.85;
      ctx.fillRect(XP(e.ph) - 1.4, y - 1.4, 2.8, 2.8);
      if (!prop) ctx.fillRect(XA(e.att) - 1.4, y - 1.4, 2.8, 2.8);
    } }
  ctx.globalAlpha = 1;
  if (S.hoverF != null) hoverLine(ctx, l, w - r, Y(S.hoverF));
  ctx.restore();
  $('#cplxSub').textContent = `μ = k·a₁ along ${C.dir}`;
}
function onHover(ev, which) {
  const g = which === 'band' ? S.bandGeom : S.cplxGeom; if (!g || !S.band) return;
  if (S.view === 'both') {
    const rect = ev.currentTarget.getBoundingClientRect(), y = ev.clientY - rect.top, [fmin, fmax] = Aids.fz('main', S.shared.fmin, S.shared.fmax);
    const f = fmin + (g.h - M.b - y) / (g.h - M.t - M.b) * (fmax - fmin);
    if (f < fmin || f > fmax) { S.hoverF = null; drawBand(); drawCplx(); return; }
    S.hoverF = f; drawBand(); drawCplx();
    const st = fam => S.fam[fam].gapsAll.some(q => f >= q.lo && f <= q.hi);
    const so = st('oop'), si = st('ip');
    $('#readout').textContent = `f = ${fmtF(f)}. Out-of-plane waves ${so ? 'are stopped' : 'propagate'}; in-plane waves ${si ? 'are stopped' : 'propagate'}.${so && si ? ' Inside a full band gap: no vibration of either polarisation propagates.' : ''}`;
    return;
  }
  const rect = ev.currentTarget.getBoundingClientRect(), y = ev.clientY - rect.top, [fmin, fmax] = Aids.fz('main', S.band.fmin || 0, S.band.fmax), F0 = S.band.fmin || 0, F1 = S.band.fmax;
  const f = fmin + (g.h - M.b - y) / (g.h - M.t - M.b) * (fmax - fmin);
  if (f < fmin || f > fmax) { S.hoverF = null; $('#readout').textContent = 'Hover over either plot to read frequency, phase and attenuation.'; drawBand(); drawCplx(); return; }
  S.hoverF = f; drawBand(); drawCplx();
  const rows = S.cplx.rows, i = Math.min(rows.length - 1, Math.max(0, Math.round((f - F0) / (F1 - F0) * rows.length) - 1)), row = rows[i];
  let txt = `f = ${fmtF(f)}`;
  if (row && row.ev) {
    const prop = row.ev.filter(e => e.att < PROP), minAtt = Math.min(...row.ev.map(e => e.att));
    txt += prop.length ? `. Propagating along ${S.cplx.dir}: ${[...new Set(prop.map(e => (e.ph / Math.PI).toFixed(2)))].join(', ')} π rad/cell` : `. No propagating wave along ${S.cplx.dir}; the least attenuated wave loses ${(8.686 * minAtt).toFixed(1)} dB per cell (${(8.686 * minAtt * 5).toFixed(0)} dB over 5 cells)`;
  }
  const inGap = S.gaps.find(gp => f >= gp.lo && f <= gp.hi); if (inGap) txt += '. Inside a complete band gap.';
  $('#readout').textContent = txt;
}
Aids.zoomable($('#cvBand'), 'main', () => { drawBand(); drawCplx(); }); Aids.zoomable($('#cvCplx'), 'main', () => { drawBand(); drawCplx(); });
['band', 'cplx'].forEach(k => { const cv = $(k === 'band' ? '#cvBand' : '#cvCplx'); cv.addEventListener('mousemove', e => onHover(e, k)); cv.addEventListener('mouseleave', () => { S.hoverF = null; drawBand(); drawCplx(); }); });
$('#cvBand').addEventListener('click', e => {
  const g = S.bandGeom, B = S.band; if (!g || !B) return;
  if (S.view === 'both') {
    const rect = e.currentTarget.getBoundingClientRect(), mx = e.clientX - rect.left, my = e.clientY - rect.top;
    let best = null, bd = 18 * 18;
    for (const f of FAMS) { const F = S.fam[f].F; B.path.xs.forEach((x, i) => { for (let b = 0; b < F[i].length; b++) { const d = (g.X(x) - mx) ** 2 + (g.Y(F[i][b]) - my) ** 2; if (d < bd) { bd = d; best = { ki: i, b, fam: f }; } } }); }
    if (!best) return;
    S.sel = best; const k = B.path.ks[best.ki], md = S.fam[best.fam].md;
    const r = Eng.solveK(md, k[0], k[1], best.b + 1, true, true);
    S.mode = { u: r.vecs[best.b], f: r.f[best.b], k, frac: r.frac[best.b], md };
    const brA = Aids.branchOf(S.fam[best.fam].F, B.path.ks, null, best.fam === 'ip' ? 2 : 1, best.ki, best.b);
    $('#modeSub').textContent = `${FAM_LABEL[best.fam]} band ${best.b + 1}${brA ? ' (' + brA + ')' : ''}, ${fmtF(S.mode.f)}, ${labelK(best.ki)}, ${Math.round(100 * S.mode.frac)} % of kinetic energy in resonators`;
    drawBand(); return;
  }
  const rect = e.currentTarget.getBoundingClientRect(), mx = e.clientX - rect.left, my = e.clientY - rect.top;
  let best = null, bd = 18 * 18;
  B.path.xs.forEach((x, i) => { for (let b = 0; b < B.nb; b++) { const d = (g.X(x) - mx) ** 2 + (g.Y(B.F[i][b]) - my) ** 2; if (d < bd) { bd = d; best = { ki: i, b }; } } });
  if (!best) return;
  S.sel = best; const k = B.path.ks[best.ki];
  const r = Eng.solveK(S.md, k[0], k[1], best.b + 1, true, true);
  S.mode = { u: r.vecs[best.b], f: r.f[best.b], k, frac: r.frac[best.b], md: S.md };
  const lab = labelK(best.ki);
  const brB = Aids.branchOf(B.F, B.path.ks, null, (B.family || 'oop') === 'ip' ? 2 : 1, best.ki, best.b);
  $('#modeSub').textContent = `${FAM_LABEL[S.md.family]} band ${best.b + 1}${brB ? ' (' + brB + ')' : ''}, ${fmtF(S.mode.f)}, ${lab}, ${Math.round(100 * S.mode.frac)} % of kinetic energy in resonators`;
  drawBand();
});
function labelK(ki) { const B = S.band, [sg, u] = B.path.seg[ki], L = B.path.labels; if (u < 1e-9) return L[sg]; if (u > 1 - 1e-9) return L[sg + 1]; return `${Math.round(100 * u)} % from ${L[sg]} to ${L[sg + 1]}`; }

// ---------- mode animation ----------
function drawMode(cvSel = '#cvMode', mode = S.mode, mdIn = (S.mode && S.mode.md) || S.md, cellsIn = null) {
  const cv = $(cvSel); const { ctx, w, h } = prep(cv);
  if (!mode || !mdIn) { ctx.fillStyle = css('--ink-2'); ctx.font = '15px system-ui, sans-serif'; ctx.fillText('Click a point on the band structure to see its mode shape.', 18, 36); return; }
  if (mdIn.family === 'ip') return drawModeIP(ctx, w, h, mode, mdIn, cellsIn);
  const md = mdIn, L = md.Lat, u = mode.u, k = mode.k, th = S.phase;
  const ca = Math.cos(-0.6), sa = Math.sin(-0.6);
  const proj = (x, y) => [x * ca - y * sa, (x * sa + y * ca) * 0.5];
  let wmax = 1e-30; for (let i = 0; i < md.N; i++) wmax = Math.max(wmax, Math.hypot(u.re[3 * i], u.im[3 * i]));
  for (const sp of md.springs) wmax = Math.max(wmax, Math.hypot(u.re[sp.dof], u.im[sp.dof]));
  const cells = cellsIn || Aids.tiles();
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (const [i, j] of cells) for (const nd of md.nodes) { const p = proj(nd.x + i * L.a1[0] + j * L.a2[0], nd.y + i * L.a1[1] + j * L.a2[1]); x0 = Math.min(x0, p[0]); x1 = Math.max(x1, p[0]); y0 = Math.min(y0, p[1]); y1 = Math.max(y1, p[1]); }
  const amp = 0.16 * (h) * parseFloat($('#amp').value), padY = amp + 12, sc = Math.min((w - 40) / (x1 - x0), (h - 2 * padY) / (y1 - y0));
  const ox = w / 2 - sc * (x0 + x1) / 2, oy = h / 2 + sc * (y0 + y1) / 2;
  const disp = (dofIdx, R) => { const ph = k[0] * R[0] + k[1] * R[1] + th; return (u.re[dofIdx] * Math.cos(ph) - u.im[dofIdx] * Math.sin(ph)) / wmax; };
  const scr = (x, y, wv) => { const p = proj(x, y); return [ox + sc * p[0], oy - sc * p[1] - amp * wv]; };
  const host = css('--host'), res = css('--res'), rule = css('--rule');
  // reference plane (cell 0)
  ctx.strokeStyle = rule; ctx.lineWidth = 1; ctx.setLineDash([3, 4]);
  for (const [i, j] of cells) for (const ch of md.chains) { const R = [i * L.a1[0] + j * L.a2[0], i * L.a1[1] + j * L.a2[1]]; const a = scr(ch.p0[0] + R[0], ch.p0[1] + R[1], 0), b = scr(ch.p1[0] + R[0], ch.p1[1] + R[1], 0); ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke(); }
  ctx.setLineDash([]);
  const segs = [];
  for (const [i, j] of cells) {
    const R = [i * L.a1[0] + j * L.a2[0], i * L.a1[1] + j * L.a2[1]];
    for (const e of md.elems) {
      const Rq = [R[0] + e.n[0] * L.a1[0] + e.n[1] * L.a2[0], R[1] + e.n[0] * L.a1[1] + e.n[1] * L.a2[1]];
      const p = md.nodes[e.p], q = md.nodes[e.q], wp = disp(3 * e.p, R), wq = disp(3 * e.q, Rq);
      const a = scr(p.x + R[0], p.y + R[1], wp), b = scr(q.x + Rq[0], q.y + Rq[1], wq);
      segs.push({ a, b, z: (a[1] + b[1]) / 2, res: e.res, wv: (wp + wq) / 2 });
    }
  }
  segs.sort((s1, s2) => s1.z - s2.z);
  for (const s of segs) { ctx.strokeStyle = s.res ? res : host; ctx.lineWidth = s.res ? 2.2 : 3; ctx.lineCap = 'round'; ctx.beginPath(); ctx.moveTo(s.a[0], s.a[1]); ctx.lineTo(s.b[0], s.b[1]); ctx.stroke(); }
  for (const [i, j] of cells) {
    const R = [i * L.a1[0] + j * L.a2[0], i * L.a1[1] + j * L.a2[1]];
    for (const c of md.massNodes) { const nd = md.nodes[c.node], wv = disp(3 * c.node, R), p = scr(nd.x + R[0], nd.y + R[1], wv); ctx.fillStyle = res; ctx.fillRect(p[0] - 4, p[1] - 4, 8, 8); }
    for (const sp of md.springs) {
      const nd = md.nodes[sp.node], wb = disp(3 * sp.node, R), wr = disp(sp.dof, R);
      const base = scr(nd.x + R[0], nd.y + R[1], wb), top = scr(nd.x + R[0], nd.y + R[1], wr); top[1] -= 22;
      ctx.strokeStyle = res; ctx.lineWidth = 1.4; ctx.beginPath(); ctx.moveTo(base[0], base[1]);
      const n = 6; for (let t = 1; t < n; t++) { const yy = base[1] + (top[1] - base[1]) * t / n; ctx.lineTo(base[0] + (t % 2 ? 4 : -4), yy); } ctx.lineTo(top[0], top[1]); ctx.stroke();
      ctx.fillStyle = res; ctx.beginPath(); ctx.arc(top[0], top[1] - 4, 5, 0, TAU); ctx.fill();
    }
  }
  // motion arrows: from the rest position to the moved position, up / down (out-of-plane)
  if (Aids.st.arrows) {
    const cellL = Math.hypot(L.a1[0], L.a1[1]), many = cells.length > 25, ak = (many ? 'm' : 'f') + md.N;
    if (!mode._arw || mode._arw.k !== ak) { const mags = md.nodes.map((_, i) => Math.hypot(u.re[3 * i], u.im[3 * i])); mode._arw = { k: ak, ids: Aids.pickNodes(md.nodes.map(n => [n.x, n.y]), mags, cellL / (many ? 2 : 4), md.massNodes.map(c => c.node)) }; }
    const col = Aids.colOOP();
    for (const [i, j] of cells) { const R = [i * L.a1[0] + j * L.a2[0], i * L.a1[1] + j * L.a2[1]];
      for (const id of mode._arw.ids) { const nd = md.nodes[id], wv = disp(3 * id, R), a = scr(nd.x + R[0], nd.y + R[1], 0), b = scr(nd.x + R[0], nd.y + R[1], wv); Aids.arrow(ctx, a[0], a[1], b[0], b[1], col, many ? 1.2 : 1.7); }
      for (const sp of md.springs) { const nd = md.nodes[sp.node], wr = disp(sp.dof, R), a = scr(nd.x + R[0], nd.y + R[1], 0), b = scr(nd.x + R[0], nd.y + R[1], wr); Aids.arrow(ctx, a[0] + 9, a[1] - 26, b[0] + 9, b[1] - 26, col, many ? 1.2 : 1.7); } }
    Aids.arrowKey(ctx, w, ['oop']);
  }
  ctx.fillStyle = css('--ink-2'); ctx.font = '12.5px system-ui, sans-serif'; ctx.textAlign = 'left'; ctx.fillText(`${cellsIn ? cells.length + ' cells' : Aids.label() + ' cells'}; dashed: plane at rest`, 12, h - 10);
}
function drawModeIP(ctx, w, h, mode, md, cellsIn) {
  const L = md.Lat, u = mode.u, k = mode.k, th = S.phase;
  const cells = cellsIn || Aids.tiles();
  let umax = 1e-30; for (let i = 0; i < md.N; i++) umax = Math.max(umax, Math.hypot(u.re[3 * i], u.im[3 * i], u.re[3 * i + 1], u.im[3 * i + 1]));
  for (const sp of md.springs) umax = Math.max(umax, Math.hypot(u.re[sp.dof], u.im[sp.dof]));
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (const [i, j] of cells) for (const nd of md.nodes) { const x = nd.x + i * L.a1[0] + j * L.a2[0], y = nd.y + i * L.a1[1] + j * L.a2[1]; x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, y); y1 = Math.max(y1, y); }
  const cell = Math.hypot(L.a1[0], L.a1[1]), pad = 0.12 * cell * parseFloat($('#amp').value) + 0.02 * cell;
  const sc = Math.min((w - 30) / (x1 - x0 + 2 * pad), (h - 30) / (y1 - y0 + 2 * pad)), ox = w / 2 - sc * (x0 + x1) / 2, oy = h / 2 + sc * (y0 + y1) / 2;
  const dsc = 0.12 * cell * parseFloat($('#amp').value) / umax;
  const val = (dof, R) => { const ph = k[0] * R[0] + k[1] * R[1] + th; return u.re[dof] * Math.cos(ph) - u.im[dof] * Math.sin(ph); };
  const P = (nd, R, deformed) => { const x = nd.x + R[0] + (deformed ? dsc * val(3 * nd.id, R) : 0), y = nd.y + R[1] + (deformed ? dsc * val(3 * nd.id + 1, R) : 0); return [ox + sc * x, oy - sc * y]; };
  md.nodes.forEach((nd, i) => nd.id = i);
  const host = css('--host'), res = css('--res');
  for (const pass of [false, true]) for (const [i, j] of cells) {
    const R = [i * L.a1[0] + j * L.a2[0], i * L.a1[1] + j * L.a2[1]];
    for (const e of md.elems) {
      const Rq = [R[0] + e.n[0] * L.a1[0] + e.n[1] * L.a2[0], R[1] + e.n[0] * L.a1[1] + e.n[1] * L.a2[1]];
      const a = P(md.nodes[e.p], R, pass), b = P(md.nodes[e.q], Rq, pass);
      ctx.strokeStyle = pass ? (e.res ? res : host) : css('--rule'); ctx.lineWidth = pass ? (e.res ? 2.2 : 2.6) : 1; ctx.setLineDash(pass ? [] : [3, 4]); ctx.lineCap = 'round';
      ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
    }
  }
  ctx.setLineDash([]);
  for (const [i, j] of cells) {
    const R = [i * L.a1[0] + j * L.a2[0], i * L.a1[1] + j * L.a2[1]];
    for (const c of md.massNodes) { const p = P(md.nodes[c.node], R, true); ctx.fillStyle = res; ctx.fillRect(p[0] - 4, p[1] - 4, 8, 8); }
    const byNode = {}; for (const sp of md.springs) (byNode[sp.node] = byNode[sp.node] || [])[sp.c] = sp.dof;
    for (const nid in byNode) { const nd = md.nodes[nid], d = byNode[nid], base = P(nd, R, true);
      const mx = nd.x + R[0] + dsc * val(d[0], R), my = nd.y + R[1] + dsc * val(d[1], R), q = [ox + sc * mx, oy - sc * my];
      ctx.strokeStyle = res; ctx.lineWidth = 1.2; ctx.beginPath(); ctx.moveTo(base[0], base[1]); ctx.lineTo(q[0], q[1]); ctx.stroke();
      ctx.fillStyle = res; ctx.beginPath(); ctx.arc(q[0], q[1], 4.5, 0, TAU); ctx.fill(); }
  }
  // motion arrows: from the rest position to the moved position, sideways (in-plane)
  if (Aids.st.arrows) {
    const many = cells.length > 25, ak = (many ? 'm' : 'f') + md.N;
    if (!mode._arw || mode._arw.k !== ak) { const mags = md.nodes.map((_, i) => Math.hypot(u.re[3 * i], u.im[3 * i], u.re[3 * i + 1], u.im[3 * i + 1])); mode._arw = { k: ak, ids: Aids.pickNodes(md.nodes.map(n => [n.x, n.y]), mags, cell / (many ? 2 : 4), md.massNodes.map(c => c.node)) }; }
    const col = Aids.colIP(), byNode = {}; for (const sp of md.springs) (byNode[sp.node] = byNode[sp.node] || [])[sp.c] = sp.dof;
    for (const [i, j] of cells) { const R = [i * L.a1[0] + j * L.a2[0], i * L.a1[1] + j * L.a2[1]];
      const G = 2.5, ar = (a, b) => Aids.arrow(ctx, a[0], a[1], a[0] + G * (b[0] - a[0]), a[1] + G * (b[1] - a[1]), col, many ? 1.3 : 1.8); // arrows drawn 2.5 × the plotted displacement
      for (const id of mode._arw.ids) ar(P(md.nodes[id], R, false), P(md.nodes[id], R, true));
      for (const nid in byNode) { const nd = md.nodes[nid], d = byNode[nid]; ar(P(nd, R, false), [ox + sc * (nd.x + R[0] + dsc * val(d[0], R)), oy - sc * (nd.y + R[1] + dsc * val(d[1], R))]); } }
    Aids.arrowKey(ctx, w, ['ip']);
  }
  ctx.fillStyle = css('--ink-2'); ctx.font = '12.5px system-ui, sans-serif'; ctx.textAlign = 'left'; ctx.fillText(`Top view, ${cellsIn ? cells.length + ' cells' : Aids.label() + ' cells'}; dashed: at rest`, 12, h - 10);
}
function animate() {
  if (S.playing && S.mode && !$('#p-disp').classList.contains('hidden')) { S.phase += 0.06; drawMode(); }
  if (S.playing && S.emode && !$('#p-edge').classList.contains('hidden')) { S.phase += 0.06; drawEdgeMode(); }
  requestAnimationFrame(animate);
}
$('#play').addEventListener('click', () => { S.playing = !S.playing; $('#play').textContent = S.playing ? 'Pause' : 'Play'; });
$('#amp').addEventListener('input', drawMode);
window.addEventListener('mmview', () => { if (S.mode) S.mode._arw = null; if (!$('#p-disp').classList.contains('hidden')) { drawBand(); drawMode(); } });
if (window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches) { S.playing = false; $('#play').textContent = 'Play'; }

// ---------- sensitivity sweep ----------
async function runSweep() {
  if (S.running) return;
  if (isDrawn()) { if (typeof DrawUI !== 'undefined') await DrawUI.runSweep(); return; }
  const P = readParams(), err = validate(P); if (err) { setStatus(err); return; }
  const key = $('#sp').value, rng = num('sr'); let n = Math.round(num('sn'));
  if (!(rng > 1)) { setStatus('Range must be greater than 1.'); return; }
  if (!(n >= 3)) n = 3; if (n % 2 === 0) n += 1; n = Math.min(n, 41);
  let spec0; try { spec0 = buildSpec(P); } catch (e) { setStatus(e.message); return; }
  const sfam = $('#sFam').value;
  const md0 = Eng.buildModel(spec0, P.mat, { ne: P.ne, path: P.path, family: sfam, beam: P.beam });
  if ((key === 'resM' || key === 'resK') && !md0.resMass) { setStatus('This design has no resonator, so there is no resonator mass or stiffness to vary.'); return; }
  startRun();
  const t0 = performance.now(), samples = [], kper = Math.max(6, Math.round(P.solver.kper / 2));
  let cap = P.solver.fmax > 0 ? P.solver.fmax : null;
  try {
    for (let s = 0; s < n; s++) {
      const fac = rng ** (2 * s / (n - 1) - 1), scale = { hostE: 1, hostRho: 1, resK: 1, resM: 1 };
      if (key !== 'a') scale[key] = fac;
      const md = Eng.buildModel(key === 'a' ? Eng.scaleSpec(spec0, fac) : spec0, P.mat, { ne: P.ne, path: P.path, scale, family: sfam, beam: P.beam }), nb = Math.min(P.solver.nb, md.n), path = Eng.kPath(md.Lat, kper), F = [];
      for (let i = 0; i < path.ks.length; i++) { F.push(Eng.solveK(md, path.ks[i][0], path.ks[i][1], nb + 1, false, false).f); if (i % 4 === 0) { setStatus(`Sweep: sample ${s + 1} of ${n}`, (s + i / path.ks.length) / n); await tick(); } }
      const top = Math.min(...F.map(r => r[r.length - 1]));
      const gaps = Eng.findGaps(F, F[0].length).filter(g => g.lo < top);
      const fr0 = resInfo(P, md0, sfam).fr; let fr = null;
      if (fr0) fr = key === 'resM' ? fr0 / Math.sqrt(fac) : key === 'resK' ? fr0 * Math.sqrt(fac) : key === 'a' ? null : fr0;
      samples.push({ fac, gaps, top, fr });
    }
    if (cap == null) { const nom = samples[(n - 1) / 2]; cap = freqCap(P, md0, nom.top, sfam); }
    const fminS = Math.min(P.solver.fmin, cap * 0.999);
    samples.forEach(s => { s.gaps = s.gaps.filter(g => g.lo < cap && g.hi > fminS); });
    S.sweep = { key, label: `${$('#sp').selectedOptions[0].textContent} (${FAM_LABEL[sfam].toLowerCase()} waves)`, samples, P, cap, fmin: fminS };
    setStatus(`Sweep done in ${((performance.now() - t0) / 1000).toFixed(1)} s`);
    drawSweep(); renderSens();
  } catch (e) { setStatus('Sweep stopped: ' + e.message); if (!/stopped by the user/.test(e.message)) console.error(e); }
  endRun();
}
$('#sweep').addEventListener('click', runSweep);
function drawSweep() {
  const cv = $('#cvSweep'); const { ctx, w, h } = prep(cv); const SW = S.sweep && !!S.sweep.drawn === isDrawn() ? S.sweep : null;
  const ink2 = css('--ink-2'), grid = css('--grid');
  if (!SW) { ctx.fillStyle = ink2; ctx.font = '15px system-ui, sans-serif'; ctx.fillText('Pick a parameter and press “Run sweep”.', 20, 40); return; }
  const sm = SW.samples, lmin = Math.log(sm[0].fac), lmax = Math.log(sm[sm.length - 1].fac);
  const fmax = SW.cap, fmin = SW.fmin || 0;
  const m = { l: 62, r: 18, t: 16, b: 44 }, U = fUnit(fmax);
  const X = lf => m.l + (lf - lmin) / (lmax - lmin) * (w - m.l - m.r), Y = f => h - m.b - (f - fmin) / (fmax - fmin) * (h - m.t - m.b);
  ctx.strokeStyle = grid; ctx.lineWidth = 1; ctx.font = '13px system-ui, sans-serif'; ctx.fillStyle = ink2;
  for (const v of niceTicks(fmin / U.d, fmax / U.d, 6)) { const y = Y(v * U.d); ctx.beginPath(); ctx.moveTo(m.l, y); ctx.lineTo(w - m.r, y); ctx.stroke(); ctx.textAlign = 'right'; ctx.fillText(+v.toFixed(3), m.l - 7, y + 4); }
  ctx.textAlign = 'center';
  sm.forEach(s => { const x = X(Math.log(s.fac)); ctx.beginPath(); ctx.moveTo(x, m.t); ctx.lineTo(x, h - m.b); ctx.stroke(); });
  [sm[0], sm[(sm.length - 1) / 2], sm[sm.length - 1]].forEach(s => ctx.fillText(`× ${s.fac < 1 ? s.fac.toFixed(2) : s.fac.toFixed(s.fac >= 10 ? 0 : 2)}`, X(Math.log(s.fac)), h - m.b + 16));
  ctx.fillText(`${SW.label}, relative to the current design (log scale)`, (m.l + w - m.r) / 2, h - 6);
  ctx.save(); ctx.translate(15, (m.t + h - m.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(`Frequency (${U.u})`, 0, 0); ctx.restore();
  const step = (lmax - lmin) / (sm.length - 1);
  sm.forEach(s => { const lf = Math.log(s.fac), xa = X(lf - step / 2), xb = X(lf + step / 2);
    for (const g of s.gaps) { if (g.hi < fmin || g.lo > fmax) continue; const y1 = Y(Math.min(g.hi, fmax)), y0 = Y(Math.max(g.lo, fmin)); ctx.fillStyle = css('--gap-tint'); ctx.fillRect(xa, y1, xb - xa, y0 - y1); ctx.globalAlpha = 0.55; ctx.fillStyle = hatch(ctx, css('--gap')); ctx.fillRect(xa, y1, xb - xa, y0 - y1); ctx.globalAlpha = 1; } });
  // edges as lines for gaps tracked by band index
  const ids = [...new Set(sm.flatMap(s => s.gaps.map(g => g.below)))];
  ctx.strokeStyle = css('--ink'); ctx.lineWidth = 1.6;
  for (const id of ids) for (const edge of ['lo', 'hi']) { ctx.beginPath(); let on = false; sm.forEach(s => { const g = s.gaps.find(x => x.below === id); if (!g || g[edge] > fmax || g[edge] < fmin) { on = false; return; } const x = X(Math.log(s.fac)), y = Y(g[edge]); on ? ctx.lineTo(x, y) : ctx.moveTo(x, y); on = true; }); ctx.stroke(); }
  ctx.save(); ctx.beginPath(); ctx.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b); ctx.clip();
  if (sm.some(s => s.fr)) { ctx.strokeStyle = css('--res'); ctx.setLineDash([6, 4]); ctx.lineWidth = 1.6; ctx.beginPath(); sm.forEach((s, i) => { const x = X(Math.log(s.fac)), y = Y(s.fr); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }); ctx.stroke(); ctx.setLineDash([]);
    const last = sm[sm.length - 1]; ctx.fillStyle = css('--res'); ctx.textAlign = 'right'; ctx.fillText('isolated resonator fᵣ', w - m.r - 4, Math.max(m.t + 12, Y(last.fr) - 6)); }
  ctx.restore();
  const xn = X(0); ctx.strokeStyle = css('--ink'); ctx.setLineDash([2, 3]); ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(xn, m.t); ctx.lineTo(xn, h - m.b); ctx.stroke(); ctx.setLineDash([]);
}
function renderSens() {
  const SW = S.sweep, sm = SW.samples, c = (sm.length - 1) / 2, a = sm[c - 1], b = sm[c + 1], dl = Math.log(b.fac) - Math.log(a.fac);
  const rows = sm[c].gaps.map(g => {
    const ov = (x, y) => Math.min(x.hi, y.hi) - Math.max(x.lo, y.lo), near = G => { let best = G.find(x => x.below === g.below && ov(x, g) > 0) || null; if (!best) for (const x of G) if (ov(x, g) > 0 && (!best || ov(x, g) > ov(best, g))) best = x; return best; };
    const ga = near(a.gaps), gb = near(b.gaps);
    const s = (k) => (ga && gb) ? ((Math.log(gb[k]) - Math.log(ga[k])) / dl) : null;
    const wa = ga && (ga.hi - ga.lo), wb = gb && (gb.hi - gb.lo);
    const sw = (ga && gb && wa > 0 && wb > 0) ? (Math.log(wb) - Math.log(wa)) / dl : null;
    const f = v => v == null || !isFinite(v) ? '–' : (v >= 0 ? '+' : '−') + Math.abs(v).toFixed(2);
    return `<tr><td>${g.below} and ${g.below + 1}</td><td>${fmtF(g.lo)}</td><td>${fmtF(g.hi)}</td><td>${f(s('lo'))}</td><td>${f(s('hi'))}</td><td>${f(sw)}</td></tr>`;
  });
  $('#sensTable').innerHTML = rows.length
    ? `<table><thead><tr><th>Gap between bands</th><th>Lower edge</th><th>Upper edge</th><th>S lower</th><th>S upper</th><th>S width</th></tr></thead><tbody>${rows.join('')}</tbody></table><p class="note" style="padding:0 12px 10px">Varied: ${SW.label.charAt(0).toLowerCase() + SW.label.slice(1)}. For an ideal mass-in-spring resonator the lower edge follows fᵣ ∝ √(k/m), so S ≈ +0.5 for stiffness and −0.5 for mass; departures from ±0.5 show how strongly the host beams take part.</p>`
    : `<p class="empty">The nominal design has no complete gap in the computed range, so there is nothing to differentiate. The gap map above still shows where gaps appear as the parameter changes.</p>`;
}

// ---------- method text ----------
$('#methodText').innerHTML = `
<h3>What the workbench solves</h3>
<p>The lattice lies in the x–y plane and every member is a slender beam. For a flat lattice whose members are symmetric about its mid-plane, in-plane motion (u, v, rotation about z) and out-of-plane motion (w, rotations about x and y) decouple exactly, so the two polarisations are solved as separate problems, each with three degrees of freedom per node: w and the two cross-section rotations for out-of-plane (flexural and flexural–torsional) waves, and u, v, θ<sub>z</sub> for in-plane waves. A full band gap is a frequency range that lies in a gap of both.</p>
<h3>Unit-cell geometry</h3>
<p>The hexagon cell has one vertical wall of length h and two inclined walls of length l at angle θ to the horizontal, with lattice vectors a₁ = (2l cos θ, 0) and a₂ = (l cos θ, h + l sin θ). θ = 30° with h = l is the regular honeycomb; other values give stretched cells, and negative θ gives the re-entrant (auxetic) cell, which stays closed while h + 2l sin θ > 0. The vertical wall can be made wider, as in expanded honeycomb cores where it is a double wall.</p>
<p>Resonators can sit at the junctions, at wall mid-spans or at any fraction t along chosen walls; the wall is split there so the arm is clamped to a real node. Arms can point into the cell, out of it, to both sides, or at a fixed angle, and can be straight or L-shaped. The central-mass option places a cube at the cell centre carried by ligaments to the cell corners or to the wall mid-spans.</p>
<p>For anything else, copy the design into the custom editor. A cell is a list of nodes (x, y in mm, optional mass in g, rotary inertia J in g·mm², role host or res), beams between nodes i and j where “shift” is the cell offset [n₁, n₂] of node j, named sections, and optional mass-on-spring resonators. Every beam is meshed in proportion to its length, and anything marked res counts as resonator in the mass ratio, the energy colouring and the sensitivity sweep.</p>
<h3>Choosing the zone path</h3>
<p>The workbench tests which point operations about the cell centre map the whole pattern onto itself. With six-fold (or four-fold) symmetry it uses the standard Γ–X–K–Γ (or Γ–X–M–Γ) path, where for the hexagonal zone X is the edge midpoint (often called M) and K the corner; this path is the default and can be forced from the Solver panel even when resonators lower the symmetry, in which case it is worth confirming the gap with the half zone; with two mirror lines it walks the quarter zone Γ–X–V–Y–Γ around the zone boundary; otherwise it walks the half zone, which is always sufficient because ω(k) = ω(−k). One-sided arms, arms off mid-span and angled arms usually break the mirrors, so expect the longer path for those cells.</p>
<h3>Beam elements: Timoshenko</h3>
<p>Each beam is split into two-node Timoshenko elements, which include shear deformation and the rotary inertia of the cross-section; St-Venant torsion is linear. The element uses interdependent interpolation, so it is exact in statics and does not lock in shear, with the consistent mass of Friedman and Kosmatka. In local coordinates (element length ℓ, DOFs v₁, θ₁, v₂, θ₂):</p>
<div class="eq">Φ = 12EI / (κGAℓ²),&nbsp;&nbsp; K<sub class="up">b</sub> = EI / ((1 + Φ)ℓ³) [12, 6ℓ, −12, 6ℓ; 6ℓ, (4 + Φ)ℓ², −6ℓ, (2 − Φ)ℓ²; …],&nbsp;&nbsp; K<sub class="up">t</sub> = GJ/ℓ [1, −1; −1, 1]</div>
<div class="eq">M<sub class="up">b</sub> = ρAℓ/(1 + Φ)² [13/35 + 7Φ/10 + Φ²/3, …] + ρI/((1 + Φ)²ℓ) [6/5, (1/10 − Φ/2)ℓ, …],&nbsp;&nbsp; M<sub class="up">t</sub> = ρI<sub class="up">p</sub>ℓ/6 [2, 1; 1, 2]</div>
<p>κ is Cowper's shear coefficient, 10(1 + ν)/(12 + 11ν) for a rectangle and 6(1 + ν)/(7 + 6ν) for a circle. With Φ = 0 and no rotary inertia the element reduces exactly to Euler–Bernoulli, which stays selectable in the Solver panel for comparison. Shear matters when a member is not slender in its bending plane: here the 5 mm tall walls bend out of the plane over only about 8 mm, and the blocks are as deep as they are long. The nodal rotations are cross-section rotations: ψ = c φ<sub>x</sub> + s φ<sub>y</sub> about the element's transverse axis and the twist χ = −s φ<sub>x</sub> + c φ<sub>y</sub>, with (c, s) the element direction. I = bh³/12 for bending out of the plane; J uses the thin-rectangle series and I<sub>p</sub> = I<sub>y</sub> + I<sub>z</sub>. Rigid joints share all three nodal values. On an infinite straight beam the element reproduces the exact Timoshenko dispersion relation.</p>
<h3>In-plane frame elements</h3>
<p>For in-plane waves each beam is a frame element: a two-node bar for stretching and a Timoshenko beam for bending in the plane, with the in-plane bending stiffness EI<sub>z</sub>, I<sub>z</sub> = hb³/12, where b is the width in the plane and h the height out of it. The bar has K = EA/ℓ [1, −1; −1, 1] and consistent mass ρAℓ/6 [2, 1; 1, 2]; the bending part uses the same Timoshenko matrices as above with EI<sub>z</sub> and rotary inertia ρI<sub>z</sub> (lumped ρI<sub>z</sub>ℓ/2 at each node in the Euler–Bernoulli option). Local axial and transverse displacements follow from u<sub>a</sub> = c u + s v and v<sub>t</sub> = −s u + c v. Point masses act on u and v, with their rotary inertia about z; a mass on a spring is taken to be isotropic in the plane, with the same k along x and y. Because walls are usually much thinner in the plane than out of it, a honeycomb is far softer in-plane and its in-plane bands sit much lower; the in-plane long-wavelength shear wave speed reproduces the Gibson–Ashby honeycomb shear modulus to within 1 %.</p>
<h3>Resonators</h3>
<p>A mass on a spring adds one degree of freedom z with energy ½k(w − z)² and kinetic energy ½mż² (two, along x and y, for in-plane waves). A cantilever arm is meshed with three beam elements clamped to the junction, with the cube as a point mass and rotary inertia m s²/6 at the tip. The isolated frequencies shown are f<sub>r</sub> = (1/2π)√(k/m) and, for the arm, (1/2π)√(3EI<sub>r</sub>/L³ ÷ (m<sub>tip</sub> + 0.24 ρA<sub>r</sub>L)), with I<sub>r</sub> the arm's out-of-plane or in-plane bending inertia.</p>
<h3>Bloch–Floquet boundary conditions</h3>
<p>For a wave vector k, a degree of freedom in the cell displaced by R = n₁a₁ + n₂a₂ equals the one in the reference cell times e<sup>i k·R</sup>. Elements that cross the cell boundary are assembled with that factor on their far node, which gives Hermitian matrices and the eigenproblem</p>
<div class="eq">[ K(k) − ω² M(k) ] u = 0</div>
<p>solved for ω(k) along Γ–X–M–Γ (square) or Γ–X–K–Γ (hexagonal lattices). The problem is mapped to a real symmetric one of twice the size, reduced with a Cholesky factor of M, and solved by Householder tridiagonalisation and implicit QL. The colour of each point is the share of kinetic energy carried by resonator degrees of freedom, so flat orange bands mark local resonance and blue bands mark host waves. A complete gap is reported where the highest point of band n lies below the lowest point of band n + 1 everywhere on the path; extrema inside the zone are not checked, which is the usual caveat.</p>
<h3>Complex band: phase and attenuation</h3>
<p>The right-hand plot inverts the question: frequency is fixed and the propagation constant μ = k·a₁ is found, complex in general, with the phase along a₂ held at zero (this is the Γ–X leg). Beams that cross the a₁ boundary are cut, the interior is condensed out of the dynamic stiffness D = K − ω²M, and the left and right faces satisfy u<sub>R</sub> = λu<sub>L</sub>, λ = e<sup>−iμ</sup>. Force balance across the cut gives the quadratic eigenproblem</p>
<div class="eq">[ D<sub class="up">RL</sub> + λ (D<sub class="up">LL</sub> + D<sub class="up">RR</sub>) + λ² D<sub class="up">LR</sub> ] u<sub class="up">L</sub> = 0</div>
<p>solved through λ = (1 + s)/(1 − s), which keeps it regular when D<sub>LR</sub> is rank-deficient (lattices where several cut beams meet at one joint). Re(μ) is the phase change per cell; Im(μ) is the attenuation in nepers per cell (× 8.686 for dB per cell). Inside a gap every wave has Im(μ) > 0: in a Bragg gap Re(μ) sticks at 0 or π, while a local-resonance gap shows a sharp attenuation peak just above f<sub>r</sub>. Standing modes in which the cut faces do not move, such as a beam vibrating between two stationary joints, produce flat bands that the transfer formulation cannot see; they appear only in the left-hand plot.</p>
<h3>Sensitivity to mass and stiffness</h3>
<p>The sweep scales one quantity, recomputes the bands and reports S = ∂ ln f<sub>edge</sub> / ∂ ln p by central differences at the nominal point. For a mass-in-spring resonator on a heavy host the lower edge sits near f<sub>r</sub>, so S ≈ −½ for resonator mass and +½ for resonator stiffness, and the relative width grows roughly like √(1 + μ) − 1 with mass ratio μ. Host stiffness and mass mainly move Bragg gaps and the frequency at which the host branch meets the resonance, which is why a local-resonance gap is far less sensitive to them.</p>
<h3>Band gaps versus edge waves</h3>
<p>This page computes bulk properties only. A band gap is a property of the infinite lattice: no real-k solution exists at that frequency in any direction, and the complex band shows every wave decaying. An edge wave is a solution of a finite or semi-infinite lattice that sits inside the bulk gap and decays away from a boundary. To separate them numerically, build a supercell ribbon that is periodic along the edge and N cells deep (free, clamped or interface termination), compute its bands against the edge wave number, and overlay the projected bulk bands. Modes that cross the bulk gap and whose energy fraction within a few cells of the boundary approaches one are edge states; bulk modes have that fraction near the geometric share, about 2/N. If the gap comes from breaking a lattice symmetry, the Zak phase or valley Chern number of the bulk bands on either side of the interface tells you whether those edge states are topologically protected. Comparing the edge-wave branch with the complex band also helps: its decay into the bulk should match the least attenuated bulk wave at that frequency.</p>
<h3>How edge waves are detected</h3>
<p>The Edge waves tab builds a strip N cells wide that repeats along a₁ (or a₂) and has free edges on both long sides, then solves its bands against the phase per cell along the edge. The bulk bands are projected onto the same edge wave number by sweeping the other wave-vector component through a full period. Each strip mode gets an edge share: the fraction of its kinetic energy in the outer cells on either side. A mode is called an edge wave when it lies outside every projected bulk band and its edge share is at least 60 %; an evenly spread mode would have about 2w/N per side. The strip uses two elements per wall, which reproduces the bulk gap edges of the default design to within 1 Hz.</p>
<h3>Drawn cells (cell type “Drawn cell”)</h3>
<p>A drawn cell is a list of parts in coordinates: beams (width b in the plane, height h, axis height z), shell walls (a vertical strip of thickness t and height H drawn as a line), shell plates (a polygon at height z with thickness t and any number of holes), point masses and mass-springs, each marked host or resonator. Building the model meshes every part at the element size: beams into two-node 3D Timoshenko elements, walls into a structured triangle strip, plates by a conforming Delaunay triangulation of the outline, holes and any points lying on them. Points that meet are welded, and a point that is a lattice translate of another is identified with it, so parts may cross the cell outline; a point that lands on another part's edge splits that edge, so T-joints connect.</p>
<p>The shell element is a flat triangle: constant-strain membrane, Discrete Kirchhoff (DKT, Batoz 1980) bending, and a drilling rotation tied to the membrane rotation by a penalty G t A (Hughes–Brezzi), which carries the in-plane moment of a beam joined to a shell edge. It passes the membrane and bending patch tests exactly (D = Et³/12(1 − ν²), D₆₆ = Gt³/12 from the homogenisation of a plain plate) and gives the plate flexural wave within 0.2 % at 8 elements per wavelength. The constant-strain membrane is stiff when a wall bends in its own plane, so for tall thin walls use an element size of a quarter of the height or less, or draw the walls as beams (for the 14 mm hexagon the shell walls give D₁₁ 31 N·m at 0.6 mm against 28.2 N·m for beams, while their twisting stiffness is higher because the joints restrain warping).</p>
<p>The full Bloch FE solves [K(k) − ω²M(k)]u = 0 for every node's six degrees of freedom by shift-invert block Krylov iteration (complex Cholesky of K − σM, M-orthonormal basis with full re-orthogonalisation, Rayleigh–Ritz), so models with a few thousand degrees of freedom stay practical. Each mode is sorted as mostly out-of-plane or mostly in-plane by the share of its translational kinetic energy along z; for a cell symmetric about its mid-plane this separation is exact. Validation: your hexagon imported as a drawing reproduces the beam-lattice Bloch FE gap to 10⁻⁵ (2091.8–2961.9 Hz).</p>
<h3>Hexagonal frame templates</h3>
<p>Three ready-made cells build an editable drawing: (a) a hexagonal frame with a central hexagonal mass on six spokes, with an optional inner ring and blocks on the spokes; (b) the same mass on two beams along one diagonal; (c) the electrospun cell, a mass carried by a thin membrane spanning the frame. Hexagons of side L tile with a₁ = (1.5L, √3L/2), a₂ = (0, √3L); each cell carries three of its six walls and the neighbours supply the rest, so shared walls are not doubled. The frame is the host and everything inside it the resonator. The membrane is modelled as a thin shell without prestress: a tensioned electrospun mat is stiffer, so raise its E or thickness to match a measured resonance.</p>
<h3>Complex band and mode shapes of drawn cells</h3>
<p>The complex band of a drawn cell uses wave finite elements along a₁ with the phase along a₂ held at zero (the Γ–X leg). The cell's dynamic stiffness is split as D(λ) = A₀ + λA₊ + λ⁻¹A₋ (λ = e<sup>iμ</sup>, μ = k·a₁) by the cell offset of each element coupling, the degrees of freedom not touched by A± are condensed at each frequency, and the quadratic is solved through λ = (1 + s)/(1 − s). For a flat cell the out-of-plane and in-plane degrees of freedom are solved separately; for a cell that is not flat all waves are shown together. Click a band point to see the deformed mesh with the undeformed shape in grey; drag to rotate.</p>
<h3>Editable cells: shapes, moving and copying</h3>
<p>Any cell can be turned into an editable cell (“Edit this cell freely”); every wall, arm and mass becomes a part. The Shape tool adds a regular polygon with any number of sides, a rectangle or a circle (a polygon with many sides), stored by its numbers: centre, size (corner radius, side length or inner radius), rotation, and its use as a shell plate, a frame of beams along its edges, shell walls along its edges, or a hole cut from the plate under its centre. Each part can be dragged, resized by its handles, typed in exactly, turned, scaled, mirrored, duplicated or copied around a point, and given its own role, material and thickness; nothing has to be symmetric.</p>
<h3>Zone path and symmetry</h3>
<p>For editable cells the workbench looks for every mirror line and rotation (about any centre) that maps the periodic drawing onto itself, comparing beams, walls, plate edges and masses with their sizes, heights, roles and materials. The standard path Γ–X–K–Γ (Γ–X–M–Γ) is used only when the cell has six-fold (four-fold) symmetry with a mirror; with two mirrors the quarter-zone boundary is used, otherwise the half-zone boundary, which is enough for any cell because ω(k) = ω(−k). The same rule now applies to the parametric cells; the forced standard path remains selectable. A custom path takes your own points, as fractions of b₁ and b₂ or as (k<sub>x</sub>, k<sub>y</sub>) in units of π/a, and the preview in the Solver panel shows it on the Brillouin zone.</p>
<h3>Sensitivity and edge waves of editable cells</h3>
<p>The sensitivity sweep rebuilds the cell for every sample: resonator or host density and modulus are scaled through the material of each part (point masses and mass-springs included), the cell size scales all in-plane coordinates and widths, and “one part” scales a chosen dimension of a chosen part. Each sample is solved with full Bloch FE at half the k-point density. Edge waves use a strip of N copies of the cell, periodic along the chosen edge and free on both long sides; it is solved with the same beam and shell elements, with the nodes ordered across the strip so that the Cholesky factor is banded (cost grows linearly with N). The bulk bands are projected by sweeping the other wave-vector component over eight values per edge wave number.</p>
<h3>Ready-made cells and Build your own</h3>
<p>Cell A is the parametric honeycomb with a cantilever and tip block, solved with the beam-lattice Bloch FE. Cells B and C, and every “build your own” cell, are generated as editable drawings from the Lattice and Resonator boxes: frame walls are host beams, spokes, beams across, arms and rings are resonator beams, centre masses and blocks are shell plates, a membrane is a shell plate with the Membrane material. They are solved with the beam and shell Bloch FE. The generator only writes a drawing, so whatever it builds can be edited part by part.</p>
<h3>Thickness reference and offsets</h3>
<p>A beam or plate whose thickness is measured upward (or downward) from its height z keeps its nodes at z, where it joins other parts, while its axis or mid-surface lies h/2 (t/2) above (below). The element matrices are transformed by the rigid link u<sub>mid</sub> = u + θ × (0, 0, d), i.e. K′ = TᵀKT and M′ = TᵀMT with u<sub>x</sub> += d θ<sub>y</sub> and u<sub>y</sub> −= d θ<sub>x</sub>. Two beams stacked this way reproduce the composite-beam bending stiffness EI about the shifted neutral axis exactly. Offsets couple in-plane and out-of-plane motion, so such cells report the complete gap only.</p>
<h3>Edge waves in a finite sample</h3>
<p>The finite sample repeats the drawn cell M × N times with every edge free and no Bloch condition. Its real symmetric matrices are ordered row by row across the shorter side, so K − ω<sub>c</sub>²M has a band of about one row of cells; it is factorised as LDLᵀ (indefinite, the shift is nudged if it hits an eigenvalue), and shift-invert block Krylov iteration returns the modes nearest ω<sub>c</sub>. Each mode is placed by the share of its kinetic energy within half a cell of the boundary and in the four corner cells; it is an edge mode when that share exceeds both 60 % and the zone's area share plus 25 %, and a corner mode when the corner cells alone carry at least 50 % (and 2.5 times their area share). Against a dense solution of the same 2 × 2 sample the modes agree to all printed digits.</p>
<h3>Driven sample and topology check</h3>
<p>The driven response of the finite sample is a modal sum over the modes found near the centre frequency, u(ω) = Σ φ<sub>n</sub> (φ<sub>n</sub><sup>H</sup> F) / (ω<sub>n</sub>² − ω² + 2iζω<sub>n</sub>ω) with ζ = 2 %, for a unit force on the node nearest the shaker cell's centre; the energy per cell is Σ m<sub>j</sub>|u<sub>j</sub>|² over its translations. The staggered rectangle repeats the cell in N columns along a₁, shifting every other column by a₂ − a₁.</p>
<p>The Chern number of a band of a passive (time-reversal symmetric) lattice is zero, because K(−k) = K(k)* and M(−k) = M(k)* make the Berry curvature odd in k. The Zak phase of the bands below a gap is the Wilson loop φ = −arg det Π U<sub>s</sub>, with links U<sub>s</sub> = Ψ(k<sub>s</sub>)<sup>H</sup> diag(m<sub>j</sub> e<sup>−iΔk·x<sub>j</sub></sup>) Ψ(k<sub>s+1</sub>) over 24 wave vectors along b₁, b₂ and b₁ + b₂ (the node positions x<sub>j</sub> carry the intra-cell phase, so the loop closes on itself), referred to the inversion centre nearest the middle of the drawn cell. With 180° symmetry it is 0 or π; shifting the cell by half a lattice vector changes a single band's phase by exactly π, which is a test the implementation passes.</p>
<h3>Limits</h3>
<p>Timoshenko beams capture shear deformation and rotary inertia, but they are still one-dimensional: the cross-section does not warp or distort, and the joint regions where walls meet, and the solid blocks, are represented by their centrelines. Joints are rigid and massless apart from the beams that meet there. No damping is included, so attenuation inside a passband is zero. For thick members, plate-like cells or 3D cube resonators, use these results to set up and check a solid-element Bloch analysis in a commercial FE code.</p>`;

// ---------- startup ----------
// ---------- edge waves: strip model ----------
const EDGE_L = 0.6;
// the parametric cell as a drawing, with its gaps, for the finite sample and the topology check
async function paramSource(famSel) {
  if (!S.band || S.bandKey !== JSON.stringify(readParams())) { await runAnalysis(false); if (!S.band) return null; }
  const fam = famSel === 'ip' ? 'ip' : 'oop', P = S.P, spec = buildSpec(P), dr = Drw.fromSpec(spec);
  const lw = (spec.walls || []).map(w => Math.hypot(w.p1[0] - w.p0[0], w.p1[1] - w.p0[1])); dr.mesh = lw.length ? Math.max(0.5, 1000 * Math.min(...lw) / P.ne) : 2;
  return { dr, gaps: famSel === 'all' ? S.complete : S.fam[fam].gaps, fmax: S.fam[fam].band.fmax, fmin: S.fam[fam].band.fmin || 0, mat: P.mat, beam: P.beam, key: S.bandKey };
}
async function runEdge() {
  if (S.running) return;
  if ($('#eKind').value === 'finite') {
    if (isDrawn()) { await DrawUI.runFinite(); return; }
    if (!S.band || S.bandKey !== JSON.stringify(readParams())) { await runAnalysis(false); if (!S.band) return; }
    const src = await paramSource($('#eFam').value); if (!src) return;
    await DrawUI.runFinite(src);
    renderAssessment(); return;
  }
  if (isDrawn()) { if (typeof DrawUI !== 'undefined') await DrawUI.runEdge(); return; }
  if (!S.band || S.bandKey !== JSON.stringify(readParams())) { await runAnalysis(false); if (!S.band) return; }
  const P = S.P, N = Math.max(3, Math.min(20, Math.round(num('eN')))), dir = $('#eDir').value, nk = Math.max(5, Math.min(61, Math.round(num('eK'))));
  const efam = $('#eFam').value, FB = S.fam[efam].band;
  const w = Math.max(0.5, Math.min(N / 2, num('eW') || 2));
  startRun();
  const t0 = performance.now();
  try {
    const spec = buildSpec(P), neR = Math.min(P.ne, 2);
    const mdB = Eng.buildModel(spec, P.mat, { ne: neR, path: 'half', family: efam, beam: P.beam });
    const rs = Eng.ribbonSpec(spec, N, dir), mdR = Eng.buildModel(rs, P.mat, { ne: neR, path: 'half', family: efam, beam: P.beam });
    const B1 = dir === 'a2' ? mdB.Lat.b2 : mdB.Lat.b1, B2 = dir === 'a2' ? mdB.Lat.b1 : mdB.Lat.b2;
    const depth = Eng.nodeDepth(mdR, rs.ribbon.A1, rs.ribbon.A2), mdiag = Eng.massDiag(mdR);
    const fmin = FB.fmin || 0, fmax = FB.fmax, nbB = Math.min(mdB.n, FB.nb + 2), NU = 12, tol = 0.003 * fmax;
    let nR = Math.min(mdR.n, (FB.nb + 2) * N);
    const mus = [], proj = [], modes = [];
    for (let i = 0; i < nk; i++) {
      const mu = Math.PI * i / (nk - 1); mus.push(mu);
      const lo = new Array(nbB).fill(Infinity), hi = new Array(nbB).fill(-Infinity);
      for (let v = 0; v < NU; v++) { const nu = TAU * v / NU, k = [(mu * B1[0] + nu * B2[0]) / TAU, (mu * B1[1] + nu * B2[1]) / TAU];
        Eng.solveK(mdB, k[0], k[1], nbB).f.forEach((x, b) => { lo[b] = Math.min(lo[b], x); hi[b] = Math.max(hi[b], x); }); }
      proj.push(lo.map((x, b) => [x, hi[b]]));
      const kr = [mu / TAU * mdR.Lat.b1[0], mu / TAU * mdR.Lat.b1[1]];
      let r = Eng.solveK(mdR, kr[0], kr[1], nR, true, false);
      while (r.f[r.f.length - 1] < fmax && nR < mdR.n) { nR = Math.min(mdR.n, Math.ceil(nR * 1.5)); r = Eng.solveK(mdR, kr[0], kr[1], nR, true, false); }
      const row = [];
      r.f.forEach((f, b) => {
        if (f > fmax * 1.02 || f < fmin * 0.98) return;
        const sh = Eng.edgeShares(mdR, r.vecs[b], depth, mdiag, w), inGap = !proj[i].some(([a, c]) => f >= a - tol && f <= c + tol);
        row.push({ f, bot: sh.bot, top: sh.top, L: Math.max(sh.bot, sh.top), inGap, mu, u: r.vecs[b], k: kr });
      });
      modes.push(row);
      setStatus(`Edge analysis: k-point ${i + 1} of ${nk}`, (i + 1) / nk); await tick();
    }
    S.edge = { mus, proj, modes, N, dir, w, fmin, fmax, mdR, key: S.bandKey, family: efam };
    S.emode = null; S.esel = null;
    setStatus(`Edge analysis done in ${((performance.now() - t0) / 1000).toFixed(1)} s`);
    renderAssessment(); showTab('edge');
  } catch (e) { setStatus('Edge analysis stopped: ' + e.message); if (!/stopped by the user/.test(e.message)) console.error(e); }
  endRun();
}
$('#edgeRun').addEventListener('click', runEdge);
function drawEdge() {
  if (S.edge && S.edge.kind === 'finite' && !!S.edge.drawn === isDrawn() && $('#eKind').value === 'finite') { DrawUI.drawFinite(); return; }
  if (typeof DrawUI !== 'undefined') DrawUI.edgeTitles(false);
  const cv = $('#cvEdge'); const { ctx, w, h } = prep(cv); const E = S.edge && S.edge.kind !== 'finite' && !!S.edge.drawn === isDrawn() && $('#eKind').value === 'strip' ? S.edge : null, ink2 = css('--ink-2');
  if (!E) { ctx.fillStyle = ink2; ctx.font = '15px system-ui, sans-serif'; ctx.fillText('Press “Run edge analysis” to compute the strip bands.', 20, 40); return; }
  const m = { l: 58, r: 14, t: 16, b: 40 }, fmin = E.fmin, fmax = E.fmax, U = fUnit(fmax);
  const X = mu => m.l + mu / Math.PI * (w - m.l - m.r), Y = f => h - m.b - (f - fmin) / (fmax - fmin) * (h - m.t - m.b);
  S.edgeGeom = { X, Y };
  ctx.save(); ctx.beginPath(); ctx.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b); ctx.clip();
  ctx.fillStyle = css('--grid');
  const nbB = Math.min(...E.proj.map(p => p.length));
  for (let b = 0; b < nbB; b++) { ctx.beginPath(); E.mus.forEach((mu, i) => { const y = Y(E.proj[i][b][1]); i ? ctx.lineTo(X(mu), y) : ctx.moveTo(X(mu), y); }); for (let i = E.mus.length - 1; i >= 0; i--) ctx.lineTo(X(E.mus[i]), Y(E.proj[i][b][0])); ctx.closePath(); ctx.fill(); }
  for (const g of (E.gaps || S.fam[E.family].gaps)) { if (g.hi < fmin || g.lo > fmax) continue; const y1 = Y(Math.min(g.hi, fmax)), y0 = Y(Math.max(g.lo, fmin)); ctx.globalAlpha = 0.45; ctx.fillStyle = hatch(ctx, css('--gap')); ctx.fillRect(m.l, y1, w - m.l - m.r, y0 - y1); ctx.globalAlpha = 1; }
  const uni = E.uni != null ? E.uni : E.w / E.N, Lthr = E.Lthr || EDGE_L, evan = css('--evan'), edge = css('--edge'), ink = css('--ink');
  E.modes.forEach((row, i) => row.forEach((md, j) => {
    const t = Math.max(0, Math.min(1, (md.L - uni) / (Math.min(1, uni + 0.45) - uni))), isEdge = md.inGap && md.L >= Lthr;
    ctx.fillStyle = mix(evan, edge, t); ctx.beginPath(); ctx.arc(X(md.mu), Y(md.f), isEdge ? 4.2 : 2.6, 0, TAU); ctx.fill();
    if (isEdge) { ctx.strokeStyle = ink; ctx.lineWidth = 1.2; ctx.stroke(); }
  }));
  if (S.esel) { const md = E.modes[S.esel.i][S.esel.j]; ctx.strokeStyle = ink; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(X(md.mu), Y(md.f), 8, 0, TAU); ctx.stroke(); }
  ctx.restore();
  ctx.font = '13px system-ui, sans-serif'; ctx.fillStyle = ink2; ctx.strokeStyle = css('--rule'); ctx.lineWidth = 1;
  for (const v of niceTicks(fmin / U.d, fmax / U.d, 6)) { const y = Y(v * U.d); ctx.textAlign = 'right'; ctx.fillText(+v.toFixed(3), m.l - 7, y + 4); }
  ctx.textAlign = 'center';
  [[0, '0'], [Math.PI / 2, 'π/2'], [Math.PI, 'π']].forEach(([mu, sx]) => ctx.fillText(sx, X(mu), h - m.b + 16));
  ctx.fillText(`Phase per cell along the edge, k∥·${E.dir === 'a2' ? 'a₂' : 'a₁'} (rad)`, (m.l + w - m.r) / 2, h - 6);
  ctx.save(); ctx.translate(15, (m.t + h - m.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(`Frequency (${U.u})`, 0, 0); ctx.restore();
  ctx.strokeRect(m.l, m.t, w - m.l - m.r, h - m.t - m.b);
  $('#edgeSub').textContent = `${(E.famLabel || FAM_LABEL[E.family]).toLowerCase()} waves, ${E.N} cells, free edges along ${E.dir === 'a2' ? 'a₂' : 'a₁'}`;
}
function drawEdgeMode() { if (S.edge && S.edge.kind === 'finite') { if (!!S.edge.drawn === isDrawn()) DrawUI.drawFiniteMode(); return; } if (S.edge && S.edge.drawn) { if (isDrawn() && typeof DrawUI !== 'undefined') DrawUI.drawEdgeMode(); return; } drawMode('#cvEdgeMode', S.emode, S.edge && S.edge.mdR, [[-1, 0], [0, 0], [1, 0]]); }
$('#cvEdge').addEventListener('click', e => {
  if (S.edge && S.edge.kind === 'finite') { DrawUI.pickFinite(e); return; }
  const E = S.edge, g = S.edgeGeom; if (!E || !g || !!E.drawn !== isDrawn()) return;
  const rect = e.currentTarget.getBoundingClientRect(), mx = e.clientX - rect.left, my = e.clientY - rect.top;
  let best = null, bd = 18 * 18;
  E.modes.forEach((row, i) => row.forEach((md, j) => { const d = (g.X(md.mu) - mx) ** 2 + (g.Y(md.f) - my) ** 2; if (d < bd) { bd = d; best = { i, j }; } }));
  if (!best) return;
  S.esel = best; const md = E.modes[best.i][best.j];
  S.emode = { u: md.u, k: md.k, f: md.f };
  if (E.drawn) { $('#edgeModeSub').textContent = `${fmtF(md.f)}, k∥·a = ${(md.mu / Math.PI).toFixed(2)}π, ${Math.round(100 * md.L)} % of energy in the ${md.top >= md.bot ? 'top' : 'bottom'} ${E.w} cells${md.inGap ? ', outside the bulk bands' : ''}`; drawEdge(); drawEdgeMode(); return; }
  const side = md.top >= md.bot ? 'top' : 'bottom';
  $('#edgeModeSub').textContent = `${fmtF(md.f)}, k∥·a = ${(md.mu / Math.PI).toFixed(2)}π, ${Math.round(100 * md.L)} % of energy in the ${side} ${E.w} cells${md.inGap ? ', outside the bulk bands' : ''}`;
  drawEdge(); drawEdgeMode();
});
function directional(cplx) {
  const iv = []; let st = null, prev = null;
  for (const r of cplx.rows) { const blocked = r.ev && Math.min(...r.ev.map(q => q.att)) > PROP; if (blocked && st == null) st = r.f; if (!blocked && st != null) { iv.push([st, prev]); st = null; } if (blocked) prev = r.f; }
  if (st != null) iv.push([st, prev]);
  return iv.filter(([a, c]) => c > a);
}
function famAssessment(f) {
  const P = S.P, x = S.fam[f], B = x.band, gaps = x.gaps, { fr } = resInfo(P, x.md, f);
  const desc = g => {
    const b = g.below - 1, share = B.R.reduce((a, r) => a + (r[b] || 0), 0) / B.R.length;
    const rows = x.cplx ? x.cplx.rows.filter(r => r.ev && r.f >= g.lo && r.f <= g.hi) : [], att = rows.length ? Math.max(...rows.map(r => Math.min(...r.ev.map(q => q.att)))) : null;
    const type = share > 0.5
      ? `a local-resonance gap: the band at its lower edge flattens out towards the zone boundary and carries ${Math.round(100 * share)} % of its kinetic energy in the resonators${fr ? ` (isolated resonator about ${fmtF(fr)})` : ''}`
      : `a Bragg-type gap: the bands that bound it are mainly host waves (resonator share ${Math.round(100 * share)} %), so it comes from the periodicity of the walls`;
    return `<li><b>${fmtF(g.lo)} to ${fmtF(g.hi)}</b> (${(200 * (g.hi - g.lo) / (g.hi + g.lo)).toFixed(1)} % wide), ${type}${att != null ? `. Peak attenuation ${(8.686 * att).toFixed(1)} dB per cell along ${x.cplx.dir}` : ''}.</li>`;
  };
  if (gaps.length) return `<p><b>${FAM_LABEL[f]} band gap: yes.</b> ${gaps.length === 1 ? 'One complete gap' : gaps.length + ' complete gaps'} for ${FAM_WAVE[f]} waves between ${fmtF(B.fmin || 0)} and ${fmtF(B.fmax)}:</p><ul>${gaps.map(desc).join('')}</ul>`;
  const ok = x.cplx ? directional(x.cplx) : [];
  return `<p><b>${FAM_LABEL[f]} band gap: no.</b> No complete gap for ${FAM_WAVE[f]} waves between ${fmtF(B.fmin || 0)} and ${fmtF(B.fmax)}.${ok.length ? ` Along ${x.cplx.dir} alone, waves are blocked at ${ok.slice(0, 3).map(([a, c]) => `${fmtF(a)} to ${fmtF(c)}`).join(', ')}${ok.length > 3 ? ' and higher' : ''}; these are directional gaps, open in that direction only.` : ''}</p>`;
}
function fullAssessment() {
  const C = S.complete, sh = S.shared;
  if (C.length) return `<p><b>Full band gap (both polarisations): yes.</b> Neither out-of-plane nor in-plane waves propagate at ${C.map(g => `<b>${fmtF(g.lo)} to ${fmtF(g.hi)}</b>`).join(', ')}. Only there does the lattice stop vibration of every polarisation.</p>`;
  const notes = [];
  for (const [f, o] of [['oop', 'ip'], ['ip', 'oop']]) for (const g of S.fam[f].gapsAll.slice(0, 2)) {
    const cross = bandRanges(S.fam[o]).filter(r => r.hi > g.lo && r.lo < g.hi).map(r => r.b);
    if (cross.length) notes.push(`${FAM_LABEL[o].toLowerCase()} band${cross.length > 1 ? 's' : ''} ${cross.length > 2 ? `${cross[0]}–${cross[cross.length - 1]}` : cross.join(' and ')} propagate${cross.length > 1 ? '' : 's'} inside the ${FAM_LABEL[f].toLowerCase()} gap at ${fmtF(g.lo)} to ${fmtF(g.hi)}`);
  }
  return `<p><b>Full band gap (both polarisations): no</b> between ${fmtF(sh.fmin)} and ${fmtF(sh.fmax)}.${notes.length ? ' ' + notes.map(t => t[0].toUpperCase() + t.slice(1)).join('; ') + ', so vibration of that polarisation still passes there.' : ''}</p>`;
}
function renderAssessment() {
  let html;
  if (!S.fam) html = '<h3>What this geometry gives</h3><p>Run the analysis to get an assessment.</p>';
  else {
    let edgeHtml, edgeYes = null, edgeWhere = ''; const E = S.edge;
    const eg = E && E.kind !== 'finite' ? S.fam[E.family].gaps : [];
    if (E && E.kind === 'finite' && !E.drawn && E.key === S.bandKey) { const fs = DrawUI.finiteSummary(); edgeHtml = fs.html; edgeYes = fs.yes; edgeWhere = fs.where; }
    else if (E && (E.kind === 'finite' || E.drawn)) edgeHtml = '<p><b>Edge waves: not analysed for this design yet.</b> Open the Edge waves tab and press “Run edge analysis”.</p>';
    else if (!E) edgeHtml = '<p><b>Edge waves: not analysed yet.</b> Open the Edge waves tab and press “Run edge analysis”.</p>';
    else if (E.key !== S.bandKey) edgeHtml = '<p><b>Edge waves: out of date.</b> The design has changed since the last edge analysis; run it again.</p>';
    else {
      const gaps = eg, fw = FAM_WAVE[E.family];
      const all = E.modes.flat().filter(q => q.f > 0.01 * E.fmax), edge = all.filter(q => q.inGap && q.L >= EDGE_L), uni = E.w / E.N, dn = E.dir === 'a2' ? 'a₂' : 'a₁';
      const leaky = all.filter(q => !q.inGap && q.L >= EDGE_L);
      if (edge.length) {
        edgeYes = true;
        const fl = Math.min(...edge.map(q => q.f)), fh = Math.max(...edge.map(q => q.f)), ml = Math.min(...edge.map(q => q.mu)) / Math.PI, mh = Math.max(...edge.map(q => q.mu)) / Math.PI;
        const nTop = edge.filter(q => q.top >= q.bot).length, side = nTop && nTop < edge.length ? 'both edges of the strip' : nTop ? 'the top edge of the strip' : 'the bottom edge of the strip';
        const avg = edge.reduce((a, q) => a + q.L, 0) / edge.length;
        let cov = '';
        for (const g of gaps) { const em = edge.filter(q => q.f >= g.lo && q.f <= g.hi); if (!em.length) continue; const c = (Math.max(...em.map(q => q.f)) - Math.min(...em.map(q => q.f))) / (g.hi - g.lo);
          cov += c > 0.8 ? ` In the ${fmtF(g.lo)} to ${fmtF(g.hi)} gap the edge branch spans ${Math.round(100 * c)} % of the gap, so this edge guides waves at frequencies the bulk blocks.` : ` In the ${fmtF(g.lo)} to ${fmtF(g.hi)} gap the edge branch covers only ${Math.round(100 * c)} % of the gap: a localised edge resonance rather than a gap-spanning edge channel.`; }
        const inComplete = edge.filter(q => gaps.some(g => q.f >= g.lo && q.f <= g.hi)).length;
        edgeWhere = inComplete ? (inComplete === edge.length ? ', inside the band gap' : ', partly inside the band gap') : (gaps.length ? ', but outside the band gap' : '');
        edgeHtml = `<p><b>Edge waves (${fw}): yes.</b> ${edge.length} strip modes lie outside the projected bulk bands with, on average, ${Math.round(100 * avg)} % of their energy in the outer ${E.w} cells (an evenly spread mode would have ${Math.round(100 * uni)} %). They run from ${fmtF(fl)} to ${fmtF(fh)} for k∥·a from ${ml.toFixed(2)}π to ${mh.toFixed(2)}π, on ${side}, with free edges along ${dn}.${cov}${inComplete < edge.length ? (inComplete ? ` The other ${edge.length - inComplete} sit in directional gaps, where the bulk still propagates in other directions.` : !gaps.length ? ` There is no complete ${fw} gap here, so these edge waves share their frequencies with bulk waves and are kept apart from them only by their wave number along the edge.` : ` None of them lies inside the complete band gap: they sit in directional gaps, below the bulk bands for that edge wave number, which is typical of free-edge waves. Inside the band gap itself the edge reflects waves rather than guiding them.`) : ''}</p>`;
      } else {
        edgeYes = false;
        edgeHtml = `<p><b>Edge waves (${fw}): no.</b> No edge-localised strip mode falls outside the bulk bands for free edges along ${dn} (${E.N}-cell strip).${gaps.length ? ' Waves in the gap are reflected at this edge rather than guided along it.' : ''}${leaky.length ? ` Edge-localised modes do appear inside the bulk bands (${fmtF(Math.min(...leaky.map(q => q.f)))} to ${fmtF(Math.max(...leaky.map(q => q.f)))}), but there they couple to bulk waves and leak away.` : ''}</p>`;
      }
    }
    if (typeof DrawUI !== 'undefined' && S.zak && !S.zak.drawn && S.zak.key === S.bandKey) edgeHtml = (edgeHtml || '') + DrawUI.zakHtml();
    const yn = v => v ? 'yes' : 'no';
    html = `<h3>What this geometry gives</h3><p class="verdict">Out-of-plane gap: ${yn(S.fam.oop.gaps.length)}. In-plane gap: ${yn(S.fam.ip.gaps.length)}. Full gap (both): ${yn(S.complete.length)}. Edge waves${E && E.kind !== 'finite' && !E.drawn ? ` (${FAM_LABEL[E.family].toLowerCase()})` : ''}: ${edgeYes === null ? 'not yet known' : edgeYes ? 'yes' + edgeWhere : 'no'}.</p>${famAssessment('oop')}${famAssessment('ip')}${fullAssessment()}${edgeHtml}<p class="note">Out-of-plane and in-plane waves are independent only for a flat lattice whose members are symmetric about its mid-plane; offsets, curvature or an asymmetric shaker load couple them, so the full gap is the safe criterion for isolating vibration of any kind. Edge waves depend on where the lattice is cut: try the other edge direction or a wider strip before concluding. Topological protection is not assessed here; that needs a Zak-phase or valley-Chern calculation.</p>`;
  }
  $('#assess1').innerHTML = html; if (!isDrawn()) $('#assess2').innerHTML = html; else if (typeof DrawUI !== 'undefined') DrawUI.edgeAssessment();
}
renderAssessment();
function redrawAll() { drawGeometry(); if (!$('#p-disp').classList.contains('hidden')) { drawBand(); drawCplx(); drawMode(); } if (!$('#p-sens').classList.contains('hidden')) drawSweep(); if (!$('#p-edge').classList.contains('hidden')) { drawEdge(); drawEdgeMode(); } if (typeof DrawUI !== 'undefined' && isDrawn()) DrawUI.redraw(); }
// ---------- Brillouin zone and the zone path (Solver panel) ----------
function drawKPath() {
  const cv = $('#cvKpath'); if (!cv || cv.offsetParent === null) return; let Lat, sym = null;
  try { if (isDrawn()) { if (typeof DrawUI === 'undefined') return; Lat = DrawUI.lattice(); } else { const P = readParams(), err = validate(P); if (err) throw new Error(err); Lat = makeModel(P, undefined, undefined, 'oop').Lat; } sym = Lat.sym; }
  catch (e) { const { ctx } = prep(cv); ctx.fillStyle = css('--danger'); ctx.font = '12px system-ui, sans-serif'; wrapText(ctx, e.message, 8, 18, cv.clientWidth - 16, 15); $('#kpNote').textContent = ''; return; }
  const { ctx, w, h } = prep(cv), bz = Lat.bz, pts = Lat.path.pts; let R = 0; for (const p of bz.concat(pts)) R = Math.max(R, Math.hypot(p[0], p[1]));
  const sc = (Math.min(w, h) / 2 - 18) / R, K = p => [w / 2 + p[0] * sc, h / 2 - p[1] * sc];
  ctx.fillStyle = css('--panel-2'); ctx.strokeStyle = css('--ink-2'); ctx.lineWidth = 1; ctx.beginPath(); bz.forEach((p, i) => { const q = K(p); i ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1]); }); ctx.closePath(); ctx.fill(); ctx.stroke();
  ctx.strokeStyle = css('--grid'); ctx.beginPath(); ctx.moveTo(8, h / 2); ctx.lineTo(w - 8, h / 2); ctx.moveTo(w / 2, 8); ctx.lineTo(w / 2, h - 8); ctx.stroke();
  ctx.strokeStyle = css('--gap'); ctx.lineWidth = 2.5; ctx.beginPath(); pts.forEach((p, i) => { const q = K(p); i ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1]); }); ctx.stroke();
  for (let i = 0; i < pts.length - 1; i++) { const a = K(pts[i]), b = K(pts[i + 1]), m = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2], an = Math.atan2(b[1] - a[1], b[0] - a[0]); if (Math.hypot(b[0] - a[0], b[1] - a[1]) < 14) continue;
    ctx.fillStyle = css('--gap'); ctx.beginPath(); ctx.moveTo(m[0] + 6 * Math.cos(an), m[1] + 6 * Math.sin(an)); ctx.lineTo(m[0] - 5 * Math.cos(an) - 4 * Math.sin(an), m[1] - 5 * Math.sin(an) + 4 * Math.cos(an)); ctx.lineTo(m[0] - 5 * Math.cos(an) + 4 * Math.sin(an), m[1] - 5 * Math.sin(an) - 4 * Math.cos(an)); ctx.closePath(); ctx.fill(); }
  ctx.fillStyle = css('--ink'); ctx.font = '600 12.5px system-ui, sans-serif'; const seen = new Set();
  pts.forEach((p, i) => { const q = K(p), key = Math.round(q[0]) + ',' + Math.round(q[1]); ctx.beginPath(); ctx.arc(q[0], q[1], 2.6, 0, TAU); ctx.fill(); if (seen.has(key)) return; seen.add(key); const r = Math.hypot(p[0], p[1]); ctx.fillText(Lat.path.labels[i], q[0] + (r ? 8 * p[0] / r : -14) - 3, q[1] + (r ? -8 * p[1] / r : 12) + 4); });
  const md = Lat.path.mode, symTxt = sym ? (['r60', 'r90', 'mx', 'my'].filter(k => sym[k]).map(k => ({ r60: '6-fold', r90: '4-fold', mx: 'mirror x', my: 'mirror y' })[k]).join(', ') || 'none besides inversion') : '';
  $('#kpNote').textContent = `Path ${Lat.path.labels.join('–')}: ${md === 'std' ? 'standard path' : md === 'custom' ? 'your custom path' : md === 'half' ? 'half-zone boundary' : md === 'quarter' ? 'quarter-zone boundary' : md}.${symTxt ? ` Cell symmetry found: ${symTxt}.` : ''}`;
}
window.MMW = { readParams, validate, buildSpec, S, setStatus, prep, css, hatch, niceTicks, fUnit, mix, fmtF, fmtM, tick, wrapText, isRunning: () => S.running, startRun, endRun, lastLat: () => S.lastLat || 'hexagon', isDrawn, paramSource, renderAssessment, showTab, drawSweep, renderSens, drawEdge, drawKPath };
$('#openEditor').addEventListener('click', () => showTab('geo'));
let rT; window.addEventListener('resize', () => { clearTimeout(rT); rT = setTimeout(redrawAll, 120); });
if (window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { hatchCache = {}; redrawAll(); });
syncVisibility();
const start = (location.hash || '').slice(1);
showTab(tabs.includes(start) ? start : 'geo');
setTimeout(() => runAnalysis(false), 60);
requestAnimationFrame(animate);
})();
