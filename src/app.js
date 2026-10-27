(() => {
const $ = id => document.getElementById(id);
const ELEMENTS = [['H', 0], ['He', 7], ['Li', 0], ['Be', 1], ['B', 2], ['C', 3], ['N', 4], ['O', 5], ['F', 6], ['Ne', 7], ['Cl', 6]];
const ROWS = [['H', null, null, null, null, null, null, 'He'], ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'], [null, null, null, null, null, null, 'Cl', null]];
const PRESETS = [
  ['Second-row series', [['H', 'H', 0], ['He', 'He', 0], ['Li', 'Li', 0], ['Be', 'Be', 0], ['B', 'B', 0], ['C', 'C', 0], ['N', 'N', 0], ['O', 'O', 0], ['F', 'F', 0], ['Ne', 'Ne', 0]]],
  ['Ions', [['H', 'H', 1, 1.06], ['He', 'He', 1, 1.081], ['O', 'O', 1, 1.116], ['O', 'O', -1, 1.35], ['O', 'O', -2, 1.49], ['N', 'O', 1, 1.063], ['C', 'N', -1, 1.16], ['O', 'H', -1, 0.964], ['He', 'H', 1, 0.774]]],
  ['Heteronuclear', [['C', 'O', 0], ['N', 'O', 0], ['B', 'N', 0], ['H', 'F', 0], ['H', 'Cl', 0], ['Li', 'H', 0], ['Cl', 'Cl', 0], ['Cl', 'F', 0]]],
];
const SUPD = { '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶' };
const HEV = 27.211386245988, KJMOL = 2625.4996;
const SPIN_NAMES = { 1: 'Singlet', 2: 'Doublet', 3: 'Triplet', 4: 'Quartet', 5: 'Quintet' };
const state = { cmp: null, act: null, engine: 'eht', A: 'N', B: 'N', q: 0, R: 1.098, mix: 1, K: 1.75, mult: 0, sel: null, view: 'slice', yaw: -0.55, pitch: 0.38, spin: true };
const reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
if (reduceMotion) state.spin = false;
let res = null, eht = null, compareRes = null, cmpReq = 0;

const chargeStr = q => (q === 0 ? '' : (Math.abs(q) > 1 ? Math.abs(q) : '') + (q > 0 ? '+' : '−'));
const formula = (A, B, q) => (A === B ? `${A}<sub>2</sub>` : A + B) + (q ? `<sup>${chargeStr(q)}</sup>` : '');
const moHTML = m => `${m.g}${m.star ? '<sup>*</sup>' : ''}${m.sub ? `<sub>${m.sub}</sub>` : ''}`;
const fmtE = e => (e < 0 ? '−' : '+') + Math.abs(e).toFixed(2);
const signed = v => (v < 0 ? '−' : '+') + Math.abs(v).toFixed(2);
const minus = s => String(s).replace('-', '−');
const nBasis = (A, B) => buildBasis(A, B, 1).length;
const KIND = { b: 'Bonding', a: 'Antibonding', n: 'Nonbonding' };
const isHF = () => state.engine === 'hf';
const hasPx = s => ATOMS[s].shells.some(sh => sh.l === 1);
const totalZ = () => HF.ZNUM[state.A] + HF.ZNUM[state.B] - state.q;
const unique = a => a.filter((v, i) => a.indexOf(v) === i);

// ---------- Hartree–Fock worker ----------
let worker = null, reqId = 0, busy = false, queued = null, curveReq = 0;
const hfAtoms = {};
function makeWorker() {
  try {
    const src = $('hf-src').textContent + '\n' + WORKER_MAIN.toString().replace(/^function WORKER_MAIN\(\) \{/, '').replace(/\}$/, '');
    const w = new Worker(URL.createObjectURL(new Blob([src], { type: 'text/javascript' })));
    w.onmessage = e => onWorker(e.data);
    w.onerror = () => { worker = null; };
    return w;
  } catch (e) { return null; }
}
// body runs inside the worker (HF is defined by the prepended source)
function WORKER_MAIN() {
  const cache = {};
  const atomInfo = s => ({ levels: HF.atomOrbitalEnergies(s), E: HF.atomEnergy(s) });
  const pack = r => { const o = Object.assign({}, r); delete o.Fa; delete o.Fb; return o; };
  self.onmessage = e => {
    const { id, type, p } = e.data, key = [p.A, p.B, p.q, p.mult].join('|');
    try {
      if (type === 'calc' || type === 'calc2') {
        const r = HF.diatomic(p.A, p.B, p.R, p.q, p.mult); // fresh start from atomic densities: lands on the ground state
        if (r.converged && type === 'calc') cache[key] = r;
        self.postMessage({ id, type, p, r: pack(r), atoms: { [p.A]: atomInfo(p.A), [p.B]: atomInfo(p.B) } });
      } else if (type === 'opt') {
        const o = HF.optimize(p.A, p.B, p.R, p.q, p.mult, cache[key] && Math.abs(cache[key].R - p.R) < 1e-6 ? cache[key] : null);
        self.postMessage({ id, type, p, R: o.R });
      } else if (type === 'curve') {
        const n = 34, lo = Math.max(0.45, 0.6 * p.R0), hi = Math.min(4.2, Math.max(2.6, 2.7 * p.R0));
        const Rs = []; for (let i = 0; i < n; i++) Rs.push(lo + (hi - lo) * i / (n - 1));
        const pts = HF.curve(p.A, p.B, p.q, p.mult, Rs, (done, tot) => { if (done % 3 === 0) self.postMessage({ id, type: 'progress', done, n: tot }); });
        const Eatoms = p.q === 0 ? HF.atomEnergy(p.A) + HF.atomEnergy(p.B) : null;
        self.postMessage({ id, type, p, pts, Eatoms });
      }
    } catch (err) { self.postMessage({ id, type, p, error: err.message }); }
  };
}
function hfParams() { return { A: state.A, B: state.B, q: state.q, R: state.R, mult: effMult() }; }
function sameP(a, b) { return a && b && a.A === b.A && a.B === b.B && a.q === b.q && a.mult === b.mult && Math.abs(a.R - b.R) < 1e-9; }
function requestHF() {
  const p = hfParams();
  if (!worker) { // fallback: run on the main thread
    try {
      const r = HF.diatomic(p.A, p.B, p.R, p.q, p.mult);
      const atoms = {}; [p.A, p.B].forEach(s => (atoms[s] = { levels: HF.atomOrbitalEnergies(s), E: HF.atomEnergy(s) }));
      onWorker({ type: 'calc', p, r, atoms });
    } catch (e) { onWorker({ type: 'calc', p, error: e.message }); }
    return;
  }
  if (busy) { queued = p; return; }
  busy = true; setStatus('Calculating…', 'busy');
  worker.postMessage({ id: ++reqId, type: 'calc', p });
}
function onWorker(d) {
  if (d.type === 'calc') {
    busy = false;
    if (queued) { const q = queued; queued = null; if (!sameP(q, d.p)) { busy = true; worker.postMessage({ id: ++reqId, type: 'calc', p: q }); } }
    if (!isHF()) return;
    if (d.error) { setStatus(d.error, 'err'); return; }
    if (d.atoms) Object.assign(hfAtoms, d.atoms);
    if (!sameP(d.p, hfParams())) return; // stale
    res = processHF(d.r, d.p);
    afterResult();
  } else if (d.type === 'calc2') {
    if (d.id !== cmpReq || !isHF()) return;
    if (d.error) { $('cmpStatus').textContent = d.error; return; }
    if (d.atoms) Object.assign(hfAtoms, d.atoms);
    compareRes = processHF(d.r, d.p); compareRes.p = d.p;
    $('cmpStatus').textContent = '';
    drawInto($('diagram'), res, true); drawCompare();
  } else if (d.type === 'opt') {
    $('optBtn').disabled = false; $('optBtn').textContent = 'Optimize bond length';
    if (d.error || !isHF() || d.p.A !== state.A || d.p.B !== state.B || d.p.q !== state.q) return;
    state.R = +d.R.toFixed(4); optNote = `Optimized: R = ${d.R.toFixed(3)} Å (HF/STO-3G minimum).`; update();
  } else if (d.type === 'progress') {
    if (d.id === curveReq) $('curveStatus').textContent = `Computing integrals for ${d.done} of ${d.n} bond lengths…`;
  } else if (d.type === 'curve') {
    if (d.id !== curveReq) return;
    if (d.error) { $('curveStatus').textContent = d.error; return; }
    curve = { p: d.p, pts: d.pts, Eatoms: d.Eatoms };
    drawCurve();
  }
}
let optNote = '';
function setStatus(t, kind) { const el = $('hfStatus'); el.innerHTML = t; el.dataset.kind = kind || ''; }

// ---------- HF result -> diagram-ready structure ----------
function coreCount(sym, block) {
  const Z = HF.ZNUM[sym];
  if (Z <= 2) return 0;
  if (Z <= 10) return block === 'σ' ? 1 : 0;
  return block === 'σ' ? 3 : 1; // Cl: 1s, 2s, 2pz | 2px | 2py
}
function processHF(r, p) {
  const basis = r.basis.map(b => Object.assign({}, b, { kind: 'gto' }));
  const S = r.S, nb = basis.length;
  const all = r.eps.map((e, k) => ({ E: e * HEV, block: r.sym[k], sym: r.sym[k] === 'σ' ? 'σ' : 'π', comp: r.sym[k] === 'πx' ? 'x' : r.sym[k] === 'πy' ? 'y' : '', c: r.C.map(row => row[k]), occ: r.occ[k] }));
  // phase convention: positive at the bond midpoint where possible
  all.forEach(m => {
    const v = gtoMO(basis, m.c, m.comp === 'y' ? 0 : (m.sym === 'π' ? 0.8 : 0), m.comp === 'y' ? 0.8 : 0, r.Rb / 2);
    let flip = Math.abs(v) > 1e-4 ? v < 0 : (() => { let bi = 0, bv = 0; m.c.forEach((c, i) => { if (basis[i].atom === 0 && Math.abs(c) > bv) { bv = Math.abs(c); bi = i; } }); return m.c[bi] < 0; })();
    if (flip) m.c = m.c.map(x => -x);
  });
  // Mulliken analysis
  all.forEach(m => {
    m.ao = basis.map((_, i) => { let g = 0; for (let j = 0; j < nb; j++) g += m.c[i] * m.c[j] * S[i][j]; return g; });
    m.wA = 0; m.wB = 0; m.q = 0;
    basis.forEach((b, i) => { if (b.atom === 0) m.wA += m.ao[i]; else m.wB += m.ao[i]; });
    for (let i = 0; i < nb; i++) for (let j = 0; j < nb; j++) if (basis[i].atom === 0 && basis[j].atom === 1) m.q += 2 * m.c[i] * m.c[j] * S[i][j];
  });
  // split off core orbitals per symmetry block
  const core = [], mos = [];
  for (const block of ['σ', 'πx', 'πy']) {
    const list = all.filter(m => m.block === block).sort((a, b) => a.E - b.E);
    const nc = coreCount(p.A, block) + coreCount(p.B, block);
    list.forEach((m, i) => (i < nc ? core : mos).push(m));
  }
  mos.sort((a, b) => a.E - b.E || (a.comp < b.comp ? -1 : 1));
  labelMOs(mos, basis, p.A, p.B, 1);
  // pair energies for split pi pairs
  // pair energy for drawing: average when both components hold the same number of electrons,
  // otherwise the energy of the more occupied one (e.g. the singly occupied π* of NO)
  mos.forEach(m => { if (m.comp === 'x') { const y = mos.find(o => o.comp === 'y' && o.name === m.name); m.pairE = !y ? m.E : y.occ === m.occ ? (m.E + y.E) / 2 : (m.occ > y.occ ? m.E : y.E); m.split = y ? Math.abs(m.E - y.E) : 0; m.yE = y ? y.E : m.E; m.unevenOcc = y ? y.occ !== m.occ : false; } });
  const unpaired = mos.reduce((s, m) => s + (m.occ === 1 ? 1 : 0), 0);
  let bonding = 0, anti = 0;
  mos.forEach(m => { if (m.kind === 'b') bonding += m.occ; else if (m.kind === 'a') anti += m.occ; });
  let popA = 0, popB = 0;
  all.forEach(m => { popA += m.occ * m.wA; popB += m.occ * m.wB; });
  const occd = mos.filter(m => m.occ > 0).sort((a, b) => a.E - b.E);
  const homo = occd.length ? occd[occd.length - 1] : null;
  const lumo = mos.filter(m => m.occ === 0).sort((a, b) => a.E - b.E)[0] || null;
  const nElec = mos.reduce((s, m) => s + m.occ, 0);
  return { engine: 'hf', symA: p.A, symB: p.B, charge: p.q, R: p.R, Rb: r.Rb, basis, S, mos, core, nElec, unpaired, bondOrder: (bonding - anti) / 2, bonding, anti,
    qA: HF.ZNUM[p.A] - popA, qB: HF.ZNUM[p.B] - popB, homo, lumo, E: r.E, Enuc: r.Enuc, iterations: r.iterations, converged: r.converged, mult: r.mult, history: r.history };
}
function gtoMO(basis, c, x, y, z) {
  let v = 0;
  basis.forEach((b, i) => {
    if (!c[i]) return;
    const dz = z - b.z0, r2 = x * x + y * y + dz * dz;
    let g = 0; for (let k = 0; k < b.exps.length; k++) g += b.coefs[k] * Math.exp(-b.exps[k] * r2);
    v += c[i] * g * (b.type === 's' ? 1 : b.type === 'px' ? x : b.type === 'py' ? y : dz);
  });
  return v;
}

// ---------- controls ----------
function buildPicker(el, which) {
  el.innerHTML = '';
  ROWS.flat().forEach(sym => {
    const b = document.createElement('button');
    b.type = 'button';
    if (!sym) { b.className = 'gap'; b.tabIndex = -1; b.setAttribute('aria-hidden', 'true'); el.appendChild(b); return; }
    b.innerHTML = `<span class="z">${ATOMS[sym].Z}</span>${sym}`;
    b.title = ATOMS[sym].name;
    b.setAttribute('aria-label', `${which === 'A' ? 'Left' : 'Right'} atom: ${ATOMS[sym].name}`);
    b.dataset.sym = sym;
    b.addEventListener('click', () => { state[which] = sym; state.R = defaultBond(state.A, state.B); state.mult = 0; clampCharge(); update(); });
    el.appendChild(b);
  });
}
function chargeOK(q) { const n = ATOMS[state.A].val + ATOMS[state.B].val - q; return n >= 0 && n <= 2 * nBasis(state.A, state.B); }
function clampCharge() { if (!chargeOK(state.q)) state.q = 0; }
function buildCharge() {
  const el = $('charge');
  [-2, -1, 0, 1, 2].forEach(q => {
    const b = document.createElement('button');
    b.type = 'button'; b.textContent = q === 0 ? '0' : (q > 0 ? '+' : '−') + Math.abs(q); b.dataset.q = q;
    b.setAttribute('aria-label', `Charge ${q}`);
    b.addEventListener('click', () => { state.q = q; state.mult = 0; update(); });
    el.appendChild(b);
  });
}
function buildPresets() {
  const el = $('presets');
  PRESETS.forEach(([title, list]) => {
    const g = document.createElement('div'); g.className = 'chip-group';
    g.innerHTML = `<span class="eyebrow">${title}</span>`;
    const wrap = document.createElement('div'); wrap.className = 'chips';
    list.forEach(([A, B, q, R]) => {
      const b = document.createElement('button');
      b.type = 'button'; b.className = 'chip'; b.innerHTML = formula(A, B, q);
      b.dataset.key = A + '|' + B + '|' + q;
      b.addEventListener('click', () => { Object.assign(state, { A, B, q, R: R ?? defaultBond(A, B), mult: 0 }); update(); });
      wrap.appendChild(b);
    });
    g.appendChild(wrap); el.appendChild(g);
  });
}
function buildParamTable() {
  let h = '<thead><tr><th>Atom</th><th>Orbital</th><th>H<sub>ii</sub> (eV)</th><th>ζ</th></tr></thead><tbody>';
  for (const [sym] of ELEMENTS) for (const sh of ATOMS[sym].shells)
    h += `<tr><td>${sym}</td><td>${sh.n}${sh.l ? 'p' : 's'}</td><td>${minus(sh.H.toFixed(1))}</td><td>${sh.z.toFixed(3)}</td></tr>`;
  $('paramTable').innerHTML = h + '</tbody>';
}
// spin multiplicity: 0 = automatic (Hund's rule from the qualitative diagram)
function autoMult() { return (eht ? eht.unpaired : 0) + 1; }
function effMult() {
  const N = totalZ();
  const m = state.mult || autoMult();
  return ((N - (m - 1)) % 2 === 0 && m - 1 <= N) ? m : (N % 2 ? 2 : 1);
}
function buildSpin() {
  const sel = $('spin'), N = totalZ(), start = N % 2 ? 2 : 1;
  const opts = [[0, `Automatic (${SPIN_NAMES[autoMult()] || autoMult()})`]];
  for (let m = start; m <= 5; m += 2) if (m - 1 <= Math.min(N, 2 * nBasis(state.A, state.B))) opts.push([m, `${SPIN_NAMES[m]} (${m - 1} unpaired)`]);
  sel.innerHTML = opts.map(([v, t]) => `<option value="${v}">${t}</option>`).join('');
  sel.value = opts.some(o => o[0] === state.mult) ? state.mult : 0;
}

$('bond').addEventListener('input', e => { state.R = +e.target.value; optNote = ''; schedule(); });
$('bondReset').addEventListener('click', () => { state.R = defaultBond(state.A, state.B); optNote = ''; update(); });
$('mix').addEventListener('input', e => { state.mix = +e.target.value / 100; schedule(); });
$('kval').addEventListener('input', e => { state.K = +e.target.value; schedule(); });
$('kReset').addEventListener('click', () => { state.K = 1.75; update(); });
$('spin').addEventListener('change', e => { state.mult = +e.target.value; update(); });
$('optBtn').addEventListener('click', () => {
  const p = hfParams();
  if (!worker) { const o = HF.optimize(p.A, p.B, p.R, p.q, p.mult); state.R = +o.R.toFixed(4); optNote = `Optimized: R = ${o.R.toFixed(3)} Å (HF/STO-3G minimum).`; update(); return; }
  $('optBtn').disabled = true; $('optBtn').textContent = 'Optimizing…';
  worker.postMessage({ id: ++reqId, type: 'opt', p });
});
document.querySelectorAll('#engineSeg button').forEach(b => b.addEventListener('click', () => { state.engine = b.dataset.engine; optNote = ''; if (!isHF() && state.view === 'curve') state.view = 'slice'; update(); setView(state.view); }));
$('tabSlice').addEventListener('click', () => setView('slice'));
$('tabCloud').addEventListener('click', () => setView('cloud'));
$('tabCurve').addEventListener('click', () => setView('curve'));
$('tabCompare').addEventListener('click', () => setView('compare'));
$('spinBtn').addEventListener('click', () => { state.spin = !state.spin; $('spinBtn').textContent = state.spin ? 'Pause spin' : 'Spin'; if (state.spin) loop(); });
function setView(v) {
  state.view = v;
  [['tabSlice', 'slice'], ['tabCloud', 'cloud'], ['tabCurve', 'curve'], ['tabCompare', 'compare']].forEach(([id, k]) => $(id).setAttribute('aria-selected', v === k));
  const orbital = v === 'slice' || v === 'cloud';
  $('cloudTools').hidden = v !== 'cloud';
  $('orb').classList.toggle('grab', v === 'cloud');
  $('orbWrap').hidden = !orbital;
  $('moDetails').hidden = !orbital;
  $('curveWrap').hidden = v !== 'curve';
  $('cmpWrap').hidden = v !== 'compare';
  $('moEyebrow').textContent = v === 'curve' ? 'Potential energy' : v === 'compare' ? 'Side by side, same energy scale' : 'Selected orbital';
  if (v === 'curve') { $('moTitle').innerHTML = 'Energy curve'; ensureCurve(); }
  else if (v === 'compare') { $('moTitle').innerHTML = 'Compare'; updateCompare(); }
  else { renderInfo(); renderViewer(); renderLegend(); }
  if (res) drawInto($('diagram'), res, true);
  if (v === 'cloud') loop();
}

let pending = false;
function schedule() { if (pending) return; pending = true; requestAnimationFrame(() => { pending = false; update(); }); }

// ---------- main update ----------
function update() {
  eht = calculate(state.A, state.B, state.q, state.R, { spMix: state.mix, K: state.K });
  renderControls();
  if (isHF()) {
    requestHF();
    if (state.view === 'curve') ensureCurve();
    if (state.view === 'compare') updateCompare();
    return; // the drawing updates when the worker answers
  }
  res = eht;
  afterResult();
  if (state.view === 'compare') updateCompare();
}
function afterResult() {
  if (!currentMO()) { const d = defaultMO(); state.sel = d ? { name: d.name } : null; }
  drawDiagram(); renderSummary(); renderNotes(); renderInfo(); renderViewer(); renderLegend(); renderHFStrip();
  if (state.view === 'curve') drawCurve();
}
function defaultMO() {
  const m = res.homo || res.mos[0];
  return m.sym === 'π' ? res.mos.find(x => x.name === m.name && x.comp === 'x') : m;
}
function currentMO() {
  if (!state.sel || !res) return null;
  return res.mos.find(m => m.name === state.sel.name && m.comp !== 'y') || null;
}

function renderControls() {
  ['A', 'B'].forEach(w => $('pick' + w).querySelectorAll('button[data-sym]').forEach(b => b.setAttribute('aria-pressed', b.dataset.sym === state[w])));
  $('charge').querySelectorAll('button').forEach(b => { const q = +b.dataset.q; b.disabled = !chargeOK(q); b.setAttribute('aria-pressed', q === state.q); });
  document.querySelectorAll('#engineSeg button').forEach(b => b.setAttribute('aria-pressed', b.dataset.engine === state.engine));
  const key = state.A + '|' + state.B + '|' + state.q;
  $('presets').querySelectorAll('.chip').forEach(b => b.setAttribute('aria-pressed', b.dataset.key === key));
  $('formula').innerHTML = formula(state.A, state.B, state.q);
  $('formulaMeta').innerHTML = `${ATOMS[state.A].name}${state.A === state.B ? '' : ' + ' + ATOMS[state.B].name}<br>${eht.nElec} valence e<sup>−</sup>`;
  $('bond').value = state.R; $('bondOut').textContent = state.R.toFixed(3) + ' Å';
  const mixable = hasPx(state.A) && hasPx(state.B);
  $('mixBlock').hidden = isHF();
  $('hfBlock').hidden = !isHF();
  $('mix').disabled = !mixable; $('mixBlock').classList.toggle('disabled', !mixable);
  $('mix').value = Math.round(state.mix * 100); $('mixOut').textContent = mixable ? Math.round(state.mix * 100) + '%' : 'n/a';
  $('kval').value = state.K; $('kOut').textContent = state.K.toFixed(2);
  $('tabCurve').hidden = !isHF();
  $('engineHint').innerHTML = isHF()
    ? 'Self-consistent field with the STO-3G minimal basis. Gives total energies, so bond lengths can be optimized.'
    : 'Fast qualitative model. Parameters come from atomic ionization energies.';
  $('optNote').textContent = optNote;
  if (isHF()) buildSpin();
}

function renderHFStrip() {
  const el = $('hfStrip');
  if (!isHF() || res.engine !== 'hf') { el.hidden = true; if (!isHF()) setStatus(`Extended Hückel · K = ${state.K.toFixed(2)}`); return; }
  el.hidden = false;
  const method = res.mult === 1 ? 'RHF' : 'ROHF';
  el.innerHTML = `
    <div><span class="eyebrow">Total energy</span><span class="mono">${minus(res.E.toFixed(6))} E<sub>h</sub></span></div>
    <div><span class="eyebrow">Nuclear repulsion</span><span class="mono">${res.Enuc.toFixed(6)} E<sub>h</sub></span></div>
    <div><span class="eyebrow">Spin state</span><span>${SPIN_NAMES[res.mult] || res.mult} · ${method}</span></div>
    <div><span class="eyebrow">SCF</span><span>${res.converged ? `Converged in ${res.iterations} cycles` : '<b style="color:var(--warn)">Not converged</b>'}</span></div>`;
  setStatus(res.converged ? `HF/STO-3G · ${minus(res.E.toFixed(6))} E<sub>h</sub>` : 'SCF did not converge', res.converged ? '' : 'err');
}

function renderSummary() {
  const bo = res.bondOrder;
  const boTxt = Number.isInteger(bo) ? String(bo) : (Math.floor(bo) ? Math.floor(bo) : '') + '½';
  const para = res.unpaired > 0;
  const homo = res.homo;
  const polar = state.A !== state.B;
  const charges = polar ? `${state.A} ${signed(res.qA)}<br>${state.B} ${signed(res.qB)}` : 'Nonpolar';
  $('stats').innerHTML = `
    <div class="stat"><span class="eyebrow">Bond order</span><span class="v">${boTxt}</span><span class="s">(${res.bonding} − ${res.anti}) ÷ 2</span></div>
    <div class="stat"><span class="eyebrow">Magnetism</span><span class="v" style="font-size:18px;padding-block:3px"><span class="pill" style="color:${para ? 'var(--neg)' : 'var(--muted)'}">${para ? 'Paramagnetic' : 'Diamagnetic'}</span></span><span class="s">${res.unpaired} unpaired e<sup>−</sup></span></div>
    <div class="stat"><span class="eyebrow">HOMO</span><span class="v">${homo ? moHTML(homo) : '—'}</span><span class="s">${homo ? fmtE(homo.E) + ' eV' : 'no electrons'}</span></div>
    <div class="stat"><span class="eyebrow">Polarity</span><span class="v" style="${polar ? 'font-size:14px;line-height:1.35;font-family:var(--f-mono);font-weight:500' : ''}">${charges}</span><span class="s">${polar ? (isHF() ? 'Mulliken charges' : 'Mulliken; EHT exaggerates') : 'Identical atoms'}</span></div>`;
  const parts = [];
  for (const m of res.mos) {
    if (m.comp === 'y') continue;
    let occ = m.occ; if (m.comp === 'x') { const y = res.mos.find(x => x.name === m.name && x.comp === 'y'); if (y) occ += y.occ; }
    if (occ) parts.push(`(${m.g}${m.star}<sub>${m.sub}</sub>)${String(occ).split('').map(c => SUPD[c]).join('')}`);
  }
  $('config').innerHTML = parts.join('') || '<span style="color:var(--muted)">No valence electrons</span>';
}

const TEXTBOOK_PI_FIRST = { B: true, C: true, N: true, O: false, F: false, Ne: false };
function renderNotes() {
  const n = [];
  const A = state.A, B = state.B;
  if (res.bondOrder <= 0) n.push(['warn', `Bond order ${res.bondOrder}: the bonding and antibonding electrons cancel, so ${formula(A, B, state.q)} is not expected to form a stable bond. ${isHF() ? '' : 'The bond length here is only for illustration.'}`]);
  // sigma/pi ordering vs textbooks (homonuclear B2..Ne2)
  if (A === B && A in TEXTBOOK_PI_FIRST && (!(!isHF()) || true)) {
    const s = res.mos.find(m => m.name === 'σ2p'), p = res.mos.find(m => m.name === 'π2p' && m.comp === 'x');
    if (s && p && (isHF() || state.mix > 0.6)) {
      const piFirst = (p.pairE ?? p.E) < s.E, tb = TEXTBOOK_PI_FIRST[A];
      if (piFirst !== tb) n.push(['warn', `Textbooks put σ<sub>2p</sub> ${tb ? 'above' : 'below'} π<sub>2p</sub> for ${A}<sub>2</sub>. ${isHF() ? 'This Hartree–Fock calculation' : 'Extended Hückel'} puts it ${tb ? 'below' : 'above'} (σ<sub>2p</sub> ${fmtE(s.E)} eV, π<sub>2p</sub> ${fmtE(p.pairE ?? p.E)} eV). The bond order and magnetism are the same either way.`]);
    }
  }
  const cap = offCap();
  const off = res.mos.filter(m => (m.comp === 'x' && m.pairE !== undefined ? m.pairE : m.E) > cap && m.comp !== 'y');
  if (off.length) n.push(['', `${unique(off.map(m => moHTML(m))).join(', ')} ${unique(off.map(m => moHTML(m))).length > 1 ? 'lie' : 'lies'} far above the other levels (${unique(off.map(m => fmtE(m.comp === 'x' && m.pairE !== undefined ? m.pairE : m.E))).join(', ')} eV), so ${unique(off.map(m => moHTML(m))).length > 1 ? 'they are' : 'it is'} drawn at the top above the axis break.${isHF() ? '' : ' Extended Hückel strongly destabilizes this antibonding orbital when overlap is large.'}`]);
  if (isHF()) {
    const split = res.mos.filter(m => m.comp === 'x' && m.split > 0.02);
    if (split.length) n.push(['', `The odd π electron sits in π<sub>x</sub>, so π<sub>x</sub> and π<sub>y</sub> end up with different energies in a single-determinant calculation (${split.map(m => `${moHTML(m)}: ${fmtE(m.E)} and ${fmtE(m.yE)} eV`).join('; ')}). The diagram draws each pair at its average, or at the occupied one's energy when only one holds the odd electron.`]);
    if (res.mult > 1) n.push(['', 'Open-shell (ROHF) orbital energies are given as α-spin values, so each occupied level approximates the energy needed to remove that electron (Koopmans’ theorem).']);
    if (res.core.length) {
      const groups = {}; res.core.forEach(m => { const k = fmtE(m.E); groups[k] = (groups[k] || 0) + 1; });
      n.push(['', `Core orbitals are included in the calculation but not drawn: ${res.core.length} orbitals between ${fmtE(Math.min(...res.core.map(m => m.E)))} and ${fmtE(Math.max(...res.core.map(m => m.E)))} eV.`]);
    }
    n.push(['', 'Dashed lines link each molecular orbital to the atomic orbitals that contribute at least 5% of it. Side columns show the free atoms’ orbital energies from the same method.']);
  } else {
    const core = [A, B].some(s => ATOMS[s].Z > 2) ? 'Core electrons (1s for Li–Ne; 1s, 2s, 2p for Cl) are not shown.' : '';
    n.push(['', `Dashed lines link each molecular orbital to the atomic orbitals that contribute at least 5% of it. Each molecule is drawn with neutral atoms on the sides. ${core}`]);
  }
  $('notes').innerHTML = n.map(([c, t]) => `<p class="note ${c}">${t}</p>`).join('');
}

// ---------- MO diagram ----------
function aoLevels(sym) {
  const base = atomLevels(sym);
  if (!isHF() || !hfAtoms[sym]) return base;
  return base.map(l => Object.assign({}, l, { H: hfAtoms[sym].levels[l.label] ?? l.H }));
}
const offCapFor = R => Math.max(...[...aoLevels(R.symA), ...aoLevels(R.symB)].map(l => l.H)) + 9;
const offCap = () => offCapFor(res);
function svgEl(tag, attrs, parent) { const e = document.createElementNS('http://www.w3.org/2000/svg', tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
function arrow(g, x, y, up, color) {
  const d = up ? `M${x},${y + 8} L${x},${y - 8} L${x - 4},${y - 3}` : `M${x},${y - 8} L${x},${y + 8} L${x + 4},${y + 3}`;
  svgEl('path', { d, stroke: color, 'stroke-width': 1.5, fill: 'none', 'stroke-linecap': 'round', 'stroke-linejoin': 'round' }, g);
}
function electrons(g, cx, y, occ, color) {
  if (occ === 2) { arrow(g, cx - 5, y, true, color); arrow(g, cx + 5, y, false, color); }
  else if (occ === 1) arrow(g, cx, y, true, color);
}
function labelText(g, x, y, m, anchor, extra) {
  const t = svgEl('text', { x, y, 'text-anchor': anchor, 'font-family': 'var(--f-body)', 'font-size': 14, fill: 'var(--ink)', 'paint-order': 'stroke', stroke: 'var(--panel)', 'stroke-width': 4, 'stroke-linejoin': 'round' }, g);
  t.appendChild(document.createTextNode(m.g));
  if (m.star) { const s = svgEl('tspan', { dy: -5, 'font-size': 11 }, t); s.textContent = '*'; }
  if (m.sub) { const s = svgEl('tspan', { dy: m.star ? 9 : 4, 'font-size': 10 }, t); s.textContent = m.sub; }
  if (extra) { const s = svgEl('tspan', { dy: -((m.star ? -5 : 0) + (m.sub ? (m.star ? 9 : 4) : 0)), dx: 4, 'font-size': 10.5, fill: 'var(--muted)', 'font-family': 'var(--f-mono)' }, t); s.textContent = extra; }
  return t;
}
function niceStep(span) { for (const s of [1, 2, 5, 10, 20]) if (span / s <= 7) return s; return 50; }
function relax(items, minGap, lo, hi) {
  items.sort((a, b) => a.y - b.y);
  for (let it = 0; it < 60; it++) {
    let moved = false;
    for (let i = 1; i < items.length; i++) { const d = items[i].y - items[i - 1].y; if (d < minGap) { const p = (minGap - d) / 2; items[i].y += p; items[i - 1].y -= p; moved = true; } }
    items.forEach(t => { t.y = Math.max(lo, Math.min(hi, t.y)); });
    if (!moved) break;
  }
}

function drawDiagram() {
  drawInto($('diagram'), res, true);
  if (state.view === 'compare') drawCompare();
}
const groupE = (R, m) => (m.comp === 'x' && m.pairE !== undefined ? m.pairE : m.E);
function levelGroups(R) {
  const groups = [];
  R.mos.forEach((m, i) => {
    if (m.comp === 'y') return;
    groups.push({ m, idx: i, E: groupE(R, m), partner: m.comp === 'x' ? R.mos.find(x => x.name === m.name && x.comp === 'y') : null });
  });
  return groups;
}
function scaleFor(list) {
  // shared energy window for one or more results
  const cap = Math.max(...list.map(R => offCapFor(R)));
  let vals = [];
  list.forEach(R => { vals = vals.concat([...aoLevels(R.symA), ...aoLevels(R.symB)].map(l => l.H), levelGroups(R).filter(g => g.E <= cap).map(g => g.E)); });
  const offMax = Math.max(...list.map(R => levelGroups(R).filter(g => g.E > cap).length));
  return { cap, lo: Math.min(...vals) - 1.5, hi: Math.max(...vals) + 1.5, offMax };
}
function drawInto(svg, R, interactive) {
  svg.innerHTML = '';
  if (!R) return;
  const W = 560, H = 600, top = 26, bot = H - 64;
  const AX = 50, ACX = 124, MCX = 280, BCX = 436;
  const aoA = aoLevels(R.symA), aoB = aoLevels(R.symB);
  const shared = state.view === 'compare' && compareRes && compareRes.engine === R.engine;
  const sc = scaleFor(shared ? [res, compareRes] : [R]);
  const cap = sc.cap, lo = sc.lo, hi = sc.hi;
  const groups = levelGroups(R);
  const on = groups.filter(g => g.E <= cap), off = groups.filter(g => g.E > cap).sort((a, b) => b.E - a.E);
  const yTop = sc.offMax ? top + 26 + sc.offMax * 22 : top;
  const Y = E => bot - (E - lo) / (hi - lo) * (bot - yTop);
  groups.forEach(g => { g.y = g.E > cap ? top + 6 + off.indexOf(g) * 22 : Y(g.E); });
  const onBars = on.map(g => ({ g, y: g.y }));
  relax(onBars, 24, yTop, bot);
  onBars.forEach(b => { b.g.y = b.y; });

  const ax = svgEl('g', { 'font-family': 'var(--f-mono)', 'font-size': 10.5, fill: 'var(--muted)' }, svg);
  svgEl('line', { x1: AX, x2: AX, y1: yTop, y2: bot, stroke: 'var(--muted)', 'stroke-width': 1 }, ax);
  const step = niceStep(hi - lo);
  for (let e = Math.ceil(lo / step) * step; e <= hi; e += step) {
    const y = Y(e);
    svgEl('line', { x1: AX - 4, x2: AX, y1: y, y2: y, stroke: 'var(--muted)' }, ax);
    const t = svgEl('text', { x: AX - 7, y: y + 3.5, 'text-anchor': 'end' }, ax); t.textContent = e === 0 ? '0' : (e < 0 ? '−' : '+') + Math.abs(e);
  }
  const at = svgEl('text', { x: 12, y: (yTop + bot) / 2, transform: `rotate(-90 12 ${(yTop + bot) / 2})`, 'text-anchor': 'middle', 'font-family': 'var(--f-display)', 'font-size': 11.5, 'letter-spacing': '0.08em', fill: 'var(--muted)' }, ax);
  at.textContent = 'ENERGY (eV)';
  if (sc.offMax) {
    svgEl('line', { x1: AX, x2: AX, y1: top - 4, y2: yTop - 16, stroke: 'var(--muted)', 'stroke-dasharray': '2 3' }, ax);
    svgEl('path', { d: `M${AX - 6},${yTop - 12} l12,-4 M${AX - 6},${yTop - 6} l12,-4`, stroke: 'var(--muted)', 'stroke-width': 1.2 }, ax);
  }

  const aoGeom = (cx, lv) => lv.deg === 1 ? [{ x: cx, w: 30 }] : [-25, 0, 25].map(o => ({ x: cx + o, w: 20 }));
  const aoLabels = [];
  const drawAtom = (lvls, cx, side) => {
    const g = svgEl('g', {}, svg);
    lvls.forEach(lv => {
      const y = Y(lv.H);
      aoGeom(cx, lv).forEach((b, k) => {
        svgEl('line', { x1: b.x - b.w / 2, x2: b.x + b.w / 2, y1: y, y2: y, stroke: 'var(--muted)', 'stroke-width': 2.5, 'stroke-linecap': 'round' }, g);
        electrons(g, b.x, y, lv.occ[k], 'var(--muted)');
      });
      aoLabels.push({ y, side, cx, text: lv.label });
    });
  };
  const lines = svgEl('g', { fill: 'none', 'stroke-dasharray': '3 3', stroke: 'var(--muted)' }, svg);
  groups.forEach(g => {
    const m = g.m, half = m.sym === 'π' ? 49 : 27;
    [[0, aoA, ACX, -1], [1, aoB, BCX, 1]].forEach(([atom, lvls, cx, dir]) => {
      lvls.forEach(lv => {
        const w = R.basis.reduce((t, b, i) => t + (b.atom === atom && b.shell === lv.label ? m.ao[i] : 0), 0);
        if (w < 0.05) return;
        const xa = cx - dir * (lv.deg === 1 ? 17 : 37), xm = MCX + dir * (half + 2);
        svgEl('line', { x1: xa, y1: Y(lv.H), x2: xm, y2: g.y, 'stroke-opacity': (0.2 + 0.7 * Math.min(1, w)).toFixed(2), 'stroke-width': 1 }, lines);
      });
    });
  });
  drawAtom(aoA, ACX, -1); drawAtom(aoB, BCX, 1);
  const al = svgEl('g', { 'font-family': 'var(--f-body)', 'font-size': 12, fill: 'var(--muted)' }, svg);
  aoLabels.forEach(it => { const t = svgEl('text', { x: it.cx + it.side * 44, y: it.y + 4, 'text-anchor': it.side < 0 ? 'end' : 'start' }, al); t.textContent = it.text; });

  const labels = groups.map(g => ({ g, y: g.y }));
  relax(labels, 16, top - 4, bot + 6);
  const sel = interactive ? currentMO() : null;
  groups.forEach(gr => {
    const m = gr.m;
    const g = svgEl('g', interactive ? { class: 'lvl' + (sel && sel.name === m.name ? ' sel' : ''), tabindex: 0, role: 'button', 'data-idx': gr.idx,
      'aria-label': `${m.name}${m.comp ? ' (degenerate pair)' : ''}, ${gr.E.toFixed(2)} electron volts, ${m.occ + (gr.partner ? gr.partner.occ : 0)} electrons` } : {}, svg);
    if (interactive) svgEl('rect', { class: 'band', x: 210, y: gr.y - 11, width: 190, height: 22, rx: 4, fill: 'transparent' }, g);
    const bars = m.sym === 'π' ? [{ x: MCX - 27, o: m.occ }, { x: MCX + 27, o: gr.partner ? gr.partner.occ : 0 }] : [{ x: MCX, o: m.occ }];
    const bw = m.sym === 'π' ? 44 : 54;
    bars.forEach(b => {
      svgEl('line', { class: 'bar', x1: b.x - bw / 2, x2: b.x + bw / 2, y1: gr.y, y2: gr.y, stroke: 'var(--ink)', 'stroke-width': 3, 'stroke-linecap': 'round' }, g);
      electrons(g, b.x, gr.y, b.o, 'var(--ink)');
    });
    const ly = labels.find(l => l.g === gr).y;
    labelText(g, 336, ly + 4.5, m, 'start', gr.E > cap ? '↑ ' + fmtE(gr.E) + ' eV' : '');
    const tag = R.homo && m.name === R.homo.name ? 'HOMO' : R.lumo && m.name === R.lumo.name ? 'LUMO' : '';
    if (tag) {
      const t = svgEl('text', { x: MCX - (m.sym === 'π' ? 54 : 32), y: gr.y + 3.5, 'text-anchor': 'end', 'font-family': 'var(--f-display)', 'font-weight': 600, 'font-size': 10, 'letter-spacing': '0.06em', fill: 'var(--accent)', 'paint-order': 'stroke', stroke: 'var(--panel)', 'stroke-width': 4 }, g);
      t.textContent = tag;
    }
    if (!interactive) return;
    g.addEventListener('click', () => selectMO(m));
    g.addEventListener('keydown', e => {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); selectMO(m); }
      if (e.key === 'ArrowUp' || e.key === 'ArrowDown') {
        e.preventDefault();
        const order = groups.slice().sort((a, b) => a.E - b.E);
        const k = order.indexOf(gr) + (e.key === 'ArrowUp' ? 1 : -1);
        if (order[k]) { selectMO(order[k].m); const nx = svg.querySelector(`.lvl[data-idx="${order[k].idx}"]`); nx && nx.focus(); }
      }
    });
  });

  const hd = svgEl('g', { 'text-anchor': 'middle' }, svg);
  const head = (x, big, small) => {
    const t = svgEl('text', { x, y: H - 30, 'font-family': 'var(--f-display)', 'font-weight': 700, 'font-size': 20, fill: 'var(--ink)' }, hd); t.textContent = big;
    const s = svgEl('text', { x, y: H - 12, 'font-family': 'var(--f-display)', 'font-size': 11, 'letter-spacing': '0.08em', fill: 'var(--muted)' }, hd); s.textContent = small;
  };
  head(ACX, R.symA, 'ATOMIC ORBITALS');
  head(BCX, R.symB, 'ATOMIC ORBITALS');
  head(MCX + 20, plainFormula(R.symA, R.symB, R.charge), 'MOLECULAR ORBITALS');
}
const plainFormula = (A, B, q) => (A === B ? A + '₂' : A + B) + (q ? ({ 1: '⁺', 2: '²⁺', '-1': '⁻', '-2': '²⁻' })[q] : '');

function selectMO(m) {
  state.sel = { name: m.name };
  drawDiagram(); renderInfo(); renderViewer(); renderLegend();
}


// ---------- compare view ----------
function defaultCmp() {
  return state.A === 'N' && state.B === 'N' && state.q === 0 ? { A: 'O', B: 'O', q: 0, R: 1.208 } : { A: 'N', B: 'N', q: 0, R: 1.098 };
}
function buildCmpChips() {
  const el = $('cmpChips');
  el.innerHTML = PRESETS.flatMap(([, list]) => list).map(([A, B, q, R]) =>
    `<button type="button" class="chip" data-a="${A}" data-b="${B}" data-q="${q}" data-r="${R ?? defaultBond(A, B)}">${formula(A, B, q)}</button>`).join('');
  el.querySelectorAll('.chip').forEach(b => b.addEventListener('click', () => {
    state.cmp = { A: b.dataset.a, B: b.dataset.b, q: +b.dataset.q, R: +b.dataset.r }; updateCompare();
  }));
  $('cmpSwap').addEventListener('click', () => {
    if (!state.cmp) return;
    const c = state.cmp; state.cmp = { A: state.A, B: state.B, q: state.q, R: state.R };
    Object.assign(state, { A: c.A, B: c.B, q: c.q, R: c.R, mult: 0, sel: null });
    update();
  });
}
function updateCompare() {
  if (state.view !== 'compare') return;
  if (!state.cmp) state.cmp = defaultCmp();
  const c = state.cmp;
  $('cmpChips').querySelectorAll('.chip').forEach(b => b.setAttribute('aria-pressed', b.dataset.a === c.A && b.dataset.b === c.B && +b.dataset.q === c.q));
  const e2 = calculate(c.A, c.B, c.q, c.R, { spMix: state.mix, K: state.K });
  if (!isHF()) { compareRes = e2; $('cmpStatus').textContent = ''; drawInto($('diagram'), res, true); drawCompare(); return; }
  const p = { A: c.A, B: c.B, q: c.q, R: c.R, mult: e2.unpaired + 1 };
  if (compareRes && compareRes.engine === 'hf' && sameP(compareRes.p, p)) { drawCompare(); return; }
  if (!worker) {
    const r = HF.diatomic(p.A, p.B, p.R, p.q, p.mult);
    [p.A, p.B].forEach(s => (hfAtoms[s] = { levels: HF.atomOrbitalEnergies(s), E: HF.atomEnergy(s) }));
    compareRes = processHF(r, p); compareRes.p = p; drawInto($('diagram'), res, true); drawCompare(); return;
  }
  $('cmpStatus').textContent = 'Calculating…';
  cmpReq = ++reqId; worker.postMessage({ id: cmpReq, type: 'calc2', p });
}
function configHTML(R) {
  const parts = [];
  for (const m of R.mos) {
    if (m.comp === 'y') continue;
    let occ = m.occ; if (m.comp === 'x') { const y = R.mos.find(x => x.name === m.name && x.comp === 'y'); if (y) occ += y.occ; }
    if (occ) parts.push(`(${m.g}${m.star}<sub>${m.sub}</sub>)${String(occ).split('').map(c => SUPD[c]).join('')}`);
  }
  return parts.join('') || '—';
}
function drawCompare() {
  const ok = compareRes && res && compareRes.engine === res.engine;
  drawInto($('cmpDiagram'), ok ? compareRes : null, false);
  if (!ok) { $('cmpTable').innerHTML = ''; return; }
  const col = R => {
    const bo = R.bondOrder, boTxt = Number.isInteger(bo) ? String(bo) : (Math.floor(bo) ? Math.floor(bo) : '') + '½';
    return { f: formula(R.symA, R.symB, R.charge), bo: boTxt, mag: R.unpaired ? `Paramagnetic (${R.unpaired} unpaired)` : 'Diamagnetic',
      homo: R.homo ? `${moHTML(R.homo)} ${fmtE(groupE(R, R.homo))} eV` : '—', cfg: configHTML(R), R: R.R.toFixed(3) + ' Å',
      E: R.engine === 'hf' ? minus(R.E.toFixed(5)) + ' E<sub>h</sub>' : null,
      qq: R.symA === R.symB ? 'Nonpolar' : `${R.symA} ${signed(R.qA)}, ${R.symB} ${signed(R.qB)}` };
  };
  const a = col(res), b = col(compareRes);
  const rows = [['Molecule', 'f'], ['Bond order', 'bo'], ['Magnetism', 'mag'], ['HOMO', 'homo'], ['Configuration', 'cfg'], ['Partial charges', 'qq'], ['Bond length used', 'R']];
  if (a.E) rows.push(['Total energy', 'E']);
  $('cmpTable').innerHTML = `<table class="cmp"><tbody>${rows.map(([t, k]) => `<tr><th scope="row">${t}</th><td>${a[k]}</td><td>${b[k]}</td></tr>`).join('')}</tbody></table>`;
}

// ---------- MO information ----------
function renderInfo() {
  const m = currentMO(); if (!m || state.view === 'curve' || state.view === 'compare') return;
  const pair = m.comp === 'x' ? res.mos.find(x => x.name === m.name && x.comp === 'y') : null;
  $('moTitle').innerHTML = moHTML(m) + (pair ? '<span style="font-size:14px;font-weight:500;color:var(--muted);font-family:var(--f-body)"> &nbsp;π<sub>x</sub> of a degenerate pair</span>' : '');
  const occ = m.occ + (pair ? pair.occ : 0);
  const homonuc = state.A === state.B;
  const share = homonuc ? 'Shared equally by both atoms' : `${Math.round(Math.max(0, m.wA) * 100)}% on ${state.A}, ${Math.round(Math.max(0, m.wB) * 100)}% on ${state.B}`;
  const eTxt = pair && Math.abs(pair.E - m.E) > 0.02 ? `E = ${fmtE(m.E)} (x), ${fmtE(pair.E)} (y) eV` : `E = ${fmtE(m.E)} eV`;
  $('moFacts').innerHTML = `<span><b>${KIND[m.kind]}</b></span><span class="mono">${eTxt}</span><span>${occ} e<sup>−</sup> ${pair ? 'in the pair' : ''}</span><span>${share}</span>`;
  const rows = res.basis.map((b, i) => ({ b, g: m.ao[i], c: m.c[i] })).filter(r => Math.abs(r.c) > (isHF() ? 0.02 : 1e-4) && r.b.type !== 'py');
  $('comp').innerHTML = `<div class="comp-row comp-head"><span>Atomic orbital</span><span>Share of this MO</span><span>%</span><span>Coeff.</span></div>` +
    rows.map(r => {
      const pct = Math.max(0, r.g) * 100;
      const lab = `${r.b.sym} ${r.b.shell[0]}${r.b.type === 's' ? 's' : 'p<sub>' + r.b.type[1] + '</sub>'}`;
      return `<div class="comp-row"><span class="lab">${lab}</span><span class="track"><span class="fill" style="display:block;width:${Math.min(100, pct).toFixed(1)}%;background:var(${r.c >= 0 ? '--pos' : '--neg'})"></span></span><span class="pct">${pct.toFixed(0)}</span><span class="cf">${signed(r.c)}</span></div>`;
    }).join('');
}

function renderLegend() {
  const m = currentMO();
  const pi = m && m.sym === 'π';
  $('legend').innerHTML = state.view === 'slice'
    ? `<span><span class="sw pos"></span>Positive phase</span><span><span class="sw neg"></span>Negative phase</span><span><span class="sw node"></span>Node</span><span>Slice through the bond axis${pi ? ' (π<sub>y</sub> is the same, rotated 90°)' : ''}</span>`
    : `<span><span class="sw pos"></span>Positive phase</span><span><span class="sw neg"></span>Negative phase</span><span>Dot density follows ψ². Drag to rotate.</span>`;
}

// ---------- orbital rendering ----------
let colors = {};
function readColors() {
  const cs = getComputedStyle(document.documentElement);
  const hex = v => { v = cs.getPropertyValue(v).trim(); if (v.length === 4) v = '#' + v[1] + v[1] + v[2] + v[2] + v[3] + v[3]; return [1, 3, 5].map(i => parseInt(v.slice(i, i + 2), 16)); };
  colors = { bg: hex('--bg'), ink: hex('--ink'), muted: hex('--muted'), pos: hex('--pos'), neg: hex('--neg'), panel: hex('--panel') };
}
function makeEval(m) {
  if (res.engine === 'hf') {
    const fs = res.basis.map((b, i) => ({ b, c: m.c[i] })).filter(f => Math.abs(f.c) > 1e-6).map(({ b, c }) => ({ z0: b.z0, t: b.type, e: b.exps, k: b.coefs.map(x => x * c) }));
    return (x, y, z) => {
      let v = 0;
      for (let i = 0; i < fs.length; i++) {
        const f = fs[i], dz = z - f.z0, r2 = x * x + y * y + dz * dz;
        let g = 0; for (let k = 0; k < f.e.length; k++) g += f.k[k] * Math.exp(-f.e[k] * r2);
        v += f.t === 's' ? g : g * (f.t === 'pz' ? dz : f.t === 'px' ? x : y);
      }
      return v;
    };
  }
  const fs = res.basis.map((b, i) => ({ b, c: m.c[i] })).filter(f => Math.abs(f.c) > 1e-6).map(({ b, c }) => ({
    z0: b.z0, n1: b.n - 1, zeta: b.zeta, k: c * b.N * (b.type === 's' ? 0.28209479177 : 0.48860251190), t: b.type,
  }));
  return (x, y, z) => {
    let v = 0;
    for (let i = 0; i < fs.length; i++) {
      const f = fs[i], dz = z - f.z0, r = Math.sqrt(x * x + y * y + dz * dz);
      let rad = f.k * Math.exp(-f.zeta * r);
      if (f.n1 === 1) rad *= r; else if (f.n1 === 2) rad *= r * r;
      if (f.t === 's') v += rad;
      else if (r > 1e-9) v += rad * (f.t === 'pz' ? dz : f.t === 'px' ? x : y) / r;
    }
    return v;
  };
}
function extent() {
  let L = 3;
  for (const s of [res.symA, res.symB]) for (const sh of ATOMS[s].shells) L = Math.max(L, (sh.n + 2.2) / sh.z);
  return Math.min(L, 10);
}
function frame(cw, ch) {
  const L = extent(), Rb = res.Rb;
  const s = Math.max((Rb + 2 * L) / cw, (2 * L) / ch);
  return { s, cz: Rb / 2, L };
}
const canvas = $('orb'), ctx = canvas.getContext('2d');
function sizeCanvas() {
  const r = canvas.getBoundingClientRect();
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  const w = Math.max(200, Math.round(r.width)), h = Math.round(w * 0.75);
  if (canvas.width !== Math.round(w * dpr)) { canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); }
  return { w, h, dpr };
}
function renderViewer() {
  if (!res || state.view === 'curve') return;
  if (state.view === 'slice') renderSlice(); else { prepCloud(); drawCloud(); }
}
const offc = document.createElement('canvas');
function renderSlice() {
  const m = currentMO(); if (!m) return;
  const { w, h, dpr } = sizeCanvas();
  const f = frame(w, h), ev = makeEval(m);
  offc.width = w; offc.height = h;
  const octx = offc.getContext('2d'), img = octx.createImageData(w, h), d = img.data;
  const vals = new Float32Array(w * h);
  let mx = 1e-12;
  for (let j = 0; j < h; j++) {
    const x = (h / 2 - j - 0.5) * f.s;
    for (let i = 0; i < w; i++) {
      const z = f.cz + (i + 0.5 - w / 2) * f.s;
      const v = ev(x, 0, z); vals[j * w + i] = v; const a = Math.abs(v); if (a > mx) mx = a;
    }
  }
  const { bg, ink, pos, neg } = colors;
  const levels = [0.06, 0.15, 0.3, 0.5, 0.75];
  for (let j = 0; j < h; j++) for (let i = 0; i < w; i++) {
    const k = j * w + i, v = vals[k] / mx, a = Math.abs(v);
    const t = Math.min(1, Math.pow(a, 0.55)) * 0.9;
    const pc = v >= 0 ? pos : neg;
    let r = bg[0] + (pc[0] - bg[0]) * t, g = bg[1] + (pc[1] - bg[1]) * t, b = bg[2] + (pc[2] - bg[2]) * t;
    const vr = i < w - 1 ? vals[k + 1] / mx : v, vd = j < h - 1 ? vals[k + w] / mx : v;
    let line = 0;
    if ((v > 0) !== (vr > 0) && Math.max(a, Math.abs(vr)) > 0.004 || (v > 0) !== (vd > 0) && Math.max(a, Math.abs(vd)) > 0.004) {
      if (((i + j) >> 2) % 2 === 0) line = 0.75;
    } else {
      const ar = Math.abs(vr), ad = Math.abs(vd);
      for (const L of levels) if ((a - L) * (ar - L) < 0 || (a - L) * (ad - L) < 0) { line = 0.28; break; }
    }
    if (line) { r += (ink[0] - r) * line; g += (ink[1] - g) * line; b += (ink[2] - b) * line; }
    d[4 * k] = r; d[4 * k + 1] = g; d[4 * k + 2] = b; d[4 * k + 3] = 255;
  }
  octx.putImageData(img, 0, 0);
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(offc, 0, 0, canvas.width, canvas.height);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const zx = z => w / 2 + (z - f.cz) / f.s;
  ctx.strokeStyle = rgba(colors.muted, 0.5); ctx.setLineDash([2, 4]); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(0, h / 2); ctx.lineTo(w, h / 2); ctx.stroke(); ctx.setLineDash([]);
  [[res.symA, 0], [res.symB, res.Rb]].forEach(([sym, z]) => nucleus(zx(z), h / 2, sym));
  scaleBar(w, h, f.s);
}
function rgba(c, a) { return `rgba(${c[0]},${c[1]},${c[2]},${a})`; }
function nucleus(x, y, sym) {
  ctx.fillStyle = rgba(colors.ink, 1); ctx.beginPath(); ctx.arc(x, y, 3, 0, 7); ctx.fill();
  ctx.font = '600 13px "IBM Plex Sans Condensed", "Arial Narrow", sans-serif'; ctx.textAlign = 'center';
  ctx.lineWidth = 3; ctx.strokeStyle = rgba(colors.panel, 0.85); ctx.strokeText(sym, x, y + 18); ctx.fillText(sym, x, y + 18);
}
function scaleBar(w, h, s) {
  const px = BOHR / s;
  ctx.strokeStyle = rgba(colors.ink, 0.8); ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(14, h - 16); ctx.lineTo(14 + px, h - 16); ctx.moveTo(14, h - 20); ctx.lineTo(14, h - 12); ctx.moveTo(14 + px, h - 20); ctx.lineTo(14 + px, h - 12); ctx.stroke();
  ctx.fillStyle = rgba(colors.ink, 0.85); ctx.font = '500 11px "IBM Plex Mono", monospace'; ctx.textAlign = 'left';
  ctx.fillText('1 Å', 18 + px, h - 12);
}

// 3D cloud
let cloud = null, cloudKey = '';
function rng(seed) { return () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
function prepCloud() {
  const m = currentMO(); if (!m) return;
  const key = [res.engine, res.symA, res.symB, res.R, res.charge, state.mix, state.K, res.mult, m.name].join('|');
  if (key === cloudKey) return;
  cloudKey = key;
  const ev = makeEval(m), L = extent(), Rb = res.Rb;
  let mx = 1e-12;
  for (let i = 0; i <= 160; i++) for (let j = 0; j <= 80; j++) { const v = Math.abs(ev(-L + 2 * L * j / 80, 0, -L + (Rb + 2 * L) * i / 160)); if (v > mx) mx = v; }
  const target = 9000, pts = new Float32Array(target * 4), rnd = rng(12345);
  let n = 0, tries = 0; const m2 = mx * mx * 1.05;
  while (n < target && tries < 2.5e6) {
    tries++;
    const x = (rnd() * 2 - 1) * L, y = (rnd() * 2 - 1) * L, z = -L + rnd() * (Rb + 2 * L);
    const v = ev(x, y, z);
    if (rnd() * m2 < v * v) { pts[4 * n] = x; pts[4 * n + 1] = y; pts[4 * n + 2] = z - Rb / 2; pts[4 * n + 3] = v >= 0 ? 1 : -1; n++; }
  }
  cloud = { pts, n, L, Rb };
}
function drawCloud() {
  if (!cloud || state.view !== 'cloud') return;
  const { w, h, dpr } = sizeCanvas();
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.fillStyle = rgba(colors.bg, 1); ctx.fillRect(0, 0, w, h);
  const f = frame(w, h), sc = 1 / f.s;
  const cy = Math.cos(state.yaw), sy = Math.sin(state.yaw), cp = Math.cos(state.pitch), sp = Math.sin(state.pitch);
  const proj = (X, Y, Z) => { const X1 = X * cy + Z * sy, Z1 = -X * sy + Z * cy; const Y1 = Y * cp - Z1 * sp, Z2 = Y * sp + Z1 * cp; return [w / 2 + X1 * sc, h / 2 - Y1 * sc, Z2]; };
  const a0 = proj(-cloud.Rb / 2 - cloud.L, 0, 0), a1 = proj(cloud.Rb / 2 + cloud.L, 0, 0);
  ctx.strokeStyle = rgba(colors.muted, 0.45); ctx.setLineDash([2, 4]); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(a0[0], a0[1]); ctx.lineTo(a1[0], a1[1]); ctx.stroke(); ctx.setLineDash([]);
  const P = cloud.pts, sz = 1.7;
  for (const sgn of [1, -1]) {
    const c = sgn > 0 ? colors.pos : colors.neg;
    ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`;
    for (let i = 0; i < cloud.n; i++) {
      if (P[4 * i + 3] !== sgn) continue;
      const [u, v, dz] = proj(P[4 * i + 2], P[4 * i], P[4 * i + 1]);
      ctx.globalAlpha = Math.max(0.18, Math.min(0.85, 0.52 + dz * 0.06 / (cloud.L / 4)));
      ctx.fillRect(u - sz / 2, v - sz / 2, sz, sz);
    }
  }
  ctx.globalAlpha = 1;
  [[res.symA, -cloud.Rb / 2], [res.symB, cloud.Rb / 2]].forEach(([sym, X]) => { const p = proj(X, 0, 0); nucleus(p[0], p[1], sym); });
  scaleBar(w, h, f.s);
}
let looping = false;
function loop() {
  if (looping) return; looping = true;
  const tick = () => {
    if (state.view !== 'cloud' || !state.spin || reduceMotion) { looping = false; return; }
    state.yaw += 0.006; drawCloud(); requestAnimationFrame(tick);
  };
  requestAnimationFrame(tick);
}
let drag = null;
canvas.addEventListener('pointerdown', e => { if (state.view !== 'cloud') return; drag = { x: e.clientX, y: e.clientY, yaw: state.yaw, pitch: state.pitch }; canvas.setPointerCapture(e.pointerId); if (state.spin) { state.spin = false; $('spinBtn').textContent = 'Spin'; } });
canvas.addEventListener('pointermove', e => { if (!drag) return; state.yaw = drag.yaw + (e.clientX - drag.x) * 0.01; state.pitch = Math.max(-1.5, Math.min(1.5, drag.pitch + (e.clientY - drag.y) * 0.01)); drawCloud(); });
canvas.addEventListener('pointerup', () => { drag = null; });
canvas.addEventListener('pointercancel', () => { drag = null; });

// ---------- potential energy curve (HF only) ----------
let curve = null;
function curveKey() { return [state.A, state.B, state.q, effMult()].join('|'); }
function ensureCurve() {
  if (!isHF()) return;
  const key = curveKey();
  if (curve && [curve.p.A, curve.p.B, curve.p.q, curve.p.mult].join('|') === key) { drawCurve(); return; }
  if (curve && curve.pendingKey === key) return;
  const p = { A: state.A, B: state.B, q: state.q, mult: effMult(), R0: defaultBond(state.A, state.B) };
  curve = { pendingKey: key, p, pts: null };
  $('curveSvg').innerHTML = ''; $('curveStatus').textContent = 'Computing the energy curve…';
  if (!worker) { $('curveStatus').textContent = 'Background calculation is not available in this browser.'; return; }
  curveReq = ++reqId;
  worker.postMessage({ id: curveReq, type: 'curve', p });
}
function drawCurve() {
  const svg = $('curveSvg');
  if (!curve || !curve.pts) return;
  const pts = curve.pts.filter(p => p.ok);
  const W = 520, H = 360, L = 58, Rr = 16, T = 16, B = 44;
  const ref = curve.Eatoms ?? Math.min(...pts.map(p => p.E));
  const ys = pts.map(p => (p.E - ref) * KJMOL);
  const emin = Math.min(...ys), far = ys[ys.length - 1];
  let ymax = Math.max(far, 0) + Math.max(150, 0.5 * (Math.max(far, 0) - emin));
  let ymin = emin - Math.max(60, 0.1 * (ymax - emin));
  const xs = pts.map(p => p.R), xmin = xs[0], xmax = xs[xs.length - 1];
  const X = r => L + (r - xmin) / (xmax - xmin) * (W - L - Rr), Y = e => T + (ymax - e) / (ymax - ymin) * (H - T - B);
  let h = '';
  const ystep = niceStep2((ymax - ymin) / 6);
  for (let e = Math.ceil(ymin / ystep) * ystep; e <= ymax; e += ystep) h += `<line x1="${L}" x2="${W - Rr}" y1="${Y(e)}" y2="${Y(e)}" stroke="var(--line)"/><text x="${L - 6}" y="${Y(e) + 3.5}" text-anchor="end" font-family="var(--f-mono)" font-size="10.5" fill="var(--muted)">${minus(Math.round(e))}</text>`;
  for (let r = Math.ceil(xmin * 2) / 2; r <= xmax; r += 0.5) h += `<line x1="${X(r)}" x2="${X(r)}" y1="${H - B}" y2="${H - B + 4}" stroke="var(--muted)"/><text x="${X(r)}" y="${H - B + 16}" text-anchor="middle" font-family="var(--f-mono)" font-size="10.5" fill="var(--muted)">${r.toFixed(1)}</text>`;
  h += `<text x="${(L + W - Rr) / 2}" y="${H - 8}" text-anchor="middle" font-family="var(--f-display)" font-size="11.5" letter-spacing="0.08em" fill="var(--muted)">BOND LENGTH (Å)</text>`;
  h += `<text x="14" y="${(T + H - B) / 2}" transform="rotate(-90 14 ${(T + H - B) / 2})" text-anchor="middle" font-family="var(--f-display)" font-size="11.5" letter-spacing="0.08em" fill="var(--muted)">${curve.Eatoms != null ? 'ENERGY vs ATOMS (kJ/mol)' : 'ENERGY vs MINIMUM (kJ/mol)'}</text>`;
  if (curve.Eatoms != null) h += `<line x1="${L}" x2="${W - Rr}" y1="${Y(0)}" y2="${Y(0)}" stroke="var(--ink)" stroke-dasharray="4 4" stroke-opacity="0.6"/><text x="${W - Rr - 4}" y="${Y(0) - 6}" text-anchor="end" font-family="var(--f-body)" font-size="11.5" fill="var(--muted)">separated atoms</text>`;
  h += `<clipPath id="cclip"><rect x="${L}" y="${T}" width="${W - L - Rr}" height="${H - T - B}"/></clipPath>`;
  const path = pts.map((p, i) => `${i ? 'L' : 'M'}${X(p.R).toFixed(1)},${Y(ys[i]).toFixed(1)}`).join('');
  h += `<path d="${path}" fill="none" stroke="var(--accent)" stroke-width="2.25" clip-path="url(#cclip)"/>`;
  pts.forEach((p, i) => { if (ys[i] <= ymax) h += `<circle cx="${X(p.R)}" cy="${Y(ys[i])}" r="2.2" fill="var(--accent)"/>`; });
  // minimum by parabola through the lowest three points
  let k = ys.indexOf(emin); let Rmin = xs[k], Emin = emin;
  if (k > 0 && k < ys.length - 1) { const [x0, x1, x2] = [xs[k - 1], xs[k], xs[k + 1]], [y0, y1, y2] = [ys[k - 1], ys[k], ys[k + 1]]; const den = (x0 - x1) * (x0 - x2) * (x1 - x2); const a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den, b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / den; if (a > 0) { Rmin = -b / (2 * a); Emin = y1 + a * (Rmin - x1) ** 2 + 0 * b; } }
  h += `<line x1="${X(Rmin)}" x2="${X(Rmin)}" y1="${Y(emin)}" y2="${H - B}" stroke="var(--muted)" stroke-dasharray="2 3"/>`;
  // current geometry marker
  if (res && res.engine === 'hf' && state.R >= xmin && state.R <= xmax) {
    const ec = (res.E - ref) * KJMOL;
    h += `<line x1="${X(state.R)}" x2="${X(state.R)}" y1="${T}" y2="${H - B}" stroke="var(--neg)" stroke-opacity="0.6"/>`;
    if (ec <= ymax && ec >= ymin) h += `<circle cx="${X(state.R)}" cy="${Y(ec)}" r="5" fill="var(--neg)" stroke="var(--panel)" stroke-width="2"/>`;
  }
  h += `<rect id="curveHit" x="${L}" y="${T}" width="${W - L - Rr}" height="${H - T - B}" fill="transparent" style="cursor:crosshair"/>`;
  svg.innerHTML = h;
  $('curveHit').addEventListener('click', e => {
    const r = svg.getBoundingClientRect(), x = (e.clientX - r.left) / r.width * W;
    state.R = +Math.max(0.5, Math.min(3.5, xmin + (x - L) / (W - L - Rr) * (xmax - xmin))).toFixed(3); optNote = ''; update();
  });
  const De = curve.Eatoms != null ? -emin : null;
  const farE = far;
  const edge = k >= ys.length - 2;
  let txt = edge ? 'No minimum: the energy keeps falling as the atoms move apart, so this method predicts no bond.' : `Minimum near <b>${Rmin.toFixed(3)} Å</b>.`;
  if (De != null && !edge) txt += De > 0 ? ` Well depth ${Math.round(De)} kJ/mol below the separated atoms.` : ` The minimum lies ${Math.round(-De)} kJ/mol <b>above</b> the separated atoms, so this method predicts no net bond.`;
  if (edge) { /* nothing more to say */ }
  else if (curve.Eatoms != null && farE > 40 && res && res.mult === 1) txt += ` At long range the curve stays ${Math.round(farE)} kJ/mol above the atoms: a single closed-shell determinant forces both electrons of each bond to stay paired, so restricted Hartree–Fock cannot break the bond correctly.`;
  else if (curve.Eatoms != null && farE > 40) txt += ` At long range the curve stays ${Math.round(farE)} kJ/mol above the atoms: one restricted determinant cannot describe the separated atoms.`;
  txt += ' A sharp corner marks a bond length where a different electron configuration becomes the lowest-energy solution. Click the curve to jump to that bond length.';
  $('curveStatus').innerHTML = txt;
}
function niceStep2(span) { const p = Math.pow(10, Math.floor(Math.log10(span))); for (const m of [1, 2, 2.5, 5, 10]) if (m * p >= span) return m * p; return 10 * p; }

// ---------- guided activities ----------
const ACTIVITIES = [
  { id: 'period2', title: 'Bond order across period 2', engine: 'eht', blurb: 'From H<sub>2</sub> to Ne<sub>2</sub>: count bonding and antibonding electrons and predict which molecules exist.', steps: [
    { title: 'Hydrogen, H<sub>2</sub>', set: { A: 'H', B: 'H', view: 'slice', sel: 'σ1s' },
      look: 'Two hydrogen 1s orbitals combine into a bonding σ<sub>1s</sub> and an antibonding σ*<sub>1s</sub>. Click each level in the diagram and compare the pictures.',
      q: 'Which orbital has a node between the nuclei? How many electrons are in bonding and in antibonding orbitals, and what is the bond order?',
      a: 'σ*<sub>1s</sub> has the node: a plane halfway between the nuclei. Both electrons sit in σ<sub>1s</sub>, so the bond order is (2 − 0) ÷ 2 = 1.' },
    { title: 'Helium, He<sub>2</sub>', set: { A: 'He', B: 'He', view: 'compare', cmp: { A: 'H', B: 'H' } },
      look: 'He<sub>2</sub> is on the left and H<sub>2</sub> on the right, drawn on the same energy scale.',
      q: 'Predict before you look: what is the bond order of He<sub>2</sub>? Compare how far σ<sub>1s</sub> drops below the He 1s level with how far σ*<sub>1s</sub> rises above it. Should He<sub>2</sub> exist?',
      a: 'Both orbitals are full, so the bond order is (2 − 2) ÷ 2 = 0. The antibonding level rises more than the bonding level falls, so filling both costs energy overall. He<sub>2</sub> is not a stable molecule.' },
    { title: 'The He<sub>2</sub><sup>+</sup> ion', set: { A: 'He', B: 'He', q: 1, R: 1.081, view: 'slice', sel: 'σ*1s' },
      look: 'One electron has been removed from He<sub>2</sub>.',
      q: 'Which orbital lost the electron, and what is the new bond order?',
      a: 'The electron comes out of σ*<sub>1s</sub>, the highest occupied orbital. The bond order becomes (2 − 1) ÷ 2 = ½. He<sub>2</sub><sup>+</sup> is a real ion that has been observed in helium discharges.' },
    { title: 'Lithium and beryllium', set: { A: 'Li', B: 'Li', view: 'compare', cmp: { A: 'Be', B: 'Be' } },
      look: 'Li<sub>2</sub> (left) and Be<sub>2</sub> (right). Core 1s electrons are not shown; they cancel out.',
      q: 'Find the bond orders of Li<sub>2</sub> and Be<sub>2</sub>. Which one would you expect to exist as a stable molecule?',
      a: 'Li<sub>2</sub> has (σ<sub>2s</sub>)<sup>2</sup>, bond order 1, and exists in lithium vapor. Be<sub>2</sub> has (σ<sub>2s</sub>)<sup>2</sup>(σ*<sub>2s</sub>)<sup>2</sup>, bond order 0, so it is not expected to be stable.' },
    { title: 'Boron, B<sub>2</sub>', set: { A: 'B', B: 'B', view: 'slice', sel: 'π2p' },
      look: 'B<sub>2</sub> has six valence electrons. Look at where the last two go.',
      q: 'Which orbitals hold the last two electrons, and how are their spins arranged? What does that predict about magnetism?',
      a: 'They go into the two π<sub>2p</sub> orbitals, one in each with parallel spins (Hund’s rule). B<sub>2</sub> has bond order 1 and is paramagnetic. Experiments confirm B<sub>2</sub> is paramagnetic, which is evidence that π<sub>2p</sub> lies below σ<sub>2p</sub> for the lighter elements.' },
    { title: 'Carbon and nitrogen', set: { A: 'C', B: 'C', view: 'compare', cmp: { A: 'N', B: 'N' } },
      look: 'C<sub>2</sub> (left) and N<sub>2</sub> (right).',
      q: 'Give the bond order, HOMO and magnetism of each.',
      a: 'C<sub>2</sub>: bond order 2, HOMO π<sub>2p</sub>, diamagnetic. N<sub>2</sub>: bond order 3, HOMO σ<sub>2p</sub>, diamagnetic. The triple bond is why N<sub>2</sub> is so unreactive.' },
    { title: 'Oxygen and fluorine', set: { A: 'O', B: 'O', view: 'compare', cmp: { A: 'F', B: 'F' } },
      look: 'O<sub>2</sub> (left) and F<sub>2</sub> (right).',
      q: 'O<sub>2</sub> has two more electrons than N<sub>2</sub>. Where do they go? What are the bond orders and magnetism of O<sub>2</sub> and F<sub>2</sub>? What would you predict for Ne<sub>2</sub>?',
      a: 'They go into π*<sub>2p</sub>, one in each orbital with parallel spins. O<sub>2</sub> has bond order 2 and is paramagnetic, which is why liquid oxygen sticks to a magnet. F<sub>2</sub> fills π*<sub>2p</sub>: bond order 1, diamagnetic. Ne<sub>2</sub> would also fill σ*<sub>2p</sub>: bond order 0, no bond.' },
  ] },
  { id: 'oxygen', title: 'Oxygen and its ions', engine: 'eht', blurb: 'Add and remove electrons from O<sub>2</sub> and connect bond order to bond length.', steps: [
    { title: 'The π* orbitals of O<sub>2</sub>', set: { A: 'O', B: 'O', view: 'slice', sel: 'π*2p' },
      look: 'The picture shows one of the two π*<sub>2p</sub> orbitals, sliced through the bond axis.',
      q: 'How many nodal planes does π*<sub>2p</sub> have, and where are they? How many electrons occupy the π* pair?',
      a: 'Two. One contains the bond axis (every π orbital has this one). The other is perpendicular to the axis, halfway between the nuclei; this one makes the orbital antibonding. The π* pair holds two electrons, one in each orbital.' },
    { title: 'O<sub>2</sub><sup>+</sup>', set: { A: 'O', B: 'O', q: 1, R: 1.116, view: 'slice', sel: 'π*2p' },
      look: 'One electron has been removed.',
      q: 'Which orbital lost the electron? Does the bond get stronger or weaker?',
      a: 'An electron leaves π*<sub>2p</sub>, so the bond order rises to 2½. Removing an antibonding electron strengthens the bond: O<sub>2</sub><sup>+</sup> is 1.12 Å long, compared with 1.21 Å for O<sub>2</sub>.' },
    { title: 'Superoxide, O<sub>2</sub><sup>−</sup>', set: { A: 'O', B: 'O', q: -1, R: 1.35, view: 'slice', sel: 'π*2p' },
      look: 'One electron has been added.',
      q: 'Where does the extra electron go? Give the bond order and magnetism.',
      a: 'It goes into π*<sub>2p</sub>. Bond order 1½, one unpaired electron, paramagnetic. The bond lengthens to about 1.35 Å.' },
    { title: 'Peroxide, O<sub>2</sub><sup>2−</sup>', set: { A: 'O', B: 'O', q: -2, R: 1.49, view: 'slice', sel: 'π*2p' },
      look: 'Two electrons have been added.',
      q: 'Give the bond order and magnetism.',
      a: 'Both π*<sub>2p</sub> orbitals are full. Bond order 1, no unpaired electrons, diamagnetic. The bond is about 1.49 Å, close to a typical O–O single bond.' },
    { title: 'Putting it together', set: { A: 'O', B: 'O', view: 'compare', cmp: { A: 'O', B: 'O', q: -2, R: 1.49 } },
      look: 'O<sub>2</sub> and peroxide side by side. Use the Compare chips to look at the other ions too.',
      q: 'Make a table of bond order and bond length for O<sub>2</sub><sup>+</sup>, O<sub>2</sub>, O<sub>2</sub><sup>−</sup> and O<sub>2</sub><sup>2−</sup>. What is the trend?',
      a: 'Bond order 2½ → 1.12 Å; 2 → 1.21 Å; 1½ → 1.35 Å; 1 → 1.49 Å. A higher bond order means a shorter, stronger bond.' },
  ] },
  { id: 'polar', title: 'Polar bonds: CO and HF', engine: 'hf', blurb: 'See how electronegativity skews orbitals toward one atom, and why CO binds metals through carbon.', steps: [
    { title: 'The bonding π orbital of CO', set: { A: 'C', B: 'O', view: 'slice', sel: 'π2p' },
      look: 'Carbon is on the left, oxygen on the right. The bars under the picture show how much each atomic orbital contributes.',
      q: 'Which atom contributes more to the bonding π<sub>2p</sub> orbital? Why?',
      a: 'Oxygen, about 70%. Oxygen is more electronegative, so its 2p orbitals start lower in energy (compare the side columns). A bonding orbital takes most of its character from the lower-energy atomic orbital.' },
    { title: 'The HOMO of CO', set: { A: 'C', B: 'O', view: 'slice', sel: 'σ2p' },
      look: 'This is the highest occupied orbital.',
      q: 'Which end of the molecule carries the large outer lobe of the HOMO?',
      a: 'The carbon end, pointing away from oxygen; about 80% of this orbital is on carbon. This lone-pair-like orbital is why CO bonds to metals through carbon, as in carbon monoxide poisoning of hemoglobin and in metal carbonyls.' },
    { title: 'The LUMO of CO', set: { A: 'C', B: 'O', view: 'slice', sel: 'π*2p' },
      look: 'This is the lowest empty orbital.',
      q: 'Where is the LUMO concentrated?',
      a: 'Mostly on carbon, about 70%. An antibonding orbital takes most of its character from the higher-energy atomic orbitals. A metal can push electron density back into this orbital (back-bonding), which weakens the C–O bond.' },
    { title: 'Hydrogen fluoride', set: { A: 'H', B: 'F', view: 'slice', sel: 'σ' },
      look: 'Hydrogen is on the left, fluorine on the right. Click each level.',
      q: 'Which molecular orbitals are bonding, and which stay on fluorine as nonbonding pairs? What is the bond order, and how does this match the Lewis structure?',
      a: 'Only σ, built from H 1s and F 2p<sub>z</sub>, is bonding. σ<sub>nb</sub> (mostly F 2s) and the two π<sub>nb</sub> orbitals stay on fluorine. That gives bond order 1 and three nonbonding pairs, matching the Lewis structure H–F with three lone pairs on fluorine.' },
    { title: 'Partial charges', set: { A: 'H', B: 'F', view: 'compare', cmp: { A: 'C', B: 'O' } },
      look: 'HF is in the main diagram and CO in the comparison. Compare their partial charges in the table.',
      q: 'Which end of HF is δ−? How do the partial charges in CO compare?',
      a: 'Fluorine is δ− and hydrogen δ+, about ±0.2 in this calculation. In CO the charges are also small. CO’s measured dipole moment is tiny (0.11 D), much smaller than the electronegativity difference alone would suggest.' },
  ] },
  { id: 'energy', title: 'Bond length and energy', engine: 'hf', blurb: 'Use Hartree–Fock energy curves to find bond lengths, and see what the method gets wrong.', steps: [
    { title: 'The N<sub>2</sub> energy curve', set: { A: 'N', B: 'N', view: 'curve' },
      look: 'The curve shows the energy of N<sub>2</sub> relative to two separate N atoms as the bond length changes.',
      q: 'Where is the minimum? Compare it with the experimental bond length, 1.098 Å.',
      a: 'Near 1.13 Å, about 3% longer than experiment. Hartree–Fock with a minimal basis typically gives bond lengths within a few percent.' },
    { title: 'How deep is the well?', set: { A: 'N', B: 'N', view: 'curve' },
      look: 'Click Optimize bond length, then read the well depth under the curve.',
      q: 'The experimental bond energy of N<sub>2</sub> is 945 kJ/mol. How close is the calculation?',
      a: 'The well is only about 160 kJ/mol deep, far too shallow. Hartree–Fock leaves out electron correlation (electrons avoiding each other moment to moment), and a single restricted wavefunction cannot describe two separated N atoms. Both errors make bonds look much weaker than they are.' },
    { title: 'He<sub>2</sub>', set: { A: 'He', B: 'He', view: 'curve' },
      look: 'The energy curve for two helium atoms.',
      q: 'What does the curve look like? Connect it to the bond order.',
      a: 'There is no minimum: the energy rises steadily as the atoms approach. With bond order 0, nothing holds the atoms together.' },
    { title: 'F<sub>2</sub>', set: { A: 'F', B: 'F', view: 'curve' },
      look: 'The energy curve for F<sub>2</sub>.',
      q: 'Read the well depth (experiment: about 155 kJ/mol). What do you notice?',
      a: 'The well is only about 20 kJ/mol deep. F<sub>2</sub> is the classic case where Hartree–Fock barely binds a molecule: many electron pairs sit close together, so the missing correlation matters a lot.' },
    { title: 'Singlet or triplet O<sub>2</sub>?', set: { A: 'O', B: 'O', mult: 1, view: 'slice' },
      look: 'This step forces O<sub>2</sub> into a singlet, with every electron paired. Note the total energy, then change Spin state to Triplet.',
      q: 'Which spin state is lower in energy, and by how much? (1 E<sub>h</sub> = 2625.5 kJ/mol.)',
      a: 'The triplet is lower, by about 210 kJ/mol here. Oxygen’s ground state is the triplet, which is why O<sub>2</sub> is paramagnetic. The real lowest singlet state (“singlet oxygen”) lies only 94 kJ/mol above the triplet; a single all-paired wavefunction describes it poorly.' },
  ] },
];

function buildActivities() {
  $('actCards').innerHTML = ACTIVITIES.map(act => `
    <button type="button" class="act-card" data-id="${act.id}">
      <span class="act-card-top"><span class="act-method">${act.engine === 'hf' ? 'Hartree–Fock' : 'Extended Hückel'}</span><span class="act-count">${act.steps.length} steps</span></span>
      <span class="act-card-title">${act.title}</span>
      <span class="act-card-blurb">${act.blurb}</span>
    </button>`).join('');
  $('actCards').querySelectorAll('.act-card').forEach(b => b.addEventListener('click', () => startActivity(b.dataset.id, 0)));
  $('actBack').addEventListener('click', () => gotoStep(state.act.step - 1));
  $('actNext').addEventListener('click', () => gotoStep(state.act.step + 1));
  $('actReveal').addEventListener('click', () => { const h = $('actA').hidden; $('actA').hidden = !h; $('actReveal').textContent = h ? 'Hide answer' : 'Show answer'; });
  $('actRedo').addEventListener('click', () => gotoStep(state.act.step));
  $('actExit').addEventListener('click', () => { state.act = null; saveAct(); renderActivity(); });
  try { const s = JSON.parse(localStorage.getItem('mo-lab-activity') || 'null'); if (s && ACTIVITIES.find(a => a.id === s.id)) state.actResume = s; } catch (e) { /* storage unavailable */ }
}
function saveAct() { try { localStorage.setItem('mo-lab-activity', JSON.stringify(state.act)); } catch (e) { /* ignore */ } }
function startActivity(id, step) { state.act = { id, step: 0 }; gotoStep(step); }
function gotoStep(i) {
  const act = ACTIVITIES.find(a => a.id === state.act.id);
  i = Math.max(0, Math.min(act.steps.length - 1, i));
  state.act.step = i; saveAct();
  const s = act.steps[i].set;
  state.engine = s.engine || act.engine;
  Object.assign(state, { A: s.A, B: s.B, q: s.q || 0, R: s.R ?? defaultBond(s.A, s.B), mult: s.mult || 0, sel: s.sel ? { name: s.sel } : null });
  if (s.cmp) state.cmp = { A: s.cmp.A, B: s.cmp.B, q: s.cmp.q || 0, R: s.cmp.R ?? defaultBond(s.cmp.A, s.cmp.B) };
  optNote = '';
  state.view = s.view || 'slice';
  update(); setView(state.view);
  renderActivity();
  $('activity').scrollIntoView({ block: 'nearest', behavior: reduceMotion ? 'auto' : 'smooth' });
}
function renderActivity() {
  const on = !!state.act;
  $('actCards').hidden = on; $('actBody').hidden = !on; $('actExit').hidden = !on;
  if (!on) {
    $('actTitle').textContent = 'Work through a topic step by step';
    const r = state.actResume, act = r && ACTIVITIES.find(a => a.id === r.id);
    $('actResume').hidden = !act;
    if (act) { $('actResume').innerHTML = `Resume “${act.title}” at step ${r.step + 1}`; $('actResume').onclick = () => startActivity(act.id, r.step); }
    return;
  }
  $('actResume').hidden = true;
  const act = ACTIVITIES.find(a => a.id === state.act.id), i = state.act.step, st = act.steps[i];
  $('actTitle').innerHTML = act.title;
  $('actStep').innerHTML = `Step ${i + 1} of ${act.steps.length}`;
  $('actDots').innerHTML = act.steps.map((_, k) => `<button type="button" class="dot${k === i ? ' on' : ''}${k < i ? ' done' : ''}" aria-label="Go to step ${k + 1}" data-k="${k}"></button>`).join('');
  $('actDots').querySelectorAll('.dot').forEach(d => d.addEventListener('click', () => gotoStep(+d.dataset.k)));
  $('actStepTitle').innerHTML = st.title;
  $('actLook').innerHTML = st.look;
  $('actQ').innerHTML = st.q;
  $('actA').innerHTML = st.a; $('actA').hidden = true; $('actReveal').textContent = 'Show answer';
  $('actBack').disabled = i === 0;
  $('actNext').disabled = i === act.steps.length - 1;
}

// theme + resize
const rerenderTheme = () => { readColors(); renderViewer(); };
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', rerenderTheme);
new MutationObserver(rerenderTheme).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
let rsT; new ResizeObserver(() => { clearTimeout(rsT); rsT = setTimeout(renderViewer, 60); }).observe(canvas);

// boot
buildPicker($('pickA'), 'A'); buildPicker($('pickB'), 'B'); buildCharge(); buildPresets(); buildParamTable(); buildCmpChips(); buildActivities(); renderActivity();
readColors();
worker = makeWorker();
update();
if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => renderViewer());
})();
