// Regression tests for both engines.  Run from the repository root:
//   python3 src/build.py && node tests/test_engines.js
// Reference total energies (hartree) were confirmed with the independent PyQuante2 program
// (STO-3G, same geometries); bond lengths are in Å.
const path = require('path');
const EHT = require(path.join(__dirname, '..', 'src', 'engine.js'));
const HF = require(path.join(__dirname, '..', 'src', 'hf.built.js'));
let fail = 0;
const check = (name, ok, detail) => { console.log(`${ok ? 'PASS' : 'FAIL'}  ${name}${detail ? '  ' + detail : ''}`); if (!ok) fail++; };

// --- Hartree–Fock total energies
const REF = [
  ['H', 'H', 0.741, 0, 1, -1.11670614], ['N', 'N', 1.098, 0, 1, -107.49597516], ['C', 'O', 1.128, 0, 1, -111.22455880],
  ['H', 'F', 0.917, 0, 1, -98.57077999], ['F', 'F', 1.412, 0, 1, -195.96743438], ['Li', 'Li', 2.673, 0, 1, -14.63872518],
  ['C', 'N', 1.16, -1, 1, -90.93765562], ['Be', 'Be', 2.45, 0, 1, -28.69877700], ['O', 'O', 1.208, 0, 3, -147.63219104],
  ['N', 'O', 1.151, 0, 2, -127.52611131], ['B', 'B', 1.59, 0, 3, -48.34749566], ['He', 'H', 0.774, 1, 1, -2.84177923],
  ['B', 'N', 1.281, 0, 1, -77.84305057], ['O', 'O', 1.208, 0, 1, -147.55112741], ['O', 'O', 1.35, -1, 2, -147.38088642],
  ['H', 'Cl', 1.275, 0, 1, -455.13483521], ['Cl', 'Cl', 1.988, 0, 1, -909.11162644],
];
for (const [A, B, R, q, m, E] of REF) {
  const r = HF.diatomic(A, B, R, q, m);
  check(`HF ${A}${B}${q ? (q > 0 ? '+' : '') + q : ''} (2S+1=${m})`, r.converged && Math.abs(r.E - E) < 1e-6, `E=${r.E.toFixed(8)} ref=${E.toFixed(8)}`);
}
// --- rotation invariance (tests every p-type integral)
{
  const Rb = 1.128 * HF.BOHR_PER_ANG, u = [1, 1, 1].map(v => v / Math.sqrt(3)), o = [0.3, -0.2, 0.1];
  const ez = HF.scf(HF.setup([{ sym: 'C', pos: [0, 0, 0] }, { sym: 'O', pos: [0, 0, Rb] }]), 14, 1).E;
  const er = HF.scf(HF.setup([{ sym: 'C', pos: o }, { sym: 'O', pos: o.map((v, k) => v + u[k] * Rb) }]), 14, 1).E;
  check('HF rotation invariance (CO)', Math.abs(ez - er) < 1e-9, `diff=${(ez - er).toExponential(2)}`);
}
// --- optimized STO-3G bond lengths
for (const [A, B, R0, m, Rexp] of [['H', 'H', 0.74, 1, 0.712], ['N', 'N', 1.1, 1, 1.134], ['C', 'O', 1.13, 1, 1.146], ['H', 'F', 0.92, 1, 0.956]]) {
  const o = HF.optimize(A, B, R0, 0, m);
  check(`HF optimized R ${A}${B}`, Math.abs(o.R - Rexp) < 0.002, `R=${o.R.toFixed(4)} expected≈${Rexp}`);
}
// --- extended Hückel: bond orders and magnetism
for (const [A, B, q, bo, unp] of [['H', 'H', 0, 1, 0], ['He', 'He', 0, 0, 0], ['Li', 'Li', 0, 1, 0], ['Be', 'Be', 0, 0, 0], ['B', 'B', 0, 1, 2],
  ['C', 'C', 0, 2, 0], ['N', 'N', 0, 3, 0], ['O', 'O', 0, 2, 2], ['F', 'F', 0, 1, 0], ['Ne', 'Ne', 0, 0, 0], ['O', 'O', 1, 2.5, 1],
  ['O', 'O', -1, 1.5, 1], ['O', 'O', -2, 1, 0], ['C', 'O', 0, 3, 0], ['N', 'O', 0, 2.5, 1], ['H', 'F', 0, 1, 0]]) {
  const r = EHT.calculate(A, B, q, EHT.defaultBond(A, B));
  check(`EHT ${A}${B}${q || ''} bond order / unpaired`, r.bondOrder === bo && r.unpaired === unp, `${r.bondOrder}, ${r.unpaired}`);
}
console.log(fail ? `\n${fail} test(s) failed` : '\nAll tests passed');
process.exit(fail ? 1 : 0);
