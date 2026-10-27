// ===== Extended Hückel engine for diatomics =====
const BOHR = 1.8897259886; // bohr per Å

// Hoffmann / Alvarez EHT parameters: Hii (eV) = -VSIP, zeta = Slater exponent
const ATOMS = {
  H:  { Z: 1,  name: 'Hydrogen',  val: 1, cov: 0.31, shells: [{ n: 1, l: 0, H: -13.6, z: 1.300 }] },
  He: { Z: 2,  name: 'Helium',    val: 2, cov: 0.28, shells: [{ n: 1, l: 0, H: -23.4, z: 1.688 }] },
  Li: { Z: 3,  name: 'Lithium',   val: 1, cov: 1.28, shells: [{ n: 2, l: 0, H: -5.4,  z: 0.650 }, { n: 2, l: 1, H: -3.5,  z: 0.650 }] },
  Be: { Z: 4,  name: 'Beryllium', val: 2, cov: 0.96, shells: [{ n: 2, l: 0, H: -10.0, z: 0.975 }, { n: 2, l: 1, H: -6.0,  z: 0.975 }] },
  B:  { Z: 5,  name: 'Boron',     val: 3, cov: 0.84, shells: [{ n: 2, l: 0, H: -15.2, z: 1.300 }, { n: 2, l: 1, H: -8.5,  z: 1.300 }] },
  C:  { Z: 6,  name: 'Carbon',    val: 4, cov: 0.76, shells: [{ n: 2, l: 0, H: -21.4, z: 1.625 }, { n: 2, l: 1, H: -11.4, z: 1.625 }] },
  N:  { Z: 7,  name: 'Nitrogen',  val: 5, cov: 0.71, shells: [{ n: 2, l: 0, H: -26.0, z: 1.950 }, { n: 2, l: 1, H: -13.4, z: 1.950 }] },
  O:  { Z: 8,  name: 'Oxygen',    val: 6, cov: 0.66, shells: [{ n: 2, l: 0, H: -32.3, z: 2.275 }, { n: 2, l: 1, H: -14.8, z: 2.275 }] },
  F:  { Z: 9,  name: 'Fluorine',  val: 7, cov: 0.57, shells: [{ n: 2, l: 0, H: -40.0, z: 2.425 }, { n: 2, l: 1, H: -18.1, z: 2.425 }] },
  Ne: { Z: 10, name: 'Neon',      val: 8, cov: 0.58, shells: [{ n: 2, l: 0, H: -43.2, z: 2.879 }, { n: 2, l: 1, H: -20.0, z: 2.879 }] },
  Cl: { Z: 17, name: 'Chlorine',  val: 7, cov: 1.02, shells: [{ n: 3, l: 0, H: -26.3, z: 2.183 }, { n: 3, l: 1, H: -14.2, z: 1.733 }] },
};

// Experimental bond lengths (Å); keys sorted by symbol order A|B as entered
const BONDS = {
  'H|H': 0.741, 'He|He': 1.50, 'He|H': 0.774, 'H|He': 0.774, 'Li|Li': 2.673, 'Be|Be': 2.45, 'B|B': 1.590,
  'C|C': 1.243, 'N|N': 1.098, 'O|O': 1.208, 'F|F': 1.412, 'Ne|Ne': 3.10, 'Cl|Cl': 1.988,
  'C|O': 1.128, 'N|O': 1.151, 'C|N': 1.172, 'H|F': 0.917, 'O|H': 0.970, 'H|O': 0.970, 'H|Cl': 1.275,
  'B|N': 1.281, 'Li|H': 1.595, 'H|Li': 1.595, 'Be|H': 1.343, 'B|H': 1.232, 'C|H': 1.120, 'H|C': 1.120,
  'N|H': 1.036, 'H|N': 1.036, 'Cl|F': 1.628, 'F|Cl': 1.628, 'B|F': 1.263, 'Li|F': 1.564, 'B|O': 1.205,
  'C|F': 1.272, 'N|F': 1.317, 'O|F': 1.354, 'C|Cl': 1.645, 'Li|O': 1.688,
};
function defaultBond(a, b) {
  const k1 = a + '|' + b, k2 = b + '|' + a;
  if (BONDS[k1]) return BONDS[k1];
  if (BONDS[k2]) return BONDS[k2];
  return +(ATOMS[a].cov + ATOMS[b].cov).toFixed(3);
}

// ---- Gauss–Legendre nodes on [-1,1]
function gauleg(n) {
  const x = new Float64Array(n), w = new Float64Array(n);
  for (let i = 0; i < Math.ceil(n / 2); i++) {
    let z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5)), pp = 0;
    for (let it = 0; it < 100; it++) {
      let p1 = 1, p2 = 0;
      for (let j = 1; j <= n; j++) { const p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
      pp = n * (z * p1 - p2) / (z * z - 1);
      const z1 = z; z = z1 - p1 / pp;
      if (Math.abs(z - z1) < 1e-15) break;
    }
    x[i] = -z; x[n - 1 - i] = z;
    w[i] = w[n - 1 - i] = 2 / ((1 - z * z) * pp * pp);
  }
  return { x, w };
}
const GL_XI = gauleg(72), GL_ETA = gauleg(56);

const FACT = n => { let f = 1; for (let i = 2; i <= n; i++) f *= i; return f; };
const Y_S = 1 / Math.sqrt(4 * Math.PI), Y_P = Math.sqrt(3 / (4 * Math.PI));

// Build basis functions. Atom A at z=0, atom B at z=Rb (bohr). Axis: z, A→B.
function buildBasis(symA, symB, Rb) {
  const basis = [];
  [[symA, 0, 0], [symB, 1, Rb]].forEach(([sym, a, z0]) => {
    for (const sh of ATOMS[sym].shells) {
      const N = Math.pow(2 * sh.z, sh.n + 0.5) / Math.sqrt(FACT(2 * sh.n));
      const shell = sh.n + (sh.l === 0 ? 's' : 'p');
      const types = sh.l === 0 ? ['s'] : ['px', 'py', 'pz'];
      for (const t of types) basis.push({ atom: a, sym, z0, n: sh.n, l: sh.l, zeta: sh.z, N, H: sh.H, type: t, shell, label: shell + (t === 's' ? '' : t.slice(1)) });
    }
  });
  return basis;
}

// Value of a basis function at Cartesian point (x,y,z)
function bf(b, x, y, z) {
  const dz = z - b.z0, r = Math.sqrt(x * x + y * y + dz * dz);
  const rad = b.N * Math.pow(r, b.n - 1) * Math.exp(-b.zeta * r);
  if (b.type === 's') return rad * Y_S;
  if (r < 1e-12) return 0;
  const c = b.type === 'px' ? x : b.type === 'py' ? y : dz;
  return rad * Y_P * c / r;
}

// Two-centre overlap by quadrature in prolate spheroidal coordinates
function overlap2c(ba, bb, Rb) {
  const piType = ba.type === 'px' || ba.type === 'py';
  if (piType !== (bb.type === 'px' || bb.type === 'py')) return 0;
  if (piType && ba.type !== bb.type) return 0;
  const h = Rb / 2, c = Math.min(4 / ((ba.zeta + bb.zeta) * h), 40);
  const ev = piType ? 'px' : null;
  const A = piType ? { ...ba, type: 'px' } : ba, B = piType ? { ...bb, type: 'px' } : bb;
  let sum = 0;
  for (let i = 0; i < GL_XI.x.length; i++) {
    const s = 0.5 * (GL_XI.x[i] + 1), ws = 0.5 * GL_XI.w[i];
    const xi = 1 + c * s / (1 - s), dxi = c / ((1 - s) * (1 - s));
    for (let j = 0; j < GL_ETA.x.length; j++) {
      const eta = GL_ETA.x[j], we = GL_ETA.w[j];
      const rho = h * Math.sqrt(Math.max(0, (xi * xi - 1) * (1 - eta * eta)));
      const z = h * (1 + xi * eta);
      sum += ws * we * dxi * (xi * xi - eta * eta) * bf(A, rho, 0, z) * bf(B, rho, 0, z);
    }
  }
  return sum * h * h * h * (piType ? Math.PI : 2 * Math.PI);
}

// ---- small symmetric eigen solver (cyclic Jacobi)
function jacobi(Ain) {
  const n = Ain.length, A = Ain.map(r => r.slice()), V = A.map((r, i) => r.map((_, j) => (i === j ? 1 : 0)));
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += A[p][q] * A[p][q];
    if (off < 1e-24) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(A[p][q]) < 1e-30) continue;
      const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
      const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const cs = 1 / Math.sqrt(t * t + 1), sn = t * cs;
      for (let k = 0; k < n; k++) { const akp = A[k][p], akq = A[k][q]; A[k][p] = cs * akp - sn * akq; A[k][q] = sn * akp + cs * akq; }
      for (let k = 0; k < n; k++) { const apk = A[p][k], aqk = A[q][k]; A[p][k] = cs * apk - sn * aqk; A[q][k] = sn * apk + cs * aqk; }
      for (let k = 0; k < n; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = cs * vkp - sn * vkq; V[k][q] = sn * vkp + cs * vkq; }
    }
  }
  return { vals: A.map((r, i) => r[i]), vecs: V };
}
const matmul = (A, B) => A.map(r => B[0].map((_, j) => r.reduce((s, v, k) => s + v * B[k][j], 0)));
const transpose = A => A[0].map((_, j) => A.map(r => r[j]));

// Solve HC = SCE for a symmetry block (Löwdin orthogonalization)
function solveGeneral(H, S) {
  const { vals, vecs } = jacobi(S);
  const X = matmul(matmul(vecs, vals.map((v, i) => vals.map((_, j) => (i === j ? 1 / Math.sqrt(v) : 0)))), transpose(vecs));
  const Hp = matmul(matmul(X, H), X);
  const e = jacobi(Hp);
  const C = matmul(X, e.vecs);
  return e.vals.map((E, k) => ({ E, c: C.map(r => r[k]) })).sort((a, b) => a.E - b.E);
}

// ---- main calculation
function calculate(symA, symB, charge, R_ang, opts = {}) {
  const K = opts.K ?? 1.75, spMix = opts.spMix ?? 1;
  const Rb = R_ang * BOHR;
  const basis = buildBasis(symA, symB, Rb), nb = basis.length;
  const S = basis.map((_, i) => basis.map((_, j) => (i === j ? 1 : 0)));
  const cache = {};
  for (let i = 0; i < nb; i++) for (let j = 0; j < nb; j++) {
    const a = basis[i], b = basis[j];
    if (a.atom !== 0 || b.atom !== 1) continue;
    const isSP = (a.type === 's' && b.type === 'pz') || (a.type === 'pz' && b.type === 's');
    const key = [a.label, a.zeta, b.label, b.zeta, a.type, b.type].join();
    let v = cache[key] ?? (cache[key] = overlap2c(a, b, Rb));
    if (isSP && hasP(symA) && hasP(symB)) v *= spMix;
    S[i][j] = S[j][i] = v;
  }
  const H = basis.map((a, i) => basis.map((b, j) => (i === j ? a.H : K * S[i][j] * (a.H + b.H) / 2)));

  // symmetry blocks: sigma (s, pz), pi-x, pi-y
  const blocks = { σ: [], πx: [], πy: [] };
  basis.forEach((b, i) => { (b.type === 'px' ? blocks.πx : b.type === 'py' ? blocks.πy : blocks.σ).push(i); });
  const mos = [];
  for (const [name, idx] of Object.entries(blocks)) {
    if (!idx.length) continue;
    const sub = M => idx.map(i => idx.map(j => M[i][j]));
    for (const r of solveGeneral(sub(H), sub(S))) {
      const c = new Array(nb).fill(0);
      idx.forEach((i, k) => (c[i] = r.c[k]));
      // phase convention: largest-magnitude coefficient on atom A positive, else first nonzero positive
      // phase convention: positive at the bond midpoint (bonding region shows the + colour);
      // if there is a node there, make atom A's largest coefficient positive
      const px = name === 'σ' ? 0 : 0.8, py = 0;
      let mid = 0; for (let i = 0; i < nb; i++) if (c[i]) mid += c[i] * bf(name === 'πy' ? { ...basis[i], type: basis[i].type === 'py' ? 'px' : basis[i].type } : basis[i], px, py, Rb / 2);
      let flip;
      if (Math.abs(mid) > 1e-4) flip = mid < 0;
      else {
        let ref = -1, big = 0;
        idx.forEach(i => { const w = Math.abs(c[i]) * (basis[i].atom === 0 ? 1.001 : 1); if (w > big) { big = w; ref = i; } });
        flip = ref >= 0 && c[ref] < 0;
      }
      if (flip) for (let i = 0; i < nb; i++) c[i] = -c[i];
      mos.push({ E: r.E, sym: name === 'σ' ? 'σ' : 'π', comp: name === 'πx' ? 'x' : name === 'πy' ? 'y' : '', c });
    }
  }
  mos.sort((a, b) => a.E - b.E || (a.comp < b.comp ? -1 : 1));

  // Mulliken analysis per MO
  for (const m of mos) {
    let wA = 0, wB = 0, q = 0;
    m.ao = basis.map((_, i) => { let g = 0; for (let j = 0; j < nb; j++) g += m.c[i] * m.c[j] * S[i][j]; return g; });
    basis.forEach((b, i) => { if (b.atom === 0) wA += m.ao[i]; else wB += m.ao[i]; });
    for (let i = 0; i < nb; i++) for (let j = 0; j < nb; j++) if (basis[i].atom === 0 && basis[j].atom === 1) q += 2 * m.c[i] * m.c[j] * S[i][j];
    m.wA = wA; m.wB = wB; m.q = q;
  }

  labelMOs(mos, basis, symA, symB, spMix);

  // Electron filling (Aufbau + Hund)
  const nElec = Math.max(0, Math.min(2 * mos.length, ATOMS[symA].val + ATOMS[symB].val - charge));
  let left = nElec;
  mos.forEach(m => (m.occ = 0));
  for (let i = 0; i < mos.length;) {
    let j = i; while (j < mos.length && Math.abs(mos[j].E - mos[i].E) < 1e-6) j++;
    const g = j - i, put = Math.min(left, 2 * g);
    for (let k = 0; k < g; k++) mos[i + k].occ = Math.min(2, Math.floor(put / g) + (k < put % g ? 1 : 0));
    // Hund: with put < g each gets at most 1; with g < put < 2g, pairs first then singles — fix ordering
    if (put <= g) for (let k = 0; k < g; k++) mos[i + k].occ = k < put ? 1 : 0;
    else for (let k = 0; k < g; k++) mos[i + k].occ = k < put - g ? 2 : 1;
    left -= put; i = j;
  }
  const unpaired = mos.reduce((s, m) => s + (m.occ === 1 ? 1 : 0), 0);
  let bonding = 0, anti = 0;
  mos.forEach(m => { if (m.kind === 'b') bonding += m.occ; else if (m.kind === 'a') anti += m.occ; });
  const bondOrder = (bonding - anti) / 2;
  let popA = 0, popB = 0;
  mos.forEach(m => { popA += m.occ * m.wA; popB += m.occ * m.wB; });
  const qA = ATOMS[symA].val - popA, qB = ATOMS[symB].val - popB;
  // effective charge distribution: assign ion charge implicitly via populations
  const occIdx = mos.map((m, i) => (m.occ > 0 ? i : -1)).filter(i => i >= 0);
  const homo = occIdx.length ? mos[occIdx[occIdx.length - 1]] : null;
  const lumo = mos.find(m => m.occ < 2 && (!homo || m.E > homo.E + 1e-6) && m.occ === 0) || null;
  return { symA, symB, charge, R: R_ang, Rb, basis, S, H, mos, nElec, unpaired, bondOrder, bonding, anti, qA, qB, homo, lumo, K, spMix };
}
const hasP = s => ATOMS[s].shells.some(sh => sh.l === 1);

function labelMOs(mos, basis, symA, symB, spMix) {
  const pA = hasP(symA), pB = hasP(symB);
  const nA = ATOMS[symA].shells[0].n, nB = ATOMS[symB].shells[0].n;
  const textbook = (pA && pB) || (!pA && !pB);
  const sig = mos.filter(m => m.sym === 'σ');
  const pix = mos.filter(m => m.sym === 'π' && m.comp === 'x');
  const piy = mos.filter(m => m.sym === 'π' && m.comp === 'y');
  const dominant = m => {
    let best = 0, bi = 0; m.ao.forEach((g, i) => { if (g > best) { best = g; bi = i; } });
    const b = basis[bi]; return b.sym + ' ' + b.shell;
  };
  mos.forEach(m => (m.where = dominant(m)));
  if (textbook) {
    const sh = nA === nB ? String(nA) : '';
    if (!pA) {
      const lab = [['σ', '', '1s', 'b'], ['σ', '*', '1s', 'a']];
      sig.forEach((m, i) => { [m.g, m.star, m.sub, m.kind] = lab[i]; });
    } else {
      const set = (m, star, s) => { m.g = 'σ'; m.star = star; m.sub = sh + s; m.kind = star ? 'a' : 'b'; };
      if (symA === symB) {
        // inversion parity: s_A→s_B, pz_A→−pz_B. g σ MOs: σs (lower), σp (upper); u: σ*s, σ*p
        // pair every function on A with its mirror image on B (same label); works with or without core functions
        const mirror = basis.map((b, i) => (b.atom === 0 ? [i, basis.findIndex(x => x.atom === 1 && x.label === b.label), b.type === 'pz' ? -1 : 1] : null)).filter(Boolean);
        const par = m => mirror.reduce((t, [i, j, s]) => t + s * m.c[i] * m.c[j], 0) >= 0 ? 'g' : 'u';
        const g = sig.filter(m => par(m) === 'g'), u = sig.filter(m => par(m) === 'u');
        if (g.length === 2 && u.length === 2) { set(g[0], '', 's'); set(g[1], '', 'p'); set(u[0], '*', 's'); set(u[1], '*', 'p'); }
        else sig.forEach((m, i) => set(m, i % 2 ? '*' : '', i < 2 ? 's' : 'p'));
        sig.forEach(m => (m.parity = par(m)));
      } else {
        // heteronuclear: the two most s-rich σ MOs form the "s" pair; within each pair lower = bonding
        const sChar = m => basis.reduce((t, b, i) => t + (b.type === 's' ? m.ao[i] : 0), 0);
        const rank = spMix === 0 ? sig.map((m, i) => [sChar(m), i]).sort((a, b) => b[0] - a[0]) : [[0, 0], [0, 1], [0, 2], [0, 3]];
        const sIdx = rank.slice(0, 2).map(r => r[1]).sort((a, b) => a - b), pIdx = rank.slice(2).map(r => r[1]).sort((a, b) => a - b);
        set(sig[sIdx[0]], '', 's'); set(sig[sIdx[1]], '*', 's'); set(sig[pIdx[0]], '', 'p'); set(sig[pIdx[1]], '*', 'p');
      }
      [pix, piy].forEach(arr => arr.forEach((m, i) => { [m.g, m.star, m.sub, m.kind] = i === 0 ? ['π', '', sh + 'p', 'b'] : ['π', '*', sh + 'p', 'a']; }));
    }
  } else {
    // One atom has only 1s (H/He), the other has s + p (textbook H–X convention)
    const X = pA ? symA : symB, xS = X + ' ' + ATOMS[X].shells[0].n + 's';
    const set = (m, k) => { m.g = m.sym; m.kind = k; m.star = k === 'a' ? '*' : ''; m.sub = k === 'n' ? 'nb' : ''; };
    [...pix, ...piy].forEach(m => set(m, 'n'));
    const kinds = sig[0].where === xS ? ['n', 'b', 'a'] : ['b', 'n', 'a'];
    sig.forEach((m, i) => set(m, kinds[i]));
  }
  mos.forEach(m => (m.name = m.g + m.star + m.sub));
}

// Atomic orbital occupancy for a neutral atom (for the side columns)
function atomLevels(sym) {
  let e = ATOMS[sym].val;
  return ATOMS[sym].shells.map(sh => {
    const deg = sh.l === 0 ? 1 : 3, put = Math.min(e, 2 * deg); e -= put;
    const occ = []; for (let k = 0; k < deg; k++) occ.push(put <= deg ? (k < put ? 1 : 0) : (k < put - deg ? 2 : 1));
    return { H: sh.H, deg, occ, label: sh.n + (sh.l === 0 ? 's' : 'p') };
  });
}

// Molecular orbital value at a point
function moValue(res, m, x, y, z) {
  let v = 0;
  for (let i = 0; i < res.basis.length; i++) if (Math.abs(m.c[i]) > 1e-9) v += m.c[i] * bf(res.basis[i], x, y, z);
  return v;
}

if (typeof module !== 'undefined') module.exports = { calculate, overlap2c, buildBasis, BOHR, ATOMS, defaultBond, atomLevels, moValue };
