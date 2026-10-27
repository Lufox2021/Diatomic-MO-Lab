// ===== Hartree–Fock (RHF / ROHF) with STO-3G for diatomics =====
// McMurchie–Davidson integrals over contracted Cartesian Gaussians (s and p only).
var HF = (function () {
  const HARTREE_EV = 27.211386245988;
  const BOHR_PER_ANG = 1.8897259886;
  const ZNUM = { H: 1, He: 2, Li: 3, Be: 4, B: 5, C: 6, N: 7, O: 8, F: 9, Ne: 10, Cl: 17 };
  /*STO3G*/
  const SQPI3 = Math.pow(Math.PI, 1.5);

  // ---------- basis ----------
  function shellsFor(sym, center, atomIndex) {
    const data = STO3G[ZNUM[sym]], out = [];
    let ns = 0, np = 1;
    for (const [L, prims] of data) {
      const l = L === 'S' ? 0 : 1, n = l === 0 ? ++ns : ++np;
      const exps = prims.map(p => p[0]);
      // fold primitive normalization (x^l component normalisation) into coefficients
      let coefs = prims.map(([a, c]) => c * Math.pow(2 * a / Math.PI, 0.75) * (l ? 2 * Math.sqrt(a) : 1));
      // renormalise contraction
      let ss = 0;
      for (let i = 0; i < exps.length; i++) for (let j = 0; j < exps.length; j++) {
        const p = exps[i] + exps[j];
        ss += coefs[i] * coefs[j] * Math.pow(Math.PI / p, 1.5) * (l ? 1 / (2 * p) : 1);
      }
      coefs = coefs.map(c => c / Math.sqrt(ss));
      out.push({ atom: atomIndex, sym, center, l, exps, coefs, shell: n + (l ? 'p' : 's') });
    }
    return out;
  }
  function buildBasis(atoms) {
    const shells = [], funcs = [];
    atoms.forEach((at, ia) => shellsFor(at.sym, at.pos, ia).forEach(sh => {
      sh.f0 = funcs.length;
      const comps = sh.l === 0 ? [[0, 0, 0, 's']] : [[1, 0, 0, 'px'], [0, 1, 0, 'py'], [0, 0, 1, 'pz']];
      sh.comps = comps;
      comps.forEach(([lx, ly, lz, type]) => funcs.push({ atom: ia, sym: at.sym, shell: sh.shell, type, lx, ly, lz, center: at.pos, exps: sh.exps, coefs: sh.coefs, label: sh.shell + (type === 's' ? '' : type.slice(1)) }));
      shells.push(sh);
    }));
    return { shells, funcs, n: funcs.length };
  }

  // ---------- Boys function ----------
  function boys(nmax, T) {
    const F = new Float64Array(nmax + 1);
    if (T < 1e-12) { for (let n = 0; n <= nmax; n++) F[n] = 1 / (2 * n + 1); return F; }
    if (T > 35) {
      F[0] = 0.5 * Math.sqrt(Math.PI / T);
      for (let n = 1; n <= nmax; n++) F[n] = F[n - 1] * (2 * n - 1) / (2 * T);
      return F;
    }
    // series for the highest order, then downward recursion (stable)
    const m = nmax, eT = Math.exp(-T);
    let term = 1 / (2 * m + 1), sum = term;
    for (let i = 1; i < 300; i++) { term *= 2 * T / (2 * m + 2 * i + 1); sum += term; if (term < 1e-17 * sum) break; }
    F[m] = eT * sum;
    for (let n = m - 1; n >= 0; n--) F[n] = (2 * T * F[n + 1] + eT) / (2 * n + 1);
    return F;
  }

  // ---------- Hermite expansion coefficients E^{ij}_t (1D) ----------
  function Etab(imax, jmax, Xab, a, b) {
    const p = a + b, mu = a * b / p, XPA = -b * Xab / p, XPB = a * Xab / p, ip2 = 1 / (2 * p);
    const T = imax + jmax + 1;
    const E = [];
    for (let i = 0; i <= imax; i++) { E.push([]); for (let j = 0; j <= jmax; j++) E[i].push(new Float64Array(T + 1)); }
    E[0][0][0] = Math.exp(-mu * Xab * Xab);
    const g = (arr, t) => (t < 0 || t > T ? 0 : arr[t]);
    for (let i = 0; i <= imax; i++) for (let j = 0; j <= jmax; j++) {
      if (i === 0 && j === 0) continue;
      const out = E[i][j];
      if (i > 0) { const s = E[i - 1][j]; for (let t = 0; t <= i + j; t++) out[t] = ip2 * g(s, t - 1) + XPA * g(s, t) + (t + 1) * g(s, t + 1); }
      else { const s = E[i][j - 1]; for (let t = 0; t <= i + j; t++) out[t] = ip2 * g(s, t - 1) + XPB * g(s, t) + (t + 1) * g(s, t + 1); }
    }
    return E;
  }

  // ---------- Hermite Coulomb integrals R_{tuv} ----------
  function Rtab(L, alpha, X, Y, Z) {
    const F = boys(L, alpha * (X * X + Y * Y + Z * Z));
    const D = L + 1;
    // R[n][t][u][v] flattened; only t+u+v <= L-n needed
    const R = new Float64Array(D * D * D * D);
    const idx = (n, t, u, v) => ((n * D + t) * D + u) * D + v;
    let f = 1;
    for (let n = 0; n <= L; n++) { R[idx(n, 0, 0, 0)] = f * F[n]; f *= -2 * alpha; }
    for (let n = L - 1; n >= 0; n--) {
      const lmax = L - n;
      for (let t = 0; t <= lmax; t++) for (let u = 0; u <= lmax - t; u++) for (let v = 0; v <= lmax - t - u; v++) {
        if (t + u + v === 0) continue;
        let val;
        if (t > 0) val = (t > 1 ? (t - 1) * R[idx(n + 1, t - 2, u, v)] : 0) + X * R[idx(n + 1, t - 1, u, v)];
        else if (u > 0) val = (u > 1 ? (u - 1) * R[idx(n + 1, t, u - 2, v)] : 0) + Y * R[idx(n + 1, t, u - 1, v)];
        else val = (v > 1 ? (v - 1) * R[idx(n + 1, t, u, v - 2)] : 0) + Z * R[idx(n + 1, t, u, v - 1)];
        R[idx(n, t, u, v)] = val;
      }
    }
    return (t, u, v) => R[idx(0, t, u, v)];
  }

  // ---------- one-electron integrals ----------
  function oneElectron(basis, atoms) {
    const n = basis.n, S = zeros(n), T = zeros(n), V = zeros(n);
    const sh = basis.shells;
    for (let I = 0; I < sh.length; I++) for (let J = 0; J <= I; J++) {
      const A = sh[I], B = sh[J];
      for (let pa = 0; pa < A.exps.length; pa++) for (let pb = 0; pb < B.exps.length; pb++) {
        const a = A.exps[pa], b = B.exps[pb], p = a + b, cc = A.coefs[pa] * B.coefs[pb];
        const Ex = Etab(1, 3, A.center[0] - B.center[0], a, b), Ey = Etab(1, 3, A.center[1] - B.center[1], a, b), Ez = Etab(1, 3, A.center[2] - B.center[2], a, b);
        const P = [0, 1, 2].map(k => (a * A.center[k] + b * B.center[k]) / p);
        const sp = Math.sqrt(Math.PI / p);
        const s1 = (E, i, j) => (j < 0 ? 0 : E[i][j][0] * sp);
        const t1 = (E, i, j) => -2 * b * b * s1(E, i, j + 2) + b * (2 * j + 1) * s1(E, i, j) - 0.5 * j * (j - 1) * s1(E, i, j - 2);
        const Rs = atoms.map(C => Rtab(2, p, P[0] - C.pos[0], P[1] - C.pos[1], P[2] - C.pos[2]));
        for (const ca of A.comps) for (const cb of B.comps) {
          const i = A.f0 + A.comps.indexOf(ca), j = B.f0 + B.comps.indexOf(cb);
          if (j > i) continue;
          const [ax, ay, az] = ca, [bx, by, bz] = cb;
          const sx = s1(Ex, ax, bx), sy = s1(Ey, ay, by), sz = s1(Ez, az, bz);
          S[i][j] += cc * sx * sy * sz;
          T[i][j] += cc * (t1(Ex, ax, bx) * sy * sz + sx * t1(Ey, ay, by) * sz + sx * sy * t1(Ez, az, bz));
          let v = 0;
          atoms.forEach((C, k) => {
            const R = Rs[k]; let s = 0;
            for (let t = 0; t <= ax + bx; t++) for (let u = 0; u <= ay + by; u++) for (let w = 0; w <= az + bz; w++)
              s += Ex[ax][bx][t] * Ey[ay][by][u] * Ez[az][bz][w] * R(t, u, w);
            v -= ZNUM[C.sym] * s;
          });
          V[i][j] += cc * 2 * Math.PI / p * v;
        }
      }
    }
    for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) { S[j][i] = S[i][j]; T[j][i] = T[i][j]; V[j][i] = V[i][j]; }
    return { S, T, V };
  }

  // ---------- two-electron integrals (full n^4 array) ----------
  function twoElectron(basis) {
    const n = basis.n, sh = basis.shells, ERI = new Float64Array(n * n * n * n);
    const id = (i, j, k, l) => ((i * n + j) * n + k) * n + l;
    // precompute shell-pair data
    const pairs = [];
    for (let I = 0; I < sh.length; I++) for (let J = 0; J <= I; J++) {
      const A = sh[I], B = sh[J], prims = [];
      for (let pa = 0; pa < A.exps.length; pa++) for (let pb = 0; pb < B.exps.length; pb++) {
        const a = A.exps[pa], b = B.exps[pb], p = a + b;
        prims.push({ p, P: [0, 1, 2].map(k => (a * A.center[k] + b * B.center[k]) / p), c: A.coefs[pa] * B.coefs[pb],
          E: [0, 1, 2].map(k => Etab(A.l, B.l, A.center[k] - B.center[k], a, b)) });
      }
      pairs.push({ I, J, A, B, prims });
    }
    const pref0 = 2 * Math.pow(Math.PI, 2.5);
    for (let P1 = 0; P1 < pairs.length; P1++) for (let P2 = 0; P2 <= P1; P2++) {
      const ab = pairs[P1], cd = pairs[P2];
      const L = ab.A.l + ab.B.l + cd.A.l + cd.B.l;
      const nA = ab.A.comps.length, nB = ab.B.comps.length, nC = cd.A.comps.length, nD = cd.B.comps.length;
      const acc = new Float64Array(nA * nB * nC * nD);
      for (const x of ab.prims) for (const y of cd.prims) {
        const p = x.p, q = y.p, alpha = p * q / (p + q);
        const R = Rtab(L, alpha, x.P[0] - y.P[0], x.P[1] - y.P[1], x.P[2] - y.P[2]);
        const pref = pref0 / (p * q * Math.sqrt(p + q)) * x.c * y.c;
        let m = 0;
        for (const ca of ab.A.comps) for (const cb of ab.B.comps) for (const cc of cd.A.comps) for (const cdd of cd.B.comps) {
          let s = 0;
          const t1m = ca[0] + cb[0], u1m = ca[1] + cb[1], v1m = ca[2] + cb[2];
          const t2m = cc[0] + cdd[0], u2m = cc[1] + cdd[1], v2m = cc[2] + cdd[2];
          for (let t = 0; t <= t1m; t++) { const e1 = x.E[0][ca[0]][cb[0]][t]; if (!e1) continue;
            for (let u = 0; u <= u1m; u++) { const e2 = e1 * x.E[1][ca[1]][cb[1]][u]; if (!e2) continue;
              for (let v = 0; v <= v1m; v++) { const e3 = e2 * x.E[2][ca[2]][cb[2]][v]; if (!e3) continue;
                for (let tt = 0; tt <= t2m; tt++) { const f1 = y.E[0][cc[0]][cdd[0]][tt]; if (!f1) continue;
                  for (let uu = 0; uu <= u2m; uu++) { const f2 = f1 * y.E[1][cc[1]][cdd[1]][uu]; if (!f2) continue;
                    for (let vv = 0; vv <= v2m; vv++) { const f3 = f2 * y.E[2][cc[2]][cdd[2]][vv]; if (!f3) continue;
                      s += e3 * f3 * ((tt + uu + vv) & 1 ? -1 : 1) * R(t + tt, u + uu, v + vv);
                    } } } } } }
          acc[m++] += pref * s;
        }
      }
      let m = 0;
      for (let a = 0; a < nA; a++) for (let b = 0; b < nB; b++) for (let c = 0; c < nC; c++) for (let d = 0; d < nD; d++) {
        const i = ab.A.f0 + a, j = ab.B.f0 + b, k = cd.A.f0 + c, l = cd.B.f0 + d, v = acc[m++];
        for (const [w, x2, y2, z2] of [[i, j, k, l], [j, i, k, l], [i, j, l, k], [j, i, l, k], [k, l, i, j], [l, k, i, j], [k, l, j, i], [l, k, j, i]]) ERI[id(w, x2, y2, z2)] = v;
      }
    }
    return ERI;
  }

  // ---------- linear algebra ----------
  function zeros(n, m) { m = m ?? n; const A = []; for (let i = 0; i < n; i++) A.push(new Array(m).fill(0)); return A; }
  function jacobi(Ain) {
    const n = Ain.length, A = Ain.map(r => r.slice()), V = A.map((r, i) => r.map((_, j) => (i === j ? 1 : 0)));
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += A[p][q] * A[p][q];
      if (off < 1e-26) break;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-300) continue;
        const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) { const x = A[k][p], y = A[k][q]; A[k][p] = c * x - s * y; A[k][q] = s * x + c * y; }
        for (let k = 0; k < n; k++) { const x = A[p][k], y = A[q][k]; A[p][k] = c * x - s * y; A[q][k] = s * x + c * y; }
        for (let k = 0; k < n; k++) { const x = V[k][p], y = V[k][q]; V[k][p] = c * x - s * y; V[k][q] = s * x + c * y; }
      }
    }
    return { vals: A.map((r, i) => r[i]), vecs: V };
  }
  const mm = (A, B) => { const n = A.length, m = B[0].length, k = B.length, C = zeros(n, m); for (let i = 0; i < n; i++) for (let l = 0; l < k; l++) { const a = A[i][l]; if (!a) continue; const Bl = B[l], Ci = C[i]; for (let j = 0; j < m; j++) Ci[j] += a * Bl[j]; } return C; };
  const tr = A => A[0].map((_, j) => A.map(r => r[j]));
  const add = (A, B, s = 1) => A.map((r, i) => r.map((v, j) => v + s * B[i][j]));
  function solveLinear(M, b) {
    const n = b.length, A = M.map((r, i) => [...r, b[i]]);
    for (let c = 0; c < n; c++) {
      let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
      [A[c], A[p]] = [A[p], A[c]];
      if (Math.abs(A[c][c]) < 1e-14) return null;
      for (let r = 0; r < n; r++) if (r !== c) { const f = A[r][c] / A[c][c]; for (let k = c; k <= n; k++) A[r][k] -= f * A[c][k]; }
    }
    return A.map((r, i) => r[n] / r[i]);
  }

  // ---------- integrals bundle for a geometry ----------
  function setup(atoms) {
    const basis = buildBasis(atoms);
    const { S, T, V } = oneElectron(basis, atoms);
    const H = add(T, V);
    const ERI = twoElectron(basis);
    let Enuc = 0;
    for (let i = 0; i < atoms.length; i++) for (let j = 0; j < i; j++) {
      const d = Math.hypot(...[0, 1, 2].map(k => atoms[i].pos[k] - atoms[j].pos[k]));
      Enuc += ZNUM[atoms[i].sym] * ZNUM[atoms[j].sym] / d;
    }
    // symmetric orthogonaliser
    const e = jacobi(S);
    const X = mm(mm(e.vecs, e.vals.map((v, i) => e.vals.map((_, j) => (i === j ? 1 / Math.sqrt(v) : 0)))), tr(e.vecs));
    // symmetry blocks for linear molecules on the z axis: sigma (s, pz), pi-x, pi-y
    const linearZ = atoms.every(a => Math.abs(a.pos[0]) < 1e-12 && Math.abs(a.pos[1]) < 1e-12);
    const blocks = linearZ
      ? [['σ', basis.funcs.map((f, i) => (f.type === 's' || f.type === 'pz' ? i : -1)).filter(i => i >= 0)],
         ['πx', basis.funcs.map((f, i) => (f.type === 'px' ? i : -1)).filter(i => i >= 0)],
         ['πy', basis.funcs.map((f, i) => (f.type === 'py' ? i : -1)).filter(i => i >= 0)]].filter(b => b[1].length)
      : [['all', basis.funcs.map((_, i) => i)]];
    return { atoms, basis, S, H, ERI, Enuc, X, blocks, n: basis.n };
  }

  function JK(ints, D) {
    const n = ints.n, E = ints.ERI, J = zeros(n), K = zeros(n);
    for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
      let jv = 0, kv = 0;
      const b1 = (i * n + j) * n * n;
      for (let k = 0; k < n; k++) for (let l = 0; l < n; l++) {
        const d = D[k][l]; if (!d) continue;
        jv += d * E[b1 + k * n + l];
        kv += d * E[((i * n + k) * n + j) * n + l];
      }
      J[i][j] = J[j][i] = jv; K[i][j] = K[j][i] = kv;
    }
    return { J, K };
  }
  const density = (C, occIdx) => { const n = C.length, D = zeros(n); for (const k of occIdx) for (let i = 0; i < n; i++) { const ci = C[i][k]; if (!ci) continue; for (let j = 0; j < n; j++) D[i][j] += ci * C[j][k]; } return D; };
  const fracDensity = (C, occ) => { const n = C.length, D = zeros(n); occ.forEach((w, k) => { if (!w) return; for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) D[i][j] += w * C[i][k] * C[j][k]; }); return D; };
  const dot = (A, B) => { let s = 0; for (let i = 0; i < A.length; i++) for (let j = 0; j < A.length; j++) s += A[i][j] * B[i][j]; return s; };

  // diagonalise F (AO) by symmetry blocks in the orthogonal basis
  function diagF(ints, F) {
    const n = ints.n, C = zeros(n), eps = [], symOf = [];
    const Xt = tr(ints.X);
    const Fp = mm(mm(Xt, F), ints.X);
    let col = 0;
    for (const [name, idx] of ints.blocks) {
      const sub = idx.map(i => idx.map(j => Fp[i][j]));
      const e = jacobi(sub);
      const order = e.vals.map((v, k) => k).sort((a, b) => e.vals[a] - e.vals[b]);
      for (const k of order) {
        const vec = new Array(n).fill(0); idx.forEach((i, r) => (vec[i] = e.vecs[r][k]));
        for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += ints.X[i][j] * vec[j]; C[i][col] = s; }
        eps.push(e.vals[k]); symOf.push(name); col++;
      }
    }
    // sort columns by energy (stable: σ, πx, πy order for ties)
    const order = eps.map((_, k) => k).sort((a, b) => (eps[a] - eps[b]) || (a - b));
    return { C: C.map(r => order.map(k => r[k])), eps: order.map(k => eps[k]), sym: order.map(k => symOf[k]) };
  }

  // Aufbau occupation that keeps degenerate πx/πy partners together.
  // Orbitals are grouped into units: each σ (or 'all') orbital alone; the k-th πx with the k-th πy.
  function aufbau(eps, nb, na, sym) {
    const units = [];
    const px = [], py = [];
    eps.forEach((e, k) => { if (sym[k] === 'πx') px.push(k); else if (sym[k] === 'πy') py.push(k); else units.push({ E: e, orbs: [k] }); });
    for (let i = 0; i < Math.max(px.length, py.length); i++) {
      const o = [px[i], py[i]].filter(k => k !== undefined);
      units.push({ E: o.reduce((s, k) => s + eps[k], 0) / o.length, orbs: o });
    }
    units.sort((a, b) => a.E - b.E);
    const flat = units.flatMap(u => u.orbs);
    return { core: flat.slice(0, nb), open: flat.slice(nb, na) };
  }

  // Maximum-overlap occupation: follow the previously occupied orbitals (keeps a scan on one state)
  function momOcc(Cold, occOld, Cnew, S, nb, na) {
    const n = Cnew.length, SC = mm(S, Cnew);
    const proj = cols => Cnew[0].map((_, k) => { let s = 0; for (const j of cols) { let d = 0; for (let i = 0; i < n; i++) d += Cold[i][j] * SC[i][k]; s += d * d; } return s; });
    const pc = proj(occOld.core), idx = pc.map((_, k) => k).sort((a, b) => pc[b] - pc[a]);
    const core = idx.slice(0, nb), rest = idx.slice(nb);
    const po = proj(occOld.open);
    rest.sort((a, b) => po[b] - po[a]);
    return { core: core.sort((a, b) => a - b), open: rest.slice(0, na - nb).sort((a, b) => a - b) };
  }

  // ---------- SCF ----------
  function scf(ints, nElec, mult, opts = {}) {
    const n = ints.n, nopen = mult - 1;
    if ((nElec - nopen) % 2 !== 0 || nopen > nElec) throw new Error('Spin multiplicity incompatible with electron count');
    const nb = (nElec - nopen) / 2, na = nb + nopen;
    if (na > n) throw new Error('Too many electrons for the basis');
    const S = ints.S, H = ints.H;
    let orb;
    if (opts.guessC) orb = { C: opts.guessC.map(r => r.slice()), eps: opts.guessEps.slice(), sym: opts.guessSym.slice() };
    else {
      // superposition of spherically averaged atomic densities (SAD)
      const D0 = zeros(n);
      ints.atoms.forEach((at, ia) => {
        const Da = atomData(at.sym).D, off = ints.basis.funcs.findIndex(f => f.atom === ia);
        for (let i = 0; i < Da.length; i++) for (let j = 0; j < Da.length; j++) D0[off + i][off + j] = Da[i][j];
      });
      const { J, K } = JK(ints, D0);
      orb = diagF(ints, add(add(H, J), K, -0.5));
    }
    const mom = !!(opts.mom && opts.guessOcc);
    let occ = mom
      ? { core: opts.guessOcc.map((o, k) => (o === 2 ? k : -1)).filter(k => k >= 0), open: opts.guessOcc.map((o, k) => (o === 1 ? k : -1)).filter(k => k >= 0) }
      : aufbau(orb.eps, nb, na, orb.sym);
    if (mom && (occ.core.length !== nb || occ.open.length !== na - nb)) occ = aufbau(orb.eps, nb, na, orb.sym);
    const diisF = [], diisE = [];
    let Eold = 0, E = 0, converged = false, it = 0, Fa, Fb, Feff, Da, Db;
    const history = [];
    for (it = 1; it <= (opts.maxIter ?? 200); it++) {
      Db = density(orb.C, occ.core); Da = add(Db, density(orb.C, occ.open));
      const Dt = add(Da, Db);
      const { J, K: Ka } = JK(ints, Da);
      const Kb = nopen ? JK(ints, Db).K : Ka;
      const { J: Jb } = nopen ? JK(ints, Db) : { J: J };
      const Jt = nopen ? add(J, Jb) : J.map(r => r.map(v => 2 * v));
      Fa = add(add(H, Jt), Ka, -1); Fb = add(add(H, Jt), Kb, -1);
      E = 0.5 * (dot(Da, add(H, Fa)) + dot(Db, add(H, Fb))) + ints.Enuc;
      if (!nopen) Feff = Fa;
      else {
        // Roothaan effective Fock (same convention as PySCF ROHF)
        const C = orb.C, Ct = tr(C);
        const fa = mm(mm(Ct, Fa), C), fb = mm(mm(Ct, Fb), C), fc = add(fa, fb).map(r => r.map(v => v / 2));
        const cls = new Array(n).fill(2); occ.core.forEach(k => (cls[k] = 0)); occ.open.forEach(k => (cls[k] = 1));
        const f = zeros(n);
        for (let p = 0; p < n; p++) for (let q = 0; q < n; q++) {
          const a = cls[p], b = cls[q];
          f[p][q] = (a === 0 && b === 1) || (a === 1 && b === 0) ? fb[p][q] : (a === 1 && b === 2) || (a === 2 && b === 1) ? fa[p][q] : fc[p][q];
        }
        const SC = mm(S, C);
        Feff = mm(mm(SC, f), tr(SC));
      }
      // DIIS error FDS - SDF in orthogonal basis
      const FDS = mm(mm(Feff, Dt), S), err = add(FDS, tr(FDS), -1);
      const errO = mm(mm(tr(ints.X), err), ints.X);
      let emax = 0; for (const r of errO) for (const v of r) emax = Math.max(emax, Math.abs(v));
      history.push(E);
      if (it > 1 && Math.abs(E - Eold) < 1e-10 && emax < 1e-7) { converged = true; break; }
      Eold = E;
      diisF.push(Feff); diisE.push(errO);
      if (diisF.length > 8) { diisF.shift(); diisE.shift(); }
      let Fuse = Feff;
      if (diisF.length >= 2 && it > 1) {
        const m = diisF.length, B = zeros(m + 1), rhs = new Array(m + 1).fill(0);
        for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) B[i][j] = dot(diisE[i], diisE[j]);
        for (let i = 0; i < m; i++) { B[i][m] = B[m][i] = -1; }
        rhs[m] = -1;
        const cvec = solveLinear(B, rhs);
        if (cvec) { Fuse = zeros(n); for (let k = 0; k < m; k++) Fuse = add(Fuse, diisF[k], cvec[k]); }
      }
      const prev = orb;
      orb = diagF(ints, Fuse);
      occ = mom && it > 1 ? momOcc(prev.C, occ, orb.C, S, nb, na) : (mom ? momOcc(prev.C, occ, orb.C, S, nb, na) : aufbau(orb.eps, nb, na, orb.sym));
    }
    // final canonical orbitals from the converged effective Fock
    { const prev = orb; orb = diagF(ints, Feff); occ = mom ? momOcc(prev.C, occ, orb.C, S, nb, na) : aufbau(orb.eps, nb, na, orb.sym); }
    if (nopen) {
      // alpha-spin semicanonicalisation: diagonalise F(alpha) within each occupation class and symmetry,
      // so orbital energies approximate alpha ionisation / attachment energies (Koopmans-like).
      const cls = new Array(n).fill('v'); occ.core.forEach(k => (cls[k] = 'c')); occ.open.forEach(k => (cls[k] = 'o'));
      const groups = {};
      orb.eps.forEach((_, k) => { const key = cls[k] + '|' + orb.sym[k]; (groups[key] = groups[key] || []).push(k); });
      const Cn = orb.C.map(r => r.slice()), en = orb.eps.slice();
      for (const key in groups) {
        const ks = groups[key];
        const Csub = orb.C.map(r => ks.map(k => r[k]));
        const fsub = mm(mm(tr(Csub), Fa), Csub);
        const e = jacobi(fsub);
        const ord = e.vals.map((_, i) => i).sort((a, b) => e.vals[a] - e.vals[b]);
        ord.forEach((src, dst) => {
          const k = ks[dst]; en[k] = e.vals[src];
          for (let i = 0; i < n; i++) { let v = 0; for (let r = 0; r < ks.length; r++) v += Csub[i][r] * e.vecs[r][src]; Cn[i][k] = v; }
        });
      }
      const order = en.map((_, k) => k).sort((a, b) => (en[a] - en[b]) || (a - b));
      const map = {}; order.forEach((k, i) => (map[k] = i));
      orb = { C: Cn.map(r => order.map(k => r[k])), eps: order.map(k => en[k]), sym: order.map(k => orb.sym[k]) };
      occ = { core: occ.core.map(k => map[k]), open: occ.open.map(k => map[k]) };
    }
    const occVec = orb.eps.map((_, k) => (occ.core.includes(k) ? 2 : occ.open.includes(k) ? 1 : 0));
    return { Fa, Fb, E, Enuc: ints.Enuc, C: orb.C, eps: orb.eps, sym: orb.sym, occ: occVec, nAlpha: na, nBeta: nb, mult, iterations: it, converged, history, Dtot: add(Da, Db) };
  }

  // spherically averaged (fractional occupation) restricted calculation for a free atom: orbital energies for diagrams
  const atomCache = {};
  function atomData(sym) {
    if (atomCache[sym]) return atomCache[sym];
    return (atomCache[sym] = atomFrac(sym));
  }
  function atomOrbitalEnergies(sym) { return atomData(sym).levels; }
  function atomFrac(sym) {
    const ints = setup([{ sym, pos: [0, 0, 0] }]);
    const N = ZNUM[sym];
    let orb = diagF(ints, ints.H), E = 0, Eold = 1, Dlast = null;
    const occFor = eps => {
      const occ = new Array(eps.length).fill(0); let left = N;
      for (let i = 0; i < eps.length && left > 0;) {
        let j = i; while (j < eps.length && Math.abs(eps[j] - eps[i]) < 1e-5) j++;
        const g = j - i, put = Math.min(left, 2 * g);
        for (let k = i; k < j; k++) occ[k] = put / g;
        left -= put; i = j;
      }
      return occ;
    };
    for (let it = 0; it < 200; it++) {
      const occ = occFor(orb.eps), D = fracDensity(orb.C, occ); Dlast = D;
      const { J, K } = JK(ints, D);
      const F = add(add(ints.H, J), K, -0.5);
      E = 0.5 * dot(D, add(ints.H, F));
      orb = diagF(ints, F);
      if (Math.abs(E - Eold) < 1e-11) break;
      Eold = E;
    }
    // map orbitals to shells by dominant function
    const out = {};
    orb.eps.forEach((e, k) => {
      let best = 0, bi = 0; for (let i = 0; i < ints.n; i++) { const v = Math.abs(orb.C[i][k]); if (v > best) { best = v; bi = i; } }
      const sh = ints.basis.funcs[bi].shell;
      if (!(sh in out)) out[sh] = e * HARTREE_EV;
    });
    // 1s vs 2s dominance can be ambiguous for the 2s orbital; assign by energy order within s shells
    const sE = orb.eps.filter((_, k) => orb.sym[k] === 'σ').map(e => e * HARTREE_EV);
    const sShells = ints.basis.shells.filter(s => s.l === 0).map(s => s.shell);
    const sOrbs = orb.eps.map((e, k) => [e, k]).filter(([e, k]) => {
      // s-type orbital: no pz character
      let pz = 0; for (let i = 0; i < ints.n; i++) if (ints.basis.funcs[i].type === 'pz') pz += orb.C[i][k] ** 2; return pz < 1e-8;
    }).map(([e]) => e * HARTREE_EV);
    sShells.forEach((sh, i) => { if (sOrbs[i] !== undefined) out[sh] = sOrbs[i]; });
    void sE;
    // Orbital energies for the diagram: alpha-spin Fock eigenvalues of the high-spin ROHF atom
    // (Koopmans-like: H 1s = -12.7 eV). p shell = average over its alpha-occupied orbitals (or all if none).
    atomCache[sym] = { D: Dlast, levels: null }; // SAD guess for the atom's own ROHF run
    const ro = scf(ints, N, HUND[sym]);
    const lev = {}, groups = {};
    ro.eps.forEach((_, k) => {
      const ck = ro.C.map(r => r[k]);
      let e = 0; for (let i = 0; i < ints.n; i++) for (let j = 0; j < ints.n; j++) e += ck[i] * ro.Fa[i][j] * ck[j];
      let best = 0, bi = 0; for (let i = 0; i < ints.n; i++) { const v = Math.abs(ck[i]); if (v > best) { best = v; bi = i; } }
      const f = ints.basis.funcs[bi];
      const isP = f.type !== 's';
      const key = isP ? f.shell : 's';
      (groups[key] = groups[key] || []).push({ e, occA: ro.occ[k] > 0 });
    });
    // s orbitals in energy order -> 1s, 2s, 3s
    const sList = (groups.s || []).sort((a, b) => a.e - b.e);
    ints.basis.shells.filter(x => x.l === 0).forEach((x, i) => { if (sList[i]) lev[x.shell] = sList[i].e * HARTREE_EV; });
    for (const key in groups) if (key !== 's') {
      const occd = groups[key].filter(g => g.occA), use = occd.length ? occd : groups[key];
      lev[key] = use.reduce((t, g) => t + g.e, 0) / use.length * HARTREE_EV;
    }
    void out;
    return { levels: lev, D: Dlast };
  }

  const HUND = { H: 2, He: 1, Li: 2, Be: 1, B: 2, C: 3, N: 4, O: 3, F: 2, Ne: 1, Cl: 2 };
  function atomEnergy(sym) {
    const ints = setup([{ sym, pos: [0, 0, 0] }]);
    return scf(ints, ZNUM[sym], HUND[sym]).E;
  }

  // ---------- public: diatomic single point ----------
  function diatomic(symA, symB, R_ang, charge, mult, guess, follow) {
    const Rb = R_ang * BOHR_PER_ANG;
    const atoms = [{ sym: symA, pos: [0, 0, 0] }, { sym: symB, pos: [0, 0, Rb] }];
    const ints = setup(atoms);
    const nElec = ZNUM[symA] + ZNUM[symB] - charge;
    let res = null;
    if (guess && guess.n === ints.n) {
      try { res = scf(ints, nElec, mult, { guessC: guess.C, guessEps: guess.eps, guessSym: guess.sym, guessOcc: guess.occ, mom: !!follow }); } catch (e) { res = null; }
      if (res && !res.converged) res = null;
    }
    if (!res) res = scf(ints, nElec, mult);
    res.basis = ints.basis.funcs.map(f => ({ atom: f.atom, sym: f.sym, shell: f.shell, type: f.type, label: f.label, z0: f.center[2], exps: f.exps, coefs: f.coefs, l: f.type === 's' ? 0 : 1 }));
    res.S = ints.S; res.Rb = Rb; res.R = R_ang; res.nElec = nElec; res.n = ints.n;
    delete res.Dtot;
    return res;
  }

  // one point of a scan: follow the previous state AND search from scratch; keep the lower converged solution
  function scanPoint(symA, symB, R_ang, charge, mult, guess) {
    const Rb = R_ang * BOHR_PER_ANG;
    const ints = setup([{ sym: symA, pos: [0, 0, 0] }, { sym: symB, pos: [0, 0, Rb] }]);
    const nElec = ZNUM[symA] + ZNUM[symB] - charge;
    const cands = [];
    if (guess && guess.n === ints.n) { try { const r = scf(ints, nElec, mult, { guessC: guess.C, guessEps: guess.eps, guessSym: guess.sym, guessOcc: guess.occ, mom: true }); if (r.converged) cands.push(r); } catch (e) { /* ignore */ } }
    try { const r = scf(ints, nElec, mult); cands.push(r); } catch (e) { /* ignore */ }
    const conv = cands.filter(r => r.converged);
    const best = (conv.length ? conv : cands).sort((a, b) => a.E - b.E)[0];
    best.n = ints.n; best.R = R_ang;
    return best;
  }

  // potential energy curve: lowest converged solution at each R from a fresh start plus sweeps in both directions
  function curve(symA, symB, charge, mult, Rs, onProgress) {
    const nElec = ZNUM[symA] + ZNUM[symB] - charge;
    const ints = Rs.map((R, i) => { const x = setup([{ sym: symA, pos: [0, 0, 0] }, { sym: symB, pos: [0, 0, R * BOHR_PER_ANG] }]); onProgress && onProgress(i + 1, Rs.length); return x; });
    const best = Rs.map(() => null);
    const keep = (i, r) => { if (r && r.converged && (!best[i] || r.E < best[i].E - 1e-9)) { best[i] = r; return true; } return false; };
    Rs.forEach((_, i) => { try { keep(i, scf(ints[i], nElec, mult)); } catch (e) { /* ignore */ } });
    const follow = (i, g) => { try { return scf(ints[i], nElec, mult, { guessC: g.C, guessEps: g.eps, guessSym: g.sym, guessOcc: g.occ, mom: true }); } catch (e) { return null; } };
    for (let pass = 0; pass < 2; pass++) {
      let changed = false;
      for (let i = 1; i < Rs.length; i++) if (best[i - 1]) changed = keep(i, follow(i, best[i - 1])) || changed;
      for (let i = Rs.length - 2; i >= 0; i--) if (best[i + 1]) changed = keep(i, follow(i, best[i + 1])) || changed;
      if (!changed) break;
    }
    return Rs.map((R, i) => ({ R, E: best[i] ? best[i].E : NaN, ok: !!best[i] }));
  }

  function optimize(symA, symB, R0, charge, mult, guess) {
    let R = R0, g = guess, last = null;
    const h = 0.005;
    for (let k = 0; k < 25; k++) {
      const e0 = diatomic(symA, symB, R, charge, mult, g, true); g = e0;
      const ep = diatomic(symA, symB, R + h, charge, mult, e0, true), em = diatomic(symA, symB, R - h, charge, mult, e0, true);
      const grad = (ep.E - em.E) / (2 * h), curv = (ep.E - 2 * e0.E + em.E) / (h * h);
      let step = curv > 1e-6 ? -grad / curv : -Math.sign(grad) * 0.1;
      step = Math.max(-0.15, Math.min(0.15, step));
      last = e0;
      R = Math.max(0.4, Math.min(4.0, R + step));
      if (Math.abs(step) < 1e-4) break;
    }
    const fin = diatomic(symA, symB, R, charge, mult, g, true);
    return { R, result: fin };
  }

  return { diatomic, scanPoint, curve, optimize, atomOrbitalEnergies, atomEnergy, setup, scf, buildBasis, HARTREE_EV, BOHR_PER_ANG, ZNUM, HUND };
})();
if (typeof module !== 'undefined') module.exports = HF;
