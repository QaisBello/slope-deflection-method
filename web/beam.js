/* Continuous-beam solver (slope-deflection method). Mirrors beam/solver.py line for line.
   Units: kN, m, kN·m. Loads downward +. End moments clockwise +. Bending moment sagging +.
   Rotations clockwise +. Settlement / deflection downward +. */
(function (root) {
  const SLICES = 1000;

  function pointLoads(span) {
    const xs = [], ps = [];
    for (const ld of span.loads || []) {
      if (ld.kind === "point") { xs.push(+ld.a); ps.push(+ld.w); continue; }
      const a = +ld.a || 0, b = ld.b == null ? span.length : +ld.b;
      const w1 = +ld.w, w2 = ld.kind === "udl" || ld.w2 == null ? +ld.w : +ld.w2;
      if (b <= a) continue;
      const dx = (b - a) / SLICES;
      for (let k = 0; k < SLICES; k++) { const m = a + (k + 0.5) * dx; xs.push(m); ps.push((w1 + (w2 - w1) * (m - a) / (b - a)) * dx); }
    }
    return { xs, ps };
  }

  function fem(span) {
    const { xs, ps } = pointLoads(span), L = span.length; let l = 0, r = 0;
    xs.forEach((x, k) => { const b = L - x; l -= ps[k] * x * b * b / (L * L); r += ps[k] * x * x * b / (L * L); });
    return [l, r];
  }

  function overhangMoment(oh, side) {
    const { xs, ps } = pointLoads(oh); let m = 0;
    xs.forEach((x, k) => { m += side === "left" ? ps[k] * (oh.length - x) : -ps[k] * x; });
    return m;
  }

  function gauss(A, b) { // solve A x = b with partial pivoting
    const n = b.length, M = A.map((r, i) => [...r, b[i]]);
    for (let c = 0; c < n; c++) {
      let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      [M[c], M[p]] = [M[p], M[c]];
      if (Math.abs(M[c][c]) < 1e-12) throw new Error("The structure is unstable: check the supports.");
      for (let r = c + 1; r < n; r++) { const f = M[r][c] / M[c][c]; for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k]; }
    }
    const x = Array(n).fill(0);
    for (let r = n - 1; r >= 0; r--) { let s = M[r][n]; for (let k = r + 1; k < n; k++) s -= M[r][k] * x[k]; x[r] = s / M[r][r]; }
    return x;
  }

  // integrate EI v'' = −M outward from a support at index s (v = v0, v' = s0); pts sorted by x
  function integrate(pts, EI, s, v0, s0) {
    const out = pts.map(p => ({ x: p.x, v: 0 })); out[s].v = v0;
    for (const dir of [1, -1]) { let v = v0, sl = s0;
      for (let q = s + dir; q >= 0 && q < pts.length; q += dir) { const h = pts[q].x - pts[q - dir].x, c0 = -pts[q - dir].M / EI, c1 = -pts[q].M / EI, s1 = sl + (c0 + c1) / 2 * h; v += (sl + s1) / 2 * h; sl = s1; out[q].v = v; } }
    return out;
  }

  function solve(beam, nPoints = 201) {
    const spans = beam.spans, nSup = spans.length + 1;
    const D = beam.settlements ? beam.settlements.map(Number) : Array(nSup).fill(0);
    const F_ = spans.map(fem), psi = spans.map((s, i) => (D[i + 1] - D[i]) / s.length);
    const free = []; for (let i = 0; i < nSup; i++) if (!((i === 0 && beam.left === "fixed") || (i === nSup - 1 && beam.right === "fixed"))) free.push(i);
    const idx = new Map(free.map((n, k) => [n, k]));
    const K = free.map(() => Array(free.length).fill(0)), F = Array(free.length).fill(0);
    const ext = Array(nSup).fill(0);
    if (beam.overhangLeft) ext[0] += overhangMoment(beam.overhangLeft, "left");
    if (beam.overhangRight) ext[nSup - 1] += overhangMoment(beam.overhangRight, "right");
    const equations = [];
    spans.forEach((s, si) => {
      const k = 2 * s.EI / s.length, i = si, j = si + 1;
      [[i, i, j, F_[si][0]], [j, j, i, F_[si][1]]].forEach(([node, near, far, f]) => {
        if (!idx.has(node)) return; const r = idx.get(node);
        K[r][idx.get(near)] += 2 * k; if (idx.has(far)) K[r][idx.get(far)] += k;
        F[r] -= f - 3 * k * psi[si];
      });
    });
    free.forEach(n => { F[idx.get(n)] -= ext[n]; });
    const theta = Array(nSup).fill(0);
    if (free.length) gauss(K, F).forEach((t, k) => { theta[free[k]] = t; });

    const endMoments = spans.map((s, si) => { const k = 2 * s.EI / s.length, i = si, j = si + 1;
      return [F_[si][0] + k * (2 * theta[i] + theta[j] - 3 * psi[si]), F_[si][1] + k * (2 * theta[j] + theta[i] - 3 * psi[si])]; });

    const reactions = Array(nSup).fill(0), pts = [], defl = [];
    const cant = (oh, side, x0) => {
      const { xs, ps } = pointLoads(oh), L = oh.length, n = Math.max(11, Math.floor(nPoints / 2)), out = [];
      for (let q = 0; q < n; q++) { const t = L * q / (n - 1); let V = 0, M = 0;
        xs.forEach((x, k) => { if (side === "left" && x < t) { V -= ps[k]; M -= ps[k] * (t - x); } if (side === "right" && x > t) { V += ps[k]; M -= ps[k] * (x - t); } });
        out.push({ x: x0 + t, V, M }); }
      return { out, P: ps.reduce((a, b) => a + b, 0) };
    };
    if (beam.overhangLeft) { const c = cant(beam.overhangLeft, "left", -beam.overhangLeft.length); reactions[0] += c.P; pts.push(...c.out); }
    let off = 0;
    spans.forEach((s, si) => {
      const { xs, ps } = pointLoads(s), L = s.length, [Mij, Mji] = endMoments[si];
      let sumM = 0, sumP = 0; xs.forEach((x, k) => { sumM += ps[k] * (L - x); sumP += ps[k]; });
      const Vi = (sumM - Mij - Mji) / L, Vj = sumP - Vi; reactions[si] += Vi; reactions[si + 1] += Vj;
      const loc = [], Ms = [];
      for (let q = 0; q < nPoints; q++) { const t = L * q / (nPoints - 1); let V = Vi, M = Mij + Vi * t;
        xs.forEach((x, k) => { if (x < t) { V -= ps[k]; M -= ps[k] * (t - x); } });
        loc.push(t); Ms.push(M); pts.push({ x: off + t, V, M, span: si }); }
      let slope = theta[si], v = D[si]; const vs = [{ x: off, v }];
      for (let q = 1; q < nPoints; q++) { const h = loc[q] - loc[q - 1], c0 = -Ms[q - 1] / s.EI, c1 = -Ms[q] / s.EI, s1 = slope + (c0 + c1) / 2 * h; v += (slope + s1) / 2 * h; slope = s1; vs.push({ x: off + loc[q], v }); }
      defl.push(vs); off += L;
    });
    if (beam.overhangRight) { const c = cant(beam.overhangRight, "right", off); reactions[nSup - 1] += c.P; pts.push(...c.out);
      defl.push(integrate(c.out, beam.overhangRight.EI, 0, D[nSup - 1], theta[nSup - 1])); }
    if (beam.overhangLeft) { const c = pts.filter(p => p.x <= 0 && p.span == null);
      defl.unshift(integrate(c, beam.overhangLeft.EI, c.length - 1, D[0], theta[0])); }
    const Ms = pts.map(p => p.M);
    return { theta, endMoments, reactions, points: pts, deflection: defl, fem: F_, psi, maxSagging: Math.max(...Ms), maxHogging: Math.min(...Ms), totalLength: off };
  }
  root.BeamSolver = { solve, fem, pointLoads };
  if (typeof module !== "undefined") module.exports = root.BeamSolver;
})(typeof window !== "undefined" ? window : globalThis);
