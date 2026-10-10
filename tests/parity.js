// node tests/parity.js  -> compares the browser solver with the Python reference on shared cases
const B = require("../web/beam.js");
const w = 10, L = 6, P = 40, EI = 20000; let fails = 0;
const close = (a, b, tol = 1e-3) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(b));
const check = (name, got, want) => { const ok = got.every((g, i) => close(g, want[i])); if (!ok) fails++; console.log(`${ok ? "PASS" : "FAIL"} ${name}: ${got.map(v => v.toFixed(3))} vs ${want.map(v => v.toFixed(3))}`); };
let r = B.solve({ spans: [{ length: L, EI, loads: [{ kind: "udl", w }] }], left: "fixed", right: "fixed" });
check("fixed-fixed UDL end moments", r.endMoments[0], [-w*L*L/12, w*L*L/12]);
r = B.solve({ spans: [{ length: L, EI, loads: [{ kind: "udl", w }] }], left: "fixed", right: "pinned" });
check("propped cantilever reactions", r.reactions, [5*w*L/8, 3*w*L/8]);
r = B.solve({ spans: [{ length: L, EI, loads: [{ kind: "udl", w }] }, { length: L, EI, loads: [{ kind: "udl", w }] }], left: "pinned", right: "pinned" });
check("two-span UDL support moment", [r.endMoments[0][1]], [w*L*L/8]);
check("two-span UDL reactions", r.reactions, [3*w*L/8, 10*w*L/8, 3*w*L/8]);
r = B.solve({ spans: [{ length: L, EI, loads: [] }], left: "fixed", right: "fixed", settlements: [0, 0.01] });
check("settlement moments", r.endMoments[0], [-6*EI*0.01/L/L, -6*EI*0.01/L/L]);
r = B.solve({ spans: [{ length: 6, EI, loads: [{ kind: "point", w: P, a: 2 }] }], left: "fixed", right: "fixed" });
check("off-centre point load FEM", r.endMoments[0], [-P*2*16/36, P*4*4/36]);
r = B.solve({ spans: [{ length: L, EI, loads: [{ kind: "point", w: P, a: L/2 }] }] });
check("simply supported deflection", [Math.max(...r.deflection[0].map(p => p.v))], [P*L**3/(48*EI)]);
r = B.solve({ spans: [{ length: L, EI, loads: [] }], overhangRight: { length: 2, EI, loads: [{ kind: "point", w: P, a: 2 }] } });
check("right overhang reactions", r.reactions, [-P*2/L, P + P*2/L]);
{ const tipD = r.deflection[r.deflection.length - 1]; check("overhang tip deflection", [tipD[tipD.length - 1].v], [P*2*2*(L+2)/(3*EI)]); }
process.exit(fails ? 1 : 0);
