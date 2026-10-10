"""
Continuous-beam analysis by the slope-deflection method.

Sign conventions
- x runs left to right from the first support (overhangs use their own local x).
- Loads act downward and are entered as positive numbers (kN, kN/m).
- Member end moments are clockwise-positive (slope-deflection convention).
- Bending moment along the beam is sagging-positive.
- Rotations are clockwise-positive; settlements and deflections are downward-positive (m).

Method
For a span i-j with length L and flexural rigidity EI:
    M_ij = FEM_ij + (2EI/L)(2θi + θj − 3ψ)
    M_ji = FEM_ji + (2EI/L)(2θj + θi − 3ψ),   ψ = (Δj − Δi)/L
Every support that can rotate (pinned ends and interior supports) contributes one
equilibrium equation  Σ member end moments = 0  (overhang moments enter as constants).
The resulting linear system is solved for the rotations, then end moments, reactions,
shear, moment and deflection follow from statics and integration of M/EI.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

SLICES = 2000  # distributed loads are integrated as this many point-load slices per span


@dataclass
class Load:
    kind: str            # "point" | "udl" | "linear"
    w: float = 0.0       # magnitude: kN for point; kN/m for distributed (start value for linear)
    a: float = 0.0       # start position from the span's left end (m); position for point loads
    b: float | None = None   # end position (m) for distributed loads; None = span end
    w2: float | None = None  # end value for linear loads (kN/m)


@dataclass
class Span:
    length: float
    EI: float = 1.0      # kN·m²
    loads: list = field(default_factory=list)


@dataclass
class Beam:
    spans: list
    left: str = "pinned"     # "fixed" | "pinned"
    right: str = "pinned"
    settlements: list | None = None          # per support, m (downward +)
    overhang_left: Span | None = None        # cantilever left of the first support (free end at far left)
    overhang_right: Span | None = None       # cantilever right of the last support


def point_loads(span: Span):
    """Return (positions, forces) with distributed loads sliced into point loads."""
    xs, ps = [], []
    for ld in span.loads:
        if ld.kind == "point":
            xs.append(ld.a); ps.append(ld.w)
        else:
            a = ld.a; b = span.length if ld.b is None else ld.b
            w1 = ld.w; w2 = ld.w if (ld.kind == "udl" or ld.w2 is None) else ld.w2
            n = SLICES; dx = (b - a) / n
            mids = a + (np.arange(n) + 0.5) * dx
            ws = w1 + (w2 - w1) * (mids - a) / (b - a) if b > a else np.zeros(n)
            xs.extend(mids.tolist()); ps.extend((ws * dx).tolist())
    return np.asarray(xs, float), np.asarray(ps, float)


def fixed_end_moments(span: Span):
    x, p = point_loads(span)
    L = span.length; b = L - x
    return float(-(p * x * b**2).sum() / L**2), float((p * x**2 * b).sum() / L**2)


def overhang_moment(oh: Span, side: str):
    """Member end moment (clockwise +) that a cantilever overhang applies at its support."""
    x, p = point_loads(oh)
    if side == "left":   # local x from the free end; support at x = length
        return float((p * (oh.length - x)).sum())
    return float(-(p * x).sum())   # right overhang: local x from the support


def solve(beam: Beam, n_points: int = 201):
    spans = beam.spans
    n_sup = len(spans) + 1
    D = np.zeros(n_sup) if beam.settlements is None else np.asarray(beam.settlements, float)
    fem = [fixed_end_moments(s) for s in spans]
    psi = [(D[i + 1] - D[i]) / s.length for i, s in enumerate(spans)]

    free = [i for i in range(n_sup) if not ((i == 0 and beam.left == "fixed") or (i == n_sup - 1 and beam.right == "fixed"))]
    idx = {node: k for k, node in enumerate(free)}
    K = np.zeros((len(free), len(free))); F = np.zeros(len(free))

    ext = np.zeros(n_sup)  # constant end moments from overhangs, added at the support
    if beam.overhang_left is not None:
        ext[0] += overhang_moment(beam.overhang_left, "left")
    if beam.overhang_right is not None:
        ext[-1] += overhang_moment(beam.overhang_right, "right")

    for s_i, s in enumerate(spans):
        k = 2 * s.EI / s.length
        i, j = s_i, s_i + 1
        # M_ij = FEM_ij + k(2θi + θj − 3ψ)  contributes to equilibrium at i
        # M_ji = FEM_ji + k(2θj + θi − 3ψ)  contributes to equilibrium at j
        for node, near, far, f in ((i, i, j, fem[s_i][0]), (j, j, i, fem[s_i][1])):
            if node not in idx:
                continue
            r = idx[node]
            K[r, idx[near]] += 2 * k
            if far in idx:
                K[r, idx[far]] += k
            F[r] -= f - 3 * k * psi[s_i]
    for node in free:
        F[idx[node]] -= ext[node]

    theta = np.zeros(n_sup)
    if free:
        theta[free] = np.linalg.solve(K, F)

    end_moments = []
    for s_i, s in enumerate(spans):
        k = 2 * s.EI / s.length; i, j = s_i, s_i + 1
        Mij = fem[s_i][0] + k * (2 * theta[i] + theta[j] - 3 * psi[s_i])
        Mji = fem[s_i][1] + k * (2 * theta[j] + theta[i] - 3 * psi[s_i])
        end_moments.append((float(Mij), float(Mji)))

    # ---- reactions, diagrams
    reactions = np.zeros(n_sup)
    xs_all, V_all, M_all, v_all = [], [], [], []
    offset = 0.0

    def cantilever(oh: Span, side: str, x0: float):
        x, p = point_loads(oh)
        L = oh.length
        loc = np.linspace(0, L, max(11, n_points // 2))
        if side == "left":   # free end at loc=0, support at loc=L
            V = np.array([-(p[x < t]).sum() for t in loc])
            M = np.array([-(p[x < t] * (t - x[x < t])).sum() for t in loc])
        else:                # support at loc=0, free end at loc=L
            V = np.array([(p[x > t]).sum() for t in loc])
            M = np.array([-(p[x > t] * (x[x > t] - t)).sum() for t in loc])
        return loc + x0, V, M, float(p.sum())

    if beam.overhang_left is not None:
        xc, Vc, Mc, P = cantilever(beam.overhang_left, "left", -beam.overhang_left.length)
        reactions[0] += P
        xs_all.append(xc); V_all.append(Vc); M_all.append(Mc)

    for s_i, s in enumerate(spans):
        x, p = point_loads(s); L = s.length
        Mij, Mji = end_moments[s_i]
        Vi = ((p * (L - x)).sum() - Mij - Mji) / L
        Vj = p.sum() - Vi
        reactions[s_i] += Vi; reactions[s_i + 1] += Vj
        loc = np.linspace(0, L, n_points)
        V = np.array([Vi - p[x < t].sum() for t in loc])
        M = np.array([Mij + Vi * t - (p[x < t] * (t - x[x < t])).sum() for t in loc])
        # deflection: EI v'' = −M, v(0)=Δi, v'(0)=θi (downward +)
        curv = -M / s.EI
        slope = theta[s_i] + np.concatenate([[0], np.cumsum((curv[1:] + curv[:-1]) / 2 * np.diff(loc))])
        v = D[s_i] + np.concatenate([[0], np.cumsum((slope[1:] + slope[:-1]) / 2 * np.diff(loc))])
        xs_all.append(loc + offset); V_all.append(V); M_all.append(M); v_all.append((loc + offset, v))
        offset += L

    if beam.overhang_right is not None:
        xc, Vc, Mc, P = cantilever(beam.overhang_right, "right", offset)
        reactions[-1] += P
        xs_all.append(xc); V_all.append(Vc); M_all.append(Mc)

    X = np.concatenate(xs_all); Vd = np.concatenate(V_all); Md = np.concatenate(M_all)
    return dict(theta=theta.tolist(), end_moments=end_moments, reactions=reactions.tolist(),
                x=X.tolist(), shear=Vd.tolist(), moment=Md.tolist(),
                deflection=[dict(x=a.tolist(), v=b.tolist()) for a, b in v_all],
                max_sagging=float(Md.max()), max_hogging=float(Md.min()))
