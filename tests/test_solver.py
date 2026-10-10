"""Checks the solver against closed-form textbook results."""
import math
import pytest
from beam.solver import Beam, Span, Load, solve

w, L, P, EI = 10.0, 6.0, 40.0, 20000.0
approx = lambda v, rel=1e-4: pytest.approx(v, rel=rel, abs=1e-6)


def test_fixed_fixed_udl():
    r = solve(Beam([Span(L, EI, [Load("udl", w)])], left="fixed", right="fixed"))
    Mij, Mji = r["end_moments"][0]
    assert Mij == approx(-w * L**2 / 12) and Mji == approx(w * L**2 / 12)
    assert r["max_sagging"] == approx(w * L**2 / 24, 1e-3)
    assert r["reactions"] == [approx(w * L / 2), approx(w * L / 2)]


def test_propped_cantilever_udl():
    r = solve(Beam([Span(L, EI, [Load("udl", w)])], left="fixed", right="pinned"))
    assert r["end_moments"][0][0] == approx(-w * L**2 / 8)
    assert r["reactions"] == [approx(5 * w * L / 8), approx(3 * w * L / 8)]


def test_simply_supported_point_load_centre():
    r = solve(Beam([Span(L, EI, [Load("point", P, a=L / 2)])]))
    assert r["max_sagging"] == approx(P * L / 4, 1e-2)
    v_max = max(r["deflection"][0]["v"])
    assert v_max == approx(P * L**3 / (48 * EI), 1e-3)


def test_two_equal_spans_udl():
    r = solve(Beam([Span(L, EI, [Load("udl", w)]), Span(L, EI, [Load("udl", w)])]))
    assert r["end_moments"][0][1] == approx(w * L**2 / 8)
    assert r["reactions"] == [approx(3 * w * L / 8), approx(10 * w * L / 8), approx(3 * w * L / 8)]


def test_fixed_fixed_settlement():
    d = 0.01
    r = solve(Beam([Span(L, EI)], left="fixed", right="fixed", settlements=[0, d]))
    m = 6 * EI * d / L**2
    assert list(r["end_moments"][0]) == [approx(-m), approx(-m)]


def test_point_load_off_centre_fem():
    a, b = 2.0, 4.0
    r = solve(Beam([Span(a + b, EI, [Load("point", P, a=a)])], left="fixed", right="fixed"))
    assert list(r["end_moments"][0]) == [approx(-P * a * b**2 / (a + b)**2), approx(P * a**2 * b / (a + b)**2)]


def test_triangular_load_fem():
    r = solve(Beam([Span(L, EI, [Load("linear", 0.0, a=0, b=L, w2=w)])], left="fixed", right="fixed"))
    assert list(r["end_moments"][0]) == [approx(-w * L**2 / 30, 1e-3), approx(w * L**2 / 20, 1e-3)]


def test_overhang_right():
    c = 2.0
    r = solve(Beam([Span(L, EI)], overhang_right=Span(c, EI, [Load("point", P, a=c)])))
    # simply supported span with a tip load on a right overhang: R_left = −P·c/L (uplift)
    assert r["reactions"][0] == approx(-P * c / L)
    assert r["reactions"][1] == approx(P + P * c / L)
    assert r["max_hogging"] == approx(-P * c)


def test_deflection_closes_at_supports():
    r = solve(Beam([Span(5, EI, [Load("udl", w)]), Span(7, EI * 2, [Load("point", P, a=3)])], left="fixed"))
    for seg in r["deflection"]:
        assert abs(seg["v"][-1]) < 1e-6
