# What changed from version 1

Version 1 (kept in [`legacy/`](../legacy/)) was a Flask app that solved the slope-deflection equations symbolically with SymPy. Writing textbook tests for version 2 turned up three bugs in it.

| # | Bug | Effect | Textbook check |
|---|---|---|---|
| 1 | The point-load position `span_a_value` was read from the input but overwritten with `1` for every span. | Every off-centre point load was treated as sitting 1 m from the left support. | Fixed-fixed 6 m span, 40 kN at 2 m: v1 gave 27.78 / 5.56 kN·m, correct is −35.56 / +17.78. |
| 2 | The left fixed-end moment for an off-centre point load had the wrong sign (every other load case was negative). | Combined with bug 1, the moments above are wrong in both size and direction. | Same case as above. |
| 3 | At a pinned end, the modified slope-deflection equation used only the near-end FEM and dropped the carry-over term (− FEM_far ⁄ 2). | Interior support moments came out too small whenever an end span was loaded. | Two equal 6 m spans, 10 kN/m: v1 gave M_B = 30 kN·m, correct is wL²/8 = 45 kN·m. |

Version 2 avoids the modified equation entirely: a pinned end simply keeps its rotation as an unknown with its own equilibrium equation (M = 0). The system is solved numerically (NumPy / Gaussian elimination) instead of symbolically, which is faster and handles any number of spans.

Other changes:

- The committed `env/` virtual environment (thousands of third-party files) and `__pycache__/` folders were removed and git-ignored.
- `EI` was hard-coded to 1, so support settlements gave meaningless moments. Spans now take E and I.
- Cantilever overhangs, partial UDLs and triangular loads were added; loads can be combined on one span.
- Shear, moment and deflection diagrams are drawn live in the browser. No server is needed.
- Nine Python tests plus a JS parity suite run on every push (GitHub Actions).
