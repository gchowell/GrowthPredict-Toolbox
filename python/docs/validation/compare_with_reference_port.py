"""Compare the Python package with the repository's reference port of fit_model.m.

docs/validation/python/gp_port.py (repository root) is the independent
re-implementation used to validate the MATLAB toolbox: its own Dormand-Prince
integrator (numba) and scipy least_squares. The package uses LSODA and SLSQP.
This script fits the same windows with both and compares the objective value at
the optimum: three windows of 20 observations (first, middle, last) of every
bundled series with unit time spacing, the seven growth models, least squares
and Poisson maximum likelihood, fixI0 = 1.

Run from the `python` folder (needs numba; about 10 minutes):

    python docs/validation/compare_with_reference_port.py

The summary is printed and written to docs/validation/reference_port_comparison.txt.
"""
import os
import sys
import tempfile
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON_DIR = os.path.dirname(os.path.dirname(HERE))
REPO = os.path.dirname(PYTHON_DIR)
sys.path.insert(0, PYTHON_DIR)

import growthpredict as gp  # noqa: E402
from growthpredict.fitting import fit_model, initialParams  # noqa: E402

sys.path.insert(0, os.path.join(REPO, "docs", "validation", "python"))
# leave no __pycache__ beside the reference port (Python bytecode and numba's compiled cache)
sys.dont_write_bytecode = True
os.environ.setdefault("NUMBA_CACHE_DIR", tempfile.mkdtemp(prefix="numba-cache-"))

import gp_port as ref  # noqa: E402

NAMES = {-1: "EXP", 0: "GGM", 1: "GLM", 2: "GRM", 3: "LM", 4: "RICH", 5: "GOM"}
W = 20
TOL = 1e-4


def windows():
    inp = gp.example_input_dir()
    for f in sorted(x for x in os.listdir(inp) if x.endswith(".txt")):
        data = gp.load_timeseries(f)[0][:, :2]
        n = data.shape[0]
        if n < W or not np.allclose(np.diff(data[:, 0]), 1):
            continue
        for s in sorted({1, max(1, (n - W + 1) // 2), n - W + 1}):
            d = data[s - 1 : s - 1 + W].copy()
            if d[0, 1] == 0:
                d[0, 1] = 0.01  # as the runners do
            if np.sum(d[:, 1]) < 20:
                continue
            yield f, s, d


def main():
    rows = []
    for f, s, d in windows():
        y = d[:, 1]
        for flag in NAMES:
            for method in (0, 1):
                if method == 1 and (np.any(y < 0) or np.any(y != np.round(y))):
                    continue
                r = ref.fit_model(d[:, 0], y, flag, 1, method, 10, "gfix" if flag == 5 else "new", np.random.default_rng(1))
                gp.rng(1)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    m = fit_model(d, initialParams(y, flag), 1, 10, 1, flag, 0, None, method1=method)
                rows.append((f, s, flag, method, m.fval, r.fval, m.P[:5], r.P))
        print(f, s, "done", flush=True)

    nofit = [row for row in rows if not (np.isfinite(row[4]) and np.isfinite(row[5]))]
    rows = [row for row in rows if np.isfinite(row[4]) and np.isfinite(row[5])]
    rel = np.array([(a - b) / max(1.0, abs(b)) for (_, _, _, _, a, b, _, _) in rows])
    lines = [
        "Package (LSODA + SLSQP) against docs/validation/python/gp_port.py (Dormand-Prince + least_squares)",
        "series %d, windows %d, fits %d (%d more without a finite objective in one of the two)"
        % (len({r[0] for r in rows}), len({(r[0], r[1]) for r in rows}), len(rows), len(nofit)),
        "objective equal within %g (relative): %d" % (TOL, int(np.sum(np.abs(rel) <= TOL))),
        "package lower: %d, package higher: %d" % (int(np.sum(rel < -TOL)), int(np.sum(rel > TOL))),
        "median |relative difference| %.2g, 95th percentile %.2g" % (np.median(np.abs(rel)), np.percentile(np.abs(rel), 95)),
        "",
        "model method   n  equal  package lower  package higher",
    ]
    for flag in NAMES:
        for method in (0, 1):
            rr = np.array([x for x, row in zip(rel, rows) if row[2] == flag and row[3] == method])
            lines.append("%-5s %-6s %3d  %5d  %13d  %14d" % (
                NAMES[flag], "LSQ" if method == 0 else "Poisson", rr.size, int(np.sum(np.abs(rr) <= TOL)),
                int(np.sum(rr < -TOL)), int(np.sum(rr > TOL))))
    lines += ["", "fits where the objectives differ:"]
    for x, (f, s, flag, method, a, b, Pm, Pr) in sorted(zip(rel, rows), key=lambda t: -abs(t[0])):
        if abs(x) > TOL:
            lines.append("  %s, start row %d, %s, %s: package %.6g, reference %.6g (relative %.2g)" % (
                f, s, NAMES[flag], "LSQ" if method == 0 else "Poisson", a, b, x))
            lines.append("      package   [r p a K I0] = %s" % np.array2string(np.asarray(Pm), precision=4))
            lines.append("      reference [r p a K I0] = %s" % np.array2string(np.asarray(Pr), precision=4))
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(HERE, "reference_port_comparison.txt"), "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
