"""Checks for three corrections to the plotting functions (after v1.1.1).

1. Figure 400 of plotForecast_GrowthModels.m: enough panels for every number of windows (forecastPanelGrid.m).
2. Rolling-parameter figure of the two forecast plotting functions: all four panels read the same arrays.
3. Last window start in the five plotting functions: the adjustment of the runners (adjustLastWindowStart.m).

The logic is checked with the ports in plot_port.py, the expected values of the two MATLAB test files are
checked against the ports, and the MATLAB sources are read to confirm the edits are in place. MATLAB itself
is not run here.

usage: python3 docs/validation/python/test_plot_windows.py            (from the repository root)
"""
import os
import re
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from plot_port import factor, panel_grid_old, panel_grid_new, adjust_last_window_start, window_mat_name

REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
CODE = os.path.join(REPO, 'forecasting_growthmodels code')
fails = 0


def check(ok, msg):
    global fails
    fails += (not ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {msg}")


def src(name):
    with open(os.path.join(CODE, name), encoding='utf-8') as fh:
        return fh.read().replace('\r\n', '\n')


def code_only(text):
    """MATLAB source without comment lines and trailing comments (good enough for these files)."""
    out = []
    for line in text.split('\n'):
        s = line.strip()
        if s.startswith('%'):
            continue
        out.append(re.sub(r"\s%[^']*$", '', line))
    return '\n'.join(out)


# ------------------------------------------------------------------------------------------------
print("=== 1. Panel of forecasts (figure 400): rows and columns for n windows")
N = 200
too_small = [n for n in range(1, N + 1) if panel_grid_old(n)[0] * panel_grid_old(n)[1] < n]
four_plus = [n for n in range(1, N + 1) if len(factor(n)) >= 4]
print(f"  former grid has fewer panels than windows for {len(too_small)} of the counts 1..{N}: {too_small[:16]} ...")
check(too_small == four_plus, "these are exactly the counts with four or more prime factors")
check(all(panel_grid_new(n)[0] * panel_grid_new(n)[1] == n for n in range(1, 1001)),
      "new grid: rows*cols equals the number of windows for every count 1..1000")
check(all(panel_grid_new(n) == panel_grid_old(n) for n in range(1, 1001) if len(factor(n)) <= 3),
      "new grid is the former grid for every count with at most three prime factors (1..1000)")
check(panel_grid_new(6) == (2, 3) and panel_grid_new(1) == (1, 1) and panel_grid_new(7) == (7, 1) and panel_grid_new(12) == (4, 3),
      "README example and other small counts unchanged: 6 -> 2x3, 1 -> 1x1, 7 -> 7x1, 12 -> 4x3")
ok = True
for n in [m for m in range(1, 1001) if len(factor(m)) >= 4]:
    r, c = panel_grid_new(n)
    best = max(d for d in range(1, n + 1) if n % d == 0 and d * d <= n)      # brute force
    ok &= (r == best and r <= c)
check(ok, "four or more prime factors: rows is the largest divisor not above the square root (brute force, 1..1000)")
print("  windows   former grid (panels)   new grid")
for n in [16, 24, 32, 36, 40, 48, 54, 56, 60, 64, 72, 80, 81, 96, 100]:
    ro, co = panel_grid_old(n)
    rn, cn = panel_grid_new(n)
    print(f"  {n:7d}   {ro}x{co} ({ro*co:2d})              {rn}x{cn}")

print("\n=== 2. The same calls as plotForecast_GrowthModels.m, with Matplotlib subplots")


def draw(n, grid):
    """subplot(rows,cols,k) for every window, then the bottom row and the first panel (axis labels)."""
    rows, cols = grid(n)
    fig = plt.figure()
    try:
        for k in range(1, n + 1):
            fig.add_subplot(rows, cols, k)
        for c in range(1, cols + 1):
            fig.add_subplot(rows, cols, (rows - 1) * cols + c)
        fig.add_subplot(rows, cols, 1)
        return len({(ax.get_subplotspec().rowspan.start, ax.get_subplotspec().colspan.start) for ax in fig.get_axes()})
    finally:
        plt.close(fig)


for n in (6, 12, 16, 24, 32):
    try:
        old = f"{draw(n, panel_grid_old)} panels"
    except ValueError as e:
        old = "stops: " + str(e).split('\n')[0][:60]
    new = draw(n, panel_grid_new)
    check(new == n, f"{n:2d} windows: former grid -> {old}; new grid -> {new} panels")

print("\n=== 3. Last window start: plotting functions against runners")
ok = True
cases = 0
for n in (5, 16, 32, 50, 60, 107):
    for w in (2, 7, 10, 20, 30, 40):
        for t in range(1, 121):
            new_t, adjusted = adjust_last_window_start(t, w, n)
            fits = [s for s in range(-200, t + 1) if s + w - 1 <= n]          # brute force: latest start <= t that fits
            ok &= (new_t == max(fits)) and (adjusted == (t + w - 1 > n))
            cases += 1
check(ok, f"adjusted value is the latest start not above the request whose window ends inside the series ({cases} cases, brute force)")
for (ts, te, w, n, label) in ((1, 1, 20, 32, 'README quick start'), (1, 6, 20, 32, 'README rolling windows'), (13, 13, 20, 32, 'README last 20 observations'),
                             (1, 80, 10, 107, 'options_fit_cholera')):
    check(adjust_last_window_start(te, w, n) == (te, False), f"not adjusted: {label} (windows {ts}-{te} of {w}, {n} rows)")
t, adj = adjust_last_window_start(20, 20, 32)
runner = window_mat_name('Most_Recent_Timeseries_US-CDC.txt', 1, 1, 0, 0, 1, t, 20, 0)
before = window_mat_name('Most_Recent_Timeseries_US-CDC.txt', 1, 1, 0, 0, 1, 20, 20, 0)
after = window_mat_name('Most_Recent_Timeseries_US-CDC.txt', 1, 1, 0, 0, 1, adjust_last_window_start(20, 20, 32)[0], 20, 0)
print("  Run_Fit_GrowthModels(1,20,20), U.S. mpox series (32 rows): the runner adjusts tend1 to", t)
print("    runner writes         ", os.path.basename(runner))
print("    plot function, before ", os.path.basename(before))
print("    plot function, after  ", os.path.basename(after))
check(before != runner and after == runner, "the plotting function now asks for the file the runner wrote")
check(adjust_last_window_start(37, 20, 32) == (13, True), "options_fit_mpox.m as shipped (row 37 of 32): adjusted to 13, below tstart1 = 37, so nothing to plot")
check(adjust_last_window_start(1, 40, 32) == (-7, True), "window longer than the series: result -7, nothing to plot")

print("\n=== 4. Expected values in the MATLAB test files agree with the ports")
t1 = src(os.path.join('tests', 'test_forecastPanelGrid.m'))
counts = [int(x) for x in re.search(r"counts = \[([^\]]*)\];", t1).group(1).split()]
expected = [tuple(int(v) for v in pair.split()) for pair in re.search(r"expected = \[([^\]]*)\];", t1).group(1).split(';')]
check(len(counts) == len(expected) == 15 and all(panel_grid_new(n) == e for n, e in zip(counts, expected)),
      f"test_forecastPanelGrid.m: the {len(counts)} layouts of testCountsThatUsedToFail")
check(all(panel_grid_old(n)[0] * panel_grid_old(n)[1] < n for n in counts), "  and the former grid was too small for each of them")
t2 = src(os.path.join('tests', 'test_adjustLastWindowStart.m'))
calls = re.findall(r"adjustLastWindowStart\((-?\d+), (\d+), (n|\d+)\);[^\n]*\n+verifyEqual\(testCase, tend1, (-?\d+)\);\nverify(True|False)\(testCase, adjusted\)", t2)
check(len(calls) >= 8 and all(adjust_last_window_start(int(a), int(w), 32 if n == 'n' else int(n)) == (int(e), flag == 'True') for a, w, n, e, flag in calls),
      f"test_adjustLastWindowStart.m: the {len(calls)} hard-coded cases")

print("\n=== 5. MATLAB sources")
helper = code_only(src('adjustLastWindowStart.m'))
check('if (tend1+windowsize1-1) > nObservations' in helper and 'tend1=nObservations-windowsize1+1;' in helper,
      "adjustLastWindowStart.m: condition and value")
for runner_file in ('Run_Fit_GrowthModels.m', 'Run_Forecasting_GrowthModels.m'):
    r = code_only(src(runner_file))
    check('if (tend1+windowsize1-1) > length(data(:,1))' in r and 'tend1= length(data(:,1))-windowsize1+1;' in r,
          f"{runner_file}: the same condition and value, with nObservations = length(data(:,1))")

OUTPUTS = {'plotFit_GrowthModels.m': ['param_doubling', 'seq_doublingtimes', 'fit_model1'],
           'plotForecast_GrowthModels.m': ['param_doubling', 'seq_doublingtimes', 'forecast_model1'],
           'plotFit_ReproductionNumber.m': ['param_R', 'param_doubling', 'seq_doublingtimes'],
           'plotForecast_ReproductionNumber.m': ['param_R', 'param_doubling', 'seq_doublingtimes'],
           'plotRt.m': ['Rts', 'params', 'performanceC', 'AICcs']}
CALL = '[tend1,tend1Adjusted]=adjustLastWindowStart(tend1,windowsize1,length(data(:,1)));'
for name, outs in OUTPUTS.items():
    c = code_only(src(name))
    declared = [o.strip() for o in re.match(r"function \[([^\]]*)\]=", c).group(1).split(',')]
    pos = c.find(CALL)
    first_use = min(p for p in (c.find('tstart1:1:tend1'), c.find('tstart1:tend1'), c.find("num2str(tend1)")) if p >= 0)
    guard = re.search(r"if tstart1>tend1\n(.*?)\n\s*return\nend", c[pos:], re.S)
    ok = (c.count(CALL) == 1 and declared == outs
          and c.find('data = load(') < pos and c.find('windowsize1=windowsize1_INP;') < pos < first_use
          and guard is not None and all(re.search(r"^\s*%s=\[\];" % o, guard.group(1), re.M) for o in outs))
    check(ok, f"{name}: one call, after the data and the window length are known and before tend1 is used; "
              f"with no window to plot it returns empty {', '.join(outs)}")

for name in ('plotForecast_GrowthModels.m', 'plotForecast_ReproductionNumber.m'):
    c = code_only(src(name))
    block = c[c.rfind('if tend1>tstart1'):]
    block = block[:block.find("\nend\n")]
    used = re.findall(r"\b(param_[A-Za-z0-9]+s2?)\(", block)
    check(len(used) == 16 and all(u.endswith('2') for u in used) and sorted(set(used)) == ['param_Ks2', 'param_as2', 'param_ps2', 'param_rs2'],
          f"{name}: rolling-parameter figure reads {sorted(set(used))} only ({len(used)} uses)")

c = code_only(src('plotForecast_GrowthModels.m'))
loop = c[c.find('for i=tstart1:1:tend1'):c.find("if length(tstart1:1:tend1)>1\n\n    figure(400)\n    for c=1:cols")]
check('[rows,cols]=forecastPanelGrid(length(tstart1:1:tend1));' in c and 'factor(' not in c, "plotForecast_GrowthModels.m: figure 400 takes its grid from forecastPanelGrid.m")
check(len(loop) > 1000 and 'close all' not in loop and "~isequal(get(openFigures(kFigure),'Number'),400)" in loop,
      "plotForecast_GrowthModels.m: inside the window loop figure 400 is no longer closed")
g = code_only(src('forecastPanelGrid.m'))
check(all(s in g for s in ('factors=factor(nWindows);', 'rows=factors(1)*factors(2);', 'if d*d>nWindows', 'if mod(nWindows,d)==0', 'cols=nWindows/rows;')),
      "forecastPanelGrid.m: the statements ported in plot_port.py")

print(f"\n{'ALL CHECKS PASSED' if fails == 0 else str(fails) + ' CHECK(S) FAILED'}")
sys.exit(1 if fails else 0)
