"""Tests of the coverage filter for Rt rows (selectRtRows.m) with the Python ports."""
import sys, os, json, warnings, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings('ignore')
from rt_port import get_Rt_new, gi_coverage, select_rt_rows
from scipy.stats import gamma
OUT = {}
GIS = {'mpox default: 8.5 d, SD 5.0 d': (1, 8.5/7, (5/7)**2), 'cholera: 5 d, SD 8 d': (1, 5/7, (8/7)**2),
       'mpox rash onset: 7.0 d, SD 4.2 d': (1, 1.0, 0.36), 'exponential, mean 1': (2, 1.0, 0), 'fixed, 1': (3, 1.0, 0)}
def true_R(r, gi, m, v):
    if gi == 1: b = v/m; return (1 + r*b)**(m/b)
    return 1 + r*m if gi == 2 else np.exp(r*m)

print("=== 1. Which rows are kept (window labelled 12..31, 95% coverage)")
for name, (gi, m, v) in GIS.items():
    t = np.arange(12, 32.0); Rs, ts = get_Rt_new(t, np.exp(-0.2*(t - 12)), gi, m, v)
    cov = gi_coverage(len(Rs), gi, m, v); idx, rel = select_rt_rows(cov)
    first = int(np.argmax(rel)) + 1
    print(f"  {name:<34} reliable from t = {ts[first]:.2f} | CSV rows: t = {ts[idx[0]]:.0f}..{ts[idx[-1]]:.0f} ({len(idx)} rows; dropped t = {[int(ts[k]) for k in range(100, len(Rs)+1, 100) if k not in idx]})")
    OUT[f"rows|{name}"] = dict(first_reliable=float(ts[first]), first_index=first, n_rows=len(idx), first_row=float(ts[idx[0]]))

print("\n=== 2. Exponential growth or decline, exact R known: relative error of the rows dropped and kept (unit-time rows)")
print(f"  {'generation interval':<34}{'r':>6} | dropped rows: errors | kept rows: largest |error|")
worst_kept = {}
for name, (gi, m, v) in GIS.items():
    for r in (-0.3, -0.2, 0.1, 0.3):
        t = np.arange(12, 32.0); Rs, ts = get_Rt_new(t, 100*np.exp(r*(t - 12)), gi, m, v)
        cov = gi_coverage(len(Rs), gi, m, v); idx, rel = select_rt_rows(cov); ex = true_R(r, gi, m, v)
        unit = np.arange(100, len(Rs) + 1, 100); dropped = [k for k in unit if k not in idx]
        ed = [Rs[k-1]/ex - 1 for k in dropped]; ek = max(abs(Rs[k-1]/ex - 1) for k in idx)
        worst_kept[name] = max(worst_kept.get(name, 0), ek)
        print(f"  {name:<34}{r:>6} | " + (", ".join(f"{e:+.0%}" if np.isfinite(e) else "none" for e in ed) or "-") + f" | {ek:.1%}")
        OUT[f"exp|{name}|{r}"] = dict(dropped=[float(e) for e in ed], kept_max=float(ek))
print("  largest |error| among kept rows, per interval:", {k: f"{v:.1%}" for k, v in worst_kept.items()})

print("\n=== 3. Choice of threshold: first exported row and largest |error| among kept unit-time rows (window 12..31)")
for gname in ('mpox default: 8.5 d, SD 5.0 d', 'cholera: 5 d, SD 8 d', 'exponential, mean 1'):
    gi, m, v = GIS[gname]
    for thr in (0.90, 0.95, 0.99):
        cells = []
        for r in (-0.3, -0.2, 0.1, 0.3):
            t = np.arange(12, 32.0)
            Rs, ts = get_Rt_new(t, 100*np.exp(r*(t - 12)), gi, m, v); idx, _ = select_rt_rows(gi_coverage(len(Rs), gi, m, v), thr)
            cells.append((r, max(abs(Rs[k-1]/true_R(r, gi, m, v) - 1) for k in idx), ts[idx[0]]))
        print(f"  {gname:<32} threshold {thr:.2f}: first row t = {cells[0][2]:.0f} | " + "  ".join(f"r={r:+.1f}: {e:.1%}" for r, e, _ in cells))
        OUT[f"thr|{gname}|{thr}"] = [dict(r=r, err=float(e), first=float(f)) for r, e, f in cells]

print("\n=== 4. Window too short for the generation interval (mean 10, SD 5, window 0..5)")
t = np.arange(0, 6.0); Rs, ts = get_Rt_new(t, np.exp(0.2*t), 1, 10.0, 25.0); cov = gi_coverage(len(Rs), 1, 10.0, 25.0); idx, rel = select_rt_rows(cov)
print(f"  coverage at the end {cov[-1]:.3f} | any reliable: {bool(rel.any())} | rows exported: t = {[float(ts[k]) for k in idx]} (the last time, flagged by its coverage)")

if len(sys.argv) > 2:
    import gp_port as gp
    print("\n=== 5. U.S. mpox series, Rt of the fitted GLM curve with the mpox default: unit-time rows before and after the filter")
    d = np.loadtxt(os.path.join(sys.argv[2], 'Most_Recent_Timeseries_US-CDC.txt')); gi, m, v = GIS['mpox default: 8.5 d, SD 5.0 d']
    for a0 in (1, 8, 13):
        t = d[a0-1:a0+19, 0]; y = d[a0-1:a0+19, 1].astype(float)
        f = gp.fit_model(t, y, gp.GLM, 1, 0, 10, 'new', np.random.default_rng(1))
        Rs, ts = get_Rt_new(t, f.curve, gi, m, v); cov = gi_coverage(len(Rs), gi, m, v); idx, rel = select_rt_rows(cov)
        unit = np.arange(100, len(Rs) + 1, 100); dropped = [k for k in unit if k not in idx]
        print(f"  rows {a0:>2}-{a0+19:<2}: dropped " + ", ".join(f"t={ts[k]:.0f}: Rt {Rs[k-1]:.2f} (coverage {cov[k-1]:.0%})" for k in dropped) + f" | kept from t={ts[idx[0]]:.0f}: Rt {Rs[idx[0]-1]:.2f} (coverage {cov[idx[0]-1]:.0%}) ... t={ts[idx[-1]]:.0f}: Rt {Rs[-1]:.3f} | {len(idx)} rows | largest kept Rt {max(Rs[k-1] for k in idx):.2f}")
        OUT[f"mpox|{a0}"] = dict(dropped=[(float(ts[k]), float(Rs[k-1]), float(cov[k-1])) for k in dropped], first=(float(ts[idx[0]]), float(Rs[idx[0]-1])), end=float(Rs[-1]), n=len(idx))
json.dump(OUT, open(sys.argv[1], 'w'), indent=1)
