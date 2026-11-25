import sys, os, json, warnings, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rt_port import get_Rt_old, get_Rt_new, DT
from scipy.stats import gamma
warnings.filterwarnings('ignore')
OUT = {}
def true_R(r, gi, mean1, var1):
    if gi == 1: b = var1/mean1; a = mean1/b; return (1 + r*b)**a
    if gi == 2: return 1 + r*mean1
    return np.exp(r*mean1)
GIS = {'default gamma (mean 5/7, var (8/7)^2)': (1, 5/7, (8/7)**2), 'gamma mean 1.4, sd 0.6': (1, 1.4, 0.36),
       'exponential mean 1': (2, 1.0, 0), 'fixed lag 1': (3, 1.0, 0)}

print("=== 1. Exponential growth, weekly samples, window of 20: estimate at the window end vs exact R")
print(f"{'generation interval':<40}{'r':>6}{'exact':>9} | " + " | ".join(f"start {t0:>2}: old / new" for t0 in (0, 1, 12)))
rows = []
for name, (gi, m, v) in GIS.items():
    for r in (-0.2, 0.1, 0.3):
        Rt = true_R(r, gi, m, v); cells = []
        for t0 in (0, 1, 12):
            t = np.arange(t0, t0 + 20.0); I = 100*np.exp(r*(t - t0))
            try: o = get_Rt_old(t, I, gi, m, v)[0][-1]
            except Exception as e: o = np.nan
            nw = get_Rt_new(t, I, gi, m, v)[0][-1]
            cells.append((o, nw)); rows.append(dict(gi=name, r=r, exact=Rt, t0=t0, old=o, new=nw))
        print(f"{name:<40}{r:>6}{Rt:>9.4f} | " + " | ".join(f"{o:>13.5g} / {n:<8.4f}" for o, n in cells))
OUT['exact'] = rows
err_new = max(abs(x['new']/x['exact'] - 1) for x in rows)
err_old0 = max(abs(x['old']/x['exact'] - 1) for x in rows if x['t0'] == 0)
err_old = [abs(x['old']/x['exact'] - 1) for x in rows if x['t0'] > 0 and np.isfinite(x['old'])]
print(f"largest relative error, new (all starts): {err_new:.3%} | old at start 0: {err_old0:.3%} | old at start 1 or 12: median {np.median(err_old):.0%}, max {max(err_old):.3g}x")

print("\n=== 2. Shift invariance: same curve, window start moved")
t = np.arange(0, 20.0); I = 50*np.exp(0.25*t) / (1 + np.exp(0.5*(t - 10)))       # rise and fall
for name, (gi, m, v) in GIS.items():
    base = get_Rt_new(t, I, gi, m, v)[0]; baseo = get_Rt_old(t, I, gi, m, v)[0]
    dn = max(np.nanmax(np.abs(get_Rt_new(t + s, I, gi, m, v)[0] / base - 1)) for s in (1, 5, 12, 50, 1000))
    do = []
    for s in (1, 5, 12):
        x = get_Rt_old(t + s, I, gi, m, v)[0]; ok = np.isfinite(x) & np.isfinite(baseo)
        do.append(np.max(np.abs(x[ok]/baseo[ok] - 1)) if ok.any() else np.nan)
    print(f"  {name:<40} new: max relative change = {dn:.1e} | old: max relative change = {np.nanmax(do):.3g}")
OUT.setdefault('shift_new_max', []).append(float(dn))

print("\n=== 3. Fixed lag: every lag k/100, k = 1..1000 (dense exponential input, r = 0.3, window 0..30)")
tf = DT*np.arange(3001); If = np.exp(0.3*tf); bad_old = []; lag_err_old = []; lag_err_new = []; fail_new = 0
for k in range(1, 1001):
    g = k/100
    try:
        Ro = get_Rt_old(tf, If, 3, g)[0][-1]; lag_err_old.append(np.log(Ro)/0.3 - g)
    except IndexError: bad_old.append(g)
    try:
        Rn = get_Rt_new(tf, If, 3, g)[0][-1]; lag_err_new.append(np.log(Rn)/0.3 - g)
    except Exception: fail_new += 1
print(f"  old: {len(bad_old)} lags raise an indexing error (e.g. {bad_old[:6]}); implied lag error where it runs: mean {np.mean(lag_err_old):+.4f}")
print(f"  new: {fail_new} failures; implied lag error: max |{np.max(np.abs(lag_err_new)):.1e}|")
for g in (5/7, 0.714, 2.5/7):
    Rn = get_Rt_new(tf, If, 3, g)[0][-1]; print(f"  new, lag {g:.4f} (not on the grid): implied lag {np.log(Rn)/0.3:.4f}")
OUT['lags'] = dict(old_errors=len(bad_old), old_lag_bias=float(np.mean(lag_err_old)), new_fail=fail_new, new_max=float(np.max(np.abs(lag_err_new))))
Rn, tsn = get_Rt_new(np.arange(0, 20.0), np.exp(0.3*np.arange(20.0)), 3, 25.0)
print(f"  new, lag longer than the window: all NaN = {bool(np.all(np.isnan(Rn)))}")
try: get_Rt_new(np.arange(0, 20.0), np.ones(20), 3, 0.004); print("  new, lag 0.004: no error (unexpected)")
except ValueError as e: print("  new, lag 0.004: rejected")

print("\n=== 4. Renewal simulation with a known, changing R(t) (gamma GI mean 1.4, sd 0.6; 40 time units)")
m, v = 1.4, 0.36; b = v/m; a = m/b
T = 40; N = int(T/DT) + 1; tt = DT*np.arange(N)
Rtrue = 0.7 + 1.9/(1 + np.exp(0.45*(tt - 14)))                     # 2.6 falling to 0.7
w = np.diff(gamma.cdf(DT*np.arange(N), a, scale=b))
I = np.zeros(N); seed = int(3/DT); I[:seed + 1] = 5*np.exp(0.5*tt[:seed + 1])
for i in range(seed + 1, N):
    I[i] = Rtrue[i]*np.dot(I[i - 1::-1][:i], w[:i])
ti = np.arange(0, T + 1.0); Ii = I[::100]
print(f"  simulated incidence: peak {Ii.max():.0f} at t = {ti[Ii.argmax()]:.0f}; R(t) from {Rtrue[seed]:.2f} to {Rtrue[-1]:.2f}")
res4 = []
for (a0, b0) in ((0, 20), (5, 25), (12, 32), (20, 40)):
    sel = (ti >= a0) & (ti <= b0 - 1)
    for nm, fn in (('old', get_Rt_old), ('new', get_Rt_new)):
        Rs, ts = fn(ti[sel], Ii[sel], 1, m, v)
        tR = ts[1:]; truth = np.interp(tR, tt, Rtrue)
        late = tR >= a0 + 6                                           # after 6 time units of history
        res4.append(dict(win=f"{a0}-{b0-1}", est=nm, end=float(Rs[-1]), truth_end=float(truth[-1]),
                         max_rel_late=float(np.nanmax(np.abs(Rs[late]/truth[late] - 1)))))
    o, n_ = res4[-2], res4[-1]
    print(f"  window {o['win']:>5}: true R at end {o['truth_end']:.3f} | old {o['end']:.4g} | new {n_['end']:.3f} | max relative error after 6 units of history: old {o['max_rel_late']:.3g}, new {n_['max_rel_late']:.3%}")
OUT['renewal'] = res4

print("\n=== 5. How much history the estimate needs (exponential growth r = 0.2, new function, relative error vs elapsed time)")
for name, (gi, m_, v_) in list(GIS.items())[:3]:
    t = np.arange(0, 30.0); Rs, ts = get_Rt_new(t, np.exp(0.2*t), gi, m_, v_); rel = Rs/true_R(0.2, gi, m_, v_) - 1; tR = ts[1:]
    within = lambda tol: next((tR[k] for k in range(len(tR)) if np.all(np.abs(rel[k:]) < tol)), np.nan)
    print(f"  {name:<40} error at elapsed 1, 2, 4, 8: " + ", ".join(f"{rel[np.searchsorted(tR, x)]:+.1%}" for x in (1, 2, 4, 8)) + f" | within 5% from {within(.05):.2f}, within 1% from {within(.01):.2f}")

print("\n=== 6. CSV export labels (window 0..19): row pairing of time label and estimate")
t = np.arange(0, 20.0); Rs, ts = get_Rt_new(t, 100*np.exp(0.1*t), 1, 5/7, (8/7)**2)
old_rows = [(ts[k], ts[k+1]) for k in range(1, len(Rs), 100)][:3]                     # MATLAB ts(2:100:end) with Rs(2:100:end)
idx = list(range(100, len(Rs) + 1, 100))
if not idx or idx[-1] != len(Rs): idx.append(len(Rs))
new_rows = [(ts[k], ts[k]) for k in idx][:3]                                          # MATLAB ts(idx+1) with Rs(idx)
print("  old (label, time the estimate refers to):", [(round(a,2), round(b,2)) for a, b in old_rows], f"... {len(range(1, len(Rs), 100)) + 1} rows")
print("  new (label, time the estimate refers to):", [(round(a,2), round(b,2)) for a, b in new_rows], f"... {len(idx)} rows, last label {ts[idx[-1]]:.2f}")

if len(sys.argv) > 2:
    print("\n=== 7. Bundled series: Rt of the fitted GLM curve at the window end (default GI), old vs new")
    import glob, gp_port as gp
    rows7 = []
    def one(path, a0, w):
        d = np.loadtxt(path); t = d[a0-1:a0-1+w, 0]; y = d[a0-1:a0-1+w, 1].astype(float).copy()
        if y[0] == 0: y[0] = 0.01
        f = gp.fit_model(t, y, gp.GLM, 1, 0, 10, 'new', np.random.default_rng(1))
        try: o = get_Rt_old(t, f.curve, 1, 5/7, (8/7)**2)[0][-1]
        except Exception: o = np.nan
        n_ = get_Rt_new(t, f.curve, 1, 5/7, (8/7)**2)[0][-1]
        rows7.append(dict(series=os.path.basename(path)[:-4], rows=f"{a0}-{a0+w-1}", t_first=float(t[0]), old=float(o), new=float(n_)))
        print(f"  {os.path.basename(path)[:44]:<44} rows {a0:>3}-{a0+w-1:<3} first time label {t[0]:>4.0f} | old {o:>11.5g} | new {n_:.3f}")
    mp = os.path.join(sys.argv[2], 'Most_Recent_Timeseries_US-CDC.txt')
    for a0 in (1, 4, 8, 13): one(mp, a0, 20)
    for path in sorted(glob.glob(os.path.join(sys.argv[2], '*.txt'))):
        d = np.loadtxt(path)
        if d[0, 0] != 0 and np.allclose(np.diff(d[:, 0]), 1) and len(d) >= 12: one(path, 1, min(20, len(d)))
    OUT['bundled'] = rows7
json.dump(OUT, open(sys.argv[1], 'w'), indent=1)
