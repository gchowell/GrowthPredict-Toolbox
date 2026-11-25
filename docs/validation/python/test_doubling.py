"""Doubling times of cumulative incidence: window-relative count (old) vs epidemic count (new)."""
import sys, os, json, warnings, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings('ignore')
from dt_port import doubling_times
OUT = {}
f = lambda a: "[" + ", ".join(f"{v:.2f}" for v in a) + "]" if len(a) else "none"

print("=== 1. Exponential growth, incidence exp(r k) for k = 0..39; windows of 20. Doubling time of the incidence is ln(2)/r.")
for r in (0.1, 0.2, 0.35):
    k = np.arange(40); I = np.exp(r*k)
    print(f"  r = {r}: ln(2)/r = {np.log(2)/r:.2f}")
    for a0 in (0, 10, 20):
        w = I[a0:a0+20]; prior = I[:a0].sum()
        old, _ = doubling_times(w); new, _ = doubling_times(w, 1.0, prior)
        # reference: the cumulative count of the whole series, followed from the window start to the window end
        ref, _ = doubling_times(np.r_[I[:a0+1].sum(), I[a0+1:a0+20]])
        print(f"    window starts at {a0:>2}: window-relative {f(old)} | epidemic count {f(new)} | max |new - whole series| = {np.max(np.abs(new-ref)) if len(new)==len(ref) and len(new) else 0:.1e}")
        OUT[f"exp|{r}|{a0}"] = dict(old=old.tolist(), new=new.tolist(), ref=ref.tolist())

print("\n=== 2. Epidemic that rises and falls (logistic cumulative, K = 20000, r = 0.5, midpoint 12); windows of 20")
k = np.arange(41); C = 20000/(1 + np.exp(-0.5*(k - 12))); I = np.r_[C[0], np.diff(C)]
for a0 in (0, 6, 12, 16, 20):
    w = I[a0:a0+20]; prior = I[:a0].sum()
    old, _ = doubling_times(w); new, _ = doubling_times(w, 1.0, prior)
    trend = 'rising' if w[3] > w[1] else 'falling'
    print(f"  window starts at {a0:>2} (incidence {trend}; {prior + w[0]:.0f} cases so far, {C[min(a0+19, 40)]:.0f} by its end): window-relative {f(old)} | epidemic count {f(new)}")
    OUT[f"logistic|{a0}"] = dict(old=old.tolist(), new=new.tolist())

if len(sys.argv) > 2:
    import gp_port as gp
    d = np.loadtxt(os.path.join(sys.argv[2], 'Most_Recent_Timeseries_US-CDC.txt')); y = d[:, 1]
    print("\n=== 3. U.S. mpox series, windows of 20: observed counts and the fitted GLM curve")
    for a0 in (1, 4, 8, 13):
        w = y[a0-1:a0+19]; prior = y[:a0-1].sum()
        fit = gp.fit_model(d[a0-1:a0+19, 0], w.copy(), gp.GLM, 1, 0, 10, 'new', np.random.default_rng(1)).curve
        o_obs, _ = doubling_times(w); n_obs, _ = doubling_times(w, 1.0, prior)
        o_fit, _ = doubling_times(np.maximum(fit, 0)); n_fit, _ = doubling_times(np.maximum(fit, 0), 1.0, prior)
        print(f"  rows {a0:>2}-{a0+19:<2} (cases before the window {prior:>6.0f}; cumulative {prior+w[0]:>6.0f} -> {prior+w.sum():>6.0f})")
        print(f"      observed counts : window-relative {f(o_obs)} | epidemic count {f(n_obs)}")
        print(f"      fitted GLM curve: window-relative {f(o_fit)} | epidemic count {f(n_fit)}")
        OUT[f"mpox|{a0}"] = dict(prior=float(prior), old_obs=o_obs.tolist(), new_obs=n_obs.tolist(), old_fit=o_fit.tolist(), new_fit=n_fit.tolist())
    print("\n=== 4. As the plotters report it: bootstrap of the GLM fit (least squares, normal noise, B = 200), mpox rows 13-32 and 8-27")
    for a0 in (13, 8):
        t = d[a0-1:a0+19, 0]; w = y[a0-1:a0+19].astype(float); prior = y[:a0-1].sum(); n = len(w)
        fobs = gp.fit_model(t, w.copy(), gp.GLM, 1, 0, 10, 'new', np.random.default_rng(1))
        sd = np.sqrt(fobs.fval/(n - 4)); rng = np.random.default_rng(42); seqs = {'old': [], 'new': []}
        for j in range(200):
            yb = fobs.curve + sd*rng.standard_normal(n); yb[0] = fobs.curve[0]
            g = gp.fit_model(t, yb, gp.GLM, 1, 0, 1, 'new', rng, params0=list(fobs.P) + [1, 1], box=(fobs.LB, fobs.UB))
            c = np.maximum(g.curve, 0)
            seqs['old'].append(doubling_times(c)[0]); seqs['new'].append(doubling_times(c, 1.0, prior)[0])
        for nm in ('old', 'new'):
            S = seqs[nm]; maxd = max(1, max(len(x) for x in S)); rows = []
            for j in range(maxd):
                v = np.array([x[j] for x in S if len(x) > j])
                rows.append((j + 1, np.mean(v) if len(v) else np.nan, np.quantile(v, .025) if len(v) else np.nan, np.quantile(v, .975) if len(v) else np.nan, len(v)/len(S)))
            means = np.array([np.mean(x) if len(x) else np.nan for x in S])
            lab = 'window-relative' if nm == 'old' else 'epidemic count '
            print(f"  rows {a0}-{a0+19} {lab}: " + "; ".join(f"doubling {r[0]}: {r[1]:.2f} ({r[2]:.2f}, {r[3]:.2f}), reached by {r[4]:.0%}" for r in rows) + f" | mean per trajectory {np.nanmean(means):.2f}")
            OUT[f"boot|{a0}|{nm}"] = dict(rows=[[float(z) for z in r] for r in rows], mean=float(np.nanmean(means)))
json.dump(OUT, open(sys.argv[1], 'w'), indent=1)
