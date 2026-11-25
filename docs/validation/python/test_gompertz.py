"""Gompertz model: dependence on the time labels and on DT, parameter recovery, bootstrap. HEAD vs adopted."""
import sys, os, json, warnings, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings('ignore')
import gp_port as gp
OUT = {}
d = np.loadtxt(os.path.join(sys.argv[2], 'Most_Recent_Timeseries_US-CDC.txt')); yall = d[:, 1].astype(float)
RNG = lambda: np.random.default_rng(1)

print("=== 1. Same 20 observations (mpox rows 1-20), time labels shifted: least-squares Gompertz fit")
y = yall[:20]
for shift in (0, 1, 7, 12, 30, 100):
    t = np.arange(20.0) + shift; cells = []
    for sch in ('old', 'gfix'):
        f = gp.fit_model(t, y, gp.GOM, 1, 0, 10, sch, RNG())
        cells.append(f"{'HEAD   ' if sch == 'old' else 'adopted'}: r = {f.P[0]:<8.4g} a = {f.P[2]:<7.4g} SSE = {f.fval:<10.5g} r on bound: {str(bool(f.atLB[0] or f.atUB[0])):<5}")
        OUT[f"shift|{shift}|{sch}"] = dict(r=float(f.P[0]), a=float(f.P[2]), sse=float(f.fval), r_on=bool(f.atLB[0] or f.atUB[0]))
    print(f"  first label {shift:>3} | " + " | ".join(cells))

print("\n=== 2. Same observations, time step DT changed (labels 0..19; times are labels*DT)")
for DT in (1.0, 7.0, 1/7, 30.0):
    t = np.arange(20.0); cells = []
    for sch in ('old', 'gfix'):
        f = gp.fit_model(t, y, gp.GOM, 1, 0, 10, sch, RNG(), DT=DT)
        cells.append(f"{'HEAD   ' if sch == 'old' else 'adopted'}: r*DT = {f.P[0]*DT:<7.4g} a*DT = {f.P[2]*DT:<7.4g} SSE = {f.fval:<10.5g} on bound: {str(bool(f.atLB[0] or f.atUB[0] or f.atLB[2] or f.atUB[2])):<5}")
        OUT[f"dt|{DT:.4g}|{sch}"] = dict(rDT=float(f.P[0]*DT), aDT=float(f.P[2]*DT), sse=float(f.fval))
    print(f"  DT = {DT:<6.4g} | " + " | ".join(cells))

print("\n=== 3. Simulated Gompertz epidemics with Poisson noise, 60 replicates, Poisson MLE (true r = 0.8, a = 0.15, C0 = 20; 25 points; 7-step forecast)")
r0, a0, C0, n, h = 0.8, 0.15, 20.0, 25, 7
tt = np.arange(n + h, dtype=float); mu = gp.incidence([r0, 1, a0, 1, C0], tt, gp.GOM)
for first in (0, 10, 40):
    res = {s: [] for s in ('old', 'gfix')}
    for rep in range(60):
        rng = np.random.default_rng(5000 + rep); yy = rng.poisson(mu[:n]).astype(float); yy[0] = C0
        for sch in res:
            f = gp.fit_model(tt[:n] + first, yy, gp.GOM, 1, 1, 10, sch, np.random.default_rng(rep))
            fc = gp.forecast(f.P, tt[:n] + first, h, gp.GOM, elapsed=sch in gp.ELAPSED)
            r_start = f.P[0]*np.exp(-f.P[2]*first) if sch == 'old' else f.P[0]      # growth rate at the first observation
            res[sch].append((f.P[0], f.P[2], r_start, bool(f.atLB[0] or f.atUB[0]), np.mean(np.abs(fc[-h:] - mu[n:])), f.fval))
    for sch in res:
        A = np.array(res[sch], float)
        print(f"  first label {first:>2} {'HEAD   ' if sch == 'old' else 'adopted'}: reported r median {np.median(A[:,0]):<9.4g} | growth rate at first observation: median {np.median(A[:,2]):.3f} (true 0.8) | a median {np.median(A[:,1]):.3f} (true 0.15) | r on bound {100*A[:,3].mean():3.0f}% | 7-step MAE / mean {A[:,4].mean()/mu[n:].mean():.3f}")
        OUT[f"sim|{first}|{sch}"] = dict(r=float(np.median(A[:,0])), r_start=float(np.median(A[:,2])), a=float(np.median(A[:,1])), r_on=float(A[:,3].mean()), err=float(A[:,4].mean()/mu[n:].mean()))

print("\n=== 4. Bootstrap as in the runners (least squares, normal noise, B = 200), mpox rows 13-32 and 8-27")
for a_row in (13, 8):
    t = d[a_row-1:a_row+19, 0]; w = yall[a_row-1:a_row+19]; nn = len(w)
    fut = yall[a_row+19:a_row+23] if len(yall) >= a_row + 23 else None
    for sch in ('old', 'gfix'):
        f = gp.fit_model(t, w.copy(), gp.GOM, 1, 0, 10, sch, RNG()); sd = np.sqrt(f.fval/(nn - 3)); rng = np.random.default_rng(42); P = []
        for j in range(200):
            yb = f.curve + sd*rng.standard_normal(nn); yb[0] = f.curve[0]
            g = gp.fit_model(t, yb, gp.GOM, 1, 0, 1, sch, rng, params0=list(f.P) + [1, 1], box=(f.LB, f.UB)); P.append(g.P)
        P = np.array(P); q = lambda k: np.quantile(P[:, k], [.025, .975])
        msg = f"  rows {a_row}-{a_row+19} {'HEAD   ' if sch == 'old' else 'adopted'}: RMSE {np.sqrt(f.fval/nn):6.1f} | r = {f.P[0]:.4g} ({q(0)[0]:.4g}, {q(0)[1]:.4g}), refits at upper bound {np.mean(P[:,0] >= f.UB[0]*(1-1e-3)):.0%} | a = {f.P[2]:.4g} ({q(2)[0]:.4g}, {q(2)[1]:.4g})"
        if fut is not None:
            fc = gp.forecast(f.P, t, 4, gp.GOM, elapsed=sch in gp.ELAPSED); msg += f" | 4-step MAE {np.mean(np.abs(fc[-4:] - fut)):.0f}"
        print(msg)
        OUT[f"boot|{a_row}|{sch}"] = dict(rmse=float(np.sqrt(f.fval/nn)), r=float(f.P[0]), rCI=q(0).tolist(), a=float(f.P[2]), aCI=q(2).tolist(), r_ub=float(np.mean(P[:,0] >= f.UB[0]*(1-1e-3))))
json.dump(OUT, open(sys.argv[1], 'w'), indent=1)
