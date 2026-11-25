import sys, time, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import gp_port as gp
from scipy.integrate import solve_ivp
d = np.loadtxt(sys.argv[1]); t, y = d[:, 0], d[:, 1]
# 1) integrator against scipy LSODA on awkward parameter sets
rng = np.random.default_rng(0); worst = 0
for flag in (0, 1, 2, 3, 4):
    for _ in range(200):
        r = 10**rng.uniform(-3, 3.5); p = rng.uniform(.01, 1) if flag in (0,1,2) else 1.0
        a = rng.uniform(.05, 5) if flag in (2,4) else 1.0; K = 10**rng.uniform(3, 6); C0 = 10**rng.uniform(0, 2.5)
        tv = np.arange(0., 25.)
        mine, ok = gp.solve_cum(flag, r, p, a, K, C0, tv, 1e-10, 1e-10, 400000)
        f = {0: lambda tt,C: r*max(C[0],0)**p, 1: lambda tt,C: r*max(C[0],0)**p*(1-C[0]/K), 2: lambda tt,C: r*max(C[0],0)**p*(1-(max(C[0],0)/K)**a),
             3: lambda tt,C: r*C[0]*(1-C[0]/K), 4: lambda tt,C: r*C[0]*(1-(max(C[0],0)/K)**a)}[flag]
        s = solve_ivp(f, (0, 24), [C0], t_eval=tv, method='LSODA', rtol=1e-11, atol=1e-11)
        if flag == 0 and (not ok or not s.success or s.y.max() > 1e60): continue
        if not ok: print('port solver failed', flag, r, p, a, K, C0); continue
        worst = max(worst, np.max(np.abs(mine - s.y[0])/np.maximum(np.abs(s.y[0]), 1)))
print(f"integrator: worst relative difference vs LSODA over 1000 random parameter sets = {worst:.2e}")
# 2) README quick start (rows 1-20, GLM, LS, fixI0=1): README figure shows r=1.7, p=0.85, K=2.9e4
t0 = time.time()
f = gp.fit_model(t[:20], y[:20], gp.GLM, 1, 0, 10, 'old', np.random.default_rng(1))
print(f"quick start (old box): r={f.P[0]:.4f} p={f.P[1]:.4f} K={f.P[3]:.1f} SSE={f.fval:.5g} agree={f.agree:.2f} box r=[{f.LB[0]:.4g},{f.UB[0]:.4g}]  ({time.time()-t0:.2f}s)")
for sch in ('new', 'i0'):
    f = gp.fit_model(t[:20], y[:20], gp.GLM, 1, 0, 10, sch, np.random.default_rng(1))
    print(f"quick start ({sch} box): r={f.P[0]:.4f} p={f.P[1]:.4f} K={f.P[3]:.1f} SSE={f.fval:.5g} agree={f.agree:.2f} box slot1=[{f.LB[0]:.4g},{f.UB[0]:.4g}]")
for a, b in ((13, 32), (8, 27)):
    for sch in ('old', 'new', 'i0'):
        f = gp.fit_model(t[a-1:b], y[a-1:b], gp.GLM, 1, 0, 10, sch, np.random.default_rng(1))
        print(f"rows {a}-{b} {sch:>3}: r={f.P[0]:.4g} p={f.P[1]:.3f} K={f.P[3]:.5g} RMSE={np.sqrt(f.fval/20):.1f} atBound(r/i0)={bool(f.atLB[0] or f.atUB[0])} widened={f.widened} agree={f.agree:.2f} valid={f.n_valid}/{f.n_starts}")
