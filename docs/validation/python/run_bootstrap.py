"""Parametric bootstrap as in the runners (LS fit, normal noise, first point fixed), old vs proposed box."""
import sys, json, warnings, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__))); import gp_port as gp
warnings.filterwarnings('ignore')
d = np.loadtxt(sys.argv[1]); B = int(sys.argv[2]); res = []
for (a, b) in ((1, 20), (8, 27), (13, 32)):
    t, y = d[a-1:b, 0], d[a-1:b, 1].astype(float)
    print(f"\nrows {a}-{b}")
    for sch in ('old', 'new'):
        f = gp.fit_model(t, y, gp.GLM, 1, 0, 10, sch, np.random.default_rng(1))
        n = len(y); sd = np.sqrt(f.fval / (n - 4))               # runner: SSE/(n - numparams), numparams = 3 + 1
        rng = np.random.default_rng(42); P = np.empty((B, 5)); fc = np.empty((B, 4))
        for j in range(B):
            yb = f.curve + sd * rng.standard_normal(n); yb[0] = f.curve[0]
            g = gp.fit_model(t, yb, gp.GLM, 1, 0, 1, sch, rng, params0=list(f.P) + [1, 1], box=(f.LB, f.UB))
            P[j] = g.P; fc[j] = gp.forecast(g.P, t, 4, gp.GLM)[-4:]
        q = lambda k: np.quantile(P[:, k], [.025, .975])
        hit = lambda k: (np.mean(P[:, k] <= f.LB[k] * (1 + 1e-3)), np.mean(P[:, k] >= f.UB[k] * (1 - 1e-3)))
        fut = d[b:b+4, 1] if len(d) >= b + 4 else None
        fcd = gp.forecast(f.P, t, 4, gp.GLM)[-4:]
        row = dict(rows=f"{a}-{b}", scheme=sch, r=f.P[0], p=f.P[1], K=f.P[3], rCI=q(0).tolist(), pCI=q(1).tolist(), KCI=q(3).tolist(),
                   r_hit=hit(0), p_hit=hit(1), K_hit=hit(3), rmse=float(np.sqrt(f.fval/n)), sd=float(sd))
        s = (f"  {sch}: RMSE {row['rmse']:.1f} | r={f.P[0]:.4g} [{q(0)[0]:.4g}, {q(0)[1]:.4g}] at LB/UB {hit(0)[0]:.0%}/{hit(0)[1]:.0%}"
             f" | p={f.P[1]:.3f} [{q(1)[0]:.3f}, {q(1)[1]:.3f}] at LB/UB {hit(1)[0]:.0%}/{hit(1)[1]:.0%}"
             f" | K={f.P[3]:.5g} [{q(3)[0]:.5g}, {q(3)[1]:.5g}] at LB/UB {hit(3)[0]:.0%}/{hit(3)[1]:.0%}")
        if fut is not None:
            lo, hi = np.quantile(fc + sd * rng.standard_normal(fc.shape), [.025, .975], axis=0); lo = np.maximum(lo, 0)
            row.update(mae4=float(np.mean(np.abs(fcd - fut))), cover=float(np.mean((fut >= lo) & (fut <= hi))))
            s += f" | 4-step MAE {row['mae4']:.0f}, 95% PI covers {int(round(4*row['cover']))}/4"
        print(s); res.append(row)
json.dump(res, open(sys.argv[3], 'w'))
