"""Parameter recovery on simulated Poisson data: old vs proposed r box (Poisson MLE, fixI0=1)."""
import sys, json, warnings, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__))); import gp_port as gp
from multiprocessing import Pool
warnings.filterwarnings('ignore')
SC = {  # name: (fit model, true r, p, K, C0, n calibration points, horizon)
 'A  GLM r=0.18 p=0.94 K=1e4  (classic, near-exponential)': (gp.GLM, 0.18, 0.94, 1e4, 1.0, 60, 7),
 'B  GLM r=40   p=0.40 K=5e5  (sub-exponential, large counts)': (gp.GLM, 40., 0.40, 5e5, 50., 30, 7),
 'C  GGM r=25   p=0.30        (polynomial growth)': (gp.GGM, 25., 0.30, 0., 20., 25, 7),
 'D  GLM r=3    p=0.70 K=2e4  (intermediate)': (gp.GLM, 3., 0.70, 2e4, 5., 40, 7),
}
def one(args):
    name, rep = args
    flag, r, p, K, C0, n, h = SC[name]
    t = np.arange(n + h, dtype=float)
    mu = gp.incidence([r, p, 1, K, C0], t, flag)
    rng = np.random.default_rng(1000 * (list(SC).index(name) + 1) + rep)
    y = rng.poisson(mu[:n]).astype(float); y[0] = max(round(C0), 1)
    out = dict(name=name, rep=rep, true_r=r, true_p=p, true_K=K)
    for sch in ('old', 'new'):
        f = gp.fit_model(t[:n], y, flag, 1, 1, 10, sch, np.random.default_rng(rep))
        fc = gp.forecast(f.P, t[:n], h, flag)
        out[sch] = dict(r=f.P[0], p=f.P[1], K=f.P[3], nll=f.fval, r_on=bool(f.atLB[0] or f.atUB[0]), p_on=bool(f.atLB[1] or f.atUB[1]),
                        widened=bool(f.widened), agree=f.agree, fmae=float(np.mean(np.abs(fc[-h:] - mu[n:]))), fmean=float(np.mean(mu[n:])))
    return out
if __name__ == '__main__':
    R = int(sys.argv[2]); jobs = [(nm, k) for nm in SC for k in range(R)]
    with Pool(2) as pool: res = pool.map(one, jobs, chunksize=4)
    json.dump(res, open(sys.argv[1], 'w'))
    for nm in SC:
        rr = [x for x in res if x['name'] == nm]; flag, r, p, K, *_ = SC[nm]
        print(f"\n{nm}   ({len(rr)} replicates)")
        for sch in ('old', 'new'):
            g = lambda k: np.array([x[sch][k] for x in rr], float)
            line = f"  {sch}: r median {np.median(g('r')):.4g} (true {r:g}) | p median {np.median(g('p')):.3f}, mean bias {np.mean(g('p'))-p:+.3f}, RMSE {np.sqrt(np.mean((g('p')-p)**2)):.3f}"
            if K: line += f" | K median {np.median(g('K')):.4g}"
            line += f" | r on bound {100*np.mean(g('r_on')):.0f}% | p on bound {100*np.mean(g('p_on')):.0f}% | 7-step MAE/mean {np.mean(g('fmae'))/np.mean(g('fmean')):.3f} | start agreement {np.mean(g('agree')):.2f}"
            print(line)
        d = np.array([x['old']['nll'] - x['new']['nll'] for x in rr])
        print(f"  NLL(old) - NLL(new): median {np.median(d):.3g}, min {d.min():.3g}, max {d.max():.3g}")
