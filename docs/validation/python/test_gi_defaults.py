"""Effect of the generation-interval default on Rt, and reference values for tests/test_optionsRt.m."""
import sys, os, json, warnings, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings('ignore')
import gp_port as gp
from rt_port import get_Rt_new
from scipy.stats import gamma
GI = {'old default / cholera: 5 d, SD 8 d': (5/7, (8/7)**2),
      'new default / mpox, US symptom onset: 8.5 d, SD 5.0 d': (8.5/7, (5.0/7)**2),
      'mpox, US rash onset: 7.0 d, SD 4.2 d': (7.0/7, (4.2/7)**2)}
out = {}
print("=== 1. Distributions (weekly units)")
for k, (m, v) in GI.items():
    b = v/m; a = m/b
    print(f"  {k:<54} mean {m:.4f}  variance {v:.4f}  shape {a:.4f}  scale {b:.4f}  P(within 1 day) {gamma.cdf(1/7, a, scale=b):.1%}  95% within {gamma.ppf(.95, a, scale=b):.2f} weeks")
    out[k] = dict(mean=m, var=v, shape=a, scale=b, p1day=float(gamma.cdf(1/7, a, scale=b)), q95=float(gamma.ppf(.95, a, scale=b)))

print("\n=== 2. Exponential growth, t = 0..19, incidence 100*exp(r t): estimate at the window end vs exact (1 + r b)^a")
for k, (m, v) in GI.items():
    b = v/m; a = m/b
    for r in (-0.2, 0.3):
        t = np.arange(20.0); Rs, ts = get_Rt_new(t, 100*np.exp(r*t), 1, m, v)
        print(f"  {k:<54} r = {r:+.1f}: Rs(end) = {Rs[-1]:.10f}  Rs(500) = {Rs[499]:.10f}  exact = {(1 + r*b)**a:.6f}  ({Rs[-1]/(1 + r*b)**a - 1:+.2%})")
        out[f"{k}|r={r}"] = dict(end=float(Rs[-1]), i500=float(Rs[499]), exact=float((1 + r*b)**a))

if len(sys.argv) > 2:
    print("\n=== 3. U.S. mpox series: Rt of the fitted GLM curve (least squares, fixI0 = 1)")
    d = np.loadtxt(os.path.join(sys.argv[2], 'Most_Recent_Timeseries_US-CDC.txt'))
    for a0 in (1, 8, 13):
        t = d[a0-1:a0+19, 0]; y = d[a0-1:a0+19, 1].astype(float)
        f = gp.fit_model(t, y, gp.GLM, 1, 0, 10, 'new', np.random.default_rng(1))
        for k, (m, v) in GI.items():
            Rs, ts = get_Rt_new(t, f.curve, 1, m, v); tR = ts[1:]
            cross = next((tR[i] for i in range(300, len(Rs)) if Rs[i] < 1), np.nan)      # first time below 1 after 3 weeks of history
            at = lambda x: Rs[int(round((x - t[0])/0.01)) - 1]
            print(f"  rows {a0:>2}-{a0+19:<2} {k:<54} Rt at start+5: {at(t[0]+5):.2f}  start+10: {at(t[0]+10):.2f}  end: {Rs[-1]:.3f}  " + (f"falls below 1 at t = {cross:.1f}" if Rs[300] >= 1 else "below 1 throughout"))
            out[f"rows{a0}|{k}"] = dict(p5=float(at(t[0]+5)), p10=float(at(t[0]+10)), end=float(Rs[-1]), cross=float(cross))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
