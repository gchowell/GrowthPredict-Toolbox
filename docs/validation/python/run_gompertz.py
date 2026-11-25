"""Gompertz fits on every bundled series: HEAD vs elapsed time vs elapsed time with data-independent box."""
import sys, os, glob, json, zlib, warnings, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gp_port as gp
from multiprocessing import Pool
warnings.filterwarnings('ignore')
INPUT, OUT = sys.argv[1], sys.argv[2]
SCHEMES = (('old', 10), ('gfix', 10), ('gnew', 10), ('gwide', 40))

def cases():
    out = []
    for path in sorted(glob.glob(os.path.join(INPUT, '*.txt'))):
        d = np.loadtxt(path); N = len(d); name = os.path.basename(path)[:-4]
        w = 20 if N >= 30 else max(8, N - 4); last = N - w + 1
        starts = sorted(set(np.round(np.linspace(1, last, min(8, last))).astype(int).tolist()))
        for (s, ww) in [(s, w) for s in starts] + [(1, N)]:
            for method in (0, 1):
                out.append((path, name, s, ww, method))
    return out

def run(case):
    path, name, s, w, method = case
    d = np.loadtxt(path); t = d[s-1:s-1+w, 0]; y = d[s-1:s-1+w, 1].astype(float).copy()
    base = dict(dataset=name, start=s, w=w, method=method, t_first=float(t[0]) if np.isfinite(t[0]) else None)
    if not (np.all(np.isfinite(y)) and np.all(np.isfinite(t)) and np.allclose(np.diff(t), 1.0)) or y.sum() < 10 or np.sum(y > 0) < 4:
        return [dict(base, scheme='skip')]
    if y[0] == 0: y[0] = 0.01
    fut = d[s-1+w:s-1+w+4, 1] if len(d) >= s-1+w+4 else None
    seed = zlib.crc32(repr(case[1:]).encode()); rows = []
    for sch, ns in SCHEMES:
        t0 = time.time()
        try:
            f = gp.fit_model(t, y, gp.GOM, 1, method, ns, sch, np.random.default_rng(seed))
        except Exception as e:
            rows.append(dict(base, scheme=sch, error=repr(e)[:120])); continue
        row = dict(base, scheme=sch, fval=f.fval, secs=round(time.time()-t0, 3), n_starts=f.n_starts, n_valid=f.n_valid, agree=f.agree)
        if f.P is not None:
            row.update(r=f.P[0], a=f.P[2], r_atLB=bool(f.atLB[0]), r_atUB=bool(f.atUB[0]), a_atLB=bool(f.atLB[2]), a_atUB=bool(f.atUB[2]),
                       LBr=float(f.LB[0]), UBr=float(f.UB[0]), LBa=float(f.LB[2]), UBa=float(f.UB[2]))
            if fut is not None:
                fc = gp.forecast(f.P, t, 4, gp.GOM, elapsed=sch in gp.ELAPSED)
                if fc is not None: row['mae4'] = float(np.mean(np.abs(fc[-4:] - fut)))
        rows.append(row)
    return rows

if __name__ == '__main__':
    cs = cases(); t0 = time.time()
    with Pool(2) as pool, open(OUT, 'w') as fh:
        for rows in pool.imap_unordered(run, cs, chunksize=4):
            for r in rows: fh.write(json.dumps(r) + '\n')
    print(len(cs), 'cases in', round(time.time()-t0), 's')
