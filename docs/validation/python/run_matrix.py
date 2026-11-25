"""Old vs proposed r box on every bundled dataset. Writes one JSON line per (case, scheme)."""
import sys, os, glob, json, time, zlib, warnings
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import gp_port as gp
from multiprocessing import Pool
warnings.filterwarnings('ignore')
INPUT, OUT = sys.argv[1], sys.argv[2]

def cases():
    out = []
    for path in sorted(glob.glob(os.path.join(INPUT, '*.txt'))):
        d = np.loadtxt(path); N = len(d); name = os.path.basename(path)[:-4]
        w = 20 if N >= 30 else max(8, N - 4)
        last = N - w + 1
        starts = sorted(set(np.round(np.linspace(1, last, min(8, last))).astype(int).tolist()))
        wins = [(s, w) for s in starts] + [(1, N)]
        for (s, ww) in wins:
            for flag in gp.FREE_P:
                for method in (0, 1):
                    out.append((path, name, s, ww, flag, method, 1))
            if ww == w and s in starts[::2]:
                out.append((path, name, s, ww, gp.GLM, 0, 0))      # fixI0 = 0 subset
    return out

def run(case):
    try:
        return _run(case)
    except Exception as e:
        path, name, s, w, flag, method, fixI0 = case
        return [dict(dataset=name, start=s, w=w, model=gp.NAMES[flag], method=method, fixI0=fixI0, scheme='error', why=repr(e)[:200])]

def _run(case):
    path, name, s, w, flag, method, fixI0 = case
    d = np.loadtxt(path); t = d[s-1:s-1+w, 0]; y = d[s-1:s-1+w, 1].astype(float).copy()
    base = dict(dataset=name, start=s, w=w, model=gp.NAMES[flag], method=method, fixI0=fixI0, ymax=float(y.max()), y1=float(y[0]), ysum=float(y.sum()))
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(t)):
        return [dict(base, scheme='skip', why='non-finite data')]
    if not np.allclose(np.diff(t), 1.0):
        return [dict(base, scheme='skip', why='time index is not unit-spaced')]
    if y.sum() < 10 or np.sum(y > 0) < 4:
        return [dict(base, scheme='skip', why='too few counts')]
    if y[0] == 0 and fixI0 == 1:
        y[0] = 0.01                                        # as the runners do
    fut = d[s-1+w:s-1+w+4, 1] if len(d) >= s-1+w+4 else None
    seed = zlib.crc32(repr(case[1:]).encode())
    rows, fits = [], {}
    for sch, ns in (('old', 10), ('new', 10), ('i0', 10), ('new40', 40), ('i040', 40)):
        t0 = time.time()
        f = gp.fit_model(t, y, flag, fixI0, method, ns, sch[:3].rstrip('4'), np.random.default_rng(seed))
        fits[sch] = f
        row = dict(base, scheme=sch, fval=f.fval, secs=round(time.time()-t0, 3), n_starts=f.n_starts, n_valid=f.n_valid, agree=f.agree, widened=bool(f.widened))
        if f.P is not None:
            row.update(r=f.P[0], p=f.P[1], a=f.P[2], K=f.P[3], I0=f.P[4], slot1=float(f.slot1),
                       LB1=float(f.LB[0]), UB1=float(f.UB[0]),
                       r_atLB=bool(f.atLB[0]), r_atUB=bool(f.atUB[0]), p_atLB=bool(f.atLB[1]), p_atUB=bool(f.atUB[1]),
                       K_atLB=bool(f.atLB[3]), K_atUB=bool(f.atUB[3]), a_atLB=bool(f.atLB[2]), a_atUB=bool(f.atUB[2]))
            if fut is not None:
                fc = gp.forecast(f.P, t, 4, flag)
                if fc is not None:
                    row['mae4'] = float(np.mean(np.abs(fc[-4:] - fut)))
                    row['fut_mean'] = float(np.mean(fut))
        rows.append(row)
    return rows

if __name__ == '__main__':
    cs = cases()
    seen = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            r = json.loads(l); seen.add((r['dataset'], r['start'], r['w'], r['model'], r['method'], r['fixI0']))
    cs = [c for c in cs if (c[1], c[2], c[3], gp.NAMES[c[4]], c[5], c[6]) not in seen]
    print(len(cs), 'cases to run', flush=True)
    t0 = time.time(); done = 0
    with Pool(2) as pool, open(OUT, 'a') as fh:
        for rows in pool.imap_unordered(run, cs, chunksize=4):
            for r in rows: fh.write(json.dumps(r) + '\n')
            fh.flush(); done += 1
            if done % 50 == 0: print(done, round(time.time()-t0), 's', flush=True)
    print('finished', round(time.time()-t0), 's', flush=True)
