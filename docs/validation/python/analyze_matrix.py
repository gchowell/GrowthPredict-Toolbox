import json, sys, collections, numpy as np
rows = [json.loads(l) for l in open(sys.argv[1])]
key = lambda r: (r['dataset'], r['start'], r['w'], r['model'], r['method'], r['fixI0'])
by = collections.defaultdict(dict)
for r in rows: by[key(r)][r['scheme']] = r
skips = [k for k, v in by.items() if 'skip' in v]; errs = [k for k, v in by.items() if 'error' in v]
cases = {k: v for k, v in by.items() if all(s in v for s in ('old', 'new', 'i0', 'new40', 'i040'))}
print(f"cases fitted: {len(cases)}  | skipped: {len(skips)} ({collections.Counter(by[k]['skip'].get('why') for k in skips)}) | errors: {len(errs)}")
print("datasets:", len({k[0] for k in cases}))
TOL = 1e-4
def rel(a, b):   # relative excess of a over b, robust for NLL (can be negative) and SSE
    return (a - b) / max(abs(b), 1.0)
def summarize(sel, label):
    ks = [k for k in cases if sel(k)]
    if not ks: return
    ref = {k: min(cases[k][s]['fval'] for s in cases[k]) for k in ks}
    print(f"\n=== {label}: {len(ks)} cases")
    for s in ('old', 'new', 'i0'):
        v = [cases[k][s] for k in ks]
        at_ref = np.mean([rel(cases[k][s]['fval'], ref[k]) <= TOL for k in ks])
        rb = np.mean([x.get('r_atLB', False) or x.get('r_atUB', False) for x in v])
        pb = np.mean([x.get('p_atLB', False) or x.get('p_atUB', False) for x in v])
        Kb = np.mean([x.get('K_atLB', False) or x.get('K_atUB', False) for x in v])
        valid = np.mean([x['n_valid'] / x['n_starts'] for x in v])
        agree = np.nanmean([x['agree'] for x in v])
        nofit = sum(x['n_valid'] == 0 for x in v)
        secs = np.mean([x['secs'] for x in v])
        print(f"  {s:>3}: reaches best-known objective {100*at_ref:5.1f}% | slot-1 (r or i0) on bound {100*rb:5.1f}% | p on bound {100*pb:5.1f}% | K on bound {100*Kb:5.1f}% | valid starts {100*valid:5.1f}% | start agreement {agree:.2f} | no fit {nofit} | {secs:.2f}s/fit")
    worse = [k for k in ks if rel(cases[k]['new']['fval'], cases[k]['old']['fval']) > TOL]
    better = [k for k in ks if rel(cases[k]['old']['fval'], cases[k]['new']['fval']) > TOL]
    print(f"  new vs old objective: better in {len(better)} ({100*len(better)/len(ks):.1f}%), same in {len(ks)-len(better)-len(worse)}, worse in {len(worse)}")
    if better:
        if ks[0][4] == 0:
            ratio = np.array([cases[k]['old']['fval'] / max(cases[k]['new']['fval'], 1e-12) for k in better])
            print(f"    where better: SSE(old)/SSE(new) median {np.median(ratio):.2f}, 90th pct {np.quantile(ratio,.9):.1f}, max {ratio.max():.1f}")
        else:
            dn = np.array([cases[k]['old']['fval'] - cases[k]['new']['fval'] for k in better])
            print(f"    where better: NLL(old)-NLL(new) median {np.median(dn):.1f}, 90th pct {np.quantile(dn,.9):.1f}, max {dn.max():.1f}")
        dp = np.array([abs(cases[k]['old']['p'] - cases[k]['new']['p']) for k in better])
        print(f"    where better: |p(old)-p(new)| median {np.median(dp):.3f}, share > 0.1: {100*np.mean(dp>0.1):.0f}%")
    for k in worse[:6]:
        o, n = cases[k]['old'], cases[k]['new']
        print(f"    WORSE {k}: old {o['fval']:.6g} new {n['fval']:.6g} ref {ref[k]:.6g}  new valid {n['n_valid']}/{n['n_starts']} agree {n['agree']:.2f}")
    miss = [k for k in ks if rel(cases[k]['new']['fval'], ref[k]) > TOL]
    print(f"  new misses best-known in {len(miss)} cases; of these, 40-start run in the same box finds it: {sum(rel(cases[k]['new40']['fval'], ref[k]) <= TOL for k in miss)}")
    fm = [(cases[k]['old']['mae4'], cases[k]['new']['mae4'], cases[k]['old'].get('fut_mean', np.nan)) for k in ks if 'mae4' in cases[k]['old'] and 'mae4' in cases[k]['new']]
    if fm:
        fm = np.array(fm); diff = fm[:, 0] != fm[:, 1]
        ch = np.abs(fm[:, 0] - fm[:, 1]) > 1e-3 * np.maximum(fm[:, 0], 1)
        print(f"  4-step forecast MAE (n={len(fm)}): changed in {ch.sum()}; among changed: new better in {np.sum(fm[ch,1] < fm[ch,0])}, worse in {np.sum(fm[ch,1] > fm[ch,0])}; median MAE old {np.median(fm[ch,0]):.1f} -> new {np.median(fm[ch,1]):.1f}; median ratio new/old {np.median(fm[ch,1]/np.maximum(fm[ch,0],1e-9)):.2f}")
    return ks
summarize(lambda k: k[5] == 1, "ALL (fixI0 = 1)")
for m in ('GGM', 'GLM', 'GRM'):
    for meth in (0, 1):
        summarize(lambda k, m=m, meth=meth: k[3] == m and k[4] == meth and k[5] == 1, f"{m}, {'least squares' if meth == 0 else 'Poisson MLE'}")
summarize(lambda k: k[5] == 0, "GLM least squares, fixI0 = 0")
print("\nnew r box hits (fixI0=1):", collections.Counter((k[3], 'UB' if cases[k]['new'].get('r_atUB') else 'LB') for k in cases if k[5]==1 and (cases[k]['new'].get('r_atUB') or cases[k]['new'].get('r_atLB'))))
print("old: widened and still on r bound:", sum(1 for k in cases if cases[k]['old']['widened'] and (cases[k]['old'].get('r_atUB') or cases[k]['old'].get('r_atLB'))), "| widened total:", sum(cases[k]['old']['widened'] for k in cases))
print("per dataset (fixI0=1): share of cases where new improves on old")
dd = collections.defaultdict(list)
for k in cases:
    if k[5] == 1: dd[k[0]].append(rel(cases[k]['old']['fval'], cases[k]['new']['fval']) > TOL)
for d, v in sorted(dd.items()): print(f"   {d[:62]:<62} {sum(v):>3}/{len(v)}")
