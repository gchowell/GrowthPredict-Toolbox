import json, sys, collections, numpy as np
rows = [json.loads(l) for l in open(sys.argv[1])]
key = lambda r: (r['dataset'], r['start'], r['w'], r['method'])
by = collections.defaultdict(dict)
for r in rows: by[key(r)][r['scheme']] = r
cases = {k: v for k, v in by.items() if all(s in v and 'fval' in v[s] for s in ('old', 'gfix', 'gnew', 'gwide'))}
print(f"Gompertz fits: {len(cases)} cases ({len({k[:3] for k in cases})} windows of {len({k[0] for k in cases})} series, least squares and Poisson); skipped {sum('skip' in v for v in by.values())}; errors {sum(any('error' in x for x in v.values()) for v in by.values())}")
rel = lambda a, b: (a - b)/max(abs(b), 1.0); TOL = 1e-4
ref = {k: min(v[s]['fval'] for s in v) for k, v in cases.items()}
def line(ks, label):
    print(f"\n=== {label}: {len(ks)} cases")
    for s, nm in (('old', 'HEAD'), ('gfix', 'elapsed time (adopted)'), ('gnew', 'elapsed time, wide box')):
        v = [cases[k][s] for k in ks]
        print(f"  {nm:<24}: best-known objective {100*np.mean([rel(cases[k][s]['fval'], ref[k]) <= TOL for k in ks]):5.1f}% | r on bound {100*np.mean([x['r_atLB'] or x['r_atUB'] for x in v]):5.1f}% (LB {sum(x['r_atLB'] for x in v)}, UB {sum(x['r_atUB'] for x in v)}) | a on bound {100*np.mean([x['a_atLB'] or x['a_atUB'] for x in v]):5.1f}% (LB {sum(x['a_atLB'] for x in v)}, UB {sum(x['a_atUB'] for x in v)}) | valid starts {100*np.mean([x['n_valid']/x['n_starts'] for x in v]):5.1f}% | agreement {np.nanmean([x['agree'] for x in v]):.2f} | agreement < 0.5: {100*np.mean([x['agree'] < 0.5 for x in v]):.1f}% | {np.mean([x['secs'] for x in v]):.2f}s")
    for a_, b_ in (('old', 'gfix'), ('gfix', 'gnew'), ('old', 'gnew')):
        better = [k for k in ks if rel(cases[k][a_]['fval'], cases[k][b_]['fval']) > TOL]; worse = [k for k in ks if rel(cases[k][b_]['fval'], cases[k][a_]['fval']) > TOL]
        msg = f"  {b_} vs {a_}: better {len(better)}, same {len(ks)-len(better)-len(worse)}, worse {len(worse)}"
        ls = [k for k in better if k[3] == 0]
        if ls:
            ra = np.array([cases[k][a_]['fval']/max(cases[k][b_]['fval'], 1e-9) for k in ls]); msg += f" | least squares, where better: SSE ratio median {np.median(ra):.2f}, >2x in {100*np.mean(ra > 2):.0f}%, >10x in {100*np.mean(ra > 10):.0f}%"
        print(msg)
        for k in worse[:4]: print("     WORSE", k, cases[k][a_]['fval'], cases[k][b_]['fval'], 'ref', ref[k], 'valid', cases[k][b_]['n_valid'], 'agree', cases[k][b_]['agree'])
    for a_, b_, lab in (('old', 'gfix', 'adopted vs HEAD'), ('gfix', 'gnew', 'wide box vs adopted')):
        fm = np.array([(cases[k][a_]['mae4'], cases[k][b_]['mae4']) for k in ks if 'mae4' in cases[k][a_] and 'mae4' in cases[k][b_]])
        if len(fm):
            ch = np.abs(fm[:, 0] - fm[:, 1]) > 1e-3*np.maximum(fm[:, 0], 1)
            print(f"  4-step forecast MAE, {lab} (n={len(fm)}): changed {ch.sum()}; better {np.sum(fm[ch,1] < fm[ch,0])}, worse {np.sum(fm[ch,1] > fm[ch,0])}; median ratio {np.median(fm[ch,1]/np.maximum(fm[ch,0],1e-9)) if ch.any() else float('nan'):.2f}")
allk = list(cases)
line(allk, "ALL")
line([k for k in allk if cases[k]['old']['t_first'] in (0.0,) and k[1] == 1], "windows whose first time label is 0")
line([k for k in allk if not (cases[k]['old']['t_first'] == 0.0)], "windows whose first time label is not 0")
miss = [k for k in allk if rel(cases[k]['gnew']['fval'], ref[k]) > TOL]
print("\nnew misses best-known:", len(miss))
for k in miss[:8]:
    v = cases[k]; print("  ", k, "new %.6g ref %.6g (gap %.2e) r=%.3g a=%.3g | ref r=%.3g a=%.3g | new agree %.2f valid %d/%d" % (v['gnew']['fval'], ref[k], rel(v['gnew']['fval'], ref[k]), v['gnew']['r'], v['gnew']['a'], v['gwide']['r'], v['gwide']['a'], v['gnew']['agree'], v['gnew']['n_valid'], v['gnew']['n_starts']))
g = [cases[k]['gnew'] for k in allk]
print("\nnew-box estimates: r range [%.3g, %.3g], 1st-99th pct [%.3g, %.3g]; a range [%.3g, %.3g], 1st-99th pct [%.3g, %.3g]" % (min(x['r'] for x in g), max(x['r'] for x in g), *np.quantile([x['r'] for x in g], [.01, .99]), min(x['a'] for x in g), max(x['a'] for x in g), *np.quantile([x['a'] for x in g], [.01, .99])))
print("HEAD box with elapsed time: r below its LB would be needed in", sum(cases[k]['gnew']['r'] < cases[k]['gfix']['LBr']*0.999 for k in allk), "cases; r above UB in", sum(cases[k]['gnew']['r'] > cases[k]['gfix']['UBr']*1.001 for k in allk), "; a below LB in", sum(cases[k]['gnew']['a'] < cases[k]['gfix']['LBa']*0.999 for k in allk), "; a above UB in", sum(cases[k]['gnew']['a'] > cases[k]['gfix']['UBa']*1.001 for k in allk))
