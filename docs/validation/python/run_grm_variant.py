import sys, json, zlib, collections, warnings, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__))); import gp_port as gp
warnings.filterwarnings('ignore')
INPUT=sys.argv[1]; RES=sys.argv[2]
rows=[json.loads(l) for l in open(RES+'/matrix.jsonl')]
key=lambda r:(r['dataset'],r['start'],r['w'],r['model'],r['method'],r['fixI0'])
by=collections.defaultdict(dict)
for r in rows: by[key(r)][r['scheme']]=r
grm={k:v for k,v in by.items() if k[3]=='GRM' and 'new' in v}
rel=lambda a,b:(a-b)/max(abs(b),1.0)
hit=[k for k,v in grm.items() if v['new']['r_atUB'] or rel(v['new']['fval'],min(v[s]['fval'] for s in v))>1e-4]
rng=np.random.default_rng(3); others=[k for k in grm if k not in hit]; idx=rng.choice(len(others),40,replace=False)
sel=hit+[others[i] for i in idx]
out=[]
for k in sel:
    name,s,w,_,method,fixI0=k
    d=np.loadtxt(f"{INPUT}/{name}.txt"); t=d[s-1:s-1+w,0]; y=d[s-1:s-1+w,1].astype(float).copy()
    if y[0]==0: y[0]=0.01
    f=gp.fit_model(t,y,gp.GRM,1,method,10,'newa',np.random.default_rng(zlib.crc32(repr(k).encode())))
    v=grm[k]; ref=min(min(v[x]['fval'] for x in v), f.fval)
    out.append(dict(k=k, group='hit' if k in hit else 'other', fa=f.fval, fnew=v['new']['fval'], ref=ref, r_on=bool(f.atLB[0] or f.atUB[0]), a_lb=bool(f.atLB[2]), a=float(f.P[2]), r=float(f.P[0]), agree=f.agree, valid=f.n_valid/f.n_starts))
json.dump(out,open(RES+'/grm_variant.json','w'))
for g in ('hit','other'):
    o=[x for x in out if x['group']==g]
    print(f"{g}: n={len(o)} | variant at best-known {np.mean([rel(x['fa'],x['ref'])<=1e-4 for x in o]):.2f} (current new: {np.mean([rel(x['fnew'],x['ref'])<=1e-4 for x in o]):.2f}) | variant better than new {sum(rel(x['fnew'],x['fa'])>1e-4 for x in o)}, worse {sum(rel(x['fa'],x['fnew'])>1e-4 for x in o)} | r on bound {sum(x['r_on'] for x in o)} | a at 0.05 {sum(x['a_lb'] for x in o)} | max gap variant {max(rel(x['fa'],x['ref']) for x in o):.3g}, new {max(rel(x['fnew'],x['ref']) for x in o):.3g} | valid starts {np.mean([x['valid'] for x in o]):.3f} | agreement {np.mean([x['agree'] for x in o]):.2f}")
