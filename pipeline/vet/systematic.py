# Field-systematic test: how many other field stars (search ckpt, FAP_r<1e-6) have P within 0.1% of the lead's P (x1, x2, x0.5),
# compared with the local density of periods (+-5% window) -> Poisson probability. Flags a shared-period systematic.
import json,pandas as pd,numpy as np
from scipy.stats import poisson
L=pd.read_csv('leads.csv',dtype={'source_id':str}); out=[]; C={}
for r in L.itertuples():
    if r.ckpt not in C:
        c=pd.DataFrame([x for x in map(json.loads,open(r.ckpt)) if 'P_zr' in x]); c['Source']=c.Source.astype(str); C[r.ckpt]=c
    c=C[r.ckpt]; me=c[c.Source==r.source_id]
    if not len(me): out.append(dict(source_id=r.source_id,sys_n=None)); continue
    P=float(me.P_zr.iloc[0]); o=c[(c.Source!=r.source_id)&(c.fap_zr<1e-6)].P_zr.values
    n=0; wide=0
    for k in (1,2,.5):
        x=np.abs(o*k/P-1); n+=int((x<0.001).sum()); wide+=int(((x<0.05)&(x>=0.001)).sum())
    exp=wide*0.002/0.098; p=float(poisson.sf(n-1,max(exp,1e-3))) if n>0 else 1.0
    out.append(dict(source_id=r.source_id,P=P,sys_n=n,sys_expected=round(exp,2),sys_p=p))
pd.DataFrame(out).to_csv('systematic.csv',index=False); print(pd.DataFrame(out).sort_values('sys_p').head(30).to_string())
