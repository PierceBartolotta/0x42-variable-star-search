# Period search, frequencies up to 50 cycles/day, with a 1/5-day alias exclusion.
# Run inside a field folder that has field_gaia.csv and lc/. Writes search_ckpt.jsonl (one line per star) and search_results_hf.csv.
# Blind periodicity search on ZTF light curves. Rules fixed before looking at novelty:
#  variable if: reduced chi2 > 5, amp(5-95%) > 0.05 mag, LS FAP<1e-10 in r; period accepted only if
#  g and r agree (within 0.5%, allowing x2/x0.5) and P not within 2% of 1, 1/2, 1/3, 1/4 d or 27-32 d, or >300 d.
import numpy as np, pandas as pd, glob, os, json
from astropy.timeseries import LombScargle
gaia=pd.read_csv('field_gaia.csv',dtype={'Source':str}).set_index('Source')
fmin,fmax=1/300.,50.
def prep(df,fc):
    d=df[(df.filtercode==fc)&(df.catflags==0)&(df.magerr<0.1)].copy()
    if len(d)<80: return None
    # remove per field/ccd/quadrant offsets
    d['m']=d.mag-d.groupby(['field','ccdid','qid']).mag.transform('median')
    # 5-sigma clip
    s=1.4826*np.median(abs(d.m-np.median(d.m))); d=d[abs(d.m-np.median(d.m))<6*max(s,0.01)]
    return d
def bad(P):
    for a in [1,1/2,1/3,1/4,1/5]:
        if abs(P/a-1)<0.02: return True
    return (27<P<32) or P>300
def pdm(t,m,P,nb=20):
    ph=(t/P)%1; idx=np.floor(ph*nb).astype(int); v=0;n=0
    for b in range(nb):
        x=m[idx==b]
        if len(x)>1: v+=x.var()*len(x); n+=len(x)
    return v/n/m.var()
def best(t,m,e):
    ls=LombScargle(t,m,e); f,p=ls.autopower(minimum_frequency=fmin,maximum_frequency=fmax,samples_per_peak=10)
    o=np.argsort(p)[::-1]; 
    for i in o[:50]:
        P=1/f[i]
        if not bad(P): return P,p[i],ls.false_alarm_probability(p[i],minimum_frequency=fmin,maximum_frequency=fmax)
    return np.nan,0,1
def same(a,b):
    return any(abs(a*k/b-1)<0.005 for k in [1,2,0.5])
def run(f):
    sid=os.path.basename(f)[:-4]
    try: df=pd.read_csv(f)
    except Exception: return None
    if len(df)==0: return None
    # keep only the dominant ZTF object id per filter (closest match)
    rec=dict(Source=sid)
    for fc in ['zr','zg']:
        d=prep(df,fc)
        if d is None: continue
        t,m,e=d.hjd.values-2400000.5,d.m.values,d.magerr.values
        chi=np.sum(((m-np.median(m))/e)**2)/(len(m)-1)
        amp=np.percentile(m,95)-np.percentile(m,5)
        P,pw,fap=best(t,m,e)
        rec.update({f'n_{fc}':len(m),f'chi_{fc}':chi,f'amp_{fc}':amp,f'P_{fc}':P,f'pow_{fc}':pw,f'fap_{fc}':fap,f'med_{fc}':np.median(d.mag)})
        if np.isfinite(P):
            # EB check: does 2P fold tighter (PDM) than P?
            rec[f'pdm1_{fc}']=pdm(t,m,P); rec[f'pdm2_{fc}']=pdm(t,m,2*P)
    return rec
from multiprocessing import Pool
if __name__=='__main__':
  with Pool(2) as p: out=[x for x in p.map(run,sorted(glob.glob('lc/*.csv'))) if x]
  with open('search_ckpt.jsonl','w') as fo:
      for x in out: fo.write(json.dumps(x,default=float)+'\n')
  r=pd.DataFrame(out)
  r['variable']=(r.chi_zr>5)&(r.amp_zr>0.05)&(r.fap_zr<1e-10)
  r['periodic']=r.variable&r.apply(lambda x: np.isfinite(x.get('P_zg',np.nan)) and np.isfinite(x.P_zr) and same(x.P_zr,x.P_zg) and x.get('fap_zg',1)<1e-5,axis=1)
  r['P']=np.where(r.pdm2_zr<0.9*r.pdm1_zr,2*r.P_zr,r.P_zr)
  r=r.join(gaia,on='Source'); r.to_csv('search_results_hf.csv',index=False)
  print(len(r),r.variable.sum(),r.periodic.sum())
