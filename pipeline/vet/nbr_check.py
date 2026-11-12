# nbr_check.py: shared-systematic test with real neighbour light curves (fetches
# ZTF LCs of up to 8 field stars within 6' and |dG|<1.0 from IRSA). Reports LS power/FAP at the lead period (+-0.5%) and
# Pearson r of nightly-binned r-band residuals. usage: python3 nbr_check.py SOURCE
import sys,os,io,json,time,warnings,numpy as np,pandas as pd,requests;warnings.filterwarnings('ignore')
from astropy.timeseries import LombScargle
sid=sys.argv[1]; L=pd.read_csv('leads.csv',dtype={'source_id':str}).set_index('source_id').loc[sid]
o=json.load(open(f'json/{sid}.json')); P=o['zr']['P_near']; ra,de,G=o['ra'],o['dec'],o['G']
fg=pd.read_csv(os.path.dirname(L.ckpt)+'/field_gaia.csv',dtype={'Source':str})
fg['sep']=60*np.hypot((fg.RA_ICRS-ra)*np.cos(np.radians(de)),fg.DE_ICRS-de)
nb=fg[(fg.sep<6)&(fg.sep>0.25)&(abs(fg.Gmag-G)<1.0)&(fg.Source!=sid)].sort_values('sep').head(8)
def lc(s,r,d):
    p=f'nbr/{s}.csv'
    if not os.path.exists(p):
        for t in range(3):
            try:
                x=requests.get(f'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20{r}%20{d}%200.0004&FORMAT=csv',timeout=180)
                if x.text.startswith('oid'): open(p,'w').write(x.text); break
            except Exception: time.sleep(10)
        else: return None
    return pd.read_csv(p)
def prep(df,fc='zr'):
    d=df[(df.filtercode==fc)&(df.catflags==0)&(df.magerr<0.1)].copy()
    if len(d)<30: return None
    d['m']=d.mag-d.groupby(['field','ccdid','qid']).mag.transform('median'); d['night']=np.floor(d.hjd).astype(int); return d
def pw(d):
    ls=LombScargle(d.hjd,d.m,d.magerr);fr=np.linspace(0.995/P,1.005/P,201);p=ls.power(fr).max()
    return float(p),float(ls.false_alarm_probability(p,minimum_frequency=1/300,maximum_frequency=50))
tl=f'lc/{sid}.csv' if os.path.exists(f'lc/{sid}.csv') else f'{L.lcdir}/{sid}.csv'
d0=prep(pd.read_csv(tl)); n0=d0.groupby('night').m.mean(); tp=pw(d0); rows=[]
for r in nb.itertuples():
    df=lc(r.Source,r.RA_ICRS,r.DE_ICRS)
    if df is None: rows.append((r.Source,'LC fetch failed')); continue
    d=prep(df)
    if d is None: rows.append((r.Source,'too few points')); continue
    p,fap=pw(d); n=d.groupby('night').m.mean(); j=n0.index.intersection(n.index)
    cc=float(np.corrcoef(n0[j],n[j])[0,1]) if len(j)>20 else None
    rows.append((r.Source,round(r.sep,2),round(r.Gmag,2),round(p,4),fap,None if cc is None else round(cc,3)))
ok=[x for x in rows if len(x)==6]
# neighbours whose own blind top peak is at P (+-0.1%, x1/x2/x0.5)
top=[]; ens=[]
for r in nb.itertuples():
    p=f'nbr/{r.Source}.csv'
    if not os.path.exists(p): continue
    d=prep(pd.read_csv(p))
    if d is None: continue
    ls=LombScargle(d.hjd,d.m,d.magerr); f,pp=ls.autopower(minimum_frequency=1/300,maximum_frequency=50,samples_per_peak=10)
    Pt=1/f[pp.argmax()]
    if any(abs(Pt*k/P-1)<0.001 for k in (1,2,.5)): top.append((r.Source,round(float(Pt),6),round(r.sep,2)))
    ens.append(d.groupby('night').m.median().rename(r.Source))
# ensemble correction: subtract per-night median residual of neighbours from target, re-measure signal at P
from numpy.linalg import lstsq
def amp(t,y,P):
    ph=2*np.pi*t/P; A=np.column_stack([np.ones_like(t)]+[fn(h*ph) for h in (1,2) for fn in (np.sin,np.cos)]); c=lstsq(A,y,rcond=None)[0]; m=A@c; return float(m.max()-m.min())
a0=amp(d0.hjd.values,d0.m.values,P); corr=None
if len(ens)>=3:
    E=pd.concat(ens,axis=1).median(axis=1); dd=d0[d0.night.isin(E.index)].copy(); dd['m2']=dd.m-dd.night.map(E)
    l2=LombScargle(dd.hjd,dd.m2,dd.magerr); fr=np.linspace(0.995/P,1.005/P,201); p2=l2.power(fr).max()
    corr=dict(n=len(dd),amp_before=a0,amp_after=amp(dd.hjd.values,dd.m2.values,P),fap_after=float(l2.false_alarm_probability(p2,minimum_frequency=1/300,maximum_frequency=50)))
out=dict(source_id=sid,P=P,target_pow_fap=tp,n_nbr=len(ok),n_nbr_sigP=sum(x[4]<1e-6 for x in ok),n_nbr_corr03=sum((x[5] or 0)>0.3 for x in ok),nbr_top_at_P=top,ensemble_corrected=corr,nbr=rows)
json.dump(out,open(f'nbr/{sid}_nbr.json','w'),indent=1); print(sid,out['n_nbr'],'top@P',len(top),'corr',corr)
