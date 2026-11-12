# vet_backlog.py: first-stage vet for one lead (ZTF r/g periods, halves, fold, aliases,
# field-systematic count from search ckpt, Gaia DR3 neighbours (IRSA TAP), catalog xmatch (VizieR CfA mirror,
# VSX live API, SIMBAD Harvard mirror). Writes json/<id>.json. Any failed query is recorded as ERR, never guessed.
import sys,os,io,json,time,warnings,numpy as np,pandas as pd,requests
warnings.filterwarnings('ignore')
from astropy.timeseries import LombScargle
L=pd.read_csv('leads.csv',dtype={'source_id':str}).set_index('source_id')
def gaia_row(sid,ck):
    d=os.path.dirname(ck); f=pd.read_csv(f'{d}/field_gaia.csv',dtype={'Source':str}).set_index('Source'); return f.loc[sid]
def get(url,params=None,data=None,tries=3,timeout=90):
    last=None
    for t in range(tries):
        try:
            r=(requests.post(url,data=data,timeout=timeout) if data else requests.get(url,params=params,timeout=timeout))
            if r.ok: return r
            last=f'HTTP {r.status_code}'
        except Exception as e: last=str(e)[:80]
        time.sleep(5*(t+1))
    raise RuntimeError(last)
def fetch_lc(sid,ra,de,lcdir):
    p=f'{lcdir}/{sid}.csv' if lcdir else ''
    if p and os.path.exists(p) and os.path.getsize(p)>0: return pd.read_csv(p),'local'
    q=f'lc/{sid}.csv'
    if os.path.exists(q) and os.path.getsize(q)>200: return pd.read_csv(q),'irsa(cached)'
    r=get(f'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20{ra}%20{de}%200.0004&FORMAT=csv&BAD_CATFLAGS_MASK=32768',timeout=180)
    if not r.text.startswith('oid'): raise RuntimeError('IRSA LC bad response: '+r.text[:80])
    open(q,'w').write(r.text); return pd.read_csv(io.StringIO(r.text)),'irsa'
def prep(df,fc):
    d=df[(df.filtercode==fc)&(df.catflags==0)&(df.magerr<0.1)].copy()
    if len(d)<30: return None
    d['m']=d.mag-d.groupby(['field','ccdid','qid']).mag.transform('median'); return d.sort_values('hjd')
def harm(t,y,P,nh=2):
    ph=2*np.pi*t/P; A=np.column_stack([np.ones_like(t)]+[fn(h*ph) for h in range(1,nh+1) for fn in (np.sin,np.cos)])
    c=np.linalg.lstsq(A,y,rcond=None)[0]; mod=A@c; return mod,float(mod.max()-mod.min())
def band(d,P0):
    t,m,e=d.hjd.values,d.m.values,d.magerr.values; ls=LombScargle(t,m,e)
    f,p=ls.autopower(minimum_frequency=1/300,maximum_frequency=50,samples_per_peak=10); i=p.argmax()
    fr=np.linspace(0.99/P0,1.01/P0,2001); pp=ls.power(fr); Pn=1/fr[pp.argmax()]
    o=dict(n=len(t),P_blind=float(1/f[i]),fap_blind=float(ls.false_alarm_probability(p[i],minimum_frequency=1/300,maximum_frequency=50)),
           P_near=float(Pn),fap_near=float(ls.false_alarm_probability(pp.max(),minimum_frequency=1/300,maximum_frequency=50)))
    # 5 strongest independent peaks, for alias reporting
    idx=np.argsort(p)[::-1]; top=[]
    for j in idx:
        if all(abs(f[j]-f[k])>0.002*f[k]+1/3000 for k in top): top.append(j)
        if len(top)==5: break
    o['top5']=[round(float(1/f[j]),6) for j in top]
    h=len(t)//2; hh=[]
    for x in (slice(0,h),slice(h,None)):
        l2=LombScargle(t[x],m[x],e[x]); hh.append(float(1/fr[l2.power(fr).argmax()]))
    o['halves']=hh; o['stable']=abs(hh[0]/hh[1]-1)<0.003
    mod,amp=harm(t,m,Pn); res=m-mod
    o['amp_model']=amp; o['R_fold']=float(1-np.var(res)/np.var(m)); o['rms']=float(np.std(m)); o['snr']=float(amp/(np.std(res)/np.sqrt(len(t))*2))
    return o
def alias(P):
    bad=[1,0.5,1/3,0.25,1.0027]; 
    if any(abs(P/b-1)<0.01 for b in bad): return 'near 1-day harmonic'
    if abs(P/1.035-1)<0.005 or abs(P/29.53-1)<0.03 or abs(P/14.77-1)<0.03: return 'lunar alias'
    if abs(P/365.25-1)<0.05 or abs(P/182.6-1)<0.05: return 'yearly alias'
    if P>100: return 'P>100 d (beyond rule)'
    return ''
def viz(cat,ra,de,rs):
    errs=[]
    for H in ('https://vizier.cfa.harvard.edu','https://vizier.idia.ac.za'):
        try: r=get(H+'/viz-bin/asu-tsv',params={'-source':cat,'-c':f'{ra} {de}','-c.rs':rs,'-out.max':50,'-out.add':'_r','-sort':'_r'},tries=2)
        except Exception as e: errs.append(H[8:]+' '+str(e)); continue
        if 'Error=' in r.text or 'QUERY_STATUS=ERROR' in r.text or '#Table' not in r.text and '#RESOURCE' not in r.text:
            errs.append(H[8:]+' db error'); continue
        lines=[l for l in r.text.splitlines() if l and not l.startswith('#')]
        if len(lines)<3: return []
        hdr=lines[0].split('\t'); rows=[dict(zip(hdr,l.split('\t'))) for l in lines[3:]]
        return [{k:v.strip() for k,v in list(z.items())[:14] if v.strip()} for z in rows]
    raise RuntimeError('; '.join(errs))
CATS={'VSX':'B/vsx/vsx','ZTF_Chen2020':'J/ApJS/249/18/table2','ATLAS':'J/AJ/156/241/table4','ASASSN':'II/366/catv2021',
      'GaiaDR3_vclass':'I/358/vclassre','GaiaDR3_varisum':'I/358/varisum','TESS_EB_Prsa2022':'J/ApJS/258/16/tess-ebs'}
def main(sid):
    out=f'json/{sid}.json'
    if os.path.exists(out): return
    l=L.loc[sid]; g=gaia_row(sid,l.ckpt); ra,de,G=float(g.RA_ICRS),float(g.DE_ICRS),float(g.Gmag)
    O=dict(source_id=sid,field=l.field,ra=ra,dec=de,G=G,bp_rp=float(g['BP-RP']) if pd.notna(g['BP-RP']) else None,plx=float(g.Plx) if pd.notna(g.Plx) else None,note=l.note if isinstance(l.note,str) else '')
    ck=[json.loads(x) for x in open(l.ckpt)]; c=pd.DataFrame([x for x in ck if 'P_zr' in x]); c['Source']=c.Source.astype(str)
    me=c[c.Source==sid]; P0=float(me.P_zr.iloc[0]) if len(me) else None
    O['P_search']=P0
    try: df,O['lc_src']=fetch_lc(sid,ra,de,l.lcdir if isinstance(l.lcdir,str) else '')
    except Exception as e: O['lc_err']=str(e); json.dump(O,open(out,'w'),indent=1); return
    if P0 is None: P0=float(LombScargle(*(lambda d:(d.hjd,d.m))(prep(df,'zr'))).autopower(minimum_frequency=1/300,maximum_frequency=50)[0][0])
    for fc in ['zr','zg']:
        d=prep(df,fc); O[fc]=band(d,P0) if d is not None else None
    O['alias']=alias(P0)
    # field systematic: other field stars whose search period (r) is within 0.3% of P0 (or x2, /2) with fap<1e-6
    cc=c[(c.Source!=sid)&(c.fap_zr<1e-6)]
    O['field_n_scored']=int(len(c)); O['field_same_P']=int(sum(any(abs(p*k/P0-1)<0.003 for k in (1,2,.5)) for p in cc.P_zr))
    # Gaia DR3 neighbours within 15" (IRSA TAP)
    try:
        q=f"SELECT source_id,phot_g_mean_mag,ra,dec FROM gaia_dr3_source WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{de},{15/3600}))=1"
        n=pd.read_csv(io.StringIO(get('https://irsa.ipac.caltech.edu/TAP/sync',data={'QUERY':q,'FORMAT':'csv'}).text),dtype={'source_id':str})
        n['sep']=3600*np.hypot((n.ra-ra)*np.cos(np.radians(de)),n.dec-de); n=n[n.source_id!=sid].sort_values('sep')
        O['neigh15']=[(r.source_id,round(r.sep,2),None if pd.isna(r.phot_g_mean_mag) else round(r.phot_g_mean_mag,2)) for r in n.itertuples()]
    except Exception as e: O['neigh15']='ERR '+str(e)
    X={}
    for k,cat in CATS.items():
        try: X[k]=viz(cat,ra,de,10)
        except Exception as e: X[k]='ERR '+str(e)
    try:
        j=get('https://vsx.aavso.org/index.php',params=dict(view='api.list',ra=ra,dec=de,radius=30/3600,format='json')).json()
        o=j.get('VSXObjects',{}); o=o.get('VSXObject',[]) if isinstance(o,dict) else []
        if isinstance(o,dict): o=[o]
        X['VSX_live30']=[(z.get('Name'),z.get('VariabilityType'),z.get('Period'),round(3600*np.hypot((float(z['RA2000'])-ra)*np.cos(np.radians(de)),float(z['Declination2000'])-de),1)) for z in o]
    except Exception as e: X['VSX_live30']='ERR '+str(e)
    try:
        Q=f"SELECT b.main_id,b.otype,DISTANCE(POINT('ICRS',b.ra,b.dec),POINT('ICRS',{ra},{de}))*3600 AS sep FROM basic b WHERE CONTAINS(POINT('ICRS',b.ra,b.dec),CIRCLE('ICRS',{ra},{de},{10/3600}))=1"
        r=get('https://simbad.harvard.edu/simbad/sim-tap/sync',data=dict(REQUEST='doQuery',LANG='ADQL',FORMAT='json',QUERY=Q)).json()
        X['SIMBAD10']=sorted([(a,b,round(s,1)) for a,b,s in r['data']],key=lambda z:z[2])
    except Exception as e: X['SIMBAD10']='ERR '+str(e)
    O['xmatch']=X
    json.dump(O,open(out,'w'),indent=1,default=str); print(sid,'ok',flush=True)
if __name__=='__main__':
    from concurrent.futures import ThreadPoolExecutor
    ids=sys.argv[1:] or list(L.index)
    def w(s):
        try: main(s)
        except Exception as e: print(s,'FAIL',repr(e)[:200],flush=True)
    with ThreadPoolExecutor(4) as ex: list(ex.map(w,ids))
