# Flag rule: same period in ZTF r and g (within 0.3%, x1/x2/x0.5), FAP<1e-6 in both, no 1-day/lunar aliases, P<100 d.
# Then drop stars Gaia DR3 already lists as variable (vari_summary or phot_variable_flag).
# usage: python3 05_flag_rule.py FIELD_NAME [FIELD_DIR]   (FIELD_DIR holds search_ckpt.jsonl and field_gaia.csv)
import json,io,sys,requests,numpy as np,pandas as pd
F=sys.argv[1]; D=sys.argv[2] if len(sys.argv)>2 else '.'
r=pd.DataFrame([json.loads(l) for l in open(f'{D}/search_ckpt.jsonl')]); r['Source']=r.Source.astype(str)
for c in ['P_zg','fap_zg','amp_zg','chi_zr','amp_zr','fap_zr','P_zr']:
    if c not in r: r[c]=np.nan
g=pd.read_csv(f'{D}/field_gaia.csv',dtype={'Source':str}); r=r.merge(g,on='Source',how='left')
ra0,de0=g.RA_ICRS.mean(),g.DE_ICRS.mean()
rad=np.max(np.hypot((g.RA_ICRS-ra0)*np.cos(np.radians(de0)),g.DE_ICRS-de0))+0.02
def same(a,b,tol): return np.isfinite(a)&np.isfinite(b)&np.any([abs(a*k/b-1)<tol for k in [1,2,.5]],axis=0)
al=np.zeros(len(r),bool)
for a in [1,0.5,1/3,1.0351,0.5087,29.53,14.77,1/(1+1/29.53)]: al|=abs(r.P_zr/a-1)<0.03
al|=abs(r.P_zr/0.2-1)<0.02
r['cand']=same(r.P_zr.values,r.P_zg.values,0.003)&(r.fap_zr<1e-6)&(r.fap_zg<1e-6)&~al&(r.P_zr<100)
gv=set(); status='gaia_vari: not reachable'
q=f"SELECT v.source_id FROM gaiadr3.vari_summary v JOIN gaiadr3.gaia_source g USING(source_id) WHERE 1=CONTAINS(POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',{ra0},{de0},{rad}))"
for i in range(3):
    try:
        x=requests.post("https://gea.esac.esa.int/tap-server/tap/sync",data={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":q},timeout=300); x.raise_for_status()
        gv=set(pd.read_csv(io.StringIO(x.text),dtype=str).iloc[:,0]); status=f'gaia_vari: {len(gv)} rows'; break
    except Exception as e: print('ERR',e)
r['in_gaia_vari']=r.Source.isin(gv); r['known']=r.in_gaia_vari|(r.VarFlag=='VARIABLE')
c=r[r.cand]; n=c[~c.known]
print(F,f'center {ra0:.3f} {de0:.3f} r {rad:.3f}',status); print('scored',len(r),'rule cands',len(c),'known(Gaia)',int(c.known.sum()),'not known to Gaia',len(n))
print(n[['Source','RA_ICRS','DE_ICRS','Gmag','P_zr','P_zg','amp_zr','fap_zr']].sort_values('P_zr').to_string())
n.to_csv(f'{D}/{F}_flags_unvetted.csv',index=False)
json.dump(dict(field=F,center=[ra0,de0,rad],scored=len(r),rule_cands=len(c),known_gaia=int(c.known.sum()),flags_not_in_gaia=len(n),catalog_status=status),open(f'{D}/{F}_flag_summary.json','w'),indent=1)
