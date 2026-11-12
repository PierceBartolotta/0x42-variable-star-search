# r/g phase consistency: phase of the 1st harmonic (and 2nd) at the common P in ZTF r and g. Real variables are ~in phase.
import os,json,numpy as np,pandas as pd
L=pd.read_csv('leads.csv',dtype={'source_id':str}).set_index('source_id'); out=[]
def fit(t,y,P):
    ph=2*np.pi*t/P; A=np.column_stack([np.ones_like(t),np.sin(ph),np.cos(ph),np.sin(2*ph),np.cos(2*ph)]); c,res,_,_=np.linalg.lstsq(A,y,rcond=None)
    s=np.sqrt(np.sum((y-A@c)**2)/(len(y)-5)); cov=s**2*np.linalg.inv(A.T@A)
    return c,cov
for sid in L.index:
    f=f'json/{sid}.json'; o=json.load(open(f))
    if 'lc_err' in o or not o.get('zr') or not o.get('zg'): out.append(dict(source_id=sid)); continue
    P=o['zr']['P_near']; p=f'lc/{sid}.csv' if os.path.exists(f'lc/{sid}.csv') else f"{L.loc[sid,'lcdir']}/{sid}.csv"; df=pd.read_csv(p); R={}
    for fc in ['zr','zg']:
        d=df[(df.filtercode==fc)&(df.catflags==0)&(df.magerr<0.1)].copy(); d['m']=d.mag-d.groupby(['field','ccdid','qid']).mag.transform('median')
        R[fc]=fit(d.hjd.values-2458000,d.m.values,P)
    row=dict(source_id=sid)
    for h,(i,j) in {1:(1,2),2:(3,4)}.items():
        ph={fc:np.arctan2(R[fc][0][j],R[fc][0][i]) for fc in R}; A={fc:np.hypot(R[fc][0][i],R[fc][0][j]) for fc in R}
        e={fc:np.sqrt(R[fc][1][i,i]+R[fc][1][j,j])/max(A[fc],1e-9) for fc in R}   # approx phase error (rad)
        dphi=(ph['zr']-ph['zg']+np.pi)%(2*np.pi)-np.pi
        row[f'h{h}_dphase_cyc']=round(float(dphi/(2*np.pi)),3); row[f'h{h}_err_cyc']=round(float(np.hypot(e['zr'],e['zg'])/(2*np.pi)),3)
        row[f'h{h}_Ar']=round(float(A['zr']),4); row[f'h{h}_Ag']=round(float(A['zg']),4)
    out.append(row)
D=pd.DataFrame(out); D.to_csv('phase_check.csv',index=False); print(D.to_string())
