# Phase-fold PNGs (ZTF r and g, P and 2P) for given ids
import sys,os,json,numpy as np,pandas as pd,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
L=pd.read_csv('leads.csv',dtype={'source_id':str}).set_index('source_id')
for sid in sys.argv[1:]:
    o=json.load(open(f'json/{sid}.json')); P=o['zr']['P_near']
    p=f'lc/{sid}.csv' if os.path.exists(f'lc/{sid}.csv') else f"{L.loc[sid,'lcdir']}/{sid}.csv"; df=pd.read_csv(p)
    fig,ax=plt.subplots(2,2,figsize=(10,6),sharex=True)
    for i,fc in enumerate(['zr','zg']):
        d=df[(df.filtercode==fc)&(df.catflags==0)&(df.magerr<0.1)].copy()
        if len(d)<30: continue
        d['m']=d.mag-d.groupby(['field','ccdid','qid']).mag.transform('median')+d.mag.median()
        for j,PP in enumerate([P,2*P]):
            ph=((d.hjd-2458000.0)/PP)%1; a=ax[i,j]; a.plot(ph,d.m,'.',ms=2,alpha=.3,color='r' if fc=='zr' else 'g')
            b=np.floor(ph*30); m=d.m.groupby(b).median(); a.plot((m.index+.5)/30,m.values,'ko-',ms=3)
            lo,hi=np.percentile(d.m,[1,99]); a.set_ylim(hi+0.02,lo-0.02); a.set_title(f'ZTF {fc} n={len(d)} fold P={PP:.6f} d',fontsize=9)
    ax[1,0].set_xlabel('phase'); ax[1,1].set_xlabel('phase')
    fig.suptitle(f"Gaia DR3 {sid} ({o['field']}) G={o['G']:.2f}  CANDIDATE ONLY, not confirmed",fontsize=10)
    fig.tight_layout(); fig.savefig(f'plots/{sid}_fold.png',dpi=90); plt.close()
    print('plot',sid)
