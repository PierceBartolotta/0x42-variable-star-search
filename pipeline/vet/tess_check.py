# tess_check.py: per-sector TESS FFI check (TESScut) plus a pixel amplitude map, to see whether the ZTF period shows up in TESS
# and whether the signal sits on the target pixel or on a neighbour. Gaia neighbours come from the IRSA TAP. Output in tess/.
# usage: python3 tess_check.py SOURCE RA DEC P_ztf [maxsec]
import sys,json,warnings,numpy as np,pandas as pd,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
from astroquery.mast import Tesscut;from astropy.coordinates import SkyCoord;from astropy.wcs import WCS
from astropy.timeseries import LombScargle;from scipy.ndimage import median_filter;
sid,ra,de,P=sys.argv[1],float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]);mx=int(sys.argv[5]) if len(sys.argv)>5 else 6
c0=SkyCoord(ra,de,unit='deg');N=15;OD='tess'
secs=sorted(int(s) for s in Tesscut.get_sectors(coordinates=c0)['sector'])
use=secs if len(secs)<=mx else [secs[int(round(i))] for i in np.linspace(0,len(secs)-1,mx)]
import io,requests
_q=f"SELECT source_id,ra,dec,phot_g_mean_mag FROM gaia_dr3_source WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{de},{8*21/3600}))=1 AND phot_g_mean_mag<18"
g=pd.read_csv(io.StringIO(requests.post('https://irsa.ipac.caltech.edu/TAP/sync',data={'QUERY':_q,'FORMAT':'csv'},timeout=300).text),dtype={'source_id':str})
g['sep']=SkyCoord(g.ra.values,g.dec.values,unit='deg').separation(c0).arcsec;g=g.sort_values('sep').reset_index(drop=True)
Gt=g.phot_g_mean_mag.iloc[0]
def detr(t,f):
    if P<2.5:
        dt=np.median(np.diff(t));w=int(max(1.0,4*P)/dt)|1;return f/median_filter(f,w,mode='nearest')-1
    out=np.empty_like(f);seg=np.r_[0,np.where(np.diff(t)>1.0)[0]+1,len(t)]
    for a,b in zip(seg[:-1],seg[1:]):
        cc=np.polyfit(t[a:b]-t[a],f[a:b],1) if b-a>10 else [0,np.median(f[a:b])];out[a:b]=f[a:b]/np.polyval(cc,t[a:b]-t[a])-1
    return out
def harmfit(t,y,P,nh=2):
    md=np.median(y);mad=1.4826*np.median(abs(y-md))+1e-12;k=abs(y-md)<6*mad;t,y=t[k],y[k]
    ph=2*np.pi*t/P;A=np.column_stack([np.ones_like(t)]+[fn(h*ph) for h in range(1,nh+1) for fn in (np.sin,np.cos)])
    cc=np.linalg.lstsq(A,y,rcond=None)[0];mod=A@cc;return mod.max()-mod.min()
res=[]
for sec in use:
    try: h=Tesscut.get_cutouts(coordinates=c0,size=N,sector=sec)[0]
    except Exception as e: res.append(dict(sector=sec,err=str(e)[:80]));continue
    d=h[1].data;q=d['QUALITY']==0;t=d['TIME'][q];F=d['FLUX'][q].astype(float);ok=np.isfinite(F).all(axis=(1,2))&np.isfinite(t);t,F=t[ok],F[ok]
    if len(t)<200: res.append(dict(sector=sec,err=f'only {len(t)} good frames'));continue
    w=WCS(h[2].header);x0,y0=[float(v) for v in w.world_to_pixel(c0)]
    fl=F.reshape(len(t),-1);bkg=np.percentile(fl,20,axis=1);F=F-bkg[:,None,None]
    iy,ix=int(round(y0)),int(round(x0));ap=F[:,iy-1:iy+2,ix-1:ix+2].sum(axis=(1,2))
    fr=detr(t,ap/np.nanmedian(ap))
    ls=LombScargle(t,fr);f,p=ls.autopower(minimum_frequency=1/13,maximum_frequency=min(50,0.5/np.median(np.diff(t))*0.98),samples_per_peak=20)
    i=p.argmax();Ptop=1/f[i];fap=float(ls.false_alarm_probability(p[i]))
    ff=np.linspace(0.98/P,1.02/P,401);pp=ls.power(ff);Pn=1/ff[pp.argmax()];fapn=float(ls.false_alarm_probability(pp.max(),minimum_frequency=1/13,maximum_frequency=f.max()))
    rank=int((p>pp.max()).sum())+1
    amp_ap=harmfit(t,fr,P)
    A=np.zeros((N,N))
    for a in range(N):
        for b in range(N):
            y=F[:,a,b];med=np.nanmedian(y)
            if not np.isfinite(med) or med<=0: continue
            A[a,b]=harmfit(t,detr(t,y/med)*med,P)
    from scipy.ndimage import uniform_filter
    B=np.clip(A-np.median(A),0,None);yy,xx=np.indices(A.shape);S=uniform_filter(B,3);pk=np.unravel_index(S.argmax(),S.shape)
    m=(abs(yy-pk[0])<=1)&(abs(xx-pk[1])<=1)
    cy,cx=(B[m]*yy[m]).sum()/B[m].sum(),(B[m]*xx[m]).sum()/B[m].sum()
    cand=[]
    for r in g[g.phot_g_mean_mag<Gt+3].itertuples():
        sx,sy=[float(v) for v in w.world_to_pixel(SkyCoord(r.ra,r.dec,unit='deg'))]
        if not(1<sx<N-2 and 1<sy<N-2): continue
        M=np.exp(-((xx-sx)**2+(yy-sy)**2)/(2*0.9**2));cr=np.corrcoef(M.ravel(),B.ravel())[0,1]
        cand.append((str(r.source_id),round(float(r.phot_g_mean_mag),2),round(float(r.sep),1),round(float(np.hypot(sx-cx,sy-cy)),2),round(float(cr),3)))
    cand.sort(key=lambda z:-z[4])
    res.append(dict(sector=sec,n=len(t),P_top=float(Ptop),fap_top=fap,P_near=float(Pn),fap_near=fapn,rank_near=rank,amp_rel_ap=float(amp_ap),
        match_top=bool(any(abs(Ptop*k/P-1)<0.01 for k in [1,2,.5])),target_px=(round(y0,2),round(x0,2)),peak3x3_px=tuple(int(v) for v in pk),
        centroid_px=(round(float(cy),2),round(float(cx),2)),centroid_off_px=round(float(np.hypot(cy-y0,cx-x0)),2),peak_over_median=round(float(A.max()/np.median(A)),1),best_psf_match=cand[:3]))
    fig,ax=plt.subplots(1,2,figsize=(9,3.6));ph=(t/P)%1;ax[0].plot(ph,fr,'.',ms=1,alpha=.3);bb=np.floor(ph*25);mm=pd.Series(fr).groupby(bb).median();ax[0].plot((mm.index+.5)/25,mm.values,'ko-',ms=3)
    ax[0].set_title(f'S{sec} fold P={P:.5f} top={Ptop:.4f}',fontsize=8);ax[1].imshow(A,origin='lower');ax[1].plot(x0,y0,'r+',ms=14);ax[1].plot(cx,cy,'wx',ms=10)
    for r in g[g.phot_g_mean_mag<Gt+2].itertuples():
        sx,sy=w.world_to_pixel(SkyCoord(r.ra,r.dec,unit='deg'));ax[1].plot(sx,sy,'c.',ms=4)
    ax[1].set_title('amp map at P (+target, x centroid, . Gaia G<Gt+2)',fontsize=7);fig.tight_layout();fig.savefig(f'{OD}/{sid}_tess_S{sec}.png',dpi=80);plt.close()
    del F,d,h
out=dict(Source=sid,P=P,sectors_available=secs,sectors_used=use,res=res)
json.dump(out,open(f'{OD}/{sid}_tess.json','w'),indent=1,default=str)
for r in res: print(r)
