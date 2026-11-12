# Gaia DR3 stars in a circular field via the IRSA TAP service.
# usage: RA=100.0 DEC=30.0 R=0.6 GR=12..15.5 python3 01_field.py   -> field_gaia.csv
import os,io,requests,pandas as pd
RA,DEC,R=[float(os.environ[k]) for k in ("RA","DEC","R")]; lo,hi=os.environ.get("GR","12..15.5").split("..")
q=f"SELECT source_id,ra,dec,phot_g_mean_mag,bp_rp,parallax,phot_variable_flag FROM gaia_dr3_source WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},{R}))=1 AND phot_g_mean_mag BETWEEN {lo} AND {hi}"
x=requests.post("https://irsa.ipac.caltech.edu/TAP/sync",data={"QUERY":q,"FORMAT":"csv"},timeout=900); x.raise_for_status()
d=pd.read_csv(io.StringIO(x.text),dtype={'source_id':str})
d.columns=['Source','RA_ICRS','DE_ICRS','Gmag','BP-RP','Plx','VarFlag']; d.to_csv('field_gaia.csv',index=False); print(len(d))
