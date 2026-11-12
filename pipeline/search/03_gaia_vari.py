# Gaia DR3 known variables in the field (vari_summary + classifier) via the ESA Gaia TAP.
# The other catalogs (VSX, Chen+2020, ATLAS, ASAS-SN, SIMBAD, TESS EB) are checked per lead later in vet/vet_backlog.py.
# usage: RA=.. DEC=.. R=.. python3 03_gaia_vari.py
import os,io,requests,pandas as pd
from astropy.table import Table
RA,DEC,R=[float(os.environ[k]) for k in ("RA","DEC","R")]; R+=0.02
q=f"SELECT v.source_id,g.ra AS _RAJ2000,g.dec AS _DEJ2000,c.best_class_name FROM gaiadr3.vari_summary v JOIN gaiadr3.gaia_source g USING(source_id) LEFT JOIN gaiadr3.vari_classifier_result c USING(source_id) WHERE 1=CONTAINS(POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',{RA},{DEC},{R}))"
for i in range(3):
  try:
    x=requests.post("https://gea.esac.esa.int/tap-server/tap/sync",data={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":q},timeout=900); x.raise_for_status(); break
  except Exception as e: print('ERR',e)
d=pd.read_csv(io.StringIO(x.text)).rename(columns={'_raj2000':'_RAJ2000','_dej2000':'_DEJ2000'}); Table.from_pandas(d).write('cat_gaia_vari.ecsv',overwrite=True); print('gaia_vari',len(d))
