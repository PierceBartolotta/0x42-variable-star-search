# Fetch ZTF public data release light curves (IRSA) within 1.5" of each Gaia star in field_gaia.csv -> lc/<source_id>.csv
import csv,os,requests,concurrent.futures as cf
os.makedirs('lc',exist_ok=True)
rows=list(csv.DictReader(open('field_gaia.csv')))
def get(r):
    f=f"lc/{r['Source']}.csv"
    if os.path.exists(f) and os.path.getsize(f)>0: return
    u=f"https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20{r['RA_ICRS']}%20{r['DE_ICRS']}%200.00042&BAD_CATFLAGS_MASK=32768&FORMAT=csv"
    for _ in range(3):
        try:
            x=requests.get(u,timeout=180); 
            if x.ok: open(f,'w').write(x.text); return
        except Exception: pass
with cf.ThreadPoolExecutor(6) as ex: list(ex.map(get,rows))
print(len(os.listdir('lc')))
