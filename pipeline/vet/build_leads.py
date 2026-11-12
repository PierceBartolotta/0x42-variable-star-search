# Build the deduplicated list of leads to vet from one or more searched fields.
# usage: python3 build_leads.py FIELD:DIR [FIELD:DIR ...]
#   DIR holds search_ckpt.jsonl, field_gaia.csv, lc/ and FIELD_flags_unvetted.csv (from search/05_flag_rule.py).
# (Simplified public version of the script used for the Nov 2026 backlog, where the flag files sat in several run folders.)
import sys,os,pandas as pd
rows=[]
for arg in sys.argv[1:]:
    field,d=arg.split(':',1); d=os.path.abspath(d)
    x=pd.read_csv(f'{d}/{field}_flags_unvetted.csv',dtype={'Source':str})
    if 'in_gaia_vari' in x: x=x[x.in_gaia_vari!=True]
    for r in x.itertuples():
        rows.append(dict(source_id=str(r.Source),field=field,src_file=f'{d}/{field}_flags_unvetted.csv',ckpt=f'{d}/search_ckpt.jsonl',lcdir=f'{d}/lc',note=''))
L=pd.DataFrame(rows); print(L.groupby('field').size()); n0=len(L); L=L.drop_duplicates('source_id'); print('total',n0,'dedup',len(L))
L.to_csv('leads.csv',index=False)
