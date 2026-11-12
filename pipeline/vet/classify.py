# classify.py: verdicts from json/<id>.json (+ tess/<id>_tess.json if present). Rules fixed before looking at results.
import json,glob,os,numpy as np,pandas as pd
VARTYPES=('V*','EB*','RR*','dS*','Ce*','LP*','Mi*','RG*','Ro*','BY*','RS*','SB*','Pu*','Er*','Or*','WU*','Al*','bL*','El*','gD*','sr*','Ir*','CV*','No*','Psr','XB*','HXB','LXB','YSO','TT*','Ae*')
def harmonic(a,b,tol=0.003): return any(abs(a*k/b-1)<tol for k in (1,2,.5))
SY=pd.read_csv('systematic.csv',dtype={'source_id':str}).set_index('source_id')
def stab(h,P): return abs(h[0]/h[1]-1)<max(0.003,P/2700)  # half-baseline (~1350 d) frequency resolution for long P
PH=pd.read_csv('phase_check.csv',dtype={'source_id':str}).set_index('source_id')
rows=[]
for f in sorted(glob.glob('json/*.json')):
    o=json.load(open(f)); sid=o['source_id']; R=[]; W=[]; rej=[]; notrun=[]; I=[]
    r=dict(source_id=sid,field=o['field'],ra=round(o['ra'],6),dec=round(o['dec'],6),G=round(o['G'],2),bp_rp=o.get('bp_rp'))
    if 'lc_err' in o:
        r.update(verdict='WEAK',reasons='ZTF light curve could not be fetched: '+o['lc_err']); rows.append(r); continue
    zr,zg=o.get('zr'),o.get('zg'); P0=o['P_search']
    P=zr['P_near'] if zr else P0; r['period_d']=round(P,6)
    r['amp_r']=round(zr['amp_model'],3) if zr else None; r['amp_g']=round(zg['amp_model'],3) if zg else None
    # --- catalogs
    X=o['xmatch']; hits=[]
    for k,v in X.items():
        if isinstance(v,str): notrun.append(f'{k} ({v[:60]})'); continue
        if k=='VSX_live30':
            for nm,typ,per,sep in v:
                if sep<=5: hits.append(f'VSX {nm} {typ} at {sep}"')
                elif per and harmonic(float(per),P,0.01): hits.append(f'VSX {nm} {typ} P={per} at {sep}" (same period: leakage)')
                elif sep<=10: W.append(f'VSX {nm} ({typ}) at {sep}"')
                else: I.append(f'VSX {nm} ({typ}, P={per}) at {sep}" (different period)')
        elif k=='SIMBAD10':
            for nm,ot,sep in v:
                isvar=ot in VARTYPES or 'V*' in ot or ot.startswith(('EB','RR','dS','Ce','LP','Mi','Ro','BY','RS','WU','Pu','Ir','sr','El','Er','Or','CV','No','YSO','TT','Y*'))
                if isvar: (hits if sep<=5 else W).append(f'SIMBAD {nm} [{ot}] at {sep}"')
        else:
            for z in v:
                sep=float(z.get('_r',99)); 
                if k.startswith('GaiaDR3') and z.get('Source')==sid: hits.append(f'{k} self'); continue
                if sep<=5: hits.append(f'{k} at {sep}"')
                else:
                    per=z.get('Period') or z.get('Per') or z.get('fp-lngfitper') or z.get('P')
                    try: same=per and harmonic(float(per),P,0.01)
                    except: same=False
                    (hits if same else W).append(f'{k} star at {sep}"'+(' same period' if same else ''))
    if hits: rej.append('catalogued: '+'; '.join(hits[:3]))
    # --- period quality
    if o.get('alias'): rej.append('alias: '+o['alias'])
    if zr is None: rej.append('too few clean ZTF r points')
    else:
        if zr['fap_near']>1e-6: rej.append(f"r-band signal at P not significant now (FAP {zr['fap_near']:.1e})")
        if not stab(zr['halves'],P): rej.append(f"r halves disagree {zr['halves'][0]:.5f}/{zr['halves'][1]:.5f}")
        if not harmonic(zr['P_blind'],P): W.append(f"r blind top peak {zr['P_blind']:.5f} not at P")
        if zr['R_fold']<0.1: W.append(f"weak r fold (explains {100*zr['R_fold']:.0f}% of variance)")
    if zg is None: W.append('too few clean ZTF g points')
    else:
        if zg['fap_near']>1e-6: W.append(f"g band not significant at P (FAP {zg['fap_near']:.1e})")
        elif not harmonic(zg['P_blind'],P): W.append(f"g blind top peak {zg['P_blind']:.5f} differs from r")
        if not stab(zg['halves'],P): W.append(f"g halves disagree {zg['halves'][0]:.5f}/{zg['halves'][1]:.5f}")
    sy=SY.loc[sid]; n=sy.sys_n
    if n==n and n>=3 and sy.sys_p<1e-3: rej.append(f'systematic: {int(n)} other field stars share the period (expect {sy.sys_expected})')
    elif n==n and n>0 and sy.sys_p<0.05: W.append(f'{int(n)} other field star(s) share the period (expect {sy.sys_expected})')
    r['field_same_P']=None if n!=n else int(n)
    # --- blending
    nb=o['neigh15']
    if isinstance(nb,str): notrun.append('Gaia neighbours ('+nb[:50]+')')
    else:
        for s2,sep,g2 in nb:
            if g2 is None: continue
            if sep<=5 and g2<o['G']+2.5: rej.append(f'blend: G={g2} star at {sep}"')
            elif sep<=5 and r.get('amp_r') and r['amp_r']*10**(0.4*(g2-o['G']))<1.0: W.append(f'neighbour G={g2} at {sep}" could produce the signal if it varied by {r["amp_r"]*10**(0.4*(g2-o["G"])):.2f} mag')
            elif sep<=10 and g2<o['G']: W.append(f'brighter star G={g2} at {sep}"')
    ph=PH.loc[sid] if sid in PH.index else None
    if ph is not None and ph.h1_dphase_cyc==ph.h1_dphase_cyc:
        r['rg_phase_diff_cyc']=ph.h1_dphase_cyc
        if abs(ph.h1_dphase_cyc)>0.1 and abs(ph.h1_dphase_cyc)>3*ph.h1_err_cyc: W.append(f'r and g out of phase by {ph.h1_dphase_cyc:+.2f} cycle')
    # --- neighbour LC / ensemble test
    nf=f'nbr/{sid}_nbr.json'; r['nbr_test']='not run'
    if os.path.exists(nf):
        N=json.load(open(nf)); t=N.get('nbr_top_at_P',[]); e=N.get('ensemble_corrected')
        r['nbr_test']=f"{N['n_nbr']} nbrs, {len(t)} with top peak at P"+(f", ensemble-corrected amp {e['amp_after']:.3f}/{e['amp_before']:.3f}" if e else ', ensemble corr not possible')
        if len(t)>=2: rej.append(f'systematic: {len(t)} nearby stars have their top peak at P')
        elif len(t)==1: W.append(f'nearby star {t[0][0]} ({t[0][2]}\') also peaks at P')
        if e and (e['fap_after']>1e-6 or e['amp_after']<0.5*e['amp_before']): rej.append('signal mostly removed by neighbour-ensemble correction')
    # --- TESS
    tf=f"tess/{sid}_tess.json"; tess='not run'
    if os.path.exists(tf):
        T=json.load(open(tf)); res=[x for x in T['res'] if 'err' not in x]
        if not res: tess='no usable sectors'
        else:
            m=[x for x in res if x['match_top'] and x['fap_top']<1e-6]
            on=[x for x in m if x['centroid_off_px']<=1.5]
            if on: tess=f"CONFIRMED in {len(on)}/{len(res)} sectors (top peak = P, signal centred on target, S"+','.join(str(x['sector']) for x in on)+")"
            elif m: tess=f"period seen in {len(m)}/{len(res)} sectors but signal centred {min(x['centroid_off_px'] for x in m):.1f} px (~{21*min(x['centroid_off_px'] for x in m):.0f}\") off target"; W.append('TESS sees P but centred on another star')
            else:
                nearP=[x for x in res if x['fap_top']<1e-6 and x['P_top']<12.5 and any(abs(x['P_top']*k/P-1)<0.03 for k in (1,2,.5)) and x['centroid_off_px']>1.5]
                tess=f"not detected in {len(res)} sectors (amplitude likely below TESS limit at this G / crowding)"
                if nearP: tess+=f"; a nearby star ~{21*nearP[0]['centroid_off_px']:.0f}\" away shows P={nearP[0]['P_top']:.4f}"; W.append('TESS shows a similar period from a nearby star')
    r['tess']=tess
    r['catalog_checks_not_run']='; '.join(notrun)
    if rej: r['verdict']='REJECTED'; reasons=rej+W
    elif W or notrun: r['verdict']='WEAK'; reasons=W+(['catalog check not run: '+'; '.join(notrun)] if notrun else [])
    else: r['verdict']='STRONG'; reasons=['r and g agree, both halves stable, clean fold, no catalog match, no blend within 5", no field systematic']
    if o.get('note'): reasons.append(o['note'])
    reasons+=I[:2]
    r['reasons']=' | '.join(reasons)
    # rough class guess (guess only)
    ratio=(zg['amp_model']/zr['amp_model']) if zg and zr and zr['amp_model']>0 else np.nan
    if P<0.1: cg='DSCT/HADS:'
    elif P<0.25: cg='DSCT: or EW: (check 2P)'
    elif P<1.0: cg='RRC:/RRAB:' if ratio>1.3 else 'EW: (2P) or ELL:'
    elif P<20: cg='ROT:/ELL:/EA:' if not ratio>1.3 else 'CEP:/pulsator:'
    else: cg='SR:/LPV: or ROT:'
    r['class_guess']=cg; r['g_over_r_amp']=round(ratio,2) if ratio==ratio else None
    rows.append(r)
D=pd.DataFrame(rows)
for i,a in D.iterrows():
    if a.verdict=='REJECTED' or a.period_d!=a.period_d: continue
    for j,b in D.iterrows():
        if i!=j and a.field==b.field and b.period_d==b.period_d and any(abs(a.period_d*k/b.period_d-1)<0.002 for k in (1,2,.5)):
            D.loc[i,'reasons']=(D.loc[i,'reasons'] if a.verdict!='STRONG' else '')+f' | another lead in {a.field} ({b.source_id}) has nearly the same period'; 
            if a.verdict=='STRONG': D.loc[i,'verdict']='WEAK'; D.loc[i,'reasons']=D.loc[i,'reasons'].lstrip(' |')
            break
D.to_csv('backlog_vet_1110.csv',index=False)
print(D.verdict.value_counts()); print(D.groupby(['field','verdict']).size().unstack(fill_value=0))
