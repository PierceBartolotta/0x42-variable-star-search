# Pipeline

Python 3, packages in `../requirements.txt`. All data comes from free public services (IRSA, MAST, ESA Gaia archive, VizieR, SIMBAD, the VSX API). No accounts or keys are needed.

## 1. Search one field (`search/`)
Run inside an empty folder for the field:

```
RA=100.0 DEC=30.0 R=0.6 GR=12..15.5 python3 01_field.py     # Gaia DR3 stars -> field_gaia.csv
python3 02_ztf.py                                         # ZTF light curves -> lc/
RA=100.0 DEC=30.0 R=0.6 python3 03_gaia_vari.py            # Gaia DR3 known variables (optional, for reference)
python3 04_search.py                                      # periods in ZTF r and g -> search_ckpt.jsonl
python3 05_flag_rule.py dH .                              # flag rule, drop Gaia-known -> dH_flags_unvetted.csv
```

## 2. Vet the flags (`vet/`)
Run in a working folder:

```
python3 build_leads.py dH:/path/to/dH dI:/path/to/dI      # leads.csv
python3 vet_backlog.py                                    # per-star checks -> json/<id>.json
python3 systematic.py                                     # shared-period field test -> systematic.csv
python3 phase_check.py                                    # r/g phase test -> phase_check.csv
python3 classify.py                                       # first pass verdicts
python3 nbr_check.py <id>                                 # neighbour light curves, for leads still standing
python3 tess_check.py <id> <ra> <dec> <P> 3               # TESS check, for leads still standing
python3 classify.py                                       # final verdicts -> backlog_vet_1110.csv
python3 plot_folds.py <id> ...                            # fold plots -> plots/
```

These are the scripts as they were run, with only comments and file paths tidied up. They are research scripts, not a polished package.
