# Looking for new variable stars in free public data

This is a hobby research project. I search free public sky survey data for stars that change in brightness but are not in any variable star catalog yet. Everything here is a **candidate**. Nothing has been confirmed by a spectrum, by another observer, or by the AAVSO VSX moderators.

## What it does

1. Pick a small patch of sky and get every Gaia DR3 star in it brighter than about G 16.
2. Download each star's ZTF light curve (Zwicky Transient Facility, public data releases, 2018 to 2025) in the r and g bands.
3. Look for a repeating signal with a Lomb-Scargle periodogram, separately in r and in g.
4. Flag a star if r and g give the same period, the signal is strong in both, and the period is not a known observing alias (1 day, the Moon's cycle, a year).
5. Throw out stars Gaia DR3 already lists as variable.
6. Run every flag that is left through the vetting checks below.

## Vetting checks

A flag only counts as a **strong candidate** if it passes all of these:

- **Period stability.** The light curve is split into a first half and a second half. Both halves must give the same period.
- **r and g agree.** The period must show up independently in both ZTF bands, and the two bands must vary in phase with each other.
- **Field-wide shared-period test.** If many other stars in the same field show the same period, it is probably an instrument or data artifact, not the star. I count how many field stars share the period and compare that with what chance would give. I also pull the light curves of up to 8 similar stars nearby and check that none of them shows the signal.
- **Blend test.** Gaia DR3 neighbours within 15 arcsec are checked. A lead is rejected if a brighter star sits within 5 arcsec, and marked weak if a fainter neighbour could produce the signal with a reasonable amount of variation.
- **Catalog cross-match** within 10 arcsec (30 arcsec for the live VSX search): AAVSO VSX, ZTF periodic variables (Chen et al. 2020), ATLAS variables (Heinze et al. 2018), ASAS-SN variables, Gaia DR3 variability tables, TESS eclipsing binaries (Prsa et al. 2022) and SIMBAD. Any match rejects the lead.
- **TESS independent check.** For strong leads, I cut out TESS full-frame images (TESScut, up to 3 sectors), look for the same period, and map which pixels carry the signal, so a nearby star can't be mistaken for the target.

## Results (November 2026 backlog)

108 uncatalogued leads from 8 fields were vetted: **22 strong candidates, 36 weak, 50 rejected.** Most rejections were stars that turned out to be in an existing catalog already (mostly ATLAS).

- **Tier A (9):** ZTF r and g agree, and TESS independently shows the same period centred on the target star.
- **Tier B (13):** passes all ZTF checks, but TESS does not detect it. For stars this faint with variations of only 0.01 to 0.05 mag, that is expected and is not evidence against them, but it does mean they rest on one survey.

| tier | field | Gaia DR3 | RA (deg) | Dec (deg) | G | P (d) | ZTF amp r / g (mag) | type guess | TESS |
|---|---|---|---|---|---|---|---|---|---|
| A | ssD | 4320783833335254528 | 290.47610 | +15.18215 | 14.73 | 0.123508 | 0.027 / 0.040 | DSCT: or EW: (check 2P) | seen in 3/3 sectors |
| A | ssF | 511010986934527872 | 19.65206 | +62.45938 | 14.95 | 0.770712 | 0.036 / 0.026 | EW: (2P) or ELL: | seen in 1/3 sectors |
| A | ssF | 511040742465502848 | 20.56147 | +62.33228 | 15.19 | 0.182418 | 0.027 / 0.026 | DSCT: or EW: (check 2P) | seen in 2/3 sectors |
| A | ssF | 511044281518550528 | 20.96206 | +62.34744 | 14.63 | 0.136443 | 0.031 / 0.048 | DSCT: or EW: (check 2P) | seen in 2/3 sectors |
| A | ssF | 522603172387615872 | 18.76512 | +61.85266 | 13.49 | 0.860790 | 0.041 / 0.047 | EW: (2P) or ELL: | seen in 3/3 sectors |
| A | ssG | 3167111877664739840 | 110.51897 | +15.20047 | 13.82 | 0.034545 | 0.020 / 0.031 | DSCT/HADS: | seen in 1/3 sectors |
| A | ssH | 365394175204442112 | 10.10199 | +34.70653 | 14.20 | 1.488666 (EB at 2P = 2.977332) | 0.059 / 0.041 | ROT:/ELL:/EA: | seen in 3/3 sectors |
| A | ssH | 365668778233507712 | 9.01758 | +35.54700 | 14.82 | 0.047268 | 0.013 / 0.016 | DSCT/HADS: | seen in 2/3 sectors |
| A | ssI | 1976511380901555328 | 330.10919 | +49.49843 | 13.76 | 0.369475 | 0.019 / 0.028 | RRC:/RRAB: | seen in 2/3 sectors |
| B | dH | 3434889688272647424 | 99.57301 | +30.07908 | 15.24 | 12.416612 | 0.027 / 0.022 | ROT:/ELL:/EA: | not detected (3 sectors) |
| B | dI | 320414047865173504 | 20.60853 | +35.05867 | 13.51 | 12.681611 | 0.014 / 0.015 | ROT:/ELL:/EA: | not detected (3 sectors) |
| B | dI | 320628109035231104 | 20.38252 | +35.43081 | 15.49 | 20.340019 | 0.019 / 0.023 | SR:/LPV: or ROT: | not detected (3 sectors) |
| B | ssD | 4319934598041831296 | 290.32041 | +14.96850 | 15.46 | 8.672827 | 0.032 / 0.028 | ROT:/ELL:/EA: | not detected (3 sectors) |
| B | ssD | 4320705080851462144 | 290.13454 | +15.22241 | 15.01 | 2.945021 | 0.026 / 0.029 | ROT:/ELL:/EA: | not detected (3 sectors) |
| B | ssD | 4320725456180455936 | 289.48228 | +15.10909 | 15.43 | 6.418089 | 0.039 / 0.029 | ROT:/ELL:/EA: | not detected (3 sectors) |
| B | ssD | 4320807541555438208 | 290.19945 | +15.50189 | 15.96 | 0.593345 | 0.018 / 0.031 | RRC:/RRAB: | not detected (3 sectors) |
| B | ssF | 510568880183507200 | 19.34866 | +61.51697 | 14.57 | 3.119449 | 0.017 / 0.013 | ROT:/ELL:/EA: | not detected (3 sectors) |
| B | ssF | 510993807065560064 | 20.34385 | +62.27850 | 15.97 | 2.511374 | 0.046 / 0.046 | ROT:/ELL:/EA: | not detected (3 sectors) |
| B | ssF | 511088777381358336 | 20.44652 | +62.38944 | 13.89 | 0.690047 | 0.014 / 0.015 | EW: (2P) or ELL: | not detected (3 sectors) |
| B | ssG | 3167013127776175616 | 110.33588 | +14.98135 | 14.59 | 1.268475 | 0.024 / 0.019 | ROT:/ELL:/EA: | not detected (3 sectors) |
| B | ssG | 3167893355555020032 | 110.42196 | +15.60812 | 15.51 | 15.977808 | 0.018 / 0.027 | CEP:/pulsator: | not detected (2 sectors) |
| B | ssI | 1979541394428535040 | 329.69246 | +49.97195 | 15.45 | 4.816824 | 0.016 / 0.024 | CEP:/pulsator: | not detected (2 sectors) |

Amplitudes are peak-to-peak of a 2-harmonic fit. Type guesses use VSX style codes with a colon, which means uncertain. Several short-period ones could be contact binaries at twice the listed period.

Files:
- `data/strong_candidates.csv`: the 22 strong candidates with every check result.
- `data/all_108_leads_vetted.csv`: all 108 leads, including the weak and rejected ones and why.
- `plots/fold/<id>_fold.png`: ZTF r and g light curves folded at P and 2P.
- `plots/tess/<id>_tess_S<sector>.png`: TESS fold and pixel amplitude map per sector.
- `figures/365394175204442112_ztf_tess_fold.png`: the best candidate (below), ZTF and TESS folded at 2.977332 d.
- `pipeline/`: the scripts. See `pipeline/README.md`.

### Best candidate: Gaia DR3 365394175204442112

G 14.2, RA 10.10199, Dec +34.70653. It looks like an eclipsing binary with a period of 2.977332 d, twice the 1.488666 d the periodogram finds, with two minima per cycle. TESS shows the same signal on this star in all 3 sectors it was observed (17, 57, 84). The two minima are nearly equal in depth. TESS hints that alternate minima differ by about 0.003 mag, but that is not certain, so a 1.49 d period with only one visible eclipse is not ruled out.

## VSX submissions

VSX accepts one pending submission per person at a time.

- **Gaia DR3 2834068121983089152** (candidate delta Scuti / SX Phe type, P 0.0466004 d): submitted 25 Oct 2026, still waiting for a moderator (status "New") as of the last check in November 2026.
- **Gaia DR3 1049612805106097408** (candidate spotted rotating star, P 2.783 d): saved as a draft, not submitted yet.

These two came from an earlier round of this project and are not in the table above.

## Dormant black hole search

I also searched the Gaia DR3 two-body orbit catalog (all 443,205 orbital solutions) for stars orbiting something heavy and dark, the way Gaia BH1 and BH2 were found. The method does recover Gaia BH1 and BH2.

**It found no black hole candidates.** A handful of binaries came out with a large unseen companion on paper, but each one has an ordinary explanation that I can't rule out: a hidden pair of normal stars, a wrong mass for the visible star, or a doubtful Gaia orbit fit. They are just stars to look at again when Gaia DR4 comes out (expected December 2026), and I'm not listing them here as black hole candidates.

## Limitations

- Every signal comes from one survey. ZTF r and g are two filters on the same telescope, not two independent surveys. Only the 9 tier A stars have an independent check from TESS.
- Most amplitudes are tiny, 0.01 to 0.05 mag, close to what ZTF can measure at these magnitudes.
- No spectra. Type guesses are from light curve shape and colour only.
- Some checks were done with backup services because the main ones were down that day: Gaia variability tables via VizieR instead of the ESA archive, SIMBAD via the Harvard mirror. A star could still be in a catalog I didn't check.
- Some thresholds (the field-wide shared-period test, the half-split tolerance, the r/g phase test) were set after looking at a few leads, though before the final classification.
- ATLAS and ASAS-SN light curves need accounts and were not used.
- One TESScut product is labelled sector "1751" in the plots. That is the label the service returned; I haven't checked what it is.
- The ssI field was only partly searched (4000 of 8405 stars).
- I'm an amateur. If any of these stars is already published, the credit belongs to the earlier work, and I'd like to hear about it.

## Data and credits

ZTF public data releases (IRSA, Caltech/IPAC), Gaia DR3 (ESA), TESS full-frame images via TESScut (MAST), AAVSO VSX, VizieR, SIMBAD. This work uses data from these free public archives. Thank you to the teams that run them.

## License

Code is MIT licensed (see `LICENSE`). Tables and plots may be reused with credit.
