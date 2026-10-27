# Findings: convergence clubs and explanatory models

Convergence and explanatory-modelling half of the analysis of domestic electricity consumption per meter across England and Wales, 2010–2024 (35,672 LSOAs, 7,264 MSOAs, 318 LADs, May 2026 vintage). Code: `projects/python/energy_inequality/analysis/b_*.py`, run end to end with

```bash
/Users/jay/anaconda3/envs/envs563/bin/python -I -W ignore projects/python/energy_inequality/analysis/b_run_models.py
```

Outputs are written to `data/data_python/energy_inequality/results/models/` (`tables/`, `figures/` as PNG and SVG, `cache/` for pickled fits, `logs/`). Paths below are relative to that folder. Heavy fits are cached, so a re-run with a warm cache takes a few minutes; `--force` re-estimates everything. Citations use keys in `references.bib`; nine methods references were added for this half (Hodrick–Prescott, Newey–West, da Silva–Fotheringham, Kelejian–Prucha, Anselin et al., Halleck Vega–Elhorst, LeSage–Pace, McCullagh, Barro–Sala-i-Martin), each with its DOI checked against Crossref.

Throughout, "consumption" means mean domestic kWh per meter. Means at MSOA and LAD level are total kWh divided by total meters, summed from LSOAs. Every association reported here is ecological: it describes places, not households [@openshaw1984modifiable; @fotheringham1991modifiable].

## Headline findings

1. **No convergence; stratified persistence.** Phillips–Sul tests reject convergence of consumption per meter at every scale and in every window: LAD log *t* = −52 for 2010–24 and −31 for 2010–19, before COVID-19. LADs sort into 12 ordered clubs, with 21 at MSOA level, and relative paths that were nearly fixed by 2010. Relative dispersion rose. The *ordering* of clubs is robust (rank correlation 0.85–0.98 across windows and transformations); their exact boundaries are not.
2. **Clubs are geographical.** Club membership is strongly spatially clustered: Moran's *I* of the club number is 0.41–0.58, and same-club neighbour links run 1.2–1.8 times their random expectation. It is profiled by electric and off-gas heating, household size, affluence and detached housing.
3. **Composition explains most of the cross-section.** Dwelling size and electric heating dominate. Census, EPC and deprivation covariates explain about 80% of the between-area variance. Spatial dependence is error-type (robust LM-error ≫ LM-lag), so it reflects unmeasured smooth context rather than spillovers.
4. **Relationships vary strongly over space.** GWR and MGWR raise *R*² to about 0.95 and cut AICc by 2,500–2,700 against the best global model. MGWR places age at the national scale, off-gas heating and social renting at the regional scale, electric heating and EPC efficiency at local scales (about 250 MSOAs), and household size and detached share at very local but mostly non-significant scales.
5. **Changes over 2010–24 and through the shocks.** Consumption per meter fell by about 19% over 2010–12 → 2022–24, and by about 11% over 2019 → 2022. In both changes, conditional on level, electric-heated, large-dwelling, rented and heat-pump areas fell least, and new-build growth and density went with larger falls. MGWR puts the heat-pump association at a uniform national scale. The level effect and new-build dilution are very local. For the long change MGWR is the best-fitting model (AICc 13,398 against 13,512 for GWR and 14,806 for SDM).
6. **Electrification is visible.** Areas with more renewable (heat-pump) heating in 2021 drifted steadily upwards relative to others, by about 1.3% per SD over 2010–24. Electric-heated areas did not cut more than others in the 2022 crisis once 2019 is the baseline.

### Hypotheses

| Hypothesis | Verdict | Evidence |
|---|---|---|
| H3: log *t* rejects global convergence | **Supported** | *t* = −27 to −95 at all scales and windows |
| H3: club membership is spatially autocorrelated | **Supported** | Same-club link ratio 1.2–1.8; Moran's *I* of club number 0.41–0.58; all *p* ≤ 0.023 |
| H3: clubs are associated with heating fuel and dwelling type | **Supported** | Ordered-logit odds ratios per SD (MSOA): electric heating 6.4, off-gas 2.1, detached 2.2, household size 3.1 |
| H4: household size and detached share act at near-global bandwidths | **Rejected** | MGWR bandwidths 44 and 53 MSOAs, significant in only 12% and 27% of MSOAs |
| H4: heating-fuel and EPC effects are local and survive spatial dependence | **Partly supported** | Electric heating local (257) and significant everywhere; off-gas regional (1,389); EPC local (241), negative where significant (40%). The error-type dependence is absorbed by SEM/SDEM without changing these signs |
| H5 (models' angle): areas with high consumption or electric heating cut most in 2022 | **Partly supported** | Over 2019→22, higher-consuming LSOAs fell proportionally more (unconditional β = −0.028, *t* = −15.7). But *conditional* on their level, electric-heated LSOAs fell *less* (+1.8% per SD, SEM). In the event study, the electric-heating gap was +0.65% per SD wider in 2022 than in 2019. Higher-consumption areas cut more, but not because they are electrically heated |

H1, H2 and the descriptive parts of H5 are assessed in `findings_dynamics.md`.

## 1. Club convergence

### Methods

We test for convergence with the log *t* regression of @phillips2007transition and group units into convergence clubs with the clustering algorithm of @phillips2009economic, following the implementation logic of the R package *ConvergenceClubs* [@sichera2019convergenceclubs]. For each unit, log kWh per meter is first smoothed with a Hodrick–Prescott filter [@hodrick1997postwar]. We use λ = 400, the value Phillips and Sul use for annual data. Relative transition paths are then *h*<sub>*it*</sub> = *X*<sub>*it*</sub> / mean<sub>*i*</sub> *X*<sub>*it*</sub>. We regress log(*H*<sub>1</sub>/*H*<sub>*t*</sub>) − 2 log log *t* on log *t*, where *H*<sub>*t*</sub> is the cross-sectional variance of *h*<sub>*it*</sub>, after discarding the first third of the sample (*r* = 1/3). Convergence is rejected when the slope's *t*-statistic falls below −1.65. The standard error is Newey–West HAC with a Bartlett kernel [@newey1987simple].

Clustering proceeds in four steps. Units are ordered by their final filtered value. A core group is formed by maximising *t*<sub>*k*</sub> over the leading units. Each remaining unit is sieved in if its log *t* with the core exceeds *c*\* = 0, raising *c*\* until the whole club passes. The procedure then recurses on the remainder. Adjacent clubs are merged whenever their joint log *t* statistic exceeds −1.65.

The implementation (`b_phillips_sul.py`) is our own vectorised Python code. Group variances are computed from running sums, which makes the algorithm feasible for all 35,636 complete LSOAs. A Monte Carlo check on a Phillips–Sul-type data-generating process recovers the planted clubs (`tables/clubs_monte_carlo.csv`). Under true convergence the test never rejected (size 0 in 200 replications, N = 300 or 2,000, T = 10 or 15). Under divergence and three-club designs it always rejected. The adjusted Rand index between true and estimated clubs was 0.985–1.000, and three clubs were recovered in 94.5–100% of replications. The test is therefore conservative at our T, as Phillips and Sul's small-T results lead one to expect.

The LAD level is the main specification. Robustness runs cover MSOAs and LSOAs and two further windows: the pre-shock window 2010–2019 (T = 10) and 2015–2024, which starts after the 2014/15 DESNZ processing break. A further variant filters kWh in levels rather than logs. In that variant *h*<sub>*it*</sub> is a pure ratio, so common proportional shocks such as COVID-19 and the 2022 price fall cancel exactly. We complement the test with σ-convergence diagnostics (`tables/clubs_sigma_dispersion.csv`).

Clubs are a-spatial by construction. We therefore test whether membership is spatially clustered, using queen contiguity (islands attached to their nearest neighbour; MSOA and LAD contiguity derived from LSOA contiguity) and three statistics, each with 999 permutations: (i) the share of neighbour links joining units in the same club, a multi-category join-count statistic; (ii) BB join counts for each club; and (iii) Moran's *I* of the ordinal club number [@moran1950notes; @rey2022pysal].

Club membership is then profiled with ordered logit [@mccullagh1980regression] on standardised 2021 covariates. Higher categories are higher-consumption clubs. At MSOA level, standard errors are clustered by LAD. As a robustness check that drops the ordering assumption, a multinomial logit is fitted on three tiers of clubs (high, middle, low), with average marginal effects.

### Key numbers

| | LAD 2010–24 | LAD 2010–19 | LAD 2015–24 | MSOA 2010–24 | LSOA 2010–24 |
|---|---|---|---|---|---|
| Global log *t* (*b̂*, *t*) | −1.08, −52.2 | −1.33, −30.9 | −1.28, −26.9 | −1.10, −67.9 | −1.12, −94.6 |
| Initial → merged clubs | 20 → 12 | 29 → 18 | 25 → 15 | 22 → 21 | 17 → 9 |
| Divergent units | 5 | 10 | 9 | 2 | 2 |
| Same-club link share, obs./random | 1.33 (*p* = 0.001) | 1.65 | 1.23 (*p* = 0.023) | 1.61 | 1.42 |
| Moran's *I* of club number | 0.41 | 0.43 | 0.41 | 0.55 | 0.51 |

- **σ-divergence.** The cross-LAD coefficient of variation rose from 0.117 (2010) to 0.128 (2019) and 0.129 (2024). Across LSOAs it rose from 0.216 to 0.246.
- **Levels variant.** Using kWh in levels gives the same verdict: LAD log *t* = −54.9, 12 clubs. Adjusted Rand index with the log specification is 0.71 and the Spearman correlation of club rank is 0.98.
- **Divergent LADs, 2010–24:** Isles of Scilly, Uttlesford, Cotswold, Elmbridge and South Tyneside.
- **Club 1** (highest steady state) contains West Berkshire, Maldon, Sevenoaks, Mid Suffolk, Mole Valley, Tandridge and Waverley. **Clubs 11–12** contain Hartlepool, Middlesbrough, Redcar and Cleveland, North Tyneside, Sunderland, Gateshead and Blaenau Gwent.
- **Stability of club rank.** Club rank is highly stable across windows: Spearman correlation 0.94 for LAD full vs pre-shock, 0.96 at MSOA and 0.85 at LSOA. Exact club boundaries are not stable: adjusted Rand 0.28, 0.37 and 0.14 respectively.
- **Spatial clustering of individual clubs.** BB join counts exceed random expectation for 19 of 21 MSOA clubs and for all 9 LSOA clubs (median ratio 2.2 and 4.7). At LAD level only 6 of 12 clubs are individually significant.
- **Two LSOA clubs with distinctive composition.** The lowest LSOA club (club 9, n = 271) fell from 3,687 to 2,304 kWh per meter. It is dominated by flats (55%), district heating (9%) and social renting (30%), and its median meter count grew by 28% over 2010–24, so it is largely new-build and regenerated stock. The two divergent LSOAs grew their meter counts more than fourfold.
- **Ordered logit, MSOA, 2010–24 clubs** (odds ratio per SD; McFadden *R*² = 0.26):
  - Higher-consumption clubs: electric heating 6.39, household size 3.13, less deprived 2.85, detached 2.16, oil/LPG/solid fuel 2.14, rural 1.47.
  - Lower-consumption clubs: private renting 0.52, higher EPC efficiency 0.56, Wales 0.71.
  - LAD ordered logit gives the same ranking. Less deprived (5.54), detached (4.93) and electric heating (4.84) lead.
- **Multinomial tiers.** A 1-SD higher electric-heating share raises the probability of the high tier by 8.5 points and lowers that of the low tier by 22.5 points. For renewable heating, the AME on the low tier is positive (+9.0 points), and the AME on the high tier is also positive (+3.0).

### Figures and tables

- Transition paths: `figures/clubs_paths_{lad,msoa,lsoa}_{full,pre,post15,levels}.png|svg`.
- Club maps: `figures/clubs_map_{lad,msoa}_{full,pre,post15,levels}`, `figures/clubs_map_lsoa_full`.
- Profile: `figures/clubs_profile_ordered_logit`.
- Tables:
  - `tables/clubs_summary.csv` (every club's *n*, *b̂*, *t* and first/last-year kWh)
  - `clubs_merge_log.csv`
  - `clubs_membership_{lad,msoa}.csv`, `clubs_membership_lsoa.parquet`
  - `clubs_join_counts.csv`
  - `clubs_profile_ordered_logit.csv`, `clubs_profile_mnlogit_ame.csv`, `clubs_profile_fit.csv`
  - `clubs_agreement.csv`, `clubs_crosstab_lad_full_pre.csv`
  - `clubs_monte_carlo.csv`, `clubs_sigma_dispersion.csv`

### Results

There is no evidence of national convergence in domestic electricity consumption per meter. The log *t* test rejects convergence decisively at every scale and in every window. The rejection holds before COVID-19 and the price crisis (2010–19), after the processing break (2015–24), and when common proportional shocks are removed by working in levels. Relative dispersion rose over the period, so the absolute fall in consumption (from about 4,240 to 3,400 kWh per LAD meter) did not narrow relative gaps.

The units instead fall into many clubs. Each club converges internally to its own relative level, and the levels are ordered. With 12 LAD clubs and 21 MSOA clubs, the picture is a graded stratification of steady states, not two or three discrete regimes. Most transition paths are nearly flat (see `clubs_paths_lad_full`): relative positions were largely fixed by 2010. The top clubs drift slowly upwards and the bottom clubs slowly downwards, which is what drives the global rejection.

The ordering of clubs is robust across windows, scales and transformations, with rank correlations of 0.85–0.98. Where an individual unit sits relative to a club boundary is not robust. Club numbers should therefore be read as ordinal strata, not as sharply bounded groups.

Membership is strongly spatially clustered at every scale. The same-club link share is 1.2–1.8 times its random expectation, and Moran's *I* of the club number is 0.41–0.58, all with *p* ≤ 0.023. This confirms the caveat in `methods_refs.md`: the a-spatial clustering recovers a spatially organised structure. High-consumption clubs cover the affluent, detached and partly off-gas commuter belts of the South East and East of England, together with rural Wales and the South West. Low-consumption clubs cover Tyneside, Teesside, the South Wales valleys and inner London.

The profiles show what the clubs proxy:

- **Heating fuel.** Electric and off-gas heating carry the largest odds ratios, which is consistent with @roberts2015fuel.
- **Dwelling and household size.**
- **Affluence.** The deprivation percentile has an odds ratio of 2.8–5.5.

Private renting and higher EPC efficiency mark the lower clubs. Areas with more renewable (heat-pump) heating are over-represented at *both* tails of the club distribution. This fits a technology that adds load where it replaces gas [@love2017addition], while new, efficient estates with heat pumps sit in low-consumption clubs.

### Discussion and limitations

H3 is **supported**. The Phillips–Sul test rejects global convergence, club membership is spatially autocorrelated, and membership is associated with heating fuel (electric and off-gas) and dwelling type. This extends to the neighbourhood scale the club evidence that @herrerias2017residential and @tastan2023club found for regions and cities. The number of clubs, however, should not be over-interpreted, for four reasons:

1. **Short and broken series.** T is short (10–15), and the first third is discarded, leaving 8–11 points in the regression. The series also contains the 2014/15 processing break.
2. **Many units.** With very large N, the sieve will find many fine strata. This is why the LSOA pre-shock run splits into 29 clubs.
3. **Smoothing over shocks.** The HP filter smooths over the 2020 and 2022 shocks rather than removing them.
4. **Inherited 2010–14 values.** Split LSOAs inherit their parent's 2010–14 values, which weakens LSOA-level paths in the early years.

The robust findings are the rejection of convergence, the stability of the club *ordering*, and the spatial concentration of clubs. Because clubs are defined on per-meter means, compositional change shows up as club movement: new-build flats with district heating, and the growth in meters, appear as the falling LSOA club 9. Club profiles are cross-sectional associations with 2021 characteristics, not explanations of the dynamics.

## 2. Cross-sectional explanatory models

### Methods

**Outcomes.** The main outcome is the log of the meter-weighted mean of 2022–2024, the recent post-shock level. The comparison outcome is 2019, the last pre-COVID year.

**Covariates.** All are standardised (*z*-scores), so each coefficient is the difference in log kWh per meter associated with a 1-SD difference in the covariate:

- Census 2021 central-heating shares, with mains gas as the reference: electric; oil/LPG/solid fuel ("off-gas"); renewable (including heat pumps); district or communal.
- Dwelling type, with semi-detached and terraced as the reference: detached share, flats share.
- Mean rooms and mean household size.
- Tenure, with owner-occupation as the reference: social rent, private rent.
- Age shares: 0–14 and 65+.
- The IMD/WIMD 2019 within-country rank percentile (higher = less deprived). The 2019 index is used because it predates both outcomes.
- Median EPC energy-efficiency score of existing dwellings (MSOA, FYE2025).
- Log population density, a rural dummy (RUC 2021) and a Wales dummy.

**Screening.** Starting from 18 candidates, VIF screening at 5 removed only the one-person-household share, leaving 17 covariates with a maximum VIF of 4.8 (`tables/xs_vif_final.csv`). LSOAs with fewer than 50 meters in any year of the window are excluded, leaving 35,662 for 2022–24 and 35,670 for 2019.

**LSOA models.** We estimate:

- OLS with White standard errors and Lagrange-multiplier diagnostics, including the robust LM tests [@anselin1996simple];
- SLX [@hallekvega2015slx];
- the spatial lag (SAR), spatial error (SEM), spatial Durbin (SDM) and spatial Durbin error (SDEM) models.

All use row-standardised queen weights. At n = 35,662, spreg's maximum-likelihood estimators need a dense (I − ρ*W*)⁻¹ for their variance, which is infeasible here. The LSOA spatial models are therefore estimated by GMM/IV: spatial two-stage least squares with White errors for the lag models, and the Kelejian–Prucha heteroskedasticity-robust GMM for the error models [@kelejian2010specification]. For SAR and SDM we report LeSage–Pace average direct, indirect and total impacts [@lesage2009introduction]. A sensitivity run drops LSOAs with any outlier flag.

**MSOA design for local models.** Local models are fitted at MSOA level (n = 7,264) on the whole national surface. The choice is principled:

1. The EPC covariate is native to MSOAs. LSOA-level local models would replicate one EPC value across about five LSOAs.
2. MSOA means rest on about 3,600 meters each, so they are less noisy.
3. A full multiscale GWR (MGWR) at n = 7,264 is computationally feasible without sampling.

The MSOA models use 13 continuous covariates (`b_models.XS_LOCAL`), the same set for the global and the local models so that AICc compares like with like. The rural and Wales dummies are excluded because they are constant inside most local kernels and would make local designs singular. District heating (near-zero in most kernels) and the child share (overlapping with household size at MSOA level) are also dropped.

At MSOA level we fit the ML global suite (OLS, SLX, SAR, SEM, SDM, SDEM), GWR [@brunsdon1996geographically] and MGWR [@fotheringham2017multiscale; @oshan2019mgwr]. Both local models use adaptive bisquare kernels with AICc-selected bandwidths, on standardised *y* and *X*.

MGWR is fitted for the main 2022–24 outcome; the 2019 outcome gets GWR only. Backfitting starts from the GWR bandwidth, with score-of-change tolerance 10⁻⁵, a minimum bandwidth of 30 MSOAs, and at most 200 iterations.

**MGWR implementation.** The mgwr package's inference step builds each kernel row in a Python loop and replays every backfitting iteration. At n = 7,264 with 14 terms, a timing test put it at more than 100 hours on the study machine. We therefore wrote `b_mgwr_fast.py`, which:

- calls mgwr's own backfitting (`multi_bw`) and golden-section bandwidth search unchanged, but computes the kernels and local fits for all locations at once;
- runs mgwr's inference recursion as dense matrix products, stopping only once further sweeps can change nothing;
- passes the result to mgwr's own `MGWRResults` for every diagnostic.

On a 500-MSOA subsample it reproduces the mgwr package exactly (`tables/mgwr_fast_validation.csv`): identical bandwidths and number of backfitting iterations, coefficients and standard errors equal to within 4 × 10⁻¹³, and the same AICc.

Local inference uses the corrected critical *t* of @dasilva2016multiple, which is based on each covariate's effective number of parameters; local coefficients that do not pass are shown in grey. The models are compared by AICc, computed in the same way for every model (−2 log *L* + 2*K* + 2*K*(*K* + 1)/(*n* − *K* − 1), with *K* = ENP + 1 for local models), and by Moran's *I* of the residuals. This follows the good-practice route of @comber2023route.

### Key numbers (LSOA, 2022–24)

- **OLS fit.** *R*² = 0.787 (2019: 0.773). Residual Moran's *I* = 0.45 (*z* = 141).
- **LM diagnostics.** LM-error 19,765 and LM-lag 6,441. Robust LM-error 13,326 (*p* < 0.001) against robust LM-lag 2.1 (*p* = 0.15). The tests therefore point to spatial error dependence, not a substantive spillover in *y*.
- **SEM (GMM).** λ = 0.72. Residual Moran's *I* falls to −0.06.
- **SAR.** Leaves residual *I* at 0.43.
- **SDM (S2SLS).** ρ = 0.98, near the boundary, and residual *I* = −0.20. With *WX* already in the model, the remaining instruments (*W*²*X*) are weak, so we do not interpret the LSOA SDM. SDEM (λ = 0.67) is the preferred Durbin-type specification.
- **SEM coefficients** (log points per SD; 2019 in brackets):

| Covariate | 2022–24 (2019) |
|---|---|
| Mean rooms | +0.117 (+0.113) |
| Electric heating | +0.113 (+0.116) |
| Off-gas heating | +0.053 (+0.052) |
| Household size | +0.038 (+0.036) |
| Renewable heating | +0.013 (+0.010) |
| Social rent | +0.010 |
| Aged 65+ | +0.010 (+0.017) |
| Detached | +0.008 (+0.003) |
| Rural | +0.007 |
| EPC efficiency | −0.004 |
| Less deprived | −0.005 |
| District heating | −0.022 (−0.026) |
| Wales | −0.025 |
| Log density | −0.040 |

- **Deprivation reverses sign between models.** In OLS, less deprived LSOAs consume more (+0.018). The SDEM separates the levels: the LSOA's own deprivation percentile has −0.009, and its neighbours' (*W*·dep) has +0.026. Within an affluent area, the more deprived LSOAs consume *more* per meter, while affluent surroundings raise consumption overall.
- **Robustness.** Dropping flagged outliers (159 LSOAs for 2022–24) changes no OLS coefficient by more than 0.008.

### Key numbers (MSOA, 2022–24; standardised *y* and *X*; 13 covariates; maximum VIF 7.4 for off-gas heating)

**Global diagnostics.** OLS *R*² = 0.80, residual Moran's *I* = 0.52. Robust LM-error = 3,803 against robust LM-lag = 1.2 (*p* = 0.27).

**Model comparison** (`tables/xs_msoa_y2224_model_comparison.csv`):

| Model | AICc | ΔAICc | Residual Moran's *I* | Pseudo-*R*² |
|---|---|---|---|---|
| OLS | 8,980 | 7,681 | 0.52 | 0.80 |
| SLX | 7,812 | 6,514 | 0.51 | 0.83 |
| SAR (ρ = 0.27) | 7,899 | 6,601 | 0.32 | 0.83 |
| SEM (λ = 0.81) | 4,463 | 3,164 | −0.08 | 0.74 |
| SDEM | 4,189 | 2,890 | −0.05 | 0.81 |
| SDM (ρ = 0.71) | 4,014 | 2,715 | −0.05 | 0.91 |
| **GWR** (bandwidth 178; ENP 1,350) | **1,299** | **0** | 0.05 | 0.956 |
| MGWR (ENP 1,321) | 1,492 | 193 | −0.03 | 0.954 |

The 2019 outcome gives the same ordering: GWR AICc 1,452 against SDM 4,280 and OLS 9,189.

**MGWR bandwidths** (nearest MSOAs out of 7,264). The backfitting converged at iteration 137 (score of change 9.9 × 10⁻⁶).

| Scale | Bandwidth | Local coefficients (mean; range) | Share significant after correction |
|---|---|---|---|
| Intercept | 41 | | 51% |
| **Very local** | | | |
| Household size | 44 | mean +0.08 | 12% (8% positive, 4% negative) |
| Detached share | 53 | mean +0.20, range −0.22 to +0.91 | 27% |
| **Local** | | | |
| Mean rooms | 168 | 0.19 to 1.28 | 99% |
| EPC efficiency | 241 | −0.22 to +0.04 | 40% (all negative) |
| Electric heating | 257 | +0.31 to +0.64 | 100% |
| Private rent | 283 | | 29% (all positive) |
| Renewable heating | 305 | −0.61 to +0.18 | 17% (15% negative) |
| **Regional** | | | |
| Deprivation | 505 | | 23%; 10% positive, 13% negative |
| Flats | 722 | | 25% |
| Social rent | 981 | | 100% |
| Off-gas heating | 1,389 | 0.27 to 0.56 | 100% |
| **National** | | | |
| Log density | 3,653 | | 100% |
| Aged 65+ | 7,263 (global) | +0.044 | 100% |

The corrected critical |*t*| ranges from 1.99 for age 65+ to 3.84 for the intercept.

### Figures and tables

- Figures: `figures/xs_lsoa_{y2224,y2019}_forest`, `figures/xs_mgwr_msoa_y2224_bandwidths`, `figures/xs_mgwr_msoa_y2224_map_<term>` (one map per covariate; grey = not significant after correction).
- Tables:
  - `tables/xs_vif_steps.csv`, `xs_vif_final.csv`, `xs_vif_msoa_local_set.csv`
  - `xs_{lsoa,msoa}_{y2224,y2019}_{coefs,fit,diag,impacts}.csv`
  - `xs_lsoa_outlier_sensitivity.csv`
  - `xs_msoa_{y2224,y2019}_model_comparison.csv`
  - `xs_{gwr,mgwr}_msoa_*_local_summary.csv`, `xs_mgwr_msoa_y2224_local.parquet`

### Results

**What the global models explain.** Small-area consumption per meter is highly predictable from housing and heating composition. Census, EPC and deprivation covariates explain about 79% of the variance between LSOAs and 80% between MSOAs. The ranking is stable between 2019 and 2022–24:

- **Dwelling size** (mean rooms) and **electric central heating** are the two dominant correlates, each worth about 11–12% more consumption per meter per SD.
- **Off-gas heating** (oil, LPG, solid fuel) follows at about 5%. In those homes, electricity drives pumps and supplementary heaters.
- **Household size** comes next at about 4%.
- **Density, district heating and Wales** carry negative associations. In Wales the effect is about 2.5% lower consumption per meter, conditional on composition.
- **EPC efficiency** has the expected negative sign but a small one, about −0.4% per SD at LSOA level and about −6% of an SD of *y* at MSOA level. This is consistent with modelled efficiency being a weak and partly mismeasured proxy for realised demand [@few2023overprediction; @summerfield2019empirical].

**Spatial dependence is a nuisance, not a spillover.** The residuals are strongly spatially autocorrelated, and the robust LM tests favour error dependence over a spatial lag at both scales. The spatial error and Durbin-type models remove almost all of the residual autocorrelation and greatly improve fit. The SAR does neither. We read this as unmeasured, spatially smooth determinants: climate, housing vintage, network and tariff history. It is not evidence that a neighbourhood's consumption causally raises its neighbours'.

**Deprivation changes sign once context is separated.** With SDEM, the LSOA's own deprivation percentile has a small negative association: within its surroundings, the less deprived LSOA consumes slightly less per meter. Living in a less deprived *neighbourhood* has a positive one. The positive raw association between affluence and consumption is therefore a contextual, regional gradient, not a within-neighbourhood one. This is a textbook case of why such coefficients must not be read at household level [@openshaw1984modifiable].

**Local models.**

- **Fit.** Local models fit far better than any global specification. GWR and MGWR reach *R*² ≈ 0.95 and cut AICc by 2,500–2,700 relative to the best global model (SDM). That is evidence of strong non-stationarity in the relationships, beyond smooth residual dependence.
- **GWR against MGWR.** GWR, with one bandwidth of 178 neighbours, has a lower AICc than MGWR (ΔAICc = 193), though MGWR leaves less residual autocorrelation. On information criteria the extra flexibility of covariate-specific bandwidths is therefore not rewarded here. The value of MGWR lies in what its bandwidths say about the scale of each process.
- **Electric heating** is positive and significant everywhere: +0.31 to +0.64 SD of *y* per SD. Its bandwidth is local (257 MSOAs, about 3.5% of the country). The effect is weakest in Surrey and around Bath (about 0.32–0.34 SD), and strongest in inner London (about 0.63–0.64 in Westminster, Southwark and the City), where electrically heated flats are concentrated. It is also strong across the East Midlands and Lincolnshire.

The deprivation coefficient changes sign across the country. In inner London (Camden, Islington, Westminster) less deprived MSOAs consume more (about +0.22). In the West Midlands (Walsall, Wolverhampton) they consume less (about −0.12). The few significant negative renewable-heating coefficients cluster in north-west London (Harrow, Brent, Barnet), where heat-pump shares are tiny and local collinearity is highest.
- **Off-gas heating, social renting and density** operate at regional to national scales.
- **Age 65+** is effectively global.
- **Household size and detached share** are the most local terms (bandwidths 44 and 53). Their local coefficients are significant in only 12% and 27% of MSOAs after the multiple-testing correction.
- **EPC efficiency** is significantly negative in 40% of MSOAs, never significantly positive.
- **Renewable heating** is mostly non-significant, and negative where it is significant.

### Discussion and limitations

H4 is **partly supported and partly rejected**:

- **Supported:** the heating-fuel effects are local (electric heating) or regional (off-gas), and EPC efficiency is local and remains so after spatial dependence is modelled.
- **Rejected:** household size and detached share do *not* operate at near-global bandwidths. They are the most local terms, but are seldom significant after correction. The likeliest reading is that, once rooms, heating fuel and tenure are held constant, their residual associations are weak and locally idiosyncratic. Collinearity also plays a part: detached share, rooms and household size overlap at MSOA level.

The near-global bandwidth for age 65+ and the regional bandwidths for off-gas heating and social renting identify the processes that a single national coefficient describes adequately.

Several limitations apply:

1. **Ecological covariates.** All covariates are area shares, so coefficients describe places [@fotheringham1991modifiable].
2. **EPC measurement.** EPC is native to MSOAs and selective of the stock that has changed hands since 2008 [@hardy2019analysis].
3. **Outcome.** Electricity per meter mixes standard and Economy 7 meters and is not weather-corrected.
4. **MGWR assumptions.** MGWR assumes the same kernel shape for all processes. Local multicollinearity is summarised in `tables/xs_mgwr_msoa_y2224_diagnostics.csv` (local condition numbers: median 7.8, 95th percentile 24, 2% of MSOAs above the usual warning level of 30, maximum 81). The condition number exceeds 30 mainly in inner London (Tower Hamlets, Westminster, Southwark, Hackney), where flats, density and social renting move together. Local coefficients there, especially for renewable heating and deprivation, are the least reliable.. The very small bandwidths for household size and detached share should be read alongside it.
5. **Effective parameters.** With ENP of about 1,300, the local models use a large share of the 7,264 degrees of freedom. Their maps are best read as descriptions of where associations differ, not as precise local estimates.
6. **Estimator choice at LSOA level.** The S2SLS SDM estimates ρ = 0.98 at LSOA level, which we treat as weak-instrument failure, so the LSOA SDM is not interpreted.

## 3. Models of change

### Methods

**Outcomes.** Two changes in log kWh per meter:

- **Long change:** meter-weighted mean of 2010–12 to that of 2022–24. Three-year averages damp annual weather and measurement noise.
- **Shock change:** 2019 (the last pre-COVID meter year) to 2022 (the first price-crisis meter year, February 2022 to January 2023).

**Regressors.**

- The 2021 baselines from Section 2. These are cross-sections, so they say which kinds of place changed more, not which within-place change caused what.
- The log change in meter counts over the same window. This is the one genuinely time-varying regressor and largely captures new-build growth.
- The initial log level. Including it makes each model a conditional β-convergence regression [@barro1992convergence].

**Exclusions.** At LSOA level we exclude split and irregular (S/X) LSOAs from the long change, because their 2010–12 values are apportioned from the parent, and any LSOA with fewer than 50 meters in a window year. This leaves 33,760 LSOAs for the long change and 35,668 for the shock.

**Estimation.** The same estimator suite is used as in Section 2: GMM at LSOA level; ML plus GWR at MSOA level, with the 13-covariate MSOA set plus meter growth and the initial level. MGWR (the vectorised implementation described in Section 2) is fitted for the long change; the shock change gets GWR.

**Solar PV.** No LSOA or MSOA data on solar PV exist in the open statistics. The final DESNZ sub-regional Feed-in Tariff table (cumulative to March 2019) gives domestic PV installations by local authority. We mapped these to the May 2026 LADs: 2019 codes directly where unchanged, abolished districts through the 2011 OA lookups, and Somerset West and Taunton to Somerset. We then expressed them per 100 households (2021) and added them to a LAD-level change model (n = 318; OLS and ML SEM). The PV file is stored in `external/`.

### Key numbers

**β-convergence** (`tables/ch_beta_convergence.csv`):

| | LSOA | MSOA | LAD |
|---|---|---|---|
| Long change, unconditional β | +0.008 (*t* = 3.5) | +0.017 (*t* = 4.1) | +0.027 (*t* = 1.9) |
| Long change, conditional β | −0.205 | −0.165 | — |
| Shock change, unconditional β | −0.028 (*t* = −15.7) | −0.013 | −0.020 |
| Shock change, conditional β | −0.126 | −0.066 | −0.036 |

- Over 2010–24, consumption per meter fell by 0.215 log points on average (about −19%) at every level.
- In the long change, a positive unconditional β means that higher-consuming places fell slightly *less* in proportion. That is weak β-divergence, matching the Phillips–Sul result.
- Over 2019→22, the mean fall was 0.113 log points (about −11%).
- The LAD "conditional" row in the table controls only for electric heating; the LAD covariate model is reported below.

**LSOA long change** (SEM, GMM, per SD; *n* = 33,760; OLS *R*² = 0.33; λ = 0.50):

| Smaller falls (positive) | Coefficient | Larger falls (negative) | Coefficient |
|---|---|---|---|
| Mean rooms | +0.037 | Initial level | −0.051 |
| Electric heating | +0.032 | Log density | −0.020 |
| Household size | +0.020 | Meter growth | −0.015 |
| Off-gas heating | +0.015 | Wales | −0.007 |
| Flats | +0.014 | District heating | −0.003 |
| Social rent | +0.012 | | |
| Private rent | +0.010 | | |
| Renewable heating | +0.010 | | |
| Detached | +0.006 | | |
| Less deprived | +0.003 | | |

**LSOA shock change** (SEM, per SD; *n* = 35,668; OLS *R*² = 0.26; λ = 0.31):

| Smaller falls (positive) | Coefficient | Larger falls (negative) | Coefficient |
|---|---|---|---|
| Mean rooms | +0.023 | Initial level | −0.031 |
| Electric heating | +0.018 | Log density | −0.010 |
| Flats | +0.013 | Meter growth | −0.010 |
| Private rent | +0.008 | Wales | −0.006 |
| Household size | +0.007 | Aged 65+ | −0.004 |
| Renewable heating | +0.005 | | |
| Off-gas heating | +0.005 | | |
| Less deprived | +0.003 | | |

**LAD change models with solar PV** (`tables/ch_lad_solar_coefs.csv`):

- PV installations per 100 households are not associated with either change once heating, housing and deprivation are controlled. Long change (SEM): −0.0005 per SD (SE 0.0017). Shock change: −0.0006 (SE 0.0009).
- PV correlates 0.68 with the detached share and 0.15 with the deprivation percentile.
- Renewable heating is the clearest positive LAD-level correlate: +0.013 per SD (SE 0.003) for the long change and +0.007 (SE 0.002) for the shock.

**MSOA long change** (standardised; 15 covariates; `tables/ch_msoa_long_model_comparison.csv`):

- **Global fit.** OLS *R*² = 0.41, with residual Moran's *z* = 40. Robust LM-error = 606 against robust LM-lag = 11.
- **SEM** (λ = 0.61), per SD of *y* per SD of *x*:

| Covariate | SEM coefficient |
|---|---|
| Initial level | −0.74 |
| Mean rooms | +0.59 |
| Electric heating | +0.50 |
| Social rent | +0.28 |
| Household size | +0.32 |
| Private rent | +0.26 |
| Flats | +0.23 |
| Renewable heating | +0.19 |
| Off-gas heating | +0.18 |
| Detached | +0.15 |
| Less deprived | +0.08 |
| EPC efficiency | −0.08 |
| Meter growth | −0.33 |
| Log density | −0.34 |

- **AICc:**

| Model | AICc | Residual *z* | Pseudo-*R*² |
|---|---|---|---|
| OLS | 16,811 | | |
| SDM | 14,806 | | |
| GWR (bandwidth 250) | 13,512 | | 0.74 |
| **MGWR** | **13,398** | −4.4 | 0.73 |

- **MGWR** (ENP 998; backfitting converged at iteration 130; local condition numbers all below 15). Bandwidths and share of MSOAs significant after correction:

| Scale | Covariate | Bandwidth | Coefficient | Significant |
|---|---|---|---|---|
| Very local | Initial level | 93 | Mean −0.91 | 99.8%, all negative |
| Very local | Meter growth | 120 | Mean −0.23 | 38%, all negative |
| Local | Aged 65+ | 165 | | 40% |
| Local | Electric heating | 250 | Mean +0.57 | 97% |
| Local | EPC efficiency | 372 | | 15%, all negative |
| Local | Private rent | 496 | | 99.6% |
| Regional | Household size | 1,123 | | 100% |
| Regional | Detached | 2,173 | | 59% |
| Regional | Flats | 3,728 | | 60% |
| Global (≥ 7,089) | Renewable heating | | +0.21 | 100% |
| Global (≥ 7,089) | Off-gas heating | | +0.31 | 100% |
| Global (≥ 7,089) | Mean rooms | | +0.60 | 100% |
| Global (≥ 7,089) | Social rent | | +0.37 | 100% |
| Global (≥ 7,089) | Less deprived | | +0.12 | 100% |
| Global (≥ 7,089) | Log density | | −0.33 | 100% |

**MSOA shock change.** OLS *R*² = 0.41. The best global model is SDM (AICc 15,712). GWR (bandwidth 353) reaches AICc 15,081 and *R*² = 0.64. No MGWR was fitted for the shock change.

### Figures and tables

- Figures: `figures/ch_lsoa_{long,shock}_forest`, `figures/ch_mgwr_msoa_{long,shock}_bandwidths`, `figures/ch_mgwr_msoa_{long,shock}_map_<term>`.
- Tables:
  - `tables/ch_beta_convergence.csv`
  - `ch_{lsoa,msoa}_{long,shock}_{coefs,fit,diag,impacts}.csv`
  - `ch_msoa_{long,shock}_model_comparison.csv`
  - `ch_{gwr,mgwr}_msoa_{long,shock}_local_summary.csv`, `ch_mgwr_msoa_*_local.parquet`
  - `ch_lad_solar_fit.csv`, `ch_lad_solar_coefs.csv`, `ch_lad_solar_fit_stats.csv`

### Results

**Long change, 2010–12 to 2022–24.** Consumption per meter fell by about a fifth almost everywhere, but not uniformly.

*Unconditional pattern.* Higher-consuming places fell slightly *less* in proportion (weak β-divergence at every level). This matches the Phillips–Sul rejection of convergence and the rise in σ.

*Conditional pattern.* Holding baseline composition constant, the initial level is strongly negative. That is the usual signature of conditional convergence, but here it is also inflated by regression to the mean in a noisy 2010–12 baseline. It should not be read as a behavioural catch-up. The places that fell *least* are:

- **electrically heated** areas;
- areas of **large dwellings and larger households**;
- **renting** areas, both social and private;
- areas with more **flats**, and, consistent across every model, those with more **renewable (heat-pump) heating**.

**Meter growth** (new-build) and **density** are associated with larger falls. New homes are more efficient, and dense urban stock uses less per meter.

**Local structure (MGWR).** The MGWR splits these processes by scale:

- Heat-pump heating, dwelling size, social renting, deprivation and density act **uniformly across the country**.
- Electric heating, private renting and EPC efficiency act **locally**.
- The level effect and new-build dilution are **very local**.

In this model MGWR improves on GWR, by ΔAICc 114, and it is better conditioned (local condition numbers all below 15).

**Shock change, 2019 to 2022.** The pattern is the same, with smaller magnitudes:

- Higher-consuming LSOAs fell proportionally more, both unconditionally and conditionally.
- Conditional on level, larger-dwelling, electric-heated, flatted and privately rented areas fell *less*.
- Areas with an older population fell slightly *more*.

The models explain much less of the shock (OLS *R*² = 0.26 at LSOA level) than of the levels. This is consistent with household-level responses to prices and weather [@zapatawebborn2024winter], which area composition captures only partly.

**Solar PV.** At LAD level, domestic PV installed under the Feed-in Tariff shows no association with either change once housing, heating and deprivation are controlled. Renewable heating remains positive at that level too.

### Discussion and limitations

**What the change models show.** The 2021 baselines say which kinds of place changed more. They do not identify mechanisms: a Census 2021 share is not a measure of when heat pumps were installed or rooms were added.

**Electrification.** The positive, nationally uniform association of renewable heating with change is the clearest electrification signal in the data. Where heat pumps had arrived by 2021, consumption per meter fell less than elsewhere. This agrees with the event study (Section 4) and with evidence that heat pumps add substantial electrical load [@love2017addition; @arvanitopoulos2025decarbonising].

**Electric-heated areas in the crisis.** Their smaller conditional fall in 2022 cuts against the expectation that they would cut most. Two explanations are plausible:

- Resistance-heated homes, often flats and rented, have little discretionary load to shed.
- Electricity is not weather-corrected, and the 2022 meter year includes the cold December 2022.

**Solar.** The null LAD result should not be read as "solar does not reduce net consumption" [@mckenna2018solar]. FiT counts stop in 2019, cover only scheme installations, are measured at a coarse scale, and are confounded with detached, owner-occupied housing (*r* = 0.68). Testing the solar effect properly needs LSOA-level MCS or smart-meter export data.

**Measurement limitations.**

- The 2014/15 processing break and the household-weighted apportionment of split LSOAs (excluded here) affect the 2010–12 baseline.
- Regression to the mean biases the conditional β downwards.
- The deprivation result differs from the descriptive one, in which the least-deprived areas cut most in 2021→22 (`findings_dynamics.md`). The difference arises because here the coefficient is conditional on level and composition, and the window starts in the pre-COVID year 2019.

## 4. Spatio-temporal model of the panel

### Why GTWR was replaced

We did not fit a geographically and temporally weighted regression (GTWR) [@huang2010geographically; @fotheringham2015geographical]. The covariates that describe dwellings, heating, tenure and deprivation are Census 2021, EPC FYE2025 and IMD 2019/2025 cross-sections. With those covariates, the temporal variation in local coefficients is not identified: a GTWR would only smooth the outcome through its temporal kernel.

The genuinely time-varying covariates at LSOA level are meter counts (2010–24) and the ONS mid-year population and its age structure (2011–24). After removing LSOA and year effects, they keep only 3–7% of their variance (`tables/st_within_variation.csv`), and log kWh keeps 2.5%. A local panel estimator on 35k × 15 observations with so little within-unit variation would be computationally heavy and statistically noisy. We fitted two identified alternatives instead.

### Methods

**(A) Two-way fixed effects and event study.** Log kWh per meter is regressed on LSOA and year fixed effects plus the time-varying covariates: log meters, log residents per meter, and the shares aged 0–15 and 65+. The sample is the balanced panel for 2011–24.

The event-study extension interacts year dummies with standardised 2021 baselines: electric, renewable and off-gas heating, deprivation, detached share, household size, child share and private renting. 2019 is the reference year. The interactions are identified because they are time-invariant traits multiplied by year. They trace how the consumption gap associated with each baseline evolved through the 2014/15 break, COVID-19 and the price crisis.

Standard errors are clustered by LAD (CR1). The sample is 33,750 LSOAs: S/X boundary changes and any LSOA-year with fewer than 50 meters are excluded.

**(B) Basis-function trend surface.** In the spirit of Fixed Rank Kriging [@cressie2008fixed; @zammitmangion2021frk], relative log consumption (log kWh minus the E&W mean for the year) is projected onto a tensor product of:

- 458 bisquare spatial basis functions on three resolutions: 100, 50 and 25 km grids, with radius 1.5 × spacing;
- 10 piecewise-linear temporal basis functions, with knots at 2010, 2012, 2014, 2015, 2017, 2019, 2020, 2021, 2022 and 2024.

The 4,580 coefficients are fitted by ridge regression. The penalty (λ = 100) was chosen by five-fold cross-validation with LAD blocks (`tables/st_trend_surface_cv.csv`). The fitted surface is the large-scale (≳ 25 km) spatio-temporal trend; the residual is fine-scale variation, which we split into a persistent LSOA component and a transitory one.

This is a Python re-implementation of the FRK idea. The R package was not used, and the model has no explicit fine-scale variance parameter. We therefore treat it as a descriptive decomposition, not a predictive or causal model.

### Key numbers

**Two-way FE, 2011–24** (472,500 LSOA-years, within *R*² = 0.074; `tables/st_fe_timevarying.csv`):

| Term | Coefficient (SE) | Reading |
|---|---|---|
| Log meters | −0.091 (0.007) | Growth in meters lowers the mean, consistent with new-build dilution |
| Log residents per meter | +0.161 (0.009) | Elasticity |
| Share aged 0–15 | −0.157 (0.022) | Per unit share |
| Share aged 65+ | −0.068 (0.021) | Per unit share |

**Event study** (% difference per SD of the baseline, relative to the 2019 gap; `tables/st_event_study.csv`):

| Baseline | Pattern |
|---|---|
| Renewable heating | Rises almost linearly: −0.85% (2010) → 0 (2019) → +0.48% (2022) → +0.45% (2024) |
| Electric heating | About −1.1% before 2015, ≈ 0 for 2017–19, then +0.65% (2022), −0.18% (2023), +0.57% (2024) |
| Less deprived | +0.83% (2020) and +0.60% (2021), back to ≈ 0 by 2022–24 |
| Detached share | +0.31% (2020) and +0.52% (2021) |

**Trend surface** (`tables/st_trend_surface_variance.csv`):

- The large-scale surface explains 17–20% of the cross-sectional variance of relative consumption in every year. Blocked-CV MSE is 0.0354 against a variance of about 0.04.
- The persistent LSOA-level fine-scale component accounts for 67–83%.
- The SD of the surface rose from 0.082 (2010) to 0.090 (2017) and was 0.088 in 2024. The total SD rose from 0.195 to 0.213.

### Figures and tables

- Figures: `figures/st_event_study`, `figures/st_trend_surface_years` (2010, 2019, 2022, 2024), `figures/st_trend_surface_change_{shock,long}`.
- Tables: `tables/st_fe_timevarying.csv`, `st_within_variation.csv`, `st_event_study.csv`, `st_trend_surface_{cv,variance}.csv`, `st_trend_surface_msoa.parquet`; `st_trend_surface_meta.json`.

### Results

The large-scale geography of relative consumption is stable. It is high in the South East outside London, East Anglia, mid and west Wales, Devon and Cornwall. It is low in Tyneside and Teesside, the Pennine conurbations, the South Wales valleys and inner London. The surface for 2024 is almost the 2010 surface with slightly larger amplitude.

Most spatial variation is fine-grained. Four-fifths of the variance between LSOAs lies below the 25 km resolution of the surface, and most of that is persistent. This matches the descriptive analysis, which puts 73–76% of meter-weighted Theil *T* within LADs (`findings_dynamics.md`), and it explains why the convergence clubs sort neighbourhoods rather than regions.

The time-varying panel identifies two compositional effects:

- **New meters lower the mean.** As meter counts grow, typically through new-build flats and estates, consumption per meter falls.
- **More residents per meter raises consumption.** The elasticity is 0.16.

Holding residents per meter constant, a larger share of children or older people slightly lowers consumption per meter.

The event study shows three distinct temporal signatures:

1. **Heat pumps.** Areas with more renewable (largely heat-pump) heating in 2021 have drifted steadily upwards relative to other areas since 2010. The drift is about 1.3% per SD over the period and continued through the price crisis. This is the expected footprint of gradual electrification of heat [@love2017addition; @arvanitopoulos2025decarbonising].
2. **COVID-19.** Less deprived areas and areas of detached housing rose most in 2020–21 and fell back by 2022. This is consistent with a home-working effect concentrated in affluent households [@zapatawebborn2023impact].
3. **Electric heating.** Electric-heated areas moved up relative to others at the 2014/15 break, which may be partly an artefact of DESNZ's processing change. They then oscillate from year to year: up in 2022 and 2024, down in 2023. That pattern is consistent with sensitivity to winter temperatures. Electricity is not weather-corrected, and the 2022 meter year includes the cold December of 2022.

Electric-heated areas therefore did *not* cut consumption more than other areas in 2022. Relative to 2019 they rose slightly.

### Discussion and limitations

Neither model is a GTWR. The event study answers the question a GTWR with time-varying coefficients would have asked (did the association between a 2021 baseline and consumption change over time?) with a transparent, identified design.

The trend surface uses a ridge penalty, not FRK's estimated covariance, so the split between large-scale and fine-scale variance depends on the resolution of the coarsest-to-finest basis (25 km). The persistent fine-scale share is an upper bound on "neighbourhood-specific" variation, because it also absorbs any small-scale processes the basis cannot represent.

Coefficients on residents per meter are partly mechanical, because population estimates and meter counts share a denominator. The electric-heating event-study coefficients cannot be separated from weather without regional degree days [@summerfield2010two], which the panel does not yet include.
