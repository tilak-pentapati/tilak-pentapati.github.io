# Findings: spatial structure, inequality and distribution dynamics

Domestic electricity consumption per meter, LSOAs and LADs in England and Wales, 2010–2024.

**Reproduce.** Run `/Users/jay/anaconda3/envs/envs563/bin/python -I -W ignore projects/python/energy_inequality/analysis/run_dynamics.py`, adding `--fresh` to rebuild the cache. A full run from an empty cache takes about 18 minutes. Outputs go to `data/data_python/energy_inequality/results/dynamics/`: `tables/` (Parquet and CSV), `figures/` (PNG at 300 dpi and SVG), `*_meta.json` and `MANIFEST.json` (package versions and the file list). Spatial weights, LISA results for every year and local tau are cached in `cache/`. All random seeds are fixed (`a_common.SEED`), so a rerun reproduces every table exactly.

**Analysis set.** The panel is balanced: 35,636 of the 35,672 2021 LSOAs have a published or converted value in all 15 years. The other 36 are mostly new-build LSOAs with no meters in 2015–17. LADs follow the May 2026 vintage (318). Relative consumption is the LSOA (or LAD) mean divided by the meter-weighted England and Wales mean of the retained units in that year, which removes common shocks. The core period for dynamics is 2015–2024, which uses one geography and one DESNZ processing system [@desnz2025lsoa]. 2010–2014 are reported separately, and the 2014→15 transition is tested as a break. The 1,906 split LSOAs (ONS change codes S and X [@ons2023lsoalookup]) inherit their parent's mean in 2010–14, so every result touching those years is also run without them.

**Weights.** The main weights are queen contiguity on the 2021 LSOA polygons, row-standardised. They have 5.84 neighbours on average and a minimum of 1. The single island (Isles of Scilly) is joined to its nearest LSOA by centroid distance. Four other components remain: the Isle of Wight, Anglesey and two small tidal islands. They are kept because each still has internal neighbours. The robustness weights are 8 nearest neighbours on population-weighted centroids. LAD weights are queen contiguity, with the Isle of Wight, Isles of Scilly and Anglesey attached to their nearest LAD.

---

## 1. Exploratory spatial data analysis (RQ1)

### Methods

For each year we summarise the LSOA distribution of mean kWh per meter: moments, quantiles, P90/P10 and Bowley skewness. Kernel densities of relative consumption are estimated on the log scale with a Gaussian kernel and Silverman's bandwidth [@silverman1986density], then mapped back to the relative scale. Logging is used because the distribution has a long upper tail, and a fixed rule matters because bandwidth governs apparent multi-modality [@quah1996twin]. Maps use seven Fisher–Jenks classes [@fisher1958grouping; @jenks1967data]. The breaks are estimated once on a 5% sample of the pooled 2015–2024 relative values (sampled Fisher–Jenks [@rey2017evaluation]) and applied to every year, so the maps are comparable over time. Global spatial autocorrelation is measured by Moran's *I* [@moran1950notes] with 999 conditional permutations. It is computed on kWh, on log kWh, with KNN-8 weights, and without split LSOAs. Local clusters are LISAs [@anselin1995local] with 9,999 permutations per year. Significance is controlled with the Benjamini–Hochberg false discovery rate at 5%, because the 35,636 local tests are many and dependent [@caldasdecastro2006controlling]. All statistics use PySAL [@rey2022pysal].

### Key numbers

| | 2010 | 2015 | 2019 | 2020 | 2022 | 2024 |
|---|---|---|---|---|---|---|
| Mean kWh/meter (meter-weighted) | 4,148 | 3,945 | 3,641 | 3,818 | 3,256 | 3,343 |
| Median LSOA | 3,970 | 3,761 | 3,444 | 3,611 | 3,056 | 3,135 |
| CV | 0.216 | 0.222 | 0.236 | 0.241 | 0.241 | 0.246 |
| Skewness (level / log) | 1.68 / 0.84 | 1.97 / 0.91 | 1.97 / 1.04 | 1.92 / 1.00 | 2.12 / 1.17 | 2.08 / 1.15 |
| Moran's *I*, queen (KNN-8) | 0.586 (0.415) | 0.573 (0.398) | 0.582 (0.401) | 0.584 (0.403) | 0.574 (0.395) | 0.575 (0.393) |
| LISA HH / LL after FDR | 3,433 / 2,972 | 3,211 / 2,607 | 2,994 / 2,653 | 3,235 / 2,922 | 2,850 / 2,713 | 3,037 / 2,666 |

- Every Moran's *I* has pseudo-*p* = 0.001, the minimum possible with 999 permutations, for all weights, transformations and years. Log kWh gives almost identical values (0.573–0.592).
- Without split LSOAs, *I* is 0.590–0.604 in every year and shows no step at 2015 (0.598 in 2014, 0.600 in 2015). With them, *I* falls from 0.587 to 0.573. The inherited parent means therefore raise 2010–14 autocorrelation by about 0.01 relative to the clean series.
- The KDE is unimodal in every year (one mode above 2% of the peak), with the mode at 0.89–0.93 of the mean. The mode drifts down (0.92 in 2015, 0.89 in 2024) while the upper tail thickens.
- Pooled Fisher–Jenks upper bounds (× E&W mean): 0.87, 1.01, 1.19, 1.44, 1.76, 2.32 and 4.93.
- About 17% of LSOAs are significant LISAs after FDR (5,800–6,700 per year). The FDR threshold is pseudo-*p* ≤ 0.008–0.009. Without the correction, at *p* < 0.05, 11,400–12,000 are significant. Outliers are rare: HL 85–109 and LH 144–203 per year.
- Queen and KNN-8 LISA labels agree for 85% of LSOAs (Cohen's κ = 0.56 in 2015 and 2024). KNN-8 finds more LL (5,126 vs 2,666 in 2024) and fewer HH (1,808 vs 3,037). With KNN-8, LL clusters in dense urban areas extend further, while rural HH clusters, whose neighbours are spread out, become less distinct.

### Figures and tables

- Figures: `esda_kde_relative`, `esda_map_fj_2015_2024`, `esda_map_fj_{2010,2015,2020,2022,2024}`, `esda_moran_series`, `esda_lisa_map_{2010,2015,2022,2024}`
- Tables: `esda_summary_stats`, `esda_kde_bandwidth`, `esda_kde_curves`, `esda_fj_breaks`, `esda_moran_global`, `esda_lisa_counts`, `esda_lisa_knn_agreement`
- Cache: `cache/lisa_queen.parquet` (LSOA × year: *I*ᵢ, quadrant, pseudo-*p*, raw and FDR cluster labels)

### Results and discussion

Per-meter electricity use is strongly and persistently clustered. Global Moran's *I* stays between 0.57 and 0.59 under contiguity in every year from 2010 to 2024. It barely changes through the pandemic and the price crisis, even though mean consumption fell by a fifth over the period. The cluster map is very stable. High–high clusters cover much of rural England and Wales: Cornwall and Devon, the Welsh uplands, Cumbria and North Yorkshire, Lincolnshire, East Anglia and the outer South East. Low–low clusters are concentrated in inner London, the Tyneside and Teesside conurbations, the South Wales valleys, and parts of West Yorkshire and the East Midlands. The pattern matches the known geography of off-gas-grid, electrically heated and larger detached dwellings in the countryside, and of flats and smaller households in cities [@roberts2015fuel]. These are place-level associations, not household behaviour [@openshaw1984modifiable].

The level of autocorrelation depends on the weights. *I* is about 0.40 with eight nearest neighbours against 0.58 with contiguity, and the agreement between the two LISA maps is moderate rather than high. Contiguity weights give rural LSOAs, which are physically large, few and similar neighbours, whereas KNN-8 reaches further and dilutes them. Both specifications reject spatial randomness overwhelmingly, so the qualitative conclusion does not depend on the choice. The size of any local cluster should still be read with the weights in mind [@fotheringham1991modifiable].

The cross-sectional distribution is unimodal and right-skewed, and it is becoming more skewed. Skewness of log consumption rose from 0.91 in 2015 to 1.15 in 2024, and the mode moved further below the mean. The kernel densities give no sign of twin peaks [@quah1996twin]. The change is a slow lengthening of the upper tail relative to a body that is drifting down. Two caveats apply. First, electricity is not weather-corrected, so year-to-year shifts include winter temperature. Second, kWh per meter is not kWh per household: Economy 7, second homes and communal supplies distort small-area means.

---

## 2. Inequality (RQ1, RQ5)

### Methods

Inequality between LSOAs is measured by the Theil index *T* (GE(1)) and the mean log deviation (GE(0)) [@theil1967economics]. Both are weighted by meters, so they describe the distribution of consumption across meters grouped by place. Both are additively decomposable [@shorrocks1980class]. We decompose *T* exactly into three nested parts: between regions (nine English regions and Wales), between LADs within regions, and within LADs. The parts sum to the total to machine precision. Region membership is read from the DESNZ local-authority table. The between-LAD component is tested against 999 random relabellings of LSOAs to LADs of the same sizes. An unweighted (place-based) Theil is cross-checked against `inequality.theil.TheilD`. The Gini coefficient is reported both place-based and meter-weighted. Following @rey2013spatial, we also compute the spatial Gini, which splits the Gini numerator into contiguous-neighbour pairs and non-neighbour pairs. Its permutation inference (999 permutations) uses a vectorised implementation that matches `inequality.gini.Gini_Spatial` exactly. We report the ratio of mean absolute differences between neighbours and between all pairs, and the polarisation index, whose expected value is 1. The spatial-Gini framing builds on @jacobson2005letting and @rey2025spatial.

For each LAD and year we compute the dispersion of its LSOA means: the CV, a meter-weighted CV and the quartile coefficient of dispersion (QCD), which is robust to the tail. Each is related to the LAD level by a Spearman correlation with a 2,000-draw bootstrap CI. Each is also regressed in OLS as log dispersion on log level, log number of LSOAs and region fixed effects, with HC3 errors. The CV is scale-invariant but not decomposable, so this is descriptive [@cowell2011measuring]. Residual spatial dependence is checked with Moran's *I* on LAD contiguity and with a GMM spatial-error model with heteroskedasticity-robust inference [@arraiz2010spatial]. The level is the LAD mean (total kWh / total meters) and, as a check, the LAD median of LSOA means. The median check guards against a mechanical link in which a heavy upper tail raises both the mean and the CV. The Isles of Scilly LAD has a single LSOA and is excluded, leaving 317 LADs. For RQ5 we compute aggregate (meter-weighted) changes in kWh per meter over each shock window. These are broken down by start-year consumption quintile and by deprivation quintile, using the 2019 IMD/WIMD rank percentile within each country.

### Key numbers

| | 2010 | 2015 | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 |
|---|---|---|---|---|---|---|---|---|
| Theil *T* (meter-weighted) | 0.0217 | 0.0222 | 0.0253 | 0.0262 | 0.0270 | 0.0259 | 0.0244 | 0.0269 |
| between regions, % | 9.2 | 9.3 | 9.8 | 10.2 | 9.3 | 9.2 | 8.7 | 8.5 |
| between LADs within region, % | 18.2 | 17.1 | 16.4 | 16.1 | 16.4 | 15.6 | 15.4 | 15.8 |
| within LADs, % | 72.6 | 73.6 | 73.7 | 73.6 | 74.3 | 75.2 | 75.8 | 75.7 |
| Gini (LSOAs) | 0.111 | 0.111 | 0.117 | 0.120 | 0.121 | 0.118 | 0.115 | 0.121 |
| Neighbour / all-pairs mean abs. difference | 0.68 | 0.70 | 0.70 | 0.70 | 0.70 | 0.70 | 0.71 | 0.70 |
| Spatial polarisation index | 1.46 | 1.43 | 1.43 | 1.43 | 1.43 | 1.42 | 1.42 | 1.42 |

- The MLD gives the same split (within-LAD share 71–74%). The unweighted Theil and its between-LAD share (25–27%) match PySAL. The observed between-LAD component (about 0.006) is roughly 25 times its permutation mean (0.0002–0.0003), *p* = 0.001.
- Spatial Gini: the non-neighbour component exceeds its permutation expectation in every year (*z* = 101–116, *p* = 0.001). Meter-weighted and place-based Gini agree to the third decimal.
- MSOA scale (`sens_theil_msoa`): 68–69% of *T* lies between MSOAs and 31–32% within them. About 42–43% lies between MSOAs within the same LAD.
- Dispersion vs level, all 15 years. Spearman ρ(CV, LAD mean) ranges from 0.68 [95% CI 0.61–0.74] in 2010 to 0.75 [0.70–0.80] in 2022–24. The OLS elasticity of CV with respect to the mean is 2.3–2.9 (HC3 SE about 0.14, *p* < 10⁻⁶). The spatial-error estimate is 2.3–2.95, with λ between 0.26 and 0.40. Other checks give the same picture: QCD, ρ 0.64–0.74 and elasticity 2.9–3.5; excluding outlier LSOAs, ρ 0.63–0.72 and elasticity 2.1–2.7; median as level, CV ρ 0.54–0.58 with elasticity 2.0–2.5, and QCD ρ 0.54–0.60. OLS residuals remain spatially autocorrelated (Moran's *I* 0.11–0.20, *p* ≤ 0.004).
- Shock response, aggregate change in kWh per meter, 2021→22 (E&W −8.8%):

  | Group | lowest / most deprived | highest / least deprived |
  |---|---|---|
  | Consumption quintile | −6.8% | −9.6% |
  | Deprivation quintile | −6.8% | −10.0% |

  The gradient is monotonic. In 2019→20 (E&W +4.9%), the least deprived quintile rose most (+5.9%) and the most deprived least (+3.5%). Over 2019→24, the gap between quintiles is small: −6.6% to −9.1% by consumption and −7.5% to −8.5% by deprivation.

### Figures and tables

- Figures: `ineq_theil_gini`, `ineq_cv_vs_mean_2024`, `ineq_cv_mean_rho_series`, `ineq_shock_2022_by_group`
- Tables: `ineq_theil`, `ineq_gini_spatial`, `ineq_lad_dispersion`, `ineq_dispersion_vs_mean`, `ineq_shock_response`, `sens_theil_msoa`

### Results and discussion

Inequality in per-meter consumption between LSOAs is modest in absolute terms (Gini about 0.11–0.12) but has risen. The Theil index grew by about a fifth between 2015 and 2021, and the Gini from 0.111 to 0.121. Most of this inequality is local. Roughly three quarters lies within local authorities, a share that rose from 73% to 76% over the period, and almost a third lies within MSOAs. Regions account for less than a tenth. The between-region and between-LAD shares have fallen slowly since about 2020, so the growth in inequality has been within places rather than between them. The between-group shares depend on how many groups there are and how large they are [@openshaw1984modifiable]. Even so, the dominance of the within-LAD component holds whether LSOAs are weighted by meters or counted equally, and whether Theil *T* or the MLD is used.

The spatial Gini confirms that the inequality is spatially structured. Contiguous LSOAs differ on average by 30% less than random pairs, and the share of total inequality carried by non-neighbour pairs is far above what spatial randomness would produce (*z* > 100). This is the inequality-space counterpart of the high Moran's *I* in Section 1.

Within-LAD dispersion rises steeply with LAD consumption. In every year, the authorities with the highest mean use per meter are also the most internally unequal. The elasticity of the CV with respect to the mean (about 2.5) holds with tail-robust dispersion, with outlier LSOAs removed, with region fixed effects and LAD size controlled, under a spatial-error specification, and with the LAD median as the level. The last check rules out a pure tail artefact. The relationship strengthened slightly after 2020. We read this as polarisation within high-consumption authorities. These are largely rural LADs in which electrically heated, off-gas-grid LSOAs sit beside gas-heated market towns, so the high mean and the high dispersion share a compositional cause. The association is descriptive and should not be read as a causal effect of consumption level on dispersion.

On RQ5, the price crisis compressed relative differences for a short time. In 2021→22, high-consumption and less deprived LSOAs cut consumption by about 10%, against about 7% in the lowest-consumption and most deprived quintiles. The Gini fell from 0.121 to 0.118 in 2022 and to 0.115 in 2023, but was back at 0.121 in 2024, and the cumulative 2019→24 changes differ little across quintiles. This pattern fits evidence that the poorest households had little room left to cut and that much of the 2022 reduction was rationing rather than efficiency [@zapatawebborn2024winter; @fawcett2024using]. Lower consumption in deprived areas should not be read as a better outcome. The pandemic year shows the opposite social gradient. Consumption rose most in the least deprived quintile, which is consistent with home working concentrated among higher-income households and with the national smart-meter estimates [@zapatawebborn2023impact]. These are area-level contrasts on non-weather-corrected data, and the meter years run February to January, so "2022" covers February 2022 to January 2023.

---

## 3. Distribution dynamics (RQ2, RQ5)

### Methods

Each year, LSOAs and LADs are assigned to five classes of relative consumption. The class bounds are pooled quintiles over 2010–2024 and fixed across all windows, so matrices and ergodic distributions are comparable. The LSOA bounds are 0.842, 0.914, 0.987 and 1.116 × E&W mean; the LAD bounds are 0.912, 0.973, 1.042 and 1.128. We estimate classic first-order Markov chains [@quah1996twin; @quah1997empirics] for 2015–24 (core), 2010–24, each sub-period (2010–14, 2015–19, 2020–24) and each annual transition. From each chain we derive ergodic distributions, mean first passage times (MFPT), sojourn times 1/(1 − *p*ᵢᵢ) and the Shorrocks mobility index (*k* − trace *P*)/(*k* − 1) [@shorrocks1978measurement].

Time homogeneity is tested with the Bickenbach–Bode *Q* and likelihood-ratio statistics [@bickenbach2003evaluating] across the three sub-periods, across the nine annual core transitions, and across 2013→14, 2014→15 and 2015→16, with 2013→14 vs 2015→16 as a no-break control. The first-order Markov property is tested by comparing t→t+1 transition matrices conditioned on the class at t−1. It is also checked by comparing observed 2- and 5-year transition diagonals with those implied by *P*² and *P*⁵.

Spatial Markov chains condition transitions on the quintile of the spatial lag [@rey2001spatial]. Spatial homogeneity is tested with *Q* and LR [@rey2016properties; @kang2018conditional] and with a permutation version of the χ² statistic. The permutation version (99 random spatial permutations) does not assume that transitions are independent.

LISA Markov chains track transitions between Moran-scatterplot quadrants [@rey2001spatial]. Following the caveat in @rey2016properties, they are run twice: on all quadrants, and conditioned on FDR significance, with five states (HH, LH, LL, HL and not significant). Independence of own and neighbourhood moves is tested with a χ² test. Directional LISA rose diagrams [@rey2011visualizing] use eight sectors over the joint movement of an LSOA and its spatial lag in relative-consumption space. Sector counts are tested against 999 random spatial permutations (two-sided).

### Key numbers

**Classic Markov, 2015–24 (320,724 LSOA transitions; 2,862 LAD transitions)**

| from \ to | Q1 | Q2 | Q3 | Q4 | Q5 | sojourn (yrs) | ergodic |
|---|---|---|---|---|---|---|---|
| LSOA Q1 | **0.883** | 0.114 | 0.002 | 0.000 | 0.000 | 8.6 | 0.264 |
| LSOA Q2 | 0.137 | **0.722** | 0.138 | 0.003 | 0.000 | 3.6 | 0.220 |
| LSOA Q3 | 0.003 | 0.171 | **0.700** | 0.125 | 0.000 | 3.3 | 0.178 |
| LSOA Q4 | 0.001 | 0.003 | 0.136 | **0.796** | 0.064 | 4.9 | 0.163 |
| LSOA Q5 | 0.000 | 0.000 | 0.001 | 0.060 | **0.939** | 16.5 | 0.174 |

- LAD diagonals are 0.957, 0.911, 0.879, 0.893 and 0.950.
- Shorrocks index, LSOA: 0.240 (2015–24), 0.220 (2010–14), 0.208 (2015–19) and 0.272 (2020–24).
- Shorrocks index, LAD: 0.103 (2015–24), 0.074 (2010–14), 0.071 (2015–19) and 0.140 (2020–24).
- Annual LSOA Shorrocks: 0.20–0.21 in 2015–19, then 0.240 (2019→20), 0.252 (2020→21), **0.307 (2021→22)**, 0.265 and 0.264. The LAD series peaks at 0.184 in 2021→22.
- MFPT, LSOA Q1→Q5 = 128 years and Q5→Q1 = 69 years. For LADs they are 230 and 181 years.
- Time homogeneity is rejected for the three sub-periods: *Q* = 1,920, df = 40, *p* < 10⁻¹⁰. The largest difference in any transition probability is 0.054, and the mean diagonal is 0.824, 0.834 and 0.782.
- For LADs time homogeneity is also rejected (*Q* = 56.7, df 16, *p* < 0.001).
- Break at 2015: the 2014→15 LSOA transition is more mobile than its neighbours, with mean diagonal 0.810 against 0.836 and 0.837 (*Q* = 590, df 40). Without split LSOAs the gap shrinks to 0.827 against 0.839 and 0.840. LADs show no break (*Q* = 16.0, df 16, *p* = 0.45).
- First-order property rejected (*Q* = 52,681, df 80, largest difference in any transition probability 0.62). The first-order chain also understates persistence. The observed 5-year diagonals for Q1–Q5 are 0.85, 0.57, 0.53, 0.66 and 0.91, against 0.63, 0.36, 0.33, 0.42 and 0.76 implied by *P*⁵.

**Spatial Markov, LSOAs 2015–24.**

- *Q* = 6,150 and LR = 5,715, df 80, *p* < 10⁻¹⁰. Permutation χ² *p* = 0.01, the minimum with 99 permutations: the observed 6,150 against a permutation maximum of 176. Each lag-class χ² is significant.
- Probability of leaving Q1 upwards: 0.084 when neighbours are in Q1, rising to 0.143, 0.166, 0.179 and 0.185 across lag quintiles 2–5.
- Probability of leaving Q5 downwards: 0.155 when neighbours are in Q1 and 0.031 when they are in Q5.
- Ergodic share in Q1: 0.59 in the lowest lag class and 0.04 in the highest.
- LADs: *Q* = 120, df 32, *p* < 0.001, but only the two lowest lag classes differ significantly from the pooled chain.

**LISA Markov, 2015–24.**

- 92.1% of LSOA-years keep their quadrant. 4.6% move through their own value only, 3.0% through the lag only, and 0.26% jointly. Independence of own and lag moves is rejected (χ² = 729,816, df 9).
- With FDR significance, sojourn times are 17.2 years for HH and 8.2 years for LL. The ergodic shares are HH 7.8%, LL 7.8% and not significant 83.7%.

**Rose diagrams (2015→2024).**

- The four co-movement sectors (own and neighbourhood moving the same way) are over-represented. For example, sector 5 (both down) has 8,596 against 7,334 expected, and sector 6 has 4,836 against 3,092.
- Opposed moves are under-represented: sector 4 has 3,819 against 6,450 and sector 8 has 3,262 against 5,914. All sectors give *p* ≤ 0.002.
- The same pattern holds in every window, including 2019→20 and 2021→22.

### Figures and tables

- Figures: `markov_classic_lsoa_core`, `markov_classic_lad_core`, `markov_spatial_lsoa_extremes`, `markov_shorrocks_annual`, `markov_lisa_fdr`, `rose_directional_lisa`
- Tables: `markov_class_bounds`, `markov_matrices` (long format: every *P* and MFPT, including conditional and LISA chains), `markov_summary`, `markov_annual`, `markov_time_homogeneity_tests`, `markov_spatial_homogeneity_tests`, `markov_spatial_by_lag`, `markov_spatial_updown`, `rose_counts`, `sens_markov`
- JSON: `markov_meta.json` (LISA Markov tests, multi-step diagonals)

### Results and discussion

Places keep their relative positions. Between 70% and 94% of LSOAs stay in the same quintile from one year to the next. The extremes are the most persistent: 0.88 for the lowest quintile and 0.94 for the highest. At LAD level the diagonals are 0.88–0.96. Moves almost always go to an adjacent class: in a typical year only 0.15–0.3% of LSOAs move two or more quintiles, and 0.46% in 2021→22. The mean first passage time from the lowest to the highest quintile is over a century.

Persistence is, if anything, understated by the first-order chain. The Markov property is decisively rejected, and the observed five-year probability of remaining in Q1 (0.85) or Q5 (0.91) is far higher than *P*⁵ implies (0.63 and 0.76). Much of the year-to-year mobility near class boundaries is short-lived back-and-forth: estimated readings, meter churn and weather all add transitory noise to annual means. The one-step chain mistakes this noise for movement. The sojourn times (8.6 years for Q1 and 16.5 years for Q5 under the first-order chain) are therefore lower bounds on how long places actually stay. The ergodic distributions are projections of a process that is neither first-order nor time-homogeneous, so they are illustrative rather than forecasts. With that caveat, the core-period chain puts more mass in the lowest class (0.26) than the pooled quintiles do (0.20). This is the dynamic counterpart of the growing skew in Section 1: the body of the distribution is drifting below the mean while a persistent group of high-consumption places stays above it.

Mobility is not constant over time. Each test rejects time homogeneity, although with 140,000 or more transitions per sub-period even small differences are significant. The effect sizes are what matter: the largest change in any transition probability between sub-periods is 0.05. Mobility was lowest in 2015–19, rose through the pandemic and peaked in 2021→22, when the Shorrocks index reached 0.31 for LSOAs and 0.18 for LADs. It has remained above its pre-2020 level since. The shocks did not affect all places uniformly in relative terms, and they re-sorted places, especially in the middle of the distribution.

The 2014→15 transition is more mobile than its neighbours. About two thirds of that excess disappears when split LSOAs are excluded, and none is visible at LAD level. Most of the break is therefore an artefact of apportioning 2011 LSOAs to 2021 LSOAs. The rest reflects the change in DESNZ processing systems, and this is why 2015–24 is the core period.

Neighbourhood context conditions movement strongly. Spatial homogeneity is rejected by the analytical tests and by a permutation test that does not assume independent transitions. The effect is large and symmetric:

- A low-consumption LSOA surrounded by low-consumption neighbours leaves the bottom quintile upwards with probability 0.08, against 0.19 when its neighbours are in the top quintile.
- A high-consumption LSOA among high-consumption neighbours drops out of the top quintile with probability 0.03, against 0.16 among low-consumption neighbours.

Places are pulled towards their neighbourhood's level. As a result, spatial clusters reproduce themselves, and the conditional ergodic distributions are almost mirror images (59% in Q1 under the lowest lag class, 59% in Q5 under the highest). The LISA Markov chain and the rose diagrams tell the same story from different angles. Quadrant membership is very stable: significant HH clusters last 17 years on average. When LSOAs move, they move with their neighbours far more often than spatial randomness predicts, and moves against the neighbourhood are rare. At LAD scale the spatial conditioning is weaker and confined to the lower lag classes, which is expected with 318 large, internally heterogeneous units. These are descriptive dependencies. They do not identify spillovers, because shared exposure to climate, housing stock and the gas grid would produce the same pattern.

---

## 4. Place mobility and rank concordance (RQ2, RQ5)

### Methods

Following @vankerm2004what, as applied to places by @rey2023place [see also @rey2010interregional; @rey2023inequality], we decompose mobility between two years into growth, dispersion and exchange components. Mobility is the mean absolute change across places, measured in kWh per meter and in log points. Three counterfactual factors are switched on in turn:

- growth: the end-year mean;
- dispersion: the end-year shape, i.e. the quantile function divided by its mean;
- exchange: the end-year ranks.

Each component is the Shapley value over all orderings of the three factors [@shorrocks2013decomposition], and the components sum exactly to the total. Windows are each annual transition plus 2010–15, 2015–19, 2019–24, 2015–24, 2010–24 and 2019–22. The decomposition is run for LSOAs (all, without outlier-flagged LSOAs, and without split LSOAs) and for LADs.

Rank concordance between years is measured by Kendall's τ_b, by spatial τ (the τ computed over contiguous pairs only, with 999 permutations) and by local neighbour τ (LIMA) [@rey2016space] with 999 conditional permutations and FDR. The LSOA LIMA is a vectorised implementation that reproduces `giddy.rank.Tau_Local_Neighbor` on the LAD data. Neighbours start with similar values, and pairs that start close together swap ranks more easily wherever they are. We therefore also compare spatial τ with a value-matched null: 199 samples of random pairs, each the size of the neighbour-pair set and stratified to the neighbour pairs' distribution of |Δx| over 20 quantile bins.

### Key numbers

**LSOAs, log-point index (× 100)**

| window | total | growth | dispersion | exchange |
|---|---|---|---|---|
| typical year 2015–19 | 2.4–3.3 | 0.6–2.0 (25–58%) | 0.04–0.22 | 1.3–1.8 (40–73%) |
| 2019→20 | 5.0 | 3.6 (71%) | 0.13 | 1.3 (27%) |
| 2020→21 | 6.9 | 5.4 (78%) | 0.12 | 1.4 (20%) |
| 2021→22 | 9.3 | 7.7 (82%) | 0.13 | 1.5 (16%) |
| 2023→24 | 2.7 | 0.2 (8%) | 0.35 | 2.2 (79%) |
| 2015→24 | 17.2 | 14.1 (82%) | 0.8 | 2.4 (14%) |
| 2010→24 | 22.4 | 18.6 (83%) | 0.9 | 2.9 (13%) |

- In kWh, 2021→22 is total 322, growth 262, dispersion 6 and exchange 54. Over 2015→24 it is total 615, growth 495, dispersion 26 and exchange 94.
- For LADs, growth is 92–95% of mobility in 2019–22 and exchange is 4–6%.
- Excluding outlier-flagged LSOAs changes no component by more than 0.0012 log points. Excluding split LSOAs leaves growth and dispersion essentially unchanged but lowers exchange: by 0.44 points in 2014→15 (2.0 against 2.5), the geography break, and by about 0.1 points a year afterwards. Split LSOAs, many of them growing new-build areas, re-rank more than average.
- Gas on the same LSOAs: 2021→22 total 14.2, of which growth is 12.9, and exchange is 0.8–2.2 per year.
- Kendall τ_b between consecutive years: 0.86–0.90 for LSOAs (lowest in 2021→22) and 0.92–0.97 for LADs. Between 2015 and 2024 it is 0.76 for LSOAs and 0.88 for LADs.
- Spatial τ, LSOAs: below global τ in every window (for example 0.697 against 0.757 for 2015–24; permutation *p* = 0.001). It is at or above the value-matched null, for example 0.697 against 0.668 [0.663, 0.672] for 2015–24. In annual windows it is within or slightly above the null.
- Spatial τ, LADs: not different from global τ (*p* 0.04–0.44) but above the value-matched null in most windows.
- LIMA, LSOAs: mean local τ is 0.69 for 2015–24, and 3.5% of LSOAs have a negative local τ. 2,409 are significant at raw *p* < 0.05, which is close to the roughly 1,800 false positives expected by chance, and none survive FDR. LADs: 8 at raw *p* < 0.05, none after FDR.

### Figures and tables

- Figures: `mobility_vankerm_annual`, `rank_tau_annual`, `rank_lima_map_2015_2024`
- Tables: `mobility_vankerm`, `rank_concordance`, `rank_lima_summary`, `rank_lima_lad_2015_2024`, `sens_gas_vankerm`
- Cache: `cache/lima_lsoa.parquet`

### Results and discussion

The decomposition separates the common shocks from the re-sorting of places. In shock years, mobility is dominated by growth. Between 2019 and 2022, 71–82% of LSOA mobility and over 90% of LAD mobility is the national level moving up (COVID-19) and then sharply down (the price crisis). In calm years, and over the 2013–15 plateau when the national mean barely moved, exchange dominates. Dispersion, meaning a change in the shape of the distribution with ranks held fixed, is consistently the smallest component, at a tenth or less of exchange. Exchange is the spatially interesting part, and it is small and roughly constant. It is about 1.3–2.2 log points a year, and over nine years it accumulates to only 2.4. If annual re-rankings were permanent they would sum to more than 13 points, so most of them reverse. This confirms from a different angle the transitory, boundary-hopping mobility found in the Markov analysis. Estimated meter readings and annual noise inflate one-year exchange, so the multi-year windows are the better guide to genuine re-ranking. Gas shows the same structure: a very large 2022 growth component and modest exchange.

Rank concordance is high and falls slowly with the window length. It falls to 0.76 (LSOAs) and 0.88 (LADs) over nine years, with the clearest single-year disruption in 2021→22. Spatial τ for LSOAs is well below global τ, which appears to say that neighbours swap ranks more than random pairs do. The value-matched null shows this is a level effect. Neighbours start with similar values, and similar places swap ranks easily wherever they are. Against equally similar random pairs, neighbouring LSOAs re-rank no more, and over longer windows somewhat less. This matches the co-movement in the rose diagrams. Local rank concordance gives no evidence of hotspots of re-ranking: no LSOA or LAD is significant after FDR correction. The LIMA map should therefore be read as descriptive. Neighbourhoods do not form local pockets of churn. Re-ranking is diffuse and small-scale, and it is superimposed on stable neighbourhood levels.

---

## 5. Sensitivity

### Methods

1. **Outliers.** Every cross-sectional and dynamic statistic was re-estimated excluding the 216 balanced-panel LSOAs flagged as outliers in any year (`flag_outlier`), leaving 35,420.
2. **Split LSOAs.** Results touching 2010–14 were also re-estimated without the S/X LSOAs (33,751 remain).
3. **Gas.** Gas, which DESNZ weather-corrects, was analysed alongside electricity on the 34,551 LSOAs with complete series for both fuels. This makes the comparison like-for-like and excludes off-gas-grid areas from both.
4. **Scale.** Theil was also decomposed by MSOA.
5. **Independence.** The spatial Markov χ² was given a permutation *p*-value that does not assume independent transitions.

### Key numbers

**Excluding outlier LSOAs**

- Moran's *I* changes by −0.003 to −0.011 (0.561–0.580).
- Gini falls by about 0.004 (0.105–0.117) and the between-LAD Theil share by 0.003–0.006.
- Markov core diagonals change by ≤ 0.003 and the Shorrocks index is 0.242 against 0.240.
- Spatial Markov: *Q* = 6,070. The upward exit from Q1 is 0.087 with low-quintile neighbours and 0.196 with top-quintile neighbours.
- The dispersion–mean elasticity is 2.1–2.8.

**Gas vs electricity on the same LSOAs**

| | Gas | Electricity |
|---|---|---|
| Moran's *I* | 0.50–0.53 | 0.51–0.53 |
| Gini | 0.122–0.132 (flat) | 0.100–0.111 (rising) |
| Between-LAD Theil share | 0.21–0.24 (rising after 2021) | 0.23–0.25 (falling) |
| Shorrocks index, 2015–24 | 0.134 | 0.254 |
| Upward exit from Q1, low vs high neighbours | 0.052 vs 0.118 | 0.086 vs 0.190 |

On these LSOAs, gas is less mobile than electricity, and its spatial conditioning is weaker but the same in direction.

**Spatial Markov permutation test.** χ² = 6,150 against a permutation mean of 118 and maximum of 176 (99 permutations, *p* = 0.01).

### Figures and tables

- Figures: `sens_elec_gas_outliers`
- Tables: `sens_cross_section`, `sens_markov`, `sens_gas_vankerm`, `sens_theil_msoa`, plus `excl_flagged` rows in `ineq_lad_dispersion`/`ineq_dispersion_vs_mean` and "LSOA excl. flagged/S/X" rows in `mobility_vankerm`, `markov_*`
- JSON: `sens_meta.json`

### Discussion

None of the substantive conclusions depends on the outlier flags, the split LSOAs, the weights or the inference method. The flagged LSOAs are real, mostly electrically heated places in the far upper tail. Removing them lowers inequality slightly but leaves clustering, persistence, spatial conditioning and the mobility decomposition essentially unchanged.

The gas comparison is informative. Restricting electricity to gas-connected LSOAs lowers its Moran's *I* from about 0.58 to 0.52. Much of electricity's clustering therefore comes from off-gas-grid areas, where electricity substitutes for heating fuel [@roberts2015fuel]. Even on gas-connected LSOAs, electricity is the more mobile fuel, and its inequality is rising while gas inequality is flat. One plausible reason is that gas is weather-corrected and electricity is not. Another is that electricity increasingly reflects uneven new loads (heat pumps, electric vehicles) and self-generation (rooftop PV) [@love2017addition; @mckenna2018solar]; the next section of the study examines these.

---

## Hypotheses: verdicts

| Hypothesis (from `literature.md`) | Verdict | Evidence |
|---|---|---|
| H1a: Moran's *I* positive and significant every year | **Supported** | *I* = 0.57–0.59 (queen), 0.39–0.42 (KNN-8), *p* = 0.001 in all 15 years and all variants |
| H1b: within-LAD component dominates Theil | **Supported** | 73–76% within LADs (MLD 71–74%), and the share is rising; 31–32% lies within MSOAs |
| H1c: within-LAD dispersion systematically related to LAD mean | **Supported** | ρ = 0.68–0.76, elasticity about 2.5; robust to QCD, outliers, median level and spatial error |
| H2a: high persistence on Markov diagonals | **Supported** | LSOA 0.70–0.94 and LAD 0.88–0.96 per year; true multi-year persistence is higher still |
| H2b: spatial Markov rejects homogeneity; low-consumption neighbourhoods reduce upward mobility | **Supported** | *Q* = 6,150 (df 80), permutation *p* = 0.01; P(up from Q1) 0.08 vs 0.19. The effect is symmetric, with high neighbourhoods also holding places at the top |
| H2c: growth dominates in shock years | **Supported** | 71–82% of LSOA mobility (over 90% for LADs) in 2019→20, 2020→21 and 2021→22 |
| H2d: exchange outweighs dispersion | **Supported** | Exchange is about 10× dispersion in most windows. Most annual exchange is transitory |
| H5a: relative inequality fell in 2022 as high-consuming areas cut most | **Partly supported** | The top consumption quintile cut 9.6% against 6.8% in the bottom quintile, and the Gini fell 0.121 → 0.118 → 0.115 (2023). The compression reversed in 2024 (0.121) |
| H5b: low-consumption, more deprived areas show the smallest reductions | **Supported (area level)** | The most deprived quintile cut 6.8% against 10.0% in the least deprived in 2021→22; the gradient is monotonic |

RQ3 (convergence clubs) and RQ4 (heterogeneous drivers) are addressed in a separate part of the study.

## Limitations

- **Meters are not households.** Means per meter mix dwelling size, heating fuel, occupancy, multi-rate tariffs and empty or second homes. All results describe places, not behaviour (ecological inference).
- **Weather is not removed.** Electricity is not weather-corrected, so some year-to-year movement, particularly the national growth component, reflects winter temperatures. The relative measures largely remove the common part of that effect.
- **Two geographies and two processing systems.** 2010–14 are converted from 2011 LSOAs and come from DESNZ's old system. 2015–24 is the core period throughout, and the 2014→15 break is documented above.
- **Formal tests are over-powered.** With 35,000 units, almost any difference is "significant", and transitions are neither independent nor first-order. We therefore report effect sizes, permutation-based checks and multi-step comparisons alongside the asymptotic *p*-values.
- **Results depend on the weights.** Moran's *I* and LISA maps differ between contiguity and KNN weights. The conclusions do not change, but cluster extents do.
- **Inequality between meters, not households.** Theil is meter-weighted and Gini is place-based (they agree closely). Neither measures inequality between households, which is much larger within LSOAs than between them.
