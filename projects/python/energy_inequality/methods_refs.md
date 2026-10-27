## Methods and their canonical references

Each planned method is listed with the citations to use and the objections a reviewer is likely to raise. Keys refer to `references.bib`.

### Data and geography

| Element | Cite | Assumptions and caveats reviewers will raise |
|---|---|---|
| DESNZ LSOA domestic electricity, 2010–latest | @desnz2025lsoa; @desnz2026subnational | Domestic means NHH profile classes 1–2 (standard and Economy 7), with automatic reallocation above 100,000 kWh. **Not weather-corrected.** The meter year runs from late January to late January, so "2020" covers Feb 2020–Jan 2021 and "2022" covers Feb 2022–Jan 2023. About 20% of readings are estimated. Meters ≠ households. **Only 2015 onwards is published on 2021 LSOAs with consistent postcodes; 2010–2014 are on 2011 LSOAs and are not revised.** Test for a 2015 break, or treat 2015+ as the core panel. Mean per meter is skewed, so report the median alongside. Metered consumption is *net import*: PV self-consumption lowers it (@mckenna2018solar), while heat pumps and EVs raise it (@love2017addition). |
| 2011→2021 LSOA harmonisation | @ons2023lsoalookup; @ons2021statgeog; @cockings2011maintaining | The exact-fit lookup codes U/S/M/X. Splits (S) can be apportioned by meter counts; merges (M) can be summed. **X (irregular) cases cannot be compared like-for-like**, so drop them or flag them. State the counts affected. |
| EPC covariates | @mhclg2025epc; @hardy2019analysis; @few2023overprediction; @summerfield2019empirical | The register covers only dwellings sold, let or built since 2008, so LSOA aggregates are biased towards the rental and turnover stock. Records contain errors. EPCs are *modelled* demand and over-predict most for bands D–G. Report EPC coverage per LSOA and use it as a weight or a sensitivity filter. |
| Census 2021 covariates | @ons2021statgeog | These are a single cross-section, so they cannot explain *within-LSOA change over time* (see the GTWR row). The inference is ecological. |
| MAUP and ecological inference | @openshaw1984modifiable; @fotheringham1991modifiable | Repeat the key results at MSOA (and LAD) scale. Phrase findings as place-level associations, not household behaviour. |

### Exploratory spatial data analysis

| Method | Cite | Caveats |
|---|---|---|
| Kernel density of the cross-sectional distribution | @silverman1986density; @quah1996twin | Bandwidth choice drives the apparent multi-modality. Report the bandwidth rule. Use log or relative consumption (LSOA ÷ national mean per year) to remove common shocks. |
| Fisher–Jenks classification | @fisher1958grouping; @jenks1967data; @rey2017evaluation | With about 35k units, use sampled Fisher–Jenks. Fix class breaks across years, or use quantiles, to keep maps comparable over time. |
| Global Moran's *I* | @moran1950notes; @rey2022pysal | State the weights (queen contiguity vs *k*-NN; islands such as Isles of Scilly and Anglesey links). Use permutation inference. Run a sensitivity check on alternative weights. |
| LISA clusters | @anselin1995local; @caldasdecastro2006controlling | Local tests are many and dependent, so apply FDR or a Bonferroni-type adjustment. Pseudo-*p* values depend on the number of permutations. |

### Inequality and dispersion

| Method | Cite | Caveats |
|---|---|---|
| Coefficient of variation vs LAD mean | @cowell2011measuring | The CV is scale-invariant but not decomposable. Its correlation with the mean is descriptive, so use it as motivation and not as a test. |
| Theil decomposition (within/between LAD, region) | @theil1967economics; @shorrocks1980class; @rey2025spatial | Additively decomposable. The between-group share depends on the number and size of groups (MAUP again). Weight by meters. |
| Gini and spatial Gini | @jacobson2005letting; @rey2013spatial | The standard Gini is not additively decomposable by subgroup. Use the spatial Gini (neighbour vs non-neighbour pairs) with permutation inference instead. |

### Distribution dynamics

| Method | Cite | Caveats |
|---|---|---|
| Classic Markov chain | @quah1996twin; @quah1997empirics | Class boundaries are arbitrary, so use quantiles of *relative* consumption. Test time homogeneity (pre-/post-2020) and the first-order Markov property (@bickenbach2003evaluating). The ergodic distribution assumes stationarity, which shocks violate. |
| Spatial Markov chain | @rey2001spatial; @rey2016properties; @kang2018conditional | Use the likelihood-ratio / χ² tests of spatial homogeneity. Report conditional *and* joint tests. Results are sensitive to how the spatial lag is classed. |
| LISA Markov and directional rose diagrams | @rey2011visualizing; @rey2016space | Transitions between quadrants mix significant and non-significant LISAs, so either condition on significance or say so. Rose-diagram inference uses permutation. |
| Place-mobility decomposition (growth / dispersion / exchange) | @vankerm2004what; @shorrocks1978measurement; @rey2023place; @rey2010interregional; @rey2023inequality | Growth captures common shocks (COVID-19, prices). Exchange (re-ranking) is the spatially interesting part. Measurement noise from estimated readings inflates apparent exchange mobility, so compare 1-year and multi-year windows. |

### Club convergence

| Method | Cite | Caveats |
|---|---|---|
| Phillips–Sul log *t* test and clustering | @phillips2007transition; @phillips2009economic; @sichera2019convergenceclubs; @bartkowska2012regional | **T is short (about 10–15 years) relative to N (about 35k).** PS recommend discarding the first third of the sample and filtering out the business cycle (HP filter). The test's size and power with short T, and the computational cost with very large N, need to be justified. Consider running at MSOA level, or within regions. Clubs are a-spatial, so test club membership for spatial autocorrelation afterwards. Shocks in 2020 and 2022 can create spurious divergence. |
| Energy applications | @herrerias2017residential; @tastan2023club | Precedents exist for regional and city electricity, but none at the neighbourhood scale. |

### Explanatory spatial and spatio-temporal models

| Method | Cite | Caveats |
|---|---|---|
| GWR / MGWR (cross-section, e.g. 2021) | @brunsdon1996geographically; @fotheringham2017multiscale; @oshan2019mgwr; @comber2023route; @feng2025examining | Cost is O(n²) per bandwidth, so ~35k LSOAs needs a parallel or chunked implementation, or MSOA level. Check local multicollinearity (local VIF/condition numbers) and correct local inference for multiple testing. Compare with a global spatial lag/error model as a baseline. |
| GTWR | @huang2010geographically; @fotheringham2015geographical | **If the covariates are Census 2021 snapshots, temporal variation in the coefficients is not identified from the covariates.** GTWR then mostly smooths the outcome. Use it only with time-varying covariates (EPC stock by year, degree days, meter counts), or model *change* on 2021 baselines instead. |
| Fixed Rank Kriging trend surface | @cressie2008fixed; @zammitmangion2021frk; @sainsburydale2024modeling | This is an R package; calling it from a Python page needs rpy2 or a separate R step. Basis-function resolution and the fine-scale variance component need reporting. It is best suited to smoothing and prediction rather than causal attribution. |

### Shock periods and covariates to add

| Issue | Cite | Implication |
|---|---|---|
| COVID-19 | @zapatawebborn2023impact; @bahmanyar2020impact | Electricity up +7.8% (year 1) and +2.2% (year 2), largest for households with children. The 2020 and 2021 meter years should carry indicators, and child share is a key interaction. |
| 2022 price crisis | @zapatawebborn2024winter; @fawcett2024using; @guan2023burden | Weather-adjusted electricity fell 8.4% in winter 2022/23. Reductions are tied to underheating, so lower consumption ≠ better outcome. |
| Weather | @summerfield2010two | Electricity is not weather-corrected. Add annual heating degree days by region or LSOA (for example from HadUK-Grid) as a covariate. |
| Electrification | @love2017addition; @arvanitopoulos2025decarbonising; @mckenna2018solar; @baltaozkan2015regional | Add heat-pump and PV install counts (MCS / feed-in tariff data) where available. They can produce convergence or divergence unrelated to need. |
| Heating fuel / off-gas grid | @roberts2015fuel | The Census 2021 central-heating type variable is essential. Electric-heating LSOAs form the upper tail. |
