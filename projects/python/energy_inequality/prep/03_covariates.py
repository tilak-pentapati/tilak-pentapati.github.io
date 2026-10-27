"""Cross-sectional covariates on 2021 LSOAs (E&W) and a mid-year population panel.

Outputs (processed/)
- covariates_lsoa21.parquet      one row per LSOA21 (35,672)
- population_lsoa_panel.parquet  ONS mid-year estimates on 2021 LSOAs, 2011-2024 (2010 not published on 2021 LSOAs)

Methods
- Census 2021 (Nomis bulk TS tables) are published on 2021 LSOAs: used directly.
  Shares use each table's own total as denominator (households, usual residents or
  residents aged 16+ as published).
- IMD 2025 (England) and WIMD 2025 (Wales) are published on 2021 LSOAs: used directly.
- IMD 2019 / WIMD 2019 scores are on 2011 LSOAs. Each LSOA21 gets the household-weighted
  average of the scores of its LSOA11 sources (pair households from 01_geography.py):
  unchanged and split LSOAs simply inherit their parent's score, merged / complex ones
  get a weighted average. Ranks and deciles are then recomputed within each country on
  the converted scores, so they are NOT the published ranks for changed areas.
- EPC (ONS energy efficiency of housing): the FYE2025 edition is on 2021 MSOAs and is
  assigned to every LSOA21 in the MSOA (LSOA21 nest exactly in MSOA21). The FYE2021 and
  FYE2022 editions (local copies of earlier releases) are on 2011 MSOAs and are moved
  to 2021 MSOAs with household-weighted averages built from the same LSOA pair weights.
  EPC values are MSOA-level medians/percentages of dwellings WITH an EPC, not LSOA values.
"""
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # keep __pycache__ out of the repo
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from _common import AUDIT, PROCESSED, RAW, log, read_census, read_lookup, safe_div  # noqa: E402

DEP = RAW / "deprivation"
EPC = RAW / "epc_ons"
DISS = RAW / "superseded_local"


# ---------------------------------------------------------------------------- Census
def share(t: pd.DataFrame, total: str, parts: dict[str, list[str]]) -> pd.DataFrame:
    out = {}
    for name, cols in parts.items():
        missing = [c for c in cols if c not in t.columns]
        if missing:
            raise KeyError(f"{name}: {missing}")
        out[name] = safe_div(t[cols].sum(axis=1), t[total])
    return pd.DataFrame(out, index=t.index)


def census() -> pd.DataFrame:
    frames = []
    t = read_census("ts001")
    frames.append(pd.DataFrame({"usual_residents_2021": t["Total"],
                                "share_communal": safe_div(t["Lives in a communal establishment"], t["Total"])}, index=t.index))
    t = read_census("ts041")
    frames.append(t.rename(columns={"Number of households": "households_2021"}))

    t = read_census("ts017")
    sizes = {f"{k} people in household": k for k in range(2, 8)} | {"1 person in household": 1, "8 or more people in household": 8}
    occupied = t[list(sizes)].sum(axis=1)
    frames.append(pd.DataFrame({
        "mean_hh_size": safe_div(sum(t[c] * k for c, k in sizes.items()), occupied),
        "share_1person_hh": safe_div(t["1 person in household"], occupied),
        "share_5plus_hh": safe_div(t[["5 people in household", "6 people in household", "7 people in household",
                                      "8 or more people in household"]].sum(axis=1), occupied),
    }, index=t.index))

    t = read_census("ts044")
    frames.append(share(t, "Total: All households", {
        "share_detached": ["Detached"], "share_semi": ["Semi-detached"], "share_terraced": ["Terraced"],
        "share_flat": ["In a purpose-built block of flats or tenement",
                       "Part of a converted or shared house, including bedsits",
                       "Part of another converted building, for example, former school, church or warehouse",
                       "In a commercial building, for example, in an office building, hotel or over a shop"],
        "share_caravan": ["A caravan or other mobile or temporary structure"],
    }))
    t = read_census("ts054")
    frames.append(share(t, "Total: All households", {
        "share_owned_outright": ["Owned: Owns outright"], "share_owned_mortgage": ["Owned: Owns with a mortgage or loan"],
        "share_shared_ownership": ["Shared ownership: Shared ownership"], "share_social_rent": ["Social rented"],
        "share_private_rent": ["Private rented", "Lives rent free"],
    }))
    t = read_census("ts046")
    frames.append(share(t, "Total: All households", {
        "share_no_central_heating": ["No central heating"], "share_ch_mains_gas": ["Mains gas only"],
        "share_ch_electric": ["Electric only"], "share_ch_oil": ["Oil only"],
        "share_ch_lpg": ["Tank or bottled gas only"], "share_ch_solid_wood": ["Wood only", "Solid fuel only"],
        "share_ch_renewable": ["Renewable energy only"], "share_ch_district": ["District or communal heat networks only"],
        "share_ch_other": ["Other central heating only"],
        "share_ch_multi": ["Two or more types of central heating (not including renewable energy)",
                           "Two or more types of central heating (including renewable energy)"],
    }))
    t = read_census("ts052")
    frames.append(share(t, "Total: All households", {
        "share_overcrowded_bedrooms": ["Occupancy rating of bedrooms: -1", "Occupancy rating of bedrooms: -2 or less"],
        "share_underoccupied_bedrooms": ["Occupancy rating of bedrooms: +2 or more"],
    }))
    t = read_census("ts053")
    frames.append(share(t, "Total: All households", {
        "share_overcrowded_rooms": ["Occupancy rating of rooms: -1", "Occupancy rating of rooms: -2 or less"],
    }))
    t = read_census("ts062")
    tot = "Total: All usual residents aged 16 years and over"
    frames.append(share(t, tot, {
        "share_nssec_higher_mgr_prof": ["L1, L2 and L3 Higher managerial, administrative and professional occupations"],
        "share_nssec_lower_mgr_prof": ["L4, L5 and L6 Lower managerial, administrative and professional occupations"],
        "share_nssec_intermediate": ["L7 Intermediate occupations"],
        "share_nssec_small_employer": ["L8 and L9 Small employers and own account workers"],
        "share_nssec_lower_supervisory": ["L10 and L11 Lower supervisory and technical occupations"],
        "share_nssec_semi_routine": ["L12 Semi-routine occupations"],
        "share_nssec_routine": ["L13 Routine occupations"],
        "share_nssec_never_worked_ltu": ["L14.1 and L14.2 Never worked and long-term unemployed"],
        "share_nssec_student": ["L15 Full-time students"],
    }))
    t = read_census("ts021")
    frames.append(share(t, "Total: All usual residents", {
        "share_eth_asian": ["Asian, Asian British or Asian Welsh"],
        "share_eth_black": ["Black, Black British, Black Welsh, Caribbean or African"],
        "share_eth_mixed": ["Mixed or Multiple ethnic groups"], "share_eth_white": ["White"],
        "share_eth_white_british": ["White: English, Welsh, Scottish, Northern Irish or British"],
        "share_eth_other": ["Other ethnic group"],
    }))
    t = read_census("ts007a")
    frames.append(share(t, "Total", {
        "share_age_0_14": ["Aged 4 years and under", "Aged 5 to 9 years", "Aged 10 to 14 years"],
        "share_age_15_24": ["Aged 15 to 19 years", "Aged 20 to 24 years"],
        "share_age_65plus": ["Aged 65 to 69 years", "Aged 70 to 74 years", "Aged 75 to 79 years",
                             "Aged 80 to 84 years", "Aged 85 years and over"],
        "share_age_75plus": ["Aged 75 to 79 years", "Aged 80 to 84 years", "Aged 85 years and over"],
    }))
    t = read_census("ts050")
    beds = {"1 bedroom": 1, "2 bedrooms": 2, "3 bedrooms": 3, "4 or more bedrooms": 4}
    frames.append(pd.DataFrame({
        "mean_bedrooms": safe_div(sum(t[c] * k for c, k in beds.items()), t["Total: All households"]),
        "share_4plus_bedrooms": safe_div(t["4 or more bedrooms"], t["Total: All households"]),
    }, index=t.index))
    t = read_census("ts051")
    rooms = {f"{k} rooms": k for k in range(2, 9)} | {"1 room": 1, "9 or more rooms": 9}
    frames.append(pd.DataFrame({"mean_rooms": safe_div(sum(t[c] * k for c, k in rooms.items()), t["Total: All households"])},
                               index=t.index))
    t = read_census("ts045")
    frames.append(share(t, "Total: All households", {
        "share_no_car": ["No cars or vans in household"],
        "share_2plus_cars": ["2 cars or vans in household", "3 or more cars or vans in household"],
    }))
    t = read_census("ts066")
    tot = "Total: All usual residents aged 16 years and over"
    emp = ["Economically active (excluding full-time students):In employment",
           "Economically active and a full-time student:In employment"]
    unemp = ["Economically active (excluding full-time students): Unemployed",
             "Economically active and a full-time student: Unemployed"]
    frames.append(share(t, tot, {
        "share_employed": emp, "share_unemployed": unemp,
        "share_inactive": ["Economically inactive"], "share_retired": ["Economically inactive: Retired"],
        "share_inactive_student": ["Economically inactive: Student"],
        "share_long_term_sick": ["Economically inactive: Long-term sick or disabled"],
        "share_looking_after_home": ["Economically inactive: Looking after home or family"],
    }))
    act = t[["Economically active (excluding full-time students)", "Economically active and a full-time student"]].sum(axis=1)
    frames.append(pd.DataFrame({"unemployment_rate": safe_div(t[unemp].sum(axis=1), act)}, index=t.index))
    c = pd.concat(frames, axis=1)
    c.index.name = "lsoa21cd"
    return c


# ---------------------------------------------------------------------- Deprivation
def rerank(score: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Rank 1 = most deprived (highest score); decile 1 = most deprived 10%."""
    r = score.rank(ascending=False, method="min")
    dec = np.ceil(r / r.max() * 10).clip(1, 10)
    return r, dec


def to_lsoa21_weighted(vals: pd.DataFrame, w: pd.DataFrame) -> pd.DataFrame:
    """Household-weighted average of LSOA11 values for each LSOA21 (vals indexed by lsoa11cd)."""
    x = w[["lsoa11cd", "lsoa21cd", "hh_pair"]].merge(vals, left_on="lsoa11cd", right_index=True, how="inner")
    x["hh_pair"] = x["hh_pair"].clip(lower=1e-9)  # keep zero-household pairs from vanishing entirely
    cols = list(vals.columns)
    num = x[cols].mul(x["hh_pair"], axis=0).groupby(x["lsoa21cd"]).sum(min_count=1)
    den = x[cols].notna().mul(x["hh_pair"], axis=0).groupby(x["lsoa21cd"]).sum()
    return num / den.replace(0, np.nan)


def imd2025() -> pd.DataFrame:
    d = pd.read_csv(DEP / "File_7_IoD2025_All_Ranks_Scores_Deciles_Population_Denominators.csv")
    keep = {
        "LSOA code (2021)": "lsoa21cd",
        "Index of Multiple Deprivation (IMD) Score": "imd25_score",
        "Index of Multiple Deprivation (IMD) Rank (where 1 is most deprived)": "imd25_rank",
        "Index of Multiple Deprivation (IMD) Decile (where 1 is most deprived 10% of LSOAs)": "imd25_decile",
        "Income Score (rate)": "imd25_income_score",
        "Employment Score (rate)": "imd25_employment_score",
        "Living Environment Score": "imd25_living_env_score",
    }
    indoors = [c for c in d.columns if c.startswith("Indoors Sub-domain Score")]
    if indoors:
        keep[indoors[0]] = "imd25_indoors_score"
    d = d[list(keep)].rename(columns=keep).set_index("lsoa21cd")
    return d


def wimd2025() -> pd.DataFrame:
    d = pd.read_csv(DEP / "WIMD2025_LSOA_ranks_statswales.csv", encoding="utf-8-sig")
    d = d[d["Area code"].str.startswith("W01")]
    piv = d.pivot_table(index="Area code", columns=["Data description", "Domain"], values="Data values", aggfunc="first")
    out = pd.DataFrame(index=piv.index)
    out["wimd25_rank"] = pd.to_numeric(piv[("Rank", "WIMD")])
    out["wimd25_decile"] = pd.to_numeric(piv[("Decile", "WIMD")], errors="coerce")
    for dom, name in [("Income", "income"), ("Employment", "employment"), ("Housing", "housing"),
                      ("Physical environment", "physical_env")]:
        out[f"wimd25_{name}_rank"] = pd.to_numeric(piv[("Rank", dom)])
    out.index.name = "lsoa21cd"
    return out


def imd2019_converted(w: pd.DataFrame) -> pd.DataFrame:
    d = pd.read_csv(DEP / "File_7_-_All_IoD2019_Scores__Ranks__Deciles_and_Population_Denominators_3.csv")
    d = d.set_index("LSOA code (2011)")[["Index of Multiple Deprivation (IMD) Score", "Income Score (rate)"]]
    d.columns = ["imd19_score", "imd19_income_score"]
    e = to_lsoa21_weighted(d, w)
    e = e[e.index.str.startswith("E")]
    e["imd19_rank_reranked"], e["imd19_decile_reranked"] = rerank(e["imd19_score"])
    return e


def wimd2019_converted(w: pd.DataFrame) -> pd.DataFrame:
    s = pd.read_excel(DEP / "wimd-2019-index-and-domain-scores-by-small-area.ods", sheet_name="Data", header=None,
                      engine="calamine")
    h = s.index[s[0].astype(str).str.strip() == "LSOA code"][0]
    s.columns = [str(c).strip() for c in s.iloc[h]]
    s = s.iloc[h + 1:]
    s = s[s["LSOA code"].astype(str).str.startswith("W01")].set_index("LSOA code")
    d = pd.DataFrame({"wimd19_score": pd.to_numeric(s["WIMD 2019"]), "wimd19_income_score": pd.to_numeric(s["Income"])})
    e = to_lsoa21_weighted(d, w)
    e["wimd19_rank_reranked"], e["wimd19_decile_reranked"] = rerank(e["wimd19_score"])
    return e


# ------------------------------------------------------------------------------- EPC
def read_epc_sheet(path: Path, sheet: str, prefix: str) -> pd.DataFrame:
    s = pd.read_excel(path, sheet_name=sheet, header=None, engine="calamine")
    is_hdr = s.astype(str).apply(lambda c: c.str.contains("MSOA") & c.str.contains("code|CD", regex=True)).any(axis=1)
    h = is_hdr.idxmax()
    hdr = [" ".join(str(v).split()) for v in s.iloc[h]]
    mcol = next(i for i, c in enumerate(hdr) if "MSOA" in c and ("code" in c or c.endswith("CD")))
    body = s.iloc[h + 1:].copy()
    body = body[body[mcol].astype(str).str.match(r"[EW]02")]
    vals = {}
    for i, c in enumerate(hdr):
        if i <= mcol + 1 or c in ("nan", "None") or not c:
            continue
        name = prefix + "_" + "".join(ch if ch.isalnum() else "_" for ch in c.lower()).strip("_")
        while "__" in name:
            name = name.replace("__", "_")
        vals[name] = pd.to_numeric(body[i], errors="coerce").values
    return pd.DataFrame(vals, index=body[mcol].astype(str).values)


def msoa11_to_msoa21_weights(w: pd.DataFrame) -> pd.DataFrame:
    l11 = read_lookup("OA11_LSOA11_MSOA11_LAD11_EW_LU.csv")[["LSOA11CD", "MSOA11CD"]].drop_duplicates()
    l21 = read_lookup("OA21_LSOA21_MSOA21_LAD22_EW_LU_v3.csv")[["LSOA21CD", "MSOA21CD"]].drop_duplicates()
    x = w.merge(l11, left_on="lsoa11cd", right_on="LSOA11CD").merge(l21, left_on="lsoa21cd", right_on="LSOA21CD")
    m = x.groupby(["MSOA11CD", "MSOA21CD"], as_index=False)["hh_pair"].sum()
    return m.rename(columns={"MSOA11CD": "lsoa11cd", "MSOA21CD": "lsoa21cd"})  # reuse to_lsoa21_weighted


def epc() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (MSOA21-level EPC table, audit of coverage)."""
    f25 = [
        read_epc_sheet(EPC / "median_ee_score_FYE2025.xlsx", "3a", "epc25_median_ee"),
        read_epc_sheet(EPC / "median_ee_score_FYE2025.xlsx", "3b", "epc25_median_ee"),
        read_epc_sheet(EPC / "epc_band_c_or_above_FYE2025.xlsx", "3a", "epc25_pct_band_c_plus"),
        read_epc_sheet(EPC / "main_fuel_type_FYE2025.xlsx", "3a", "epc25_pct_fuel"),
        read_epc_sheet(EPC / "median_co2_FYE2025.xlsx", "3a", "epc25_median_co2_t"),
    ]
    e25 = pd.concat(f25, axis=1)
    e25 = e25.loc[:, ~e25.columns.duplicated()]

    mw = None
    olds = []
    for yr, f in [("21", "Energy_Eff_2021.xlsx"), ("22", "Energy_Eff_2022.xlsx")]:
        t = pd.concat([read_epc_sheet(DISS / f, "1a", f"epc{yr}_median_ee"),
                       read_epc_sheet(DISS / f, "1b", f"epc{yr}_pct_band_c_plus")], axis=1)
        t.columns = [f"epc{yr}_median_ee_all", f"epc{yr}_pct_band_c_plus_all"]
        olds.append(t)
    if olds:
        w = pd.read_parquet(PROCESSED / "lsoa11_lsoa21_weights.parquet")
        mw = msoa11_to_msoa21_weights(w)
        old = to_lsoa21_weighted(pd.concat(olds, axis=1), mw)
        e = e25.join(old, how="outer")
    else:
        e = e25
    e.index.name = "msoa21cd"
    audit = pd.DataFrame({"non_missing": e.notna().sum(), "msoa_rows": len(e)})
    return e, audit


# ------------------------------------------------------------------------ Population
def population() -> pd.DataFrame:
    out = []
    for f in ["sapelsoabroadage20112022.xlsx", "sapelsoabroadage20222024.xlsx"]:
        sheets = pd.read_excel(RAW / "population" / f, sheet_name=None, header=None, engine="calamine")
        for name, s in sheets.items():
            if not (name.startswith("Mid-") and name.endswith("LSOA 2021")):
                continue
            year = int(name[4:8])
            h = s.index[s[2].astype(str).str.strip() == "LSOA 2021 Code"][0]
            s.columns = [str(c).strip() for c in s.iloc[h]]
            b = s.iloc[h + 1:]
            b = b[b["LSOA 2021 Code"].astype(str).str.match(r"[EW]01")]
            agecols = [c for c in b.columns if c[:1] in ("F", "M") and c[1:2].isdigit()]
            v = b[agecols].apply(pd.to_numeric, errors="coerce")
            grp = lambda pat: v[[c for c in agecols if c[1:] == pat]].sum(axis=1)  # noqa: E731
            out.append(pd.DataFrame({
                "lsoa21cd": b["LSOA 2021 Code"].values, "year": year, "source_file": f,
                "pop_total": v.sum(axis=1).values, "pop_0_15": grp("0 to 15").values,
                "pop_65plus": grp("65 and over").values,
            }))
    p = pd.concat(out, ignore_index=True)
    # Mid-2022 appears in both files; keep the Nov 2025 revision.
    p = p.sort_values("source_file").drop_duplicates(["lsoa21cd", "year"], keep="last")
    p["pop_16_64"] = p["pop_total"] - p["pop_0_15"] - p["pop_65plus"]
    return p.drop(columns="source_file").sort_values(["lsoa21cd", "year"]).reset_index(drop=True)


def main() -> None:
    lsoa = pd.read_parquet(PROCESSED / "lsoa21_ew.parquet",
                           columns=["lsoa21cd", "msoa21cd", "lad_cd", "country", "area_km2"]).set_index("lsoa21cd")
    w = pd.read_parquet(PROCESSED / "lsoa11_lsoa21_weights.parquet")

    cen = census()
    assert set(cen.index) >= set(lsoa.index), "Census missing LSOAs"
    cov = lsoa.join(cen, how="left")
    cov["pop_density_2021_per_km2"] = cov["usual_residents_2021"] / cov["area_km2"]

    ruc = read_lookup("LSOA21_RUC21_EW_LU.csv").set_index("LSOA21CD")[["RUC21CD", "RUC21NM", "Urban_rural_flag"]]
    cov = cov.join(ruc.rename(columns={"RUC21CD": "ruc21_cd", "RUC21NM": "ruc21_nm", "Urban_rural_flag": "urban_rural"}))

    cov = cov.join(imd2025()).join(wimd2025()).join(imd2019_converted(w)).join(wimd2019_converted(w))
    # Within-country deprivation percentile (0 = most deprived) to pool England and Wales cautiously.
    n = cov.groupby("country")["country"].transform("size")
    for tag, e, wl in [("2025", "imd25_rank", "wimd25_rank"), ("2019", "imd19_rank_reranked", "wimd19_rank_reranked")]:
        r = cov[e].fillna(cov[wl])
        cov[f"dep{tag}_pct_within_country"] = (r - 1) / (n - 1)
    chg = w.drop_duplicates("lsoa21cd").set_index("lsoa21cd")["lsoa21_change"]
    cov["lsoa21_change"] = chg
    cov["dep2019_converted_approx"] = cov["lsoa21_change"] != "U"

    e, epc_audit = epc()
    cov = cov.join(e, on="msoa21cd")
    epc_audit.to_csv(AUDIT / "epc_coverage.csv")

    cov = cov.reset_index()
    cov.to_parquet(PROCESSED / "covariates_lsoa21.parquet", index=False)
    miss = cov.isna().sum()
    miss.rename("missing").to_csv(AUDIT / "covariates_missing.csv", header=True)
    log(f"covariates_lsoa21: {cov.shape[0]:,} rows x {cov.shape[1]} cols")
    log("columns with missing values:\n" + miss[miss > 0].to_string())

    pop = population()
    assert pop["lsoa21cd"].nunique() == len(lsoa), "population LSOAs != 35,672"
    pop.to_parquet(PROCESSED / "population_lsoa_panel.parquet", index=False)
    log(f"population_lsoa_panel: {len(pop):,} rows, years {pop['year'].min()}-{pop['year'].max()}")
    log(pop.groupby("year")["pop_total"].sum().map("{:,.0f}".format).to_string())


if __name__ == "__main__":
    main()
