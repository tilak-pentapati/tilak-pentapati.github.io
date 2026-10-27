"""Copy the figures used on the project page into ``projects/python/energy_inequality/figures/``.

The analysis writes 300 dpi PNGs to the (gitignored) results folder. The page only
needs a web-sized copy of a few of them, so each one is downscaled to at most
1,600 px wide and palette-quantised. Re-run after ``run_dynamics.py`` /
``b_run_models.py`` refresh their figures:

    /Users/jay/anaconda3/envs/envs563/bin/python -I projects/python/energy_inequality/analysis/c_page_figures.py
"""

from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[4]
RESULTS = ROOT / "data" / "data_python" / "energy_inequality" / "results"
OUT = Path(__file__).resolve().parents[1] / "figures"
MAX_WIDTH = 1600

FIGURES = {
    "dynamics": [
        "esda_map_fj_2015_2024", "esda_kde_relative", "esda_moran_series", "esda_lisa_map_2024",
        "ineq_theil_gini", "ineq_cv_vs_mean_2024", "ineq_shock_2022_by_group",
        "markov_shorrocks_annual", "markov_spatial_lsoa_extremes",
        "rose_directional_lisa", "mobility_vankerm_annual", "rank_tau_annual", "sens_elec_gas_outliers",
    ],
    "models": [
        "clubs_paths_lad_full", "clubs_map_msoa_full", "clubs_profile_ordered_logit",
        "xs_lsoa_y2224_forest", "xs_mgwr_msoa_y2224_bandwidths",
        "xs_mgwr_msoa_y2224_map_share_ch_electric", "xs_mgwr_msoa_y2224_map_epc_ee",
        "ch_lsoa_long_forest", "ch_mgwr_msoa_long_bandwidths",
        "st_event_study", "st_trend_surface_years",
    ],
}


def thumbnail():
    """Listing thumbnail: LAD convergence clubs beside the core LSOA Markov matrix."""
    import geopandas as gpd
    import pandas as pd

    lad = gpd.read_parquet(RESULTS.parent / "processed" / "lad21_ew.parquet")[["lad_cd", "geometry"]]
    clubs = pd.read_csv(RESULTS / "models" / "tables" / "clubs_membership_lad.csv")[["lad_cd", "lad_full"]]
    lad = lad.merge(clubs, on="lad_cd")
    lad["geometry"] = lad.geometry.simplify(1500, preserve_topology=True)
    x0, y0, x1, y1 = lad.total_bounds
    scale = 520 / (y1 - y0)
    # 12 clubs, high = dark blue, low = pale; divergent (0) in amber
    ramp = ["#08306b", "#0b4a8c", "#1563a6", "#2479b8", "#3a8fc7", "#55a2d3",
            "#73b5dc", "#93c6e4", "#b2d6ec", "#cbe3f2", "#dfedf7", "#eef5fb"]

    def path(geom):
        polys = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
        d = []
        for p in polys:
            pts = [(40 + (x - x0) * scale, 60 + (y1 - y) * scale) for x, y in p.exterior.coords]
            d.append("M" + "L".join(f"{x:.1f},{y:.1f}" for x, y in pts) + "Z")
        return "".join(d)

    shapes = "\n".join(
        f'<path d="{path(g)}" fill="{"#f59e0b" if c == 0 else ramp[c - 1]}"/>'
        for g, c in zip(lad.geometry, lad.lad_full))
    mm = pd.read_parquet(RESULTS / "dynamics" / "tables" / "markov_matrices.parquet")
    P = mm.query("chain == 'lsoa 2015-24 (core)' and matrix == 'P'").pivot(index="from", columns="to", values="value")
    order = sorted(P.index, key=lambda s: s[:2])
    cells = []
    for i, r in enumerate(order):
        for j, c in enumerate(order):
            v = P.loc[r, c]
            cells.append(f'<rect x="{560 + j * 72}" y="{170 + i * 72}" width="68" height="68" rx="6" '
                         f'fill="#1563a6" fill-opacity="{0.06 + 0.94 * v:.2f}"/>')
            if v >= 0.05:
                cells.append(f'<text x="{594 + j * 72}" y="{210 + i * 72}" text-anchor="middle" font-size="17" '
                             f'fill="{"#ffffff" if v > 0.5 else "#0f172a"}">{v:.2f}</text>')
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 650" width="100%" height="100%" style="background-color: #f8fafc; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;">
<rect width="1000" height="650" fill="#f8fafc"/>
<g stroke="#f8fafc" stroke-width="0.6" stroke-linejoin="round">
{shapes}
</g>
<text x="560" y="120" font-size="26" font-weight="600" fill="#0f172a">Persistence of place</text>
<text x="560" y="150" font-size="16" fill="#475569">LSOA quintile transitions, 2015–24</text>
{"".join(cells)}
<text x="40" y="620" font-size="16" fill="#475569">Convergence clubs of 318 local authorities, 2010–24</text>
</svg>
"""
    dest = OUT.parent / "thumbnail.svg"
    dest.write_text(svg)
    print(f"{dest.stat().st_size / 1024:7.0f} KB  thumbnail.svg")


def main():
    OUT.mkdir(exist_ok=True)
    thumbnail()
    total = 0
    for part, names in FIGURES.items():
        for name in names:
            img = Image.open(RESULTS / part / "figures" / f"{name}.png").convert("RGB")
            if img.width > MAX_WIDTH:
                img = img.resize((MAX_WIDTH, round(img.height * MAX_WIDTH / img.width)), Image.LANCZOS)
            img = img.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
            dest = OUT / f"{name}.png"
            img.save(dest, optimize=True)
            total += dest.stat().st_size
            print(f"{dest.stat().st_size / 1024:7.0f} KB  {name}")
    print(f"{total / 1024 ** 2:.1f} MB in {OUT}")


if __name__ == "__main__":
    main()
