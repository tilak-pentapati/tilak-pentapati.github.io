# Reproducibility

This file is not rendered on the site (Quarto renders only `*.qmd`). It explains how
each project page can be re-run, and how to archive the repository on Zenodo for a
citable DOI.

## What is in the repository

| Layer | Where |
|---|---|
| Frozen computational results (the site renders from these; no data or envs needed) | `_freeze/` |
| Per-project environment, pinned to the local envs the freezes were made with | `projects/*/*/environment.yml` (conda-forge); `projects/python/satellites_overhead/requirements.txt`; `projects/r/*/cran-packages.txt` |
| Data sources, licences, how to fetch, expected paths | each page's "Data availability" section, and `data/README.md` |
| Dev container (GitHub Codespaces / VS Code) | `.devcontainer/devcontainer.json` |
| Binder | `binder/environment.yml` + `binder/postBuild` |
| Citation metadata | `CITATION.cff`, `.zenodo.json` |

Each project page has a **Reproduce** row in its summary card: *Codespaces*, *Binder*,
*Environment* (the pinned file on GitHub) and *Data* (jumps to the page's
"Data availability" section). The badge links point at the `main` branch, so they work
once this branch is merged.

## Environments

- **Python** (5 pages): `environment.yml` lists Python 3.14.6 plus only the packages the
  page imports, with the versions from the local `envs563` env, and the Quarto Jupyter
  toolchain (`ipykernel`, `nbclient`, `nbformat`, `pyyaml`). Packages that came from PyPI
  locally are in its `pip:` section. The satellites page keeps its `requirements.txt`,
  which the nightly CI job installs.
- **R** (6 pages): `environment.yml` pins `r-base=4.5.3` and the conda-forge build of each
  package the page loads, at the version in the local `r_env_gds` env, plus `knitr` and
  `rmarkdown`. Three packages are not on conda-forge (`FRK`, `spgwr`, `tidyterra`); they are
  pinned in `cran-packages.txt` and installed from the CRAN archive by
  `scripts/repro/install_cran.R`.

  Why conda rather than `renv.lock`: the R that rendered the freezes is itself a conda env
  (conda-forge `r-base` with conda-forge builds of most packages, plus a few CRAN installs
  into the same library). A conda spec reproduces that, including the compiled GDAL / GEOS /
  PROJ / udunits stack that `sf`, `terra` and `tmap` link against, as binaries on Linux,
  macOS and Windows, and the same micromamba/conda tooling builds the Python and R pages in
  Codespaces and Binder. `renv` is not installed in `r_env_gds`, a `renv.lock` would record
  conda builds as CRAN sources (restoring them means compiling from source against system
  libraries), and an renv project changes `.Rprofile` start-up behaviour site-wide. The files
  are still a plain pinned list (`Package=version` plus the R version), so non-conda users
  can install the same versions with `remotes::install_version()`.
- Four pages have no reproducibility layer (no environment file, badges or data notes):
  `bd_gps` (Microsoft Research GeoLife licence, non-commercial research use) and
  `StrokePredictionModel` (Kaggle data with no clear licence) because of their data
  licences, and `GCH_census` and `aps_eda` because they need `ENVS450`, an unpublished
  University of Liverpool course package, and course data owned by others.

Regenerate the specs after changing a page's imports or upgrading the local envs:

```bash
/Users/jay/anaconda3/envs/envs563/bin/python scripts/build_env_specs.py
/Users/jay/anaconda3/envs/envs563/bin/python scripts/build_data_readme.py
```

(edit the import lists at the top of `build_env_specs.py`; data notes live in each
project's `_metadata.yml` under `reproduce:`).

## Re-running a page

In a Codespace (or locally with micromamba/mamba/conda and Quarto 1.10.18):

```bash
scripts/repro/setup-env.sh projects/r/stm_covid
```

It builds the env under `~/.repro-envs/<project>` and prints the render command, e.g.
`QUARTO_R=~/.repro-envs/stm_covid/bin quarto render projects/r/stm_covid/index.qmd --execute`.
Fetch the page's data first (see its "Data availability" section) into the expected path
under `data/`; the pages with Reproduce badges read it relative to the project root.

On Binder every page's packages are preinstalled in one environment (the first build takes
a while); open a terminal and run `quarto render <page> --execute`.

## Citable releases with Zenodo (manual steps)

Zenodo archives a GitHub release and gives it a DOI. These steps need your accounts, so do
them yourself once this branch is merged to `main`:

1. **Choose a licence first** (the TODO in `projects/_metadata.yml` and `CITATION.cff`);
   Zenodo shows it on the record. Add it as `"license": "cc-by-4.0"` (or your choice) to
   `.zenodo.json`, `license:` (SPDX id) to `CITATION.cff`, and `license:` to
   `projects/_metadata.yml`.
2. **ORCID (optional, recommended):** once you have one, add
   `"orcid": "0000-0000-0000-0000"` to the creator in `.zenodo.json`, uncomment `orcid:` in
   `CITATION.cff`, and uncomment `orcid:` in `projects/_metadata.yml`.
3. Log in to <https://zenodo.org> with **"Log in with GitHub"** and authorise Zenodo.
4. Open **Account → GitHub** (<https://zenodo.org/account/settings/github/>), click
   **Sync now**, and switch on `tilak-pentapati/tilak-pentapati.github.io`.
5. On GitHub, **create a release**: Releases → *Draft a new release* → tag e.g. `v1.0.0` on
   `main` → *Publish release*. Zenodo archives it within a few minutes (it reads
   `.zenodo.json`). Check the release appears under Account → GitHub.
6. Zenodo shows two DOIs: a **concept DOI** (always resolves to the latest version) and a
   version DOI. Use the **concept DOI** below.

### Where to paste the DOI

- `projects/_metadata.yml`: add `doi: 10.5281/zenodo.NNNNNNN`. Quarto puts it in each
  project page's "Cite this" appendix (BibTeX + citation) and Google Scholar tags, and the
  summary-card filter then shows a **DOI** badge in every page's Reproduce row
  automatically. `_metadata.yml` is outside the freeze hash, so a plain `quarto render`
  picks this up without re-executing anything.
- `CITATION.cff`: uncomment `doi:` and paste the same DOI.
- `README` / GitHub repo description (optional): the Zenodo badge Markdown from the record page.

Never type a DOI before Zenodo has issued it.
