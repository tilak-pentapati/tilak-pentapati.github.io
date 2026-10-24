#!/usr/bin/env bash
# Build the environment for one project page and print the command to render it.
#
#   scripts/repro/setup-env.sh projects/python/bd_gps
#   scripts/repro/setup-env.sh projects/r/stm_covid
#
# Uses the project's environment.yml (conda-forge, pinned) or, for the satellites
# page, requirements.txt; R pages with a cran-packages.txt also get those packages
# from the CRAN archive. Envs go to $REPRO_ENVS (default ~/.repro-envs/<project>).
# Works with micromamba, mamba or conda, whichever is on PATH. Raw data is not in
# the repository: see data/README.md for what to fetch and where to put it.
set -euo pipefail

proj="${1:?usage: $0 projects/<python|r>/<project>}"
proj="${proj%/}"
root="$(cd "$(dirname "$0")/../.." && pwd)"
dir="$root/$proj"
[ -f "$dir/index.qmd" ] || { echo "No index.qmd in $proj" >&2; exit 1; }

name="$(basename "$proj" | tr '[:upper:]' '[:lower:]')"
prefix="${REPRO_ENVS:-$HOME/.repro-envs}/$name"
solver="$(command -v micromamba || command -v mamba || command -v conda)" \
  || { echo "Need micromamba, mamba or conda on PATH" >&2; exit 1; }

if [ -f "$dir/environment.yml" ]; then
  "$solver" env create -y -p "$prefix" -f "$dir/environment.yml" 2>/dev/null \
    || "$solver" create -y -p "$prefix" -f "$dir/environment.yml"
elif [ -f "$dir/requirements.txt" ]; then
  "$solver" create -y -p "$prefix" -c conda-forge python=3.14.6 pip
  "$prefix/bin/python" -m pip install -r "$dir/requirements.txt"
else
  echo "No environment.yml or requirements.txt in $proj" >&2; exit 1
fi

if [ -f "$dir/cran-packages.txt" ]; then
  "$prefix/bin/Rscript" "$root/scripts/repro/install_cran.R" "$dir/cran-packages.txt"
fi

echo
echo "Environment ready: $prefix"
if [ -x "$prefix/bin/Rscript" ]; then
  echo "Render: QUARTO_R=$prefix/bin quarto render $proj/index.qmd --execute"
else
  echo "Render: QUARTO_PYTHON=$prefix/bin/python quarto render $proj/index.qmd --execute"
fi
echo "Data:   see the 'Data availability' section of the page and data/README.md"
