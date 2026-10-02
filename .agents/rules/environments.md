---
trigger: always_on
---

# Environment Execution Rules

## Script Execution (Direct Binaries)
- For Python scripts: `/Users/jay/anaconda3/envs/envs563/bin/python <script.py>`
- For R scripts: `/Users/jay/anaconda3/envs/r_env_gds/bin/Rscript <script.R>`

## Quarto Rendering (Explicit Overrides)
- To render Python Quarto projects: `QUARTO_PYTHON=/Users/jay/anaconda3/envs/envs563/bin/python quarto render <target>`
- To render R Quarto projects: `QUARTO_R=/Users/jay/anaconda3/envs/r_env_gds/bin quarto render <target>`

## Constraints
- Never run bare `python`, `Rscript`, or `quarto render` commands without absolute paths or variables.
- Never use `conda run` or `conda activate`.