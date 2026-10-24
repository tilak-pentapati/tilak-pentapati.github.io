# Install pinned CRAN packages (lines "Package==version") that are not on
# conda-forge into the R library of the env running this script.
#   Rscript scripts/repro/install_cran.R projects/r/stm_covid/cran-packages.txt [...]
# Packages already installed at the pinned version are skipped.
files <- commandArgs(trailingOnly = TRUE)
if (!length(files)) stop("usage: Rscript install_cran.R <cran-packages.txt> [...]")
for (f in files) {
  lines <- trimws(sub("#.*", "", readLines(f)))
  for (spec in lines[nzchar(lines)]) {
    parts <- strsplit(spec, "==", fixed = TRUE)[[1]]
    pkg <- parts[1]; ver <- parts[2]
    have <- tryCatch(packageDescription(pkg)$Version, warning = function(w) NA)
    if (identical(have, ver)) { message(pkg, " ", ver, " already installed"); next }
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    message("Installing ", pkg, " ", ver, " from the CRAN archive")
    remotes::install_version(pkg, version = ver, repos = "https://cloud.r-project.org",
                             upgrade = "never", quiet = TRUE)
  }
}
