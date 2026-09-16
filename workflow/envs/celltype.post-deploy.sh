#!/usr/bin/env bash
# Runs once, inside the new env, after Snakemake creates workflow/envs/celltype.yaml.
# Needs internet access, which the login node (where the controller runs) has.
set -euo pipefail

Rscript - <<'RS'
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 3600)
gh <- c("immunogenomics/presto", "satijalab/seurat-data",
        "mojaveazure/seurat-disk", "satijalab/azimuth")
for (r in gh) remotes::install_github(r, upgrade = "never", dependencies = FALSE)

# Azimuth reference (Allen human motor cortex), ~0.5 GB, into the env's library.
SeuratData::InstallData("humancortexref")

for (p in c("Azimuth", "SeuratData", "humancortexref.SeuratData", "scDblFinder", "harmony")) {
  stopifnot(requireNamespace(p, quietly = TRUE))
}
RS
