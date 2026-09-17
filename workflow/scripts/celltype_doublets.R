#!/usr/bin/env Rscript
# =============================================================================
# Cell typing, stage 2 of 3: scDblFinder on ONE sublibrary.
#
#   <prep>/counts/<sublibrary>.rds  ->  <doublets>/<sublibrary>.tsv
#                                       (bc_wells, doublet_score, doublet_class)
#
# One sublibrary per job, because that is where Parse doublets live: two nuclei
# become one barcode only if they share a round 1-3 combination AND land in the
# same sublibrary (the round-4 index then tells the rest apart). scDblFinder
# builds its artificial doublets from cells of the same batch, and each
# sublibrary also has its own depth (ParseFANS15/16 are much shallower).
#
# Narrow on purpose: 16 small jobs queue far faster than one wide one. BLAS /
# OpenMP (xgboost) are pinned to the job's slots so it uses what it reserved.
# =============================================================================

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  dirname(normalizePath(f))
}
source(file.path(script_dir(), "celltype_common.R"))

need_args(c("counts", "out", "threads", "seed", "dbr", "dbr_sd"))
pin_threads(as.integer(num("threads")))
seed <- as.integer(num("seed"))
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)

m <- readRDS(opt$counts)
step("read ", opt$counts, ": ", nrow(m), " genes x ", ncol(m), " cells")

# scDblFinder cannot fit a model on a handful of cells (only possible in a
# very small smoke test); leave those uncalled rather than fail the job.
min_cells <- 200
if (ncol(m) < min_cells) {
  warning(ncol(m), " cells is too few for scDblFinder (< ", min_cells, "); doublets not called")
  res <- data.table(bc_wells = colnames(m), doublet_score = NA_real_, doublet_class = NA_character_)
} else {
  # scDblFinder's own rate (0.8% per 1000 cells) is a 10x figure; the
  # configured rate is the Parse collision rate -- see config.yaml.
  set.seed(seed)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
  rm(m); invisible(gc())
  sce <- scDblFinder::scDblFinder(sce, dbr = num_or_null("dbr"), dbr.sd = num_or_null("dbr_sd"),
                                  BPPARAM = BiocParallel::SerialParam(RNGseed = seed),
                                  verbose = TRUE)
  res <- data.table(bc_wells      = colnames(sce),
                    doublet_score = sce$scDblFinder.score,
                    doublet_class = as.character(sce$scDblFinder.class))
}

fwrite(res, opt$out, sep = "\t")
step(sum(res$doublet_class %in% "doublet"), " doublets called (",
     round(100 * mean(res$doublet_class %in% "doublet"), 2), "%)")
log_peak()
