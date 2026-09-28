#!/usr/bin/env Rscript
# =============================================================================
# Ambient RNA, step 3: SoupX-corrected counts for ONE sublibrary.
#
#   <run>/prep/counts/<sublib>.rds     raw counts (genes x cells), every cell
#                                      that passed the low-quality filters
#   <run>/<labels>                     frozen labels: sort + cluster per cell
#   <run>/soupx/rho/<sort>.tsv         one rho per sort (soupx_rho.R)
#   soup/<sort>/soup_profile.tsv.gz    the measured soup
#     -> per sort within the sublibrary: SoupChannel + setSoupProfile +
#        setClusters + setContaminationFraction + adjustCounts(roundToInt)
#     -> <run>/soupx/counts/<sublib>.rds   corrected counts, same layout and
#                                          cell order as prep/counts
#        <run>/soupx/cells/<sublib>.tsv    per cell: sort, rho, cluster used,
#                                          UMIs before / after
#
# Every prep cell is corrected, so the cluster stage can apply its QC and
# doublet filters to the corrected counts exactly as it does to the raw ones.
# Cells that were clustered in the frozen run are corrected at cluster level
# (SoupX's recommended mode). Cells that were not -- unassigned donors, old
# doublet calls, upper-bound QC fails -- have no cluster, and are corrected
# at cell level, ADJ_CHUNK cells per call (see the loop for why).
#
# One SoupChannel per sort: each sort has its own soup and its own rho.
# roundToInt uses stochastic rounding, so the seed is set per sort.
# =============================================================================

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  dirname(normalizePath(f))
}
source(file.path(script_dir(), "celltype_common.R"))
suppressPackageStartupMessages({
  library(Matrix)
  library(SoupX)
})

need_args(c("counts", "labels", "sublib", "rho_dir", "soup_dir", "rho_override",
            "out_counts", "out_cells", "threads", "seed"))
setup_threads(as.integer(num("threads")))
seed <- as.integer(num("seed"))
dir.create(dirname(opt$out_counts), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$out_cells), recursive = TRUE, showWarnings = FALSE)
step("SoupX ", as.character(packageVersion("SoupX")), "  sublibrary=", opt$sublib)

m <- readRDS(opt$counts)
step("read ", opt$counts, ": ", nrow(m), " genes x ", ncol(m), " cells")

lab <- fread(opt$labels, select = c("bc_wells", "sample", "clustered", "seurat_clusters"),
             colClasses = c(seurat_clusters = "character"))
idx <- match(colnames(m), lab$bc_wells)
if (anyNA(idx)) stop(sum(is.na(idx)), " cells in the counts are missing from ", opt$labels)
cell_sort <- lab$sample[idx]
clustered <- lab$clustered[idx] %in% TRUE
cluster   <- ifelse(clustered, paste0("c", lab$seurat_clusters[idx]), paste0("cell_", colnames(m)))
rm(lab); invisible(gc())

# rho per sort: soupx_rho's estimate, unless overridden in config.
sorts <- unique(cell_sort)
rho <- vapply(sorts, function(s) fread(file.path(opt$rho_dir, paste0(s, ".tsv")))$rho, numeric(1))
if (!tolower(opt$rho_override) %in% c("", "none", "null")) {
  ov <- strsplit(csv("rho_override"), ":", fixed = TRUE)
  ov <- setNames(as.numeric(vapply(ov, `[`, "", 2)), vapply(ov, `[`, "", 1))
  step("rho overridden for: ", paste(names(ov), collapse = ", "))
  rho[intersect(names(ov), sorts)] <- ov[intersect(names(ov), sorts)]
}
print(rho)
stopifnot(!anyNA(rho))

# =============================================================================
# SoupX, one channel per sort
# =============================================================================
ADJ_CHUNK <- 500L  # unclustered cells per cell-level adjustCounts call

# One SoupChannel + adjustCounts. toc stands in for tod: the soup is
# supplied, not estimated from it. clusters = FALSE -> cell level.
soupx_run <- function(toc, soup_df, rho, clusters) {
  sc <- SoupChannel(toc, toc, calcSoupProfile = FALSE)
  sc <- setSoupProfile(sc, soup_df)
  if (!isFALSE(clusters)) sc <- setClusters(sc, clusters)
  sc <- setContaminationFraction(sc, rho, forceAccept = TRUE)
  adj <- adjustCounts(sc, clusters = if (isFALSE(clusters)) FALSE else NULL,
                      roundToInt = TRUE, verbose = 0)
  as(adj, "CsparseMatrix")
}

pieces <- list()
for (s in sorts) {
  k <- which(cell_sort == s)
  toc <- m[, k, drop = FALSE]

  prof <- fread(file.path(opt$soup_dir, s, "soup_profile.tsv.gz"))
  soup_df <- data.frame(est = prof$est, counts = prof$counts, row.names = prof$gene)
  stopifnot(identical(rownames(soup_df), rownames(toc)))

  step(s, ": ", length(k), " cells (", sum(clustered[k]), " clustered), rho = ", signif(rho[[s]], 4))
  set.seed(seed)
  # Clustered cells in one call at cluster level; unclustered cells at cell
  # level in chunks. adjustCounts' subtraction pulls each column out of a
  # dgTMatrix (O(nnz) per column), so one call is quadratic in its columns:
  # positive's ~7.6k unclustered cells as singleton clusters ran for hours.
  kc <- colnames(toc)[clustered[k]]
  ku <- colnames(toc)[!clustered[k]]
  parts <- list()
  if (length(kc)) {
    parts[["clustered"]] <- soupx_run(toc[, kc, drop = FALSE], soup_df, rho[[s]],
                                      clusters = setNames(cluster[k][clustered[k]], kc))
  }
  for (ch in split(ku, ceiling(seq_along(ku) / ADJ_CHUNK))) {
    parts[[length(parts) + 1]] <- soupx_run(toc[, ch, drop = FALSE], soup_df, rho[[s]],
                                            clusters = FALSE)
  }
  if (length(ku)) step(s, ": ", length(ku), " unclustered cells corrected at cell level")
  adj <- do.call(cbind, unname(parts))[, colnames(toc), drop = FALSE]
  stopifnot(identical(dim(adj), dim(toc)), identical(colnames(adj), colnames(toc)))
  pieces[[s]] <- adj
  rm(toc, adj, parts); invisible(gc())
}

out <- do.call(cbind, pieces)[, colnames(m), drop = FALSE]
stopifnot(identical(rownames(out), rownames(m)), identical(colnames(out), colnames(m)))
rm(pieces); invisible(gc())

cells <- data.table(bc_wells = colnames(m), sort = cell_sort, rho = unname(rho[cell_sort]),
                    cluster_level = clustered,
                    n_umis_raw = colSums(m), n_umis_corrected = colSums(out))
cells[, frac_removed := 1 - n_umis_corrected / n_umis_raw]
print(cells[, .(cells = .N, rho = rho[1],
                frac_removed_total = 1 - sum(n_umis_corrected) / sum(n_umis_raw),
                frac_removed_median = median(frac_removed)), by = .(sort, cluster_level)])

saveRDS(out, opt$out_counts, compress = FALSE)
fwrite(cells, opt$out_cells, sep = "\t")
step("wrote ", opt$out_counts, " and ", opt$out_cells)
log_peak()
