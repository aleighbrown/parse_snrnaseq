#!/usr/bin/env Rscript
# =============================================================================
# Ambient RNA, step 2: the contamination fraction (rho) of ONE FANS sort, by
# SoupX's own marker-gene method.
#
#   <run>/<labels> (a frozen cell_annotations.csv)   cluster + cell-type labels
#   <run>/prep/counts/*.rds                          raw counts
#   <run>/doublets_by_sort/*.tsv                     per-sort doublet calls
#   soup/<sort>/soup_profile.tsv.gz                  the measured soup
#     -> glia of this sort, confidently labelled, not doublets, not OPC
#     -> SoupChannel + setSoupProfile + setClusters
#     -> estimateNonExpressingCells (per cluster AND per cell)
#     -> calculateContaminationFraction   (Poisson GLM, counts ~ offset(expected))
#     -> <out_dir>/<sort>.tsv              rho + 95% CI, both exclusion modes
#        <out_dir>/<sort>_breakdown.tsv    the same estimate by subclass /
#                                          sublibrary / marker set
#        <out_dir>/<sort>_report.pdf
#
# Why glia only. The markers are neuronal genes, and the estimate uses only
# cells that should not express them; neurons would be excluded by SoupX
# anyway. OPCs genuinely express some of these synaptic genes (54-77% of them
# exceeded the soup expectation in glia_gene_exploration.R), so they are left
# out explicitly. Cells come from every sublibrary: rho belongs to the sort,
# and pooling gives the steadiest estimate.
#
# Exclusion. SoupX drops a whole cluster if any of its cells expresses a
# marker beyond what 100% soup could give (clusters = default), or only
# those cells (clusters = FALSE). Large Harmony clusters with a few missed
# doublets can lose every cluster, so both are fitted; the cluster mode is
# used when it keeps >= min_est_cells cells, else the cell mode.
#
# The breakdown uses the same estimator as SoupX (for an intercept-only
# Poisson GLM with an offset, the MLE is sum(observed) / sum(expected)), so
# its pooled value reproduces SoupX's rho exactly -- checked below.
#
# Next: soupx_adjust.R. See workflow/rules/soupx.smk.
# =============================================================================

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  dirname(normalizePath(f))
}
source(file.path(script_dir(), "celltype_common.R"))
suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(SoupX)
})

need_args(c("sort", "labels", "prep_dir", "doublet_dir", "sublibraries", "soup_profile",
            "barcode_summary", "marker_sets", "estimate_class", "exclude_subclasses",
            "min_label_score", "max_cells", "min_est_cells", "out_dir", "threads", "seed"))
setup_threads(as.integer(num("threads")))
seed <- as.integer(num("seed"))
this_sort <- opt$sort
out <- function(f) file.path(opt$out_dir, paste0(this_sort, f))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
step("SoupX ", as.character(packageVersion("SoupX")), "  sort=", this_sort)

# "neuronal=SYT1|SNAP25;other=A|B" -> list(neuronal = c(...), other = c(...))
ms <- strsplit(strsplit(opt$marker_sets, ";", fixed = TRUE)[[1]], "=", fixed = TRUE)
marker_sets <- setNames(lapply(ms, function(x) strsplit(x[2], "|", fixed = TRUE)[[1]]),
                        vapply(ms, `[`, "", 1))
print(marker_sets)

# =============================================================================
# 1. Which cells
# =============================================================================
lab <- fread(opt$labels,
             select = c("bc_wells", "sublibrary", "sample", "clustered", "seurat_clusters",
                        "predicted.class", "predicted.subclass", "predicted.subclass.score"),
             colClasses = c(seurat_clusters = "character"))
lab <- lab[sample == this_sort]
step(nrow(lab), " split-pipe cells in ", this_sort)

sublibs <- read.delim(opt$sublibraries, comment.char = "#")$sublibrary
dbl <- rbindlist(lapply(sublibs, function(s)
  fread(doublet_file(opt$doublet_dir, s), select = c("bc_wells", "doublet_class"),
        colClasses = c(doublet_class = "character"))))
lab[, by_sort_doublet := bc_wells %chin% dbl[doublet_class %in% "doublet", bc_wells]]
rm(dbl)

excl <- csv("exclude_subclasses")
lab[, why := fcase(!(clustered %in% TRUE),                               "not clustered",
                   by_sort_doublet,                                     "doublet (by sort)",
                   !predicted.class %in% csv("estimate_class"),         "not an estimate class",
                   predicted.subclass %in% excl,                        "excluded subclass",
                   is.na(predicted.subclass.score) |
                     predicted.subclass.score < num("min_label_score"),   "low label score",
                   default = "candidate")]
print(lab[, .N, by = why][order(-N)])

cand <- lab[why == "candidate"]
if (nrow(cand) > num("max_cells")) {
  set.seed(seed)
  keep <- sort(base::sample(nrow(cand), num("max_cells")))
  cand <- cand[keep]
}
step(nrow(cand), " candidate cells; by subclass:")
print(cand[, .N, by = predicted.subclass][order(-N)])
if (nrow(cand) < num("min_est_cells")) stop("too few candidate cells to estimate rho")

# =============================================================================
# 2. Their raw counts, one sublibrary at a time
# =============================================================================
want <- cand$bc_wells
pieces <- lapply(sublibs, function(s) {
  m <- readRDS(counts_file(opt$prep_dir, s))
  m[, colnames(m) %chin% want, drop = FALSE]
})
toc <- do.call(cbind, pieces)
rm(pieces); invisible(gc())
stopifnot(setequal(colnames(toc), want))
cand <- cand[match(colnames(toc), bc_wells)]
step("counts: ", nrow(toc), " genes x ", ncol(toc), " cells")

miss <- setdiff(unlist(marker_sets), rownames(toc))
if (length(miss)) stop("marker genes not in the count matrix: ", paste(miss, collapse = ", "))

# The soup profile was written in the same gene order and with the same
# make.unique() names as prep/counts (both from all-sample all_genes.csv).
prof <- fread(opt$soup_profile)
soup_df <- data.frame(est = prof$est, counts = prof$counts, row.names = prof$gene)
stopifnot(identical(rownames(soup_df), rownames(toc)))

# =============================================================================
# 3. SoupX
# =============================================================================
# tod is only used to estimate the soup, which is supplied instead
# (calcSoupProfile = FALSE + setSoupProfile), so toc stands in for it.
sc <- SoupChannel(toc, toc, calcSoupProfile = FALSE)
sc <- setSoupProfile(sc, soup_df)
sc <- setClusters(sc, setNames(cand$seurat_clusters, cand$bc_wells))

fit_mode <- function(clusters_arg) {
  use <- if (is.null(clusters_arg)) estimateNonExpressingCells(sc, marker_sets)
         else estimateNonExpressingCells(sc, marker_sets, clusters = clusters_arg)
  n_use <- sum(rowSums(use) > 0)
  if (n_use == 0) return(list(use = use, n = 0L, rho = NA_real_, low = NA_real_, high = NA_real_))
  fit <- calculateContaminationFraction(sc, marker_sets, useToEst = use, forceAccept = TRUE)
  list(use = use, n = n_use, rho = fit$metaData$rho[1],
       low = fit$metaData$rhoLow[1], high = fit$metaData$rhoHigh[1])
}
by_cluster <- fit_mode(NULL)
step(sprintf("cluster exclusion: %d cells kept, rho = %.4f (%.4f-%.4f)",
             by_cluster$n, by_cluster$rho, by_cluster$low, by_cluster$high))
by_cell <- fit_mode(FALSE)
step(sprintf("cell exclusion:    %d cells kept, rho = %.4f (%.4f-%.4f)",
             by_cell$n, by_cell$rho, by_cell$low, by_cell$high))

mode <- if (by_cluster$n >= num("min_est_cells")) "cluster" else "cell"
chosen <- if (mode == "cluster") by_cluster else by_cell
step("using ", mode, " exclusion: rho = ", signif(chosen$rho, 4))
if (chosen$rho > 0.3) warning("rho above 0.3 -- check the report before adjusting")

# =============================================================================
# 4. Breakdown with the same estimator, and the floor for comparison
# =============================================================================
n_umis <- sc$metaData$nUMIs
long <- rbindlist(lapply(names(marker_sets), function(nm) {
  g <- marker_sets[[nm]]
  data.table(bc_wells = colnames(toc), set = nm,
             obs = colSums(toc[g, , drop = FALSE]),
             exp = n_umis * sum(soup_df[g, "est"]),
             use = chosen$use[colnames(toc), nm])
}))
long <- merge(long, cand[, .(bc_wells, sublibrary, subclass = predicted.subclass)], by = "bc_wells")
used <- long[use == TRUE]

pooled <- used[, sum(obs) / sum(exp)]
step(sprintf("pooled sum(obs)/sum(exp) = %.4f (SoupX: %.4f)", pooled, chosen$rho))
if (abs(pooled - chosen$rho) > 1e-3 * max(chosen$rho, 1e-6)) {
  warning("breakdown estimator does not reproduce SoupX's rho: read the breakdown as relative only")
}

ci <- function(o, e) {
  r <- poisson.test(round(o), T = e)$conf.int
  list(low = r[1], high = r[2])
}
breakdown <- rbind(
  used[, c(.(level = "subclass", n_cells = uniqueN(bc_wells), obs = sum(obs), exp = sum(exp),
             rho = sum(obs) / sum(exp)), ci(sum(obs), sum(exp))), by = .(group = subclass)],
  used[, c(.(level = "sublibrary", n_cells = uniqueN(bc_wells), obs = sum(obs), exp = sum(exp),
             rho = sum(obs) / sum(exp)), ci(sum(obs), sum(exp))), by = .(group = sublibrary)],
  used[, c(.(level = "marker_set", n_cells = uniqueN(bc_wells), obs = sum(obs), exp = sum(exp),
             rho = sum(obs) / sum(exp)), ci(sum(obs), sum(exp))), by = .(group = set)])
setcolorder(breakdown, c("level", "group"))
breakdown[, sort := this_sort]
print(breakdown)
fwrite(breakdown, out("_breakdown.tsv"), sep = "\t")

# What the measured background floor alone would give these cells: the share
# of rho that a fixed count per barcode explains (glia_gene_exploration.R
# found ~70-75% in positive, ~30-50% in the negative sorts).
fl <- fread(opt$barcode_summary)
floor_rho <- mean(fl$uncalled_median_umis[match(cand$sublibrary, fl$sublibrary)] / n_umis)
step(sprintf("floor-only rho for the same cells (median floor / UMIs, mean) = %.4f", floor_rho))

res <- data.table(sort = this_sort, soupx_version = as.character(packageVersion("SoupX")),
                  n_candidates = nrow(cand), est_mode = mode,
                  rho = chosen$rho, rho_low = chosen$low, rho_high = chosen$high,
                  n_est = chosen$n,
                  rho_cluster = by_cluster$rho, rho_cluster_low = by_cluster$low,
                  rho_cluster_high = by_cluster$high, n_est_cluster = by_cluster$n,
                  rho_cell = by_cell$rho, rho_cell_low = by_cell$low,
                  rho_cell_high = by_cell$high, n_est_cell = by_cell$n,
                  floor_only_rho = floor_rho,
                  median_umis = median(n_umis),
                  marker_sets = opt$marker_sets)
print(t(res))
fwrite(res, out(".tsv"), sep = "\t")

# =============================================================================
# 5. Report
# =============================================================================
pdf(out("_report.pdf"), width = 9, height = 6)
for (lv in c("subclass", "sublibrary", "marker_set")) {
  print(ggplot(breakdown[level == lv], aes(reorder(group, rho), rho)) +
          geom_hline(yintercept = chosen$rho, colour = "#A50F15") +
          geom_hline(yintercept = floor_rho, colour = "grey50", linetype = "dashed") +
          geom_pointrange(aes(ymin = low, ymax = high)) +
          coord_flip() +
          scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
          labs(x = NULL, y = "rho = observed / expected soup counts (95% CI)",
               title = sprintf("%s: rho by %s", this_sort, lv),
               subtitle = sprintf("red = SoupX rho (%s exclusion, %.2f%%); dashed = floor only (%.2f%%)",
                                  mode, 100 * chosen$rho, 100 * floor_rho)) +
          theme_minimal(base_size = 10))
}

# Per cell: observed marker counts against the soup expectation at rho = 1.
# Excluded cells (real expressers, doublets) sit above the diagonal.
cell_tot <- long[, .(obs = sum(obs), exp = sum(exp), use = any(use)), by = .(bc_wells, subclass)]
print(ggplot(cell_tot, aes(exp, obs + 0.5, colour = use)) +
        geom_point(size = 0.3, alpha = 0.3) +
        geom_abline(slope = chosen$rho, intercept = 0, colour = "#A50F15") +
        geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
        scale_x_log10() + scale_y_log10() +
        scale_colour_manual(values = c(`TRUE` = "grey40", `FALSE` = "#D98C00"),
                            labels = c(`TRUE` = "used", `FALSE` = "excluded"), name = NULL) +
        facet_wrap(~subclass) +
        labs(x = "expected marker counts if the nucleus were all soup",
             y = "observed marker counts + 0.5",
             title = paste0(this_sort, ": marker counts per cell"),
             subtitle = "red line = rho; dotted = 100% soup") +
        theme_minimal(base_size = 9))
invisible(dev.off())

log_peak()
