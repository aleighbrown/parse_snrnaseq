#!/usr/bin/env Rscript
# =============================================================================
# Compare the two doublet arms, cell by cell:
#   default   celltype_doublets.R   one fixed dbr per sublibrary
#   by_sort   doublets_by_sort.R    one scDblFinder call per sort per
#                                   sublibrary, dbr from barcode occupancy
#
#   <prep>/cell_meta.rds + <default_dir>/*.tsv + <by_sort_dir>/*.tsv
#   + cell_annotations.csv (labels, UMAP; from the default arm's clustering)
#     -> cell_doublet_calls.tsv.gz        every tested cell, both calls, and
#                                         qc_pass under each arm (a whitelist
#                                         source for downstream work)
#        rates_by_sublibrary_sort.tsv     expected vs called rate, both arms
#        summary_by_sort.tsv              overlap of the two doublet sets
#        summary_by_sort_celltype.tsv     where the calls move, by cell type
#        doublets_compare.pdf
#
# Cell-type labels exist only for cells the default arm kept, so a cell the
# default arm removed as a doublet and by_sort keeps is shown as unlabelled.
# qc_pass_* assume celltype.doublets.remove: true, as in the default arm.
# =============================================================================

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  dirname(normalizePath(f))
}
source(file.path(script_dir(), "celltype_common.R"))
suppressPackageStartupMessages(library(ggplot2))

need_args(c("prep_dir", "default_dir", "by_sort_dir", "annotations", "sublibraries",
            "default_dbr", "out_dir", "threads"))
setDTthreads(as.integer(num("threads")))
out <- function(f) file.path(opt$out_dir, f)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
sublibs <- read.delim(opt$sublibraries, comment.char = "#")$sublibrary

# =============================================================================
# 1. Join
# =============================================================================
meta <- readRDS(file.path(opt$prep_dir, "cell_meta.rds"))
fail_cols <- grep("^fail_", names(meta), value = TRUE)          # prep's QC flags only
meta[, qc_no_doublet := !Reduce(`|`, .SD), .SDcols = fail_cols]
meta <- meta[pre_doublet_pass == TRUE,
             .(bc_wells, sublibrary, sample, donor_id, assigned, n_transcripts, n_genes, qc_no_doublet)]
step(nrow(meta), " cells tested for doublets")

read_arm <- function(dir, pattern) rbindlist(lapply(sublibs, function(s) {
  fread(file.path(dir, sprintf(pattern, s)), colClasses = c(doublet_class = "character"))
}), fill = TRUE)
dflt <- read_arm(opt$default_dir, "%s.tsv")
bys  <- read_arm(opt$by_sort_dir, "%s.tsv")
rates <- read_arm(opt$by_sort_dir, "%s.rates.tsv")
setnames(dflt, c("doublet_score", "doublet_class"), c("score_default", "class_default"))
setnames(bys, c("doublet_score", "doublet_class"), c("score_by_sort", "class_by_sort"))

d <- merge(meta, dflt, by = "bc_wells", all.x = TRUE)
d <- merge(d, bys, by = "bc_wells", all.x = TRUE)
stopifnot(nrow(d) == nrow(meta), all(d$sort == d$sample, na.rm = TRUE))
d[, sort := NULL]

ann_hdr  <- names(fread(opt$annotations, nrows = 0))
ann_cols <- intersect(c("bc_wells", "predicted.class", "predicted.subclass",
                        "umap_harmony_1", "umap_harmony_2"), ann_hdr)
ann <- fread(opt$annotations, select = ann_cols)
d <- merge(d, ann, by = "bc_wells", all.x = TRUE)
label_col <- intersect(c("predicted.subclass", "predicted.class"), names(d))[1]
d[, label := if (is.na(label_col)) NA_character_ else as.character(get(label_col))]
rm(meta, dflt, bys, ann); invisible(gc())

# Both arms called = comparable. A sort left uncalled by by_sort (too few
# cells, test run only) is kept in the table but out of the comparisons.
d[, `:=`(dbl_default = class_default %in% "doublet", dbl_by_sort = class_by_sort %in% "doublet")]
d[, compared := !is.na(class_default) & !is.na(class_by_sort)]
d[, call := fifelse(!compared, "not compared",
            fifelse(dbl_default & dbl_by_sort, "doublet in both",
            fifelse(dbl_default, "default only",
            fifelse(dbl_by_sort, "by_sort only", "singlet in both"))))]
d[, `:=`(qc_pass_default = qc_no_doublet & !dbl_default,
         qc_pass_by_sort = qc_no_doublet & !dbl_by_sort)]
d[, label := fifelse(is.na(label) & call %in% c("default only", "doublet in both"),
                     "(removed by default arm)", fifelse(is.na(label), "(unlabelled)", label))]
step(sum(d$compared), " cells called by both arms")

keep <- c("bc_wells", "sublibrary", "sample", "donor_id", "assigned", "n_transcripts", "n_genes",
          "label", "score_default", "class_default", "score_by_sort", "class_by_sort", "dbr_used",
          "call", "qc_pass_default", "qc_pass_by_sort")
fwrite(d[, ..keep], out("cell_doublet_calls.tsv.gz"), sep = "\t")

# =============================================================================
# 2. Summaries
# =============================================================================
default_dbr <- num("default_dbr")
rate_tab <- merge(
  rates[, .(sublibrary, sort, n_called_barcodes, barcode_space, lambda,
            expected_by_sort = dbr, n_tested, called_by_sort = called_rate)],
  d[compared == TRUE, .(called_default = mean(dbl_default)), by = .(sublibrary, sort = sample)],
  by = c("sublibrary", "sort"), all.x = TRUE)
rate_tab[, expected_default := default_dbr]
fwrite(rate_tab, out("rates_by_sublibrary_sort.tsv"), sep = "\t")

overlap <- function(x) x[, .(
  n_compared      = .N,
  doublets_default = sum(dbl_default),
  doublets_by_sort = sum(dbl_by_sort),
  in_both         = sum(dbl_default & dbl_by_sort),
  default_only    = sum(dbl_default & !dbl_by_sort),
  by_sort_only    = sum(!dbl_default & dbl_by_sort),
  rate_default    = mean(dbl_default),
  rate_by_sort    = mean(dbl_by_sort),
  jaccard         = sum(dbl_default & dbl_by_sort) / max(1, sum(dbl_default | dbl_by_sort)),
  score_spearman  = suppressWarnings(cor(score_default, score_by_sort, method = "spearman"))
)]
by_sort <- rbind(d[compared == TRUE, overlap(.SD), by = .(sort = sample)],
                 overlap(d[compared == TRUE])[, sort := "TOTAL"], use.names = TRUE)
fwrite(by_sort, out("summary_by_sort.tsv"), sep = "\t")
print(by_sort)

by_type <- d[compared == TRUE, .(
  n_compared   = .N,
  default_only = sum(dbl_default & !dbl_by_sort),
  by_sort_only = sum(!dbl_default & dbl_by_sort),
  in_both      = sum(dbl_default & dbl_by_sort),
  frac_by_sort_only = mean(!dbl_default & dbl_by_sort),
  median_score_default = median(score_default),
  median_score_by_sort = median(score_by_sort)
), by = .(sort = sample, label)][order(sort, -n_compared)]
fwrite(by_type, out("summary_by_sort_celltype.tsv"), sep = "\t")

# =============================================================================
# 3. Report
# =============================================================================
call_cols <- c("singlet in both" = "grey75", "doublet in both" = "black",
               "default only" = "#1b9e77", "by_sort only" = "#d95f02")
dc <- d[compared == TRUE]
dc[, call := factor(call, levels = names(call_cols))]

pdf(out("doublets_compare.pdf"), width = 11, height = 7)

rl <- melt(rate_tab, id.vars = c("sublibrary", "sort"),
           measure.vars = c("expected_by_sort", "called_by_sort", "called_default"),
           variable.name = "measure", value.name = "rate")
print(ggplot(rl[!is.na(rate)], aes(sort, pmax(rate, 1e-4), colour = measure)) +
        geom_hline(yintercept = default_dbr, linetype = 2, colour = "grey50") +
        geom_point(position = position_dodge(width = 0.6), alpha = 0.7) +
        scale_y_log10(labels = scales::percent) +
        labs(x = NULL, y = "doublet rate (log scale; 0 shown at 0.01%)",
             title = "Doublet rate per sublibrary x sort",
             subtitle = sprintf("dashed = default arm's dbr (%.2f); one point per sublibrary", default_dbr)) +
        theme_bw())

cnt <- dc[call != "singlet in both", .N, by = .(sort = sample, call)]
print(ggplot(cnt, aes(sort, N, fill = call)) + geom_col() +
        scale_fill_manual(values = call_cols) +
        labs(x = NULL, y = "cells called doublet", title = "Doublet calls: overlap of the two arms") +
        theme_bw())

set.seed(1)
pts <- if (nrow(dc) > 60000) dc[sample(.N, 60000)] else dc
pts <- rbind(pts, dc[call != "singlet in both" & !bc_wells %in% pts$bc_wells])   # keep every doublet
print(ggplot(pts[order(call)], aes(score_default, score_by_sort, colour = call)) +
        geom_point(size = 0.3, alpha = 0.5) +
        scale_colour_manual(values = call_cols) +
        facet_wrap(~ sample) +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        labs(x = "scDblFinder score, default arm", y = "scDblFinder score, by_sort arm",
             title = "Scores, cell by cell (singlets subsampled)") +
        theme_bw())

print(ggplot(dc, aes(call, n_transcripts, fill = call)) +
        geom_boxplot(outlier.size = 0.2) +
        scale_y_log10() + scale_fill_manual(values = call_cols, guide = "none") +
        facet_wrap(~ sample) +
        labs(x = NULL, y = "UMIs per cell", title = "Depth by call",
             subtitle = "real doublets should sit above the singlets") +
        theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1)))

bt <- by_type[!label %in% c("(unlabelled)", "(removed by default arm)") & n_compared >= 20]
if (nrow(bt)) {
  print(ggplot(bt, aes(reorder(label, frac_by_sort_only), frac_by_sort_only)) +
          geom_col(fill = call_cols[["by_sort only"]]) +
          coord_flip() + facet_wrap(~ sort, scales = "free_x") +
          scale_y_continuous(labels = scales::percent) +
          labs(x = NULL, y = "share of the default arm's kept cells newly called doublet",
               title = "Where the by_sort arm removes extra cells",
               subtitle = paste0("label: ", label_col, "; types with >= 20 cells")) +
          theme_bw())
}

if (all(c("umap_harmony_1", "umap_harmony_2") %in% names(dc))) {
  u <- dc[!is.na(umap_harmony_1)]
  print(ggplot(u, aes(umap_harmony_1, umap_harmony_2)) +
          geom_point(data = u[call == "singlet in both"], colour = "grey85", size = 0.1) +
          geom_point(data = u[call == "by_sort only"], colour = call_cols[["by_sort only"]], size = 0.3) +
          facet_wrap(~ sample) + coord_equal() +
          labs(title = "Cells the default arm kept and by_sort calls doublet",
               subtitle = "Harmony UMAP from celltype_cluster") +
          theme_void() + theme(strip.text = element_text(size = 11)))
}

invisible(dev.off())
step("wrote ", out("doublets_compare.pdf"))
log_peak()
