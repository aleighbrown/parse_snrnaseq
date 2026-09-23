#!/usr/bin/env Rscript
# =============================================================================
# Ambient RNA, step 1: the soup profile of ONE FANS sort.
#
#   splitpipe/combined/<sort>/DGE_unfiltered/  +  all-sample called cells
#     -> empty barcodes: not a called cell, UMIs within [umi_min, umi_max]
#     -> their counts summed per gene (and per sublibrary), streamed with awk
#     -> <out_dir>/soup_profile.tsv.gz         pooled: gene, counts, est
#        <out_dir>/soup_by_sublibrary.tsv.gz   the same split by sublibrary
#        <out_dir>/barcode_summary.tsv         barcodes / empties / soup UMIs
#        <out_dir>/top_soup_genes.tsv          for choosing rho marker sets
#        <out_dir>/soup_report.pdf             barcode ranks, consistency
#
# Why per sort. All 36 donors were pooled within each sort before loading, so
# each sort's nuclei sat in one suspension with one pool of free RNA. Ambient
# molecules are reverse-transcribed in the round-1 wells, so they carry the
# sort's round-1 barcode -- which is how split-pipe already separated them
# into combined/<sort>/. Sublibraries are split AFTER pooling, so every
# sublibrary holds a share of all three soups; the per-sublibrary split is
# kept only to check that they agree.
#
# Cells are excluded using the ALL-SAMPLE called cells, because those are the
# cells the rest of the pipeline uses.
#
# Next (not yet written): soupx_rho.R, then soupx_adjust.R.
# See workflow/rules/soupx.smk.
# =============================================================================

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  dirname(normalizePath(f))
}
source(file.path(script_dir(), "celltype_common.R"))
suppressPackageStartupMessages({
  library(ggplot2)
})

need_args(c("dge_dir", "called_cells", "ref_genes", "sublibraries", "sort", "out_dir",
            "threads", "umi_min", "umi_max", "n_top"))
setup_threads(as.integer(num("threads")))
umi_min <- num("umi_min")
umi_max <- num("umi_max")
stopifnot(umi_min >= 1, umi_max > umi_min)
out <- function(f) file.path(opt$out_dir, f)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
step("sort=", opt$sort, "  empty range=[", umi_min, ", ", umi_max, "]  out_dir=", opt$out_dir)

# =============================================================================
# 1. Genes: must line up with the all-sample matrix the counts come from
# =============================================================================
genes <- fread(file.path(opt$dge_dir, "all_genes.csv"))
ref   <- fread(opt$ref_genes)
if (!identical(genes$gene_id, ref$gene_id)) {
  stop("all_genes.csv of ", opt$sort, " differs from all-sample's: the soup would ",
       "not line up with the count matrices")
}
genes[, gene := make.unique(gene_name)]          # same rownames as celltype_prep.R
rm(ref)

# =============================================================================
# 2. Barcodes: called cells, empties, everything else
# =============================================================================
cells <- fread(file.path(opt$dge_dir, "cell_metadata.csv"), select = c("bc_wells", "tscp_count"))
cells[, mtx_row := .I]                                   # row in count_matrix.mtx
cells[, sub_n   := as.integer(sub(".*__s", "", bc_wells))]
stopifnot(!anyNA(cells$sub_n), !anyDuplicated(cells$bc_wells))

# s<N> -> sublibrary name by row of sublibraries.tsv. celltype_prep.R checks
# that assumption against barcode overlap; its log says if it ever fails.
sublibs <- read.delim(opt$sublibraries, comment.char = "#")$sublibrary
cells[, sublibrary := sublibs[sub_n]]
stopifnot(!anyNA(cells$sublibrary))

called <- fread(opt$called_cells, select = "bc_wells")$bc_wells
cells[, status := fifelse(bc_wells %chin% called, "cell",
                  fifelse(tscp_count >= umi_min & tscp_count <= umi_max, "empty", "other"))]
rm(called)
if (!any(cells$status == "cell")) {
  stop("none of ", opt$sort, "'s barcodes are all-sample called cells -- bc_wells formats differ?")
}
if (!any(cells$status == "empty")) stop("no barcodes in the empty range")

n_cell_in_range <- cells[status == "cell", sum(tscp_count <= umi_max)]
step(nrow(cells), " barcodes: ", sum(cells$status == "cell"), " called cells (min ",
     cells[status == "cell", min(tscp_count)], " UMIs), ", sum(cells$status == "empty"),
     " empties, ", sum(cells$status == "other"), " other")
if (n_cell_in_range > 0) {
  warning(n_cell_in_range, " called cells have <= ", umi_max, " UMIs. They are excluded from ",
          "the soup regardless, but the range reaches into cell territory: see the rank plot")
}

# The background floor: what a barcode with no nucleus picks up. Taken over
# ALL uncalled barcodes, not just the empty window, so it does not depend on
# where the window was put. Nearly every possible barcode is seen in these
# sublibraries, so the uncalled ones are a near-complete sample of it. The
# floor / cell ratio is the per-barcode contamination a median cell would
# carry if the background is a fixed count per barcode.
bc_summary <- cells[, .(n_barcodes = .N,
                        n_cells    = sum(status == "cell"),
                        n_empty    = sum(status == "empty"),
                        empty_umis = sum(as.numeric(tscp_count[status == "empty"])),
                        median_cell_umis     = as.numeric(median(tscp_count[status == "cell"])),
                        uncalled_q25_umis    = as.numeric(quantile(tscp_count[status != "cell"], 0.25)),
                        uncalled_median_umis = as.numeric(median(tscp_count[status != "cell"])),
                        uncalled_q75_umis    = as.numeric(quantile(tscp_count[status != "cell"], 0.75)),
                        uncalled_mean_umis   = mean(tscp_count[status != "cell"])),
                    by = .(sub_n, sublibrary)][order(sub_n)]
bc_summary[, floor_to_cell := uncalled_median_umis / median_cell_umis]
bc_summary[, sort := opt$sort]
fwrite(bc_summary, out("barcode_summary.tsv"), sep = "\t")
print(bc_summary)

# =============================================================================
# 3. Stream the matrix, summing the empties' counts per sublibrary x gene
# =============================================================================
# count_matrix.mtx is CELLS x GENES, rows in cell_metadata.csv order (as in
# celltype_prep.R). Nothing is held but the lookup of empty rows and the
# running sums (<= 16 x 63k keys), so memory stays small whatever the file.
mtx <- file.path(opt$dge_dir, "count_matrix.mtx")
hdr <- readLines(mtx, n = 50)
n_skip <- sum(startsWith(hdr, "%")) + 1L                 # comments + dims line
dims <- as.numeric(strsplit(trimws(hdr[n_skip]), "\\s+")[[1]])
if (dims[1] != nrow(cells) || dims[2] != nrow(genes)) {
  stop(sprintf("matrix is %s x %s; expected cells x genes = %d x %d",
               dims[1], dims[2], nrow(cells), nrow(genes)))
}

keep_file <- tempfile(fileext = ".txt")
fwrite(cells[status == "empty", .(mtx_row, sub_n)], keep_file, sep = " ", col.names = FALSE)
awk <- if (nzchar(Sys.which("mawk"))) "mawk" else "awk"
# printf, not print: awk would write a large sum as %.6g.
prog <- sprintf(paste0("NR == FNR { k[$1] = $2; next } ",
                       "FNR > %d && ($1 in k) { s[k[$1] \" \" $2] += $3 } ",
                       "END { for (x in s) printf \"%%s %%.0f\\n\", x, s[x] }"), n_skip)
cmd <- sprintf("%s '%s' %s %s", awk, prog, shQuote(keep_file), shQuote(mtx))
step("streaming ", round(file.size(mtx) / 1e9, 1), " GB with ", awk)
soup <- fread(cmd = cmd, sep = " ", header = FALSE,
              col.names = c("sub_n", "gene_idx", "count"),
              colClasses = c("integer", "integer", "numeric"))
unlink(keep_file)
step("read ", nrow(soup), " sublibrary x gene sums")

# Guard against a row-order mistake: the summed matrix counts must match
# split-pipe's own UMI totals for the same barcodes.
expected <- sum(as.numeric(cells[status == "empty", tscp_count]))
ratio <- sum(soup$count) / expected
step(sprintf("soup UMIs %.0f vs tscp_count of the empties %.0f (ratio %.3f)",
             sum(soup$count), expected, ratio))
if (ratio < 0.5 || ratio > 2) stop("matrix totals do not match tscp_count: wrong rows selected?")
if (abs(ratio - 1) > 0.05) warning("soup total differs from tscp_count by more than 5%")

# =============================================================================
# 4. Profiles
# =============================================================================
soup[, `:=`(gene = genes$gene[gene_idx], gene_id = genes$gene_id[gene_idx],
            sublibrary = sublibs[sub_n])]
soup[, est := count / sum(count), by = sub_n]
setorder(soup, sub_n, gene_idx)
fwrite(soup[, .(sort = opt$sort, sub_n, sublibrary, gene, gene_id, count, est)],
       out("soup_by_sublibrary.tsv.gz"), sep = "\t")

# Pooled over sublibraries, every gene (zeros included), in matrix row order.
# `counts` and `est` are the two columns SoupX's setSoupProfile() wants.
pooled <- data.table(gene = genes$gene, gene_id = genes$gene_id,
                     counts = 0)
tot <- soup[, .(counts = sum(count)), by = gene_idx]
pooled[tot$gene_idx, counts := tot$counts]
pooled[, est := counts / sum(counts)]
fwrite(pooled, out("soup_profile.tsv.gz"), sep = "\t")
step(sum(pooled$counts > 0), " of ", nrow(pooled), " genes seen in the soup")

pooled[, is_mt := grepl("^MT-", gene)]
step(sprintf("%.1f%% of soup UMIs are mitochondrial", 100 * sum(pooled$est[pooled$is_mt])))

top <- pooled[order(-counts)][seq_len(min(num("n_top"), .N))]
top[, `:=`(rank = .I, cum_est = cumsum(est))]
fwrite(top[, .(rank, gene, gene_id, counts, est, cum_est, is_mt)], out("top_soup_genes.tsv"), sep = "\t")
print(head(top[, .(rank, gene, counts, est = signif(est, 3))], 30))

# =============================================================================
# 5. Report
# =============================================================================
pdf(out("soup_report.pdf"), width = 9, height = 6)

# Barcode ranks, thinned to log-spaced ranks: millions of points add nothing.
rk <- cells[order(-tscp_count), .(tscp_count, status)]
rk[, rank := .I]
show <- unique(round(10^seq(0, log10(nrow(rk)), length.out = 4000)))
print(ggplot(rk[show], aes(rank, tscp_count, colour = status)) +
        annotate("rect", xmin = 1, xmax = Inf, ymin = umi_min, ymax = umi_max,
                 fill = "grey85", alpha = 0.5) +
        geom_point(size = 0.5) +
        scale_x_log10() + scale_y_log10() +
        scale_colour_manual(values = c(cell = "#1F4E79", empty = "#A50F15", other = "grey60")) +
        labs(x = "barcode rank", y = "UMIs (tscp_count)", colour = NULL,
             title = paste0(opt$sort, ": barcode ranks"),
             subtitle = sprintf("shaded = empty range [%g, %g]; it should sit well below the cells",
                                umi_min, umi_max)) +
        theme_minimal(base_size = 10))

# Floor per sublibrary against its cells. A flat ratio across sublibraries
# means the background scales with depth like the cells do.
print(ggplot(bc_summary, aes(factor(sublibrary, levels = sublibs), floor_to_cell)) +
        geom_col(fill = "grey40") +
        scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
        labs(x = NULL, y = "median uncalled UMIs / median cell UMIs",
             title = paste0(opt$sort, ": background floor relative to the cells"),
             subtitle = "the contamination a median cell carries if background is a fixed count per barcode") +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

print(ggplot(bc_summary, aes(factor(sublibrary, levels = sublibs), empty_umis)) +
        geom_col(fill = "#A50F15") +
        labs(x = NULL, y = "soup UMIs",
             title = paste0(opt$sort, ": soup UMIs per sublibrary")) +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Consistency: each sublibrary's soup against the pooled one, over genes that
# are well measured. Sublibraries share one soup, so these should all be high;
# a low one is a sublibrary with too few empty UMIs or a problem of its own.
well <- pooled[order(-counts)][seq_len(min(2000, sum(pooled$counts > 0))), gene]
cons <- soup[gene %chin% well,
             .(r = cor(log10(est), log10(pooled$est[match(gene, pooled$gene)]))),
             by = .(sublibrary)]
fwrite(cons, out("sublibrary_consistency.tsv"), sep = "\t")
print(ggplot(cons, aes(factor(sublibrary, levels = sublibs), r)) +
        geom_point(size = 2.5) +
        labs(x = NULL, y = "Pearson r, log10 soup fraction",
             title = paste0(opt$sort, ": each sublibrary's soup vs the pooled soup"),
             subtitle = "top 2000 soup genes; should be uniformly high") +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

print(ggplot(head(top, 40), aes(reorder(gene, est), est, fill = is_mt)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c(`FALSE` = "grey40", `TRUE` = "#D98C00"),
                          labels = c("nuclear", "mitochondrial"), name = NULL) +
        labs(x = NULL, y = "fraction of soup UMIs",
             title = paste0(opt$sort, ": top 40 soup genes")) +
        theme_minimal(base_size = 9))

invisible(dev.off())
step("report written")
log_peak()
