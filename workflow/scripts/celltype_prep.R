#!/usr/bin/env Rscript
# =============================================================================
# Cell typing, stage 1 of 3: load + QC.
#
#   count_matrix.mtx (cells x genes) + cell_metadata.csv + donor table
#     -> donor join (s<N> -> sublibrary name, checked on barcode overlap)
#     -> per-cell QC metrics and filter flags for EVERY split-pipe cell
#     -> <out_dir>/cell_meta.rds             all cells, QC columns
#        <out_dir>/counts/<sublibrary>.rds   raw counts (genes x cells) of the
#                                            cells that pass the low-quality
#                                            filters, one file per sublibrary
#
# Next: celltype_doublets.R (one job per sublibrary), then celltype_cluster.R.
# See workflow/rules/celltype.smk.
# =============================================================================

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  dirname(normalizePath(f))
}
source(file.path(script_dir(), "celltype_common.R"))
suppressPackageStartupMessages(library(Matrix))

need_args(c("dge_dir", "donor_table", "sublibraries", "out_dir", "threads", "seed",
            "n_cells", "only_assigned",
            "min_genes", "max_genes", "min_counts", "max_counts", "max_pct_mt"))
threads <- as.integer(num("threads"))
seed    <- as.integer(num("seed"))
setup_threads(threads)
dir.create(file.path(opt$out_dir, "counts"), recursive = TRUE, showWarnings = FALSE)
step("threads=", threads, "  out_dir=", opt$out_dir)

# =============================================================================
# 1. Cell metadata + donor join  (metadata only, no matrix yet)
# =============================================================================
cells <- fread(file.path(opt$dge_dir, "cell_metadata.csv"))
genes <- fread(file.path(opt$dge_dir, "all_genes.csv"))

cells[, mtx_row := .I]                                   # row in count_matrix.mtx
cells[, sub_n   := as.integer(sub(".*__s", "", bc_wells))]
cells[, barcode := sub("__s.*", "", bc_wells)]
stopifnot(!anyNA(cells$sub_n), !anyDuplicated(cells$bc_wells))

sublibs <- read.delim(opt$sublibraries, comment.char = "#")$sublibrary
step(nrow(cells), " cells, ", nrow(genes), " genes, ", length(sublibs), " sublibraries")

dx <- fread(opt$donor_table)
# The collaborator's raw table works too: barcode = "<bc_wells>__s<sublibrary>".
if (!"donor_id" %in% names(dx) && "donor" %in% names(dx)) {
  setnames(dx, "donor", "donor_id")
  dx[, sublibrary := sub(".*__s", "", barcode)]
  dx[, barcode    := sub("__s.*", "", barcode)]
}
stopifnot(all(c("barcode", "sublibrary", "donor_id") %in% names(dx)))
stopifnot(!anyDuplicated(dx, by = c("barcode", "sublibrary")))
# Keep the collaborator's QC / labels apart from the ones made here.
collab <- intersect(c("cell_type", "cluster_label", "n_genes", "n_counts", "percent_mito"), names(dx))
setnames(dx, collab, paste0("collab_", collab))

# bc_wells ends in the sublibrary NUMBER (s1) while the donor table uses the
# NAME (ParseFANS1). s<N> should be row N of sublibraries.tsv -- check it against
# barcode overlap rather than assume it (barcodes repeat across sublibraries, so
# the true pairing has to stand far above the chance matches).
hits <- merge(cells[, .(barcode, sub_n)], dx[, .(barcode, sublibrary)],
              by = "barcode", allow.cartesian = TRUE)[, .N, by = .(sub_n, sublibrary)]
best <- hits[order(sub_n, -N)][, .(best = sublibrary[1], n_best = N[1], n_runner_up = N[2]), by = sub_n]
best[, expected := sublibs[sub_n]]
print(best)

sub_map <- setNames(sublibs, seq_along(sublibs))
if (!all(best$best == best$expected)) {
  stopifnot(!anyDuplicated(best$best))
  warning("s<N> does NOT follow sublibraries.tsv order -- using the barcode-overlap mapping")
  sub_map[as.character(best$sub_n)] <- best$best
}
cells[, sublibrary := sub_map[as.character(sub_n)]]

drop <- c("cell_id", grep("^sub_lib", names(dx), value = TRUE))
cols <- setdiff(names(dx), drop)
n_before <- nrow(cells)
cells <- merge(cells, dx[, ..cols], by = c("barcode", "sublibrary"), all.x = TRUE, sort = FALSE)
setorder(cells, mtx_row)
stopifnot(nrow(cells) == n_before, identical(cells$mtx_row, seq_len(n_before)))

cells[, assigned := !is.na(donor_id)]
cells[is.na(donor_id), donor_id := "unassigned"]
if (!any(cells$assigned)) stop("No donor-table barcodes matched -- check the barcode / sublibrary formats.")
step(sprintf("%d of %d donor-table cells found; %d of %d split-pipe cells have a donor (%.1f%%)",
             sum(cells$assigned), nrow(dx), sum(cells$assigned), nrow(cells), 100 * mean(cells$assigned)))
rm(dx, hits); invisible(gc())

# =============================================================================
# 2. Count matrix
# =============================================================================
only_assigned <- flag("only_assigned")
n_cells <- as.integer(num("n_cells"))

# n_cells > 0 is the smoke-test path: read a random subset only.
use_rows <- NULL
if (n_cells > 0) {
  pool <- if (only_assigned) cells[assigned == TRUE, mtx_row] else cells$mtx_row
  if (n_cells < length(pool)) {
    set.seed(seed)
    use_rows <- sort(sample(pool, n_cells))
  }
}

# count_matrix.mtx is CELLS x GENES, rows in cell_metadata.csv order. Returns
# genes x cells. 1.1e9 entries: read as three int columns (~13 GB), then --
# since split-pipe writes it sorted by cell -- the column pointers are just
# cumulative counts, and the dgCMatrix is built without sparseMatrix()'s
# triplet sort (which would roughly double the peak).
read_parse_mtx <- function(path, n_cells, n_genes, rows = NULL) {
  hdr <- readLines(path, n = 50)
  n_skip <- sum(startsWith(hdr, "%")) + 1L                 # comments + dims line
  dims <- as.numeric(strsplit(trimws(hdr[n_skip]), "\\s+")[[1]])
  if (dims[1] != n_cells || dims[2] != n_genes) {
    stop(sprintf("matrix is %s x %s; expected cells x genes = %d x %d",
                 dims[1], dims[2], n_cells, n_genes))
  }
  if (dims[3] > .Machine$integer.max) {
    stop("more than 2^31-1 non-zeros: too many for one dgCMatrix; split the run")
  }

  if (is.null(rows)) {
    trip <- fread(path, skip = n_skip, sep = " ", header = FALSE,
                  col.names = c("cell", "gene", "x"), colClasses = "integer",
                  showProgress = FALSE)
    if (nrow(trip) != dims[3]) stop("read ", nrow(trip), " entries; header says ", dims[3])
    rows <- seq_len(n_cells)
  } else {
    keep_file <- tempfile(fileext = ".txt")
    fwrite(list(rows), keep_file, col.names = FALSE)
    awk <- if (nzchar(Sys.which("mawk"))) "mawk" else "awk"
    cmd <- sprintf("%s 'NR == FNR { k[$1]; next } FNR > %d && ($1 in k)' %s %s",
                   awk, n_skip, shQuote(keep_file), shQuote(path))
    trip <- fread(cmd = cmd, sep = " ", header = FALSE,
                  col.names = c("cell", "gene", "x"), colClasses = "integer")
    unlink(keep_file)
    new_col <- integer(n_cells)
    new_col[rows] <- seq_along(rows)
    trip[, cell := new_col[cell]]
  }
  step("read ", nrow(trip), " matrix entries")

  n_out <- length(rows)
  if (!is.unsorted(trip$cell)) {
    p <- c(0L, cumsum(tabulate(trip$cell, nbins = n_out)))
    i <- trip$gene - 1L
    x <- as.numeric(trip$x)
    rm(trip); invisible(gc())
    # new() validates, which includes gene order within each cell.
    m <- tryCatch(new("dgCMatrix", i = i, p = p, x = x, Dim = c(n_genes, n_out)),
                  error = function(e) NULL)
    if (is.null(m)) {
      step("gene order within cells is not sorted; falling back to sparseMatrix()")
      m <- sparseMatrix(i = i + 1L, p = p, x = x, dims = c(n_genes, n_out), repr = "C")
    }
  } else {
    step("matrix is not sorted by cell; using sparseMatrix()")
    m <- sparseMatrix(i = trip$gene, j = trip$cell, x = as.numeric(trip$x),
                      dims = c(n_genes, n_out), repr = "C")
  }
  m
}

counts <- read_parse_mtx(file.path(opt$dge_dir, "count_matrix.mtx"),
                         n_cells = nrow(cells), n_genes = nrow(genes), rows = use_rows)
meta <- if (is.null(use_rows)) cells else cells[use_rows]
rm(cells); invisible(gc())
rownames(counts) <- make.unique(genes$gene_name)           # 358 duplicated names
colnames(counts) <- meta$bc_wells
step("count matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")

# =============================================================================
# 3. Per-cell QC
# =============================================================================
mt_genes   <- grepl("^MT-", rownames(counts))
ribo_genes <- grepl("^RP[SL]", rownames(counts))
if (!any(mt_genes)) warning("No ^MT- genes found: percent_mt will be 0 for every cell")

meta[, n_transcripts := colSums(counts)]
meta[, n_genes       := diff(counts@p)]
meta[, pct_mt        := 100 * colSums(counts[mt_genes, , drop = FALSE])   / pmax(n_transcripts, 1)]
meta[, pct_ribo      := 100 * colSums(counts[ribo_genes, , drop = FALSE]) / pmax(n_transcripts, 1)]

# Guards against a row-order or orientation mistake: our totals must track
# split-pipe's own per-cell counts.
tscp_cor <- cor(meta$n_transcripts, meta$tscp_count)
step("cor(our UMI totals, split-pipe tscp_count) = ", round(tscp_cor, 4))
stopifnot(tscp_cor > 0.99)

max_counts <- num("max_counts")
meta[, `:=`(
  fail_low_genes        = n_genes < num("min_genes"),
  fail_high_genes       = n_genes > num("max_genes"),
  fail_low_transcripts  = n_transcripts < num("min_counts"),
  fail_high_transcripts = max_counts > 0 & n_transcripts > max_counts,
  fail_high_mt          = pct_mt > num("max_pct_mt")
)]
# scDblFinder wants low-quality cells gone but NOT high-count ones -- those are
# exactly the doublet-like cells it learns from. So the upper bounds are only
# applied in the cluster stage.
meta[, pre_doublet_pass := !(fail_low_genes | fail_low_transcripts | fail_high_mt)]
step(sum(meta$pre_doublet_pass), " of ", nrow(meta), " cells pass the low-quality filters")

# =============================================================================
# 4. Write
# =============================================================================
# One file per sublibrary in sublibraries.tsv -- even an empty one, so the
# outputs Snakemake expects always exist. Uncompressed: these are scratch
# intermediates (~13 GB in total for the full run) and gzip would add many
# minutes on each read and write.
for (s in sublibs) {
  k <- which(meta$pre_doublet_pass & meta$sublibrary == s)
  saveRDS(counts[, k, drop = FALSE], counts_file(opt$out_dir, s), compress = FALSE)
  step("wrote ", s, ": ", length(k), " cells")
}
saveRDS(meta, file.path(opt$out_dir, "cell_meta.rds"))
step("wrote cell_meta.rds: ", nrow(meta), " cells")

log_peak()
