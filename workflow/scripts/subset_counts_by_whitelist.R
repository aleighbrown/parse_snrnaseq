#!/usr/bin/env Rscript
# Subset the per-sublibrary count matrices to a whitelist of cells.
#
#   Rscript subset_counts_by_whitelist.R <whitelist.csv> <prep/counts dir> <out.rds> [sublibraries.tsv]
#
# prep/counts/*.rds are bare dgCMatrix, genes x cells, with globally unique
# column names of the form `01_02_02__s16`. One file is held in memory at a
# time, so peak usage is about one sublibrary plus the accumulating subset --
# nowhere near the 15 GB combined matrix, and no big-memory job needed.
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

args     <- commandArgs(trailingOnly = TRUE)
wl_path  <- args[1]
counts_d <- args[2]
out_rds  <- args[3]
sublib_f <- if (length(args) >= 4) args[4] else NA_character_

wl <- fread(wl_path)
stopifnot("bc_wells" %in% names(wl))
cat(nrow(wl), "cells on the whitelist\n")

files <- list.files(counts_d, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)

# Orientation is decided by which axis carries barcodes, NOT by which is
# longer: 62,696 genes vs 24,652 cells means the gene axis is the long one.
is_barcode <- function(x) mean(grepl("__s", x, fixed = TRUE)) > 0.9

as_genes_by_cells <- function(x) {
  if (inherits(x, "Seurat")) x <- Seurat::GetAssayData(x, layer = "counts")
  if (is.list(x) && !isS4(x)) {
    x <- x[[which(vapply(x, function(e) inherits(e, "Matrix"), logical(1)))[1]]]
  }
  stopifnot(!is.null(rownames(x)), !is.null(colnames(x)))
  if (is_barcode(rownames(x)) && !is_barcode(colnames(x))) t(x) else x
}

# The combined run names sublibraries by number (`__s16`); other tables use the
# name (`__sParseFANS16`). Only convert if the whitelist is in the other one.
probe <- colnames(as_genes_by_cells(readRDS(files[1])))
if (!any(wl$bc_wells %in% probe)) {
  cat("no direct match on bc_wells:\n  counts:    ", probe[1],
      "\n  whitelist: ", wl$bc_wells[1], "\n")
  if (!is.na(sublib_f) && file.exists(sublib_f)) {
    subs <- fread(sublib_f)$sublibrary
    map  <- setNames(seq_along(subs), subs)          # name -> number
    wl[, bc_wells := paste0(sub("__s.*", "", bc_wells),
                            "__s", map[sub(".*__s", "", bc_wells)])]
    cat("converted whitelist to number form, e.g.", wl$bc_wells[1], "\n")
    stopifnot(any(wl$bc_wells %in% probe))
  } else {
    stop("pass config/sublibraries.tsv as the 4th argument to convert")
  }
}
rm(probe); gc()

want   <- unique(wl$bc_wells)
pieces <- list()
for (f in files) {
  m   <- as_genes_by_cells(readRDS(f))
  hit <- intersect(colnames(m), want)
  cat(sprintf("%-14s %6d cells in file, %5d kept\n",
              sub("\\.rds$", "", basename(f)), ncol(m), length(hit)))
  if (length(hit)) pieces[[basename(f)]] <- m[, hit, drop = FALSE]
  rm(m); gc()
}
stopifnot(length(pieces) > 0)

# Same reference throughout, but make.unique() on the 358 duplicated gene names
# could in principle have run per file -- check rather than assume.
genes <- rownames(pieces[[1]])
ok    <- vapply(pieces, function(p) identical(rownames(p), genes), logical(1))
if (!all(ok)) {
  cat("gene names differ in:", paste(names(ok)[!ok], collapse = ", "), "\n")
  genes  <- Reduce(intersect, lapply(pieces, rownames))
  cat("falling back to", length(genes), "shared genes\n")
  pieces <- lapply(pieces, function(p) p[genes, , drop = FALSE])
}

counts <- do.call(cbind, pieces)
cat("subset:", nrow(counts), "genes x", ncol(counts), "cells,",
    length(counts@x), "non-zero entries\n")

missed <- setdiff(want, colnames(counts))
if (length(missed)) {
  cat(length(missed), "whitelisted cells not found, e.g.:",
      paste(head(missed, 3), collapse = ", "), "\n")
}

# Metadata in the column order of the matrix, so the two line up without a
# later join.
meta <- as.data.frame(wl[match(colnames(counts), wl$bc_wells), ])
rownames(meta) <- colnames(counts)
stopifnot(identical(rownames(meta), colnames(counts)))

saveRDS(list(counts = counts, meta = meta), out_rds, compress = "gzip")
cat("wrote", out_rds, "--", round(file.size(out_rds) / 1e6), "MB\n")
