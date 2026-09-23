#!/usr/bin/env Rscript
# =============================================================================
# Alternative doublet arm: scDblFinder on ONE sublibrary, once per FANS sort,
# each at the collision rate of that sort's own barcode space.
#
#   <prep>/counts/<sublibrary>.rds + split-pipe cell_metadata.csv
#     -> <out>        bc_wells, sort, doublet_score, doublet_class, dbr_used
#        <rates_out>  per sort: called barcodes, barcode space, occupancy,
#                     dbr used, cells tested, doublets called
#
# Sits beside celltype_doublets.R (fixed rate for the whole sublibrary); see
# workflow/rules/doublets_by_sort.smk for why. Same job shape and threading.
# =============================================================================

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  dirname(normalizePath(f))
}
source(file.path(script_dir(), "celltype_common.R"))

need_args(c("counts", "called_cells", "sublibraries", "sublib", "sort_wells",
            "barcodes_per_r1_well", "dbr_sd_rel", "min_cells",
            "out", "rates_out", "threads", "seed"))
pin_threads(as.integer(num("threads")))
seed      <- as.integer(num("seed"))
min_cells <- as.integer(num("min_cells"))
sd_rel    <- num("dbr_sd_rel")
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)

# "K_negative:43,TDP_negative:43,positive:10" -> named vector of round-1 wells
sw <- strsplit(csv("sort_wells"), ":", fixed = TRUE)
sort_wells <- setNames(as.integer(vapply(sw, `[`, "", 2)), vapply(sw, `[`, "", 1))
space <- sort_wells * num("barcodes_per_r1_well")

m <- readRDS(opt$counts)
step("read ", opt$counts, ": ", nrow(m), " genes x ", ncol(m), " cells")

# =============================================================================
# 1. Barcode occupancy per sort, from every split-pipe called barcode
# =============================================================================
# bc_wells ends in __s<N>. Take N from this sublibrary's own counts, which
# celltype_prep has already matched to the name; fall back to the row of
# sublibraries.tsv only if no cell survived prep.
hdr    <- names(fread(opt$called_cells, nrows = 0))
called <- fread(opt$called_cells, select = intersect(c("bc_wells", "sample", "bc1_well"), hdr))
if (!"bc1_well" %in% names(called)) called[, bc1_well := NA_character_]
suffix <- unique(sub(".*__s", "", colnames(m)))
if (length(suffix) == 0) {
  suffix <- as.character(match(opt$sublib, read.delim(opt$sublibraries, comment.char = "#")$sublibrary))
}
stopifnot(length(suffix) == 1)
called <- called[sub(".*__s", "", bc_wells) == suffix]
stopifnot(all(colnames(m) %in% called$bc_wells))
step(nrow(called), " split-pipe called barcodes in ", opt$sublib, " (__s", suffix, ")")

unknown <- setdiff(unique(called$sample), names(sort_wells))
if (length(unknown)) stop("Sorts in cell_metadata.csv not in the sample list: ", paste(unknown, collapse = ", "))

# Occupancy lambda from the fraction of barcodes seen; dbr = share of occupied
# barcodes holding more than one nucleus (zero-truncated Poisson).
rates <- called[, .(n_called_barcodes = .N, r1_wells_seen = uniqueN(bc1_well)), by = .(sort = sample)]
rates[, `:=`(r1_wells = sort_wells[sort], barcode_space = space[sort])]
rates[, lambda := -log(1 - n_called_barcodes / barcode_space)]
rates[, dbr := 1 - lambda * exp(-lambda) / (1 - exp(-lambda))]
rates[, dbr_sd := sd_rel * dbr]
rates[all(is.na(called$bc1_well)), r1_wells_seen := NA_integer_]
if (any(rates$r1_wells_seen > rates$r1_wells, na.rm = TRUE)) {
  stop("More round-1 wells seen than the sample list gives a sort -- check sort_wells")
}
print(rates)

# =============================================================================
# 2. scDblFinder, one call per sort
# =============================================================================
cell_sort <- called$sample[match(colnames(m), called$bc_wells)]
res <- rbindlist(lapply(rates$sort, function(s) {
  k <- which(cell_sort == s)
  r <- rates[sort == s]
  if (length(k) < min_cells) {
    warning(s, ": ", length(k), " cells is too few for scDblFinder (< ", min_cells, "); doublets not called")
    return(data.table(bc_wells = colnames(m)[k], sort = s, doublet_score = NA_real_,
                      doublet_class = NA_character_, dbr_used = r$dbr))
  }
  step(s, ": ", length(k), " cells, dbr = ", signif(r$dbr, 3), " +/- ", signif(r$dbr_sd, 3))
  set.seed(seed)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m[, k, drop = FALSE]))
  sce <- scDblFinder::scDblFinder(sce, dbr = r$dbr, dbr.sd = r$dbr_sd,
                                  BPPARAM = BiocParallel::SerialParam(RNGseed = seed),
                                  verbose = TRUE)
  out <- data.table(bc_wells      = colnames(sce),
                    sort          = s,
                    doublet_score = sce$scDblFinder.score,
                    doublet_class = as.character(sce$scDblFinder.class),
                    dbr_used      = r$dbr)
  rm(sce); invisible(gc())
  out
}))
stopifnot(nrow(res) == ncol(m), setequal(res$bc_wells, colnames(m)))

rates <- merge(rates,
               res[, .(n_tested = .N,
                       n_doublets = sum(doublet_class %in% "doublet"),
                       called_rate = if (all(is.na(doublet_class))) NA_real_
                                     else mean(doublet_class %in% "doublet")), by = sort],
               by = "sort", all.x = TRUE)
rates[is.na(n_tested), n_tested := 0L]
rates[, sublibrary := opt$sublib]
setcolorder(rates, "sublibrary")

fwrite(res, opt$out, sep = "\t")
fwrite(rates, opt$rates_out, sep = "\t")
print(rates[, .(sort, n_called_barcodes, dbr = signif(dbr, 3), n_tested, n_doublets,
                called_rate = signif(called_rate, 3))])
log_peak()
