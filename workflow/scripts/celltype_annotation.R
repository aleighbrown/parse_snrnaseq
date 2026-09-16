#!/usr/bin/env Rscript
# =============================================================================
# First-pass QC, Harmony clustering and cell type annotation of the combined
# split-pipe DGE (all 16 sublibraries, all-sample).
#
#   count_matrix.mtx (cells x genes) + donor assignment table
#     -> per-cell QC: genes, transcripts, % mitochondrial
#     -> scDblFinder, one sublibrary at a time
#     -> Seurat log-norm / HVG / PCA -> Harmony -> UMAP + Louvain clusters
#     -> Azimuth (humancortexref), one sublibrary at a time
#     -> marker-gene broad class per cluster
#     -> cell_annotations.csv: ONE ROW PER SPLIT-PIPE CELL, including the ones
#        QC removed, so every barcode can be traced to why it is (not) labelled
#
# Pipeline version of 01_qc_harmony_celltyping.Rmd. That notebook works on a
# 50k-cell subsample. This script takes all 426k cells in one large-memory job
# (see workflow/rules/celltype.smk for the sizing). Called as
#   Rscript celltype_annotation.R --key=value ...
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

# --- arguments ---------------------------------------------------------------
parse_args <- function(args) {
  bad <- !grepl("^--[^=]+=", args)
  if (any(bad)) stop("Arguments must look like --key=value, got: ", paste(args[bad], collapse = " "))
  setNames(as.list(sub("^--[^=]+=", "", args)), sub("^--([^=]+)=.*", "\\1", args))
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))

required <- c(
  "dge_dir", "donor_table", "sublibraries", "out_dir", "threads", "seed",
  "n_cells", "only_assigned",
  "min_genes", "max_genes", "min_counts", "max_counts", "max_pct_mt", "min_cells_per_gene",
  "doublet_dbr", "doublet_dbr_sd", "remove_doublets",
  "n_hvg", "n_pcs", "harmony_vars", "cluster_resolution", "extra_resolutions",
  "run_azimuth", "azimuth_reference", "azimuth_workers"
)
if (length(miss <- setdiff(required, names(opt)))) stop("Missing arguments: ", paste(miss, collapse = ", "))

num  <- function(k) as.numeric(opt[[k]])
flag <- function(k) tolower(opt[[k]]) %in% c("true", "1", "yes")
csv  <- function(k) { v <- trimws(strsplit(opt[[k]], ",")[[1]]); v[nzchar(v)] }
num_or_null <- function(k) if (tolower(opt[[k]]) %in% c("", "none", "null", "na")) NULL else num(k)

threads <- as.integer(num("threads"))
seed    <- as.integer(num("seed"))
n_pcs   <- as.integer(num("n_pcs"))
harmony_vars <- csv("harmony_vars")
resolutions  <- unique(c(as.numeric(csv("extra_resolutions")), num("cluster_resolution")))  # primary LAST
out <- function(f) file.path(opt$out_dir, f)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# --- threading ---------------------------------------------------------------
# BLAS stays single-threaded except around PCA: OpenBLAS thread pools and
# fork() (scDblFinder, future::multicore, mclapply) are a known deadlock mix.
setDTthreads(threads)
set_blas <- function(n) if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(n)
set_blas(1)
options(future.globals.maxSize = Inf)

# Seurat's RunUMAP / ScaleData take their thread count from the future plan.
with_workers <- function(expr) {
  if (threads > 1 && future::supportsMulticore()) {
    old <- future::plan(future::multicore, workers = threads)
    on.exit(future::plan(old), add = TRUE)
  }
  expr
}

# --- logging -----------------------------------------------------------------
# The GB figure is R's own heap after a gc -- use it (and the peak printed at the
# end) to trim config/cluster.yaml once a full run has reported back.
t0 <- Sys.time()
step <- function(...) {
  mem <- sum(gc(verbose = FALSE)[, 2]) / 1024
  message(sprintf("[%s | %6.1f min | %6.1f GB] ", format(Sys.time(), "%H:%M:%S"),
                  as.numeric(difftime(Sys.time(), t0, units = "mins")), mem), ...)
}
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
  fail_low_genes       = n_genes < num("min_genes"),
  fail_high_genes      = n_genes > num("max_genes"),
  fail_low_transcripts = n_transcripts < num("min_counts"),
  fail_high_transcripts = max_counts > 0 & n_transcripts > max_counts,
  fail_high_mt         = pct_mt > num("max_pct_mt")
)]
# scDblFinder wants low-quality cells gone but NOT high-count ones -- those are
# exactly the doublet-like cells it learns from. So the upper bounds wait.
meta[, pre_doublet_pass := !(fail_low_genes | fail_low_transcripts | fail_high_mt)]

keep <- meta$pre_doublet_pass
counts <- counts[, keep]
invisible(gc())
step(sum(keep), " of ", length(keep), " cells pass the low-quality filters")

# =============================================================================
# 4. Doublets (scDblFinder, split by sublibrary)
# =============================================================================
# Doublets in Parse are barcode collisions, which happen within a sublibrary
# (the sublibrary index is added after round 3), so that is the `samples` unit.
# scDblFinder's default rate (0.8% per 1000 cells, a 10x figure) would give
# ~20% here; the configured rate is used instead.
run_doublets <- function(counts, samples) {
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  n_workers <- min(threads, length(unique(samples)))
  bp <- if (n_workers > 1) {
    BiocParallel::MulticoreParam(workers = n_workers, RNGseed = seed)
  } else {
    BiocParallel::SerialParam(RNGseed = seed)
  }
  sce <- scDblFinder::scDblFinder(sce, samples = samples,
                                  dbr = num_or_null("doublet_dbr"),
                                  dbr.sd = num_or_null("doublet_dbr_sd"),
                                  BPPARAM = bp, verbose = FALSE)
  data.table(doublet_score = sce$scDblFinder.score,
             doublet_class = as.character(sce$scDblFinder.class))
}

set.seed(seed)
dbl <- run_doublets(counts, samples = as.character(meta$sublibrary[keep]))
meta[keep, `:=`(doublet_score = dbl$doublet_score, doublet_class = dbl$doublet_class)]
meta[, fail_doublet := flag("remove_doublets") & doublet_class %in% "doublet"]
rm(dbl); invisible(gc())
step("doublets: ", sum(meta$doublet_class %in% "doublet"), " called (",
     round(100 * mean(meta$doublet_class[keep] == "doublet"), 2), "% of tested cells)")

# --- final QC call + what goes into clustering --------------------------------
fail_cols <- grep("^fail_", names(meta), value = TRUE)
meta[, qc_pass := !Reduce(`|`, .SD), .SDcols = fail_cols]
meta[, clustered := qc_pass & (assigned | !only_assigned)]

reason <- character(nrow(meta))
for (f in fail_cols) reason <- paste0(reason, ifelse(meta[[f]], paste0(sub("^fail_", "", f), ";"), ""))
reason <- paste0(reason, ifelse(meta$qc_pass & !meta$clustered, "no_donor;", ""))
meta[, exclude_reason := sub(";$", "", reason)]
rm(reason)

qc_summary <- meta[, c(list(n_cells = .N), lapply(.SD, sum),
                       list(qc_pass = sum(qc_pass), clustered = sum(clustered),
                            median_transcripts = as.numeric(median(n_transcripts)),
                            median_genes = as.numeric(median(n_genes)),
                            median_pct_mt = median(pct_mt))),
                   by = sublibrary, .SDcols = fail_cols]
qc_summary <- rbind(qc_summary,
                    qc_summary[, c(list(sublibrary = "TOTAL"), lapply(.SD, sum)),
                               .SDcols = c("n_cells", fail_cols, "qc_pass", "clustered")],
                    fill = TRUE)
fwrite(qc_summary, out("qc_summary.tsv"), sep = "\t")
print(qc_summary)

# counts currently holds the pre_doublet_pass cells.
counts <- counts[, meta$clustered[keep]]
invisible(gc())
step(ncol(counts), " cells go into clustering")

# =============================================================================
# 5. Seurat -> Harmony -> clusters
# =============================================================================
md <- as.data.frame(meta[clustered == TRUE])
rownames(md) <- md$bc_wells
md$sublibrary <- factor(md$sublibrary, levels = sublibs)
obj <- CreateSeuratObject(counts, meta.data = md, project = "parse_celltype",
                          min.cells = as.integer(num("min_cells_per_gene")), min.features = 0)
rm(counts, md); invisible(gc())
step("Seurat object: ", nrow(obj), " genes x ", ncol(obj), " cells")

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = as.integer(num("n_hvg")), verbose = FALSE)
obj <- with_workers(ScaleData(obj, verbose = FALSE))
set_blas(threads)
obj <- RunPCA(obj, npcs = max(50, n_pcs), seed.use = seed, verbose = FALSE)
set_blas(1)
step("PCA done")

obj <- with_workers(RunUMAP(obj, reduction = "pca", dims = 1:n_pcs, seed.use = seed,
                            reduction.name = "umap_uncorrected", reduction.key = "uUMAP_", verbose = FALSE))
step("uncorrected UMAP done")

set.seed(seed)
obj <- RunHarmony(obj, group.by.vars = harmony_vars, reduction.use = "pca",
                  dims.use = 1:n_pcs, reduction.save = "harmony",
                  ncores = threads, verbose = FALSE)
step("Harmony done (", paste(harmony_vars, collapse = " + "), ")")

obj <- with_workers(RunUMAP(obj, reduction = "harmony", dims = 1:n_pcs, seed.use = seed,
                            reduction.name = "umap_harmony", reduction.key = "hUMAP_", verbose = FALSE))
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
obj <- FindClusters(obj, resolution = resolutions, random.seed = seed, verbose = FALSE)
Idents(obj) <- "seurat_clusters"                            # = the primary resolution
step(nlevels(obj$seurat_clusters), " clusters at resolution ", num("cluster_resolution"))

# Share of each embedding dimension explained by a variable (between-group / total
# variance, i.e. the R^2 of a one-way ANOVA).
embedding_r2 <- function(obj, reduction, vars, dims) {
  emb <- Embeddings(obj, reduction)[, dims, drop = FALSE]
  rbindlist(lapply(vars, function(v) {
    g <- as.character(obj@meta.data[[v]])
    data.table(embedding = reduction, variable = v, dim = dims,
               r2 = apply(emb, 2, function(p) sum((ave(p, g) - mean(p))^2) / sum((p - mean(p))^2)))
  }))
}
r2_vars <- intersect(c("sublibrary", "donor_id", "sample"), colnames(obj@meta.data))
r2 <- rbind(embedding_r2(obj, "pca", r2_vars, 1:n_pcs),
            embedding_r2(obj, "harmony", r2_vars, 1:n_pcs))

# =============================================================================
# 6. Azimuth reference mapping, one sublibrary per call
# =============================================================================
# Mapping all cells at once would build a dense reference-feature residual
# matrix for every cell; per sublibrary it stays at ~26k cells a go, and the
# chunks run `azimuth_workers` at a time.
pred_cols <- character(0)
if (flag("run_azimuth")) {
  cnt <- LayerData(obj, assay = "RNA", layer = "counts")
  chunks <- split(seq_len(ncol(cnt)), droplevels(obj$sublibrary))
  map_chunk <- function(k) {
    q <- CreateSeuratObject(cnt[, k, drop = FALSE], min.cells = 0, min.features = 0)
    q <- suppressWarnings(Azimuth::RunAzimuth(q, reference = opt$azimuth_reference, verbose = FALSE))
    qm <- q@meta.data
    keep_cols <- grep("^predicted\\.|^mapping\\.score$", colnames(qm), value = TRUE)
    data.table(bc_wells = rownames(qm), qm[, keep_cols, drop = FALSE])
  }
  workers <- min(as.integer(num("azimuth_workers")), length(chunks))
  step("Azimuth: ", length(chunks), " chunks, ", workers, " at a time, reference ", opt$azimuth_reference)
  az <- if (workers > 1) {
    parallel::mclapply(chunks, map_chunk, mc.cores = workers, mc.preschedule = FALSE)
  } else {
    lapply(chunks, map_chunk)
  }
  # A worker killed by the OOM killer comes back as NULL, not as an error.
  failed <- vapply(az, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
  if (any(failed)) {
    stop("Azimuth failed for: ", paste(names(az)[failed], collapse = ", "), "\n",
         paste(unique(vapply(az[failed], function(r) if (is.null(r)) "worker died (out of memory?)" else as.character(r), "")), collapse = "\n"))
  }
  az <- rbindlist(az, use.names = TRUE, fill = TRUE)
  rm(cnt); invisible(gc())

  pred_cols <- setdiff(names(az), "bc_wells")
  idx <- match(colnames(obj), az$bc_wells)
  stopifnot(!anyNA(idx))
  for (p in pred_cols) obj@meta.data[[p]] <- az[[p]][idx]
  rm(az); invisible(gc())
  step("Azimuth done: ", paste(pred_cols, collapse = ", "))
}

# =============================================================================
# 7. Cluster-level labels
# =============================================================================
cl_meta <- as.data.table(obj@meta.data, keep.rownames = "cell")

# 7a. Majority Azimuth label per cluster. Per-cell labels are noisy at Parse
# depth; the cluster vote is the steadier call.
majority <- function(dt, label) {
  dt[, .N, by = c("seurat_clusters", label)][
    , frac := N / sum(N), by = seurat_clusters][
    order(-N), .SD[1], by = seurat_clusters][
    , .(seurat_clusters, label = get(label), frac)]
}
cluster_summary <- cl_meta[, .(n_cells = .N,
                               top_donor = names(which.max(table(donor_id))),
                               top_donor_frac = max(table(donor_id)) / .N),
                           by = seurat_clusters]
for (lvl in intersect(c("predicted.class", "predicted.subclass"), pred_cols)) {
  short <- sub("^predicted\\.", "", lvl)
  m <- majority(cl_meta, lvl)
  setnames(m, c("label", "frac"), paste0("azimuth_", short, c("", "_frac")))
  cluster_summary <- merge(cluster_summary, m, by = "seurat_clusters")
}

# 7b. Broad class from marker genes (offline; no reference needed).
markers <- list(
  Exc_neuron      = c("SLC17A7", "SATB2", "CAMK2A"),
  Inh_neuron      = c("GAD1", "GAD2", "SLC32A1"),
  Astrocyte       = c("AQP4", "GFAP", "SLC1A2", "GJA1"),
  Oligodendrocyte = c("MBP", "MOBP", "PLP1", "ST18"),
  OPC             = c("PDGFRA", "VCAN", "CSPG4"),
  Microglia       = c("CSF1R", "P2RY12", "CX3CR1", "C3"),
  Endothelial     = c("CLDN5", "FLT1"),
  Mural_VLMC      = c("PDGFRB", "DCN")
)
markers <- Filter(length, lapply(markers, intersect, rownames(obj)))
obj <- AddModuleScore(obj, features = markers, name = "ms_", seed = seed)
score_cols <- paste0("ms_", seq_along(markers))
sc <- as.data.table(obj@meta.data)[, lapply(.SD, mean), by = seurat_clusters, .SDcols = score_cols]
score_mat <- as.matrix(sc[, ..score_cols])
top2 <- t(apply(score_mat, 1, sort, decreasing = TRUE))[, 1:2, drop = FALSE]
cluster_summary <- merge(cluster_summary,
                         data.table(seurat_clusters = sc$seurat_clusters,
                                    marker_class = names(markers)[max.col(score_mat, ties.method = "first")],
                                    marker_margin = top2[, 1] - top2[, 2]),
                         by = "seurat_clusters")
obj@meta.data[score_cols] <- NULL

# 7c. Top markers per cluster (presto is used automatically when installed).
cl_markers <- FindAllMarkers(obj, only.pos = TRUE, max.cells.per.ident = 300,
                             logfc.threshold = 0.5, min.pct = 0.25,
                             random.seed = seed, verbose = FALSE)
fwrite(cl_markers, out("cluster_markers.tsv"), sep = "\t")
if (nrow(cl_markers)) {
  top <- as.data.table(cl_markers)[order(-avg_log2FC)][
    , .(top_markers = paste(head(gene, 8), collapse = ", ")), by = .(seurat_clusters = cluster)]
  cluster_summary <- merge(cluster_summary, top, by = "seurat_clusters", all.x = TRUE)
}
setorder(cluster_summary, seurat_clusters)
fwrite(cluster_summary, out("cluster_summary.tsv"), sep = "\t")
print(cluster_summary)
step("cluster labels done")

# Cluster labels back onto cells.
lab_cols <- setdiff(names(cluster_summary), c("seurat_clusters", "n_cells", "top_donor", "top_donor_frac", "top_markers"))
idx <- match(obj$seurat_clusters, cluster_summary$seurat_clusters)
for (l in lab_cols) obj@meta.data[[paste0("cluster_", l)]] <- cluster_summary[[l]][idx]

# =============================================================================
# 8. Report
# =============================================================================
pdf(out("celltype_report.pdf"), width = 14, height = 8)
theme_set(theme_minimal(base_size = 11))
rot_x <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

qc_long <- melt(meta[, .(sublibrary, n_genes, n_transcripts, pct_mt)],
                id.vars = "sublibrary", variable.factor = FALSE)
qc_lines <- data.table(variable = c("n_genes", "n_genes", "n_transcripts", "n_transcripts", "pct_mt"),
                       y = c(num("min_genes"), num("max_genes"), num("min_counts"),
                             if (max_counts > 0) max_counts else NA, num("max_pct_mt")))
print(ggplot(qc_long, aes(factor(sublibrary, levels = sublibs), value)) +
        geom_violin(scale = "width", linewidth = 0.2, fill = "grey80") +
        geom_hline(data = qc_lines[!is.na(y)], aes(yintercept = y), linetype = "dashed", colour = "red") +
        scale_y_continuous(trans = "log1p") +
        facet_wrap(~ variable, ncol = 1, scales = "free_y") +
        labs(x = NULL, y = NULL, title = "Per-cell QC, all split-pipe cells (dashed = thresholds)") + rot_x)
rm(qc_long)

fail_long <- melt(qc_summary[sublibrary != "TOTAL", c("sublibrary", fail_cols), with = FALSE],
                  id.vars = "sublibrary", variable.factor = FALSE)
print(ggplot(fail_long, aes(factor(sublibrary, levels = sublibs), value, fill = variable)) +
        geom_col(position = "dodge") +
        labs(x = NULL, y = "cells failing", fill = NULL,
             title = "Cells failing each filter (a cell can fail several)") + rot_x)

print(ggplot(meta[!is.na(doublet_score)], aes(factor(sublibrary, levels = sublibs), doublet_score)) +
        geom_violin(scale = "width", linewidth = 0.2, fill = "grey80") +
        labs(x = NULL, y = "scDblFinder score", title = "Doublet scores by sublibrary") + rot_x)

print(ggplot(r2, aes(dim, r2, colour = variable)) +
        geom_line() + geom_point(size = 1) + facet_wrap(~ embedding) +
        labs(x = "dimension", y = expression(R^2), title = "Variance per dimension explained, before vs after Harmony"))

umap_pair <- function(var, legend = TRUE) {
  p <- lapply(c("umap_uncorrected", "umap_harmony"), function(red) {
    DimPlot(obj, reduction = red, group.by = var, shuffle = TRUE, raster = TRUE) +
      ggtitle(paste(sub("umap_", "", red), ":", var))
  })
  p[[1]] <- p[[1]] + NoLegend()
  if (!legend) p[[2]] <- p[[2]] + NoLegend()
  p[[1]] | p[[2]]
}
print(umap_pair("sublibrary"))
print(umap_pair("donor_id", legend = length(unique(obj$donor_id)) <= 40))
print(umap_pair("sample"))

print(DimPlot(obj, reduction = "umap_harmony", group.by = "seurat_clusters", label = TRUE, raster = TRUE) +
        NoLegend() + ggtitle("Harmony: clusters") |
        DimPlot(obj, reduction = "umap_harmony", group.by = "cluster_marker_class",
                label = TRUE, repel = TRUE, raster = TRUE) + ggtitle("Harmony: marker class per cluster"))

if ("predicted.subclass" %in% pred_cols) {
  print(DimPlot(obj, reduction = "umap_harmony", group.by = "predicted.subclass",
                label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() +
          ggtitle("Harmony: Azimuth subclass (per cell)") |
          FeaturePlot(obj, reduction = "umap_harmony", features = "mapping.score", raster = TRUE))
  print(DimPlot(obj, reduction = "umap_harmony", group.by = "cluster_azimuth_subclass",
                label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() +
          ggtitle("Harmony: Azimuth subclass (cluster majority)"))
}

print(DotPlot(obj, features = markers, group.by = "seurat_clusters", cluster.idents = TRUE) +
        RotatedAxis() + labs(x = NULL, y = "cluster"))

comp <- as.data.table(obj@meta.data)[, .N, by = .(seurat_clusters, donor_id)][, frac := N / sum(N), by = seurat_clusters]
print(ggplot(comp, aes(seurat_clusters, frac, fill = donor_id)) +
        geom_col(width = 0.9) +
        labs(x = "cluster", y = "fraction of cells", title = "Donor make-up of each cluster") +
        guides(fill = guide_legend(ncol = 2)))
invisible(dev.off())
step("report written")

# =============================================================================
# 9. Outputs
# =============================================================================
# One row per split-pipe cell. QC-failed / unclustered cells keep their QC
# columns and get NA annotation, with exclude_reason saying why.
obj_meta <- as.data.table(obj@meta.data, keep.rownames = "bc_wells")
ann_cols <- c(grep("^RNA_snn_res\\.", names(obj_meta), value = TRUE), "seurat_clusters",
              pred_cols, grep("^cluster_", names(obj_meta), value = TRUE))
emb <- cbind(Embeddings(obj, "umap_harmony"), Embeddings(obj, "umap_uncorrected"))
colnames(emb) <- c("umap_harmony_1", "umap_harmony_2", "umap_uncorrected_1", "umap_uncorrected_2")
ann <- cbind(obj_meta[, c("bc_wells", ann_cols), with = FALSE], as.data.table(emb))
ann[, seurat_clusters := as.character(seurat_clusters)]

meta[, cell_id := paste0(sublibrary, "|", barcode)]            # = demux table's cell_id
meta[, collab_barcode := paste0(barcode, "__s", sublibrary)]   # = collaborator's barcode
ann <- merge(meta, ann, by = "bc_wells", all.x = TRUE, sort = FALSE)
setorder(ann, mtx_row)

id_cols <- c("bc_wells", "cell_id", "collab_barcode", "barcode", "sublibrary", "sub_n", "mtx_row",
             "bc1_well", "bc2_well", "bc3_well", "sample", "donor_id", "assigned")
qc_out  <- c("n_transcripts", "n_genes", "pct_mt", "pct_ribo", "doublet_score", "doublet_class",
             "qc_pass", "clustered", "exclude_reason")
first <- intersect(c(id_cols, qc_out, ann_cols, colnames(emb)), names(ann))
setcolorder(ann, first)
ann[, c("pre_doublet_pass", fail_cols) := NULL]
fwrite(ann, out("cell_annotations.csv"))
step("wrote cell_annotations.csv: ", nrow(ann), " cells x ", ncol(ann), " columns")

# The Harmony embedding is for clustering and labels only -- pseudobulk from the
# raw counts layer. scale.data and the SNN graphs are dropped to shrink the file.
obj[["RNA"]]$scale.data <- NULL
obj@graphs <- list()
saveRDS(obj, out("seurat_harmony.rds"), compress = FALSE)
step("wrote seurat_harmony.rds")

peak <- sum(gc(verbose = FALSE)[, 6]) / 1024
step(sprintf("done. Peak R heap %.1f GB (the job's RSS will be higher: forks + fread buffers)", peak))
print(sessionInfo())
