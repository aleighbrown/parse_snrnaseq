#!/usr/bin/env Rscript
# =============================================================================
# Cell typing, stage 3 of 3: final QC call, clustering, labels, outputs.
#
#   <prep>/cell_meta.rds + <prep>/counts/*.rds + <doublets>/*.tsv
#     -> doublet + upper-bound filters -> qc_summary.tsv
#     -> Seurat log-norm / HVG / PCA -> Harmony -> UMAP + Louvain clusters
#     -> Azimuth (humancortexref), one sublibrary per worker process
#     -> marker-gene broad class per cluster
#     -> cell_annotations.csv: ONE ROW PER SPLIT-PIPE CELL, including the ones
#        QC removed, so every barcode can be traced to why it is (not) labelled
#
# Optional arguments, used by the SoupX re-clustering (celltype_cluster_soupx
# in soupx.smk); without them this script behaves exactly as before:
#   --counts_dir=     read counts/<sublib>.rds from here instead of prep_dir
#                     (SoupX-corrected counts). QC metrics stay on the raw
#                     counts in cell_meta.rds.
#   --floor_files=    comma-separated soup/<sort>/barcode_summary.tsv
#   --min_floor_k=    fail cells whose raw UMIs are < k x the background floor
#                     of their sort x sublibrary (fail_below_floor)
#   --soupx_cells_dir= soupx/cells/<sublib>.tsv, merged into the output
#
# Pipeline version of 01_qc_harmony_celltyping.Rmd. Clustering is one
# in-memory computation over every cell, so this stage is the one large job;
# the slow, per-sublibrary doublet calling happens before it in small jobs.
# See workflow/rules/celltype.smk.
# =============================================================================

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  dirname(normalizePath(f))
}
source(file.path(script_dir(), "celltype_common.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

need_args(c("prep_dir", "doublet_dir", "sublibraries", "out_dir", "threads", "seed",
            "only_assigned", "min_genes", "max_genes", "min_counts", "max_counts", "max_pct_mt",
            "min_cells_per_gene", "remove_doublets",
            "n_hvg", "n_pcs", "harmony_vars", "cluster_resolution", "extra_resolutions",
            "run_azimuth", "azimuth_reference", "azimuth_workers"))
threads <- as.integer(num("threads"))
seed    <- as.integer(num("seed"))
n_pcs   <- as.integer(num("n_pcs"))
harmony_vars  <- csv("harmony_vars")
resolutions   <- unique(c(as.numeric(csv("extra_resolutions")), num("cluster_resolution")))  # primary LAST
only_assigned <- flag("only_assigned")
max_counts    <- num("max_counts")
setup_threads(threads)
out <- function(f) file.path(opt$out_dir, f)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
step("threads=", threads, "  out_dir=", opt$out_dir)

sublibs <- read.delim(opt$sublibraries, comment.char = "#")$sublibrary
meta <- readRDS(file.path(opt$prep_dir, "cell_meta.rds"))
step("cell_meta.rds: ", nrow(meta), " cells")

given <- function(k) !is.null(opt[[k]]) && !tolower(opt[[k]]) %in% c("", "none", "null")
counts_dir <- if (given("counts_dir")) opt$counts_dir else opt$prep_dir
step("counts from ", counts_dir)

# SoupX per-cell summary (rho, UMIs removed), carried into the output.
if (given("soupx_cells_dir")) {
  sx <- rbindlist(lapply(sublibs, function(s) fread(file.path(opt$soupx_cells_dir, paste0(s, ".tsv")))))
  meta[sx, `:=`(soupx_rho = i.rho, n_transcripts_corrected = i.n_umis_corrected,
                soupx_frac_removed = i.frac_removed), on = "bc_wells"]
  rm(sx)
  step("SoupX: ", sum(!is.na(meta$soupx_rho)), " cells with corrected counts")
}

# =============================================================================
# 4. Doublet calls + final QC
# =============================================================================
dbl <- rbindlist(lapply(sublibs, function(s) {
  fread(doublet_file(opt$doublet_dir, s), colClasses = c(doublet_class = "character"))
}))
idx <- match(meta$bc_wells, dbl$bc_wells)
stopifnot(identical(sort(which(!is.na(idx))), which(meta$pre_doublet_pass)))
meta[, `:=`(doublet_score = dbl$doublet_score[idx], doublet_class = dbl$doublet_class[idx])]
meta[, fail_doublet := flag("remove_doublets") & doublet_class %in% "doublet"]
rm(dbl, idx); invisible(gc())
step("doublets: ", sum(meta$doublet_class %in% "doublet"), " called (",
     round(100 * mean(meta$doublet_class[meta$pre_doublet_pass] %in% "doublet"), 2), "% of tested cells)")

# Minimum size relative to the background floor: a cell with k x its floor
# carries at most 1/k background from the per-barcode floor alone.
if (given("floor_files") && given("min_floor_k")) {
  fl <- rbindlist(lapply(csv("floor_files"), fread))[
    , .(sample = sort, sublibrary, floor_umis = uncalled_median_umis)]
  meta[fl, floor_umis := i.floor_umis, on = .(sample, sublibrary)]
  if (anyNA(meta$floor_umis)) stop(sum(is.na(meta$floor_umis)), " cells have no floor for their sort x sublibrary")
  meta[, floor_k := n_transcripts / floor_umis]
  meta[, fail_below_floor := floor_k < num("min_floor_k")]
  step("below ", num("min_floor_k"), " x floor: ", sum(meta$fail_below_floor), " cells")
  rm(fl)
}

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
                            doublet_rate = mean(doublet_class[pre_doublet_pass] %in% "doublet"),
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

# --- counts of the clustered cells -------------------------------------------
# Stitch the per-sublibrary matrices together by hand: Matrix's cbind() on
# many pieces works pairwise and copies the growing result each time.
cbind_csc <- function(ms) {
  stopifnot(all(vapply(ms, function(m) identical(rownames(m), rownames(ms[[1]])), logical(1))))
  nnz <- vapply(ms, function(m) length(m@x), numeric(1))
  if (sum(nnz) > .Machine$integer.max) stop("more than 2^31-1 non-zeros in the clustered cells")
  p <- c(0L, unlist(lapply(seq_along(ms), function(j) {
    ms[[j]]@p[-1] + as.integer(sum(nnz[seq_len(j - 1)]))
  })))
  new("dgCMatrix",
      i = unlist(lapply(ms, slot, "i")), p = p, x = unlist(lapply(ms, slot, "x")),
      Dim = c(nrow(ms[[1]]), sum(vapply(ms, ncol, integer(1)))),
      Dimnames = list(rownames(ms[[1]]), unlist(lapply(ms, colnames))))
}

clustered_cells <- meta$bc_wells[meta$clustered]
pieces <- lapply(sublibs, function(s) {
  m <- readRDS(counts_file(counts_dir, s))
  m[, colnames(m) %in% clustered_cells, drop = FALSE]
})
counts <- cbind_csc(pieces)
rm(pieces); invisible(gc())
step(ncol(counts), " cells go into clustering")

# =============================================================================
# 5. Seurat -> Harmony -> clusters
# =============================================================================
md <- as.data.frame(meta[match(colnames(counts), bc_wells)])
rownames(md) <- md$bc_wells
md$sublibrary <- factor(md$sublibrary, levels = sublibs)
obj <- CreateSeuratObject(counts, meta.data = md, project = "parse_celltype",
                          min.cells = as.integer(num("min_cells_per_gene")), min.features = 0)
rm(counts, md); invisible(gc())
step("Seurat object: ", nrow(obj), " genes x ", ncol(obj), " cells")

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = as.integer(num("n_hvg")), verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)            # sequential: see with_workers() in celltype_common.R
set_blas(threads)
obj <- RunPCA(obj, npcs = max(50, n_pcs), seed.use = seed, verbose = FALSE)
set_blas(1)
step("PCA done")

obj <- with_workers(threads, RunUMAP(obj, reduction = "pca", dims = 1:n_pcs, seed.use = seed,
                            reduction.name = "umap_uncorrected", reduction.key = "uUMAP_", verbose = FALSE))
step("uncorrected UMAP done")

set.seed(seed)
obj <- RunHarmony(obj, group.by.vars = harmony_vars, reduction.use = "pca",
                  dims.use = 1:n_pcs, reduction.save = "harmony",
                  ncores = threads, verbose = FALSE)
step("Harmony done (", paste(harmony_vars, collapse = " + "), ")")

obj <- with_workers(threads, RunUMAP(obj, reduction = "harmony", dims = 1:n_pcs, seed.use = seed,
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
# chunks run `azimuth_workers` at a time, each in a fresh R process (see
# "threading" in celltype_common.R). Each worker reads its own sublibrary's
# counts file, so nothing large is sent to it. A worker that dies (e.g. out of
# memory) makes bplapply stop with an error naming the lost worker.
# Self-contained: a Snow worker is a separate R session that sees none of this
# script's functions, so the thread pinning from pin_threads() is repeated.
azimuth_worker <- function(file, cells, reference) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
  suppressPackageStartupMessages(library(Seurat))
  m <- readRDS(file)
  q <- CreateSeuratObject(m[, cells, drop = FALSE], min.cells = 0, min.features = 0)
  rm(m)
  q <- suppressWarnings(Azimuth::RunAzimuth(q, reference = reference, verbose = FALSE))
  qm <- q@meta.data
  keep_cols <- grep("^predicted\\.|^mapping\\.score$", colnames(qm), value = TRUE)
  data.frame(bc_wells = rownames(qm), qm[, keep_cols, drop = FALSE], check.names = FALSE)
}

pred_cols <- character(0)
if (flag("run_azimuth")) {
  cells_by_sub <- split(colnames(obj), as.character(obj$sublibrary))
  cells_by_sub <- cells_by_sub[intersect(sublibs, names(cells_by_sub))]
  workers <- min(as.integer(num("azimuth_workers")), threads, length(cells_by_sub))
  step("Azimuth: ", length(cells_by_sub), " sublibraries, ", workers, " at a time, reference ", opt$azimuth_reference)
  az <- BiocParallel::bpmapply(
    azimuth_worker,
    file = vapply(names(cells_by_sub), function(s) counts_file(counts_dir, s), ""),
    cells = cells_by_sub,
    MoreArgs = list(reference = opt$azimuth_reference),
    SIMPLIFY = FALSE,
    BPPARAM = snow_param(workers, seed)
  )
  az <- rbindlist(az, use.names = TRUE, fill = TRUE)
  rm(cells_by_sub); invisible(gc())

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

log_peak()
print(sessionInfo())
