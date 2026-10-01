# =============================================================================
# Canonical marker dot plots by Azimuth predicted.subclass, split by FANS sort,
# on the SoupX-corrected clustering (celltype_soupx/seurat_harmony.rds).
#
#   Rscript marker_dotplot.R --rds=<seurat_harmony.rds> --out_dir=<dir> \
#       [--min_score=0.8] [--min_cells=10]
#
# Two cell sets, each written as <out_dir>/marker_dotplot_<set>.{pdf,tsv}:
#   all_qc      qc_pass & doublet_class == "singlet"
#   score_gt08  the above & predicted.subclass.score > min_score
#
# Dot size  = fraction of cells with >= 1 (SoupX-corrected) count of the gene.
# Dot color = mean log-normalised expression, z-scored per gene across every
#             subclass x sort group in that plot, so sorts share one scale.
# Groups with fewer than min_cells cells are dropped (listed in the log).
#
# Each PDF: page 1 all three sorts stacked (shared gene axis), then one page
# per sort.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

parse_args <- function(args) {
  bad <- !grepl("^--[^=]+=", args)
  if (any(bad)) stop("Arguments must look like --key=value, got: ", paste(args[bad], collapse = " "))
  setNames(as.list(sub("^--[^=]+=", "", args)), sub("^--([^=]+)=.*", "\\1", args))
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("rds", "out_dir")) if (is.null(opt[[k]])) stop("Missing --", k)
min_score <- as.numeric(if (is.null(opt$min_score)) 0.8 else opt$min_score)
min_cells <- as.integer(if (is.null(opt$min_cells)) 10 else opt$min_cells)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()
step <- function(...) message(sprintf("[%s | %5.1f min] ", format(Sys.time(), "%H:%M:%S"),
                                      as.numeric(difftime(Sys.time(), t0, units = "mins"))), ...)

# --- markers -----------------------------------------------------------------
# Grouped by the cell type they mark; plotted in this order. RBFOX3 = NeuN,
# the sorting antigen.
markers <- list(
  Neuron       = c("RBFOX3", "SNAP25", "SYT1"),
  Excitatory   = c("SLC17A7", "SATB2", "CAMK2A"),
  `Exc layer`  = c("CUX2", "RORB", "FEZF2", "THEMIS", "FOXP2", "NXPH4"),
  Inhibitory   = c("GAD1", "GAD2", "SLC32A1"),
  `Inh class`  = c("LAMP5", "VIP", "SST", "PVALB"),
  Astrocyte    = c("AQP4", "GFAP", "SLC1A2"),
  Oligo        = c("MBP", "MOBP", "PLP1"),
  OPC          = c("PDGFRA", "VCAN"),
  Microglia    = c("CSF1R", "P2RY12", "CX3CR1"),
  Endothelial  = c("CLDN5", "FLT1"),
  `Mural/VLMC` = c("PDGFRB", "RGS5", "DCN", "COL1A2")
)
gene_group <- setNames(rep(names(markers), lengths(markers)), unlist(markers))

# Azimuth humancortexref subclasses, neurons then glia; anything else is
# appended alphabetically.
subclass_order <- c("L2/3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5/6 NP",
                    "L6 CT", "L6b", "Lamp5", "Sncg", "Vip", "Sst Chodl", "Sst",
                    "Pvalb", "Chandelier", "Astro", "Oligo", "OPC", "Micro-PVM",
                    "Endo", "VLMC", "Pericyte", "SMC")
sort_levels <- c("positive", "TDP_negative", "K_negative")

# --- load, keep only what the plots need --------------------------------------
step("reading ", opt$rds)
obj <- readRDS(opt$rds)
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj, search = "counts")) > 1) obj <- JoinLayers(obj)

genes <- intersect(names(gene_group), rownames(obj))
missing <- setdiff(names(gene_group), genes)
if (length(missing)) step("markers not in the object, skipped: ", paste(missing, collapse = ", "))

need <- c("sample", "qc_pass", "doublet_class", "predicted.subclass", "predicted.subclass.score")
miss_cols <- setdiff(need, colnames(obj@meta.data))
if (length(miss_cols)) stop("metadata columns missing: ", paste(miss_cols, collapse = ", "))

meta <- as.data.table(obj@meta.data[, need], keep.rownames = "cell")
counts <- LayerData(obj, assay = "RNA", layer = "counts")[genes, , drop = FALSE]
lognorm <- LayerData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
rm(obj); invisible(gc())
step(nrow(meta), " cells in the object, ", length(genes), " marker genes")

base <- meta$qc_pass %in% TRUE & meta$doublet_class %in% "singlet" & !is.na(meta$predicted.subclass)
sets <- list(
  all_qc     = list(keep = base,
                    title = "QC-pass singlets"),
  score_gt08 = list(keep = base & meta$predicted.subclass.score > min_score,
                    title = sprintf("QC-pass singlets, predicted.subclass.score > %s", min_score))
)

# --- per group summary ---------------------------------------------------------
summarise_set <- function(keep) {
  m <- meta[keep]
  grp <- paste(m$sample, m$predicted.subclass, sep = "\r")
  ids <- split(match(m$cell, colnames(counts)), grp)
  rbindlist(lapply(names(ids), function(g) {
    j <- ids[[g]]
    data.table(sort = sub("\r.*", "", g), subclass = sub(".*\r", "", g),
               n_cells = length(j), gene = genes,
               pct_expr = Matrix::rowSums(counts[, j, drop = FALSE] > 0) / length(j),
               mean_expr = Matrix::rowMeans(lognorm[, j, drop = FALSE]))
  }))
}

# One hue light -> dark for the z-scored mean (sequential, not diverging: low
# expression is "absent", not an opposite pole).
seq_cols <- c("#E8EEF7", "#9DB6DB", "#4A78B8", "#1F4A8A", "#0B2350")

dot_plot <- function(d, title) {
  ggplot(d, aes(gene, subclass_lab)) +
    geom_point(aes(size = pct_expr, colour = z), shape = 16) +
    scale_size_area(name = "Fraction of cells\nexpressing", max_size = 5,
                    limits = c(0, 1), breaks = c(0.1, 0.25, 0.5, 0.75, 1)) +
    scale_colour_gradientn(name = "Mean expression\n(z-score per gene)",
                           colours = seq_cols, limits = c(-2.5, 2.5), oob = scales::squish) +
    facet_grid(sort ~ group, scales = "free", space = "free") +
    labs(title = title, x = NULL, y = "Azimuth predicted.subclass (n cells)") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
          strip.text.x = element_text(angle = 90, hjust = 0, size = 8),
          strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(2, "pt"),
          plot.title = element_text(face = "bold"),
          legend.title = element_text(size = 8))
}

for (nm in names(sets)) {
  s <- sets[[nm]]
  step(nm, ": ", sum(s$keep), " cells")
  d <- summarise_set(s$keep)

  small <- unique(d[n_cells < min_cells, .(sort, subclass, n_cells)])
  if (nrow(small)) step(nm, ": dropping groups with < ", min_cells, " cells: ",
                        paste(sprintf("%s/%s (%d)", small$sort, small$subclass, small$n_cells), collapse = ", "))
  d <- d[n_cells >= min_cells]

  d[, z := { sdv <- sd(mean_expr); if (is.na(sdv) || sdv == 0) 0 else (mean_expr - mean(mean_expr)) / sdv }, by = gene]
  d[, group := factor(gene_group[gene], levels = names(markers))]
  d[, gene := factor(gene, levels = genes)]
  d[, sort := factor(sort, levels = c(intersect(sort_levels, unique(sort)), setdiff(sort(unique(sort)), sort_levels)))]
  sub_lv <- c(intersect(subclass_order, unique(d$subclass)), sort(setdiff(unique(d$subclass), subclass_order)))
  d[, subclass := factor(subclass, levels = sub_lv)]
  # y label carries n per sort; order rows by subclass, top to bottom
  d[, subclass_lab := sprintf("%s (%s)", subclass, format(n_cells, big.mark = ","))]
  lab_lv <- unique(d[order(sort, -as.integer(subclass)), subclass_lab])
  d[, subclass_lab := factor(subclass_lab, levels = lab_lv)]

  fwrite(d[order(sort, subclass, gene), .(sort, subclass, n_cells, group, gene, pct_expr, mean_expr, z)],
         file.path(opt$out_dir, paste0("marker_dotplot_", nm, ".tsv")), sep = "\t")

  n_rows <- uniqueN(d[, .(sort, subclass)])
  pdf(file.path(opt$out_dir, paste0("marker_dotplot_", nm, ".pdf")),
      width = 3 + 0.17 * length(genes), height = 2.5 + 0.17 * n_rows)
  print(dot_plot(d, s$title))
  for (so in levels(d$sort)) {
    print(dot_plot(d[sort == so], paste0(s$title, " - ", so, " sort")))
  }
  dev.off()
  step(nm, ": wrote marker_dotplot_", nm, ".pdf / .tsv")
}

step("done")
print(sessionInfo())
