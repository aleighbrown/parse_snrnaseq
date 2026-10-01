# soupx_vs_raw_celltype.R -------------------------------------------------
# How did SoupX ambient-RNA correction change the Azimuth cell typing?
# Compares the per-cell annotations from the uncorrected run (`celltype`)
# against the SoupX-corrected run (`celltype_soupx`), cell by cell.

library(data.table)   # fread only
library(tidyverse)
library(pheatmap)

base <- "/Users/annaleigh/cluster/neuroscience/direct_nanopore_rna/ariana_parse/spipe_jh_donor_annotations"
full_cellmeta <- fread(file.path(base, "celltype/cell_annotations.csv")) %>% as_tibble()
soup_meta     <- fread(file.path(base, "celltype_soupx/cell_annotations.csv")) %>% as_tibble()

# columns we compare; per-cell Azimuth calls + scores, and the cluster-level call
label_cols <- c("predicted.class", "predicted.subclass", "predicted.cluster",
                "cluster_azimuth_class", "cluster_azimuth_subclass")
score_cols <- c("predicted.class.score", "predicted.subclass.score",
                "predicted.cluster.score", "mapping.score")

# a cell with no Azimuth call has an empty string (or NA)
no_class <- function(x) is.na(x) | x == ""


# which cells are in each run ----------------------------------------------
cat("raw cells:", nrow(full_cellmeta), " soupx cells:", nrow(soup_meta), "\n")
cat("shared:", length(intersect(full_cellmeta$cell_id, soup_meta$cell_id)),
    " raw only:", length(setdiff(full_cellmeta$cell_id, soup_meta$cell_id)),
    " soupx only:", length(setdiff(soup_meta$cell_id, full_cellmeta$cell_id)), "\n")

# join on cell_id; only cells that were annotated in both runs are comparable
keep <- c("cell_id", "sample", "donor_id", "n_transcripts", "n_genes", "qc_pass",
          "seurat_clusters", intersect(c(label_cols, score_cols), names(full_cellmeta)))

cmp <- full_cellmeta %>%
    select(any_of(keep)) %>%
    inner_join(soup_meta %>% select(any_of(keep)),
               by = "cell_id", suffix = c(".raw", ".soup"))
cat("cells compared:", nrow(cmp), "\n")


# cells lost to SoupX: fail qc_pass and/or have no predicted.class ----------
cmp <- cmp %>%
    mutate(lost_qc    = qc_pass.raw %in% TRUE & !(qc_pass.soup %in% TRUE),
           lost_class = !no_class(predicted.class.raw) & no_class(predicted.class.soup))

cmp %>%
    summarise(passed_qc_raw      = sum(qc_pass.raw %in% TRUE),
              passed_qc_soupx    = sum(qc_pass.soup %in% TRUE),
              lost_qc            = sum(lost_qc),
              lost_class         = sum(lost_class),
              lost_both          = sum(lost_qc & lost_class),
              lost_either        = sum(lost_qc | lost_class),
              # the reverse direction, in case soupx rescued cells
              rescued_qc         = sum(!(qc_pass.raw %in% TRUE) & qc_pass.soup %in% TRUE),
              rescued_class      = sum(no_class(predicted.class.raw) & !no_class(predicted.class.soup))) %>%
    pivot_longer(everything(), names_to = "check", values_to = "n_cells") %>%
    print()

# raw cells that don't appear in the soupx table at all
full_cellmeta %>%
    filter(qc_pass %in% TRUE, !cell_id %in% soup_meta$cell_id) %>%
    nrow() %>%
    cat("raw qc-passing cells absent from soupx table:", ., "\n")

# cross-tab of the two outcomes, and which raw classes / samples the losses come from
cmp %>% count(lost_qc, lost_class)

cmp %>%
    filter(lost_qc | lost_class) %>%
    count(predicted.class.raw, sort = TRUE)

cmp %>%
    filter(lost_qc | lost_class) %>%
    count(sample.raw, sort = TRUE)

# if soupx also writes why a cell was excluded, say so
if ("exclude_reason" %in% names(soup_meta)) {
    soup_meta %>%
        filter(cell_id %in% (cmp %>% filter(lost_qc) %>% pull(cell_id))) %>%
        count(exclude_reason, sort = TRUE) %>%
        print()
}


# how many cells changed label at each level --------------------------------
changed <- map_dfr(intersect(label_cols, names(full_cellmeta)), function(col) {
    a <- cmp[[paste0(col, ".raw")]]
    b <- cmp[[paste0(col, ".soup")]]
    tibble(level = col,
           n_changed = sum(a != b, na.rm = TRUE),
           pct_changed = 100 * mean(a != b, na.rm = TRUE))
})
print(changed)


# confusion matrix: raw label (rows) vs soupx label (cols) -------------------
confusion <- function(col, normalise = TRUE) {
    tab <- table(raw = cmp[[paste0(col, ".raw")]], soupx = cmp[[paste0(col, ".soup")]])
    if (normalise) tab <- tab / rowSums(tab)   # fraction of each raw type going where
    tab
}

pheatmap(confusion("predicted.class"), cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "Azimuth class: raw (rows) -> SoupX (cols), row-normalised")
pheatmap(confusion("predicted.subclass"), cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Azimuth subclass: raw (rows) -> SoupX (cols), row-normalised")

# the actual moves, ranked
confusion("predicted.subclass", normalise = FALSE) %>%
    as.data.frame() %>%
    filter(raw != soupx, Freq > 0) %>%
    arrange(desc(Freq)) %>%
    head(20)


# composition shift ---------------------------------------------------------
comp <- bind_rows(
    cmp %>% count(predicted.subclass = predicted.subclass.raw)  %>% mutate(run = "raw"),
    cmp %>% count(predicted.subclass = predicted.subclass.soup) %>% mutate(run = "soupx")
) %>%
    group_by(run) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()

ggplot(comp, aes(x = predicted.subclass, y = frac, fill = run)) +
    geom_col(position = "dodge") +
    coord_flip() +
    theme_classic()


# did confidence improve? ---------------------------------------------------
scores <- map_dfr(intersect(score_cols, names(full_cellmeta)), function(col) {
    tibble(score = col,
           raw  = cmp[[paste0(col, ".raw")]],
           soup = cmp[[paste0(col, ".soup")]])
})

ggplot(scores, aes(x = raw, y = soup)) +
    geom_hex(bins = 60) +
    geom_abline(colour = "red") +
    facet_wrap(~score) +
    theme_classic()

scores %>%
    group_by(score) %>%
    summarise(median_raw   = median(raw, na.rm = TRUE),
              median_soup  = median(soup, na.rm = TRUE),
              median_delta = median(soup - raw, na.rm = TRUE))


# which cells flipped, and are they low-count / low-confidence? -------------
cmp <- cmp %>%
    mutate(flipped = predicted.subclass.raw != predicted.subclass.soup)

cmp %>%
    group_by(flipped) %>%
    summarise(n = n(),
              median_umis = median(n_transcripts.raw, na.rm = TRUE),
              median_score_raw = median(predicted.subclass.score.raw, na.rm = TRUE))

# per-sample flip rate (is the correction sample-specific, i.e. ambient-driven?)
cmp %>%
    group_by(sample.raw) %>%
    summarise(n = n(), pct_flipped = 100 * mean(flipped, na.rm = TRUE)) %>%
    arrange(desc(pct_flipped))

# seurat cluster membership: how well do raw clusters map to soupx clusters?
table(raw = cmp$seurat_clusters.raw, soupx = cmp$seurat_clusters.soup) %>%
    {. / rowSums(.)} %>%
    pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,
             main = "Seurat clusters: raw -> SoupX")


# cell type composition per donor and per sort, raw vs soupx -----------------
# Confident passing cells only: qc_pass and predicted.subclass.score >= min_score
# (0.8 matches `min_label_score` in config.yaml / soupx_rho.R). Each run is
# filtered on its own calls, so cells SoupX pushed below the line drop out of
# the soupx bar; the denominator is therefore the confident cells in that run.
min_score <- 0.8

ct_palette <- c(
    # excitatory - blues, upper -> deep layers
    "L2/3 IT"    = "#c6dbef",
    "L5 IT"      = "#9ecae1",
    "L6 IT"      = "#6baed6",
    "L6 IT Car3" = "#4292c6",
    "L5 ET"      = "#2171b5",
    "L5/6 NP"    = "#08519c",
    "L6 CT"      = "#08306b",
    "L6b"        = "#041f4a",
    # inhibitory - MGE reds, CGE oranges
    "Pvalb"      = "#99000d",
    "Sst"        = "#cb181d",
    "Sst Chodl"  = "#ef3b2c",
    "Vip"        = "#fd8d3c",
    "Lamp5"      = "#fdae6b",
    "Sncg"       = "#fdd0a2",
    # oligodendrocyte lineage - greens
    "Oligo"      = "#238b45",
    "OPC"        = "#74c476",
    # other glia
    "Astro"      = "#6a51a3",
    "Micro-PVM"  = "#8c6d31",
    # vascular - greys
    "Endo"       = "#737373",
    "VLMC"       = "#bdbdbd"
)

confident_raw <- cmp %>%
    filter(qc_pass.raw %in% TRUE, predicted.subclass.score.raw >= min_score) %>%
    transmute(cell_id, donor_id = donor_id.raw, sort = sample.raw,
              subclass = predicted.subclass.raw, run = "raw")

confident_soup <- cmp %>%
    filter(qc_pass.soup %in% TRUE, predicted.subclass.score.soup >= min_score) %>%
    transmute(cell_id, donor_id = donor_id.raw, sort = sample.raw,
              subclass = predicted.subclass.soup, run = "soupx")

confident <- bind_rows(confident_raw, confident_soup) %>%
    mutate(run = factor(run, c("raw", "soupx")))

confident %>% count(run) %>% print()

# any Azimuth subclass missing from the palette would plot as grey NA
setdiff(unique(confident$subclass), names(ct_palette)) %>% print()

confident <- confident %>%
    filter(subclass %in% names(ct_palette)) %>%
    mutate(subclass = factor(subclass, names(ct_palette)))

plot_composition <- function(by, ncol = 4) {
    confident %>%
        count(.data[[by]], run, subclass) %>%
        group_by(.data[[by]], run) %>%
        mutate(frac = n / sum(n)) %>%
        ungroup() %>%
        ggplot(aes(x = run, y = frac, fill = subclass)) +
        geom_col(width = 0.85) +
        facet_wrap(vars(.data[[by]]), ncol = ncol) +
        scale_fill_manual(values = ct_palette, drop = FALSE) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(x = NULL, y = "fraction of confident passing cells", fill = "Azimuth subclass",
             title = paste("Cell type composition by", by, ": raw vs SoupX")) +
        theme_classic() +
        theme(strip.background = element_blank())
}

plot_composition("donor_id")
plot_composition("sort", ncol = 3)

# the same, as the change in fraction (soupx - raw) so small shifts are visible
plot_composition_delta <- function(by, ncol = 4) {
    confident %>%
        count(.data[[by]], run, subclass, .drop = FALSE) %>%
        group_by(.data[[by]], run) %>%
        mutate(frac = n / sum(n)) %>%
        ungroup() %>%
        select(-n) %>%
        pivot_wider(names_from = run, values_from = frac) %>%
        mutate(delta = soupx - raw) %>%
        ggplot(aes(x = subclass, y = delta, fill = subclass)) +
        geom_col() +
        geom_hline(yintercept = 0) +
        facet_wrap(vars(.data[[by]]), ncol = ncol) +
        scale_fill_manual(values = ct_palette, drop = FALSE, guide = "none") +
        coord_flip() +
        labs(x = NULL, y = "change in fraction (soupx - raw)",
             title = paste("SoupX shift in composition by", by)) +
        theme_classic() +
        theme(strip.background = element_blank())
}

plot_composition_delta("donor_id")
plot_composition_delta("sort", ncol = 3)


# SoupX run only: confident cell counts and NeuN+ (neuron) yield -------------
# "Expected NeuN positive" = neuronal subclasses (excitatory + inhibitory, the
# first 14 entries of ct_palette). Confident passing cells as defined above.
neuron_subclasses <- names(ct_palette)[1:14]
sort_levels <- c("positive", "TDP_negative", "K_negative")
sort_cols   <- c(positive = "blue", TDP_negative = "green", K_negative = "orange")

soup_conf <- confident %>%
    filter(run == "soupx") %>%
    mutate(neun_pos = subclass %in% neuron_subclasses,
           sort = factor(sort, intersect(sort_levels, unique(sort))))

summarise_by <- function(...) {
    soup_conf %>%
        group_by(...) %>%
        summarise(n_confident   = n(),
                  n_neun_pos    = sum(neun_pos),
                  frac_neun_pos = mean(neun_pos),
                  .groups = "drop")
}
by_donor      <- summarise_by(donor_id)
by_sort       <- summarise_by(sort)
by_donor_sort <- summarise_by(donor_id, sort)

print(by_donor)
print(by_sort)

bar_theme <- theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# raw number of confident cells
ggplot(by_donor, aes(x = reorder(donor_id, -n_confident), y = n_confident)) +
    geom_col() + labs(x = NULL, y = "confident passing cells (SoupX)") + bar_theme
ggplot(by_sort, aes(x = sort, y = n_confident, fill = sort)) +
    geom_col() + scale_fill_manual(values = sort_cols, guide = "none") +
    labs(x = NULL, y = "confident passing cells (SoupX)") + bar_theme
ggplot(by_donor_sort, aes(x = donor_id, y = n_confident, fill = sort)) +
    geom_col(position = "dodge") + scale_fill_manual(values = sort_cols) +
    labs(x = NULL, y = "confident passing cells (SoupX)") + bar_theme

# raw number of NeuN+ (neuronal) cells
ggplot(by_donor, aes(x = reorder(donor_id, -n_neun_pos), y = n_neun_pos)) +
    geom_col() + labs(x = NULL, y = "NeuN+ (neuronal) confident cells (SoupX)") + bar_theme
ggplot(by_sort, aes(x = sort, y = n_neun_pos, fill = sort)) +
    geom_col() + scale_fill_manual(values = sort_cols, guide = "none") +
    labs(x = NULL, y = "NeuN+ (neuronal) confident cells (SoupX)") + bar_theme
ggplot(by_donor_sort, aes(x = donor_id, y = n_neun_pos, fill = sort)) +
    geom_col(position = "dodge") + scale_fill_manual(values = sort_cols) +
    labs(x = NULL, y = "NeuN+ (neuronal) confident cells (SoupX)") + bar_theme

# fraction of confident cells that are NeuN+
ggplot(by_donor, aes(x = reorder(donor_id, -frac_neun_pos), y = frac_neun_pos)) +
    geom_col() + labs(x = NULL, y = "fraction NeuN+ of confident cells (SoupX)") + bar_theme
ggplot(by_sort, aes(x = sort, y = frac_neun_pos, fill = sort)) +
    geom_col() + scale_fill_manual(values = sort_cols, guide = "none") +
    labs(x = NULL, y = "fraction NeuN+ of confident cells (SoupX)") + bar_theme
ggplot(by_donor_sort, aes(x = donor_id, y = frac_neun_pos, fill = sort)) +
    geom_col(position = "dodge") + scale_fill_manual(values = sort_cols) +
    labs(x = NULL, y = "fraction NeuN+ of confident cells (SoupX)") + bar_theme


# K mislocalisation score vs fraction NeuN+, by sort (SoupX, no controls) ----
# Donor-level K score is carried on every cell row of the annotation table.
# The raw header isn't snake_case (sorts_checking.R needs janitor::clean_names()
# to reach `k_mislocalisation_score`), so match on the cleaned names.
k_col <- names(soup_meta)[janitor::make_clean_names(names(soup_meta)) == "k_mislocalisation_score"]
if (length(k_col) != 1) {
    print(grep("k|mislocal", names(soup_meta), ignore.case = TRUE, value = TRUE))
    stop("couldn't find a single K mislocalisation score column in soup_meta; set k_col by hand")
}

donor_info <- soup_meta %>%
    group_by(donor_id) %>%
    summarise(k_score   = as.numeric(first(.data[[k_col]])),
              path_diag = first(path_diag),
              .groups = "drop")

k_vs_neun <- by_donor_sort %>%
    left_join(donor_info, by = "donor_id") %>%
    filter(path_diag != "control", !is.na(k_score))

# Spearman per sort, printed and used as the facet label
k_cor <- k_vs_neun %>%
    group_by(sort) %>%
    summarise(n_donors = n(),
              rho = cor(k_score, frac_neun_pos, method = "spearman"),
              p   = if (n() > 2) cor.test(k_score, frac_neun_pos, method = "spearman", exact = FALSE)$p.value else NA_real_,
              .groups = "drop") %>%
    mutate(label = sprintf("rho = %.2f, p = %.2g, n = %d", rho, p, n_donors))
print(k_cor)

ggplot(k_vs_neun, aes(x = k_score, y = frac_neun_pos)) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey40", linewidth = 0.5) +
    geom_point(aes(colour = sort), size = 2) +
    geom_text(data = k_cor, aes(label = label), x = -Inf, y = Inf,
              hjust = -0.05, vjust = 1.5, size = 3, inherit.aes = FALSE) +
    facet_wrap(~sort) +
    scale_colour_manual(values = sort_cols, guide = "none") +
    labs(x = "K mislocalisation score", y = "fraction NeuN+ of confident cells (SoupX)",
         title = "K score vs NeuN+ fraction by sort (controls excluded)") +
    theme_classic()


# which patients to keep? NeuN+ cells in the K_negative sort ------------------
# Candidate rule: keep patients with >= min_k_neun confident NeuN+ cells in
# their K_negative sort (SoupX run, controls excluded). Change the threshold
# and re-run to see who falls out.
min_k_neun <- 300

# one row per patient x sort, zero-filled so a missing sort shows up as 0
neun_per_donor_sort <- by_donor_sort %>%
    complete(donor_id, sort, fill = list(n_confident = 0, n_neun_pos = 0)) %>%
    left_join(donor_info, by = "donor_id") %>%
    filter(path_diag != "control")

k_neun <- neun_per_donor_sort %>%
    filter(sort == "K_negative") %>%
    mutate(keep = n_neun_pos >= min_k_neun)

k_neun %>% count(keep) %>% print()

keep_donors <- k_neun %>% filter(keep) %>% pull(donor_id)
drop_donors <- k_neun %>% filter(!keep) %>% pull(donor_id)
cat("keep:", length(keep_donors), " drop:", length(drop_donors), "\n")
cat("dropped:", paste(drop_donors, collapse = ", "), "\n")

# 1. every patient's NeuN+ cell count in each sort, K_negative with the cutoff;
#    patients are ordered by their K_negative count and coloured by keep/drop
donor_order <- k_neun %>% arrange(n_neun_pos) %>% pull(donor_id)

neun_per_donor_sort %>%
    left_join(k_neun %>% select(donor_id, keep), by = "donor_id") %>%
    mutate(donor_id = factor(donor_id, donor_order)) %>%
    ggplot(aes(x = n_neun_pos, y = donor_id, fill = keep)) +
    geom_col() +
    geom_vline(xintercept = min_k_neun, linetype = "dashed") +
    facet_wrap(~sort, nrow = 1) +
    scale_fill_manual(values = c("TRUE" = "grey30", "FALSE" = "firebrick"),
                      labels = c("TRUE" = "keep", "FALSE" = "drop"), name = NULL) +
    labs(x = "confident NeuN+ cells (SoupX)", y = NULL,
         title = paste0("NeuN+ cells per patient and sort; cutoff = ", min_k_neun, " in K_negative")) +
    theme_classic()

# 2. how many patients survive at each possible cutoff
thresholds <- tibble(cutoff = seq(0, max(k_neun$n_neun_pos), length.out = 200)) %>%
    mutate(n_kept = map_int(cutoff, ~ sum(k_neun$n_neun_pos >= .x)))

ggplot(thresholds, aes(x = cutoff, y = n_kept)) +
    geom_step() +
    geom_vline(xintercept = min_k_neun, linetype = "dashed") +
    labs(x = "min NeuN+ cells in K_negative", y = "patients kept",
         title = "Patients retained vs K_negative NeuN+ cutoff") +
    theme_classic()

# 3. does dropping patients bias the K score story? kept vs dropped K scores
ggplot(k_neun, aes(x = k_score, y = n_neun_pos, colour = keep)) +
    geom_point(size = 2) +
    geom_hline(yintercept = min_k_neun, linetype = "dashed") +
    scale_colour_manual(values = c("TRUE" = "grey30", "FALSE" = "firebrick")) +
    labs(x = "K mislocalisation score", y = "confident NeuN+ cells in K_negative",
         title = "Who the cutoff removes, by K score") +
    theme_classic()


# What population of cells is relatively enriched in the K_negative sorts? ----
# SoupX run, confident cells, controls excluded. Composition is computed per
# patient x sort, then K_negative is compared with another sort *within the
# same patient*, so donor differences don't masquerade as sort differences.
# Patients need >= min_cells_per_sort confident cells in a sort to count there,
# otherwise the fractions are too noisy to compare.
min_cells_per_sort <- 50

sub_frac <- soup_conf %>%
    left_join(donor_info, by = "donor_id") %>%
    filter(path_diag != "control") %>%
    mutate(subclass = as.character(subclass), sort = as.character(sort)) %>%
    count(donor_id, sort, subclass) %>%
    complete(nesting(donor_id, sort), subclass = names(ct_palette), fill = list(n = 0)) %>%
    group_by(donor_id, sort) %>%
    mutate(n_total = sum(n), frac = n / n_total) %>%
    ungroup() %>%
    filter(n_total >= min_cells_per_sort) %>%
    mutate(subclass = factor(subclass, names(ct_palette)))

# mean per-patient composition in each sort
sub_frac %>%
    group_by(sort, subclass) %>%
    summarise(mean_frac = mean(frac), n_donors = n(), .groups = "drop") %>%
    ggplot(aes(x = factor(sort, sort_levels), y = mean_frac, fill = subclass)) +
    geom_col(width = 0.85) +
    scale_fill_manual(values = ct_palette, drop = FALSE) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "mean per-patient fraction of confident cells", fill = "Azimuth subclass",
         title = "Average composition by sort (SoupX, controls excluded)") +
    theme_classic()

# paired K_negative vs reference-sort enrichment per subclass
k_enrichment <- function(reference, pseudo = 0.001) {
    paired <- sub_frac %>%
        filter(sort %in% c("K_negative", reference)) %>%
        select(donor_id, subclass, sort, frac) %>%
        pivot_wider(names_from = sort, values_from = frac) %>%
        drop_na() %>%                      # only patients with both sorts
        mutate(log2_enrich = log2((K_negative + pseudo) / (.data[[reference]] + pseudo)),
               diff_frac   = K_negative - .data[[reference]])

    stats <- paired %>%
        group_by(subclass) %>%
        summarise(n_donors     = n(),
                  mean_k       = mean(K_negative),
                  mean_ref     = mean(.data[[reference]]),
                  median_log2  = median(log2_enrich),
                  p            = if (n() >= 3) wilcox.test(K_negative, .data[[reference]], paired = TRUE, exact = FALSE)$p.value else NA_real_,
                  .groups = "drop") %>%
        mutate(fdr = p.adjust(p, method = "BH"),
               reference = reference) %>%
        arrange(desc(median_log2))

    list(paired = paired, stats = stats)
}

plot_k_enrichment <- function(res) {
    ref <- unique(res$stats$reference)
    ggplot(res$paired, aes(x = subclass, y = log2_enrich, fill = subclass)) +
        geom_hline(yintercept = 0) +
        geom_boxplot(outlier.shape = NA, alpha = 0.8) +
        geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
        scale_fill_manual(values = ct_palette, drop = FALSE, guide = "none") +
        coord_flip() +
        labs(x = NULL, y = paste0("log2 (K_negative / ", ref, ") fraction, per patient"),
             title = paste0("Subclass enrichment in K_negative vs ", ref)) +
        theme_classic()
}

k_vs_positive <- k_enrichment("positive")
k_vs_tdp      <- k_enrichment("TDP_negative")

print(k_vs_positive$stats, n = Inf)
print(k_vs_tdp$stats, n = Inf)

plot_k_enrichment(k_vs_positive)
plot_k_enrichment(k_vs_tdp)

# same question in absolute terms: change in fraction of the sort, per subclass
bind_rows(k_vs_positive$paired %>% mutate(reference = "positive") %>% rename(ref_frac = positive),
          k_vs_tdp$paired      %>% mutate(reference = "TDP_negative") %>% rename(ref_frac = TDP_negative)) %>%
    ggplot(aes(x = subclass, y = diff_frac, fill = subclass)) +
    geom_hline(yintercept = 0) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    facet_wrap(~reference) +
    scale_fill_manual(values = ct_palette, drop = FALSE, guide = "none") +
    coord_flip() +
    labs(x = NULL, y = "K_negative fraction - reference fraction, per patient",
         title = "Which subclasses gain share in K_negative?") +
    theme_classic()


# Relative likelihood of each cell type, per sort ----------------------------
# SoupX run, confident cells, controls INCLUDED, all patients pooled.
# relative likelihood = (fraction of this sort's cells that are subclass X) /
#                       (fraction of ALL confident cells that are subclass X)
# > 1 means a cell from this sort is more likely to be X than a cell drawn at
# random from the whole dataset; < 1 less likely. Plotted on a log2 axis so
# 2x enriched and 2x depleted are the same distance from the baseline.
baseline <- soup_conf %>%
    filter(!is.na(sort)) %>%
    count(subclass) %>%
    mutate(baseline_frac = n / sum(n)) %>%
    select(subclass, baseline_frac)

rel_likelihood <- soup_conf %>%
    filter(!is.na(sort)) %>%
    count(sort, subclass, .drop = FALSE) %>%
    group_by(sort) %>%
    mutate(n_sort = sum(n), sort_frac = n / n_sort) %>%
    ungroup() %>%
    left_join(baseline, by = "subclass") %>%
    mutate(rel_likelihood = sort_frac / baseline_frac,
           log2_rel = log2(rel_likelihood))

rel_likelihood %>%
    select(sort, subclass, n, sort_frac, baseline_frac, rel_likelihood) %>%
    arrange(sort, desc(rel_likelihood)) %>%
    print(n = Inf)

ggplot(rel_likelihood, aes(x = subclass, y = log2_rel, fill = sort)) +
    geom_hline(yintercept = 0) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = sort_cols) +
    coord_flip() +
    labs(x = NULL, y = "log2 relative likelihood vs all-cell baseline",
         fill = "sort",
         title = "How much more (or less) likely is each cell type in each sort?") +
    theme_classic()

# the same as a heatmap, easier to read across sorts
ggplot(rel_likelihood, aes(x = sort, y = subclass, fill = log2_rel)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1f", rel_likelihood)), size = 3) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0,
                         limits = c(-3, 3), oob = scales::squish) +
    labs(x = NULL, y = NULL, fill = "log2\nrel. likelihood",
         title = "Relative likelihood of each subclass by sort (numbers = fold vs baseline)") +
    theme_classic()


# L5 ET per patient and sort ---------------------------------------------------
# L5 ET is rare overall but more likely than baseline in K_negative. SoupX run,
# confident cells, controls included (shown as their own facet). Fractions are
# of all confident cells in that patient x sort; the neuron-only fraction
# removes the effect of glia dilution when sorts differ in neuron content.
target <- "L5 ET"

l5et <- soup_conf %>%
    filter(!is.na(sort)) %>%
    left_join(donor_info, by = "donor_id") %>%
    group_by(donor_id, path_diag, sort) %>%
    summarise(n_total    = n(),
              n_neurons  = sum(neun_pos),
              n_target   = sum(subclass == target),
              frac_all     = n_target / n_total,
              frac_neurons = n_target / n_neurons,
              .groups = "drop")

l5et %>% arrange(desc(n_target)) %>% print(n = Inf)

# total L5 ET found in each sort
l5et %>% group_by(sort) %>%
    summarise(n_target = sum(n_target), n_total = sum(n_total),
              n_donors_with_any = sum(n_target > 0), n_donors = n(),
              frac_all = n_target / n_total) %>%
    print()

l5et_bars <- function(y, ylab) {
    ggplot(l5et, aes(x = donor_id, y = .data[[y]], fill = sort)) +
        geom_col(position = position_dodge2(preserve = "single")) +
        facet_grid(~path_diag, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = sort_cols) +
        labs(x = NULL, y = ylab, title = paste(target, "per patient and sort (SoupX)")) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              strip.background = element_blank())
}

l5et_bars("n_target",     paste(target, "confident cells (raw number)"))
l5et_bars("frac_all",     paste("fraction of confident cells that are", target))
l5et_bars("frac_neurons", paste("fraction of NeuN+ cells that are", target))

# per-patient trajectory across sorts: is the K_negative enrichment consistent?
l5et %>%
    ggplot(aes(x = factor(sort, sort_levels), y = frac_all, group = donor_id, colour = path_diag)) +
    geom_line(alpha = 0.6) +
    geom_point(aes(size = n_target)) +
    labs(x = NULL, y = paste("fraction of confident cells that are", target),
         size = paste("n", target), colour = "diagnosis",
         title = paste(target, "fraction across sorts, one line per patient")) +
    theme_classic()

# is it just depth? raw L5 ET count against the number of confident cells in the sort
ggplot(l5et, aes(x = n_total, y = n_target, colour = sort)) +
    geom_point(size = 2) +
    scale_colour_manual(values = sort_cols) +
    labs(x = "confident cells in patient x sort", y = paste(target, "cells"),
         title = paste(target, "count vs sort size")) +
    theme_classic()

# K_negative vs positive within patient, on the raw L5 ET fraction
l5et %>%
    filter(sort %in% c("K_negative", "positive")) %>%
    select(donor_id, path_diag, sort, frac_all) %>%
    pivot_wider(names_from = sort, values_from = frac_all) %>%
    drop_na() %>%
    ggplot(aes(x = positive, y = K_negative, colour = path_diag)) +
    geom_abline(colour = "grey60") +
    geom_point(size = 2) +
    labs(x = paste(target, "fraction in positive sort"), y = paste(target, "fraction in K_negative sort"),
         colour = "diagnosis", title = paste(target, ": K_negative vs positive, per patient")) +
    theme_classic()


# UMAP coloured by Azimuth mapping score, singlets only ----------------------
# sort (sample) across the columns, Azimuth predicted.class down the rows.
plot_umap_mapping <- function(meta, run_label) {
    meta %>%
        filter(doublet_class == "singlet", !is.na(umap_harmony_1), !is.na(mapping.score)) %>%
        mutate(predicted.class = if_else(no_class(predicted.class), "no call", predicted.class)) %>%
        arrange(mapping.score) %>%          # best-mapped cells drawn last, on top
        ggplot(aes(x = umap_harmony_1, y = umap_harmony_2, colour = mapping.score)) +
        geom_point(size = 0.1, stroke = 0) +
        facet_grid(rows = vars(predicted.class), cols = vars(sample)) +
        scale_colour_viridis_c(limits = c(0, 1)) +
        coord_fixed() +
        labs(colour = "mapping\nscore",
             title = paste("Azimuth mapping score,", run_label, "(singlets)")) +
        theme_classic() +
        theme(strip.background = element_blank())
}

plot_umap_mapping(soup_meta, "SoupX")
plot_umap_mapping(full_cellmeta, "raw")


# Mapping score by sort and predicted.class, singlets only --------------------
# Two views of the same question: the score itself (boxplot, Wilcoxon), and the
# fraction of cells scoring >= min_score (0.8; proportion test). With this many
# cells almost everything is "significant", so read the effect sizes
# (median difference, difference in fraction) alongside the p-values.
sort_pairs <- tibble(a = c("K_negative", "K_negative", "TDP_negative"),
                     b = c("positive",   "TDP_negative", "positive"))

compare_mapping <- function(meta, run_label) {
    cells <- meta %>%
        filter(doublet_class == "singlet", !is.na(mapping.score), sample %in% sort_levels) %>%
        transmute(sort = factor(sample, sort_levels),
                  class = if_else(no_class(predicted.class), "no call", predicted.class),
                  mapping.score,
                  high = mapping.score >= min_score)

    # -- boxplot of the score -------------------------------------------------
    p_box <- ggplot(cells, aes(x = class, y = mapping.score, fill = sort)) +
        geom_hline(yintercept = min_score, linetype = "dashed") +
        geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), width = 0.7) +
        scale_fill_manual(values = sort_cols) +
        labs(x = NULL, y = "Azimuth mapping score",
             title = paste("Mapping score by sort and class,", run_label, "(singlets)")) +
        theme_classic()

    # -- fraction of cells >= min_score (approximate 95% CI) -------------------
    frac <- cells %>%
        group_by(class, sort) %>%
        summarise(n = n(), n_high = sum(high), frac_high = mean(high), .groups = "drop") %>%
        mutate(se = sqrt(frac_high * (1 - frac_high) / n),
               lo = pmax(0, frac_high - 1.96 * se),
               hi = pmin(1, frac_high + 1.96 * se))

    p_frac <- ggplot(frac, aes(x = class, y = frac_high, fill = sort)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7) +
        geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.8), width = 0.25) +
        scale_fill_manual(values = sort_cols) +
        labs(x = NULL, y = paste("fraction of cells with mapping score >=", min_score),
             title = paste("High-confidence mapping by sort and class,", run_label, "(singlets)")) +
        theme_classic()

    # -- stats: pairwise between sorts, within each class ---------------------
    pairwise <- crossing(class = unique(cells$class), sort_pairs) %>%
        mutate(res = pmap(list(class, a, b), function(cl, a, b) {
            x <- cells %>% filter(class == cl, sort == a)
            y <- cells %>% filter(class == cl, sort == b)
            if (nrow(x) < 10 || nrow(y) < 10) return(NULL)
            tibble(n_a = nrow(x), n_b = nrow(y),
                   median_a = median(x$mapping.score), median_b = median(y$mapping.score),
                   median_diff = median_a - median_b,
                   p_score = wilcox.test(x$mapping.score, y$mapping.score)$p.value,
                   frac_high_a = mean(x$high), frac_high_b = mean(y$high),
                   frac_diff = frac_high_a - frac_high_b,
                   p_frac = suppressWarnings(
                       prop.test(c(sum(x$high), sum(y$high)), c(nrow(x), nrow(y)))$p.value))
        })) %>%
        unnest(res) %>%
        mutate(fdr_score = p.adjust(p_score, method = "BH"),
               fdr_frac  = p.adjust(p_frac,  method = "BH")) %>%
        arrange(class, a, b)

    # -- stats: any difference across the three sorts, within each class -------
    overall <- cells %>%
        group_by(class) %>%
        filter(n_distinct(sort) == 3) %>%
        summarise(p_score = kruskal.test(mapping.score ~ sort)$p.value,
                  p_frac  = suppressWarnings(chisq.test(table(sort, high))$p.value),
                  .groups = "drop") %>%
        mutate(fdr_score = p.adjust(p_score, method = "BH"),
               fdr_frac  = p.adjust(p_frac,  method = "BH"))

    list(box = p_box, frac_plot = p_frac, frac = frac, pairwise = pairwise, overall = overall)
}

map_soup <- compare_mapping(soup_meta, "SoupX")
map_raw  <- compare_mapping(full_cellmeta, "raw")

map_soup$box
map_soup$frac_plot
print(map_soup$frac, n = Inf)
print(map_soup$overall, n = Inf)
print(map_soup$pairwise, n = Inf)

map_raw$box
map_raw$frac_plot
print(map_raw$overall, n = Inf)
print(map_raw$pairwise, n = Inf)


# What drives mapping score: N genes, N UMIs, or both? -------------------------
# Singlets only. n_genes and n_transcripts are strongly correlated, so the
# question is whether either carries information the other doesn't. Answered
# three ways: depth by sort/class (same boxplots + stats as above), score
# against depth, and models comparing genes alone, UMIs alone and both.
depth_cells <- function(meta) {
    meta %>%
        filter(doublet_class == "singlet", !is.na(mapping.score), sample %in% sort_levels,
               n_genes > 0, n_transcripts > 0) %>%
        transmute(sort = factor(sample, sort_levels),
                  class = if_else(no_class(predicted.class), "no call", predicted.class),
                  mapping.score,
                  high = mapping.score >= min_score,
                  n_genes, n_transcripts,
                  log_genes = log10(n_genes),
                  log_umis  = log10(n_transcripts))
}

# -- 1. depth itself by sort and class, with pairwise stats ---------------------
compare_depth <- function(cells, var, run_label) {
    p <- ggplot(cells, aes(x = class, y = .data[[var]], fill = sort)) +
        geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), width = 0.7) +
        scale_y_log10() +
        scale_fill_manual(values = sort_cols) +
        labs(x = NULL, y = var, title = paste(var, "by sort and class,", run_label, "(singlets)")) +
        theme_classic()

    stats <- crossing(class = unique(cells$class), sort_pairs) %>%
        mutate(res = pmap(list(class, a, b), function(cl, a, b) {
            x <- cells %>% filter(class == cl, sort == a) %>% pull(!!sym(var))
            y <- cells %>% filter(class == cl, sort == b) %>% pull(!!sym(var))
            if (length(x) < 10 || length(y) < 10) return(NULL)
            tibble(n_a = length(x), n_b = length(y),
                   median_a = median(x), median_b = median(y),
                   fold_diff = median_a / median_b,
                   p = wilcox.test(x, y)$p.value)
        })) %>%
        unnest(res) %>%
        mutate(fdr = p.adjust(p, method = "BH"), variable = var) %>%
        arrange(class, a, b)
    list(plot = p, stats = stats)
}

# -- 2. score against depth: binned trend, per sort ------------------------------
depth_trend <- function(cells, run_label) {
    cells %>%
        pivot_longer(c(log_genes, log_umis), names_to = "metric", values_to = "depth") %>%
        group_by(metric) %>%
        mutate(depth_bin = ntile(depth, 20)) %>%
        group_by(metric, sort, depth_bin) %>%
        summarise(depth = median(depth), mean_score = mean(mapping.score),
                  frac_high = mean(high), .groups = "drop") %>%
        pivot_longer(c(mean_score, frac_high), names_to = "measure", values_to = "value") %>%
        ggplot(aes(x = depth, y = value, colour = sort)) +
        geom_line() + geom_point(size = 1) +
        facet_grid(measure ~ metric, scales = "free") +
        scale_colour_manual(values = sort_cols) +
        labs(x = "log10 depth (bin median)", y = NULL,
             title = paste("Mapping score vs depth in 20 quantile bins,", run_label)) +
        theme_classic() +
        theme(strip.background = element_blank())
}

depth_hex <- function(cells, run_label) {
    cells %>%
        pivot_longer(c(log_genes, log_umis), names_to = "metric", values_to = "depth") %>%
        ggplot(aes(x = depth, y = mapping.score)) +
        geom_hex(bins = 60) +
        facet_grid(sort ~ metric, scales = "free_x") +
        scale_fill_viridis_c(trans = "log10") +
        labs(x = "log10 depth", title = paste("Mapping score vs depth,", run_label)) +
        theme_classic() +
        theme(strip.background = element_blank())
}

# -- 3. which one matters? models --------------------------------------------------
# R2 of each model, and the unique R2 each variable adds on top of the other
# (full model minus the model without it). Unique R2 near zero for one of them
# means the other already explains what it explains.
r2 <- function(formula, data) summary(lm(formula, data = data))$r.squared

driver_models <- function(d) {
    r2_genes <- r2(mapping.score ~ log_genes, d)
    r2_umis  <- r2(mapping.score ~ log_umis, d)
    r2_both  <- r2(mapping.score ~ log_genes + log_umis, d)
    tibble(n = nrow(d),
           cor_genes_umis = cor(d$log_genes, d$log_umis),
           spearman_genes = cor(d$mapping.score, d$log_genes, method = "spearman"),
           spearman_umis  = cor(d$mapping.score, d$log_umis,  method = "spearman"),
           r2_genes_only  = r2_genes,
           r2_umis_only   = r2_umis,
           r2_both        = r2_both,
           unique_r2_genes = r2_both - r2_umis,    # what genes add beyond UMIs
           unique_r2_umis  = r2_both - r2_genes)   # what UMIs add beyond genes
}

# standardised effects in the joint model, adjusting for class and sort
joint_effects <- function(cells) {
    fit <- lm(scale(mapping.score) ~ scale(log_genes) + scale(log_umis) + class + sort, data = cells)
    summary(fit)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("term") %>%
        filter(grepl("^scale", term)) %>%
        as_tibble()
}

run_drivers <- function(meta, run_label) {
    cells <- depth_cells(meta)

    list(cells = cells,
         genes  = compare_depth(cells, "n_genes", run_label),
         umis   = compare_depth(cells, "n_transcripts", run_label),
         trend  = depth_trend(cells, run_label),
         hex    = depth_hex(cells, run_label),
         overall    = driver_models(cells),
         by_class   = cells %>% group_by(class) %>% filter(n() >= 100) %>%
                          group_modify(~ driver_models(.x)) %>% ungroup(),
         by_sort    = cells %>% group_by(sort) %>% group_modify(~ driver_models(.x)) %>% ungroup(),
         joint      = joint_effects(cells))
}

drv_soup <- run_drivers(soup_meta, "SoupX")
drv_raw  <- run_drivers(full_cellmeta, "raw")

# depth by sort and class
drv_soup$genes$plot;  print(drv_soup$genes$stats, n = Inf)
drv_soup$umis$plot;   print(drv_soup$umis$stats, n = Inf)

# score vs depth
drv_soup$hex
drv_soup$trend

# which one drives it
print(drv_soup$overall)
print(drv_soup$by_sort)
print(drv_soup$by_class, n = Inf)
print(drv_soup$joint)

# raw run for comparison (SoupX removes ambient counts, so UMIs shift)
drv_raw$trend
print(drv_raw$overall)
print(drv_raw$joint)
