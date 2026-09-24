# -----------------------------------------------------------------------------
# Ambient RNA (SoupX), in three steps:
#
#   1. soup_profile   3 jobs    one per FANS sort: sum the empty barcodes of
#                               splitpipe/combined/<sort>/DGE_unfiltered
#   2. soupx_rho      3 jobs    one per sort: rho by SoupX's marker method
#                               (estimateNonExpressingCells +
#                               calculateContaminationFraction) on that sort's
#                               glia, pooled over all sublibraries
#   3. soupx_adjust   16 jobs   one per sublibrary: adjustCounts per sort
#   4. celltype_cluster_soupx  ONE big job  celltype_cluster.R again, on the
#                               corrected counts, with the per-sort doublet
#                               calls and a minimum size of min_floor_k x the
#                               background floor -> <results_dir>/celltype_soupx/
#
#   source submit.sh soup_profiles
#   source submit.sh soupx_rhos        # step 2 only: check rho before adjusting
#   source submit.sh soupx             # step 3 (runs step 2 first if needed)
#   source submit.sh celltype_soupx    # step 4 (runs 2-3 first if needed)
#   (soupx_rhos_test / soupx_test: the same on the celltype_test cells)
#
# Steps 2-3 sit between celltype_prep and clustering: they read prep/counts
# and correct every prep cell, so the cluster stage can filter the corrected
# counts exactly as it filters the raw ones. Cluster labels come from a
# FROZEN copy of a finished run's cell_annotations.csv (soupx.labels_file,
# in the run's directory), not from the DAG -- otherwise clustering would
# need SoupX and SoupX would need clustering. Doublets for the rho estimate
# come from doublets_by_sort.smk.
#
# Why per sort: all donors were pooled within each sort before loading, so each
# sort had one suspension and one soup. Ambient RNA is reverse-transcribed in
# the round-1 wells, so it carries the sort's round-1 barcode, which is how
# split-pipe already split it into combined/<sort>/. Sublibraries are split
# after pooling and each holds all three soups.
#
# soup_profile depends only on split-pipe output, not on the celltype stages.
#
# Outputs of step 1, under <results_dir>/soup/<sort>/:
#   soup_profile.tsv.gz         pooled per gene: counts, est (SoupX's format)
#   soup_by_sublibrary.tsv.gz   the same per sublibrary
#   barcode_summary.tsv, sublibrary_consistency.tsv, top_soup_genes.tsv
#   soup_report.pdf             barcode ranks with the empty range shaded
# Steps 2-3, under <results_dir>/<run>/soupx/:
#   rho/<sort>.tsv, rho/<sort>_breakdown.tsv, rho/<sort>_sensitivity.tsv,
#   rho/<sort>_report.pdf
#   counts/<sublib>.rds         corrected counts, same layout as prep/counts
#   cells/<sublib>.tsv          per cell: sort, rho, UMIs before / after
#
# Uses CT, CT_DGE, ct_fmt, _CT_SHELL and _CT_IMAGE from celltype.smk, so it is
# included after it.
# -----------------------------------------------------------------------------

SX = CT["soupx"]
SOUP_DIR = RESULTS + "/soup"
SPIPE_COMBINED = os.path.join(RESULTS, "splitpipe", "combined")

# The sorts are the split-pipe samples: first column of the Round-1 sample list.
with open(config["splitpipe"]["sample_list"]) as fh:
    SOUP_SORTS = [ln.split()[0] for ln in fh if ln.strip() and not ln.startswith("#")]

wildcard_constraints:
    soup_sort="|".join(SOUP_SORTS),

localrules: soup_profiles, soupx_rhos, soupx_rhos_test, soupx, soupx_test, celltype_soupx


def soup_umi_range(sort):
    """[min, max] UMIs of an empty barcode for this sort: the sort's own entry
    in soupx.empty_umi_range, else its `default` (a bare [min, max] applies
    to every sort)."""
    r = SX["empty_umi_range"]
    return r.get(sort, r["default"]) if isinstance(r, dict) else r


def soup_profile_args(wc):
    umi_min, umi_max = soup_umi_range(wc.soup_sort)
    return ct_fmt({
        "dge_dir": os.path.join(SPIPE_COMBINED, wc.soup_sort, "DGE_unfiltered"),
        "called_cells": CT_DGE + "/cell_metadata.csv",
        "ref_genes": CT_DGE + "/all_genes.csv",
        "sublibraries": config["sublibraries"],
        "sort": wc.soup_sort,
        "out_dir": SOUP_DIR + "/" + wc.soup_sort,
        "umi_min": umi_min,
        "umi_max": umi_max,
        "n_top": SX["n_top_genes"],
    })


rule soup_profile:
    input:
        # ancient(): split-pipe output; a touched timestamp must never pull
        # split-pipe back into the DAG from here.
        mtx=ancient(SPIPE_COMBINED + "/{soup_sort}/DGE_unfiltered/count_matrix.mtx"),
        cells=ancient(SPIPE_COMBINED + "/{soup_sort}/DGE_unfiltered/cell_metadata.csv"),
        genes=ancient(SPIPE_COMBINED + "/{soup_sort}/DGE_unfiltered/all_genes.csv"),
        called=ancient(CT_DGE + "/cell_metadata.csv"),
        ref_genes=ancient(CT_DGE + "/all_genes.csv"),
        sublibs=config["sublibraries"],
    output:
        profile=SOUP_DIR + "/{soup_sort}/soup_profile.tsv.gz",
        by_sublib=SOUP_DIR + "/{soup_sort}/soup_by_sublibrary.tsv.gz",
        barcodes=SOUP_DIR + "/{soup_sort}/barcode_summary.tsv",
        consistency=SOUP_DIR + "/{soup_sort}/sublibrary_consistency.tsv",
        top=SOUP_DIR + "/{soup_sort}/top_soup_genes.tsv",
        report=SOUP_DIR + "/{soup_sort}/soup_report.pdf",
    params:
        script=os.path.join(CT_SCRIPTS, "soup_profile.R"),
        args=soup_profile_args,
        **_CT_IMAGE,
    threads: nthreads("soup_profile")
    resources:
        mem_mb=16000,
        runtime=360,
    log:
        "logs/soupx/soup_profile_{soup_sort}.log",
    shell:
        _CT_SHELL


# --- 2. rho per sort ---------------------------------------------------------
def soupx_dir(run):
    return ct_dir(run) + "/soupx"


def soupx_labels(run):
    return ct_dir(run) + "/" + SX["labels_file"]


def soupx_rho_args(wc):
    return ct_fmt({
        "sort": wc.soup_sort,
        "labels": soupx_labels(wc.ct_run),
        "prep_dir": ct_dir(wc.ct_run) + "/prep",
        "doublet_dir": ct_dir(wc.ct_run) + "/doublets_by_sort",
        "sublibraries": config["sublibraries"],
        "soup_profile": SOUP_DIR + "/" + wc.soup_sort + "/soup_profile.tsv.gz",
        "barcode_summary": SOUP_DIR + "/" + wc.soup_sort + "/barcode_summary.tsv",
        "marker_sets": ";".join("%s=%s" % (k, "|".join(v))
                                for k, v in SX["marker_sets"].items()),
        "estimate_class": ",".join(SX["estimate_class"]),
        "exclude_subclasses": ",".join(SX.get("exclude_subclasses") or []),
        "min_label_score": SX["min_label_score"],
        "max_cells": SX["max_cells"],
        "min_est_cells": SX["min_est_cells"],
        "min_est_frac": SX["min_est_frac"],
        "max_contamination": SX["max_contamination"],
        "max_contamination_check": ",".join(str(v) for v in SX["max_contamination_check"]),
        "out_dir": soupx_dir(wc.ct_run) + "/rho",
        "seed": CT["seed"],
    })


rule soupx_rho:
    input:
        # No rule makes this: copy a finished run's cell_annotations.csv here.
        labels=RESULTS + "/{ct_run}/" + SX["labels_file"],
        # ancient(): upstream products; their timestamps must not pull prep or
        # the doublet jobs back into the DAG from here.
        counts=ancient(expand(RESULTS + "/{{ct_run}}/prep/counts/{sublib}.rds", sublib=SUBLIBS)),
        doublets=ancient(expand(RESULTS + "/{{ct_run}}/doublets_by_sort/{sublib}.tsv", sublib=SUBLIBS)),
        profile=SOUP_DIR + "/{soup_sort}/soup_profile.tsv.gz",
        floor=SOUP_DIR + "/{soup_sort}/barcode_summary.tsv",
    output:
        rho=RESULTS + "/{ct_run}/soupx/rho/{soup_sort}.tsv",
        breakdown=RESULTS + "/{ct_run}/soupx/rho/{soup_sort}_breakdown.tsv",
        sensitivity=RESULTS + "/{ct_run}/soupx/rho/{soup_sort}_sensitivity.tsv",
        report=RESULTS + "/{ct_run}/soupx/rho/{soup_sort}_report.pdf",
    params:
        script=os.path.join(CT_SCRIPTS, "soupx_rho.R"),
        args=soupx_rho_args,
        **_CT_IMAGE,
    threads: nthreads("soupx_rho")
    resources:
        mem_mb=32000,
        runtime=240,
    log:
        "logs/soupx/{ct_run}/rho_{soup_sort}.log",
    shell:
        _CT_SHELL


# --- 3. corrected counts per sublibrary --------------------------------------
def soupx_rho_override():
    """soupx.rho: null -> use the estimates; {sort: value} -> override those."""
    r = SX.get("rho")
    return ",".join("%s:%s" % kv for kv in r.items()) if r else "none"


def soupx_adjust_args(wc):
    return ct_fmt({
        "counts": ct_dir(wc.ct_run) + "/prep/counts/" + wc.sublib + ".rds",
        "labels": soupx_labels(wc.ct_run),
        "sublib": wc.sublib,
        "rho_dir": soupx_dir(wc.ct_run) + "/rho",
        "soup_dir": SOUP_DIR,
        "rho_override": soupx_rho_override(),
        "out_counts": soupx_dir(wc.ct_run) + "/counts/" + wc.sublib + ".rds",
        "out_cells": soupx_dir(wc.ct_run) + "/cells/" + wc.sublib + ".tsv",
        "seed": CT["seed"],
    })


rule soupx_adjust:
    input:
        counts=ancient(RESULTS + "/{ct_run}/prep/counts/{sublib}.rds"),
        labels=RESULTS + "/{ct_run}/" + SX["labels_file"],
        rho=expand(RESULTS + "/{{ct_run}}/soupx/rho/{soup_sort}.tsv", soup_sort=SOUP_SORTS),
        profiles=expand(SOUP_DIR + "/{soup_sort}/soup_profile.tsv.gz", soup_sort=SOUP_SORTS),
    output:
        counts=RESULTS + "/{ct_run}/soupx/counts/{sublib}.rds",
        cells=RESULTS + "/{ct_run}/soupx/cells/{sublib}.tsv",
    params:
        script=os.path.join(CT_SCRIPTS, "soupx_adjust.R"),
        args=soupx_adjust_args,
        **_CT_IMAGE,
    threads: nthreads("soupx_adjust")
    resources:
        mem_mb=24000,
        runtime=360,
    log:
        "logs/soupx/{ct_run}/adjust_{sublib}.log",
    shell:
        _CT_SHELL


# --- 4. re-cluster on the corrected counts ------------------------------------
# Same script and settings as celltype_cluster, with four differences, all
# passed as arguments (see the header of celltype_cluster.R):
#   counts        soupx/counts/ instead of prep/counts/ (QC metrics stay raw)
#   doublets      doublets_by_sort/ (per-sort rates) instead of doublets/
#   minimum size  fail_below_floor: raw UMIs < min_floor_k x the floor
#   output        <results_dir>/celltype_soupx/, beside the first pass
# Full run only: it reuses the full run's prep/, doublets_by_sort/ and soupx/.
SOUPX_RUN = "celltype_soupx"

rule celltype_cluster_soupx:
    input:
        meta=ancient(RESULTS + "/celltype/prep/cell_meta.rds"),
        counts=expand(soupx_dir("celltype") + "/counts/{sublib}.rds", sublib=SUBLIBS),
        cells=expand(soupx_dir("celltype") + "/cells/{sublib}.tsv", sublib=SUBLIBS),
        doublets=ancient(expand(RESULTS + "/celltype/doublets_by_sort/{sublib}.tsv", sublib=SUBLIBS)),
        floors=expand(SOUP_DIR + "/{soup_sort}/barcode_summary.tsv", soup_sort=SOUP_SORTS),
    output:
        **ct_cluster_outputs(SOUPX_RUN),
    params:
        script=os.path.join(CT_SCRIPTS, "celltype_cluster.R"),
        args=ct_cluster_args(
            "celltype",
            out_dir=ct_dir(SOUPX_RUN),
            counts_dir=soupx_dir("celltype"),
            doublet_dir=ct_dir("celltype") + "/doublets_by_sort",
            floor_files=",".join(expand(SOUP_DIR + "/{s}/barcode_summary.tsv", s=SOUP_SORTS)),
            min_floor_k=SX["min_floor_k"],
            soupx_cells_dir=soupx_dir("celltype") + "/cells",
        ),
        **_CT_IMAGE,
    threads: nthreads("celltype_cluster_soupx")
    resources:
        mem_mb=128000,
        runtime=720,
    log:
        "logs/celltype/celltype_soupx/cluster.log",
    shell:
        _CT_SHELL


# --- targets -----------------------------------------------------------------
rule soup_profiles:
    input:
        expand(SOUP_DIR + "/{soup_sort}/soup_report.pdf", soup_sort=SOUP_SORTS),


rule soupx_rhos:
    input:
        expand(soupx_dir("celltype") + "/rho/{soup_sort}_report.pdf", soup_sort=SOUP_SORTS),


rule soupx_rhos_test:
    input:
        expand(soupx_dir("celltype_test") + "/rho/{soup_sort}_report.pdf", soup_sort=SOUP_SORTS),


rule soupx:
    input:
        expand(soupx_dir("celltype") + "/cells/{sublib}.tsv", sublib=SUBLIBS),


rule soupx_test:
    input:
        expand(soupx_dir("celltype_test") + "/cells/{sublib}.tsv", sublib=SUBLIBS),


rule celltype_soupx:
    input:
        ct_dir(SOUPX_RUN) + "/cell_annotations.csv",
