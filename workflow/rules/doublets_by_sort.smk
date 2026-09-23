# -----------------------------------------------------------------------------
# Alternative doublet arm: scDblFinder per FANS sort within each sublibrary, at
# a rate worked out from that sort's own barcode occupancy. Runs BESIDE the
# default arm (celltype_doublets, one fixed rate per sublibrary) and never
# feeds celltype_cluster, so the existing annotations are untouched.
#
#   1. celltype_doublets_by_sort   16 jobs   <run>/prep/counts/<sublib>.rds
#                                  -> <run>/doublets_by_sort/<sublib>.tsv
#                                     <run>/doublets_by_sort/<sublib>.rates.tsv
#   2. doublets_compare            ONE job   both arms + cell_meta.rds + cell_annotations.csv
#                                  -> <run>/doublets_compare/
#
#   source submit.sh doublets_by_sort_test   # on the celltype_test cells
#   source submit.sh doublets_by_sort        # all cells
#
# Why per sort. Round-1 wells hold one sort each (config/parse_sample_list.txt),
# and two nuclei share a barcode only if they share round-1, -2 and -3 wells in
# the same sublibrary. So K_negative and TDP_negative (43 wells, ~2-5k cells
# per sublibrary) almost never collide, while positive (10 wells, ~20k cells)
# does so often. The default arm's 0.05 +/- 0.03 assumes the full 96^3 space.
# Running each sort separately also builds scDblFinder's artificial doublets
# only from pairs that can really occur.
#
# The occupancy is counted from split-pipe's own cell_metadata.csv (every
# called barcode), not from the prep output, which is QC-filtered and, in the
# test run, a random subset.
#
# In the test run each sublibrary holds ~1.2k sampled cells, so the negative
# sorts fall below celltype.doublets_by_sort.min_cells and stay uncalled; only
# positive is really exercised. That is expected -- the test is for plumbing.
#
# Uses CT, CT_DGE, CT_RUNS, ct_dir, ct_fmt, _CT_SHELL and _CT_IMAGE from
# celltype.smk, so it is included after it.
# -----------------------------------------------------------------------------

DBS = CT["doublets_by_sort"]


def _plate_wells(rng):
    """'A1-D7' -> round-1 wells covered, filling the 8 x 12 plate row-major
    (A1..A12, B1..); comma-separated ranges are summed."""
    def idx(w):
        return (ord(w[0].upper()) - ord("A")) * 12 + int(w[1:])
    n = 0
    for part in rng.split(","):
        a, _, b = part.strip().partition("-")
        n += idx(b or a) - idx(a) + 1
    return n


with open(config["splitpipe"]["sample_list"]) as fh:
    DBS_SORT_WELLS = {
        ln.split()[0]: sum(_plate_wells(r) for r in ln.split()[1:])
        for ln in fh if ln.strip() and not ln.startswith("#")
    }

localrules: doublets_by_sort, doublets_by_sort_test


def dbs_dir(run):
    return ct_dir(run) + "/doublets_by_sort"


def dbs_args(wc):
    return ct_fmt({
        "counts": ct_dir(wc.ct_run) + "/prep/counts/" + wc.sublib + ".rds",
        "called_cells": CT_DGE + "/cell_metadata.csv",
        "sublibraries": config["sublibraries"],
        "sublib": wc.sublib,
        "sort_wells": ",".join("%s:%d" % kv for kv in DBS_SORT_WELLS.items()),
        "barcodes_per_r1_well": DBS["barcodes_per_r1_well"],
        "dbr_sd_rel": DBS["dbr_sd_rel"],
        "min_cells": DBS["min_cells"],
        "out": dbs_dir(wc.ct_run) + "/" + wc.sublib + ".tsv",
        "rates_out": dbs_dir(wc.ct_run) + "/" + wc.sublib + ".rates.tsv",
        "seed": CT["seed"],
    })


def dbs_compare_args(wc):
    return ct_fmt({
        "prep_dir": ct_dir(wc.ct_run) + "/prep",
        "default_dir": ct_dir(wc.ct_run) + "/doublets",
        "by_sort_dir": dbs_dir(wc.ct_run),
        "annotations": ct_dir(wc.ct_run) + "/cell_annotations.csv",
        "sublibraries": config["sublibraries"],
        "default_dbr": CT["doublets"]["dbr"],
        "out_dir": ct_dir(wc.ct_run) + "/doublets_compare",
    })


# --- 1. scDblFinder per sort, one job per sublibrary -------------------------
rule celltype_doublets_by_sort:
    input:
        counts=RESULTS + "/{ct_run}/prep/counts/{sublib}.rds",
        called=ancient(CT_DGE + "/cell_metadata.csv"),
        samples=config["splitpipe"]["sample_list"],
    output:
        calls=RESULTS + "/{ct_run}/doublets_by_sort/{sublib}.tsv",
        rates=RESULTS + "/{ct_run}/doublets_by_sort/{sublib}.rates.tsv",
    params:
        script=os.path.join(CT_SCRIPTS, "doublets_by_sort.R"),
        args=dbs_args,
        **_CT_IMAGE,
    threads: nthreads("celltype_doublets_by_sort")
    resources:
        mem_mb=16000,
        runtime=480,
    log:
        "logs/celltype/{ct_run}/doublets_by_sort_{sublib}.log",
    shell:
        _CT_SHELL


# --- 2. compare with the default arm -----------------------------------------
# Needs cell_annotations.csv for cell-type labels and the UMAP, so it runs
# after celltype_cluster (which reads only the default arm).
rule doublets_compare:
    input:
        meta=RESULTS + "/{ct_run}/prep/cell_meta.rds",
        default=expand(RESULTS + "/{{ct_run}}/doublets/{sublib}.tsv", sublib=SUBLIBS),
        by_sort=expand(RESULTS + "/{{ct_run}}/doublets_by_sort/{sublib}.tsv", sublib=SUBLIBS),
        rates=expand(RESULTS + "/{{ct_run}}/doublets_by_sort/{sublib}.rates.tsv", sublib=SUBLIBS),
        annotations=RESULTS + "/{ct_run}/cell_annotations.csv",
    output:
        calls=RESULTS + "/{ct_run}/doublets_compare/cell_doublet_calls.tsv.gz",
        rates=RESULTS + "/{ct_run}/doublets_compare/rates_by_sublibrary_sort.tsv",
        by_sort=RESULTS + "/{ct_run}/doublets_compare/summary_by_sort.tsv",
        by_type=RESULTS + "/{ct_run}/doublets_compare/summary_by_sort_celltype.tsv",
        report=RESULTS + "/{ct_run}/doublets_compare/doublets_compare.pdf",
    params:
        script=os.path.join(CT_SCRIPTS, "doublets_compare.R"),
        args=dbs_compare_args,
        **_CT_IMAGE,
    threads: nthreads("doublets_compare")
    resources:
        mem_mb=24000,
        runtime=120,
    log:
        "logs/celltype/{ct_run}/doublets_compare.log",
    shell:
        _CT_SHELL


# --- targets -----------------------------------------------------------------
rule doublets_by_sort:
    input:
        ct_dir("celltype") + "/doublets_compare/doublets_compare.pdf",


rule doublets_by_sort_test:
    input:
        ct_dir("celltype_test") + "/doublets_compare/doublets_compare.pdf",
