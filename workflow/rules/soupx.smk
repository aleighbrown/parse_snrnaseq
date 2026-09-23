# -----------------------------------------------------------------------------
# Ambient RNA (SoupX), built up in steps. So far:
#
#   1. soup_profile   3 jobs    one per FANS sort: sum the empty barcodes of
#                               splitpipe/combined/<sort>/DGE_unfiltered
#
#   source submit.sh soup_profiles
#
# Planned: soupx_rho (one contamination fraction per sort, pooled over
# sublibraries, from marker genes and the celltype labels), then soupx_adjust
# (16 jobs, corrected counts per sublibrary).
#
# Why per sort: all donors were pooled within each sort before loading, so each
# sort had one suspension and one soup. Ambient RNA is reverse-transcribed in
# the round-1 wells, so it carries the sort's round-1 barcode, which is how
# split-pipe already split it into combined/<sort>/. Sublibraries are split
# after pooling and each holds all three soups.
#
# soup_profile depends only on split-pipe output, not on the celltype stages.
#
# Outputs, under <results_dir>/soup/<sort>/:
#   soup_profile.tsv.gz         pooled per gene: counts, est (SoupX's format)
#   soup_by_sublibrary.tsv.gz   the same per sublibrary
#   barcode_summary.tsv, sublibrary_consistency.tsv, top_soup_genes.tsv
#   soup_report.pdf             barcode ranks with the empty range shaded
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

localrules: soup_profiles


def soup_profile_args(wc):
    return ct_fmt({
        "dge_dir": os.path.join(SPIPE_COMBINED, wc.soup_sort, "DGE_unfiltered"),
        "called_cells": CT_DGE + "/cell_metadata.csv",
        "ref_genes": CT_DGE + "/all_genes.csv",
        "sublibraries": config["sublibraries"],
        "sort": wc.soup_sort,
        "out_dir": SOUP_DIR + "/" + wc.soup_sort,
        "umi_min": SX["empty_umi_range"][0],
        "umi_max": SX["empty_umi_range"][1],
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


# --- targets -----------------------------------------------------------------
rule soup_profiles:
    input:
        expand(SOUP_DIR + "/{soup_sort}/soup_report.pdf", soup_sort=SOUP_SORTS),
