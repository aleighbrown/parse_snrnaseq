# -----------------------------------------------------------------------------
# Cell type annotation of the combined split-pipe DGE, in three stages.
#
#   1. celltype_prep       ONE job      mtx + donor table -> QC flags for every
#                                       cell, raw counts split by sublibrary
#   2. celltype_doublets   16 jobs      scDblFinder, one sublibrary each
#   3. celltype_cluster    ONE big job  final QC -> Seurat -> Harmony -> UMAP +
#                                       clusters -> Azimuth + marker labels
#                                       -> cell_annotations.csv (every cell)
#
#   source submit.sh celltype_test   # 20k random cells, small jobs throughout
#   source submit.sh celltype        # all 426k cells
#
# Why three stages. The first full-size design was one job. The 20k-cell test
# showed doublet calling taking ~60% of the time, at a few GB per sublibrary,
# with nothing shared between sublibraries -- so it runs as 16 narrow jobs,
# which start quickly, and the one large reservation (stage 3) is held for a
# few hours instead of most of a day. It also means changing clustering
# settings reruns stage 3 only.
#
# Outputs, under <results_dir>/<run>/ with run = celltype | celltype_test:
#   prep/cell_meta.rds, prep/counts/<sublib>.rds   (~13 GB for the full run)
#   doublets/<sublib>.tsv
#   cell_annotations.csv, cluster_summary.tsv, cluster_markers.tsv,
#   qc_summary.tsv, celltype_report.pdf, seurat_harmony.rds
#
# Memory, full run (426,041 cells x 62,696 genes, 1.12e9 non-zeros). The 20k
# test peaked at 5.1 GB of R heap in the single-job version; scaled by
# non-zeros that is ~90-100 GB, almost all of it in stage 3.
#   prep:     3 int columns (~13 GB) -> dgCMatrix (~13 GB); ~30 GB peak,
#             plus fread's memory map of the 16 GB file (counts to h_vmem)
#   doublets: one ~26k-cell sublibrary + ~20k artificial doublets; a few GB
#   cluster:  counts + log-norm layers (~26 GB), scale.data, SNN graph,
#             Harmony, UMAP, markers; ~90-100 GB, plus ~10 GB per Azimuth
#             worker. Every stage logs R's heap at each step and its peak --
#             trim config/cluster.yaml from those.
#
# Not part of `rule all`: that target is split-pipe + demux only.
# -----------------------------------------------------------------------------
import os
import shlex

CT = config["celltype"]
CT_DGE = os.path.join(
    RESULTS, "splitpipe", "combined", os.path.dirname(config["splitpipe"]["cell_meta"])
)
CT_SCRIPTS = os.path.join(workflow.basedir, "scripts")
CT_DONORS = RESULTS + "/demux/donor_assignments.all_sublibraries.tsv"

# R runs inside a Singularity image, not a conda env: the login node's glibc
# is too old for current conda-forge R builds. The image (satijalab/seurat)
# supplies R 4.5 / Seurat 5 / harmony; Azimuth, humancortexref, scDblFinder
# and RhpcBLASctl live in a separate library on SAN, put first on the
# library path through R_LIBS. Called with `singularity exec` in the shell
# block, not `container:`, so submit.sh needs no --use-singularity and the
# other rules are untouched.
CT_SIF = CT["container"]
CT_RLIB = CT["r_libs"]
CT_BINDS = ",".join(CT.get("binds") or ["/SAN"])

# Output directory name -> cells to sample (0 = all).
CT_RUNS = {"celltype": 0, "celltype_test": CT["test_n_cells"]}

wildcard_constraints:
    ct_run="|".join(CT_RUNS),

localrules: celltype, celltype_test


def ct_dir(run):
    return RESULTS + "/" + run


def ct_prep_outputs(run):
    return {
        "meta": ct_dir(run) + "/prep/cell_meta.rds",
        "counts": expand(ct_dir(run) + "/prep/counts/{sublib}.rds", sublib=SUBLIBS),
    }


def ct_cluster_inputs(run):
    return {
        **ct_prep_outputs(run),
        "doublets": expand(ct_dir(run) + "/doublets/{sublib}.tsv", sublib=SUBLIBS),
    }


CT_OUTPUTS = {
    "annotations": "cell_annotations.csv",
    "clusters": "cluster_summary.tsv",
    "markers": "cluster_markers.tsv",
    "qc": "qc_summary.tsv",
    "rds": "seurat_harmony.rds",
    "report": "celltype_report.pdf",
}


def ct_cluster_outputs(run):
    return {k: ct_dir(run) + "/" + v for k, v in CT_OUTPUTS.items()}


def ct_fmt(args):
    """--key=value arguments for the celltype_*.R scripts (--threads is
    added in the shell block from {threads}, so it always matches the job).

    %-format, not an f-string: Snakemake 7 on Python 3.12 re-tokenises
    f-strings in .smk files and pads them with spaces.
    """
    return " ".join(shlex.quote("--%s=%s" % (k, v)) for k, v in args.items())


def ct_qc_args():
    qc = CT["qc"]
    return {
        "min_genes": qc["min_genes"],
        "max_genes": qc["max_genes"],
        "min_counts": qc["min_counts"],
        "max_counts": qc["max_counts"] or 0,
        "max_pct_mt": qc["max_pct_mt"],
    }


def ct_prep_args(run):
    return ct_fmt({
        "dge_dir": CT_DGE,
        "donor_table": CT_DONORS,
        "sublibraries": config["sublibraries"],
        "out_dir": ct_dir(run) + "/prep",
        "seed": CT["seed"],
        "n_cells": CT_RUNS[run],
        "only_assigned": CT["only_assigned"],
        **ct_qc_args(),
    })


def ct_doublet_args(wc):
    dbl = CT["doublets"]
    return ct_fmt({
        "counts": ct_dir(wc.ct_run) + "/prep/counts/" + wc.sublib + ".rds",
        "out": ct_dir(wc.ct_run) + "/doublets/" + wc.sublib + ".tsv",
        "seed": CT["seed"],
        "dbr": dbl["dbr"],
        "dbr_sd": dbl["dbr_sd"],
    })


def ct_cluster_args(run):
    cl, az = CT["clustering"], CT["azimuth"]
    return ct_fmt({
        "prep_dir": ct_dir(run) + "/prep",
        "doublet_dir": ct_dir(run) + "/doublets",
        "sublibraries": config["sublibraries"],
        "out_dir": ct_dir(run),
        "seed": CT["seed"],
        "only_assigned": CT["only_assigned"],
        **ct_qc_args(),
        "min_cells_per_gene": CT["qc"]["min_cells_per_gene"],
        "remove_doublets": CT["doublets"]["remove"],
        "n_hvg": cl["n_hvg"],
        "n_pcs": cl["n_pcs"],
        "harmony_vars": ",".join(cl["harmony_vars"]),
        "cluster_resolution": cl["resolution"],
        "extra_resolutions": ",".join(str(r) for r in cl.get("extra_resolutions") or []),
        "run_azimuth": az["run"],
        "azimuth_reference": az["reference"],
        "azimuth_workers": az["workers"],
    })


_CT_SHELL = r"""
        set -euo pipefail
        command -v singularity > /dev/null || {{ echo "singularity not on PATH" > {log}; exit 1; }}

        # fread(cmd=) spools the smoke test's awk output (~1 GB) to R's
        # tempdir; use the job's scratch if SGE set one, else the log dir.
        # --cleanenv drops the host environment, so TMPDIR is passed (and,
        # when it is not already under a bound path, bound) explicitly.
        TMPDIR="${{TMPDIR:-$(cd "$(dirname {log})" && pwd)}}"
        BINDS="{params.binds}"
        case "$TMPDIR" in
            /SAN/*) ;;
            *) BINDS="$BINDS,$TMPDIR" ;;
        esac

        # R_LIBS: the image's Renviron appends its own libraries after it, so
        # the SAN library comes first and the image's Seurat is still used.
        # OMP / OpenBLAS threads = the job's slots (fresh worker processes
        # would otherwise start one thread per core on the node).
        singularity exec --cleanenv \
            -B "$BINDS" \
            --env R_LIBS={params.rlib},TMPDIR="$TMPDIR",OMP_NUM_THREADS={threads},OPENBLAS_NUM_THREADS={threads} \
            {params.sif} \
            Rscript {params.script} --threads={threads} {params.args} > {log} 2>&1
"""

_CT_IMAGE = dict(sif=CT_SIF, rlib=CT_RLIB, binds=CT_BINDS)


# --- stage 1: load + QC ------------------------------------------------------
rule celltype_prep:
    input:
        # ancient(): these are the product of a ~15 h split-pipe comb; a touched
        # timestamp must never make Snakemake want to redo it from here.
        mtx=ancient(CT_DGE + "/count_matrix.mtx"),
        genes=ancient(CT_DGE + "/all_genes.csv"),
        cells=ancient(CT_DGE + "/cell_metadata.csv"),
        donors=CT_DONORS,
        sublibs=config["sublibraries"],
    output:
        **ct_prep_outputs("celltype"),
    params:
        script=os.path.join(CT_SCRIPTS, "celltype_prep.R"),
        args=ct_prep_args("celltype"),
        **_CT_IMAGE,
    threads: nthreads("celltype_prep")
    resources:
        mem_mb=64000,
        runtime=240,
    log:
        "logs/celltype/celltype/prep.log",
    shell:
        _CT_SHELL


# Smoke test: same script on CT["test_n_cells"] random cells. Reading the
# subset still streams the whole 16 GB mtx once through awk (~5-7 min).
use rule celltype_prep as celltype_prep_test with:
    output:
        **ct_prep_outputs("celltype_test"),
    params:
        script=os.path.join(CT_SCRIPTS, "celltype_prep.R"),
        args=ct_prep_args("celltype_test"),
        **_CT_IMAGE,
    threads: nthreads("celltype_prep_test")
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        "logs/celltype/celltype_test/prep.log",


# --- stage 2: doublets, one job per sublibrary (full run and test alike) -----
rule celltype_doublets:
    input:
        RESULTS + "/{ct_run}/prep/counts/{sublib}.rds",
    output:
        RESULTS + "/{ct_run}/doublets/{sublib}.tsv",
    params:
        script=os.path.join(CT_SCRIPTS, "celltype_doublets.R"),
        args=ct_doublet_args,
        **_CT_IMAGE,
    threads: nthreads("celltype_doublets")
    resources:
        mem_mb=16000,
        runtime=480,
    log:
        "logs/celltype/{ct_run}/doublets_{sublib}.log",
    shell:
        _CT_SHELL


# --- stage 3: cluster + annotate ---------------------------------------------
rule celltype_cluster:
    input:
        **ct_cluster_inputs("celltype"),
    output:
        **ct_cluster_outputs("celltype"),
    params:
        script=os.path.join(CT_SCRIPTS, "celltype_cluster.R"),
        args=ct_cluster_args("celltype"),
        **_CT_IMAGE,
    threads: nthreads("celltype_cluster")
    resources:
        mem_mb=128000,
        runtime=720,
    log:
        "logs/celltype/celltype/cluster.log",
    shell:
        _CT_SHELL


use rule celltype_cluster as celltype_cluster_test with:
    input:
        **ct_cluster_inputs("celltype_test"),
    output:
        **ct_cluster_outputs("celltype_test"),
    params:
        script=os.path.join(CT_SCRIPTS, "celltype_cluster.R"),
        args=ct_cluster_args("celltype_test"),
        **_CT_IMAGE,
    threads: nthreads("celltype_cluster_test")
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        "logs/celltype/celltype_test/cluster.log",


# --- targets -----------------------------------------------------------------
rule celltype:
    input:
        ct_dir("celltype") + "/cell_annotations.csv",


rule celltype_test:
    input:
        ct_dir("celltype_test") + "/cell_annotations.csv",
