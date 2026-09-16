# -----------------------------------------------------------------------------
# Cell type annotation of the combined split-pipe DGE.
#
#   splitpipe/combined/all-sample/DGE_filtered/{count_matrix.mtx, all_genes.csv,
#   cell_metadata.csv}  +  demux/donor_assignments.all_sublibraries.tsv
#     -> QC (genes / transcripts / % mito) -> scDblFinder -> Seurat -> Harmony
#     -> UMAP + clusters -> Azimuth (humancortexref) + marker-gene class
#     -> celltype/cell_annotations.csv  (one row per split-pipe cell)
#
#   source submit.sh celltype_test   # 20k-cell smoke test, small job
#   source submit.sh celltype        # all 426k cells, ONE big job
#
# ONE BIG JOB, deliberately. Clustering is a single shared-memory computation
# over every cell, so it cannot be fanned out, and on this cluster one large
# reservation (-R y) starts sooner than a chain of mid-sized ones each queueing
# in turn. Everything runs in workflow/scripts/celltype_annotation.R.
#
# Memory, for 426,041 cells x 62,696 genes, 1.12e9 non-zeros:
#   mtx read as 3 int columns          ~13 GB  (transient)
#   raw counts dgCMatrix               ~13 GB
#   + log-normalised layer             ~13 GB
#   scale.data, 2000 HVG x 426k        ~7 GB
#   SNN graph, Harmony, UMAP           ~10 GB
#   Azimuth, per sublibrary x workers  ~8-10 GB each
# -> ~60-90 GB steady, spikes above that at subsetting. cluster.yaml reserves
# 8 x 24G = 192G. The log prints R's heap at every step and the peak at the
# end -- trim the reservation once a full run has reported back.
#
# Not part of `rule all`: that target is split-pipe + demux only.
# -----------------------------------------------------------------------------
import os
import shlex

CT = config["celltype"]
CT_DGE = os.path.join(
    RESULTS, "splitpipe", "combined", os.path.dirname(config["splitpipe"]["cell_meta"])
)
CT_SCRIPT = os.path.join(workflow.basedir, "scripts", "celltype_annotation.R")

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

localrules: celltype, celltype_test

CT_OUTPUTS = {
    "annotations": "cell_annotations.csv",
    "clusters": "cluster_summary.tsv",
    "markers": "cluster_markers.tsv",
    "qc": "qc_summary.tsv",
    "rds": "seurat_harmony.rds",
    "report": "celltype_report.pdf",
}


def ct_outputs(out_dir):
    return {k: os.path.join(out_dir, v) for k, v in CT_OUTPUTS.items()}


def ct_args(out_dir, n_cells, threads):
    """--key=value arguments for celltype_annotation.R."""
    qc, dbl, cl, az = CT["qc"], CT["doublets"], CT["clustering"], CT["azimuth"]
    args = {
        "dge_dir": CT_DGE,
        "donor_table": RESULTS + "/demux/donor_assignments.all_sublibraries.tsv",
        "sublibraries": config["sublibraries"],
        "out_dir": out_dir,
        "threads": threads,
        "seed": CT["seed"],
        "n_cells": n_cells,
        "only_assigned": CT["only_assigned"],
        "min_genes": qc["min_genes"],
        "max_genes": qc["max_genes"],
        "min_counts": qc["min_counts"],
        "max_counts": qc["max_counts"] or 0,
        "max_pct_mt": qc["max_pct_mt"],
        "min_cells_per_gene": qc["min_cells_per_gene"],
        "doublet_dbr": dbl["dbr"],
        "doublet_dbr_sd": dbl["dbr_sd"],
        "remove_doublets": dbl["remove"],
        "n_hvg": cl["n_hvg"],
        "n_pcs": cl["n_pcs"],
        "harmony_vars": ",".join(cl["harmony_vars"]),
        "cluster_resolution": cl["resolution"],
        "extra_resolutions": ",".join(str(r) for r in cl.get("extra_resolutions") or []),
        "run_azimuth": az["run"],
        "azimuth_reference": az["reference"],
        "azimuth_workers": az["workers"],
    }
    # %-format, not an f-string: Snakemake 7 on Python 3.12 re-tokenises
    # f-strings in .smk files and pads them with spaces.
    return " ".join(shlex.quote("--%s=%s" % (k, v)) for k, v in args.items())


_CT_SHELL = r"""
        set -euo pipefail
        command -v singularity > /dev/null || {{ echo "singularity not on PATH" > {log}; exit 1; }}

        # fread(cmd=) spools the smoke test's awk output (~1 GB) to R's
        # tempdir; use the job's scratch if SGE set one, else the output dir.
        # --cleanenv drops the host environment, so TMPDIR is passed (and,
        # when it is not already under a bound path, bound) explicitly.
        TMPDIR="${{TMPDIR:-$(dirname {output.annotations})}}"
        BINDS="{params.binds}"
        case "$TMPDIR" in
            /SAN/*) ;;
            *) BINDS="$BINDS,$TMPDIR" ;;
        esac

        # R_LIBS: the image's Renviron appends its own libraries after it, so
        # the SAN library comes first and the image's Seurat is still used.
        singularity exec --cleanenv \
            -B "$BINDS" \
            --env R_LIBS={params.rlib},TMPDIR="$TMPDIR" \
            {params.sif} \
            Rscript {params.script} {params.args} > {log} 2>&1
"""


rule celltype_annotation:
    input:
        # ancient(): these are the product of a ~15 h split-pipe comb; a touched
        # timestamp must never make Snakemake want to redo it from here.
        mtx=ancient(CT_DGE + "/count_matrix.mtx"),
        genes=ancient(CT_DGE + "/all_genes.csv"),
        cells=ancient(CT_DGE + "/cell_metadata.csv"),
        donors=RESULTS + "/demux/donor_assignments.all_sublibraries.tsv",
        sublibs=config["sublibraries"],
    output:
        **ct_outputs(RESULTS + "/celltype"),
    params:
        script=CT_SCRIPT,
        sif=CT_SIF,
        rlib=CT_RLIB,
        binds=CT_BINDS,
        args=ct_args(RESULTS + "/celltype", 0, nthreads("celltype_annotation")),
    threads: nthreads("celltype_annotation")
    resources:
        mem_mb=192000,
        runtime=2880,
    log:
        "logs/celltype/celltype_annotation.log",
    shell:
        _CT_SHELL


# Same script on a random subset of cells, into its own directory, with its own
# (small) entry in cluster.yaml. Proves the image + library, the inputs and every step
# before the big reservation is spent. Reading the subset still streams the
# whole 16 GB mtx once through awk, so it is not instant.
use rule celltype_annotation as celltype_annotation_test with:
    output:
        **ct_outputs(RESULTS + "/celltype_test"),
    params:
        script=CT_SCRIPT,
        sif=CT_SIF,
        rlib=CT_RLIB,
        binds=CT_BINDS,
        args=ct_args(
            RESULTS + "/celltype_test",
            CT["test_n_cells"],
            nthreads("celltype_annotation_test"),
        ),
    threads: nthreads("celltype_annotation_test")
    resources:
        mem_mb=32000,
        runtime=720,
    log:
        "logs/celltype/celltype_annotation_test.log",


rule celltype:
    input:
        RESULTS + "/celltype/cell_annotations.csv",


rule celltype_test:
    input:
        RESULTS + "/celltype_test/cell_annotations.csv",
