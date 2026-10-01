#$ -S /bin/bash
#$ -cwd
#$ -N marker_dotplot
#$ -j y
#$ -o logs/marker_dotplot/
#$ -pe smp 2
#$ -R y
#$ -l h_vmem=14G,tmem=14G,h_rt=2:00:00,tscratch=2G
# 2 x 14G = 28G. Peak is readRDS of the whole object: it was saved with
# compress = FALSE, so it needs about its file size in RAM (15 GB for
# celltype_soupx), plus a few GB for the marker-gene subset. Everything after
# that is small. If the rds grows, re-check with
#   ls -lh <run_dir>/seurat_harmony.rds
# and keep 2 x h_vmem >= ~1.4 x that size. The script is single-threaded;
# the 2 slots are only there to split the memory request.
#
# Canonical marker dot plots by Azimuth predicted.subclass x FANS sort, on the
# SoupX-corrected clustering. Same Singularity image + R library as the
# celltype rules (config/config.yaml: celltype.container / r_libs).
#
#   mkdir -p logs/marker_dotplot && qsub marker_dotplot_qsub.sh
#
# Optional overrides: qsub -v MIN_SCORE=0.9,MIN_CELLS=20 marker_dotplot_qsub.sh

set -euo pipefail

RUN_DIR=/SAN/neuroscience/direct_nanopore_rna/ariana_parse/spipe_jh_donor_annotations/celltype_soupx
RDS=${RUN_DIR}/seurat_harmony.rds
OUT_DIR=${RUN_DIR}/marker_dotplots
SIF=/SAN/vyplab/alb_projects/tools/seurat.sif
RLIB=/SAN/vyplab/alb_projects/tools/R_celltype_lib
SCRIPT=$(pwd)/workflow/scripts/marker_dotplot.R

echo "job ${JOB_ID:-local} on $(hostname)"
command -v singularity > /dev/null || { echo "singularity not on PATH"; exit 1; }
[ -f "$RDS" ] || { echo "missing $RDS"; exit 1; }

# Same bind logic as celltype.smk: TMPDIR must be visible inside the image.
TMPDIR="${TMPDIR:-$(pwd)/logs/marker_dotplot}"
BINDS=/SAN
case "$TMPDIR" in
    /SAN/*) ;;
    *) BINDS="$BINDS,$TMPDIR" ;;
esac
# The script lives in the repo; bind it if the checkout is not under /SAN.
case "$SCRIPT" in
    /SAN/*) ;;
    *) BINDS="$BINDS,$(dirname "$SCRIPT")" ;;
esac

singularity exec --cleanenv \
    -B "$BINDS" \
    --env R_LIBS=${RLIB},TMPDIR="$TMPDIR",OMP_NUM_THREADS=2,OPENBLAS_NUM_THREADS=2 \
    "$SIF" \
    Rscript "$SCRIPT" \
        --rds="$RDS" \
        --out_dir="$OUT_DIR" \
        --min_score="${MIN_SCORE:-0.8}" \
        --min_cells="${MIN_CELLS:-10}"
