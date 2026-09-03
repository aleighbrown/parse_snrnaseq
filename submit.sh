#!/bin/bash
# Submit the pipeline to the UCL CS cluster (SGE).
#
#   qsub submit.sh                    # full run (rule all)
#   qsub submit.sh reference          # build the shared Parse reference only
#   qsub submit.sh reference gencode  # ... and tag this submission "gencode"
#
# split-pipe does NOT need to be on PATH: the rules that call it declare
# `conda: <path to the spipe env>` (config.conda_envs.spipe), which Snakemake
# activates per job. Only snakemake itself has to be available here.

#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=48:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y

set -euo pipefail

# Env providing snakemake (7.x -- --cluster-config is gone in 8). The spipe env
# is not activated here; the split-pipe rules bring it in themselves.
# source activate /path/to/snakemake_env

TARGET="${1:-all}"

if [ "${2:-}" != "" ]; then
    RUN_NAME="${TARGET}_${2}"
else
    RUN_NAME="${TARGET}"
fi

# Snapshot the config that produced this run, next to its job logs.
FOLDER="submissions/$(date +"%Y%m%d%H%M")_${RUN_NAME}"
mkdir -p "${FOLDER}"
cp config/config.yaml "${FOLDER}/${RUN_NAME}_config.yaml"
cp config/cluster.yaml "${FOLDER}/${RUN_NAME}_cluster.yaml"

snakemake "${TARGET}" \
-s workflow/Snakefile \
--conda-prefix "/SAN/vyplab/vyplab_reference_genomes/conda_envs/" \
--use-conda \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o ${FOLDER} {cluster.submission_string}" \
-j 40 \
--nolock \
--rerun-incomplete \
--latency-wait 100
