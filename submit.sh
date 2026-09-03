#!/bin/bash
# Submit the pipeline to the UCL CS cluster (SGE).
#
#   qsub submit.sh                    # full run (rule all)
#   qsub submit.sh reference          # build the shared Parse reference only
#   qsub submit.sh reference gencode  # ... and tag this submission "gencode"
#
# split-pipe must be on PATH inside the JOBS, not just here. Both this script
# and cluster_qsub.sh use `#$ -V`, so the environment active at qsub time is
# exported to the jobs -- activate the spipe env below (or before qsub) and
# split-pipe is available throughout.

#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=48:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y

set -euo pipefail

# Env providing split-pipe (v1.9.0). snakemake must also be importable -- either
# from this env or from one already active when you qsub.
SPIPE_ENV="/SAN/vyplab/FUS_RNA/conda_envs/spipe"
# source activate "${SPIPE_ENV}"

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
