# Run the pipeline from a login node, inside a screen (or tmux) session:
#
#   screen -S parse
#   source submit.sh                    # full run (rule all)
#   source submit.sh reference          # build the shared Parse reference only
#   source submit.sh all rerun          # tag this submission "all_rerun"
#   # detach with ctrl-a d ; reattach with `screen -r parse`
#
# UCL does not allow a cluster job to submit further jobs, so the snakemake
# controller runs HERE, in the foreground, and fans each rule out as its own
# qsub (see config/cluster.yaml). The screen session must stay alive for the
# duration of the run.
#
# This file is meant to be SOURCED, so it deliberately avoids `set -e`, `set -u`
# and `exit` -- in a sourced script each of those acts on your interactive
# shell, and an error would kill the screen session mid-run.
#
# split-pipe does not need to be on PATH: the rules that call it declare
# `conda: <path to the spipe env>` (config.conda_envs.spipe) and snakemake
# activates it per job. Only snakemake itself is needed here.
#
# Env providing snakemake (7.x -- --cluster-config was removed in 8).
# source activate /path/to/snakemake_env

SM_TARGET="${1:-all}"

if [ -n "${2:-}" ]; then
    SM_RUN_NAME="${SM_TARGET}_${2}"
else
    SM_RUN_NAME="${SM_TARGET}"
fi

# Snapshot the config that produced this run, next to its job logs.
SM_FOLDER="submissions/$(date +"%Y%m%d%H%M")_${SM_RUN_NAME}"
mkdir -p "${SM_FOLDER}"
cp config/config.yaml  "${SM_FOLDER}/${SM_RUN_NAME}_config.yaml"
cp config/cluster.yaml "${SM_FOLDER}/${SM_RUN_NAME}_cluster.yaml"

echo "[submit] target=${SM_TARGET}  logs=${SM_FOLDER}"

# -j 6 caps CONCURRENT cluster jobs, not total ones -- all 16 sublibraries go in
# with one submission and snakemake feeds them through six at a time. split-pipe
# writes ~900 GB per sublibrary (measured on ParseFANS1), so running all 16
# abreast would put ~10 TB of concurrent writes on shared /SAN and inflate every
# job's runtime through contention. Six keeps peak intermediates near 1-1.5 TB
# and per-job runtimes near their measured baseline. Queue time does not count
# against SGE h_rt -- it starts when the job starts -- so throttling costs
# wall-clock but puts no job at risk of being killed.

snakemake "${SM_TARGET}" \
-s workflow/Snakefile \
--conda-prefix "/SAN/vyplab/vyplab_reference_genomes/conda_envs/" \
--use-conda \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o ${SM_FOLDER} {cluster.submission_string}" \
-j 6 \
--nolock \
--rerun-incomplete \
--latency-wait 100 \
2>&1 | tee "${SM_FOLDER}/snakemake.log"
