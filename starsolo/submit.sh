# Run from starsolo/ on a login node, inside screen (same pattern as ../submit.sh):
#
#   screen -S starsolo
#   source submit.sh tools           # check python/pysam/sinto/pigz/STAR (seconds)
#   source submit.sh                 # pilot (rule all)
#   source submit.sh layout          # barcode read layout only
#   source submit.sh reference       # STAR index only
#   source submit.sh split           # donor x sort BAMs
#   source submit.sh all rerun       # tag this submission "all_rerun"
#
# SOURCED, not qsub-ed, so no `set -e` / `set -u` / `exit` (they would act on
# your interactive shell). Needs snakemake 7.x on PATH.

SM_TARGET="${1:-all}"

if [ -n "${2:-}" ]; then
    SM_RUN_NAME="${SM_TARGET}_${2}"
else
    SM_RUN_NAME="${SM_TARGET}"
fi

SM_FOLDER="submissions/$(date +"%Y%m%d%H%M")_${SM_RUN_NAME}"
mkdir -p "${SM_FOLDER}"
cp config/config.yaml  "${SM_FOLDER}/${SM_RUN_NAME}_config.yaml"
cp config/cluster.yaml "${SM_FOLDER}/${SM_RUN_NAME}_cluster.yaml"

echo "[submit] target=${SM_TARGET}  logs=${SM_FOLDER}"

# No conda: jobs use config tools.bin_dirs and tools.star (see config.yaml).
# `source submit.sh tools` checks them on the login node in seconds.
if [ "${SM_TARGET}" = "tools" ]; then
    snakemake tools -s Snakefile -j 1 --nolock
    return 2>/dev/null || true
fi

snakemake "${SM_TARGET}" \
-s Snakefile \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o ${SM_FOLDER} {cluster.submission_string}" \
-j 4 \
--nolock \
--rerun-incomplete \
--latency-wait 100 \
2>&1 | tee "${SM_FOLDER}/snakemake.log"
