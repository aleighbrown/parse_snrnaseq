# -----------------------------------------------------------------------------
# split-pipe, the Parse pipeline proper:
#   --mode all   one job per sublibrary (UDI) -> aligned/annotated BAM +
#                filtered cell metadata, with Round-1 wells split into the
#                three FANS sorts via config/parse_sample_list.txt
#   --mode comb  merge all 16 sublibraries into one combined DGE matrix
#
# NOTE: split-pipe is Parse's proprietary binary (register + download from
# Parse Biosciences) and has no conda recipe, but it IS already installed in an
# environment. Both rules activate config.conda_envs.spipe in their own shell
# block -- NOT via Snakemake's `conda:` directive, which breaks on a path (see
# the comment inside the rule). split-pipe need not be on PATH at submission.
# -----------------------------------------------------------------------------

_chem = str(config["splitpipe"].get("chemistry") or "auto").strip()
# "auto"/blank -> omit the flag entirely and let split-pipe detect the kit.
CHEM_FLAG = "" if _chem.lower() in ("auto", "none", "") else f"--chemistry {_chem}"


rule split_pipe:
    input:
        fq1=fq1,
        fq2=fq2,
        genome=GENOME_DIR,
        sample_list=config["splitpipe"]["sample_list"],
    output:
        bam="results/splitpipe/{sublib}/" + config["splitpipe"]["bam_name"],
        cell_meta="results/splitpipe/{sublib}/" + config["splitpipe"]["cell_meta"],
    params:
        outdir="results/splitpipe/{sublib}",
        chem=CHEM_FLAG,
        extra=config["splitpipe"].get("all_extra", ""),
        spipe=config["conda_envs"]["spipe"],
    threads: nthreads("split_pipe")
    resources:
        mem_mb=64000,
        runtime=1440,
    log:
        "logs/splitpipe/{sublib}.log",
    shell:
        r"""
        set -euo pipefail
        # Activate the existing split-pipe env EXPLICITLY rather than via
        # Snakemake's `conda:` directive. Snakemake 7 treats any non-YAML
        # conda value as an environment NAME, and when recording provenance
        # runs `conda env export --name <value>`; conda cannot resolve a
        # path as a name, so that fails AFTER the rule's work is done, marks
        # the job failed, and deletes the output. Doing it here keeps
        # snakemake's conda bookkeeping out of the picture entirely.
        set +u
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "{params.spipe}"
        set -u

        split-pipe \
            --mode all \
            {params.chem} \
            --genome_dir {input.genome} \
            --fq1 {input.fq1} \
            --fq2 {input.fq2} \
            --samp_list {input.sample_list} \
            --output_dir {params.outdir} \
            --nthreads {threads} \
            {params.extra} > {log} 2>&1
        """


# Combine all sublibraries into a single DGE matrix (this is the object that
# goes on to Pegasus for QC / clustering).
rule split_pipe_comb:
    input:
        cell_meta=expand(
            "results/splitpipe/{sublib}/" + config["splitpipe"]["cell_meta"],
            sublib=SUBLIBS,
        ),
    output:
        cell_meta="results/splitpipe/combined/" + config["splitpipe"]["cell_meta"],
    params:
        sublib_dirs=expand("results/splitpipe/{sublib}", sublib=SUBLIBS),
        outdir="results/splitpipe/combined",
        spipe=config["conda_envs"]["spipe"],
    threads: nthreads("split_pipe_comb")
    resources:
        mem_mb=64000,
        runtime=480,
    log:
        "logs/splitpipe/combine.log",
    shell:
        r"""
        set -euo pipefail
        # Activate the existing split-pipe env EXPLICITLY rather than via
        # Snakemake's `conda:` directive. Snakemake 7 treats any non-YAML
        # conda value as an environment NAME, and when recording provenance
        # runs `conda env export --name <value>`; conda cannot resolve a
        # path as a name, so that fails AFTER the rule's work is done, marks
        # the job failed, and deletes the output. Doing it here keeps
        # snakemake's conda bookkeeping out of the picture entirely.
        set +u
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "{params.spipe}"
        set -u

        split-pipe \
            --mode comb \
            --sublibraries {params.sublib_dirs} \
            --output_dir {params.outdir} \
            --nthreads {threads} > {log} 2>&1
        """


# Extract the retained cell-barcode list in the exact string form that appears
# in the BAM's cell tag, so cellsnp-lite can match cells. The barcode column
# name is set in config (splitpipe.barcode_column) -- verify it once.
rule barcodes_from_meta:
    input:
        cell_meta="results/splitpipe/{sublib}/" + config["splitpipe"]["cell_meta"],
    output:
        barcodes="results/splitpipe/{sublib}/barcodes.tsv",
    params:
        col=config["splitpipe"]["barcode_column"],
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/extract_barcodes.py"


# Index the BAM if split-pipe did not.
rule index_bam:
    input:
        bam="results/splitpipe/{sublib}/" + config["splitpipe"]["bam_name"],
    output:
        bai="results/splitpipe/{sublib}/" + config["splitpipe"]["bam_name"] + ".bai",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "samtools index {input.bam}"
