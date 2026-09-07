# -----------------------------------------------------------------------------
# split-pipe, the Parse pipeline proper.
#
# In split-pipe the PIPELINE STEPS ARE MODES (run_modes.py:RUN_MODE_DEFS):
#
#     inqc  pre  align  post  mol  dge  sub  ana  rep  clean
#
# `--mode all` is just a wrapper that walks that list, and the help text
# (pb_help.py) documents running a subset as "start at <mode1> and use
# --until_step <mode2>". So the per-sublibrary work can be submitted as
# several right-sized cluster jobs instead of one 4-slot job that is
# CPU-starved for half its life and idle for the other half.
#
# Measured on ParseFANS1 (1.69e9 reads, 4 threads; process/mem_profile.txt):
#
#     step         hours   mean cores   peak RSS
#     InQC          0.03      0.86        11.2G
#     Preprocess    5.73      2.77        11.3G
#     Align (STAR)  9.71      3.87        17.9G   <- CPU-bound, saturates 4
#     PostProcess   3.28      1.87         4.3G   <- I/O-bound
#     Molinfo       1.09+     1.20         3.7G   <- I/O-bound
#
# Hence the three rules below: `pre` and `align` get threads (STAR is ~all of
# the align step and scales near-linearly), everything from `post` onward runs
# narrow. Set splitpipe.split_steps: false in config.yaml to fall back to the
# original single `--mode all` job.
#
#   --mode comb  merges all 16 sublibraries into one combined DGE matrix
#
# EVERY step invocation gets the SAME flags. This matters: this run's
# process/run_proc_def.json has an empty "parent" object, so split-pipe's
# init_parent() (spc_init.py:1280) cannot recover chemistry/kit from a parent
# dir -- they must be given again on each call, and sample definitions come
# from --samp_list every time. Only --mode / --until_step / --nthreads differ.
#
# NOTE: split-pipe is Parse's proprietary binary (register + download from
# Parse Biosciences) and has no conda recipe, but it IS already installed in an
# environment. Every rule activates config.conda_envs.spipe in its own shell
# block -- NOT via Snakemake's `conda:` directive, which breaks on a path (see
# the comment inside the first rule). split-pipe need not be on PATH.
# -----------------------------------------------------------------------------

# split-pipe v1.9.0 does NOT auto-detect chemistry for --mode all: with no
# --chemistry it exits instantly with "Problem: No chemistry set". Fail here,
# while building the DAG, rather than letting 16 jobs each discover it.
_chem = str(config["splitpipe"].get("chemistry") or "").strip()
if _chem.lower() in ("", "auto", "none", "null"):
    raise WorkflowError(
        "splitpipe.chemistry is not set.\n"
        "split-pipe v1.9.0 has no auto-detection for --mode all; it exits with "
        "'Problem: No chemistry set'.\n"
        "Run `split-pipe --explain` (spipe env active) to list the chemistry and "
        "kit values your build accepts, then set splitpipe.chemistry in "
        "config/config.yaml."
    )
CHEM_FLAG = f"--chemistry {_chem}"

# --kit is required too. ('kit_source: automatic' in split-pipe's parameter
# dump applies to mkref, where the kit is forced to 'custom' -- it does not mean
# --mode all can infer it.)
_kit = str(config["splitpipe"].get("kit") or "").strip()
if not _kit:
    raise WorkflowError(
        "splitpipe.kit is not set.\n"
        "split-pipe v1.9.0 exits with 'Problem: Kit must be specified'.\n"
        "Run `split-pipe --kit_list` (spipe env active) to list the kits your "
        "build accepts, then set splitpipe.kit in config/config.yaml."
    )
KIT_FLAG = f"--kit {_kit}"

# True  -> the three-rule step split below (default)
# False -> the original single `--mode all` job per sublibrary
SPLIT_STEPS = bool(config["splitpipe"].get("split_steps", True))

SP_DIR = RESULTS + "/splitpipe/{sublib}"


if SPLIT_STEPS:

    # -- inqc + pre ------------------------------------------------------------
    # Barcode correction. Spawns one process per thread ("Multiprocessing N
    # parts" in the log), so it does scale, but it also wrote ~674 GB for a
    # 100 GB sublibrary and is partly I/O-bound -- expect sub-linear returns.
    rule split_pipe_pre:
        input:
            fq1=fq1,
            fq2=fq2,
            genome=GENOME_DIR,
            sample_list=config["splitpipe"]["sample_list"],
        output:
            # A SENTINEL, not a real file. split-pipe's `clean` step deletes the
            # process/*.bam intermediates at the end of the run; if a rule
            # declared one of those as its output, Snakemake would see it
            # missing on the next invocation and re-run the whole chain.
            flag=touch(SP_DIR + "/.spipe_pre.done"),
        params:
            outdir=SP_DIR,
            chem=CHEM_FLAG,
            kit=KIT_FLAG,
            extra=config["splitpipe"].get("all_extra", ""),
            spipe=config["conda_envs"]["spipe"],
        threads: nthreads("split_pipe_pre")
        resources:
            mem_mb=32000,
            runtime=1440,
        log:
            "logs/splitpipe/{sublib}.pre.log",
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
                --mode inqc \
                --until_step pre \
                {params.chem} \
                {params.kit} \
                --genome_dir {input.genome} \
                --fq1 {input.fq1} \
                --fq2 {input.fq2} \
                --samp_list {input.sample_list} \
                --output_dir {params.outdir} \
                --nthreads {threads} \
                {params.extra} > {log} 2>&1
            """

    # -- align -----------------------------------------------------------------
    # Effectively pure STAR: the step took 9.71 h and STAR's own log accounted
    # for 9 h 40 m of it, at 3.87 of 4 cores. This is the one genuinely
    # CPU-bound step, and the reason the split exists.
    rule split_pipe_align:
        input:
            pre=SP_DIR + "/.spipe_pre.done",
            fq1=fq1,
            fq2=fq2,
            genome=GENOME_DIR,
            sample_list=config["splitpipe"]["sample_list"],
        output:
            flag=touch(SP_DIR + "/.spipe_align.done"),
        params:
            outdir=SP_DIR,
            chem=CHEM_FLAG,
            kit=KIT_FLAG,
            extra=config["splitpipe"].get("all_extra", ""),
            spipe=config["conda_envs"]["spipe"],
        threads: nthreads("split_pipe_align")
        resources:
            mem_mb=48000,
            runtime=2880,
        log:
            "logs/splitpipe/{sublib}.align.log",
        shell:
            r"""
            set -euo pipefail
            set +u
            source "$(conda info --base)/etc/profile.d/conda.sh"
            conda activate "{params.spipe}"
            set -u

            split-pipe \
                --mode align \
                --until_step align \
                {params.chem} \
                {params.kit} \
                --genome_dir {input.genome} \
                --fq1 {input.fq1} \
                --fq2 {input.fq2} \
                --samp_list {input.sample_list} \
                --output_dir {params.outdir} \
                --nthreads {threads} \
                {params.extra} > {log} 2>&1
            """

    # -- post -> clean ---------------------------------------------------------
    # PostProcess, Molinfo, DGE, SubSample, Analysis, Report, Clean. Measured at
    # 1.2-1.9 cores, so this runs narrow -- but DGE builds the cell x gene
    # matrices and has NOT yet been profiled, so keep the memory reservation
    # generous until one full run reports its peak.
    rule split_pipe_post:
        input:
            align=SP_DIR + "/.spipe_align.done",
            fq1=fq1,
            fq2=fq2,
            genome=GENOME_DIR,
            sample_list=config["splitpipe"]["sample_list"],
        output:
            cell_meta=SP_DIR + "/" + config["splitpipe"]["cell_meta"],
            flag=touch(SP_DIR + "/.spipe_post.done"),
        params:
            outdir=SP_DIR,
            chem=CHEM_FLAG,
            kit=KIT_FLAG,
            extra=config["splitpipe"].get("all_extra", ""),
            spipe=config["conda_envs"]["spipe"],
        threads: nthreads("split_pipe_post")
        resources:
            mem_mb=64000,
            runtime=4320,
        log:
            "logs/splitpipe/{sublib}.post.log",
        shell:
            r"""
            set -euo pipefail
            set +u
            source "$(conda info --base)/etc/profile.d/conda.sh"
            conda activate "{params.spipe}"
            set -u

            split-pipe \
                --mode post \
                --until_step clean \
                {params.chem} \
                {params.kit} \
                --genome_dir {input.genome} \
                --fq1 {input.fq1} \
                --fq2 {input.fq2} \
                --samp_list {input.sample_list} \
                --output_dir {params.outdir} \
                --nthreads {threads} \
                {params.extra} > {log} 2>&1
            """


else:

    # -- the original single-job form ------------------------------------------
    # Fallback: set splitpipe.split_steps: false to use this instead. Sized for
    # the whole run, so it holds its widest reservation for the entire ~28 h.
    rule split_pipe:
        input:
            fq1=fq1,
            fq2=fq2,
            genome=GENOME_DIR,
            sample_list=config["splitpipe"]["sample_list"],
        output:
            cell_meta=SP_DIR + "/" + config["splitpipe"]["cell_meta"],
        params:
            outdir=SP_DIR,
            chem=CHEM_FLAG,
            kit=KIT_FLAG,
            extra=config["splitpipe"].get("all_extra", ""),
            spipe=config["conda_envs"]["spipe"],
        threads: nthreads("split_pipe")
        resources:
            mem_mb=64000,
            runtime=7200,
        log:
            "logs/splitpipe/{sublib}.log",
        shell:
            r"""
            set -euo pipefail
            set +u
            source "$(conda info --base)/etc/profile.d/conda.sh"
            conda activate "{params.spipe}"
            set -u

            split-pipe \
                --mode all \
                {params.chem} \
                {params.kit} \
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
            RESULTS + "/splitpipe/{sublib}/" + config["splitpipe"]["cell_meta"],
            sublib=SUBLIBS,
        ),
    output:
        cell_meta=RESULTS + "/splitpipe/combined/" + config["splitpipe"]["cell_meta"],
    params:
        sublib_dirs=expand(RESULTS + "/splitpipe/{sublib}", sublib=SUBLIBS),
        outdir=RESULTS + "/splitpipe/combined",
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
        # See split_pipe_pre above for why the env is activated here and not
        # via Snakemake's `conda:` directive.
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
        cell_meta=RESULTS + "/splitpipe/{sublib}/" + config["splitpipe"]["cell_meta"],
    output:
        barcodes=RESULTS + "/splitpipe/{sublib}/barcodes.tsv",
    params:
        col=config["splitpipe"]["barcode_column"],
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/extract_barcodes.py"


# Index the BAM if split-pipe did not.
rule index_bam:
    input:
        bam=RESULTS + "/splitpipe/{sublib}/" + config["splitpipe"]["bam_name"],
    output:
        bai=RESULTS + "/splitpipe/{sublib}/" + config["splitpipe"]["bam_name"] + ".bai",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "samtools index {input.bam}"
