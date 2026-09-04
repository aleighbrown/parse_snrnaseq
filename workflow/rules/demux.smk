# -----------------------------------------------------------------------------
# Genotype demultiplexing, per sublibrary:
#   cellsnp-lite : pileup each cell at the reference SNP sites (region mode -R)
#   vireo        : assign each cell to a NAMED donor using the donor genotypes
#
# Both tools are in the Demuxafy container, or install via conda (envs/demux.yaml).
# -----------------------------------------------------------------------------

rule cellsnp:
    input:
        bam=RESULTS + "/splitpipe/{sublib}/" + config["splitpipe"]["bam_name"],
        bai=RESULTS + "/splitpipe/{sublib}/" + config["splitpipe"]["bam_name"] + ".bai",
        barcodes=RESULTS + "/splitpipe/{sublib}/barcodes.tsv",
        ref=REF_VCF,
        ref_tbi=REF_VCF + ".tbi",
    output:
        # cellsnp-lite writes a directory of outputs; this is its base VCF.
        vcf=RESULTS + "/cellsnp/{sublib}/cellSNP.base.vcf.gz",
    params:
        outdir=RESULTS + "/cellsnp/{sublib}",
        min_maf=config["cellsnp"]["min_maf"],
        min_count=config["cellsnp"]["min_count"],
        cell_tag=config["cell_tag"],
        umi_tag=config["umi_tag"],
    threads: nthreads("cellsnp")
    conda:
        "../envs/demux.yaml"
    shell:
        r"""
        set -euo pipefail
        cellsnp-lite \
            -s {input.bam} \
            -b {input.barcodes} \
            -O {params.outdir} \
            -R {input.ref} \
            -p {threads} \
            --minMAF {params.min_maf} \
            --minCOUNT {params.min_count} \
            --cellTAG {params.cell_tag} \
            --UMItag {params.umi_tag} \
            --gzip
        """


rule vireo:
    input:
        cellsnp=RESULTS + "/cellsnp/{sublib}/cellSNP.base.vcf.gz",
        ref=REF_VCF,
    output:
        assign=RESULTS + "/vireo/{sublib}/donor_ids.tsv",
        summary=RESULTS + "/vireo/{sublib}/summary.tsv",
    params:
        indir=RESULTS + "/cellsnp/{sublib}",
        outdir=RESULTS + "/vireo/{sublib}",
        gt_tag=config["vireo"]["genotype_tag"],
    threads: nthreads("vireo")
    conda:
        "../envs/demux.yaml"
    shell:
        r"""
        set -euo pipefail
        # With donor genotypes (-d), the donor count comes from the VCF samples.
        vireo \
            -c {params.indir} \
            -d {input.ref} \
            -t {params.gt_tag} \
            -o {params.outdir} \
            -p {threads}
        """
