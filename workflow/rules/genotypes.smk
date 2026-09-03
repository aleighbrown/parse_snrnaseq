# -----------------------------------------------------------------------------
# Build the merged, normalised multi-sample GENOTYPE PANEL used by cellsnp+vireo.
# (Not to be confused with the GENOME reference built in reference.smk.)
# Merges the per-donor VCFs listed in the `donor_vcfs` TSV, keeps biallelic SNPs,
# and (optionally) harmonises contig names to match the split-pipe BAM.
# Skip entirely by setting build_reference_vcf: false and pointing
# reference_vcf at an existing joint multi-sample VCF. NOTE: no VCFs exist for
# this cohort, so this file is currently unreachable (demux.mode: provided).
# -----------------------------------------------------------------------------

REF_VCF = config["reference_vcf"]

if config.get("build_reference_vcf", True):
    # Per-donor VCFs to merge, from the donor_vcfs TSV (donor_id, vcf).
    _dv = pd.read_csv(config["donor_vcfs"], sep="\t", comment="#")
    DONOR_VCFS = _dv["vcf"].tolist()


    rule prepare_reference_vcf:
        input:
            vcfs=DONOR_VCFS,
        output:
            vcf=REF_VCF,
            tbi=REF_VCF + ".tbi",
        params:
            rename=config.get("rename_chrs") or "",
        threads: nthreads("prepare_reference_vcf")
        conda:
            "../envs/bcftools.yaml"
        shell:
            r"""
            set -euo pipefail
            # index any per-donor VCFs that lack a tabix index
            for v in {input.vcfs}; do
                [ -f "$v.tbi" ] || [ -f "$v.csi" ] || bcftools index -t "$v"
            done

            # merge -> keep biallelic SNPs only -> left-align/normalise
            bcftools merge --threads {threads} -m none -Ou {input.vcfs} \
              | bcftools view -m2 -M2 -v snps -Ou \
              | bcftools norm -m -any -Oz -o {output.vcf}.tmp.vcf.gz

            if [ -n "{params.rename}" ]; then
                bcftools annotate --rename-chrs {params.rename} \
                    -Oz -o {output.vcf} {output.vcf}.tmp.vcf.gz
                rm -f {output.vcf}.tmp.vcf.gz
            else
                mv {output.vcf}.tmp.vcf.gz {output.vcf}
            fi
            bcftools index -t {output.vcf}
            """
