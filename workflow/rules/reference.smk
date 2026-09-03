# -----------------------------------------------------------------------------
# Genome reference: build the split-pipe reference with `split-pipe --mode mkref`
# as a REUSABLE, SHARED resource under genome.ref_root -- built once, used by
# any project that needs a Parse human reference.
#
# Each entry in `genome.assemblies` gives the FASTA and GTF either as a LOCAL
# PATH (`fasta:` / `gtf:`) or as a URL to download (`fasta_url:` / `gtf_url:`).
#
#   local + uncompressed -> read in place, wrapped in ancient() so its mtime can
#                           never trigger a rebuild
#   local + .gz          -> decompressed once into <ref_root>/_sources/ (mkref
#                           needs plain text), also ancient() thereafter
#   url                  -> downloaded + decompressed into <ref_root>/_sources/
#
# Layout under ref_root:
#   <ref_root>/<TAG>/                  the split-pipe reference itself
#   <ref_root>/<TAG>.mkref.done        provenance: what it was built from, when
#   <ref_root>/_sources/<assembly>/    decompressed FASTA/GTF, shared too
#
# Once <ref_root>/<TAG>/ exists, Snakemake treats it as up to date and will not
# rebuild it -- every input is ancient(), so no timestamp can invalidate it.
# To force a rebuild, delete the directory.
#
# Two or more entries build ONE mixed-species reference: split-pipe mkref takes
# parallel --genome_name / --fasta / --genes lists.
#
# NOTE: split-pipe is Parse's proprietary binary and must be on PATH -- not
# conda-installable, so no `conda:` directive. A human reference needs a lot of
# RAM (~64 GB) and an hour or two.
# -----------------------------------------------------------------------------

_gcfg = config["genome"]
ASSEMBLIES = {a["name"]: a for a in _gcfg["assemblies"]}
# Reference dir is named after the assemblies it contains ("GRCh38_gencode_v42",
# "hg38_mm10"), so builds never collide within the shared tree.
GENOME_TAG = "_".join(ASSEMBLIES)

# Shared reference tree. Project-local default if unset.
REF_ROOT = (_gcfg.get("ref_root") or "resources/genomes").rstrip("/")

# local-path key -> url key, per file kind that mkref consumes
_SOURCE_KEYS = {"genome.fa": ("fasta", "fasta_url"), "genes.gtf": ("gtf", "gtf_url")}


def _source(assembly, kind):
    """Return (local_path_or_None, url_or_None) for one assembly's FASTA/GTF."""
    spec = ASSEMBLIES[assembly]
    local_key, url_key = _SOURCE_KEYS[kind]
    local, url = spec.get(local_key), spec.get(url_key)
    if not local and not url:
        raise WorkflowError(
            f"genome.assemblies['{assembly}'] needs either '{local_key}' "
            f"(a path) or '{url_key}' (a URL)."
        )
    return local, url


def genome_input(assembly, kind):
    """Path mkref should read. Uncompressed local files are used where they lie;
    anything else is staged into the shared _sources tree. ancient() throughout,
    so no source timestamp can ever trigger a rebuild of the reference."""
    local, _ = _source(assembly, kind)
    if local and not local.endswith(".gz"):
        return ancient(local)
    return ancient(f"{REF_ROOT}/_sources/{assembly}/{kind}")


if _gcfg.get("prebuilt_dir"):
    # Already-built split-pipe reference: use as-is, build nothing.
    GENOME_DIR = _gcfg["prebuilt_dir"]

else:
    GENOME_DIR = f"{REF_ROOT}/{GENOME_TAG}"

    wildcard_constraints:
        assembly="|".join(ASSEMBLIES),
        kind="genome\\.fa|genes\\.gtf",

    rule stage_genome_file:
        """Decompress a local .gz (or download) into the shared _sources tree."""
        output:
            f"{REF_ROOT}/_sources/{{assembly}}/{{kind}}",
        params:
            local=lambda wc: _source(wc.assembly, wc.kind)[0] or "",
            url=lambda wc: _source(wc.assembly, wc.kind)[1] or "",
        log:
            "logs/reference/stage_{assembly}_{kind}.log",
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output})"
            if [ -n "{params.local}" ]; then
                gzip -cd "{params.local}" > {output} 2> {log}
            else
                curl -fsSL "{params.url}" 2> {log} | gzip -cdf > {output}
            fi
            """

    rule mkref:
        input:
            fastas=[genome_input(a, "genome.fa") for a in ASSEMBLIES],
            gtfs=[genome_input(a, "genes.gtf") for a in ASSEMBLIES],
        output:
            genome_dir=directory(GENOME_DIR),
            done=f"{REF_ROOT}/{GENOME_TAG}.mkref.done",
        params:
            names=" ".join(ASSEMBLIES),
            extra=_gcfg.get("mkref_extra", ""),
        threads: nthreads("mkref")
        resources:
            mem_mb=_gcfg.get("mkref_mem_mb", 64000),
            runtime=_gcfg.get("mkref_runtime", 480),   # minutes
        log:
            "logs/reference/mkref.log",
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.genome_dir})"
            split-pipe \
                --mode mkref \
                --genome_name {params.names} \
                --fasta {input.fastas} \
                --genes {input.gtfs} \
                --output_dir {output.genome_dir} \
                --nthreads {threads} \
                {params.extra} > {log} 2>&1

            # Provenance stamp: this reference is shared, so record what built it.
            {{
                echo "built:       $(date -Iseconds)"
                echo "by:          ${{USER:-unknown}}@$(hostname)"
                echo "split-pipe:  $(split-pipe --version 2>&1 | tr '\n' ' ')"
                echo "genome_name: {params.names}"
                for f in {input.fastas}; do echo "fasta:       $f"; done
                for g in {input.gtfs};   do echo "genes:       $g"; done
                echo "extra:       {params.extra}"
            }} > {output.done}
            """


# Convenience target: `snakemake reference -c 16` builds only the reference.
rule reference:
    input:
        GENOME_DIR,
