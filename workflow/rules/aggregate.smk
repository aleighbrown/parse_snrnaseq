# -----------------------------------------------------------------------------
# Produce the final per-cell donor assignment table + per-donor summary.
#
# demux.mode == "vireo"     combine the per-sublibrary Vireo calls
# demux.mode == "provided"  import an assignment table produced elsewhere
#                           (the collaborator's demuxalot + Vireo result)
#
# Both branches write the same two files, so everything downstream is identical.
# -----------------------------------------------------------------------------

# DONOR_META is optional: [] means "no metadata file, skip the join".
_meta_in = [DONOR_META] if DONOR_META else []


if DEMUX_MODE == "vireo":

    rule aggregate_assignments:
        input:
            assigns=expand(RESULTS + "/vireo/{sublib}/donor_ids.tsv", sublib=SUBLIBS),
            donor_meta=_meta_in,
        output:
            table=RESULTS + "/demux/donor_assignments.all_sublibraries.tsv",
            summary=RESULTS + "/demux/demux_summary.tsv",
        params:
            sublibs=SUBLIBS,
        conda:
            "../envs/py.yaml"
        script:
            "../scripts/aggregate_assignments.py"

else:

    rule import_assignments:
        input:
            table=config["demux"]["provided"]["table"],
            donor_meta=_meta_in,
        output:
            table=RESULTS + "/demux/donor_assignments.all_sublibraries.tsv",
            summary=RESULTS + "/demux/demux_summary.tsv",
            manifest=RESULTS + "/demux/donor_manifest.tsv",
        params:
            provided=config["demux"]["provided"],
            sublibs=SUBLIBS,
        log:
            "logs/demux/import_assignments.log",
        conda:
            "../envs/py.yaml"
        script:
            "../scripts/import_assignments.py"


# Regenerate the derived donor table: `snakemake donor_manifest -c1`
rule donor_manifest:
    input:
        RESULTS + "/demux/donor_manifest.tsv",
