"""Import an existing per-barcode donor assignment table.

Used when `demux.mode: provided` -- i.e. the demultiplexing was already done
elsewhere (here: the collaborator's demuxalot + Vireo result) and we only need
it normalised into the same shape the Vireo branch produces, so everything
downstream is identical either way.

Input barcodes look like "01_01_67__sParseFANS1", i.e.
<bc_wells>__s<sublibrary>. The part before "__s" is the per-sublibrary barcode
that matches split-pipe's `bc_wells`; the part after is the sublibrary.
"""
import sys

import pandas as pd

cfg = snakemake.params.provided            # noqa: F821
in_path = snakemake.input.table            # noqa: F821
meta_in = list(snakemake.input.donor_meta)  # noqa: F821  (empty = no metadata)
sublibs = list(snakemake.params.sublibs)   # noqa: F821
out_table = snakemake.output.table         # noqa: F821
out_summary = snakemake.output.summary     # noqa: F821
out_manifest = snakemake.output.manifest   # noqa: F821

bc_col = cfg["barcode_column"]
donor_col = cfg["donor_column"]
sublib_col = cfg.get("sublibrary_column")

df = pd.read_csv(in_path, sep="\t", dtype={bc_col: str, donor_col: str})

for col in (bc_col, donor_col):
    if col not in df.columns:
        raise KeyError(
            f"Column '{col}' not in {in_path}. Available: {list(df.columns)}. "
            "Fix demux.provided.* in config/config.yaml."
        )

# --- split "<bc_wells>__s<sublibrary>" ---------------------------------------
split = df[bc_col].str.split("__s", n=1, expand=True)
if split.shape[1] != 2 or split[1].isna().any():
    raise ValueError(
        f"Barcodes in '{bc_col}' are not of the form <barcode>__s<sublibrary>; "
        f"first value seen: {df[bc_col].iloc[0]!r}"
    )
df["barcode"] = split[0]
df["sublibrary"] = split[1]

# Cross-check against the table's own sublibrary column, if there is one.
if sublib_col and sublib_col in df.columns:
    bad = int((df["sublibrary"] != df[sublib_col].astype(str)).sum())
    if bad:
        raise ValueError(
            f"{bad} rows disagree between the barcode suffix and '{sublib_col}'."
        )

unknown = sorted(set(df["sublibrary"]) - set(sublibs))
if unknown:
    raise ValueError(
        f"Sublibraries in {in_path} that are not in config/sublibraries.tsv: "
        f"{unknown}"
    )
missing = sorted(set(sublibs) - set(df["sublibrary"]))
if missing:
    print(
        f"[import_assignments] WARNING: no assigned cells for {missing}",
        file=sys.stderr,
    )

# --- optionally drop doublets -------------------------------------------------
if cfg.get("drop_doublets") and "dmux_type" in df.columns:
    before = len(df)
    df = df[df["dmux_type"] != "DBL"]
    print(f"[import_assignments] dropped {before - len(df)} doublet rows")

# --- normalise to the Vireo-branch schema ------------------------------------
df = df.rename(columns={donor_col: "donor_id"})
df.insert(0, "cell_id", df["sublibrary"] + "|" + df["barcode"])

keep = ["cell_id", "sublibrary", "barcode", "donor_id"]
keep += [c for c in cfg.get("keep_columns", []) if c in df.columns and c not in keep]
skipped = [c for c in cfg.get("keep_columns", []) if c not in df.columns]
if skipped:
    print(f"[import_assignments] keep_columns not in input, skipped: {skipped}")

out = df[keep].copy()

# --- optional donor metadata --------------------------------------------------
# Only the columns a human maintains are joined; the count columns in
# donor_metadata.tsv are derived and would go stale, so they are recomputed into
# the manifest below instead of being carried onto every cell.
_DERIVED = {
    "n_cells", "n_sublibraries", "high_conf_donor", "only_demuxalot",
    "likely_missing_donor", "K_negative", "TDP_negative", "positive",
}
if meta_in:
    meta = pd.read_csv(meta_in[0], sep="\t", comment="#", dtype=str)
    if "donor_id" not in meta.columns:
        raise KeyError(f"{meta_in[0]} has no 'donor_id' column.")
    join_cols = ["donor_id"] + [
        c for c in meta.columns if c != "donor_id" and c not in _DERIVED
    ]
    overlap = len(set(out["donor_id"]) & set(meta["donor_id"]))
    if not overlap:
        print(
            f"[import_assignments] WARNING: no donor_id overlap with "
            f"{meta_in[0]} -- metadata not joined",
            file=sys.stderr,
        )
    else:
        out = out.merge(meta[join_cols], how="left", on="donor_id")
        missing = sorted(set(out["donor_id"]) - set(meta["donor_id"]))
        if missing:
            print(
                f"[import_assignments] WARNING: no metadata row for {missing}",
                file=sys.stderr,
            )
        print(f"[import_assignments] joined metadata for {overlap} donors")
else:
    print("[import_assignments] donor_metadata not set -- no metadata joined")

out.to_csv(out_table, sep="\t", index=False)

# --- summary: cells per donor per sublibrary ---------------------------------
summary = (
    out.groupby(["sublibrary", "donor_id"], dropna=False)
    .size()
    .reset_index(name="n_cells")
    .pivot(index="donor_id", columns="sublibrary", values="n_cells")
    .reindex(columns=[s for s in sublibs if s in set(out["sublibrary"])])
    .fillna(0)
    .astype(int)
)
summary["total"] = summary.sum(axis=1)
summary = summary.sort_values("total", ascending=False)
summary.to_csv(out_summary, sep="\t")

# --- donor manifest: one row per donor, all counts recomputed from the data ---
man = pd.DataFrame({"donor_id": sorted(out["donor_id"].unique())})
grp = out.groupby("donor_id")
man["n_cells"] = man["donor_id"].map(grp.size())
man["n_sublibraries"] = man["donor_id"].map(grp["sublibrary"].nunique())
for col, values in (
    ("demux_group", ["high_conf_donor", "only_demuxalot", "likely_missing_donor"]),
    ("condition", ["K_negative", "TDP_negative", "positive"]),
):
    if col in out.columns:
        counts = out.pivot_table(
            index="donor_id", columns=col, values="cell_id", aggfunc="count"
        ).fillna(0).astype(int)
        for v in values:
            man[v] = man["donor_id"].map(counts[v]) if v in counts.columns else 0
# A donorNN id means Vireo reconstructed the donor but never matched it to a
# real identity -- there is no genotype for that person in the panel.
man["source"] = [
    "vireo_virtual" if str(d).startswith("donor") else "genotyped"
    for d in man["donor_id"]
]
man = man.sort_values("n_cells", ascending=False)
man.to_csv(out_manifest, sep="\t", index=False)

print(
    f"[import_assignments] {len(out)} cells, "
    f"{out['donor_id'].nunique()} donors, "
    f"{out['sublibrary'].nunique()} sublibraries -> {out_table}"
)
