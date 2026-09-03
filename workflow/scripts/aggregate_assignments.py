"""Merge per-sublibrary Vireo `donor_ids.tsv` files into one table and a summary.

Vireo's donor_ids.tsv columns: cell, donor_id, prob_max, prob_doublet,
n_vars, best_singlet, best_doublet, doublet_logLikRatio (version-dependent).
`donor_id` is 'donorNN' / the VCF sample name, or 'doublet' / 'unassigned'.
"""
import pandas as pd

assign_paths = snakemake.input.assigns          # noqa: F821
meta_in = list(snakemake.input.donor_meta)      # noqa: F821  (empty = none)
sublibs = list(snakemake.params.sublibs)        # noqa: F821
out_table = snakemake.output.table              # noqa: F821
out_summary = snakemake.output.summary          # noqa: F821

frames = []
for sublib, path in zip(sublibs, assign_paths):
    df = pd.read_csv(path, sep="\t")
    df.insert(0, "sublibrary", sublib)
    # globally unique cell id across sublibraries
    df.insert(1, "cell_id", sublib + "|" + df["cell"].astype(str))
    frames.append(df)

allcells = pd.concat(frames, ignore_index=True)

# attach donor metadata (path_diag, mutation, ...) where donor_id matches
if meta_in:
    donors = pd.read_csv(meta_in[0], sep="\t", comment="#", dtype=str)
    meta_cols = [c for c in donors.columns if c != "vcf"]
    allcells = allcells.merge(donors[meta_cols], how="left", on="donor_id")

allcells.to_csv(out_table, sep="\t", index=False)

# summary: cells per donor per sublibrary (incl. doublet / unassigned)
summary = (
    allcells.groupby(["sublibrary", "donor_id"], dropna=False)
    .size()
    .reset_index(name="n_cells")
    .pivot(index="donor_id", columns="sublibrary", values="n_cells")
    .fillna(0)
    .astype(int)
)
summary["total"] = summary.sum(axis=1)
summary = summary.sort_values("total", ascending=False)
summary.to_csv(out_summary, sep="\t")
print(f"[aggregate] {len(allcells)} cells across {len(sublibs)} sublibraries")
