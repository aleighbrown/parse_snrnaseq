"""Check the extraction against split-pipe: reads we kept per cell vs
split-pipe's mread_count for the same cell.

If barcode handling matches split-pipe, the two track tightly at a roughly
constant ratio (ours counts all reads, split-pipe counts mapped ones). Things
this catches:
  * wrong sublibrary / cells_tag        -> weak correlation everywhere
  * round-1 T/R not merged correctly    -> odd reads_R fraction in some wells
  * a well-specific barcode problem     -> one well's ratio far from the rest

Writes a per-cell table, a per-round-1-well table and a short verdict.
"""
import argparse
import json

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True)
ap.add_argument("--cells", required=True)
ap.add_argument("--stats", required=True)
ap.add_argument("--read-count-column", default="mread_count")
ap.add_argument("--sort-column", default="sample")
ap.add_argument("--out-cells", required=True)
ap.add_argument("--out-wells", required=True)
ap.add_argument("--out-summary", required=True)
a = ap.parse_args()

counts = pd.read_csv(a.counts, sep="\t")
cells = pd.read_csv(a.cells, sep="\t")
df = cells.merge(counts[["cb", "reads_T", "reads_R", "reads"]], on="cb", how="left")
ref = a.read_count_column
df["ratio"] = df["reads"] / df[ref].replace(0, np.nan)
df.to_csv(a.out_cells, sep="\t", index=False)

group_cols = ["w1"] + ([a.sort_column] if a.sort_column in df.columns else [])
wells = df.groupby(group_cols).agg(
    cells=("cb", "size"),
    reads=("reads", "sum"),
    reads_R=("reads_R", "sum"),
    splitpipe=(ref, "sum"),
).reset_index()
wells["ratio"] = wells["reads"] / wells["splitpipe"]
wells["frac_R"] = wells["reads_R"] / wells["reads"]
med_ratio = wells["ratio"].median()
med_r = wells["frac_R"].median()
wells["flag"] = np.where(
    (wells["ratio"] / med_ratio).between(0.5, 2) & ((wells["frac_R"] - med_r).abs() < 0.2),
    "", "CHECK")
wells.to_csv(a.out_wells, sep="\t", index=False)

with open(a.stats) as fh:
    stats = json.load(fh)

rho = df["reads"].rank().corr(df[ref].rank())      # Spearman, without scipy
r_log = np.log1p(df["reads"]).corr(np.log1p(df[ref]))
with_reads = (df["reads"] > 0).mean()
flagged = int((wells["flag"] == "CHECK").sum())
ok = rho >= 0.9 and with_reads >= 0.99 and flagged == 0

lines = [
    f"verdict: {'PASS' if ok else 'CHECK'}",
    f"cells in whitelist:                 {len(df)}",
    f"cells with >=1 extracted read:      {with_reads:.2%}",
    f"spearman(reads, {ref}):     {rho:.3f}",
    f"pearson(log reads, log {ref}): {r_log:.3f}",
    f"total ratio reads / {ref}:  {df['reads'].sum() / df[ref].sum():.3f}",
    f"median per-cell ratio:              {df['ratio'].median():.3f}",
    f"random-hexamer fraction (median well): {med_r:.3f}",
    f"round-1 wells flagged CHECK:        {flagged} of {len(wells)}",
    "",
    "extraction:",
] + [f"  {k}: {v}" for k, v in stats.items()]
if a.sort_column in df.columns:
    by_sort = df.groupby(a.sort_column).agg(cells=("cb", "size"), reads=("reads", "sum"),
                                           splitpipe=(ref, "sum"))
    by_sort["ratio"] = (by_sort["reads"] / by_sort["splitpipe"]).round(3)
    lines += ["", "by sort:", by_sort.to_string()]
if flagged:
    lines += ["", "flagged wells:", wells[wells["flag"] == "CHECK"].to_string(index=False)]

text = "\n".join(lines) + "\n"
with open(a.out_summary, "w") as fh:
    fh.write(text)
print(text)
