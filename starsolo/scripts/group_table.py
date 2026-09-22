"""Map every cell barcode to a donor x sort group, across all sublibraries.

Input is the annotation table: one row per cell, with a barcode column (usually
`bc_wells`, with or without the `__s<sublibrary>` suffix), a donor column and a
sort/condition column. Only cells present in that table are grouped, so any
filtering done there (QC, cell typing, dropping controls) carries through.

Runs once for all sublibraries, so `--min-cells` counts a donor x sort group
over the whole experiment rather than one sublibrary at a time.

Writes:
  <sublib>/cell_groups.tsv   headerless cb<TAB>group, for `sinto filterbarcodes`
  groups.tsv                 one row per kept group: group, donor, sort, cells,
                             sublibraries
  cells_per_group.tsv        donor x sort cell counts, including dropped groups
"""
import argparse
import os
import re
import sys

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--cells", nargs="+", required=True,
                help="whitelist/<sublib>/cells.tsv, one per sublibrary")
ap.add_argument("--sublibraries", nargs="+", required=True,
                help="sublibrary names, in the same order as --cells")
ap.add_argument("--annotation", required=True)
ap.add_argument("--donor-column", default="donor_id")
ap.add_argument("--barcode-column", default="bc_wells")
ap.add_argument("--condition-column", default="sample")
ap.add_argument("--sort-column", default="sample", help="the sort column in cells.tsv")
ap.add_argument("--min-cells", type=int, default=0,
                help="drop donor x sort groups with fewer cells than this, in total")
ap.add_argument("--outdir", required=True)
ap.add_argument("--groups", required=True)
ap.add_argument("--summary", required=True)
a = ap.parse_args()

if len(a.cells) != len(a.sublibraries):
    sys.exit("--cells and --sublibraries must have the same length")

sep = "," if a.annotation.endswith((".csv", ".csv.gz")) else "\t"
ann = pd.read_csv(a.annotation, sep=sep, dtype=str)
for col in (a.donor_column, a.barcode_column, a.condition_column):
    if col not in ann.columns:
        sys.exit(f"'{col}' not in {a.annotation}: {list(ann.columns)}")

bc = ann[a.barcode_column].astype(str)
if bc.str.contains("__s", regex=False).any():
    parts = bc.str.split("__s", n=1, expand=True)
    ann["bc_wells"], ann["_suffix"] = parts[0], parts[1]
else:
    sys.exit(f"{a.annotation}: '{a.barcode_column}' has no __s<sublibrary> suffix, so "
             "rows cannot be assigned to a sublibrary")

ann = ann[[a.donor_column, a.condition_column, "bc_wells", "_suffix"]].copy()
dups = int(ann.duplicated(["bc_wells", "_suffix"]).sum())
if dups:
    sys.exit(f"{a.annotation}: {dups} duplicated barcodes -- one row per cell is expected")


def safe(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "-", str(s)).strip("-")


frames = []
for path, sublib in zip(a.cells, a.sublibraries):
    cells = pd.read_csv(path, sep="\t", dtype=str)
    tag = None
    if "barcode_original" in cells.columns and cells["barcode_original"].iloc[0].count("__s"):
        tag = cells["barcode_original"].iloc[0].split("__s", 1)[1]
    mine = ann[ann["_suffix"].isin({sublib} | ({tag} if tag else set()))]
    if mine.empty:
        sys.exit(f"{a.annotation}: no cells with suffix for {sublib} "
                 f"(looked for {sorted({sublib} | ({tag} if tag else set()))})")
    # split-pipe's sort column is renamed: it usually has the same name as the
    # annotation's condition column, and pandas would suffix both.
    have_sort = a.sort_column in cells.columns
    keep_cols = ["cb", "bc_wells"] + ([a.sort_column] if have_sort else [])
    m = (mine.merge(cells[keep_cols].rename(columns={a.sort_column: "_sp_sort"}),
                    on="bc_wells", how="left"))
    missing = int(m["cb"].isna().sum())
    if missing:
        # The whitelist is built from the same annotation table, so this means
        # the two runs disagree; grouping cells we never aligned is pointless.
        sys.exit(f"{sublib}: {missing} annotated cells are not in {path}. "
                 "Is cells.subset in config.yaml the same file as split.annotation?")
    if have_sort:
        bad = int((m[a.condition_column] != m["_sp_sort"]).sum())
        if bad:
            print(f"[groups] WARNING: {sublib}: {bad} cells where "
                  f"'{a.condition_column}' disagrees with split-pipe's "
                  f"'{a.sort_column}'; using split-pipe's", file=sys.stderr)
        m["sort"] = m["_sp_sort"]
    else:
        m["sort"] = m[a.condition_column]
    m["sublibrary"] = sublib
    m["group"] = m[a.donor_column].map(safe) + "__" + m["sort"].map(safe)
    frames.append(m[["cb", "group", "sublibrary", a.donor_column, "sort"]])

all_cells = pd.concat(frames, ignore_index=True)

counts = (all_cells.groupby(["group", a.donor_column, "sort"])
          .agg(cells=("cb", "size"), sublibraries=("sublibrary", "nunique"))
          .reset_index())
counts["kept"] = counts["cells"] >= a.min_cells
counts.sort_values("cells", ascending=False).to_csv(a.summary, sep="\t", index=False)

dropped = counts[~counts["kept"]]
if len(dropped):
    print(f"[groups] {len(dropped)} groups below --min-cells {a.min_cells}, dropped:\n"
          + dropped[["group", "cells"]].to_string(index=False), file=sys.stderr)

keep = set(counts.loc[counts["kept"], "group"])
all_cells = all_cells[all_cells["group"].isin(keep)]

counts[counts["kept"]].drop(columns="kept").sort_values("group").to_csv(
    a.groups, sep="\t", index=False)

for sublib, part in all_cells.groupby("sublibrary"):
    out = os.path.join(a.outdir, sublib, "cell_groups.tsv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    part[["cb", "group"]].to_csv(out, sep="\t", index=False, header=False)
    print(f"[groups] {sublib}: {len(part)} cells in {part['group'].nunique()} groups")

# A sublibrary with no kept cells still needs its (empty) file, or the split
# rule has no input.
for sublib in a.sublibraries:
    out = os.path.join(a.outdir, sublib, "cell_groups.tsv")
    if not os.path.exists(out):
        os.makedirs(os.path.dirname(out), exist_ok=True)
        open(out, "w").close()
        print(f"[groups] WARNING: {sublib}: no cells in any kept group", file=sys.stderr)

print(f"[groups] {len(all_cells)} cells, {len(keep)} donor x sort groups "
      f"(min {a.min_cells} cells)")
