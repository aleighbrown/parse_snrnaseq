"""Turn split-pipe's cell list into a STARsolo whitelist for one sublibrary.

Input is split-pipe's cell_metadata.csv. The combined file suffixes barcodes
with the sublibrary ("01_01_67__s1"); a per-sublibrary file has no suffix.
Rows for this sublibrary are kept when the suffix equals either the sublibrary's
`cells_tag` (e.g. "1") or its name (e.g. "ParseFANS1").

Each cell's 24-nt barcode is bc1(T) + bc2 + bc3 sequence. Writes:
  cells.tsv      cb, bc_wells, sublibrary, w1, w2, w3, + split-pipe's columns
  whitelist.txt  one cb per line, for --soloCBwhitelist
"""
import argparse
import sys

import pandas as pd

from parse_barcodes import describe_table, load_table, well_name, well_sequences

ap = argparse.ArgumentParser()
ap.add_argument("--cells-file", required=True)
ap.add_argument("--barcode-column", required=True)
ap.add_argument("--sublibrary", required=True)
ap.add_argument("--tag", default="")
ap.add_argument("--bc1", required=True)
ap.add_argument("--bc2", required=True)
ap.add_argument("--bc3", required=True)
ap.add_argument("--out-cells", required=True)
ap.add_argument("--out-whitelist", required=True)
a = ap.parse_args()

sep = "," if a.cells_file.endswith((".csv", ".csv.gz")) else "\t"
df = pd.read_csv(a.cells_file, sep=sep, dtype={a.barcode_column: str})
if a.barcode_column not in df.columns:
    sys.exit(f"Column '{a.barcode_column}' not in {a.cells_file}: {list(df.columns)}")

bc = df[a.barcode_column].astype(str)
suffixed = bc.str.contains("__s", regex=False)
if suffixed.any():
    if not suffixed.all():
        sys.exit(f"{a.cells_file}: some barcodes carry a __s<sublibrary> suffix and some do not")
    parts = bc.str.split("__s", n=1, expand=True)
    wanted = {a.sublibrary} | ({a.tag} if a.tag else set())
    keep = parts[1].isin(wanted)
    if not keep.any():
        seen = sorted(parts[1].unique())[:20]
        sys.exit(f"No rows with suffix in {sorted(wanted)}; suffixes present: {seen}. "
                 "Set cells_tag for this sublibrary in config/sublibraries.tsv.")
    df = df[keep].copy()
    df["barcode_original"] = bc[keep]
    df["bc_wells"] = parts[0][keep]
else:
    print(f"[whitelist] no __s suffix: treating every row as {a.sublibrary}")
    df = df.copy()
    df["barcode_original"] = bc
    df["bc_wells"] = bc

wells = df["bc_wells"].str.split("_", expand=True)
if wells.shape[1] != 3:
    sys.exit(f"bc_wells not of the form NN_NN_NN, e.g. {df['bc_wells'].iloc[0]!r}")
for i, col in enumerate(("w1", "w2", "w3")):
    df[col] = wells[i].astype(int)
    bad = ~df[col].between(1, 96)
    if bad.any():
        sys.exit(f"{int(bad.sum())} cells have {col} outside 1..96")

# Cross-check the well numbering against split-pipe's own well-name columns.
for i in (1, 2, 3):
    col = f"bc{i}_well"
    if col in df.columns:
        n_bad = int((df[f"w{i}"].map(well_name) != df[col].astype(str)).sum())
        if n_bad:
            sys.exit(f"{n_bad} cells disagree between bc_wells and {col}: "
                     "well numbering is not 8x12 row-major")

seq = {}
for i, path in ((1, a.bc1), (2, a.bc2), (3, a.bc3)):
    table, _ = load_table(path, round1=(i == 1))
    print(f"[whitelist] round {i}: " + describe_table(table, path))
    seq[i] = well_sequences(table, path)
    missing = sorted(set(df[f"w{i}"]) - set(seq[i]))
    if missing:
        sys.exit(f"{path}: no sequence for round-{i} wells "
                 f"{[well_name(w) for w in missing][:10]}")

df["cb"] = (df["w1"].map(seq[1]) + df["w2"].map(seq[2]) + df["w3"].map(seq[3]))
df["sublibrary"] = a.sublibrary

for col in ("bc_wells", "cb"):
    dup = int(df[col].duplicated().sum())
    if dup:
        sys.exit(f"{dup} duplicated {col} values for {a.sublibrary}")

front = ["cb", "bc_wells", "sublibrary", "w1", "w2", "w3"]
drop = {a.barcode_column} - {"bc_wells"}
rest = [c for c in df.columns if c not in front and c not in drop]
df[front + rest].to_csv(a.out_cells, sep="\t", index=False)
df["cb"].to_csv(a.out_whitelist, index=False, header=False)

print(f"[whitelist] {a.sublibrary}: {len(df)} cells, cb length {len(df['cb'].iloc[0])}")
if "sample" in df.columns:
    print(df["sample"].value_counts().to_string())
