"""Map each STARsolo cell barcode to a donor x sort group, for sinto.

Donor comes from the donor assignment table (the main pipeline's
results/demux/donor_assignments.all_sublibraries.tsv); sort comes from
split-pipe's own `sample` column in cells.tsv, cross-checked against the
table's `condition` column when present. Cells without a donor call are left
out. Group names are made filesystem-safe (pathology numbers contain '/').

Writes a headerless `cb<TAB>group` file and a cells-per-group summary.
"""
import argparse
import re
import sys

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--cells", required=True)
ap.add_argument("--donor-table", required=True)
ap.add_argument("--sublibrary", required=True)
ap.add_argument("--donor-column", default="donor_id")
ap.add_argument("--barcode-column", default="barcode")
ap.add_argument("--sublibrary-column", default="sublibrary")
ap.add_argument("--condition-column", default="condition")
ap.add_argument("--sort-column", default="sample")
ap.add_argument("--out", required=True)
ap.add_argument("--summary", required=True)
a = ap.parse_args()

cells = pd.read_csv(a.cells, sep="\t", dtype=str)
donors = pd.read_csv(a.donor_table, sep="\t", dtype=str)
for col in (a.donor_column, a.barcode_column, a.sublibrary_column):
    if col not in donors.columns:
        sys.exit(f"'{col}' not in {a.donor_table}: {list(donors.columns)}")
donors = donors[donors[a.sublibrary_column] == a.sublibrary]
if donors.empty:
    sys.exit(f"No rows for {a.sublibrary} in {a.donor_table}")

m = donors.merge(cells[["cb", "bc_wells", a.sort_column]],
                 left_on=a.barcode_column, right_on="bc_wells", how="left")
missing = int(m["cb"].isna().sum())
if missing:
    print(f"[groups] WARNING: {missing} donor-assigned barcodes are not in the "
          "split-pipe cell list and are skipped", file=sys.stderr)
    m = m[m["cb"].notna()]

if a.condition_column in m.columns:
    bad = int((m[a.condition_column] != m[a.sort_column]).sum())
    if bad:
        print(f"[groups] WARNING: {bad} cells where '{a.condition_column}' "
              f"disagrees with split-pipe's '{a.sort_column}'; using the latter",
              file=sys.stderr)


def safe(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "-", str(s)).strip("-")


m["group"] = m[a.donor_column].map(safe) + "__" + m[a.sort_column].map(safe)
m[["cb", "group"]].to_csv(a.out, sep="\t", index=False, header=False)

summary = (m.groupby([a.donor_column, a.sort_column]).size()
           .unstack(fill_value=0))
summary["total"] = summary.sum(axis=1)
summary.sort_values("total", ascending=False).to_csv(a.summary, sep="\t")
print(f"[groups] {a.sublibrary}: {len(m)} cells in {m['group'].nunique()} groups")
