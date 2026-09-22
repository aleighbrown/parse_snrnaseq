"""Sum STARsolo's per-cell SJ counts into one junction x (donor x sort) matrix.

Reads each sublibrary's Solo.out/SJ/raw (matrix.mtx + barcodes.tsv +
features.tsv) and the cb->group table for that sublibrary, and adds up the
counts of all of a group's cells. Cells with no group (not in the annotation
table, or in a group below min-cells) are ignored.

Junctions are keyed by their features.tsv line, so the same junction found in
several sublibraries lands on one row even when their feature lists differ.

Output: a tab-separated matrix, one row per junction, one column per group,
gzipped. Junctions with no counts in any kept cell are left out.
"""
import argparse
import gzip
import sys
from collections import defaultdict

import numpy as np


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_lines(path):
    with opener(path) as fh:
        return [line.rstrip("\n") for line in fh]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solo", nargs="+", required=True,
                    help="Solo.out/SJ/raw directories, one per sublibrary")
    ap.add_argument("--groups", nargs="+", required=True,
                    help="cell_groups.tsv (cb<TAB>group), same order as --solo")
    ap.add_argument("--sublibraries", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()

    if not len(a.solo) == len(a.groups) == len(a.sublibraries):
        sys.exit("--solo, --groups and --sublibraries must have the same length")

    # One pass to fix the column order, so every sublibrary writes into the
    # same columns.
    group_names = sorted({line.split("\t")[1]
                          for path in a.groups for line in read_lines(path) if line})
    if not group_names:
        sys.exit("no cell -> group assignments in " + ", ".join(a.groups))
    col = {g: i for i, g in enumerate(group_names)}

    counts = defaultdict(lambda: np.zeros(len(group_names), dtype=np.int64))
    report = []

    for solo, groups, sublib in zip(a.solo, a.groups, a.sublibraries):
        cb_group = {}
        for line in read_lines(groups):
            if line:
                cb, g = line.split("\t")[:2]
                cb_group[cb] = col[g]

        barcodes = read_lines(solo + "/barcodes.tsv")
        features = read_lines(solo + "/features.tsv")
        # barcode index (1-based, as in the mtx) -> column, or -1 to skip
        bc_col = np.full(len(barcodes) + 1, -1, dtype=np.int64)
        for i, cb in enumerate(barcodes, start=1):
            bc_col[i] = cb_group.get(cb.split("\t")[0], -1)
        used = int((bc_col > -1).sum())

        rows = cols = nnz = kept = total = 0
        with opener(solo + "/matrix.mtx") as fh:
            for line in fh:
                if line.startswith("%"):
                    continue
                rows, cols, nnz = (int(x) for x in line.split())
                break
            if cols != len(barcodes) or rows != len(features):
                sys.exit(f"{solo}: matrix is {rows}x{cols} but there are "
                         f"{len(features)} features and {len(barcodes)} barcodes")
            for line in fh:
                r, c, v = line.split()
                j = bc_col[int(c)]
                if j < 0:
                    continue
                counts[features[int(r) - 1]][j] += int(v)
                kept += 1
                total += int(v)

        report.append(f"{sublib}\t{len(barcodes)}\t{used}\t{nnz}\t{kept}\t{total}")
        print(f"[sj] {sublib}: {used}/{len(barcodes)} cells grouped, "
              f"{kept}/{nnz} matrix entries kept, {total} counts", file=sys.stderr)

    # features.tsv for SJ is chr, intron start, intron end (+ any extra
    # columns STAR writes); it is passed through as-is.
    n_extra = max(len(k.split("\t")) for k in counts) if counts else 0
    head = ["chr", "start", "end"] + [f"feature_col{i}" for i in range(4, n_extra + 1)]
    with gzip.open(a.out, "wt") as out:
        out.write("\t".join(head[:n_extra] + group_names) + "\n")
        for key in sorted(counts, key=lambda k: (k.split("\t")[0],
                                                 int(k.split("\t")[1]),
                                                 int(k.split("\t")[2]))):
            row = counts[key]
            if not row.any():
                continue
            out.write(key + "\t" + "\t".join(map(str, row.tolist())) + "\n")

    with open(a.report, "w") as fh:
        fh.write("sublibrary\tcells\tcells_grouped\tmatrix_entries\tentries_kept\tcounts\n")
        fh.write("\n".join(report) + "\n")
        fh.write(f"\njunctions\t{len(counts)}\ngroups\t{len(group_names)}\n")

    print(f"[sj] {len(counts)} junctions x {len(group_names)} groups -> {a.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
