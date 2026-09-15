"""Pool first-pass STAR junctions across sublibraries into one sjdb file.

Every sublibrary's second pass inserts the SAME junction set, so a junction is
not detectable in one sample merely because that sample found it in pass 1.

Keeps junctions that are unannotated (annotated ones are already in the index),
canonical-motif (unless --keep-noncanonical), off chrM, and supported by at
least --min-unique uniquely mapped reads summed over all sublibraries.

SJ.out.tab columns: chr start end strand(0/1/2) motif annotated unique multi overhang
"""
import argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("sj", nargs="+")
ap.add_argument("--min-unique", type=int, default=3)
ap.add_argument("--keep-noncanonical", action="store_true")
ap.add_argument("--out", required=True)
ap.add_argument("--report", required=True)
a = ap.parse_args()

STRAND = {"0": ".", "1": "+", "2": "-"}
unique = defaultdict(int)
motif = {}
seen = annotated = 0
for path in a.sj:
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            seen += 1
            if f[5] == "1":
                annotated += 1
                continue
            key = (f[0], int(f[1]), int(f[2]), STRAND[f[3]])
            unique[key] += int(f[6])
            motif[key] = f[4]

kept = sorted(
    k for k, n in unique.items()
    if n >= a.min_unique
    and (a.keep_noncanonical or motif[k] != "0")
    and not k[0].endswith(("chrM", "MT"))
)
with open(a.out, "w") as fh:
    for chrom, start, end, strand in kept:
        fh.write(f"{chrom}\t{start}\t{end}\t{strand}\n")

report = (
    f"SJ.out.tab files:            {len(a.sj)}\n"
    f"rows read:                   {seen}\n"
    f"annotated (skipped):         {annotated}\n"
    f"distinct unannotated:        {len(unique)}\n"
    f"kept (>= {a.min_unique} unique reads):   {len(kept)}\n"
)
with open(a.report, "w") as fh:
    fh.write(report)
print(report)
