"""Pick --soloStrand from two short STARsolo runs (Forward vs Reverse).

Compares the fraction of reads assigned to a gene. GeneFull (exons + introns)
is used because nuclei are mostly intronic; Gene is reported alongside. The
right strand should win by a wide margin; if it doesn't, stop and look.
"""
import argparse
import os
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--forward", required=True, help="STAR output dir of the Forward run")
ap.add_argument("--reverse", required=True, help="STAR output dir of the Reverse run")
ap.add_argument("--min-fold", type=float, default=2.0)
ap.add_argument("--out", required=True)
ap.add_argument("--report", required=True)
a = ap.parse_args()


def unique_fraction(run_dir, feature):
    path = os.path.join(run_dir, "Solo.out", feature, "Summary.csv")
    with open(path) as fh:
        rows = dict(line.rstrip("\n").rsplit(",", 1) for line in fh if "," in line)
    key = f"Reads Mapped to {feature}: Unique {feature}"
    if key not in rows:
        sys.exit(f"'{key}' not in {path}; rows: {list(rows)}")
    return float(rows[key])


res = {s: {f: unique_fraction(d, f) for f in ("Gene", "GeneFull")}
       for s, d in (("Forward", a.forward), ("Reverse", a.reverse))}
fwd, rev = res["Forward"]["GeneFull"], res["Reverse"]["GeneFull"]
best, other = ("Forward", rev) if fwd >= rev else ("Reverse", fwd)
fold = max(fwd, rev) / max(other, 1e-9)

with open(a.report, "w") as fh:
    fh.write("strand\tGene\tGeneFull\n")
    for s, v in res.items():
        fh.write(f"{s}\t{v['Gene']:.4f}\t{v['GeneFull']:.4f}\n")
    fh.write(f"# chosen {best}, {fold:.1f}-fold over the other strand\n")
print(open(a.report).read())

if fold < a.min_fold:
    sys.exit(f"Strands differ only {fold:.2f}-fold on GeneFull; not choosing automatically. "
             "Set star.strand in config.yaml once you've looked at the report.")
with open(a.out, "w") as fh:
    fh.write(best + "\n")
