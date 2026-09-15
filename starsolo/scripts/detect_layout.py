"""Find where the UMI and the three barcodes sit in the barcode read (fq2).

Rather than hard-code Parse's v3 read structure, scan a sample of reads for
exact matches to each round's barcode table at every offset, keep the top few
peaks per round as candidates, and choose the combination of offsets that puts
the most reads into split-pipe's called cells. Any offset given on the command
line (>= 0) is used as-is instead.

Peaks alone are not enough: a barcode that happens to equal part of a constant
linker "matches" nearly every read at the wrong offset, and when rounds 2 and 3
share a table their peaks are indistinguishable by sequence. Neither can
produce reads from real cells, so the whitelist settles both.

Also checks: fq1/fq2 not swapped, UMI not overlapping a barcode, fq1 read length
against the STAR index's sjdbOverhang.
"""
import argparse
import json
import sys
from itertools import islice, product

import dnaio

from parse_barcodes import (cell_key, correction_map, load_table, read_cells,
                            round1_value, well_value)

ap = argparse.ArgumentParser()
ap.add_argument("--fq1", required=True)
ap.add_argument("--fq2", required=True)
ap.add_argument("--cells", required=True)
ap.add_argument("--bc1", required=True)
ap.add_argument("--bc2", required=True)
ap.add_argument("--bc3", required=True)
ap.add_argument("--max-mismatch", type=int, default=1)
ap.add_argument("--sample-reads", type=int, default=200000)
ap.add_argument("--umi-start", type=int, default=0)
ap.add_argument("--umi-len", type=int, default=10)
ap.add_argument("--bc1-start", type=int, default=-1)
ap.add_argument("--bc2-start", type=int, default=-1)
ap.add_argument("--bc3-start", type=int, default=-1)
ap.add_argument("--candidates", type=int, default=3, help="peaks tried per round")
ap.add_argument("--min-rate", type=float, default=0.3)
ap.add_argument("--sjdb-overhang", type=int, default=0)
ap.add_argument("--out", required=True)
a = ap.parse_args()

warnings = []


def warn(msg):
    warnings.append(msg)
    print(f"[layout] WARNING: {msg}", file=sys.stderr)


with dnaio.open(a.fq1, a.fq2) as reader:
    pairs = [(r1.sequence, r2.sequence) for r1, r2 in islice(reader, a.sample_reads)]
n = len(pairs)
if n == 0:
    sys.exit("no reads")
seq1 = [p[0] for p in pairs]
seq2 = [p[1] for p in pairs]


def hist(seqs, keys, length):
    counts = [0] * max(0, max(map(len, seqs)) - length + 1)
    for s in seqs:
        for o in range(len(s) - length + 1):
            if s[o:o + length] in keys:
                counts[o] += 1
    return counts


def peaks(counts, length, k):
    """Top-k offsets, at least `length` apart."""
    chosen = []
    for o in sorted(range(len(counts)), key=counts.__getitem__, reverse=True):
        if counts[o] == 0 or len(chosen) == k:
            break
        if all(abs(o - c) >= length for c in chosen):
            chosen.append(o)
    return chosen


def top(counts, k=5):
    order = sorted(range(len(counts)), key=counts.__getitem__, reverse=True)[:k]
    return {o: round(counts[o] / n, 4) for o in order}


t1, L1 = load_table(a.bc1, round1=True)
t2, L2 = load_table(a.bc2)
t3, L3 = load_table(a.bc3)
same23 = set(t2) == set(t3)
k = a.candidates

# --- candidate offsets per round ---------------------------------------------
h1 = hist(seq2, t1, L1)
if max(h1 or [0]) / n < a.min_rate:
    swapped = max(hist(seq1, t1, L1) or [0]) / n
    if swapped >= a.min_rate:
        sys.exit(f"Round-1 barcodes match fq1 ({swapped:.1%}) better than fq2 "
                 f"({max(h1 or [0]) / n:.1%}): fq1/fq2 look swapped. fq1 must be cDNA, "
                 "fq2 the barcode read.")
    sys.exit(f"Round-1 barcodes match at most {max(h1 or [0]) / n:.1%} of fq2 reads "
             f"(top offsets {top(h1)}). Wrong barcode table for this kit/chemistry?")

h2 = hist(seq2, t2, L2)
h3 = h2 if same23 else hist(seq2, t3, L3)
histograms = {"bc1": top(h1), "bc2": top(h2), "bc3": top(h3)}

cand1 = [a.bc1_start] if a.bc1_start >= 0 else peaks(h1, L1, k)
n23 = k + 1 if same23 else k          # shared table: both rounds' peaks in one list
cand2 = [a.bc2_start] if a.bc2_start >= 0 else peaks(h2, L2, n23)
cand3 = [a.bc3_start] if a.bc3_start >= 0 else peaks(h3, L3, n23)

# --- score every non-overlapping combination against the called cells --------
c1 = correction_map(t1, round1_value, a.max_mismatch)
c2 = correction_map(t2, well_value, a.max_mismatch)
c3 = correction_map(t3, well_value, a.max_mismatch)
cells, _, _ = read_cells(a.cells)


def decode(cmap, offset, length):
    get = cmap.get
    return [get(s[offset:offset + length]) for s in seq2]


dec1 = {o: decode(c1, o, L1) for o in cand1}
dec2 = {o: decode(c2, o, L2) for o in cand2}
dec3 = {o: dec2[o] if same23 and o in dec2 else decode(c3, o, L3) for o in cand3}

umi = (a.umi_start, a.umi_start + a.umi_len)


def overlaps(spans):
    spans = sorted(spans)
    return any(s2 < e1 for (_, e1), (s2, _) in zip(spans, spans[1:]))


scored = []
for b1, b2, b3 in product(cand1, cand2, cand3):
    if overlaps([(b1, b1 + L1), (b2, b2 + L2), (b3, b3 + L3), umi]):
        continue
    decoded = hits = 0
    for v1, v2, v3 in zip(dec1[b1], dec2[b2], dec3[b3]):
        if v1 is None or v2 is None or v3 is None:
            continue
        decoded += 1
        hits += cell_key(v1 >> 1, v2, v3) in cells
    scored.append({"bc1_start": b1, "bc2_start": b2, "bc3_start": b3,
                   "decoded_rate": round(decoded / n, 4), "called_cell_rate": round(hits / n, 4)})

if not scored:
    sys.exit(f"No non-overlapping offset combination (candidates bc1 {cand1}, bc2 {cand2}, "
             f"bc3 {cand3}, UMI {umi}). Check umi_start / umi_len.")
scored.sort(key=lambda x: x["called_cell_rate"], reverse=True)
best = scored[0]
b1, b2, b3 = best["bc1_start"], best["bc2_start"], best["bc3_start"]
print(f"[layout] candidates bc1 {cand1}, bc2 {cand2}, bc3 {cand3}; best combinations:")
for s in scored[:3]:
    print(f"[layout]   {s}")

if best["called_cell_rate"] < 0.05:
    sys.exit(f"Best combination puts only {best['called_cell_rate']:.2%} of sampled reads in "
             "called cells. The cell list does not belong to this sublibrary, or the "
             "barcode tables are wrong.")
if len(scored) > 1 and scored[1]["called_cell_rate"] > 0.5 * best["called_cell_rate"]:
    warn(f"runner-up offsets {scored[1]} are close to the chosen ones")

rate1 = h1[b1] / n
if rate1 < a.min_rate:
    warn(f"round-1 exact-match rate at the chosen offset is only {rate1:.1%}")
if min(b1, b2, b3) != umi[1]:
    warn(f"first barcode starts at {min(b1, b2, b3)}, not right after the UMI ({umi[1]}); "
         "check umi_start / umi_len")

need = max(b1 + L1, b2 + L2, b3 + L3, umi[1])
short = sum(len(s) < need for s in seq2) / n
if short > 0.01:
    warn(f"{short:.1%} of fq2 reads are shorter than the {need} nt the layout needs")
r1_len = max(map(len, seq1))
if a.sjdb_overhang and r1_len - 1 > a.sjdb_overhang:
    warn(f"fq1 reads are up to {r1_len} nt but the index has sjdbOverhang "
         f"{a.sjdb_overhang}; fine, but {r1_len - 1} would be ideal")

layout = {
    "umi_start": a.umi_start, "umi_len": a.umi_len,
    "bc1_start": b1, "bc1_len": L1,
    "bc2_start": b2, "bc2_len": L2,
    "bc3_start": b3, "bc3_len": L3,
    "max_mismatch": a.max_mismatch,
    "sampled_reads": n,
    "bc1_exact_rate": round(rate1, 4),
    "all_rounds_decoded_rate": best["decoded_rate"],
    "called_cell_rate": best["called_cell_rate"],
    "fq1_max_len": r1_len,
    "fq2_max_len": max(map(len, seq2)),
    "bc2_bc3_share_table": same23,
    "combinations_top3": scored[:3],
    "offset_histograms_top5": histograms,
    "warnings": warnings,
}
with open(a.out, "w") as fh:
    json.dump(layout, fh, indent=2)
print(json.dumps(layout, indent=2))
