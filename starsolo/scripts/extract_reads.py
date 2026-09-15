"""Keep only read pairs from split-pipe's called cells, and rewrite the barcode
read into a form STARsolo parses trivially.

In:   fq1 = cDNA, fq2 = Parse barcode read (layout from detect_layout.py)
Out:  R1  = the cDNA read, untouched
      R2  = 24-nt corrected cell barcode + UMI   (bc1(T) + bc2 + bc3 + UMI)

Barcodes are corrected per round (exact, else a unique 1-substitution
neighbour). Random-hexamer round-1 reads are folded onto their well's oligo-dT
barcode, as split-pipe does. Pairs whose barcode is not a called cell are
dropped.

Also writes per-cell read counts (split by T/R primer) for validation against
split-pipe's mread_count, and a JSON of where every read went.
"""
import argparse
import json
import sys
import time

import dnaio

from parse_barcodes import (cell_key, correction_map, load_table, read_cells,
                            round1_value, well_value)

ap = argparse.ArgumentParser()
ap.add_argument("--fq1", required=True)
ap.add_argument("--fq2", required=True)
ap.add_argument("--cells", required=True)
ap.add_argument("--layout", required=True)
ap.add_argument("--bc1", required=True)
ap.add_argument("--bc2", required=True)
ap.add_argument("--bc3", required=True)
ap.add_argument("--out1", required=True)
ap.add_argument("--out2", required=True)
ap.add_argument("--counts", required=True)
ap.add_argument("--stats", required=True)
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--max-reads", type=int, default=0, help="stop after N pairs (testing)")
a = ap.parse_args()

with open(a.layout) as fh:
    lay = json.load(fh)

t1, L1 = load_table(a.bc1, round1=True)
t2, L2 = load_table(a.bc2)
t3, L3 = load_table(a.bc3)
if (L1, L2, L3) != (lay["bc1_len"], lay["bc2_len"], lay["bc3_len"]):
    sys.exit("barcode lengths differ from read_layout.json -- rerun layout detection")
get1 = correction_map(t1, round1_value, lay["max_mismatch"]).get
get2 = correction_map(t2, well_value, lay["max_mismatch"]).get
get3 = correction_map(t3, well_value, lay["max_mismatch"]).get
cells, cbs, names = read_cells(a.cells)
getcell = cells.get

s1, e1 = lay["bc1_start"], lay["bc1_start"] + L1
s2, e2 = lay["bc2_start"], lay["bc2_start"] + L2
s3, e3 = lay["bc3_start"], lay["bc3_start"] + L3
su, eu = lay["umi_start"], lay["umi_start"] + lay["umi_len"]
cb_qual = "F" * len(cbs[0])

n_t = [0] * len(cbs)
n_r = [0] * len(cbs)
total = fail1 = fail2 = fail3 = not_cell = 0
Record = dnaio.SequenceRecord
t0 = time.time()

# Reading and writing gzip both go through xopen's external (i)gzip processes.
rthreads = max(1, a.threads // 2)
with dnaio.open(a.fq1, a.fq2, open_threads=rthreads) as reader, \
        dnaio.open(a.out1, a.out2, mode="w", compression_level=1,
                   open_threads=max(1, a.threads - rthreads)) as writer:
    for r1, r2 in reader:
        total += 1
        s = r2.sequence
        v1 = get1(s[s1:e1])
        if v1 is None:
            fail1 += 1
            continue
        v2 = get2(s[s2:e2])
        if v2 is None:
            fail2 += 1
            continue
        v3 = get3(s[s3:e3])
        if v3 is None:
            fail3 += 1
            continue
        idx = getcell(cell_key(v1 >> 1, v2, v3))
        if idx is None:
            not_cell += 1
            continue
        if v1 & 1:
            n_r[idx] += 1
        else:
            n_t[idx] += 1
        writer.write(r1, Record(r2.name, cbs[idx] + s[su:eu],
                                cb_qual + r2.qualities[su:eu]))
        if total % 50_000_000 == 0:
            rate = total / (time.time() - t0)
            print(f"[extract] {total:,} pairs, {rate:,.0f}/s", file=sys.stderr, flush=True)
        if a.max_reads and total >= a.max_reads:
            break

elapsed = time.time() - t0
kept = sum(n_t) + sum(n_r)

with open(a.counts, "w") as fh:
    fh.write("cb\tbc_wells\treads_T\treads_R\treads\n")
    for i, cb in enumerate(cbs):
        fh.write(f"{cb}\t{names[i]}\t{n_t[i]}\t{n_r[i]}\t{n_t[i] + n_r[i]}\n")

stats = {
    "pairs_total": total,
    "fail_bc1": fail1,
    "fail_bc2": fail2,
    "fail_bc3": fail3,
    "valid_barcode_not_called_cell": not_cell,
    "kept": kept,
    "kept_fraction": round(kept / total, 4) if total else 0,
    "kept_random_hexamer_fraction": round(sum(n_r) / kept, 4) if kept else 0,
    "cells_with_reads": sum(1 for i in range(len(cbs)) if n_t[i] or n_r[i]),
    "cells_in_whitelist": len(cbs),
    "seconds": round(elapsed),
    "pairs_per_second": round(total / elapsed) if elapsed else 0,
}
with open(a.stats, "w") as fh:
    json.dump(stats, fh, indent=2)
print(json.dumps(stats, indent=2))
