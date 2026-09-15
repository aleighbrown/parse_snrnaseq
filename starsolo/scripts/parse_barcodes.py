"""Parse Evercode barcode tables, well numbering and 1-mismatch correction.

Shared by make_whitelist.py, detect_layout.py and extract_reads.py.

Conventions
-----------
* Wells are numbered 1..96 row-major on the 8x12 plate (A1=1 .. A12=12,
  B1=13 .. H12=96). This is the numbering in split-pipe's `bc_wells`:
  "01_01_67" = bc1 A1, bc2 A1, bc3 F7.
* Round 1 has two primers per well: oligo-dT ("T") and random hexamer ("R").
  split-pipe merges both onto one cell and so do we -- every cell barcode we
  write uses the round-1 T sequence, whichever primer the read came from.
* Cell keys are packed as ((w1 * 97) + w2) * 97 + w3, so a lookup is one int.
"""
import csv
import re

ROWS = "ABCDEFGH"
_WELL_RE = re.compile(r"^([A-H])(\d{1,2})$")
_SEQ_RE = re.compile(r"^[ACGT]+$")


def well_index(name):
    """'F7' -> 67. None for anything that is not a well of an 8x12 plate."""
    m = _WELL_RE.match(str(name).strip().upper())
    if not m or not 1 <= int(m.group(2)) <= 12:
        return None
    return ROWS.index(m.group(1)) * 12 + int(m.group(2))


def well_name(idx):
    r, c = divmod(int(idx) - 1, 12)
    return f"{ROWS[r]}{c + 1}"


def cell_key(w1, w2, w3):
    return (w1 * 97 + w2) * 97 + w3


def _pick_column(path, fields, rows, names, test, what, required=True):
    """Column by name if present, else the single column whose values all pass."""
    lower = {f.lower(): f for f in fields}
    for n in names:
        if n in lower:
            return lower[n]
    hits = [f for f in fields if rows and all(test(r[f] or "") for r in rows)]
    if len(hits) == 1:
        return hits[0]
    if not hits and not required:
        return None
    raise ValueError(
        f"{path}: cannot tell which column holds the {what} "
        f"(candidates: {hits or 'none'}; columns: {fields})"
    )


def load_table(path, round1=False):
    """Read one split-pipe barcode CSV.

    Returns (table, length) where table maps sequence ->
    (well_index, is_R, primer_type).
    well_index is 0 for rows that are not wells of the 96-well plate; those
    sequences are still used for correction (so a read is not "corrected" onto
    an in-plate well it does not belong to) but never match a cell.
    """
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fields = reader.fieldnames or []
    if not rows:
        raise ValueError(f"{path}: empty barcode table")

    seq_col = _pick_column(path, fields, rows, ["sequence", "seq"],
                           lambda v: bool(_SEQ_RE.match(v.strip().upper())),
                           "barcode sequence")
    well_col = _pick_column(path, fields, rows, ["well"],
                            lambda v: well_index(v) is not None, "plate well")
    type_col = _pick_column(path, fields, rows, ["stype", "type"],
                            lambda v: v.strip().upper() in ("T", "R", "L"),
                            "primer type", required=False)

    table, lengths, per_well = {}, set(), {}
    for r in rows:
        seq = r[seq_col].strip().upper()
        idx = well_index(r[well_col]) or 0
        stype = r[type_col].strip().upper() if type_col else ""
        if seq in table and table[seq][0] != idx:
            raise ValueError(f"{path}: sequence {seq} listed for two wells")
        table[seq] = (idx, stype == "R", stype)
        lengths.add(len(seq))
        per_well[idx] = per_well.get(idx, 0) + 1

    if len(lengths) != 1:
        raise ValueError(f"{path}: barcodes of mixed lengths {sorted(lengths)}")
    if round1 and not type_col and any(n > 1 for w, n in per_well.items() if w):
        raise ValueError(
            f"{path}: wells appear more than once but there is no primer-type "
            "column, so oligo-dT and random-hexamer rows cannot be told apart"
        )
    return table, lengths.pop()


def well_sequences(table, path):
    """well_index -> the sequence that names that well.

    The oligo-dT ("T") row when there is one; otherwise the well's single
    non-R row. Every other row of the well still maps reads to it.
    """
    by_well = {}
    for seq, (idx, _, stype) in table.items():
        if idx and stype != "R":
            by_well.setdefault(idx, []).append((stype, seq))
    out = {}
    for idx, rows in by_well.items():
        t_rows = [s for t, s in rows if t == "T"]
        pick = t_rows if t_rows else [s for _, s in rows]
        if len(pick) != 1:
            raise ValueError(f"{path}: well {well_name(idx)} has no single naming "
                             f"sequence (rows by type: {sorted(rows)})")
        out[idx] = pick[0]
    return out


def describe_table(table, path):
    """One line per primer type: rows, and how many fall on the 96-well plate."""
    counts = {}
    for idx, _, stype in table.values():
        c = counts.setdefault(stype or "-", [0, 0, set()])
        c[0] += 1
        if idx:
            c[1] += 1
            c[2].add(idx)
    lines = [f"{path}: {len(table)} sequences"]
    for t, (n, on_plate, wells) in sorted(counts.items()):
        lines.append(f"  type {t}: {n} rows, {on_plate} on plate, {len(wells)} wells")
    return "\n".join(lines)


def correction_map(table, value, max_mismatch=1):
    """sequence -> value(entry), plus every single-substitution neighbour
    (N included) when max_mismatch is 1.

    A neighbour reachable from two different wells is ambiguous and dropped.
    Exact sequences always win. Round 1's T and R rows of the SAME well do not
    count as a clash -- they are the same cell.
    """
    if max_mismatch not in (0, 1):
        raise ValueError("max_mismatch must be 0 or 1")
    out = {}
    if max_mismatch:
        near, clash = {}, set()
        for seq, entry in table.items():
            for i, base in enumerate(seq):
                for b in "ACGTN":
                    if b == base:
                        continue
                    n = seq[:i] + b + seq[i + 1:]
                    if n in table:
                        continue
                    prev = near.get(n)
                    if prev is None:
                        near[n] = entry
                    elif prev[0] != entry[0]:
                        clash.add(n)
        out = {n: value(e) for n, e in near.items() if n not in clash}
    for seq, entry in table.items():
        out[seq] = value(entry)
    return out


def round1_value(entry):
    """Round-1 map value: well * 2 + is_R, so the primer type rides along."""
    return entry[0] * 2 + int(entry[1])


def well_value(entry):
    return entry[0]


def read_cells(path):
    """cells.tsv from make_whitelist.py -> (keys->row index, cb list, bc_wells list)."""
    keys, cbs, names = {}, [], []
    with open(path, newline="") as fh:
        for i, r in enumerate(csv.DictReader(fh, delimiter="\t")):
            keys[cell_key(int(r["w1"]), int(r["w2"]), int(r["w3"]))] = i
            cbs.append(r["cb"])
            names.append(r["bc_wells"])
    return keys, cbs, names
