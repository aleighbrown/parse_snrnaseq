"""Pull the cell-barcode column out of split-pipe's cell_metadata.csv and write
a plain one-per-line barcodes.tsv for cellsnp-lite.

The values MUST match the cell-barcode tag inside the BAM exactly. Parse's
metadata may store the barcode split across columns or with a well suffix; if
cellsnp reports ~0 cells, inspect the BAM tag (`samtools view <bam> | head`)
and set splitpipe.barcode_column (and, if needed, adjust the transform below).
"""
import pandas as pd

meta_path = snakemake.input.cell_meta          # noqa: F821
out_path = snakemake.output.barcodes           # noqa: F821
col = snakemake.params.col                      # noqa: F821

meta = pd.read_csv(meta_path)

if col not in meta.columns:
    raise KeyError(
        f"Column '{col}' not in {meta_path}. Available: {list(meta.columns)}. "
        "Set splitpipe.barcode_column in config/config.yaml."
    )

barcodes = meta[col].astype(str).dropna().drop_duplicates()
barcodes.to_csv(out_path, index=False, header=False)
print(f"[extract_barcodes] wrote {len(barcodes)} barcodes -> {out_path}")
