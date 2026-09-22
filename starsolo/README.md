# STARsolo re-alignment of split-pipe's cells — BAMs for splicing

The split-pipe run deleted its BAMs (`clean` removes `process/*.bam`). Instead of
rerunning split-pipe (~60 core-hours per sublibrary, mostly re-aligning reads
from barcodes that were never cells), this workflow keeps only the read pairs
from cells split-pipe **already called**, and aligns those with STARsolo.

```
split-pipe cell_metadata.csv ─► whitelist: 24-nt barcode per called cell
fq2 (200k reads) ─► read_layout.json: where UMI / bc1 / bc2 / bc3 sit
fq1 + fq2 ─► extract: called-cell pairs only; R2 := corrected barcode + UMI
     └─► validation: reads per cell vs split-pipe's mread_count
GENCODE v42 ─► STAR index  (shared, built once)
2M reads ─► strand test (Forward vs Reverse)
R1 ─► STAR pass 1 ─► novel junctions, pooled over all sublibraries
R1 + R2 ─► STARsolo pass 2 ─► sorted BAM (CB/UB) + Gene / GeneFull / SJ matrices
     └─► target `split`: one BAM per donor × sort
```

Pilot: **ParseFANS1 only** (`config/sublibraries.tsv`). It exists to answer two
questions before scaling: does our barcode handling agree with split-pipe's, and
how long does each step actually take?

## Before the first run

### 0. Tools (no conda)

The cluster's conda couldn't build envs for this workflow (glibc too old for
current conda-forge builds, then a broken channel index), so the workflow
doesn't use conda at all. `tools` in `config/config.yaml` sets what every job
uses:

- **`bin_dirs`** are put at the front of `PATH` in every job. Between them they
  must provide `python` 3 with pandas and numpy, plus `pigz`, `sinto` and
  `samtools`. Your base conda env's `bin/` is the intended choice. Any
  samtools with `index -@` works (1.9 does), since it only sorts and indexes.
  Without pigz, extraction falls back to gzip and is slower.
- **`star`** is the path to one STAR binary. The static Linux build from the
  STAR GitHub release (`STAR_2.7.11b.zip` → `Linux_x86_64_static/STAR`) runs on
  any glibc.
- **`star_version`** must match what that binary reports, and the STAR that built
  `reference.dir`. Any job that runs STAR stops with an error if the version
  differs.

Check all of this on the login node in seconds:

```bash
source submit.sh tools
```

Snakemake wrappers don't get around the conda problem. Each wrapper ships its
own env file, and Snakemake builds that env with the same conda.

### 1. Parse barcode tables (set)

`parse_barcodes` in `config/config.yaml` points at the tables split-pipe
recorded for this run, in ParseFANS1's `process/run_proc_def.json`
(`barcode.bc_round_set`). All three live in the spipe env's
`site-packages/splitpipe/barcodes/`:

| round | table | split-pipe's description |
|---|---|---|
| bc1 | `bc_data_n302_R1_v3_8.csv` | 195 barcodes, 96 wells (original count 302) |
| bc2 | `bc_data_v1.csv` | 96 barcodes, 96 wells |
| bc3 | `bc_data_R3_v3.csv` | 96 barcodes, 96 wells |

The loader looks for `sequence`, `well` and `stype` columns and fails with the
column list if it can't find them. What each table contains is logged in
`logs/whitelist/<sublib>.log` (rows per primer type, on-plate rows, wells):

- **bc1 should show 96 wells of `T` and 96 of `R`.** Each well is named by its
  `T` sequence, and all its on-plate rows count as that well.
- **Off-plate rows are used only for error correction.** They never match a cell.
- **Check the 195-vs-302 difference.** If the log shows more on-plate round-1
  rows than split-pipe's 195, say so before extracting.

Nothing about the v3 read structure is hard-coded. `read_layout` takes the top
few offsets where each round's barcodes match, then picks the combination that
puts the most sampled reads into called cells. That handles two traps: a
barcode that happens to equal part of a constant linker (a huge but useless
match), and rounds that share a table. The top three combinations are in
`read_layout.json`. The UMI is assumed to be the first 10 nt of fq2
(`layout.umi_*`).

### 2. The cell list (the whitelist)

`cells.file` is split-pipe's **combined**
`splitpipe/combined/all-sample/DGE_filtered/cell_metadata.csv`: every cell
split-pipe called, all sublibraries, barcodes like `01_01_67__s1`. This is the
same file as `../data/DGE_filtered/cell_metadata.csv`. Rows are kept per
sublibrary by `cells_tag` (`s1` → ParseFANS1).

`cells.subset` narrows that list to the cells you actually want aligned: the
QC-passed, cell-typed, non-control table from R (`bc_wells`, with the `__s1`
suffix). Everything else still comes from split-pipe's `cell_metadata.csv`,
including `mread_count`, which the extraction is validated against. The job
stops if a cell in your table is not in split-pipe's list, since that would
mean the two came from different runs.

- **Fewer reads to align**, in proportion to the cells you drop.
- **But the cell list is now baked into the alignment.** Changing it means
  re-extracting and re-aligning, so settle QC and cell typing first. Set
  `subset.file: null` to align every called cell instead.
- **`called_cell_rate` in `read_layout.json` drops with it.** With ~180k of
  426k cells, expect roughly half the ~0.4 the layout step reports for the
  full list.

Use the same file for `split.annotation`, so every aligned cell is also
grouped.

### 3. The STAR reference

Built at `/SAN/vyplab/vyplab_reference_genomes/STAR/human_v42_STAR2.7.11b`: the
same GENCODE v42 FASTA + GTF as the split-pipe reference, STAR 2.7.11b
(`tools.star_version`), `sjdbOverhang 149`, with a
`build_provenance.txt` inside. Contigs are plain `chr1`, **not** split-pipe's
`GRCh38-gencode-v42_chr1`.

If the directory already holds a finished index (`SA` + `genomeParameters.txt`)
there is **no build rule at all**, so no run can delete and rebuild the shared
index. `STAR/human_v42` already exists; if its `genomeParameters.txt` shows
the same FASTA/GTF, you can set `reference.dir` to it instead.

## Run

From this directory (`starsolo/`), in screen on the login node:

```bash
snakemake -n                      # dry run
source submit.sh tools            # tool check (seconds)
source submit.sh reference        # STAR index (can go first; needs no barcode tables)
source submit.sh layout           # whitelist + read layout: minutes
source submit.sh extract          # + extraction + validation vs split-pipe
source submit.sh                  # everything through the STARsolo BAM
source submit.sh split            # per-sublibrary donor × sort BAMs
source submit.sh splicing         # merged donor × sort BAMs + junction matrix
```

Going step by step is the point of the pilot. **Read each output before
starting the next step.**

| after | read | looking for |
|---|---|---|
| `layout` | `extract/ParseFANS1/read_layout.json` | `bc1_exact_rate` high (≳0.8), `called_cell_rate` ≈ fraction of reads in cells (~0.4), empty `warnings` |
| `extract` | `extract/ParseFANS1/validation/summary.txt` | `verdict: PASS`: Spearman ≥ 0.9, ~100% cells with reads, no round-1 wells flagged |
| `extract` | `extract/ParseFANS1/extract_stats.json` | `pairs_per_second`, which sizes the other 15 |
| strand | `strand/ParseFANS1/strand_test.tsv` | one strand clearly wins |
| `all` | `solo/ParseFANS1/Log.final.out`, `Solo.out/GeneFull/Summary.csv` | unique mapping rate; `Reads With Valid Barcodes` ≈ 1.0 |
| `splicing` | `split/cells_per_group.tsv` | cells per donor × sort, and which groups fell below `split.min_cells` |
| `splicing` | `junctions/sj_counts.report.tsv` | cells grouped per sublibrary (should be all of them), junctions found |

Timings for every job are in `submissions/<run>/` logs.

## Outputs (`results_dir`)

```
whitelist/{sublib}/cells.tsv, whitelist.txt
extract/{sublib}/read_layout.json
extract/{sublib}/R1.cdna.fastq.gz, R2.cb_umi.fastq.gz     (the filtered reads)
extract/{sublib}/cell_read_counts.tsv, extract_stats.json
extract/{sublib}/validation/summary.txt, per_cell.tsv, per_round1_well.tsv
strand/{sublib}/strand.txt, strand_test.tsv
pass1/{sublib}/SJ.out.tab
pass1/pooled_junctions.sjdb.tsv, pooled_junctions.report.txt
solo/{sublib}/Aligned.sortedByCoord.out.bam(.bai)
solo/{sublib}/Solo.out/{Gene,GeneFull,SJ}/raw/             (per-cell matrices)
split/groups.tsv, cells_per_group.tsv                      (donor × sort groups)
split/{sublib}/cell_groups.tsv, bams/<donor>__<sort>.bam   (per sublibrary)
groups/<donor>__<sort>.bam(.bai)          (merged over sublibraries: target `bams`)
junctions/sj_counts.by_group.tsv.gz       (junction × group counts: `junctions`)
junctions/sj_counts.report.tsv
```

- **Cell barcodes.** The BAM `CB` tag is the 24-nt `bc1(T)+bc2+bc3` sequence;
  `whitelist/{sublib}/cells.tsv` maps it back to `bc_wells` and split-pipe's
  metadata. CBs are only unique **within** a sublibrary.
- **Group names.** Donor × sort names are made filesystem-safe: `P89/08` →
  `P89-08`, giving `P89-08__K_negative`.
- **Two junction routes.** `junctions/sj_counts.by_group.tsv.gz` comes straight
  from STARsolo's per-cell `SJ` matrices summed over each group's cells: one
  row per junction, one column per donor × sort, `chr`/`start`/`end` first.
  The merged BAMs in `groups/` are for anything that needs reads rather than
  counts — regtools + leafcutter, MAJIQ, or IGV.
- **Groups are global.** `split.min_cells` counts a donor × sort over all 16
  sublibraries, so a group is kept or dropped everywhere at once.

## Things to know

- **Duplicates are not removed from the BAMs.** About 74% of reads are PCR
  duplicates. The `Solo.out/SJ` matrix is already UMI-collapsed. For junction
  counts from BAMs, count distinct `CB`+`UB` per junction rather than
  deduplicating by read position: fragmentation gives copies of one molecule
  different start positions.
- **Adding sublibraries reruns pass 2 everywhere.** It changes
  `pooled_junctions.sjdb.tsv`, so ParseFANS1's BAM will be regenerated when the
  other 15 are added. That is intentional: every sample aligns against the same
  junction set.
- **Cell barcodes are only unique within a sublibrary.** The merged group BAMs
  carry `RG` per sublibrary (`samtools merge -r`), so a nucleus is `RG` + `CB`,
  and a molecule is `RG` + `CB` + `UB`.
- **Scratch space.** STAR's BAM sort spills to `tmp_dir` (default
  `<results_dir>/tmp`). The 50 G `tscratch` used for split-pipe is too small.
- **Resources are estimates.** `config/cluster.yaml` walltimes are generous;
  trim them once the pilot reports real numbers.
