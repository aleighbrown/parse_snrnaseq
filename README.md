# Gatt FTLD-TDP FANS — Parse WT Evercode pipeline (Snakemake)

Nuclei from 36 FTLD-TDP donors, FANS-sorted into three populations, pooled and
run on one **Parse WT Evercode Mega** kit → 16 sublibraries (`ParseFANS1..16`).
This pipeline builds the genome reference, runs split-pipe, and assigns every
nucleus to a **named donor** using reference genotypes.

```
reference: FASTA + GTF ──► split-pipe --mode mkref ──► genome_dir
                                                          │
FASTQ (16 sublibraries: ParseFANS1..16) ──────────────────┘
  └─ split-pipe --mode all ──► per-sublibrary BAM + filtered cell barcodes
       ├─ split-pipe --mode comb ──► one combined DGE matrix (→ Pegasus)
       └─ cellsnp-lite ──► pileup of each cell at the reference SNP sites
            └─ vireo (-d donors.vcf) ──► each cell → named donor (+ doublet)
                 └─ aggregate ──► one assignment table + per-donor summary

combined DGE + assignment table
  └─ celltype ──► QC + doublets ──► Harmony clusters ──► Azimuth + marker labels
                   └─► celltype/cell_annotations.csv   (target: celltype)
```

## The experiment

Ariana Gatt sorted NeuN+ nuclei by TDP-43 and K status:

| Sort | Phenotype | Round-1 wells |
|---|---|---|
| `K_negative`   | NeuN+ TDP+ K− | 1–43 (A1–D7) |
| `TDP_negative` | NeuN+ TDP− K+ | 44–86 (D8–H2) |
| `positive`     | NeuN+ TDP+ K+ | 87–96 (H3–H12) |

All 36 donors were pooled **within** each sort, so the Round-1 barcode gives the
sort for free (via `config/parse_sample_list.txt`) and genotype demux is what
separates donors. Target was 1M nuclei; the collaborator's run recovered 425k
split-pipe cells → 326k passing QC → 250k assignable to a unique donor.

## Layout

```
config/
  config.yaml            main config (genome, paths, params, BAM tags)
  sublibraries.tsv       16 sublibraries → FASTQ pairs (relative to fastq_dir)
  donor_metadata.tsv     34 observed donors, derived from the assignment table
  donors_nygc_wgs.tsv    NYGC WGS metadata (not wired in — no crosswalk)
  parse_sample_list.txt  Round-1 well → sort map for split-pipe
  cluster.yaml           per-rule SGE resources (tmem / h_vmem / h_rt / smp)
submit.sh                qsub wrapper — runs snakemake as a cluster job
cluster_qsub.sh          SGE jobscript template ({exec_job})
workflow/
  Snakefile
  rules/  reference.smk   genome FASTA/GTF → split-pipe mkref
          genotypes.smk   per-donor VCFs → merged SNP panel
          parse.smk       split-pipe --mode all / --mode comb
          demux.smk       cellsnp-lite + vireo
          aggregate.smk   combine assignments across sublibraries
          celltype.smk    QC / doublets / Harmony / Azimuth on the combined DGE
  scripts/ extract_barcodes.py  aggregate_assignments.py  import_assignments.py
           celltype_annotation.R
  envs/   bcftools.yaml  demux.yaml  py.yaml  celltype.yaml (+ .post-deploy.sh)
```

## Stage 0 — the shared Parse human reference

This stage runs `split-pipe --mode mkref` into the shared reference tree so any project
needing a Parse human reference can point at the same build:

```
/SAN/vyplab/vyplab_reference_genomes/parse_ref/human/
  GRCh38_gencode_v42/            the split-pipe reference
  GRCh38_gencode_v42.mkref.done  provenance: sources, split-pipe version, who/when
  _sources/                      decompressed FASTA/GTF, shared too
```

Configured by `genome.ref_root` plus one entry per species in
`genome.assemblies`. Each gives its FASTA and GTF as **either** a local path
(`fasta:` / `gtf:`) **or** a URL (`fasta_url:` / `gtf_url:`). Currently pointed
at the vyplab GENCODE v42 build:

```yaml
genome:
  ref_root: /SAN/vyplab/vyplab_reference_genomes/parse_ref/human
  assemblies:
    - name: GRCh38_gencode_v42
      fasta: /SAN/vyplab/vyplab_reference_genomes/sequence/human/gencode/GRCh38.primary_assembly.genome.fa
      gtf:   /SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v42.annotation.sorted.gtf.gz
```

```bash
# build the shared reference (~64 GB RAM, an hour or two for human)
snakemake reference -c 16
```

**It builds once.** Every input is wrapped in `ancient()`, so no source
timestamp can invalidate the reference, and once `GRCh38_gencode_v42/` exists
Snakemake treats it as up to date. Delete the directory to force a rebuild.

Reusing it from another project: set `genome.prebuilt_dir` to
`.../parse_ref/human/GRCh38_gencode_v42` and no download or mkref rule is
created at all. Set `ref_root: null` to build into this project's own
`resources/genomes/` instead of the shared tree.

Other notes:

- `mkref` wants the **genome** FASTA plus the GTF.
  `gencode.v42.transcripts.fa` is a transcriptome FASTA and is not used.
- The GTF is decompressed once into `_sources/` because `mkref` needs plain
  text; the uncompressed FASTA is read where it lies, no copy.
- Two or more `assemblies` entries build a single **mixed-species** reference
  (parallel `--genome_name` / `--fasta` / `--genes` lists).
- **Contigs get the genome name as a prefix.** split-pipe's
  `mkref_chr_geno_prefix` is on, so this reference's contigs are
  `GRCh38-gencode-v42_chr1`, **not** `chr1` — and BAMs aligned to it will be
  too. Anything that matches on contig name (donor VCFs for the vireo branch,
  regtools for junction extraction, IGV) has to account for that prefix. On top
  of the usual GENCODE `chr1` vs Ensembl `1` question.
- **`--genome_name` cannot contain underscores.** split-pipe silently rewrites
  them to hyphens (`GRCh38_gencode_v42` → `GRCh38-gencode-v42`), which then
  becomes the contig prefix. The configured name is already hyphenated so the
  directory matches the reference's internal name.
- Measured for this build: **1h15m wall, 26 GB peak RSS**, `cpu_eff 0.08` at 16
  threads — STAR's `genomeGenerate` is close to single-threaded here, so extra
  slots mostly buy queue time.

## Tools

- **split-pipe** — Parse Biosciences' proprietary binary (v1.9.0). No conda
  recipe exists; it lives in the existing env
  `/SAN/vyplab/FUS_RNA/conda_envs/spipe` (`conda_envs.spipe` in config), which
  the split-pipe rules activate **in their own shell block**. Deliberately not
  via Snakemake's `conda:` directive — see below.
- **cellsnp-lite** + **vireo** — conda (`workflow/envs/demux.yaml`) or the
  [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/)
  container.
- **bcftools / samtools** — genotype panel prep and BAM indexing.

## Run — UCL CS cluster (SGE)

Follows the [frattalab/rna_seq_snakemake](https://github.com/frattalab/rna_seq_snakemake)
pattern, with `config/cluster.yaml` giving per-rule resources and
`cluster_qsub.sh` as the jobscript template.

UCL does not allow a cluster job to submit further jobs, so the snakemake
controller runs **on the login node, inside a screen session**, and fans each
rule out as its own `qsub`. `submit.sh` is `source`d, not qsub-ed:

```bash
snakemake reference -n -p       # dry run first

screen -S parse
source submit.sh                # full run (rule all)
source submit.sh splitpipe_all  # all 16 sublibraries + comb, WITHOUT demux
source submit.sh reference      # build the shared Parse reference only
source submit.sh all rerun      # tag the submission "all_rerun"
# ctrl-a d to detach; `screen -r parse` to come back
```

The screen session must stay alive for the whole run — it *is* the controller.
`submit.sh` deliberately contains no `set -e` / `set -u` / `exit`: sourced, each
of those would act on your interactive shell and kill the session on any error.

Each submission gets `submissions/<timestamp>_<name>/` holding the per-rule job
logs, `snakemake.log` (the controller's own output, tee'd so it survives a lost
screen session), and a copy of the `config.yaml` and `cluster.yaml` that
produced it.

`split-pipe` does **not** need to be on PATH: `mkref`, `split_pipe` and
`split_pipe_comb` activate the spipe env themselves. Only snakemake has to be
available when you submit.

> **Why not the `conda:` directive.** Snakemake 7 activates a path fine, but
> treats any non-YAML `conda:` value as an environment *name* when recording
> provenance, running `conda env export --name <value>`. Conda cannot resolve a
> path as a name, so it exits 1 — *after* the rule's work is done. That throws
> in `persistence.finished()`, the job is marked failed, and Snakemake deletes
> its output. With a `directory()` output that means a completed 75-minute
> reference build gets erased. The rules activate the env in the shell block
> instead, which keeps Snakemake's conda bookkeeping out of it entirely.

### Resources

`config/cluster.yaml` sets `tmem` / `h_vmem` / `h_rt` / `submission_string` per
rule. **`tmem` and `h_vmem` are per slot** on this cluster, so `-pe smp 4` with
`tmem: 16G` reserves 64 GB total:

| rule | slots | tmem | total | h_rt |
|---|---|---|---|---|
| `mkref` | 4 | 16G | 64G | 24:00:00 |
| `split_pipe_pre` | 8 | 4G | 32G | 24:00:00 |
| `split_pipe_align` | 8 | 6G | 48G | 48:00:00 |
| `split_pipe_post` | 4 | 16G | 64G | 72:00:00 |
| `split_pipe` *(fallback)* | 4 | 16G | 64G | 120:00:00 |
| `split_pipe_comb` | 8 | 8G | 64G | 24:00:00 |
| `cellsnp` | 8 | 4G | 32G | 24:00:00 |
| `vireo` | 4 | 8G | 32G | 12:00:00 |
| `import_assignments` | 1 | 16G | 16G | 2:00:00 |
| `celltype_annotation` | 8 | 24G | 192G | 48:00:00 |
| `celltype_annotation_test` | 4 | 8G | 32G | 12:00:00 |

### Why split-pipe is three rules

In split-pipe the **pipeline steps are modes** (`run_modes.py:RUN_MODE_DEFS`):
`inqc pre align post mol dge sub ana rep clean`. `--mode all` just walks that
list, and the help documents running a subset as *"start at `<mode1>` and use
`--until_step <mode2>`"*. Profiling one sublibrary (`process/mem_profile.txt`,
ParseFANS1, 1.69e9 reads, 4 threads) showed a single job is the wrong shape:

| step | hours | mean cores | peak RSS |
|---|---|---|---|
| InQC | 0.03 | 0.86 | 11.2G |
| Preprocess | 5.73 | 2.77 | 11.3G |
| **Align (STAR)** | **9.71** | **3.87** | **17.9G** |
| PostProcess | 3.28 | 1.87 | 4.3G |
| Molinfo | 1.09+ | 1.20 | 3.7G |

Alignment is the only CPU-bound step and it saturates its 4 cores; everything
from `post` onward idles at 1.2--1.9 cores while still holding a 64 GB
reservation. `align` gets 8 slots rather than 16 deliberately -- 16 would need
hundreds of slots once several sublibraries run at once, and a job that queues
finishes later than a narrower one that runs. So `split_pipe_pre` (inqc+pre), `split_pipe_align` (align alone)
and `split_pipe_post` (post→clean) are submitted separately, each sized to what
it actually uses. Set `splitpipe.split_steps: false` to go back to one
`--mode all` job per sublibrary.

Two things make this safe:

- **Sentinel outputs.** `pre` and `align` declare `.spipe_*.done` files, not
  real ones. split-pipe's `clean` step deletes the `process/*.bam`
  intermediates at the end of the run; a rule declaring one of those as its
  output would see it missing on the next invocation and re-run the whole chain.
- **Every step gets the full flag block.** This run's `run_proc_def.json` has an
  empty `parent` object, so `init_parent()` (`spc_init.py:1280`) cannot recover
  chemistry/kit from a parent dir. `--chemistry`, `--kit`, `--genome_dir` and
  `--samp_list` are passed on every call; only `--mode`, `--until_step` and
  `--nthreads` differ.

**Walltimes are sized for the largest sublibrary, not the one that was tested.**
ParseFANS1 is the 3rd *smallest* of the 16 (102 GB of FASTQ vs a 122 GB mean and
ParseFANS11's 196 GB), so its timings scale by up to 1.92x — before any I/O
contention from running all 16 at once.

The slot counts are mirrored by the per-rule `threads:` block in
`config.yaml` — **keep the two in step**. A rule that asks snakemake for more
threads than it reserved slots gets throttled or killed for overrunning its
memory reservation.

`h_rt` values are quoted in the YAML on purpose: unquoted, YAML 1.1 reads
`72:00:00` as sexagesimal and passes qsub `259200`. SGE takes that as seconds so
it still works, but the quoted form is what's meant.

`cluster_qsub.sh` requests `tscratch=50G` for split-pipe's temporary files.

Note that `--cluster-config` is a Snakemake 7-and-earlier interface (removed in
8), which is what this submission pattern relies on.

`splitpipe_all` exists because `all` also depends on
`demux.provided.table` (the collaborator's assignment table). Snakemake builds
the whole DAG up front, so if that table is missing `all` fails with a
`MissingInputException` before any split-pipe job is submitted -- even though
nothing about split-pipe is blocked. `splitpipe_all` stops at the combined DGE
matrix, so the multi-day sublibrary run can start while the demux inputs are
still being sorted out.

Outputs land in `results/`:
- `splitpipe/combined/all-sample/DGE_filtered/` — combined DGE matrix. **Not
  `all-well`**: with samples given via `--samp_list` and `make_allwell` at its
  default of False, split-pipe never creates an `all-well` sample
  (`spc_setup.py:209` only adds one when no samples were given, or when
  `--yes_allwell` is passed). This run sets up four: `K_negative`,
  `TDP_negative`, `positive`, `all-sample`. `all-sample` covers all 96 wells,
  so it is the equivalent file.
- `demux/donor_assignments.all_sublibraries.tsv` — one row per cell
  (sublibrary, cell id, donor, doublet prob, joined donor metadata).
- `demux/demux_summary.tsv` — cells per donor × sublibrary.

## Before you scale to 16 sublibraries — verify these

1. **Chemistry.** `splitpipe.chemistry: auto` lets split-pipe detect the kit.
   Pin it to `v1`/`v2`/`v3` only if you know the version.
2. **BAM cell-barcode tag.** `samtools view <a_sublib_bam> | head -1` and
   confirm which tag holds the cell barcode and which the UMI. Set `cell_tag` /
   `umi_tag`. If cells have no UMI tag, set `umi_tag: "None"`.
3. **Barcode-column match.** Strings in `barcodes.tsv` (from `cell_metadata.csv`
   via `splitpipe.barcode_column`) must equal the BAM tag values exactly.
4. **Contig naming** — see Stage 0 above.

Validate the whole chain on **one** sublibrary before launching all 16.

## Demultiplexing — two modes

`config.yaml → demux.mode` picks the branch. Both write the same two files, so
everything downstream is identical.

- **`provided`** (current setting) — import an assignment table made elsewhere.
  Points at the collaborator's demuxalot + Vireo result: **250,219 barcodes,
  34 donors** (31 GSA-genotyped FTD donors + 3 Vireo virtual donors),
  16/16 sublibraries. Barcodes are `<bc_wells>__s<sublibrary>` (e.g.
  `01_01_67__sParseFANS1`), split back into sublibrary + barcode on import.
  cellsnp/vireo/genotype-panel rules are not even created in this mode.
- **`vireo`** — call donors here: merge per-donor VCFs into a panel, run
  cellsnp-lite per sublibrary, then vireo `-d`.

**There are no VCFs for this cohort on this system**, so `vireo` cannot run.
Selecting it fails immediately at parse time with an explanatory error rather
than part-way through a long run. Re-enable it by setting `donor_vcfs` (a TSV of
`donor_id` + per-donor VCF) or `reference_vcf` (one multi-sample VCF) in
`config.yaml`.

Two things to know about the provided table:

- 49,648 of the 250,219 rows are `dmux_type == DBL` / `likely_missing_donor`.
  **These are not doublets.** They are exactly the three Vireo virtual donors —
  `donor3` (18,602) + `donor33` (18,088) + `donor30` (12,958) = 49,648 — i.e.
  cells demuxalot called doublet because their true donor had no genotype in
  the panel, which Vireo then reconstructed as unnamed donors. Setting
  `drop_doublets: true` would delete three entire donors. Leave it `false`.
- Only 12 of the 34 donors have any `high_conf_donor` cells; the other 22 are
  entirely `only_demuxalot`, i.e. supported by one caller. Worth weighing before
  treating all 34 as equally solid.

## Donor tables — what we actually have

No VCFs, and no direct genotype data. The donor information available here is
whatever can be read out of the collaborator's assignment table.

**`config/donor_metadata.tsv`** — the working donor table, *derived* from that
assignment table: 34 donors with cell counts, sublibrary coverage, demux-group
and sort breakdowns, and `source` (`genotyped` = a pathology number,
`vireo_virtual` = Vireo reconstructed the donor but never matched an identity).
`nygc_id` / `path_diag` / `mutation` are deliberately **empty** — fill them in as
the information arrives. Regenerate the count columns any time with:

```bash
snakemake donor_manifest -c1   # -> results/demux/donor_manifest.tsv
```

Only the hand-maintained columns are joined onto cells; the derived counts are
always recomputed, so they cannot go stale.

**`config/donors_nygc_wgs.tsv`** — NYGC WGS metadata (23 donors: `path_diag`,
`mutation`, `PF-UCL-*` subject IDs). Not wired into the pipeline and cannot be:
these `CGND-HDA-*` IDs have **no crosswalk** to the pathology numbers the
demultiplexing used. Kept so the diagnosis information isn't lost. Its
provenance is unverified — check it against the original source metadata.

Known gaps:

- 34 donors appear in the data; 36 went into the experiment. Two are unaccounted
  for, on top of the 3 that are virtual.
- The single highest-value missing piece is a **CGND-HDA ↔ pathology-number
  crosswalk**. With it, `donors_nygc_wgs.tsv` folds into `donor_metadata.tsv`
  via `nygc_id` and every cell gains a diagnosis. Without it, no subtype
  analysis (TDP-A/B/C) is possible.
- A previous version of this config carried 13 placeholder donor rows
  (`noDNA-TDPA-01`, …) with invented IDs and fake VCF paths. They were not real
  records and have been removed; what they encoded was only a count —
  5 Control, 5 FTLD-TDPA, 2 FTLD-TDPB, 1 FTLD-TDPC had no DNA.
- The collaborator's best result combined **Vireo** (reference-free, n=36, then
  Plink IBS matching to the VCF) with **demuxalot** (reference-based). That
  result is what `demux.mode: provided` imports. demuxalot and the Vireo/IBS
  reconciliation are not yet rules here — only the reference-based Vireo branch
  is implemented natively.

## Stage 4 — cell type annotation

A separate target, not part of `rule all`. It needs the combined DGE
(`splitpipe/combined/all-sample/DGE_filtered/`) and
`demux/donor_assignments.all_sublibraries.tsv`, and writes to
`<results_dir>/celltype/`. It is the pipeline version of
`01_qc_harmony_celltyping.Rmd`, but runs on every cell instead of a 50k
subsample.

```bash
source submit.sh celltype_test   # 20k random cells, 32 GB job, output in celltype_test/
source submit.sh celltype        # all 426k cells, one 192 GB job
```

Everything is done by `workflow/scripts/celltype_annotation.R`, in one job:

1. **Donor join.** `s<N>` in `bc_wells` is mapped to a sublibrary name and
   checked against barcode overlap with the donor table.
2. **QC** on every cell: genes (200–10,000), transcripts (300–100,000) and
   % mitochondrial (≤ 5). Thresholds are in `celltype.qc`.
3. **Doublets.** scDblFinder runs per sublibrary on cells that pass the
   low-quality filters. The expected rate is set in `celltype.doublets`, because
   scDblFinder's own default is a 10x figure. The DBL calls in the donor table
   are *not* used (see Demultiplexing).
4. **Clustering.** Log-normalise, 2,000 HVGs, PCA, then Harmony on `donor_id`,
   UMAP and Louvain clusters at 0.5 (plus 0.2 and 1.0).
5. **Labels.** Azimuth (`humancortexref`) runs one sublibrary at a time, four in
   parallel, and gives per-cell class / subclass / cluster plus a mapping score.
   Each Seurat cluster also gets a majority Azimuth label and a marker-gene
   broad class.

| output | contents |
|---|---|
| `cell_annotations.csv` | **One row per split-pipe cell (426k).** Identifiers: `bc_wells`, `cell_id` (demux table form), `collab_barcode` (`<bc>__s<sublibrary>`), `barcode`, `sublibrary`, `mtx_row`, `bc1/2/3_well`. Also sort, donor, QC metrics, doublet score/class, `qc_pass`, `clustered`, `exclude_reason`, clusters, Azimuth `predicted.*`, `cluster_*` labels, UMAP coordinates, and the donor metadata. Cells that were not clustered have NA labels. |
| `cluster_summary.tsv` | Per cluster: size, dominant donor, majority Azimuth class/subclass, marker class and margin, top markers |
| `cluster_markers.tsv` | `FindAllMarkers`, up to 300 cells per cluster |
| `qc_summary.tsv` | Per sublibrary: cells failing each filter, how many pass |
| `celltype_report.pdf` | QC violins, doublet scores, PC R² before/after Harmony, UMAPs, marker dot plot, donor make-up per cluster |
| `seurat_harmony.rds` | The Seurat object (counts + data, reductions; no scale.data or graphs) |

**Why one big job.** Clustering is one in-memory computation over all cells,
so it cannot be split across jobs. On this cluster a single large `-R y`
reservation also starts sooner than a chain of mid-sized jobs that each wait in
the queue. The memory estimate is in `workflow/rules/celltype.smk`. The log
prints R's heap at every step and the peak at the end, so the 192 GB request
can be trimmed after the first full run.

**R environment.** By default Snakemake builds `workflow/envs/celltype.yaml`.
Its post-deploy script then installs Azimuth, SeuratData and `humancortexref`
from GitHub, which needs internet access on the login node. If an R env with
these packages already exists, set `conda_envs.celltype` to its path and the
rule activates it directly, like split-pipe.

## Not yet in the pipeline

Downstream steps the collaborator ran by hand, worth adding next:
- **regtools junction extraction** per sublibrary on UMI-deduplicated BAMs, to
  pull cryptic splicing events (69 junctions in 54 genes were recoverable).
  split-pipe's BAMs were not kept; `starsolo/` is a separate workflow that
  regenerates them from split-pipe's called cells (see `starsolo/README.md`).
- ~~Pegasus QC / harmony clustering~~ → done in Seurat as Stage 4 (`celltype`).
