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
  scripts/ extract_barcodes.py  aggregate_assignments.py  import_assignments.py
  envs/   bcftools.yaml  demux.yaml  py.yaml
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
- **Contig naming is decided here.** GENCODE gives `chr1,chr2,…,chrM`; Ensembl
  gives `1,2,…,MT`. Whatever you pick must match any donor VCFs, or set
  `rename_chrs`. Mismatch is the usual reason cellsnp finds ~0 cells.

## Tools

- **split-pipe** — Parse Biosciences' proprietary binary (v1.9.0 here), on
  `PATH`. Not conda-installable, so it has no conda env.
- **cellsnp-lite** + **vireo** — conda (`workflow/envs/demux.yaml`) or the
  [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/)
  container.
- **bcftools / samtools** — genotype panel prep and BAM indexing.

## Run — UCL CS cluster (SGE)

Follows the [frattalab/rna_seq_snakemake](https://github.com/frattalab/rna_seq_snakemake)
pattern: `submit.sh` is itself qsub-ed and runs snakemake, which fans each rule
out as its own `qsub` using `config/cluster.yaml`; `cluster_qsub.sh` is the
jobscript template.

```bash
conda activate /SAN/vyplab/FUS_RNA/conda_envs/spipe   # split-pipe on PATH
snakemake -n                    # dry run first, from the login node

qsub submit.sh                  # full run (rule all)
qsub submit.sh reference        # build the shared Parse reference only
qsub submit.sh all rerun        # tag the submission "all_rerun"
```

Each submission gets `submissions/<timestamp>_<name>/` holding the job logs plus
a copy of the `config.yaml` and `cluster.yaml` that produced it.

`split-pipe` must be on PATH **inside the jobs**, not just at the login node.
Both `submit.sh` and `cluster_qsub.sh` use `#$ -V`, so the environment active at
qsub time is exported — activate the spipe env before submitting (or uncomment
the `source activate` line in `submit.sh`).

### Resources

`config/cluster.yaml` sets `tmem` / `h_vmem` / `h_rt` / `submission_string` per
rule. **`tmem` and `h_vmem` are per slot** on this cluster, so `-pe smp 16` with
`tmem: 4G` reserves 64 GB total:

| rule | slots | tmem | total | h_rt |
|---|---|---|---|---|
| `mkref` | 16 | 4G | 64G | 24:00:00 |
| `split_pipe` | 16 | 4G | 64G | 72:00:00 |
| `split_pipe_comb` | 8 | 8G | 64G | 24:00:00 |
| `cellsnp` | 8 | 4G | 32G | 24:00:00 |
| `vireo` | 4 | 8G | 32G | 12:00:00 |
| `import_assignments` | 1 | 16G | 16G | 2:00:00 |

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

Outputs land in `results/`:
- `splitpipe/combined/all-well/DGE_filtered/` — combined DGE matrix.
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

## Not yet in the pipeline

Downstream steps the collaborator ran by hand, worth adding next:
- **regtools junction extraction** per sublibrary on UMI-deduplicated BAMs, to
  pull cryptic splicing events (69 junctions in 54 genes were recoverable).
- **Pegasus** QC / harmony clustering with `--demux`, using barcode-level
  metadata of donor × sort.
