# trackclusterTU

Fast interval similarity and scalable clustering for bacterial transcript units (TUs) from mapped long reads.

This repository ships the Rust `trackclustertu` CLI.

## Project Metadata

- Repository: [lrslab/trackclusterTU](https://github.com/lrslab/trackclusterTU)
- License: MIT

## What Is A TU Here?

For bacteria, each mapped read is treated primarily as a single genomic interval `[start, end)` using 0-based, half-open coordinates.
Reads are clustered into candidate transcript units using overlap-based similarity.
When the input is BED12, clustering intentionally uses only the outer transcript span `[tx_start, tx_end)`; exon blocks do not affect `score1` or `score2`.

### Similarity Metrics

- `score1` (Jaccard-like): `overlap / union`
- `score2` (length-penalized overlap): `overlap / max(lenA, lenB)`

By default, `trackclustertu cluster` and `trackclustertu run` use `score1` to form seed clusters and then use `score2` in a second pass to merge only very similar seed clusters without treating strong short/long containment as a perfect match.
That second pass allows a strand-aware 3 prime mismatch of up to `12 bp` by default; adjust it with `--three-prime-tolerance-bp`.
If you need to relax 5 prime fragmentation for near-matching reads, you can also set `--max-5p-delta` to allow merges within an explicit strand-aware 5 prime delta.
If you want to keep the score1 seed clusters as the final TUs, pass `--skip-score2-attachment`.

## Install

### Option 1: Download A GitHub Release Binary

Prebuilt release tarballs are published from GitHub Actions on tagged releases:

- Releases: [lrslab/trackclusterTU/releases](https://github.com/lrslab/trackclusterTU/releases)
- Artifact naming: `trackclustertu-<tag>-<target>.tar.gz`

Download the archive for your platform, extract it, and place `trackclustertu` somewhere on your `PATH`.

### Option 2: Clone And Build From Source

```bash
git clone https://github.com/lrslab/trackclusterTU.git
cd trackclusterTU
cargo build --release --bin trackclustertu
```

The built executable will be:

- `target/release/trackclustertu`

To install it into Cargo's bin directory instead:

```bash
cargo install --path . --bin trackclustertu
```

## External Mapping Tools

`trackclustertu map` and `trackclustertu run` require these tools to be installed and available on `PATH`:

- `minimap2`
- `samtools`

Tested versions:

- `minimap2 2.30-r1287`
- `samtools 1.22.1`

If you only use `bam-to-bed`, `gff-to-bed`, `cluster`, or `recount`, these external mapping tools are not required.

## Usage

Help:

```bash
trackclustertu --help
trackclustertu run --help
trackclustertu map --help
trackclustertu cluster --help
trackclustertu recount --help
trackclustertu diagnose-missed-tus --help
trackclustertu rescue-missed-tus --help
trackclustertu bam-to-bed --help
trackclustertu gff-to-bed --help
```

Docs and examples:

- `doc/README.md`
- `examples/README.md`

Quick examples:

```bash
trackclustertu cluster \
  --in reads.bed \
  --format bed6 \
  --out-dir results
```

```bash
trackclustertu cluster \
  --manifest samples.bed.tsv \
  --format bed6 \
  --annotation-bed genes.bed \
  --out-dir results
```

```bash
trackclustertu recount \
  --manifest samples.bed.tsv \
  --pooled-membership results/membership.tsv \
  --out-dir recount_results
```

```bash
trackclustertu gff-to-bed \
  --annotation-gff genes.gff3 \
  --out-bed genes.bed
```

```bash
trackclustertu bam-to-bed \
  --in-bam sample.sorted.bam \
  --out-bed sample.bed
```

```bash
trackclustertu diagnose-missed-tus \
  --in sample.bed \
  --existing-tu results/tus.bed \
  --annotation-bed genes.bed \
  --out-tsv results/missed_tus.tsv \
  --out-bed results/missed_tus.bed
```

`diagnose-missed-tus` groups reads into strand-aware 3 prime end families, finds high-support 5 prime modes within each family, and reports candidate intervals that are not represented by the current TU BED.

```bash
trackclustertu rescue-missed-tus \
  --in sample.bed \
  --existing-tu results/tus.bed \
  --existing-membership results/membership.tsv \
  --annotation-bed genes.bed \
  --out-tu rescue/rescued.tus.bed \
  --out-membership rescue/rescued.membership.tsv \
  --out-tu-count rescue/rescued.tu_count.csv \
  --out-candidates-tsv rescue/rescued_candidates.tsv
```

`rescue-missed-tus` uses the same boundary-mode detector to add rescued TU intervals for abundant unsupported end modes and rewrites membership assignments only for reads that match a rescued mode better than their current TU.

## Recommended Default Workflow: Cluster Then Rescue

`trackclustertu run` stops after the main clustering and count outputs. The missed-TU rescue stage is a separate post-clustering pass.

The current built-in defaults in code are:

- clustering: `--score1-threshold 0.95`, `--score2-threshold 0.80`, `--three-prime-tolerance-bp 12`
- diagnose/rescue: `--three-prime-window-bp 12`, `--five-prime-window-bp 10`, `--min-family-support 20`, `--min-mode-support 20`, `--min-mode-fraction 0.02`

For a single BED6 input, the recommended default walkthrough is:

```bash
trackclustertu cluster \
  --in reads.bed \
  --format bed6 \
  --annotation-bed genes.bed \
  --out-dir results
```

```bash
trackclustertu diagnose-missed-tus \
  --in reads.bed \
  --existing-tu results/tus.bed \
  --annotation-bed genes.bed \
  --out-tsv rescue/missed_candidates.tsv \
  --out-bed rescue/missed_candidates.bed
```

```bash
trackclustertu rescue-missed-tus \
  --in reads.bed \
  --existing-tu results/tus.bed \
  --existing-membership results/membership.tsv \
  --annotation-bed genes.bed \
  --out-tu rescue/rescued.tus.bed \
  --out-membership rescue/rescued.membership.tsv \
  --out-tu-count rescue/rescued.tu_count.csv \
  --out-candidates-tsv rescue/rescued_candidates.tsv \
  --out-candidates-bed rescue/rescued_candidates.bed
```

For manifest-based clustering, use `results/pooled.bed` as the `--in` file for `diagnose-missed-tus` and `rescue-missed-tus`, because that is the BED6 track that was actually clustered.

If you need updated sample/group count tables after rescue, rerun `recount` from the rescued membership TSV:

```bash
trackclustertu recount \
  --manifest mapped/samples.bed.tsv \
  --pooled-membership rescue/rescued.membership.tsv \
  --out-dir rescue/recount
```

## Full Pipeline

Supported workflow:

`FASTQ -> sorted BAM -> BED6 -> TU clustering -> TU/gene counts`

### Step 1: Prepare A Sample Manifest

```tsv
sample	group	reads
sampleA	control	data/sampleA.fastq.gz
sampleB	treated	data/sampleB.fastq.gz
```

- `sample` is required
- `reads` is required
- `group` is optional

Relative `reads` paths are resolved relative to the manifest file.
Sample names must remain unique after replacing non-filename characters with `_`, because output BAM/BED/log filenames are derived from them.

### Step 2: Run The Whole Pipeline

```bash
trackclustertu run \
  --manifest samples.fastq.tsv \
  --reference-fasta ref.fa \
  --annotation-gff genes.gff3 \
  --out-dir results
```

`trackclustertu run` accepts the same clustering controls as `trackclustertu cluster`, including `--score1-threshold`, `--score2-threshold`, `--three-prime-tolerance-bp`, `--max-5p-delta`, and `--skip-score2-attachment`.
It does not automatically run `diagnose-missed-tus` or `rescue-missed-tus`.

This writes:

- `results/bam/<sample>.sorted.bam`
- `results/bam/<sample>.sorted.bam.bai`
- `results/bed/<sample>.bed`
- `results/logs/<sample>.*.log`
- `results/samples.bam.tsv`
- `results/samples.bed.tsv`
- `results/annotation.bed` when `--annotation-gff` is used
- `results/tus.bed`
- `results/membership.tsv`
- `results/pooled.bed` in manifest mode
- `results/tu_count.csv`
- `results/tu_sample_long.tsv` in manifest mode
- `results/tu_sample_matrix.tsv` in manifest mode
- `results/tu_group_matrix.tsv` in manifest mode when groups exist
- `results/gene_count.csv` with annotation
- `results/gene_sample_matrix.tsv` in manifest mode with annotation
- `results/gene_group_matrix.tsv` in manifest mode with annotation when groups exist

### Step 3: Optionally Run Stages Separately

You do not need to run these commands after Step 2.
Use this staged workflow only if you want to run each step separately instead of `trackclustertu run`.

```bash
trackclustertu map \
  --manifest samples.fastq.tsv \
  --reference-fasta ref.fa \
  --out-dir mapped
```

```bash
trackclustertu bam-to-bed \
  --manifest mapped/samples.bam.tsv \
  --out-dir mapped
```

```bash
trackclustertu gff-to-bed \
  --annotation-gff genes.gff3 \
  --out-bed genes.bed
```

```bash
trackclustertu cluster \
  --manifest mapped/samples.bed.tsv \
  --format bed6 \
  --score1-threshold 0.95 \
  --score2-threshold 0.80 \
  --three-prime-tolerance-bp 12 \
  --max-5p-delta 50 \
  --annotation-bed genes.bed \
  --out-dir results
```

The default clustering thresholds are `--score1-threshold 0.95` and `--score2-threshold 0.80`.
The default second-pass 3 prime allowance is `--three-prime-tolerance-bp 12`.
`--max-5p-delta` is optional and disabled unless you set it.

```bash
trackclustertu recount \
  --manifest mapped/samples.bed.tsv \
  --pooled-membership results/membership.tsv \
  --out-dir recount_results
```

If you want the default rescue pass after clustering, use:

```bash
trackclustertu diagnose-missed-tus \
  --in results/pooled.bed \
  --existing-tu results/tus.bed \
  --annotation-bed genes.bed \
  --out-tsv rescue/missed_candidates.tsv \
  --out-bed rescue/missed_candidates.bed
```

```bash
trackclustertu rescue-missed-tus \
  --in results/pooled.bed \
  --existing-tu results/tus.bed \
  --existing-membership results/membership.tsv \
  --annotation-bed genes.bed \
  --out-tu rescue/rescued.tus.bed \
  --out-membership rescue/rescued.membership.tsv \
  --out-tu-count rescue/rescued.tu_count.csv \
  --out-candidates-tsv rescue/rescued_candidates.tsv \
  --out-candidates-bed rescue/rescued_candidates.bed
```

For updated sample/group matrices after rescue, recount from `rescue/rescued.membership.tsv`.

## Main Outputs

- `tus.bed`
- `membership.tsv`
- `pooled.bed` in manifest mode
- `tu_count.csv`
- `tu_sample_long.tsv` in manifest mode
- `tu_sample_matrix.tsv` in manifest mode
- `tu_group_matrix.tsv` in manifest mode when groups exist
- `tu_gene.tsv` with annotation
- `tus.anchored.bed12` with annotation
- `gene_count.csv` with annotation
- `gene_sample_matrix.tsv` in manifest mode with annotation
- `gene_group_matrix.tsv` in manifest mode with annotation when groups exist

## Input Formats

- `bed6`: `chrom  start  end  name  score  strand`
- `bed12`: transcript span `[tx_start, tx_end)` is used for TU clustering; BED12 block structure is ignored by clustering
- `tsv`: `contig  start  end  id  strand`
