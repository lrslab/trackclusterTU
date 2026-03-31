# Output directory layout

This page describes the files written by:

- `trackclustertu run` (full pipeline from FASTQ manifest to TU/gene outputs),
- `trackclustertu map` (FASTQ -> minimap2 -> sorted BAM -> BED manifests),
- `trackclustertu bam-to-bed` (BAM -> BED6), and
- `trackclustertu cluster` / `trackclustertu recount` (TU clustering and recounting)

## Common outputs (from `trackclustertu cluster`)

Given an input BED/TSV (usually a sorted BED6 track), `trackclustertu cluster` writes:

- `*.tus.bed` (BED6): one line per TU interval
  - columns: `chrom  start  end  tu_id  0  strand`
- `*.membership.tsv` (TSV): one line per read assignment
  - columns: `read_id  tu_id  score1_to_tu  score2_to_tu`
- `*.tu_count.csv` (CSV): TU counts table
  - header: `tu_id,count`

By default, `trackclustertu cluster` forms score1 seed clusters and then runs a second-pass score2 merge that penalizes large short/long differences.
If `--skip-score2-attachment` is used, the score1 seed clusters are kept as the final TUs.
Clustering is span-based: when the input is BED12, only the outer transcript interval `[tx_start, tx_end)` participates in `score1` and `score2`.

When `--manifest` is used, `trackclustertu cluster` can also write:

- `*.pooled.bed` (BED6): pooled reads used for clustering
  - read IDs are rewritten as `<sample>::<read_id>`
- `*.tu_sample_long.tsv` (TSV): one row per non-zero `(tu_id, sample)` pair
  - header: `tu_id  sample  count`
- `*.tu_sample_matrix.tsv` (TSV): TU-by-sample count matrix
  - header: `tu_id` plus samples in manifest order
- `*.tu_group_matrix.tsv` (TSV): TU-by-group count matrix
  - header: `tu_id` plus groups in first-seen manifest order
  - written only when the manifest has a non-empty `group` column

If `--annotation-bed` is provided, it can also write:

- `*.tu_gene.tsv` (TSV): TU×gene overlaps
  - header: `#contig  strand  tu_id  tu_start  tu_end  gene_id  gene_start  gene_end  overlap_bp`
- `*.tus.anchored.bed12` (BED12 + extra columns): TU “blocks” clipped to overlapping genes
  - extra column `name2`: comma-separated member read IDs, with `|<read_count>` suffix
  - extra column `gene_list`: comma-separated overlapping gene IDs (or `.`)
- `*.gene_count.csv` (CSV): gene counts table
  - header: `gene_id,count`
- `*.gene_sample_matrix.tsv` (TSV): gene-by-sample count matrix
  - header: `gene_id` plus samples in manifest order
- `*.gene_group_matrix.tsv` (TSV): gene-by-group count matrix
  - header: `gene_id` plus groups in first-seen manifest order
  - written only when the manifest has a non-empty `group` column

## Count-only outputs (from `trackclustertu recount`)

`trackclustertu recount` requires both `--manifest` and `--pooled-membership`.
It only writes TU count tables and does not write TU interval, membership, pooled BED, annotation, or gene outputs.

With `--out-dir`, the default outputs are:

- `tu_count.csv`
- `tu_sample_long.tsv`
- `tu_sample_matrix.tsv`
- `tu_group_matrix.tsv`

At least one of these count outputs must be requested, either explicitly or through `--out-dir`.

## Mapping outputs (from `trackclustertu map`)

The mapping stage writes:

- `bam/<sample>.sorted.bam`: sorted BAM for each FASTQ
- `bam/<sample>.sorted.bam.bai`: BAM index
- `logs/<sample>.minimap2.log`: minimap2 stderr
- `logs/<sample>.samtools_view.log`: `samtools view` stderr
- `logs/<sample>.samtools_sort.log`: `samtools sort` stderr
- `logs/<sample>.samtools_index.log`: `samtools index` stderr
- `samples.bam.tsv`: BAM manifest for later conversion/reuse
- `bed/<sample>.bed`: BED6 reads track used for clustering
- `samples.bed.tsv`: BED manifest for `trackclustertu cluster`

## Full-pipeline extras (from `trackclustertu run`)

The `run` subcommand writes the same mapping and clustering outputs as running `trackclustertu map` followed by `trackclustertu cluster`, and it may also write:

- `annotation.bed`: converted annotation BED when `--annotation-gff` is used

## Formats / conventions

- Coordinates are **0-based, half-open** intervals: `[start, end)` (BED style).
- BED6 outputs such as `bed/<sample>.bed` and `pooled.bed` use:
  `chrom  start  end  name  score  strand`
- Membership TSV columns: `read_id  tu_id  score1_to_tu  score2_to_tu`
- In pooled mode, membership read IDs are tagged as `<sample>::<read_id>`

## Default output directories

If `--out-dir` is omitted:

- single-input `trackclustertu cluster --in reads.bed ...` writes to `reads.trackclustertu/`
- manifest-based `trackclustertu cluster --manifest <manifest> ...` writes to `<manifest-stem>.trackclustertu/`
- `trackclustertu recount --manifest <manifest> ...` uses the same default directory naming

Examples:

- `samples.tsv` -> `samples.trackclustertu/`
- `samples.bed.tsv` -> `samples.bed.trackclustertu/`

## Re-running TU clustering

Once you have `samples.bed.tsv`, you can re-run clustering without re-mapping:

```bash
trackclustertu cluster \
  --manifest samples.bed.tsv \
  --format bed6 \
  --score1-threshold 0.95 \
  --score2-threshold 0.6 \
  --annotation-bed gene.bed \
  --out-dir results
```

## Multi-sample pooled clustering

Given a manifest TSV:

```text
sample  reads  group
sampleA sampleA.bed control
sampleB sampleB.bed treated
```

you can cluster once across all samples and then write shared TU and count outputs:

```bash
trackclustertu cluster \
  --manifest samples.tsv \
  --format bed6 \
  --out-tu shared.tus.bed \
  --out-membership pooled.membership.tsv \
  --out-pooled-reads pooled.bed \
  --out-tu-count shared.tu_count.csv \
  --out-tu-sample-count-long shared.tu_sample_long.tsv \
  --out-tu-sample-count-matrix shared.tu_sample_matrix.tsv \
  --out-tu-group-count-matrix shared.tu_group_matrix.tsv
```

If you already have `pooled.membership.tsv`, you can recompute the count tables without reclustering:

```bash
trackclustertu recount \
  --manifest samples.tsv \
  --pooled-membership pooled.membership.tsv \
  --out-tu-count shared.tu_count.csv \
  --out-tu-sample-count-matrix shared.tu_sample_matrix.tsv \
  --out-tu-group-count-matrix shared.tu_group_matrix.tsv
```
