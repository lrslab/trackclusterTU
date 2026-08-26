# Output directory layout

This page describes the files written by:

- `trackclustertu run` (full pipeline from FASTQ manifest to TU/gene outputs),
- `trackclustertu map` (FASTQ -> minimap2 -> sorted BAM -> BED manifests),
- `trackclustertu bam-to-bed` (BAM -> BED6), and
- `trackclustertu cluster` / `trackclustertu recount` (TU clustering and recounting), and
- `trackclustertu diagnose-missed-tus` / `trackclustertu rescue-missed-tus` (post-clustering missed-TU review and rescue)

## Common outputs (from `trackclustertu cluster`)

Given an input BED/TSV (usually a sorted BED6 track), `trackclustertu cluster` writes:

- `*.tus.bed` (BED6): one line per TU interval
  - columns: `chrom  start  end  tu_id  0  strand`
  - default IDs encode contig bytes, coordinates, and strand, so filtering another TU does not renumber this call
  - coordinate-identical final consensus families are coalesced before output, and stable-ID uniqueness is enforced
- `*.tu_id_map.tsv`: emitted, stable, and sequential compatibility IDs with genomic identity
  - pass `--tu-id-style sequential` only when an older downstream workflow requires serial `TU000001` names
- `*.membership.tsv` (TSV): one line per read assignment
  - schema marker: `#trackclustertu_membership_schema=v2`
  - cluster writes the canonical 20-column v2 form
  - the first four columns remain `read_id  tu_id  score1  score2`; `score1` is `overlap / union` and `score2` is `overlap / max(length)`
  - appended fields record assignment status, best/second candidates and scores, both endpoint deltas, assignment-score margin, primary/fractional weights, and full-length-evidence state
- `*.tu_count.csv` (CSV): TU counts table
  - header: `tu_id,count,unique_count,full_length_evidence_count,total_count,fractional_count`
  - `count` is a compatibility alias for `total_count`
- `*.tu_endpoint_stats.tsv` (TSV): auditable consensus evidence
  - support, strand-aware consensus endpoints, endpoint min/max/spreads, and the longest example molecule
- `read_rejections.tsv` (TSV, schema `trackclustertu.read-rejections.v1`): every read excluded during parsing, validation, deduplication, or length filtering
  - columns: `sample  source_path  input_format  record  line  read_id  stage  reason  detail`
  - the file is written atomically and remains header-only when every read is accepted; override it with `--out-read-rejections`

By default, `trackclustertu cluster` forms `score1` seed clusters and then runs a second-pass `score2` attachment that pools truncated molecules with their parent cluster.
Set the thresholds with `--score1-threshold` and `--score2-threshold`.
The second-pass 3 prime gate is one-sided: a shorter cluster may terminate anywhere inside its parent but may not overhang the parent's strand-aware 3 prime end by more than `12 bp` by default. Configure the limit with `--three-prime-tolerance-bp`.
If you need to relax 5 prime fragmentation for near-matching reads, you can also set `--max-5p-delta`.
If `--skip-score2-attachment` is used, the `score1` seed clusters are kept as the final TUs.
Clustering is span-based: when the input is BED12, only the outer transcript interval `[tx_start, tx_end)` participates in `score1` and `score2`.

Individual read-record errors are recoverable by default: the bad read is
quarantined, its reason is recorded, and parsing continues with the next
record. Duplicate IDs quarantine every occurrence within the same input/sample
instead of selecting one by input order. An entirely rejected sample remains a
zero-count sample, and an entirely rejected input produces a valid empty output
set plus the rejection report. `--strict-read-errors` evaluates the input,
then converts any accumulated read rejection into a command error before the
output set is published; it is atomic failure-on-any-rejection, not a promise
to stop decoding at the first bad record. Input/manifest structure and I/O
failures, output failures, and cross-output invariants are never downgraded.

When `--manifest` is used, `trackclustertu cluster` can also write:

- `*.pooled.bed` (BED6): pooled reads used for clustering
  - read IDs are rewritten as `<sample>::<read_id>`
  - `::` is the reserved sample/read separator, so manifest sample names cannot contain it
- `*.tu_sample_long.tsv` (TSV): one row per non-zero `(tu_id, sample)` pair
  - columns: `tu_id  sample  unique_count  full_length_evidence_count  total_count  fractional_count`
- `*.tu_sample_matrix.tsv` (TSV): TU-by-sample count matrix
  - header: `tu_id` plus `<sample>.<metric>` columns in manifest order
- `*.tu_group_matrix.tsv` (TSV): TU-by-group count matrix
  - header: `tu_id` plus `<group>.<metric>` columns in first-seen manifest order
  - written only when the manifest has a non-empty `group` column

Assignment status and count behavior are explicit:

- `unique`: one hard TU; contributes one to `unique_count`, `total_count`, and `fractional_count`.
- `partial`: one qualifying contained/partial hard TU; contributes one to `total_count` and `fractional_count`, but not `unique_count`.
- `ambiguous`: no hard TU (`tu_id=.`). With `--fractional-assignment`, its exactly normalized weight is propagated to each candidate's `fractional_count`; without that flag it contributes to no TU count.
- `unassigned`: no hard or fractional contribution.
- `full_length_evidence_count` is incremented only for a hard-assigned read whose manifest points to a matching `trackclustertu.bam-evidence.v1` sidecar and whose `effective_full_length` value is true; raw BED input is never labeled full-length by assumption.

The same rules are applied independently to TU, sample, group, and downstream gene matrices. Therefore `fractional_count` always includes weight `1` for every hard assignment and, when fractional assignment is enabled, also includes ambiguous fractional weights. A sample/group `total_count` remains the hard-assignment count.

`--min-tu-count` is evaluated at different stages for the two relevant
workflows. `cluster` and `run` retain a TU when its pre-assignment clustering
family has at least that many members. `recount` has no clustering families, so
it retains a row when the summed final hard `total_count` reaches the threshold;
fractional-only support does not satisfy it.

If `--annotation-bed` is provided, it can also write:

- `*.tu_gene.tsv` (TSV, schema `trackclustertu_tu_gene_schema=v2`): qualifying same-strand and antisense TU×gene relationships
  - columns: `contig  strand  tu_id  tu_start  tu_end  relation  context_order  gene_id  gene_start  gene_end  overlap_bp  tu_fraction  gene_fraction`
  - `context_order` follows transcription direction (low-to-high coordinates on `+`, high-to-low on `-`)
- `*.tu_semantics.tsv` (TSV, schema `trackclustertu_tu_semantics_schema=v1`): one row per TU
  - includes the ordered same-strand gene-context signature, ordered antisense context, and fixed-order multi-label classifications
  - criteria are embedded as metadata: same-context shared 3-prime/different 5-prime for `alternative_start`; same-context shared 5-prime/different 3-prime for `alternative_termination`; two or more same-strand genes for `readthrough`; opposite-strand overlap for `antisense`; no qualifying relationship for `intergenic`; and strict same-context containment with both boundaries internal for `probable_processing_product`
- `*.tus.gff3` (GFF3): 1-based inclusive `transcript` features plus child `gene_overlap` relationship features linked with `Parent`
  - sequence IDs and attributes are UTF-8 byte percent-escaped where required by GFF3
- `*.tus.anchored.bed12` (BED12 + extra columns): TU blocks clipped to qualifying same-strand genes
  - extra column `name2`: comma-separated member read IDs, with `|<read_count>` suffix
  - extra column `gene_list`: comma-separated qualifying same-strand gene IDs
  - when there is no qualifying same-strand gene, the record has one full-TU block and `gene_list=.`; antisense-only relationships do not create blocks or gene-list entries
- `*.gene_count.csv` (CSV): gene counts table
  - same metric columns as TU counts: compatibility `count`, unique, full-length-evidence, hard total, and fractional
- `*.gene_sample_matrix.tsv` (TSV): gene-by-sample count matrix
  - header: `gene_id` plus `<sample>.<metric>` columns in manifest order
- `*.gene_group_matrix.tsv` (TSV): gene-by-group count matrix
  - header: `gene_id` plus `<group>.<metric>` columns in first-seen manifest order
  - written only when the manifest has a non-empty `group` column

The default relationship threshold is one overlapping base. Configure all three conjunctive gates with `--gene-min-overlap-bp`, `--gene-min-tu-fraction`, and `--gene-min-gene-fraction`; both fractions must be finite values in `[0,1]`.

Gene counts use `#count_semantics=nonexclusive_same_strand_tu_overlap;...` metadata. Every hard or fractional TU assignment is copied to every qualifying same-strand gene, while antisense relationships are not counted. A polycistronic TU can therefore increment multiple genes and summed gene totals can exceed read/TU totals.

## Missed-TU diagnostic outputs (from `trackclustertu diagnose-missed-tus`)

This stage is meant to run after clustering. Pass the BED6 track that was clustered:

- single-input workflow: the original BED6 file used with `trackclustertu cluster --in`
- manifest workflow: the pooled BED6 file, usually `results/pooled.bed`

It writes:

- `--out-tsv`: TSV report of candidate boundary modes not represented by the current TU set
- `--out-bed`: optional BED6 track for quick visualization in IGV or genome browsers
- `--out-read-rejections`: recoverable read errors (defaults to `diagnose_read_rejections.tsv` beside `--out-tsv`)

Diagnostic candidate IDs use `MISS####`. They identify rows in this report and
its optional BED track; they are not TU IDs. `family_support` is the size of the
3-prime family, while `mode_support` is the number of detector-supporting reads
for that 5-prime mode. Neither value is a final rescued hard-assignment count.

Diagnosis and rescue are independent commands over the same clustering inputs.
`rescue-missed-tus` does not read either diagnostic output; it reruns the same
detector and then performs rescue-specific deduplication, support filtering,
and read competition.

The current default detector settings are:

- `--three-prime-window-bp 12`
- `--max-three-prime-family-diameter-bp 12` (when omitted, follows `--three-prime-window-bp`)
- `--five-prime-window-bp 10`
- `--min-family-support 20`
- `--min-mode-support 20`
- `--min-mode-fraction 0.02`
- `--max-candidates-per-family 3`
- `--min-read-len` unset

## Missed-TU rescue outputs (from `trackclustertu rescue-missed-tus`)

This stage uses the same clustered BED6 input plus the original TU BED and membership TSV. It writes rescued TU calls without modifying the original clustering outputs in place.

It writes:

- `--out-tu`: rescued TU BED6
- `--out-membership`: rescued membership TSV
- `--out-tu-count`: optional rescued hard-count CSV with header `tu_id,count`
- `--out-candidates-tsv`: optional TSV report of rescued candidate modes
- `--out-candidates-bed`: optional BED6 track of rescued candidate modes
- `--out-read-rejections`: recoverable read errors (defaults to `rescue_read_rejections.tsv` beside `--out-tu`)

Rescue first deduplicates nearby detected modes. Each retained read selects its
best deduplicated candidate by smaller strand-aware 3-prime delta, smaller
5-prime delta, higher mode support, and finally deterministic candidate order.
Modes whose provisional winning-read count is below `--min-mode-support` are
removed; affected reads are not reconsidered against a next-best mode. A
surviving selected candidate is then compared only with the read's current hard
TU assignment using the same ordering; the hard TU's comparison support is
`0`, so equal boundary deltas favor the supported mode. Modes are filtered
again on their final promoted-read count. A read with no hard assignment has no
existing comparison target and can be promoted if its candidate survives that
final filter. Ambiguous candidates or fractional assignments are not treated
as a current hard TU.

Rescued TU IDs use `RESC####` by default; change the prefix with
`--rescue-prefix`. Candidate TSV/BED reports still use report-only `MISS####`
IDs. In the rescue candidate TSV, `rescued_tu_id` gives the direct mapping and
`rescued_hard_count` gives the final number of promoted hard assignments; the
BED record name contains both the `MISS####` report ID and rescued TU ID.
`mode_support` remains detector support before final hard-assignment
competition and need not equal `rescued_hard_count`.

The rescued membership follows the input membership schema. Four-column legacy
input remains legacy. Canonical v2, the backward-compatible earlier six-column
v2 form, and mixed legacy/v2 input are accepted. For retained reads, rows not
promoted by rescue preserve their representation, assignment status, candidate
rankings, ambiguity/fractional weights, and full-length-evidence flag. A
promoted v2 row is rewritten as a canonical 20-column v2 `unique` hard
assignment to its rescued TU; its prior candidate ranking is replaced because
rescue has adjudicated that assignment, while its full-length-evidence flag is
retained. A homogeneous v2-bearing output declares its one layout with
`#columns`; a homogeneous legacy-only output does not gain a synthesized schema
header. A heterogeneous output declares the widths actually present with
`#row_widths=4,6,20` (omitting absent widths) and the applicable
`#legacy_columns`, `#compatible_v2_columns`, and `#canonical_v2_columns`
declarations. `#contains_legacy_rows=true` and
`#contains_compatible_v2_rows=true` explicitly mark those compatibility rows.
For v2-bearing output, rescue also writes
`#rescue_unmodified_rows=preserved_from_input`,
`#rescue_promoted_rows=canonical_unique_hard_assignment`, and
`#rescue_count_semantics=hard_assignments`.

The existing TU BED must have unique, non-empty IDs. Membership read IDs must
be unique and must occur in the BED6 input; every hard, best, second, and
fractional TU reference is a foreign key into the existing TU BED. Violations
are fatal. Invalid BED reads and every copy of a duplicate BED read ID are
quarantined as described in `--out-read-rejections`; membership rows are copied
only for reads retained after parsing, validation, deduplication, and optional
length filtering.

The rescued TU BED is the complete existing TU set plus rescued TUs. Existing
TUs are not removed merely because they have no final hard assignment. The
optional `--out-tu-count` file remains a hard-assignment count table, so such
TUs are present with count `0`; ambiguous fractional contributions are retained
in membership v2 and are included when that file is passed to `recount`.
Unlike cluster/recount `tu_count.csv`, which has six metric columns, this rescue
count sidecar deliberately has only `tu_id,count`.

`trackclustertu rescue-missed-tus` uses the same default detector settings as `diagnose-missed-tus`:

- `--three-prime-window-bp 12`
- `--max-three-prime-family-diameter-bp 12` (when omitted, follows `--three-prime-window-bp`)
- `--five-prime-window-bp 10`
- `--min-family-support 20`
- `--min-mode-support 20`
- `--min-mode-fraction 0.02`
- `--max-candidates-per-family 3`
- `--min-read-len` unset
- `--rescue-prefix RESC` (rescue only)

If you need updated sample/group count matrices after rescue, rerun `trackclustertu recount` with the rescued membership TSV.

## Count-only outputs (from `trackclustertu recount`)

`trackclustertu recount` requires both `--manifest` and `--pooled-membership`.
It only writes TU count tables and does not write TU interval, membership, pooled BED, annotation, or gene outputs.

Recount builds its TU universe from the union of hard, best, second, and
fractional TU references in the membership TSV; it has no TU BED input. Thus a
TU with zero hard count can still appear when a candidate/fractional field
references it, but an existing TU absent from every membership field cannot
appear in recount outputs even if it is retained in `rescued.tus.bed` and the
two-column rescue count sidecar.

Recount assigns all of the following output paths automatically. `--out-dir`
places defaults in the named directory; if it is omitted, the manifest-derived
default directory described below is used. Individual `--out-*` options
override their corresponding paths, so no explicit output option is required.

- `tu_count.csv`
- `tu_sample_long.tsv`
- `tu_sample_matrix.tsv`
- `tu_group_matrix.tsv` when the manifest has a non-empty `group` column

For four-column legacy membership, each hard-assignment row is projected as
`count=total_count=1`, `fractional_count=1`, `unique_count=0`, and
`full_length_evidence_count=0`. The legacy schema cannot establish uniqueness
or full-length evidence.

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
- `run_manifest.json`: atomic provenance record with effective mapping configuration, package/git state, canonical input SHA-256 values, library profile, and resolved external tool versions

`bam-to-bed --out-evidence <path>` (or manifest mode `--emit-evidence`) adds a `trackclustertu.bam-evidence.v1` TSV. It contains every BAM record in BAM order with MAPQ, full CIGAR, strand-aware 5-prime/3-prime soft clipping, retention/filter reason, library profile, and optional poly(A)/adapter/full-length/library-preparation evidence. Manifest mode also places each sidecar path in the `evidence` column of `samples.bed.tsv`; clustering validates retained sidecar rows against BED identity and boundaries before carrying `effective_full_length` into membership and count outputs. A complete length-delimited BAM block with an invalid record payload is represented by a `malformed_record` placeholder row and conversion continues; truncation, BGZF/I/O corruption, and invalid coordinate ordering remain fatal. Conversion prints deterministic category counts to stderr. Coordinate-sorted BAMs stream with actual order validation; `--min-boundary-support > 1` uses a bounded two-pass count and compares stable SHA-256 fingerprints of both decoded streams before publication. Inputs without `@HD SO:coordinate` use an explicitly warned in-memory sorting fallback.

## Full-pipeline extras (from `trackclustertu run`)

The `run` subcommand writes the same mapping and clustering outputs as running `trackclustertu map` followed by `trackclustertu cluster`, and it may also write:

- `annotation.bed`: converted annotation BED when `--annotation-gff` is used
- `run_manifest.json`: published only after the complete `run` pipeline succeeds

`trackclustertu run` also forwards the clustering controls used by `trackclustertu cluster`, including `--score1-threshold`, `--score2-threshold`, `--three-prime-tolerance-bp`, `--max-5p-delta`, and `--skip-score2-attachment`.
It does not automatically run the missed-TU diagnose/rescue stages.

Mapping always supplies minimap2 with `-ax map-ont`. It then appends repeatable `--minimap2-arg <OS_VALUE>` options without whitespace reparsing and accepts explicit `--minimap2` / `--samtools` paths. The old `--minimap2-args "..."` form appends its whitespace-split values, is deprecated, and cannot be combined with the repeatable form.

## Formats / conventions

- Coordinates are **0-based, half-open** intervals: `[start, end)` (BED style).
- BED6 outputs such as `bed/<sample>.bed` and `pooled.bed` use:
  `chrom  start  end  name  score  strand`
- Cluster-generated membership TSV keeps `read_id  tu_id  score1  score2` as its first four columns and appends the v2 audit fields described above. Rescue may preserve legacy rows, earlier six-column v2 rows, or mixed compatibility rows from its input; consult that file's metadata and per-row schema field.
- In pooled mode, membership read IDs are tagged as `<sample>::<read_id>`; `::` is reserved and is rejected in manifest sample names

`cluster --format auto` is the default. Suffix matching is case-insensitive:
`.bed12` selects BED12, `.bed` selects BED6, and `.tsv` or `.txt` selects TSV.
FASTQ suffixes (`.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz`) are recognized so
that `cluster` can direct the user to `map` or `run`, but FASTQ is not a direct
cluster input. Unrecognized suffixes fall back to BED6. In manifest mode,
recognized suffixes must all infer the same format; conflicting formats fail.
Passing `--format` explicitly applies one parser to every manifest row rather
than enabling heterogeneous per-row formats.

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
  --score2-threshold 0.80 \
  --three-prime-tolerance-bp 12 \
  --max-5p-delta 50 \
  --annotation-bed gene.bed \
  --out-dir results
```

## Recommended default rescue pass

After clustering, the default post-clustering review/rescue flow is:

```bash
trackclustertu diagnose-missed-tus \
  --in results/pooled.bed \
  --existing-tu results/tus.bed \
  --annotation-bed gene.bed \
  --out-tsv rescue/missed_candidates.tsv \
  --out-bed rescue/missed_candidates.bed
```

```bash
trackclustertu rescue-missed-tus \
  --in results/pooled.bed \
  --existing-tu results/tus.bed \
  --existing-membership results/membership.tsv \
  --annotation-bed gene.bed \
  --out-tu rescue/rescued.tus.bed \
  --out-membership rescue/rescued.membership.tsv \
  --out-tu-count rescue/rescued.tu_count.csv \
  --out-candidates-tsv rescue/rescued_candidates.tsv \
  --out-candidates-bed rescue/rescued_candidates.bed
```

If you started from a single BED input instead of a manifest, replace `results/pooled.bed` with the original BED6 reads file that you clustered.

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
