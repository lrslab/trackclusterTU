# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project follows SemVer.

## [Unreleased]

## [0.2.2] - 2026-08-26

This patch restores truncated-molecule pooling that was lost in `0.2.0`, while
preserving the final-consensus guard against daisy-chain merging introduced in
the `0.2.x` series. It can materially change TU membership and counts at
affected loci, but it does not change CLI option names or output schemas.

### Fixed

- Restored the one-sided strand-aware 3-prime attachment gate used by `0.1.x`.
  A shorter cluster may terminate anywhere inside its parent's span; it is
  rejected only when it overhangs the parent's 3-prime end by more than
  `--three-prime-tolerance-bp`. This reconnects 3-prime-truncated degradation
  ladders that the symmetric `0.2.0`/`0.2.1` gate incorrectly split into
  low-support TUs.
- Final consensus families now retain contained fragments even when both
  similarity thresholds fail. Only anchor-quality members vote on TU
  boundaries, so absorbed fragments increase support without moving a
  boundary, extending a TU span, or bridging two endpoint modes.
- Read assignment applies the same containment rule. A uniquely explained
  contained fragment is reported as `partial` and contributes to hard and
  fractional counts, including when its original low-support family was
  removed by `--min-tu-count`.

### Changed

- Exact boundary pairs are processed by descending support, ensuring that the
  best-supported pair anchors the first consensus family rather than whichever
  read sorts first genomically.
- The consensus-retention and assignment jitter window is now
  `max(--three-prime-tolerance-bp,
  floor((1 - --span-jaccard-threshold) * read_length))`. The second-pass
  attachment gate itself continues to use the fixed CLI tolerance.
- Seed components are pooled before the final consensus split. The direct
  consensus invariant is enforced once on final families, so membership still
  cannot be justified solely through a chain of intermediate reads.
- `tu_endpoint_stats.tsv` support and endpoint ranges include all retained
  members, while consensus coordinates are determined only by anchor-quality
  boundary voters.

### Compatibility

- The primary CLI flags remain `--span-jaccard-threshold`,
  `--overlap-over-longer-threshold`, and
  `--skip-overlap-over-longer-attachment`. The historical
  `--score1-threshold`, `--score2-threshold`, and `--skip-score2-attachment`
  spellings remain visible compatibility aliases.
- Skipping overlap-over-longer attachment also disables contained-fragment
  absorption and partial assignment. Final direct-consensus refinement still
  applies to span-Jaccard seed components.
- Existing BED/TSV inputs and v2 output schemas remain compatible. Because
  assignments and candidate fields can change, rerun clustering and downstream
  counting rather than combining `0.2.2` outputs with older results.

### Validation

On the retained primary, non-spliced *E. coli* RNA002 replicate 3 benchmark
(6,486 reads, default clustering parameters):

| Version | TUs | Hard assignments | Ambiguous |
| --- | ---: | ---: | ---: |
| `0.1.4` | 913 | 6,486 (100%) | not represented by the v1 schema |
| `0.2.1` | 1,587 | 4,995 (77.0%) | 1,491 |
| `0.2.2` | **888** | **6,269 (96.7%)** | **217** |

The remaining `0.2.2` reads are explicit near-tie ambiguous assignments under
the default `--ambiguity-margin 0.02`; `--fractional-assignment` can retain
their count mass without forcing a hard biological choice. Shuffling input
order produced byte-identical `tus.bed` output.

**Full commit comparison:**
[`v0.2.1...v0.2.2`](https://github.com/lrslab/trackclusterTU/compare/v0.2.1...v0.2.2)

## [0.2.1] - 2026-07-14

### Changed

- Clustering now maintains endpoint consensus incrementally and read assignment
  uses per-reference, per-strand interval indexes. In the measured synthetic
  duplicate-heavy and fixed-21-endpoint-jitter dense-locus cases, doubling reads
  from 5,000 to 10,000 took about 2x; the separated-TU assignment benchmark also
  avoids the former all-reads-by-all-TUs scan. These cases do not establish a
  general worst-case linearity guarantee for every dense locus.
- Clustering runs with at least one retained read report structured clustering
  and assignment phase boundaries, counts, and completed-stage elapsed times
  even without `--timings`. Progress read counts are the reads retained after
  parsing, validation, deduplication, and length filtering; runs with no
  retained reads publish empty results without entering those two phases.

### Fixed

- Coordinate-identical final consensus families are coalesced before assignment
  and output, preventing duplicate stable TU IDs and artificial exact-tie
  ambiguity. Stable-ID uniqueness is also enforced before files are published.
- GFF3 annotation parsing now stops at the standard `##FASTA` directive, so
  Prokka files with embedded reference sequences are accepted. Before that
  directive, every feature row still requires nine columns. Consumed `gene`
  rows retain strict coordinate and strand validation, and percent escapes in
  parsed attribute key/value pairs are validated.
- Direct `gff-to-bed` conversion now reports the requested final BED path after
  its output transaction commits instead of exposing a temporary staging path.

## [0.2.0] - 2026-07-13

This release changes TU boundary inference, assignment semantics, identifiers,
and output schemas. Existing `0.1.x` workflows should review the breaking
changes before adopting these outputs.

### Breaking Changes

- Final TU coordinates now come from strand-aware endpoint support modes rather
  than the longest read. The longest molecule is retained only as an auditable
  example, so TU boundaries can differ from `0.1.x` on the same reads.
- Cluster membership now uses schema v2 with assignment status, best and
  runner-up evidence, endpoint deltas, ambiguity margins, and optional
  fractional weights. Count outputs distinguish unique, full-length-evidence,
  hard-total, and fractional counts; parsers written for the legacy four-column
  membership file must be updated.
- Coordinate-derived TU IDs are now the default and remain stable when other
  TUs are filtered. Use `--tu-id-style sequential` for historical serial IDs;
  `tu_id_map.tsv` records both identities.
- The Cargo package is now named `trackclustertu`. Invariant-bearing model
  fields are private, CLI implementation modules are private, and supported
  Rust callers must use checked constructors, getters, typed request APIs, and
  the crate-level CLI facade.
- Gene counts use nonexclusive qualifying same-strand overlaps: a readthrough
  TU can contribute to every overlapping gene. Antisense relationships are
  reported but are not added to gene counts.

### Added

- Added transactional output publication for `cluster`, `recount`,
  `diagnose-missed-tus`, `rescue-missed-tus`, `bam-to-bed`, and `gff-to-bed`,
  including input/output alias checks and rollback-safe temporary files.
- Added BAM evidence schema v1 with MAPQ, CIGAR, strand-aware terminal clipping,
  filter reasons, deterministic conversion counters, chemistry profiles, and
  optional poly(A), adapter, full-length, and library-preparation evidence.
- Added `tu_endpoint_stats.tsv` with endpoint support, consensus, ranges,
  spread, and example-molecule information for each TU.
- Added configurable gene-overlap policy, ordered bacterial gene-context
  signatures, multi-label TU semantics, `tu_semantics.tsv`, and `tus.gff3`.
- Added atomic `read_rejections.tsv` diagnostics with source, sample, and line
  context, stable reason codes, deterministic summaries, and
  `--strict-read-errors` for atomic failure-on-any-rejection quality-control workflows.
- Added atomic `run_manifest.json` provenance for `map` and `run`, including
  effective configuration, input SHA-256 values, package and Git state,
  chemistry, thread allocation, and external tool versions.
- Added collision-safe missed-TU rescue validation, typed assignment origins,
  bounded endpoint families, and indexed endpoint-mode detection.
- Added Rust 1.88 minimum-version coverage, strict rustdoc checks, endpoint-mode
  benchmarks, and scheduled dependency auditing.

### Changed

- Missed-TU rescue now preserves unpromoted membership rows for retained reads,
  including legacy, canonical v2, earlier six-column v2, and mixed inputs;
  emits canonical v2 assignments for promoted v2 reads; declares heterogeneous
  row widths in metadata; retains full-length evidence; maps `MISS` report IDs
  to rescued TU IDs and final hard counts; and keeps every existing TU even when
  its final hard count is zero.
- Span-Jaccard seed components are rechecked against the final endpoint
  consensus, preventing single-linkage daisy chains and low-support bridges.
- Second-pass attachment is symmetric across strands, uses absolute endpoint
  deltas, and never allows the 5-prime override to bypass the 3-prime tolerance.
- Canonical score flags are `--span-jaccard-threshold` and
  `--overlap-over-longer-threshold`; `--score1-threshold` and
  `--score2-threshold` remain compatibility aliases.
- Mapping accepts explicit minimap2 and samtools paths plus repeatable,
  boundary-preserving `--minimap2-arg` values appended after the required `-ax
  map-ont` SAM preset. Mapper and sorter allocations sum to the requested thread
  budget for `--threads >= 2`; at `--threads 1`, both receive their required
  one-thread minimum. `samtools view` remains a separate single-threaded
  process. Effective values are recorded in the run manifest.
- BED manifests now retain emitted BAM evidence sidecars. Clustering validates
  each retained row against the BED read and propagates effective full-length
  evidence into membership, TU, sample, group, and downstream gene counts.
- Run manifests use Git metadata embedded at compile time; GitHub release builds
  inject the exact tagged commit and verify that its revision bytes are embedded
  in each binary instead of querying a build-runner path after distribution.
- Coordinate-sorted BAM conversion streams output while validating actual
  record order. Exact-boundary support uses a bounded two-pass implementation;
  unsorted input uses an explicit warned fallback.
- Recoverable per-read BED6, BED12, TSV, and BAM decoder errors are quarantined
  while processing continues. Structural input, I/O, ordering, output, and
  consistency errors remain fatal.
- CSV and TSV application tables use schema-aware readers and writers that
  preserve identifiers and correctly handle delimiters, quotes, comments,
  CRLF, embedded newlines, and trailing empty fields.
- BED, BED12, GFF3, manifest, membership, and clustering validation is stricter
  and reports source context where available.
- GitHub tag `v0.2.0` publishes Cargo release-mode binaries for
  `x86_64-unknown-linux-musl`, `aarch64-unknown-linux-gnu`, and
  `aarch64-apple-darwin`. Each archive includes `LICENSE` and `README.md`, and
  the release provides one required `SHA256SUMS` manifest. Crates.io publication
  remains disabled.

### Fixed

- Restored the strict GitHub release gate on Rust 1.97 by updating source that
  triggered newly enabled Clippy lints, without changing CLI behavior.
- Pinned tagged release gates and binary builds to Rust 1.97.0 so rerunning the
  release uses the same compiler and lint set.

### Security

- Updated `anyhow`, `crossbeam-epoch`, and development `rand` to resolve the
  applicable RustSec vulnerability and unsoundness advisories.

## [0.1.4] - 2026-04-09

### Added

- Added `trackclustertu diagnose-missed-tus` to report high-support boundary modes that are not represented by the current TU calls.
- Added `trackclustertu rescue-missed-tus` to promote supported missed-TU boundary modes into rescued TU and membership outputs without overwriting the original clustering results.

### Changed

- Updated docs, examples, and regression tests for the post-clustering diagnose/rescue workflow.

## [0.1.3] - 2026-04-08

### Added

- Added `--max-5p-delta` to `trackclustertu cluster` and `trackclustertu run` to optionally relax strand-aware 5 prime fragmentation during second-pass TU attachment.

### Changed

- Updated the default clustering overlap-over-longer (legacy `score2`) threshold to `0.80`.
- Updated second-pass TU attachment docs, examples, and regression tests for the optional 5 prime override.

## [0.1.2] - 2026-04-02

### Added

- Added strand-aware 3 prime tolerance for second-pass TU attachment.

### Changed

- Clarified testing requirements for TU tolerance behavior.

## [0.1.1] - 2026-03-24

### Changed

- Updated the pure Rust `trackclustertu` CLI behavior, docs, tests, and release packaging for the next patch release.

## [0.1.0] - 2026-01-09

### Added

- Interval and interval-list operations (intersection/union).
- Similarity scoring (`score1`, `score2`) for intervals and transcript-like tracks.
- TU clustering (`cluster_tus`) with sweep-line candidate generation + deterministic output.
- `trackclustertu` CLI for clustering BED6/BED12/TSV inputs into TU BED6 + membership TSV.
- Criterion benchmarks (`cargo bench`) and a baseline performance note in `doc/performance.md`.

[Unreleased]: https://github.com/lrslab/trackclusterTU/compare/v0.2.2...HEAD
[0.2.2]: https://github.com/lrslab/trackclusterTU/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/lrslab/trackclusterTU/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/lrslab/trackclusterTU/compare/v0.1.4...v0.2.0
[0.1.4]: https://github.com/lrslab/trackclusterTU/compare/v0.1.3...v0.1.4
[0.1.3]: https://github.com/lrslab/trackclusterTU/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/lrslab/trackclusterTU/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/lrslab/trackclusterTU/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/lrslab/trackclusterTU/releases/tag/v0.1.0
