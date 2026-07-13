# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project follows SemVer.

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
  `--strict-read-errors` for fail-fast quality-control workflows.
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
  map-ont` SAM preset. Mapper and sorter threads share the requested thread
  budget, and effective values are recorded in the run manifest.
- BED manifests now retain emitted BAM evidence sidecars. Clustering validates
  each retained row against the BED read and propagates effective full-length
  evidence into membership, TU, sample, group, and downstream gene counts.
- Run manifests use Git metadata embedded at compile time; GitHub release builds
  inject and verify the exact tagged commit instead of querying a build-runner
  path after the binary has been distributed.
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
  the release provides one optional `SHA256SUMS` manifest. Crates.io
  publication remains disabled.

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

[0.2.0]: https://github.com/lrslab/trackclusterTU/compare/v0.1.4...v0.2.0
[0.1.4]: https://github.com/lrslab/trackclusterTU/compare/v0.1.3...v0.1.4
[0.1.3]: https://github.com/lrslab/trackclusterTU/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/lrslab/trackclusterTU/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/lrslab/trackclusterTU/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/lrslab/trackclusterTU/releases/tag/v0.1.0
