# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project follows SemVer.

## [0.1.3] - 2026-04-08

### Added

- Added `--max-5p-delta` to `trackclustertu cluster` and `trackclustertu run` to optionally relax strand-aware 5 prime fragmentation during second-pass TU attachment.

### Changed

- Updated second-pass TU attachment docs and regression tests for the new 3 prime tolerance and optional 5 prime override behavior.
- Removed release-noise files that should not be kept in the GitHub tree or release package.

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
