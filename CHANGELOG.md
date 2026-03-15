# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project follows SemVer.

## [0.1.0] - 2026-01-09

### Added

- Interval and interval-list operations (intersection/union).
- Similarity scoring (`score1`, `score2`) for intervals and transcript-like tracks.
- TU clustering (`cluster_tus`) with sweep-line candidate generation + deterministic output.
- `trackclustertu` CLI for clustering BED6/BED12/TSV inputs into TU BED6 + membership TSV.
- Criterion benchmarks (`cargo bench`) and a baseline performance note in `doc/performance.md`.
