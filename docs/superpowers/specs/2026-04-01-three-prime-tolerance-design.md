# Strand-Aware 3 Prime Tolerance For Second-Pass TU Attachment

## Summary

The current second-pass TU attachment step allows a shorter seed cluster to attach to a longer seed cluster only when:

- the candidate TU is longer than the child TU
- the candidate TU satisfies `score2 >= score2_threshold`
- the candidate TU fully reaches the child's downstream end

That last requirement is too strict for biologically common alternative-end cases. A representative example on the `+` strand is a longer TU with an earlier start but a slightly earlier end than a shorter TU. These pairs can have strong overlap and pass `score2`, yet they stay split because the current logic requires `candidate.end >= child.end`.

This change replaces the hard downstream-end containment check with a strand-aware 3 prime tolerance measured in base pairs. The tolerance will be exposed on the CLI and default to `12`.

## Goals

- Merge near-identical seed clusters that differ only by a small strand-aware 3 prime offset.
- Preserve the current seed clustering and second-pass `score2` semantics.
- Make the tolerance explicit and user-configurable from the CLI.
- Keep the default behavior bounded so clearly distinct TU ends still remain separate.

## Non-Goals

- Changing the `score1` seed-clustering pass.
- Changing how the best parent candidate is chosen after eligibility filtering.
- Introducing exon-aware or transcript-structure-aware clustering.
- Changing behavior when `--skip-score2-attachment` is used.

## User-Facing Behavior

### New CLI Option

Add a new option to both `trackclustertu cluster` and `trackclustertu run`:

- `--three-prime-tolerance-bp <INT>`

Default:

- `12`

Meaning:

- During second-pass score2 attachment, allow a candidate TU to miss the child's 3 prime boundary by up to this many base pairs, using strand-aware coordinates.

### Attachment Rule

The current second-pass eligibility requires:

- candidate length strictly greater than child length
- candidate interval reaches at least the same downstream end as the child
- `score2(child, candidate) >= score2_threshold`

The revised rule will require:

- candidate length strictly greater than child length
- candidate interval is within `three_prime_tolerance_bp` of the child's 3 prime boundary
- `score2(child, candidate) >= score2_threshold`

Strand-aware interpretation of the 3 prime boundary:

- `+` strand: 3 prime boundary is the interval end
- `-` strand: 3 prime boundary is the interval start

Eligibility formulas:

- `+` strand: `candidate.end + tolerance >= child.end`
- `-` strand: `candidate.start <= child.start + tolerance`

This keeps the existing "shorter to longer" directionality while allowing a bounded downstream mismatch.

## Implementation Design

### Clustering Options

Extend `TuClusteringOptions` with:

- `three_prime_tolerance_bp: u32`

Default values:

- `attach_contained_reads: true`
- `three_prime_tolerance_bp: 12`

This keeps the tolerance visible in library code and avoids introducing a hidden magic constant.

### Second-Pass Eligibility Helper

Add a small helper in `src/tu/mod.rs` for the strand-aware 3 prime tolerance check. The helper should:

- accept the child and candidate `ReadRecord`s
- read the strand from the records
- return whether the candidate is within the configured 3 prime tolerance

Behavior:

- `+` strand compares interval ends
- `-` strand compares interval starts
- zero tolerance reproduces the current hard-boundary behavior

The second-pass loop should continue to reject:

- equal-length candidates
- shorter candidates
- candidates failing `score2`

The only eligibility change is replacing the hard end-containment test with the helper.

### CLI Plumbing

Add `three_prime_tolerance_bp: u32` with default `12` to:

- `RunCli` in `src/tools/trackclustertu.rs`
- `ClusterCli` in `src/tools/cluster_pipeline.rs`

When building `TuClusteringOptions`, pass:

- `attach_contained_reads: !cli.skip_score2_attachment`
- `three_prime_tolerance_bp: cli.three_prime_tolerance_bp`

The `run` command already forwards cluster-related flags to `cluster`; it must also forward `--three-prime-tolerance-bp`.

## Testing

Add or update unit tests in `src/tu/mod.rs` to cover:

1. `+` strand, longer candidate starts earlier and ends up to `12 bp` before the child's 3 prime end.
   Expected result: attach into one TU when `score2` passes.
2. `-` strand mirror of the same case.
   Expected result: attach into one TU when the candidate starts up to `12 bp` downstream of the child's 3 prime end.
3. A case where the 3 prime gap exceeds `12 bp`.
   Expected result: remain separate TUs even when `score2` passes.
4. Existing `--skip-score2-attachment` behavior remains unchanged.

CLI tests should also cover the new flag surface where appropriate:

- `trackclustertu cluster --help` mentions `--three-prime-tolerance-bp`
- `trackclustertu run` forwards the option into the cluster subcommand arguments

## Documentation

Update:

- `README.md`
- any command or output docs that describe second-pass score2 attachment

The docs should explain that second-pass score2 attachment is strand-aware and allows a configurable 3 prime mismatch, default `12 bp`.

## Risks And Mitigations

### Risk: Over-merging distinct TU ends

Mitigation:

- keep the tolerance small and bounded
- require candidate length to still be strictly greater
- keep the existing `score2` threshold gate
- provide a CLI override for stricter or looser behavior

### Risk: Incorrect minus-strand interpretation

Mitigation:

- implement the tolerance check in a dedicated helper
- add an explicit mirrored minus-strand regression test

### Risk: Hidden behavior change in default clustering

Mitigation:

- expose the value on the CLI
- document the default
- keep `0` as an implicit compatibility mode equivalent to the old hard rule

## Acceptance Criteria

- The example class of `+`-strand near-matching TUs with `score2` above threshold and a 3 prime gap of at most `12 bp` merges in the second pass.
- The mirrored `-`-strand case also merges.
- Cases exceeding the configured tolerance remain split.
- `trackclustertu cluster` and `trackclustertu run` both accept `--three-prime-tolerance-bp`.
- The default CLI behavior uses `12 bp`.
- Existing behavior with `--skip-score2-attachment` is unchanged.
