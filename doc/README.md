# Documentation

These docs describe `trackclustertu` **0.2.2**, the version in
[`Cargo.toml`](../Cargo.toml). See the [changelog](../CHANGELOG.md) for changes
between versions; check an installed binary with `trackclustertu --version`.

- [Workflow diagram](pipeline.svg)
- [Overview, installation, and usage](../README.md)
- [Output files, schemas, and count semantics](output_directory.md)
- [Runnable examples with expected outputs](../examples/README.md)
- [Performance notes and historical benchmarks](performance.md)
- [Release procedure and artifact verification](release.md)

## Commands

Use `trackclustertu <command> --help` for the options and defaults of each command.

| Command | Input | Result |
| --- | --- | --- |
| `run` | FASTQ sample manifest + reference FASTA | Mapping, BED conversion, shared TU discovery, assignments, and counts; gene outputs with optional annotation |
| `map` | FASTQ sample manifest + reference FASTA | Sorted/indexed BAMs, BED6 reads, BAM/BED manifests, logs, and run provenance |
| `bam-to-bed` | One BAM file or a BAM sample manifest | BED6 reads; optional evidence sidecars |
| `cluster` | BED6, BED12, or TSV reads, as one file or a sample manifest | TU catalogue, read membership, counts, and audit tables; gene outputs with optional BED annotation |
| `recount` | Sample manifest + existing membership TSV | TU totals and sample/group counts |
| `diagnose-missed-tus` | Clustered BED6 reads + existing TU BED | Report of supported boundary modes missing from the TU catalogue |
| `rescue-missed-tus` | Clustered BED6 reads + existing TU BED + membership TSV | Existing and rescued TUs, updated membership, and optional hard counts/candidate reports |
| `gff-to-bed` | GFF3 gene annotation | BED6 gene annotation |

`run` covers steps 1–3 of the diagram. `map` already performs BAM-to-BED
conversion. Diagnosis and rescue are independent post-clustering commands;
rescue uses the original inputs, not the diagnostic report. Recounting rescued
membership updates TU counts; it does not regenerate gene outputs.

TU discovery uses genomic spans and strand without requiring gene annotation.
In 0.2.2, contained fragments can add support without moving consensus
boundaries. Disabling the attachment pass also disables fragment absorption
and partial assignment, while final direct-consensus refinement still applies.
See [similarity scores and CLI aliases](../README.md#similarity-scores).
