# Multi-sample TU example

Inputs:

- `samples.tsv` manifest with `sample`, `reads`, and `group`
- `sampleA.bed`
- `sampleB.bed`

Selected expected outputs are under `expected/`. They cover stable
coordinate-derived TU IDs, canonical v2 membership rows, explicit count
semantics in `tu_count.csv`, and count-semantics columns in the sample/group
matrices.

The fixtures match version 0.2.2. The contained read `sampleB::r1` also records
the covering TU as a runner-up candidate while retaining its unique direct
assignment.

Run:

```bash
trackclustertu cluster \
  --manifest samples.tsv \
  --format bed6 \
  --out-dir actual \
  --out-tu-sample-count-long actual/sample_long.tsv \
  --out-tu-sample-count-matrix actual/sample_matrix.tsv \
  --out-tu-group-count-matrix actual/group_matrix.tsv
```

Compare the same-named files in `actual/` and `expected/`. The command also
writes the current audit/default sidecars in `actual/`; those are intentionally
outside this compact fixture set.
