# Minimal TU example (BED6)

Inputs:

- `reads.bed` (BED6)

Selected expected outputs are under `expected/`. They cover the three core
files: stable coordinate-derived TU IDs, canonical v2 membership rows, and
explicit count semantics in `tu_count.csv`.

The fixtures match version 0.2.2. Contained reads `r3` and `r6` also record a
covering TU as a runner-up candidate; their stronger direct matches remain
unique assignments.

Run:

```bash
trackclustertu cluster \
  --in reads.bed \
  --format bed6 \
  --span-jaccard-threshold 0.95 \
  --overlap-over-longer-threshold 0.80 \
  --out-dir actual
```

Compare `actual/tus.bed`, `actual/membership.tsv`, and `actual/tu_count.csv`
with their counterparts in `expected/`. The command also writes the current
audit/default sidecars in `actual/`; those are intentionally outside this small
fixture set.
