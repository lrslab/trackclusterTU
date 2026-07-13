# Minimal TU example (BED6)

Inputs:

- `reads.bed` (BED6)

Expected outputs are under `expected/`. They reflect the current defaults: stable
coordinate-derived TU IDs, canonical v2 membership rows, and explicit count
semantics in `tu_count.csv`.

Run:

```bash
trackclustertu cluster \
  --in reads.bed \
  --format bed6 \
  --span-jaccard-threshold 0.95 \
  --overlap-over-longer-threshold 0.80 \
  --out-tu tus.bed \
  --out-membership membership.tsv \
  --out-tu-count tu_count.csv
```
Compare the generated files with `expected/`.
