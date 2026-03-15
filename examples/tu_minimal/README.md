# Minimal TU example (BED6)

Inputs:

- `reads.bed` (BED6)

Expected outputs are under `expected/`.

Run:

```bash
trackclustertu cluster \
  --in reads.bed \
  --format bed6 \
  --score1-threshold 0.95 \
  --score2-threshold 0.99 \
  --out-tu tus.bed \
  --out-membership membership.tsv \
  --out-tu-count tu_count.csv
```
Compare the generated files with `expected/`.
