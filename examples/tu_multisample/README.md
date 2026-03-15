# Multi-sample TU example

Inputs:

- `samples.tsv` manifest with `sample`, `reads`, and `group`
- `sampleA.bed`
- `sampleB.bed`

Expected outputs are under `expected/`.

Run:

```bash
trackclustertu cluster \
  --manifest samples.tsv \
  --format bed6 \
  --out-tu tus.bed \
  --out-membership membership.tsv \
  --out-pooled-reads pooled.bed \
  --out-tu-count tu_count.csv \
  --out-tu-sample-count-long sample_long.tsv \
  --out-tu-sample-count-matrix sample_matrix.tsv \
  --out-tu-group-count-matrix group_matrix.tsv
```
Compare the generated files with `expected/`.
