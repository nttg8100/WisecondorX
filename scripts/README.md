# Cross-call audit

`cross_call_audit.py` is a benchmark utility for `wisecondorx predict` output.
It reads per-sample `*_statistics.txt` files and counts primary and secondary
target-chromosome calls at a chosen z threshold.

The utility exists because a naive pandas row sum over a mixed object/boolean
DataFrame can perform a logical OR instead of counting calls.  A sample with two
calls is then counted as one, and cross-call counts incorrectly appear as zero.
This script uses explicit per-chromosome booleans.

Usage:

```bash
python scripts/cross_call_audit.py \
  --predict-dir path/to/predict \
  --labels labels.tsv \
  --outdir cross_call_audit \
  --z-threshold 3
```

`labels.{tsv,csv}` needs `sample` and `expected` columns.  `expected` may be
`control`, `T13`, `T18`, `T21`, or a `+`-joined pattern for research mixtures.

Outputs:

- `cross_call_audit.tsv`: per-sample z-scores, primary chromosome, calls, and
  whether more than one target chromosome was called.
- `cross_call_metrics.json`: counts of samples with any call and with multiple
  calls.

This is a benchmark utility only; it does not change the WisecondorX algorithm
or its clinical interpretation.
