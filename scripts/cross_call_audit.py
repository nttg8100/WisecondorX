#!/usr/bin/env python3
# SPDX-License-Identifier: CC-BY-NC-SA-4.0
"""Audit chromosome-level calls and cross-calls from WisecondorX output.

This is a benchmark utility, not part of the WisecondorX algorithm.  It reads
the per-sample ``*_statistics.txt`` files produced by ``wisecondorx predict``
and counts primary and secondary target-chromosome calls.

The important detail is that cross-calls are counted with explicit per-
chromosome booleans.  A mixed object/boolean DataFrame row must not be summed
directly: pandas then performs a logical OR and a two-call sample is reported
as one call, making cross-call counts appear to be zero.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

TARGETS = ("13", "18", "21")


def read_statistics(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path, sep="\t")
    frame["chr"] = frame["chr"].astype(str).str.replace("chr", "", regex=False)
    values = frame.set_index("chr")["zscore"]
    return {chrom: float(values.get(chrom, np.nan)) for chrom in TARGETS}


def read_labels(path: Path) -> dict[str, str]:
    first = path.read_text(encoding="utf-8").splitlines()[0]
    separator = "\t" if "\t" in first else ","
    frame = pd.read_csv(path, sep=separator, dtype=str)
    missing = {"sample", "expected"}.difference(frame.columns)
    if missing:
        raise SystemExit(f"labels table missing columns: {sorted(missing)}")
    return dict(zip(frame["sample"], frame["expected"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predict-dir", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--z-threshold", type=float, default=3.0)
    parser.add_argument(
        "--secondary-threshold",
        type=float,
        default=None,
        help="Optional separate threshold for secondary calls; defaults to z-threshold.",
    )
    args = parser.parse_args()

    labels = read_labels(args.labels)
    secondary_threshold = (
        args.z_threshold if args.secondary_threshold is None else args.secondary_threshold
    )
    rows = []
    for sample, expected in labels.items():
        path = args.predict_dir / f"{sample}_statistics.txt"
        if not path.is_file():
            raise SystemExit(f"missing {path}")
        z = read_statistics(path)
        calls = {chrom: bool(np.isfinite(z[chrom]) and z[chrom] >= args.z_threshold) for chrom in TARGETS}
        secondary = {
            chrom: bool(np.isfinite(z[chrom]) and z[chrom] >= secondary_threshold)
            for chrom in TARGETS
        }
        primary = None
        finite = {chrom: value for chrom, value in z.items() if np.isfinite(value)}
        if finite:
            primary = max(finite, key=finite.get)
        n_calls = sum(calls.values())
        rows.append(
            {
                "sample": sample,
                "expected": expected,
                **{f"z{chrom}": z[chrom] for chrom in TARGETS},
                "primary": primary or "none",
                "calls": "+".join(chrom for chrom in TARGETS if calls[chrom]) or "none",
                "secondary_calls": "+".join(
                    chrom for chrom in TARGETS if secondary[chrom]
                )
                or "none",
                "cross_call": n_calls > 1,
                "n_calls": n_calls,
            }
        )

    frame = pd.DataFrame(rows)
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(outdir / "cross_call_audit.tsv", sep="\t", index=False)

    expected_calls = frame["calls"].str.contains("+", regex=False)
    metrics = {
        "n_samples": int(len(frame)),
        "threshold": args.z_threshold,
        "secondary_threshold": secondary_threshold,
        "samples_any_call": int((frame["n_calls"] > 0).sum()),
        "samples_multiple_calls": int((frame["n_calls"] > 1).sum()),
        "primary_matches_expected": int(
            (
                frame.apply(
                    lambda row: row["expected"] != "control"
                    and row["primary"] == str(row["expected"]).removeprefix("T"),
                    axis=1,
                )
            ).sum()
        ),
        "expected_chromosome_called": int(
            (
                frame.apply(
                    lambda row: row["expected"] != "control"
                    and str(row["expected"]).removeprefix("T") in str(row["calls"]).split("+"),
                    axis=1,
                )
            ).sum()
        ),
        "note": (
            "Cross-call means more than one target chromosome passed the "
            "threshold. It is not a diagnosis; WisecondorX also emits secondary "
            "signals in the raw statistics."
        ),
    }
    (outdir / "cross_call_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
