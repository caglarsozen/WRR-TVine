#!/usr/bin/env python3
"""Compare regenerated aggregate diagnostics with the reference outputs.

Usage:
    python analysis/compare_reference.py analysis_output/derived

This compares numerical summaries, not fitted model objects or raw data.
"""
from pathlib import Path
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("derived", type=Path, help="Output from the audit and supplementary diagnostics")
parser.add_argument("--reference", type=Path, default=Path(__file__).resolve().parents[1] / "reference_results")
parser.add_argument("--rtol", type=float, default=1e-8)
parser.add_argument("--atol", type=float, default=1e-10)
args = parser.parse_args()
files = sorted(args.reference.glob("*.csv"))
if not files:
    parser.error("The reference directory contains no CSV files.")

failed = []
for expected in files:
    actual = args.derived / expected.name
    try:
        if not actual.is_file():
            raise FileNotFoundError(f"Missing regenerated file: {actual}")
        a = pd.read_csv(expected)
        b = pd.read_csv(actual)
        if expected.name == "pit_calibration_checked.csv":
            # Compare the documented tests across both reporting conventions.
            raw_p_columns = ["calibration_hac_p", "calibration_holm_p"]
            a = a.drop(columns=raw_p_columns)
            b = b.drop(columns=raw_p_columns)
        pd.testing.assert_frame_equal(
            a, b, check_dtype=False, rtol=args.rtol, atol=args.atol
        )
        print(f"PASS {expected.name}")
    except (AssertionError, FileNotFoundError, ValueError) as exc:
        failed.append(expected.name)
        print(f"FAIL {expected.name}: {exc}")
print(f"{len(files) - len(failed)} of {len(files)} reference tables match.")
raise SystemExit(1 if failed else 0)
