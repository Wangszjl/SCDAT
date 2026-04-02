#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import pathlib
from typing import Dict, List


def read_csv(path: pathlib.Path) -> List[Dict[str, float]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows: List[Dict[str, float]] = []
        for row in reader:
            parsed: Dict[str, float] = {}
            for key, value in row.items():
                if value is None or value == "":
                    continue
                try:
                    parsed[key] = float(value)
                except ValueError:
                    continue
            rows.append(parsed)
        return rows


def rmse(lhs: List[float], rhs: List[float]) -> float:
    count = min(len(lhs), len(rhs))
    if count == 0:
        return math.nan
    return math.sqrt(sum((lhs[i] - rhs[i]) ** 2 for i in range(count)) / count)


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare two SCDAT benchmark/result CSV files.")
    parser.add_argument("actual")
    parser.add_argument("reference")
    parser.add_argument("--columns", nargs="+", default=["surface_potential_v", "total_current_density_a_per_m2"])
    args = parser.parse_args()

    actual_rows = read_csv(pathlib.Path(args.actual))
    reference_rows = read_csv(pathlib.Path(args.reference))
    print(f"actual_rows={len(actual_rows)}")
    print(f"reference_rows={len(reference_rows)}")

    for column in args.columns:
        actual_values = [row[column] for row in actual_rows if column in row]
        reference_values = [row[column] for row in reference_rows if column in row]
        value_rmse = rmse(actual_values, reference_values)
        actual_end = actual_values[min(len(actual_values), len(reference_values)) - 1] if actual_values and reference_values else math.nan
        reference_end = reference_values[min(len(actual_values), len(reference_values)) - 1] if actual_values and reference_values else math.nan
        end_delta = actual_end - reference_end if not math.isnan(actual_end) and not math.isnan(reference_end) else math.nan
        print(f"{column}: rmse={value_rmse:.12g} end_delta={end_delta:.12g}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
