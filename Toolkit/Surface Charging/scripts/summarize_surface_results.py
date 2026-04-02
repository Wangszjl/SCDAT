#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import pathlib


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize a SCDAT surface charging CSV export.")
    parser.add_argument("csv_path")
    args = parser.parse_args()

    path = pathlib.Path(args.csv_path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    print(f"path={path}")
    print(f"rows={len(rows)}")
    if not rows:
        return 0

    last = rows[-1]
    for key in [
        "surface_potential_v",
        "body_potential_v",
        "patch_potential_v",
        "total_current_density_a_per_m2",
        "floating_equilibrium_potential_v",
        "equilibrium_error_v",
    ]:
        value = last.get(key, "")
        if value != "":
            print(f"last_{key}={value}")

    node_columns = [key for key in reader.fieldnames or [] if key.startswith("surface_node_") and key.endswith("_potential_v")]
    interface_columns = [key for key in reader.fieldnames or [] if key.startswith("surface_interface_") and key.endswith("_current_a")]
    print(f"surface_node_potential_columns={len(node_columns)}")
    print(f"surface_interface_current_columns={len(interface_columns)}")

    finite_equilibrium_errors = []
    for row in rows:
        try:
            finite_equilibrium_errors.append(abs(float(row["equilibrium_error_v"])))
        except Exception:
            continue
    if finite_equilibrium_errors:
        print(f"max_abs_equilibrium_error_v={max(finite_equilibrium_errors):.12g}")
        print(f"min_abs_equilibrium_error_v={min(finite_equilibrium_errors):.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
