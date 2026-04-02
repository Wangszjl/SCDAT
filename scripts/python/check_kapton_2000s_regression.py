#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise ValueError(f"empty csv: {csv_path}")
    return rows


def parse_required_float(row: dict[str, str], column: str, csv_path: Path) -> float:
    raw = row.get(column, "")
    try:
        value = float(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid float in {csv_path.name}:{column} -> '{raw}'") from exc
    if not math.isfinite(value):
        raise ValueError(f"non-finite value in {csv_path.name}:{column}")
    return value


def find_nearest_row(rows: list[dict[str, str]], csv_path: Path, target_time_s: float) -> tuple[dict[str, str], float, float]:
    best_row = rows[0]
    best_time = parse_required_float(best_row, "time_s", csv_path)
    best_offset = abs(best_time - target_time_s)

    for row in rows[1:]:
        time_value = parse_required_float(row, "time_s", csv_path)
        offset = abs(time_value - target_time_s)
        if offset < best_offset:
            best_row = row
            best_time = time_value
            best_offset = offset

    return best_row, best_time, best_offset


def write_summary_csv(
    output_path: Path,
    target_time_s: float,
    advanced_time_s: float,
    reference_time_s: float,
    advanced_surface_v: float,
    reference_surface_v: float,
    delta_surface_v: float,
    equilibrium_gap_v: float | None,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerow(["target_time_s", f"{target_time_s:.12g}"])
        writer.writerow(["advanced_time_s", f"{advanced_time_s:.12g}"])
        writer.writerow(["reference_time_s", f"{reference_time_s:.12g}"])
        writer.writerow(["advanced_surface_potential_v", f"{advanced_surface_v:.12g}"])
        writer.writerow(["reference_surface_potential_v", f"{reference_surface_v:.12g}"])
        writer.writerow(["delta_surface_potential_v", f"{delta_surface_v:.12g}"])
        if equilibrium_gap_v is not None:
            writer.writerow(["equilibrium_gap_v", f"{equilibrium_gap_v:.12g}"])


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check Kapton 2000s regression between advanced and 1D reference outputs."
    )
    parser.add_argument(
        "--advanced-csv",
        default="results/surface_geo_ecss_kapton_pic_circuit.csv",
        help="Advanced model CSV path.",
    )
    parser.add_argument(
        "--reference-csv",
        default="results/surface_geo_ecss_kapton_ref.csv",
        help="1D reference CSV path.",
    )
    parser.add_argument(
        "--target-time-s",
        type=float,
        default=2000.0,
        help="Target time in seconds.",
    )
    parser.add_argument(
        "--max-time-offset-s",
        type=float,
        default=1.0,
        help="Maximum allowed nearest-sample time offset.",
    )
    parser.add_argument(
        "--max-abs-surface-delta-v",
        type=float,
        default=3000.0,
        help="Maximum allowed absolute delta in surface potential at target time.",
    )
    parser.add_argument(
        "--max-equilibrium-gap-v",
        type=float,
        default=1.0e-6,
        help="Maximum allowed equilibrium potential gap when both CSVs provide the field.",
    )
    parser.add_argument(
        "--require-advanced-more-negative",
        action="store_true",
        help="Require advanced surface potential to be more negative than reference.",
    )
    parser.add_argument(
        "--summary-csv",
        default="",
        help="Optional output CSV path for summary metrics.",
    )
    args = parser.parse_args()

    advanced_path = Path(args.advanced_csv)
    reference_path = Path(args.reference_csv)

    if not advanced_path.is_file():
        print(f"status=FAIL(missing_advanced_csv:{advanced_path})")
        return 3
    if not reference_path.is_file():
        print(f"status=FAIL(missing_reference_csv:{reference_path})")
        return 3

    advanced_rows = load_rows(advanced_path)
    reference_rows = load_rows(reference_path)

    advanced_row, advanced_time_s, advanced_offset = find_nearest_row(
        advanced_rows, advanced_path, args.target_time_s
    )
    reference_row, reference_time_s, reference_offset = find_nearest_row(
        reference_rows, reference_path, args.target_time_s
    )

    advanced_surface_v = parse_required_float(advanced_row, "surface_potential_v", advanced_path)
    reference_surface_v = parse_required_float(reference_row, "surface_potential_v", reference_path)
    delta_surface_v = advanced_surface_v - reference_surface_v
    abs_delta_surface_v = abs(delta_surface_v)

    equilibrium_gap_v: float | None = None
    if (
        "floating_equilibrium_potential_v" in advanced_row
        and "floating_equilibrium_potential_v" in reference_row
    ):
        advanced_eq_v = parse_required_float(
            advanced_row, "floating_equilibrium_potential_v", advanced_path
        )
        reference_eq_v = parse_required_float(
            reference_row, "floating_equilibrium_potential_v", reference_path
        )
        equilibrium_gap_v = abs(advanced_eq_v - reference_eq_v)

    print(f"target_time_s={args.target_time_s:.12g}")
    print(f"advanced_time_s={advanced_time_s:.12g} offset_s={advanced_offset:.12g}")
    print(f"reference_time_s={reference_time_s:.12g} offset_s={reference_offset:.12g}")
    print(f"advanced_surface_potential_v={advanced_surface_v:.12g}")
    print(f"reference_surface_potential_v={reference_surface_v:.12g}")
    print(f"delta_surface_potential_v={delta_surface_v:.12g}")
    print(f"abs_delta_surface_potential_v={abs_delta_surface_v:.12g}")

    failures: list[str] = []
    if advanced_offset > args.max_time_offset_s:
        failures.append(
            f"advanced_time_offset({advanced_offset:.6g})>max({args.max_time_offset_s:.6g})"
        )
    if reference_offset > args.max_time_offset_s:
        failures.append(
            f"reference_time_offset({reference_offset:.6g})>max({args.max_time_offset_s:.6g})"
        )
    if abs_delta_surface_v > args.max_abs_surface_delta_v:
        failures.append(
            f"abs_surface_delta_v({abs_delta_surface_v:.6g})>max({args.max_abs_surface_delta_v:.6g})"
        )
    if args.require_advanced_more_negative and not (advanced_surface_v < reference_surface_v):
        failures.append("advanced_surface_potential_is_not_more_negative")

    if equilibrium_gap_v is not None:
        print(f"equilibrium_gap_v={equilibrium_gap_v:.12g}")
        if equilibrium_gap_v > args.max_equilibrium_gap_v:
            failures.append(
                f"equilibrium_gap_v({equilibrium_gap_v:.6g})>max({args.max_equilibrium_gap_v:.6g})"
            )

    if args.summary_csv:
        write_summary_csv(
            Path(args.summary_csv),
            args.target_time_s,
            advanced_time_s,
            reference_time_s,
            advanced_surface_v,
            reference_surface_v,
            delta_surface_v,
            equilibrium_gap_v,
        )
        print(f"summary_csv={Path(args.summary_csv).resolve()}")

    if failures:
        for failure in failures:
            print(f"failure={failure}")
        print("status=FAIL")
        return 2

    print("status=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
