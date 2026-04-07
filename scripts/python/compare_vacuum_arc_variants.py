#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

EPSILON = 1.0e-12


def parse_float(raw: str, csv_path: Path, row_index: int, column: str) -> float:
    try:
        value = float(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"invalid float in {csv_path.name}:{row_index}:{column} -> '{raw}'"
        ) from exc
    if not math.isfinite(value):
        raise ValueError(
            f"non-finite float in {csv_path.name}:{row_index}:{column} -> '{raw}'"
        )
    return value


def load_rows(csv_path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    if not fieldnames:
        raise ValueError(f"missing CSV header: {csv_path}")
    if not rows:
        raise ValueError(f"empty CSV: {csv_path}")
    if "time_s" not in fieldnames:
        raise ValueError(f"missing required column 'time_s' in {csv_path}")
    return fieldnames, rows


def parse_dataset(
    rows: List[Dict[str, str]],
    csv_path: Path,
    columns: List[str],
) -> Tuple[List[float], Dict[str, List[float]]]:
    records: List[Tuple[float, Dict[str, float]]] = []
    for row_index, row in enumerate(rows, start=2):
        time_value = parse_float(row.get("time_s", ""), csv_path, row_index, "time_s")
        values: Dict[str, float] = {}
        for column in columns:
            values[column] = parse_float(row.get(column, ""), csv_path, row_index, column)
        records.append((time_value, values))

    records.sort(key=lambda item: item[0])

    deduped: List[Tuple[float, Dict[str, float]]] = []
    for time_value, values in records:
        if deduped and abs(time_value - deduped[-1][0]) <= EPSILON:
            deduped[-1] = (time_value, values)
            continue
        deduped.append((time_value, values))

    times = [item[0] for item in deduped]
    series_by_column: Dict[str, List[float]] = {column: [] for column in columns}
    for _, values in deduped:
        for column in columns:
            series_by_column[column].append(values[column])

    return times, series_by_column


def interpolate_series(
    target_times: Iterable[float],
    source_times: List[float],
    source_values: List[float],
) -> List[float]:
    if not source_times:
        raise ValueError("source_times must not be empty")
    if len(source_times) != len(source_values):
        raise ValueError("source_times/source_values length mismatch")
    if len(source_times) == 1:
        constant = source_values[0]
        return [constant for _ in target_times]

    output: List[float] = []
    left_index = 0
    for target_time in target_times:
        if target_time <= source_times[0]:
            output.append(source_values[0])
            continue
        if target_time >= source_times[-1]:
            output.append(source_values[-1])
            continue

        while (
            left_index + 1 < len(source_times)
            and source_times[left_index + 1] < target_time
        ):
            left_index += 1

        t0 = source_times[left_index]
        t1 = source_times[left_index + 1]
        v0 = source_values[left_index]
        v1 = source_values[left_index + 1]

        if abs(target_time - t0) <= EPSILON:
            output.append(v0)
            continue
        if abs(target_time - t1) <= EPSILON:
            output.append(v1)
            continue

        if abs(t1 - t0) <= EPSILON:
            output.append(v1)
            continue

        weight = (target_time - t0) / (t1 - t0)
        output.append(v0 + weight * (v1 - v0))

    return output


def rms(values: List[float]) -> float:
    if not values:
        return 0.0
    return math.sqrt(sum(value * value for value in values) / len(values))


def trapz_integral(times: List[float], values: List[float]) -> float:
    if len(times) < 2 or len(values) < 2:
        return 0.0
    total = 0.0
    for index in range(1, len(times)):
        dt = times[index] - times[index - 1]
        if dt < 0.0:
            raise ValueError("time axis must be non-decreasing")
        total += 0.5 * (values[index] + values[index - 1]) * dt
    return total


def first_crossing_time(
    times: List[float],
    values: List[float],
    threshold: float,
) -> float | None:
    if len(times) != len(values):
        raise ValueError("times/values length mismatch")
    if not times:
        return None

    previous_time = times[0]
    previous_value = abs(values[0])
    if previous_value >= threshold:
        return previous_time

    for index in range(1, len(times)):
        current_time = times[index]
        current_value = abs(values[index])
        if current_value >= threshold:
            if abs(current_value - previous_value) <= EPSILON:
                return current_time
            weight = (threshold - previous_value) / (current_value - previous_value)
            weight = min(1.0, max(0.0, weight))
            return previous_time + weight * (current_time - previous_time)
        previous_time = current_time
        previous_value = current_value

    return None


def relative_delta(candidate: float, reference: float) -> float:
    return abs(candidate - reference) / max(abs(reference), EPSILON)


def finite_or_inf(value: float) -> float:
    return value if math.isfinite(value) else float("inf")


def parse_tail_fraction(raw: float) -> float:
    if raw <= 0.0:
        return EPSILON
    if raw > 1.0:
        return 1.0
    return raw


def write_summary_markdown(
    output_path: Path,
    rows: List[Dict[str, str]],
    breakdown_summary: Dict[str, str] | None,
    status: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    lines.append("# Vacuum Arc Comparison Summary")
    lines.append("")
    lines.append(f"- overall_status: {status}")
    lines.append("")
    lines.append(
        "| column | samples | max_abs_delta | rel_rmse | rel_peak_delta | rel_final_delta | rel_integral_delta | rel_steady_mean_delta | status |"
    )
    lines.append(
        "|---|---:|---:|---:|---:|---:|---:|---:|---|"
    )
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["column"],
                    row["samples"],
                    row["max_abs_delta"],
                    row["rel_rmse"],
                    row["rel_peak_delta"],
                    row["rel_final_delta"],
                    row["rel_integral_delta"],
                    row["rel_steady_mean_delta"],
                    row["status"],
                ]
            )
            + " |"
        )

    if breakdown_summary is not None:
        lines.append("")
        lines.append("## Breakdown Timing")
        lines.append("")
        for key in [
            "column",
            "threshold",
            "reference_breakdown_time_s",
            "candidate_breakdown_time_s",
            "breakdown_time_delta_s",
            "status",
        ]:
            if key in breakdown_summary:
                lines.append(f"- {key}: {breakdown_summary[key]}")

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_summary_csv(output_path: Path, rows: List[Dict[str, str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "column",
                "samples",
                "max_abs_delta",
                "mean_abs_delta",
                "rmse",
                "rel_rmse",
                "rel_peak_delta",
                "rel_final_delta",
                "rel_integral_delta",
                "rel_steady_mean_delta",
                "status",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare two Vacuum Arc CSV outputs with time alignment and per-metric deltas."
    )
    parser.add_argument(
        "--reference-csv",
        required=True,
        help="Reference CSV path (e.g., ArcPIC baseline).",
    )
    parser.add_argument(
        "--candidate-csv",
        required=True,
        help="Candidate CSV path produced by current implementation.",
    )
    parser.add_argument(
        "--columns",
        nargs="+",
        default=[
            "discharge_current_a",
            "current_density_a_per_m2",
            "surface_potential_v",
            "cathode_temperature_k",
            "pic_boundary_net_current_a",
        ],
        help="Columns to compare.",
    )
    parser.add_argument(
        "--allow-missing-columns",
        action="store_true",
        help="When set, only compare columns that exist in both CSVs.",
    )
    parser.add_argument(
        "--max-rel-rmse",
        type=float,
        default=float("inf"),
        help="Fail when rel_rmse exceeds this value.",
    )
    parser.add_argument(
        "--max-rel-peak-delta",
        type=float,
        default=float("inf"),
        help="Fail when relative peak delta exceeds this value.",
    )
    parser.add_argument(
        "--max-rel-final-delta",
        type=float,
        default=float("inf"),
        help="Fail when relative final-value delta exceeds this value.",
    )
    parser.add_argument(
        "--max-rel-integral-delta",
        type=float,
        default=float("inf"),
        help="Fail when relative integral delta exceeds this value.",
    )
    parser.add_argument(
        "--max-rel-steady-mean-delta",
        type=float,
        default=float("inf"),
        help="Fail when relative steady-segment mean delta exceeds this value.",
    )
    parser.add_argument(
        "--max-abs-delta",
        type=float,
        default=float("inf"),
        help="Fail when absolute max delta exceeds this value.",
    )
    parser.add_argument(
        "--breakdown-column",
        default="discharge_current_a",
        help="Column used to estimate first breakdown time using threshold crossing.",
    )
    parser.add_argument(
        "--breakdown-threshold",
        type=float,
        default=float("nan"),
        help="Absolute threshold for first-breakdown detection; unset to disable.",
    )
    parser.add_argument(
        "--max-breakdown-time-delta-s",
        type=float,
        default=float("inf"),
        help="Fail when absolute breakdown timing delta exceeds this value.",
    )
    parser.add_argument(
        "--fail-on-missing-breakdown",
        action="store_true",
        help="Fail if breakdown threshold is set but either side never crosses it.",
    )
    parser.add_argument(
        "--steady-tail-fraction",
        type=float,
        default=0.2,
        help="Tail window fraction used for steady-segment mean calculation.",
    )
    parser.add_argument(
        "--summary-json",
        default="",
        help="Optional output JSON path for detailed summary.",
    )
    parser.add_argument(
        "--summary-csv",
        default="",
        help="Optional output CSV path for per-column summary.",
    )
    parser.add_argument(
        "--summary-md",
        default="",
        help="Optional output Markdown path for human-readable summary.",
    )
    args = parser.parse_args()

    tail_fraction = parse_tail_fraction(args.steady_tail_fraction)

    reference_path = Path(args.reference_csv)
    candidate_path = Path(args.candidate_csv)

    if not reference_path.is_file():
        print(f"status=FAIL(missing_reference_csv:{reference_path})")
        return 3
    if not candidate_path.is_file():
        print(f"status=FAIL(missing_candidate_csv:{candidate_path})")
        return 3

    try:
        reference_fields, reference_rows = load_rows(reference_path)
        candidate_fields, candidate_rows = load_rows(candidate_path)
    except ValueError as exc:
        print(f"status=FAIL({exc})")
        return 3

    requested_columns = list(dict.fromkeys(args.columns))
    reference_missing = [column for column in requested_columns if column not in reference_fields]
    candidate_missing = [column for column in requested_columns if column not in candidate_fields]

    if args.allow_missing_columns:
        compared_columns = [
            column
            for column in requested_columns
            if column in reference_fields and column in candidate_fields
        ]
    else:
        compared_columns = requested_columns

    if (reference_missing or candidate_missing) and not args.allow_missing_columns:
        if reference_missing:
            print(f"missing_in_reference={','.join(reference_missing)}")
        if candidate_missing:
            print(f"missing_in_candidate={','.join(candidate_missing)}")
        print("status=FAIL(missing_columns)")
        return 2

    if not compared_columns:
        print("status=FAIL(no_common_columns_to_compare)")
        return 2

    if args.allow_missing_columns:
        skipped = sorted(set(reference_missing + candidate_missing))
        if skipped:
            print(f"skipped_columns={','.join(skipped)}")

    try:
        reference_times, reference_series = parse_dataset(
            reference_rows, reference_path, compared_columns
        )
        candidate_times, candidate_series = parse_dataset(
            candidate_rows, candidate_path, compared_columns
        )
    except ValueError as exc:
        print(f"status=FAIL({exc})")
        return 3

    summary_rows: List[Dict[str, str]] = []
    failures: List[str] = []
    breakdown_summary: Dict[str, str] | None = None
    aligned_candidate_series: Dict[str, List[float]] = {}

    for column in compared_columns:
        reference_values = reference_series[column]
        candidate_values = interpolate_series(
            reference_times, candidate_times, candidate_series[column]
        )
        aligned_candidate_series[column] = candidate_values

        deltas = [candidate - reference for candidate, reference in zip(candidate_values, reference_values)]
        abs_deltas = [abs(value) for value in deltas]

        max_abs_delta = max(abs_deltas) if abs_deltas else 0.0
        mean_abs_delta = sum(abs_deltas) / len(abs_deltas) if abs_deltas else 0.0
        rmse_value = rms(deltas)

        reference_rms = rms(reference_values)
        rel_rmse = rmse_value / max(reference_rms, EPSILON)

        reference_peak = max(abs(value) for value in reference_values)
        candidate_peak = max(abs(value) for value in candidate_values)
        rel_peak_delta = abs(candidate_peak - reference_peak) / max(reference_peak, EPSILON)

        reference_final = reference_values[-1]
        candidate_final = candidate_values[-1]
        rel_final_delta = abs(candidate_final - reference_final) / max(abs(reference_final), EPSILON)

        reference_integral = trapz_integral(reference_times, reference_values)
        candidate_integral = trapz_integral(reference_times, candidate_values)
        rel_integral_delta = relative_delta(candidate_integral, reference_integral)

        tail_count = max(1, int(math.ceil(len(reference_values) * tail_fraction)))
        steady_reference_values = reference_values[-tail_count:]
        steady_candidate_values = candidate_values[-tail_count:]
        steady_reference_mean = statistics.fmean(steady_reference_values)
        steady_candidate_mean = statistics.fmean(steady_candidate_values)
        rel_steady_mean_delta = relative_delta(steady_candidate_mean, steady_reference_mean)

        status = "PASS"
        if max_abs_delta > args.max_abs_delta:
            status = "FAIL"
            failures.append(
                f"{column}:max_abs_delta({max_abs_delta:.6g})>{args.max_abs_delta:.6g}"
            )
        if rel_rmse > args.max_rel_rmse:
            status = "FAIL"
            failures.append(f"{column}:rel_rmse({rel_rmse:.6g})>{args.max_rel_rmse:.6g}")
        if rel_peak_delta > args.max_rel_peak_delta:
            status = "FAIL"
            failures.append(
                f"{column}:rel_peak_delta({rel_peak_delta:.6g})>{args.max_rel_peak_delta:.6g}"
            )
        if rel_final_delta > args.max_rel_final_delta:
            status = "FAIL"
            failures.append(
                f"{column}:rel_final_delta({rel_final_delta:.6g})>{args.max_rel_final_delta:.6g}"
            )
        if rel_integral_delta > args.max_rel_integral_delta:
            status = "FAIL"
            failures.append(
                f"{column}:rel_integral_delta({rel_integral_delta:.6g})>{args.max_rel_integral_delta:.6g}"
            )
        if rel_steady_mean_delta > args.max_rel_steady_mean_delta:
            status = "FAIL"
            failures.append(
                f"{column}:rel_steady_mean_delta({rel_steady_mean_delta:.6g})>{args.max_rel_steady_mean_delta:.6g}"
            )

        print(
            " ".join(
                [
                    f"column={column}",
                    f"samples={len(reference_values)}",
                    f"max_abs_delta={max_abs_delta:.12g}",
                    f"mean_abs_delta={mean_abs_delta:.12g}",
                    f"rmse={rmse_value:.12g}",
                    f"rel_rmse={rel_rmse:.12g}",
                    f"rel_peak_delta={rel_peak_delta:.12g}",
                    f"rel_final_delta={rel_final_delta:.12g}",
                    f"rel_integral_delta={rel_integral_delta:.12g}",
                    f"rel_steady_mean_delta={rel_steady_mean_delta:.12g}",
                    f"status={status}",
                ]
            )
        )

        summary_rows.append(
            {
                "column": column,
                "samples": str(len(reference_values)),
                "max_abs_delta": f"{max_abs_delta:.12g}",
                "mean_abs_delta": f"{mean_abs_delta:.12g}",
                "rmse": f"{rmse_value:.12g}",
                "rel_rmse": f"{rel_rmse:.12g}",
                "rel_peak_delta": f"{rel_peak_delta:.12g}",
                "rel_final_delta": f"{rel_final_delta:.12g}",
                "rel_integral_delta": f"{rel_integral_delta:.12g}",
                "rel_steady_mean_delta": f"{rel_steady_mean_delta:.12g}",
                "status": status,
            }
        )

    if math.isfinite(args.breakdown_threshold):
        if args.breakdown_column not in compared_columns:
            print(f"breakdown_column_not_compared={args.breakdown_column}")
            failures.append(f"breakdown_column_missing({args.breakdown_column})")
        else:
            reference_breakdown_time = first_crossing_time(
                reference_times,
                reference_series[args.breakdown_column],
                args.breakdown_threshold,
            )
            candidate_breakdown_time = first_crossing_time(
                reference_times,
                aligned_candidate_series[args.breakdown_column],
                args.breakdown_threshold,
            )

            breakdown_status = "PASS"
            breakdown_time_delta = None
            if reference_breakdown_time is None or candidate_breakdown_time is None:
                if args.fail_on_missing_breakdown:
                    breakdown_status = "FAIL"
                    failures.append("missing_breakdown_crossing")
            else:
                breakdown_time_delta = abs(candidate_breakdown_time - reference_breakdown_time)
                if breakdown_time_delta > args.max_breakdown_time_delta_s:
                    breakdown_status = "FAIL"
                    failures.append(
                        f"breakdown_time_delta_s({breakdown_time_delta:.6g})>{args.max_breakdown_time_delta_s:.6g}"
                    )

            breakdown_summary = {
                "column": args.breakdown_column,
                "threshold": f"{args.breakdown_threshold:.12g}",
                "reference_breakdown_time_s": (
                    "none"
                    if reference_breakdown_time is None
                    else f"{reference_breakdown_time:.12g}"
                ),
                "candidate_breakdown_time_s": (
                    "none"
                    if candidate_breakdown_time is None
                    else f"{candidate_breakdown_time:.12g}"
                ),
                "breakdown_time_delta_s": (
                    "none"
                    if breakdown_time_delta is None
                    else f"{breakdown_time_delta:.12g}"
                ),
                "status": breakdown_status,
            }

            print(
                " ".join(
                    [
                        f"breakdown_column={args.breakdown_column}",
                        f"threshold={args.breakdown_threshold:.12g}",
                        f"reference_breakdown_time_s={breakdown_summary['reference_breakdown_time_s']}",
                        f"candidate_breakdown_time_s={breakdown_summary['candidate_breakdown_time_s']}",
                        f"breakdown_time_delta_s={breakdown_summary['breakdown_time_delta_s']}",
                        f"status={breakdown_status}",
                    ]
                )
            )

    if args.summary_csv:
        summary_csv_path = Path(args.summary_csv)
        write_summary_csv(summary_csv_path, summary_rows)
        print(f"summary_csv={summary_csv_path.resolve()}")

    overall_status = "FAIL" if failures else "PASS"

    if args.summary_md:
        summary_md_path = Path(args.summary_md)
        write_summary_markdown(summary_md_path, summary_rows, breakdown_summary, overall_status)
        print(f"summary_md={summary_md_path.resolve()}")

    if args.summary_json:
        summary_json_path = Path(args.summary_json)
        summary_json_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "reference_csv": str(reference_path),
            "candidate_csv": str(candidate_path),
            "compared_columns": compared_columns,
            "thresholds": {
                "max_abs_delta": args.max_abs_delta,
                "max_rel_rmse": args.max_rel_rmse,
                "max_rel_peak_delta": args.max_rel_peak_delta,
                "max_rel_final_delta": args.max_rel_final_delta,
                "max_rel_integral_delta": args.max_rel_integral_delta,
                "max_rel_steady_mean_delta": args.max_rel_steady_mean_delta,
                "breakdown_column": args.breakdown_column,
                "breakdown_threshold": finite_or_inf(args.breakdown_threshold),
                "max_breakdown_time_delta_s": args.max_breakdown_time_delta_s,
            },
            "results": summary_rows,
            "breakdown": breakdown_summary,
            "status": overall_status,
            "failures": failures,
        }
        summary_json_path.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        print(f"summary_json={summary_json_path.resolve()}")

    if failures:
        for failure in failures:
            print(f"failure={failure}")
        print("status=FAIL")
        return 2

    print("status=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
