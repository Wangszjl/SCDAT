#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import glob
import math
from dataclasses import dataclass
from pathlib import Path

DEFAULT_CSV_PATTERNS = [
    "results/tmp_geo_legacy.csv",
    "results/tmp_leo_ram_legacy.csv",
    "results/tmp_leo_wake_legacy.csv",
]

DEFAULT_PAIRS = [
    ("Vs", "benchmark_vs_v"),
    ("Time_ms", "benchmark_time_ms"),
    ("Jb", "benchmark_jb_na_per_m2"),
    ("Jcond", "benchmark_jcond_na_per_m2"),
    ("Je", "benchmark_je_na_per_m2"),
    ("Ji", "benchmark_ji_na_per_m2"),
    ("Jnet", "benchmark_jnet_na_per_m2"),
    ("Jph", "benchmark_jph_na_per_m2"),
    ("Jse", "benchmark_jse_na_per_m2"),
    ("Jsi", "benchmark_jsi_na_per_m2"),
]


@dataclass
class PairMetric:
    left: str
    right: str
    sample_count: int
    max_abs_delta: float
    mean_abs_delta: float
    status: str


@dataclass
class FileMetric:
    path: Path
    row_count: int
    metrics: list[PairMetric]
    status: str
    note: str


def resolve_csv_paths(patterns: list[str]) -> list[Path]:
    found: set[Path] = set()
    for pattern in patterns:
        matches = [Path(match) for match in glob.glob(pattern)]
        if not matches:
            candidate = Path(pattern)
            if candidate.is_file():
                matches = [candidate]
        for match in matches:
            resolved = match.resolve()
            if resolved.is_file():
                found.add(resolved)
    return sorted(found)


def parse_pair_specs(pair_specs: list[str]) -> list[tuple[str, str]]:
    parsed: list[tuple[str, str]] = []
    for raw_pair in pair_specs:
        if ":" not in raw_pair:
            raise ValueError(f"invalid pair spec '{raw_pair}', expected format lhs:rhs")
        left, right = raw_pair.split(":", 1)
        left = left.strip()
        right = right.strip()
        if not left or not right:
            raise ValueError(f"invalid pair spec '{raw_pair}', empty lhs/rhs")
        parsed.append((left, right))
    return parsed


def try_parse_float(raw_value: str) -> float | None:
    try:
        value = float(raw_value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def evaluate_csv(
    csv_path: Path,
    pairs: list[tuple[str, str]],
    max_abs_delta: float,
    strict_schema: bool,
) -> FileMetric:
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return FileMetric(csv_path, 0, [], "FAIL", "missing header")

        fieldnames = set(reader.fieldnames)
        available_pairs = [pair for pair in pairs if pair[0] in fieldnames and pair[1] in fieldnames]
        missing_pairs = [pair for pair in pairs if pair not in available_pairs]

        if strict_schema and missing_pairs:
            missing = ", ".join(f"{left}:{right}" for left, right in missing_pairs)
            return FileMetric(csv_path, 0, [], "FAIL", f"missing columns for pairs: {missing}")

        if not available_pairs:
            return FileMetric(csv_path, 0, [], "FAIL", "no comparable column pairs")

        accumulators: dict[tuple[str, str], dict[str, float | int]] = {
            pair: {
                "sample_count": 0,
                "sum_abs_delta": 0.0,
                "max_abs_delta": 0.0,
            }
            for pair in available_pairs
        }

        row_count = 0
        for row in reader:
            row_count += 1
            for pair in available_pairs:
                left_raw = row.get(pair[0], "")
                right_raw = row.get(pair[1], "")
                left_value = try_parse_float(left_raw)
                right_value = try_parse_float(right_raw)
                if left_value is None or right_value is None:
                    continue

                delta = abs(left_value - right_value)
                bucket = accumulators[pair]
                bucket["sample_count"] = int(bucket["sample_count"]) + 1
                bucket["sum_abs_delta"] = float(bucket["sum_abs_delta"]) + delta
                bucket["max_abs_delta"] = max(float(bucket["max_abs_delta"]), delta)

    pair_metrics: list[PairMetric] = []
    has_failure = False
    for left, right in available_pairs:
        bucket = accumulators[(left, right)]
        sample_count = int(bucket["sample_count"])
        if sample_count == 0:
            pair_metrics.append(
                PairMetric(left, right, 0, math.inf, math.inf, "FAIL")
            )
            has_failure = True
            continue

        max_delta = float(bucket["max_abs_delta"])
        mean_delta = float(bucket["sum_abs_delta"]) / sample_count
        status = "PASS" if max_delta <= max_abs_delta else "FAIL"
        if status == "FAIL":
            has_failure = True
        pair_metrics.append(
            PairMetric(left, right, sample_count, max_delta, mean_delta, status)
        )

    note = ""
    if missing_pairs and not strict_schema:
        note = "ignored missing pairs: " + ", ".join(
            f"{left}:{right}" for left, right in missing_pairs
        )

    return FileMetric(csv_path, row_count, pair_metrics, "FAIL" if has_failure else "PASS", note)


def render_markdown_report(results: list[FileMetric], max_abs_delta: float) -> str:
    lines: list[str] = []
    lines.append("# Legacy Replay CSV Consistency Gate")
    lines.append("")
    lines.append(f"- csv_count: `{len(results)}`")
    lines.append(f"- max_abs_delta: `{max_abs_delta:.3e}`")
    lines.append("")
    lines.append("| csv | row_count | pair | max_abs_delta | mean_abs_delta | status |")
    lines.append("|---|---:|---|---:|---:|---|")

    for file_result in results:
        if not file_result.metrics:
            lines.append(
                f"| {file_result.path.name} | {file_result.row_count} | n/a | n/a | n/a | {file_result.status} |"
            )
            continue

        for metric in file_result.metrics:
            lines.append(
                "| "
                f"{file_result.path.name} | "
                f"{file_result.row_count} | "
                f"{metric.left}->{metric.right} | "
                f"{metric.max_abs_delta:.6e} | "
                f"{metric.mean_abs_delta:.6e} | "
                f"{metric.status} |"
            )
        if file_result.note:
            lines.append(
                f"| {file_result.path.name} | {file_result.row_count} | note | n/a | n/a | {file_result.note} |"
            )

    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate legacy replay CSV consistency against embedded benchmark columns."
    )
    parser.add_argument(
        "--csv",
        nargs="+",
        default=DEFAULT_CSV_PATTERNS,
        help="CSV path(s) or glob pattern(s) to validate.",
    )
    parser.add_argument(
        "--pair",
        action="append",
        default=[],
        help="Additional comparison pair in lhs:rhs format.",
    )
    parser.add_argument(
        "--max-abs-delta",
        type=float,
        default=1.0e-12,
        help="Maximum allowed absolute delta for each compared pair.",
    )
    parser.add_argument(
        "--strict-schema",
        action="store_true",
        help="Fail when any requested comparison pair is missing from a CSV.",
    )
    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="Return success when no CSV files are discovered.",
    )
    parser.add_argument(
        "--report-md",
        default="",
        help="Optional markdown report output path.",
    )
    args = parser.parse_args()

    pairs = list(DEFAULT_PAIRS)
    if args.pair:
        pairs.extend(parse_pair_specs(args.pair))

    csv_paths = resolve_csv_paths(args.csv)
    if not csv_paths:
        print("csv_count=0")
        if args.allow_empty:
            print("status=PASS(empty)")
            return 0
        print("status=FAIL(no_csv)")
        return 3

    print(f"csv_count={len(csv_paths)}")

    results: list[FileMetric] = []
    failures: list[Path] = []

    for csv_path in csv_paths:
        file_result = evaluate_csv(csv_path, pairs, args.max_abs_delta, args.strict_schema)
        results.append(file_result)

        print(
            f"csv={csv_path.name} row_count={file_result.row_count} "
            f"pair_count={len(file_result.metrics)} status={file_result.status}"
        )
        if file_result.note:
            print(f"note={file_result.note}")

        for metric in file_result.metrics:
            print(
                f"pair={metric.left}->{metric.right} "
                f"samples={metric.sample_count} "
                f"max_abs_delta={metric.max_abs_delta:.6e} "
                f"mean_abs_delta={metric.mean_abs_delta:.6e} "
                f"status={metric.status}"
            )

        if file_result.status != "PASS":
            failures.append(csv_path)

    if args.report_md:
        report_path = Path(args.report_md)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(
            render_markdown_report(results, args.max_abs_delta),
            encoding="utf-8",
        )
        print(f"markdown_report={report_path.resolve()}")

    if failures:
        print(f"status=FAIL(failed_csv={len(failures)})")
        return 2

    print("status=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
