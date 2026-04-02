#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
from pathlib import Path

from summarize_benchmark_report import evaluate_acceptance, parse_report


def resolve_report_paths(patterns: list[str]) -> list[Path]:
    found: set[Path] = set()
    for pattern in patterns:
        for match in glob.glob(pattern):
            path = Path(match).resolve()
            if path.is_file():
                found.add(path)
    return sorted(found)


def render_markdown(
    reports: list[tuple[Path, dict[str, str], dict[str, object]]],
    failures: list[Path],
    strict_not_applicable: bool,
) -> str:
    lines: list[str] = []
    lines.append("# Legacy Benchmark Regression Gate")
    lines.append("")
    lines.append(f"- report_count: `{len(reports)}`")
    lines.append(f"- strict_not_applicable: `{int(strict_not_applicable)}`")
    lines.append(f"- failures: `{len(failures)}`")
    lines.append("")
    lines.append("| report | source | mode | status | checks_failed |")
    lines.append("|---|---|---|---|---:|")
    for path, values, acceptance in reports:
        lines.append(
            "| "
            f"{path.name} | "
            f"{values.get('benchmark_source', 'unknown')} | "
            f"{values.get('execution_mode', 'unknown')} | "
            f"{acceptance['status']} | "
            f"{acceptance['checks_failed']} |"
        )

    if failures:
        lines.append("")
        lines.append("## Failed Reports")
        lines.append("")
        for path in failures:
            lines.append(f"- `{path}`")

    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run acceptance-gate checks for one or many legacy benchmark sidecar reports."
    )
    parser.add_argument(
        "--glob",
        nargs="+",
        default=["results/*.benchmark.txt"],
        help="Glob pattern(s) to locate benchmark report files.",
    )
    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="Exit success when no report files are found.",
    )
    parser.add_argument(
        "--strict-not-applicable",
        action="store_true",
        help="Treat NOT_APPLICABLE reports as failures.",
    )
    parser.add_argument(
        "--report-md",
        default="",
        help="Optional markdown output path for CI artifact.",
    )
    args = parser.parse_args()

    report_paths = resolve_report_paths(args.glob)
    if not report_paths:
        print("report_count=0")
        if args.allow_empty:
            print("status=PASS(empty)")
            return 0
        print("status=FAIL(no_reports)")
        return 3

    records: list[tuple[Path, dict[str, str], dict[str, object]]] = []
    failures: list[Path] = []

    print(f"report_count={len(report_paths)}")
    for report_path in report_paths:
        values = parse_report(report_path)
        acceptance = evaluate_acceptance(values)
        records.append((report_path, values, acceptance))

        status = str(acceptance["status"])
        checks_failed = int(acceptance["checks_failed"])
        benchmark_source = values.get("benchmark_source", "unknown")
        execution_mode = values.get("execution_mode", "unknown")
        print(
            f"report={report_path.name} "
            f"source={benchmark_source} mode={execution_mode} "
            f"status={status} checks_failed={checks_failed}"
        )

        if status == "FAIL":
            failures.append(report_path)
        elif status == "NOT_APPLICABLE" and args.strict_not_applicable:
            failures.append(report_path)

    if args.report_md:
        report_md_path = Path(args.report_md)
        report_md_path.parent.mkdir(parents=True, exist_ok=True)
        report_md_path.write_text(
            render_markdown(records, failures, args.strict_not_applicable),
            encoding="utf-8",
        )
        print(f"markdown_report={report_md_path}")

    if failures:
        print(f"status=FAIL(failed_reports={len(failures)})")
        return 2

    print("status=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
