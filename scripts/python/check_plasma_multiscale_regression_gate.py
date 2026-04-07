#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

SUMMARY_RE = re.compile(
    r"PLASMA_BENCHMARK_SUMMARY\s+case_id=(?P<case_id>\S+)\s+fluid=(?P<fluid>\d+)\s+hybrid=(?P<hybrid>\d+)\s+pic=(?P<pic>\d+)\s+signature=(?P<signature>\S+)"
)


def load_json(path: Path) -> Dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid json at {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"expected object json at {path}")
    return payload


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Plasma Multiscale Regression Gate")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- matrix_json: {report['matrix_json']}")
    lines.append(f"- baseline_json: {report['baseline_json']}")
    lines.append("")
    lines.append("| case_id | fluid | hybrid | pic | signature | status |")
    lines.append("|---|---:|---:|---:|---|---|")
    for case in report["cases"]:
        lines.append(
            f"| {case['case_id']} | {case.get('fluid_macro_cells', '-') } | {case.get('hybrid_transition_cells', '-')} | {case.get('local_pic_cells', '-')} | {case.get('reproducibility_signature', '-')} | {case['status']} |"
        )

    if report["failures"]:
        lines.append("")
        lines.append("## Failures")
        lines.append("")
        for failure in report["failures"]:
            lines.append(f"- {failure}")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_case(exe_path: Path, case_id: str, gtest_filter: str, cwd: Path) -> Dict[str, Any]:
    command = [str(exe_path), f"--gtest_filter={gtest_filter}", "--gtest_color=no"]
    completed = subprocess.run(
        command,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )

    combined_output = (completed.stdout or "") + "\n" + (completed.stderr or "")
    summary_match = SUMMARY_RE.search(combined_output)
    result: Dict[str, Any] = {
        "case_id": case_id,
        "gtest_filter": gtest_filter,
        "return_code": completed.returncode,
        "status": "passed" if completed.returncode == 0 else "failed",
    }

    if summary_match:
        result["reported_case_id"] = summary_match.group("case_id")
        result["fluid_macro_cells"] = int(summary_match.group("fluid"))
        result["hybrid_transition_cells"] = int(summary_match.group("hybrid"))
        result["local_pic_cells"] = int(summary_match.group("pic"))
        result["reproducibility_signature"] = summary_match.group("signature")

    result["stdout"] = completed.stdout
    result["stderr"] = completed.stderr
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="PA-008: Plasma multi-scale decomposition regression gate"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--plasma-test-exe", required=True)
    parser.add_argument("--matrix-json", required=True)
    parser.add_argument("--baseline-json", required=True)
    parser.add_argument("--summary-json", default="")
    parser.add_argument("--summary-md", default="")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()

    test_exe = Path(args.plasma_test_exe)
    if not test_exe.is_absolute():
        test_exe = (project_root / test_exe).resolve()

    matrix_json = Path(args.matrix_json)
    if not matrix_json.is_absolute():
        matrix_json = (project_root / matrix_json).resolve()

    baseline_json = Path(args.baseline_json)
    if not baseline_json.is_absolute():
        baseline_json = (project_root / baseline_json).resolve()

    summary_json = Path(args.summary_json) if args.summary_json else project_root / "build/plasma_multiscale_regression_gate.json"
    if not summary_json.is_absolute():
        summary_json = (project_root / summary_json).resolve()

    summary_md = Path(args.summary_md) if args.summary_md else project_root / "build/plasma_multiscale_regression_gate.md"
    if not summary_md.is_absolute():
        summary_md = (project_root / summary_md).resolve()

    matrix = load_json(matrix_json)
    baseline = load_json(baseline_json)

    benchmarks = matrix.get("benchmarks", [])
    if not isinstance(benchmarks, list) or not benchmarks:
        raise ValueError("matrix_json must contain non-empty 'benchmarks' array")

    baseline_cases = baseline.get("cases", {})
    if not isinstance(baseline_cases, dict):
        raise ValueError("baseline_json must contain object field 'cases'")

    failures: List[str] = []
    case_reports: List[Dict[str, Any]] = []

    for benchmark in benchmarks:
        if not isinstance(benchmark, dict):
            failures.append("invalid benchmark entry: expected object")
            continue

        case_id = str(benchmark.get("case_id", "")).strip()
        gtest_filter = str(benchmark.get("gtest_filter", "")).strip()
        if not case_id or not gtest_filter:
            failures.append("benchmark entry requires case_id and gtest_filter")
            continue

        report = run_case(test_exe, case_id, gtest_filter, project_root)
        case_reports.append(report)

        if report["return_code"] != 0:
            failures.append(f"{case_id}: test command failed with code {report['return_code']}")
            continue

        if "reported_case_id" not in report:
            failures.append(f"{case_id}: missing PLASMA_BENCHMARK_SUMMARY output")
            continue

        if report["reported_case_id"] != case_id:
            failures.append(
                f"{case_id}: reported case_id {report['reported_case_id']} mismatches benchmark id"
            )
            continue

        expected = baseline_cases.get(case_id)
        if not isinstance(expected, dict):
            failures.append(f"{case_id}: missing baseline entry")
            continue

        expected_fluid = int(expected.get("fluid_macro_cells", -1))
        expected_hybrid = int(expected.get("hybrid_transition_cells", -1))
        expected_pic = int(expected.get("local_pic_cells", -1))
        expected_signature = str(expected.get("reproducibility_signature", ""))

        if report["fluid_macro_cells"] != expected_fluid:
            failures.append(
                f"{case_id}: fluid mismatch {report['fluid_macro_cells']} != {expected_fluid}"
            )
        if report["hybrid_transition_cells"] != expected_hybrid:
            failures.append(
                f"{case_id}: hybrid mismatch {report['hybrid_transition_cells']} != {expected_hybrid}"
            )
        if report["local_pic_cells"] != expected_pic:
            failures.append(
                f"{case_id}: pic mismatch {report['local_pic_cells']} != {expected_pic}"
            )
        if report["reproducibility_signature"] != expected_signature:
            failures.append(
                f"{case_id}: signature mismatch {report['reproducibility_signature']} != {expected_signature}"
            )

    status = "PASS" if not failures else "FAIL"
    output: Dict[str, Any] = {
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "matrix_json": str(matrix_json),
        "baseline_json": str(baseline_json),
        "cases": case_reports,
        "failures": failures,
    }

    write_json(summary_json, output)
    write_markdown(output, summary_md)

    print(f"[PlasmaMultiscaleRegressionGate] status={status}")
    print(f"[PlasmaMultiscaleRegressionGate] summary_json={summary_json}")
    print(f"[PlasmaMultiscaleRegressionGate] summary_md={summary_md}")

    if failures:
        for failure in failures:
            print(f"[PlasmaMultiscaleRegressionGate] failure: {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[PlasmaMultiscaleRegressionGate] fatal: {exc}", file=sys.stderr)
        raise SystemExit(2)
