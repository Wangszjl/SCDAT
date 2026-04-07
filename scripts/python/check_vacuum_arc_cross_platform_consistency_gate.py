#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


DEFAULT_COLUMNS = [
    "discharge_current_a",
    "current_density_a_per_m2",
    "channel_conductivity_s_per_m",
    "cathode_temperature_k",
    "anode_temperature_k",
    "surface_potential_v",
    "surface_charge_density_c_per_m2",
    "stability_substeps_used",
    "stability_rollbacks",
]


@dataclass
class Thresholds:
    max_rel_rmse_delta: float
    max_rel_peak_delta_delta: float
    max_rel_final_delta_delta: float
    max_rel_integral_delta_delta: float
    max_rel_steady_mean_delta_delta: float
    max_breakdown_time_delta_delta_s: float


@dataclass
class PlatformCase:
    case_id: str
    status: str
    compare_json: Path
    compare_payload: Dict[str, Any]


def run_command(command: List[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )


def print_command_output(prefix: str, completed: subprocess.CompletedProcess[str]) -> None:
    if completed.stdout.strip():
        print(f"[{prefix}] stdout:")
        print(completed.stdout.strip())
    if completed.stderr.strip():
        print(f"[{prefix}] stderr:")
        print(completed.stderr.strip())


def ensure_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"missing {description}: {path}")


def load_json(path: Path) -> Dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid json at {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"expected object json at {path}")
    return payload


def parse_float_or_none(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"", "none", "nan", "null"}:
            return None
        try:
            return float(text)
        except ValueError:
            return None
    return None


def resolve_compare_json(summary_path: Path, raw_compare_json: str) -> Path:
    compare_path = Path(raw_compare_json)
    if compare_path.is_absolute():
        return compare_path

    candidate_1 = (summary_path.parent / compare_path).resolve()
    if candidate_1.is_file():
        return candidate_1

    candidate_2 = compare_path.resolve()
    return candidate_2


def load_platform_cases(summary_json: Path, platform_label: str) -> Dict[str, PlatformCase]:
    payload = load_json(summary_json)
    if str(payload.get("status", "")).upper() != "PASS":
        raise ValueError(f"{platform_label} summary status is not PASS: {summary_json}")

    cases_raw = payload.get("cases", [])
    if not isinstance(cases_raw, list):
        raise ValueError(f"{platform_label} summary cases is not list: {summary_json}")

    result: Dict[str, PlatformCase] = {}
    for case in cases_raw:
        if not isinstance(case, dict):
            continue
        case_id = str(case.get("case_id", "")).strip()
        status = str(case.get("status", "")).strip().upper()
        compare_json_raw = str(case.get("compare_json", "")).strip()
        if not case_id or not compare_json_raw:
            continue

        compare_json = resolve_compare_json(summary_json, compare_json_raw)
        ensure_file(compare_json, f"{platform_label} compare json for case {case_id}")
        compare_payload = load_json(compare_json)

        result[case_id] = PlatformCase(
            case_id=case_id,
            status=status,
            compare_json=compare_json,
            compare_payload=compare_payload,
        )

    if not result:
        raise ValueError(f"no valid cases found in {platform_label} summary: {summary_json}")
    return result


def extract_column_metrics(compare_payload: Dict[str, Any]) -> Dict[str, Dict[str, float]]:
    results = compare_payload.get("results", [])
    if not isinstance(results, list):
        raise ValueError("compare payload results must be list")

    metrics_by_column: Dict[str, Dict[str, float]] = {}
    for item in results:
        if not isinstance(item, dict):
            continue
        column = str(item.get("column", "")).strip()
        if not column:
            continue

        column_metrics: Dict[str, float] = {}
        for key in [
            "rel_rmse",
            "rel_peak_delta",
            "rel_final_delta",
            "rel_integral_delta",
            "rel_steady_mean_delta",
        ]:
            value = parse_float_or_none(item.get(key))
            if value is not None:
                column_metrics[key] = value

        metrics_by_column[column] = column_metrics

    return metrics_by_column


def compare_case_metrics(
    case_id: str,
    a_payload: Dict[str, Any],
    b_payload: Dict[str, Any],
    thresholds: Thresholds,
) -> Tuple[List[str], List[Dict[str, Any]]]:
    failures: List[str] = []
    column_reports: List[Dict[str, Any]] = []

    a_metrics = extract_column_metrics(a_payload)
    b_metrics = extract_column_metrics(b_payload)
    common_columns = sorted(set(a_metrics.keys()) & set(b_metrics.keys()))
    missing_columns = sorted((set(a_metrics.keys()) ^ set(b_metrics.keys())))
    for column in missing_columns:
        failures.append(f"{case_id}:missing_column_across_platforms:{column}")

    limits = {
        "rel_rmse": thresholds.max_rel_rmse_delta,
        "rel_peak_delta": thresholds.max_rel_peak_delta_delta,
        "rel_final_delta": thresholds.max_rel_final_delta_delta,
        "rel_integral_delta": thresholds.max_rel_integral_delta_delta,
        "rel_steady_mean_delta": thresholds.max_rel_steady_mean_delta_delta,
    }

    for column in common_columns:
        row = {
            "column": column,
            "deltas": {},
            "status": "PASS",
        }
        for metric_key, max_allowed in limits.items():
            a_value = a_metrics[column].get(metric_key)
            b_value = b_metrics[column].get(metric_key)
            if a_value is None or b_value is None:
                continue
            delta = abs(a_value - b_value)
            row["deltas"][metric_key] = delta
            if delta > max_allowed:
                row["status"] = "FAIL"
                failures.append(
                    f"{case_id}:{column}:{metric_key}_delta({delta:.6e})>{max_allowed:.6e}"
                )
        column_reports.append(row)

    a_breakdown = a_payload.get("breakdown", {})
    b_breakdown = b_payload.get("breakdown", {})
    if isinstance(a_breakdown, dict) and isinstance(b_breakdown, dict):
        a_time_delta = parse_float_or_none(a_breakdown.get("breakdown_time_delta_s"))
        b_time_delta = parse_float_or_none(b_breakdown.get("breakdown_time_delta_s"))
        if a_time_delta is not None and b_time_delta is not None:
            dd = abs(a_time_delta - b_time_delta)
            if dd > thresholds.max_breakdown_time_delta_delta_s:
                failures.append(
                    f"{case_id}:breakdown_time_delta_delta_s({dd:.6e})>{thresholds.max_breakdown_time_delta_delta_s:.6e}"
                )

    return failures, column_reports


def write_markdown(report: Dict[str, Any], output_md: Path) -> None:
    output_md.parent.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    lines.append("# Vacuum Arc Cross-Platform Consistency Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- platform_a: {report['platform_a']['label']}")
    lines.append(f"- platform_b: {report['platform_b']['label']}")
    lines.append(f"- platform_a_summary: {report['platform_a']['summary_json']}")
    lines.append(f"- platform_b_summary: {report['platform_b']['summary_json']}")
    lines.append("")

    for case in report["cases"]:
        lines.append(f"## Case: {case['case_id']}")
        lines.append("")
        lines.append(f"- status: {case['status']}")
        lines.append("| column | rel_rmse_delta | rel_peak_delta_delta | rel_final_delta_delta | rel_integral_delta_delta | rel_steady_mean_delta_delta | status |")
        lines.append("|---|---:|---:|---:|---:|---:|---|")
        for row in case["columns"]:
            deltas = row.get("deltas", {})
            lines.append(
                "| "
                + " | ".join(
                    [
                        row.get("column", ""),
                        f"{deltas.get('rel_rmse', 0.0):.6e}",
                        f"{deltas.get('rel_peak_delta', 0.0):.6e}",
                        f"{deltas.get('rel_final_delta', 0.0):.6e}",
                        f"{deltas.get('rel_integral_delta', 0.0):.6e}",
                        f"{deltas.get('rel_steady_mean_delta', 0.0):.6e}",
                        row.get("status", "UNKNOWN"),
                    ]
                )
                + " |"
            )
        lines.append("")

    if report["failures"]:
        lines.append("## Failures")
        lines.append("")
        for failure in report["failures"]:
            lines.append(f"- {failure}")
        lines.append("")

    output_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def run_matrix_gate(
    project_root: Path,
    scdat_exe: Path,
    matrix_json: Path,
    output_root: Path,
    columns: List[str],
) -> Path:
    gate_script = project_root / "scripts" / "python" / "check_vacuum_arc_matrix_gate.py"
    ensure_file(gate_script, "check_vacuum_arc_matrix_gate.py")

    summary_json = output_root / "arc_compare_matrix_gate.json"
    summary_md = output_root / "arc_compare_matrix_gate.md"

    command = [
        sys.executable,
        str(gate_script),
        "--project-root",
        str(project_root),
        "--scdat-exe",
        str(scdat_exe),
        "--matrix-json",
        str(matrix_json),
        "--output-root",
        str(output_root),
        "--columns",
        *columns,
        "--breakdown-column",
        "discharge_current_a",
        "--breakdown-threshold",
        "1e-9",
        "--max-abs-delta",
        "1e-2",
        "--max-rel-rmse",
        "1e-3",
        "--max-rel-peak-delta",
        "1e-3",
        "--max-rel-final-delta",
        "1e-3",
        "--max-rel-integral-delta",
        "1e-3",
        "--max-rel-steady-mean-delta",
        "1e-3",
        "--max-breakdown-time-delta-s",
        "1e-8",
        "--summary-json",
        str(summary_json),
        "--summary-md",
        str(summary_md),
    ]

    completed = run_command(command, project_root)
    print_command_output(f"matrix-gate:{output_root.name}", completed)
    if completed.returncode != 0:
        raise RuntimeError(f"matrix gate failed for {output_root} with code {completed.returncode}")

    ensure_file(summary_json, f"matrix summary json at {output_root}")
    return summary_json


def main() -> int:
    parser = argparse.ArgumentParser(
        description="VA-010 cross-platform consistency gate for Vacuum Arc matrix reports"
    )
    parser.add_argument("--platform-a-label", default="windows")
    parser.add_argument("--platform-b-label", default="linux")
    parser.add_argument("--platform-a-summary-json", default="")
    parser.add_argument("--platform-b-summary-json", default="")

    parser.add_argument("--synthesize-platform-summaries", action="store_true")
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--scdat-exe", default="")
    parser.add_argument("--matrix-json", default="")
    parser.add_argument("--synth-output-root", default="build/arc_cross_platform_consistency_gate")
    parser.add_argument("--columns", nargs="+", default=DEFAULT_COLUMNS)

    parser.add_argument("--max-rel-rmse-delta", type=float, default=1e-3)
    parser.add_argument("--max-rel-peak-delta-delta", type=float, default=1e-3)
    parser.add_argument("--max-rel-final-delta-delta", type=float, default=1e-3)
    parser.add_argument("--max-rel-integral-delta-delta", type=float, default=1e-3)
    parser.add_argument("--max-rel-steady-mean-delta-delta", type=float, default=1e-3)
    parser.add_argument("--max-breakdown-time-delta-delta-s", type=float, default=1e-8)

    parser.add_argument("--report-json", default="")
    parser.add_argument("--report-md", default="")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    synth_output_root = Path(args.synth_output_root)
    if not synth_output_root.is_absolute():
        synth_output_root = (project_root / synth_output_root).resolve()

    platform_a_summary = Path(args.platform_a_summary_json) if args.platform_a_summary_json else None
    platform_b_summary = Path(args.platform_b_summary_json) if args.platform_b_summary_json else None

    if args.synthesize_platform_summaries:
        if not args.scdat_exe or not args.matrix_json:
            print("status=FAIL(missing_synthesis_inputs:scdat_exe/matrix_json)")
            return 3

        scdat_exe = Path(args.scdat_exe)
        if not scdat_exe.is_absolute():
            scdat_exe = (project_root / scdat_exe).resolve()
        matrix_json = Path(args.matrix_json)
        if not matrix_json.is_absolute():
            matrix_json = (project_root / matrix_json).resolve()

        try:
            ensure_file(scdat_exe, "SCDAT executable")
            ensure_file(matrix_json, "matrix json")
            platform_a_summary = run_matrix_gate(
                project_root,
                scdat_exe,
                matrix_json,
                synth_output_root / args.platform_a_label,
                args.columns,
            )
            platform_b_summary = run_matrix_gate(
                project_root,
                scdat_exe,
                matrix_json,
                synth_output_root / args.platform_b_label,
                args.columns,
            )
        except Exception as exc:  # pylint: disable=broad-except
            print(f"status=FAIL({exc})")
            return 3

    if platform_a_summary is None or platform_b_summary is None:
        print("status=FAIL(missing_platform_summary_json)")
        return 3

    if not platform_a_summary.is_absolute():
        platform_a_summary = (project_root / platform_a_summary).resolve()
    if not platform_b_summary.is_absolute():
        platform_b_summary = (project_root / platform_b_summary).resolve()

    report_json = Path(args.report_json) if args.report_json else synth_output_root / "vacuum_arc_cross_platform_consistency.json"
    if not report_json.is_absolute():
        report_json = (project_root / report_json).resolve()
    report_md = Path(args.report_md) if args.report_md else synth_output_root / "vacuum_arc_cross_platform_consistency.md"
    if not report_md.is_absolute():
        report_md = (project_root / report_md).resolve()

    thresholds = Thresholds(
        max_rel_rmse_delta=args.max_rel_rmse_delta,
        max_rel_peak_delta_delta=args.max_rel_peak_delta_delta,
        max_rel_final_delta_delta=args.max_rel_final_delta_delta,
        max_rel_integral_delta_delta=args.max_rel_integral_delta_delta,
        max_rel_steady_mean_delta_delta=args.max_rel_steady_mean_delta_delta,
        max_breakdown_time_delta_delta_s=args.max_breakdown_time_delta_delta_s,
    )

    try:
        ensure_file(platform_a_summary, f"{args.platform_a_label} summary json")
        ensure_file(platform_b_summary, f"{args.platform_b_label} summary json")

        cases_a = load_platform_cases(platform_a_summary, args.platform_a_label)
        cases_b = load_platform_cases(platform_b_summary, args.platform_b_label)

        failures: List[str] = []
        case_reports: List[Dict[str, Any]] = []

        missing_in_b = sorted(set(cases_a.keys()) - set(cases_b.keys()))
        missing_in_a = sorted(set(cases_b.keys()) - set(cases_a.keys()))
        for case_id in missing_in_b:
            failures.append(f"missing_case_in_{args.platform_b_label}:{case_id}")
        for case_id in missing_in_a:
            failures.append(f"missing_case_in_{args.platform_a_label}:{case_id}")

        common_cases = sorted(set(cases_a.keys()) & set(cases_b.keys()))
        for case_id in common_cases:
            case_a = cases_a[case_id]
            case_b = cases_b[case_id]

            case_failure_list: List[str] = []
            if case_a.status != "PASS":
                case_failure_list.append(
                    f"{case_id}:{args.platform_a_label}_status_not_pass:{case_a.status}"
                )
            if case_b.status != "PASS":
                case_failure_list.append(
                    f"{case_id}:{args.platform_b_label}_status_not_pass:{case_b.status}"
                )

            metric_failures, column_reports = compare_case_metrics(
                case_id,
                case_a.compare_payload,
                case_b.compare_payload,
                thresholds,
            )
            case_failure_list.extend(metric_failures)

            status = "PASS" if not case_failure_list else "FAIL"
            failures.extend(case_failure_list)
            case_reports.append(
                {
                    "case_id": case_id,
                    "status": status,
                    "platform_a_compare_json": str(case_a.compare_json),
                    "platform_b_compare_json": str(case_b.compare_json),
                    "columns": column_reports,
                    "failures": case_failure_list,
                }
            )

        overall_status = "PASS" if not failures else "FAIL"
        report = {
            "status": overall_status,
            "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "platform_a": {
                "label": args.platform_a_label,
                "summary_json": str(platform_a_summary),
            },
            "platform_b": {
                "label": args.platform_b_label,
                "summary_json": str(platform_b_summary),
            },
            "thresholds": {
                "max_rel_rmse_delta": thresholds.max_rel_rmse_delta,
                "max_rel_peak_delta_delta": thresholds.max_rel_peak_delta_delta,
                "max_rel_final_delta_delta": thresholds.max_rel_final_delta_delta,
                "max_rel_integral_delta_delta": thresholds.max_rel_integral_delta_delta,
                "max_rel_steady_mean_delta_delta": thresholds.max_rel_steady_mean_delta_delta,
                "max_breakdown_time_delta_delta_s": thresholds.max_breakdown_time_delta_delta_s,
            },
            "cases": case_reports,
            "failures": failures,
        }

        write_json(report_json, report)
        write_markdown(report, report_md)

        print(f"report_json={report_json}")
        print(f"report_md={report_md}")
        print(f"status={overall_status}")
        return 0 if overall_status == "PASS" else 2

    except Exception as exc:  # pylint: disable=broad-except
        print(f"status=FAIL({exc})")
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
