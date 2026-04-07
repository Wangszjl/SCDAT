#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple


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


def parse_kv_lines(stdout_text: str) -> Dict[str, str]:
    pairs: Dict[str, str] = {}
    for raw_line in stdout_text.splitlines():
        line = raw_line.strip()
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        pairs[key.strip()] = value.strip()
    return pairs


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
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid json at {path}: {exc}") from exc


def collect_case_variant_outputs(
    manifest: Dict[str, Any],
    aligned_variant_id: str,
    legacy_variant_id: str,
) -> Tuple[Dict[str, Dict[str, Path]], List[str]]:
    cases: Dict[str, Dict[str, Path]] = {}
    failures: List[str] = []

    runs = manifest.get("runs", [])
    if not isinstance(runs, list):
        return cases, ["manifest_runs_not_list"]

    for run in runs:
        if not isinstance(run, dict):
            continue
        case_id = str(run.get("case_id", "")).strip()
        variant = str(run.get("variant", "")).strip()
        output_csv = str(run.get("output_csv", "")).strip()
        run_status = str(run.get("status", "")).strip().upper()

        if not case_id or not variant or not output_csv:
            continue

        if run_status not in {"PASS", "DRY_RUN"}:
            failures.append(f"run_failed:{case_id}:{variant}")

        if variant not in {aligned_variant_id, legacy_variant_id}:
            continue

        cases.setdefault(case_id, {})[variant] = Path(output_csv)

    for case_id, variant_map in cases.items():
        if aligned_variant_id not in variant_map:
            failures.append(f"missing_aligned_variant:{case_id}:{aligned_variant_id}")
        if legacy_variant_id not in variant_map:
            failures.append(f"missing_legacy_variant:{case_id}:{legacy_variant_id}")

    return cases, failures


def write_summary_markdown(summary: Dict[str, Any], output_md: Path) -> None:
    output_md.parent.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    lines.append("# Vacuum Arc Matrix Gate Summary")
    lines.append("")
    lines.append(f"- status: {summary.get('status', 'UNKNOWN')}")
    lines.append(f"- timestamp_utc: {summary.get('timestamp_utc', '')}")
    lines.append(f"- matrix_json: {summary.get('matrix_json', '')}")
    lines.append(f"- manifest_json: {summary.get('manifest_json', '')}")
    lines.append("")
    lines.append("| case_id | status | compare_json | aligned_csv | legacy_csv |")
    lines.append("|---|---|---|---|---|")

    for case in summary.get("cases", []):
        lines.append(
            "| "
            + " | ".join(
                [
                    str(case.get("case_id", "")),
                    str(case.get("status", "")),
                    str(case.get("compare_json", "")),
                    str(case.get("aligned_csv", "")),
                    str(case.get("legacy_csv", "")),
                ]
            )
            + " |"
        )

    failures = summary.get("failures", [])
    if failures:
        lines.append("")
        lines.append("## Failures")
        lines.append("")
        for failure in failures:
            lines.append(f"- {failure}")

    output_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="VA-009 matrix gate: run arc case matrix and compare aligned-vs-legacy outputs."
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--scdat-exe", required=True)
    parser.add_argument("--matrix-json", required=True)
    parser.add_argument("--output-root", default="build/arc_compare_matrix_gate")
    parser.add_argument("--aligned-variant-id", default="arc_pic_aligned")
    parser.add_argument("--legacy-variant-id", default="legacy_baseline")
    parser.add_argument("--timestamp", default="")
    parser.add_argument("--columns", nargs="+", default=DEFAULT_COLUMNS)
    parser.add_argument("--allow-missing-columns", action="store_true")
    parser.add_argument("--breakdown-column", default="discharge_current_a")
    parser.add_argument("--breakdown-threshold", type=float, default=1.0e-9)
    parser.add_argument("--fail-on-missing-breakdown", action="store_true")
    parser.add_argument("--max-abs-delta", type=float, default=1.0e-2)
    parser.add_argument("--max-rel-rmse", type=float, default=1.0e-3)
    parser.add_argument("--max-rel-peak-delta", type=float, default=1.0e-3)
    parser.add_argument("--max-rel-final-delta", type=float, default=1.0e-3)
    parser.add_argument("--max-rel-integral-delta", type=float, default=1.0e-3)
    parser.add_argument("--max-rel-steady-mean-delta", type=float, default=1.0e-3)
    parser.add_argument("--max-breakdown-time-delta-s", type=float, default=1.0e-8)
    parser.add_argument("--summary-json", default="")
    parser.add_argument("--summary-md", default="")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    scdat_exe = Path(args.scdat_exe)
    if not scdat_exe.is_absolute():
        scdat_exe = (project_root / scdat_exe).resolve()
    matrix_json = Path(args.matrix_json)
    if not matrix_json.is_absolute():
        matrix_json = (project_root / matrix_json).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = (project_root / output_root).resolve()

    timestamp = args.timestamp.strip() or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    summary_json = Path(args.summary_json) if args.summary_json else output_root / "arc_compare_matrix_gate.json"
    if not summary_json.is_absolute():
        summary_json = (project_root / summary_json).resolve()
    summary_md = Path(args.summary_md) if args.summary_md else output_root / "arc_compare_matrix_gate.md"
    if not summary_md.is_absolute():
        summary_md = (project_root / summary_md).resolve()

    try:
        ensure_file(scdat_exe, "SCDAT executable")
        ensure_file(matrix_json, "matrix json")

        run_matrix_script = project_root / "scripts" / "run" / "run_vacuum_arc_case_matrix.py"
        compare_script = project_root / "scripts" / "python" / "compare_vacuum_arc_variants.py"
        ensure_file(run_matrix_script, "run_vacuum_arc_case_matrix.py")
        ensure_file(compare_script, "compare_vacuum_arc_variants.py")

        output_root.mkdir(parents=True, exist_ok=True)

        matrix_command = [
            sys.executable,
            str(run_matrix_script),
            "--matrix-json",
            str(matrix_json),
            "--project-root",
            str(project_root),
            "--scdat-exe",
            str(scdat_exe),
            "--output-root",
            str(output_root),
            "--timestamp",
            timestamp,
            "--continue-on-failure",
        ]
        matrix_result = run_command(matrix_command, project_root)
        print_command_output("matrix", matrix_result)
        if matrix_result.returncode != 0:
            raise RuntimeError(f"matrix execution failed with code {matrix_result.returncode}")

        matrix_pairs = parse_kv_lines(matrix_result.stdout)
        manifest_json_raw = matrix_pairs.get("manifest_json", "")
        if not manifest_json_raw:
            raise RuntimeError("matrix execution did not report manifest_json")
        manifest_json = Path(manifest_json_raw)
        if not manifest_json.is_absolute():
            manifest_json = (project_root / manifest_json).resolve()
        ensure_file(manifest_json, "matrix manifest")

        manifest = load_json(manifest_json)
        case_variant_outputs, failures = collect_case_variant_outputs(
            manifest,
            args.aligned_variant_id,
            args.legacy_variant_id,
        )

        case_summaries: List[Dict[str, Any]] = []
        for case_id in sorted(case_variant_outputs.keys()):
            variant_outputs = case_variant_outputs[case_id]
            aligned_csv = variant_outputs.get(args.aligned_variant_id)
            legacy_csv = variant_outputs.get(args.legacy_variant_id)
            if not aligned_csv or not legacy_csv:
                continue

            ensure_file(aligned_csv, f"aligned csv for case {case_id}")
            ensure_file(legacy_csv, f"legacy csv for case {case_id}")

            compare_json = output_root / f"compare_{case_id}.json"
            compare_csv = output_root / f"compare_{case_id}.csv"
            compare_md = output_root / f"compare_{case_id}.md"

            compare_command = [
                sys.executable,
                str(compare_script),
                "--reference-csv",
                str(legacy_csv),
                "--candidate-csv",
                str(aligned_csv),
                "--columns",
                *args.columns,
                "--breakdown-column",
                args.breakdown_column,
                "--breakdown-threshold",
                str(args.breakdown_threshold),
                "--max-abs-delta",
                str(args.max_abs_delta),
                "--max-rel-rmse",
                str(args.max_rel_rmse),
                "--max-rel-peak-delta",
                str(args.max_rel_peak_delta),
                "--max-rel-final-delta",
                str(args.max_rel_final_delta),
                "--max-rel-integral-delta",
                str(args.max_rel_integral_delta),
                "--max-rel-steady-mean-delta",
                str(args.max_rel_steady_mean_delta),
                "--max-breakdown-time-delta-s",
                str(args.max_breakdown_time_delta_s),
                "--summary-json",
                str(compare_json),
                "--summary-csv",
                str(compare_csv),
                "--summary-md",
                str(compare_md),
            ]
            if args.allow_missing_columns:
                compare_command.append("--allow-missing-columns")
            if args.fail_on_missing_breakdown:
                compare_command.append("--fail-on-missing-breakdown")

            compare_result = run_command(compare_command, project_root)
            print_command_output(f"compare:{case_id}", compare_result)

            case_status = "PASS" if compare_result.returncode == 0 else "FAIL"
            case_failures: List[str] = []
            if compare_json.is_file():
                compare_payload = load_json(compare_json)
                payload_status = str(compare_payload.get("status", case_status)).upper()
                case_status = "PASS" if payload_status == "PASS" else "FAIL"
                raw_failures = compare_payload.get("failures", [])
                if isinstance(raw_failures, list):
                    case_failures = [str(item) for item in raw_failures]

            if case_status != "PASS":
                failures.append(f"compare_failed:{case_id}")
                failures.extend([f"{case_id}:{item}" for item in case_failures])

            case_summaries.append(
                {
                    "case_id": case_id,
                    "status": case_status,
                    "aligned_csv": str(aligned_csv),
                    "legacy_csv": str(legacy_csv),
                    "compare_json": str(compare_json),
                    "compare_csv": str(compare_csv),
                    "compare_md": str(compare_md),
                    "failures": case_failures,
                }
            )

        overall_status = "PASS" if not failures else "FAIL"
        summary = {
            "status": overall_status,
            "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "project_root": str(project_root),
            "matrix_json": str(matrix_json),
            "manifest_json": str(manifest_json),
            "output_root": str(output_root),
            "thresholds": {
                "max_abs_delta": args.max_abs_delta,
                "max_rel_rmse": args.max_rel_rmse,
                "max_rel_peak_delta": args.max_rel_peak_delta,
                "max_rel_final_delta": args.max_rel_final_delta,
                "max_rel_integral_delta": args.max_rel_integral_delta,
                "max_rel_steady_mean_delta": args.max_rel_steady_mean_delta,
                "max_breakdown_time_delta_s": args.max_breakdown_time_delta_s,
            },
            "cases": case_summaries,
            "failures": failures,
        }

        summary_json.parent.mkdir(parents=True, exist_ok=True)
        summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
        write_summary_markdown(summary, summary_md)

        print(f"summary_json={summary_json}")
        print(f"summary_md={summary_md}")
        print(f"manifest_json={manifest_json}")
        print(f"status={overall_status}")
        return 0 if overall_status == "PASS" else 2

    except Exception as exc:  # pylint: disable=broad-except
        print(f"status=FAIL({exc})")
        return 3


if __name__ == "__main__":
    raise SystemExit(main())