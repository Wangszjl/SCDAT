#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


def load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected top-level object in {path}")
    return payload


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


def load_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def metadata_sidecar(path: Path) -> Path:
    return Path(str(path) + ".metadata.json")


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Vacuum Arc Benchmark Matrix Gate Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- output_root: {report['output_root']}")
    lines.append("")
    lines.append("| preset | mode | command | metrics | metadata | contracts | benchmark_artifact | status |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for item in report["cases"]:
        lines.append(
            f"| {item['preset']} | {item['alignment_mode']} | {item['command_status']} | {item['metrics_status']} | {item['metadata_status']} | {item['contracts_status']} | {item['benchmark_artifact_status']} | {item['status']} |"
        )

    if report["failures"]:
        lines.append("")
        lines.append("## Failures")
        lines.append("")
        for failure in report["failures"]:
            lines.append(f"- {failure}")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the Vacuum Arc ArcPIC benchmark matrix gate."
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--scdat-exe", required=True)
    parser.add_argument(
        "--matrix",
        default="scripts/run/vacuum_arc_benchmark_matrix.json",
    )
    parser.add_argument("--output-root", default="build/vacuum_arc_benchmark_matrix_gate")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    scdat_exe = Path(args.scdat_exe)
    if not scdat_exe.is_absolute():
        scdat_exe = (project_root / scdat_exe).resolve()
    matrix_path = Path(args.matrix)
    if not matrix_path.is_absolute():
        matrix_path = (project_root / matrix_path).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = (project_root / output_root).resolve()

    if not scdat_exe.is_file():
        raise FileNotFoundError(f"missing SCDAT executable: {scdat_exe}")
    if not matrix_path.is_file():
        raise FileNotFoundError(f"missing matrix file: {matrix_path}")

    matrix = load_json(matrix_path)
    cases = matrix.get("cases", [])
    if not isinstance(cases, list) or not cases:
        raise ValueError(f"matrix must contain non-empty cases list: {matrix_path}")

    output_root.mkdir(parents=True, exist_ok=True)
    reports: List[Dict[str, Any]] = []
    failures: List[str] = []

    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            failures.append(f"case[{index}]:invalid_case_type")
            continue
        preset = str(case.get("preset", "")).strip()
        alignment_mode = str(case.get("alignment_mode", "")).strip()
        min_peak_current_a = float(case.get("min_peak_current_a", 0.0))
        min_channel_conductivity_s_per_m = float(
            case.get("min_channel_conductivity_s_per_m", 0.0)
        )
        max_charge_conservation_error_c = float(
            case.get("max_charge_conservation_error_c", 1.0)
        )
        require_collision_fallback_clear = bool(
            case.get("require_collision_fallback_clear", False)
        )
        if not preset or not alignment_mode:
            failures.append(f"case[{index}]:missing_preset_or_alignment_mode")
            continue

        csv_path = output_root / f"{preset}_{alignment_mode}.csv"
        config_path = output_root / f"{preset}_{alignment_mode}.arc_gate.json"
        payload = {
            "schema_version": "v1",
            "module": "arc",
            "base_preset": preset,
            "run": {"steps": 16, "time_step_s": 1.0e-9, "output_csv": str(csv_path)},
            "config": {
                "alignment_mode": alignment_mode,
                "enable_pic_mcc_collisions": True,
                "enable_collision_reaction_breakdown": True
            }
        }
        config_path.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

        completed = run_command([str(scdat_exe), "arc-config", str(config_path), str(csv_path)], project_root)
        command_status = "PASS" if completed.returncode == 0 else "FAIL"
        metrics_status = "FAIL"
        metadata_status = "FAIL"
        contracts_status = "FAIL"
        benchmark_artifact_status = "FAIL"
        status = "FAIL"

        if completed.returncode != 0:
            failures.append(
                f"{preset}:{alignment_mode}:command_failed:code={completed.returncode}:stderr={completed.stderr.strip()}"
            )
        elif not csv_path.is_file():
            failures.append(f"{preset}:{alignment_mode}:missing_csv:{csv_path}")
        else:
            rows = load_csv_rows(csv_path)
            peak_current = 0.0
            max_channel_conductivity = 0.0
            max_abs_charge_conservation_error = 0.0
            final_surface_potential = 0.0
            reaction_split_ok = False
            collision_fallback_clear = True
            for row in rows:
                try:
                    current = abs(float(row.get("discharge_current_a", "0") or "0"))
                    peak_current = max(peak_current, current)
                    conductivity = float(row.get("channel_conductivity_s_per_m", "0") or "0")
                    max_channel_conductivity = max(max_channel_conductivity, conductivity)
                    charge_error = abs(
                        float(row.get("charge_conservation_error_c", "0") or "0")
                    )
                    max_abs_charge_conservation_error = max(
                        max_abs_charge_conservation_error, charge_error
                    )
                    final_surface_potential = float(row.get("surface_potential_v", "0") or "0")
                    ion = float(row.get("collision_ionization_fraction_step", "0") or "0")
                    exc = float(row.get("collision_excitation_fraction_step", "0") or "0")
                    cx = float(row.get("collision_charge_exchange_fraction_step", "0") or "0")
                    fallback_flag = float(
                        row.get("collision_stage_fallback_triggered", "0") or "0"
                    )
                    collision_fallback_clear = collision_fallback_clear and fallback_flag < 0.5
                    reaction_split_ok = reaction_split_ok or (
                        ion >= 0.0 and exc >= 0.0 and cx >= 0.0 and ion + exc + cx <= 1.0 + 1.0e-9
                    )
                except ValueError:
                    peak_current = math.nan
                    break

            if (
                math.isfinite(peak_current)
                and peak_current > min_peak_current_a
                and math.isfinite(max_channel_conductivity)
                and max_channel_conductivity >= min_channel_conductivity_s_per_m
                and math.isfinite(max_abs_charge_conservation_error)
                and max_abs_charge_conservation_error <= max_charge_conservation_error_c
                and math.isfinite(final_surface_potential)
                and reaction_split_ok
                and (
                    not require_collision_fallback_clear
                    or collision_fallback_clear
                )
            ):
                metrics_status = "PASS"
            else:
                failures.append(
                    f"{preset}:{alignment_mode}:invalid_metrics:"
                    f"peak_current={peak_current}:"
                    f"max_channel_conductivity={max_channel_conductivity}:"
                    f"max_abs_charge_conservation_error={max_abs_charge_conservation_error}:"
                    f"final_surface_potential={final_surface_potential}:"
                    f"reaction_split_ok={reaction_split_ok}:"
                    f"collision_fallback_clear={collision_fallback_clear}"
                )

            sidecar_path = metadata_sidecar(csv_path)
            if sidecar_path.is_file():
                sidecar = load_json(sidecar_path)
                metadata = sidecar.get("metadata", {})
                if (
                    metadata.get("benchmark_contract_family") == "arcpic_core_metrics_v1"
                    and metadata.get("collision_diagnostic_contract_id") == "vacuum-arc-collision-emission-channel-v1"
                ):
                    metadata_status = "PASS"
                else:
                    failures.append(f"{preset}:{alignment_mode}:metadata_contract_missing")

                if (
                    metadata.get("pipeline_contract_id")
                    and metadata.get("collision_stage_enabled") == "true"
                    and metadata.get("collision_reaction_breakdown_enabled") == "true"
                ):
                    contracts_status = "PASS"
                else:
                    failures.append(f"{preset}:{alignment_mode}:pipeline_or_collision_contract_missing")

                benchmark_contract = metadata.get("benchmark_metrics_contract_id")
                benchmark_artifact_name = str(
                    metadata.get("benchmark_metrics_artifact_path", "")
                ).strip()
                benchmark_artifact_path = (
                    csv_path.parent / benchmark_artifact_name
                    if benchmark_artifact_name
                    else csv_path.with_suffix(".benchmark_metrics.json")
                )
                if (
                    benchmark_contract == "vacuum-arc-benchmark-metrics-v1"
                    and benchmark_artifact_path.is_file()
                ):
                    benchmark_payload = load_json(benchmark_artifact_path)
                    reaction_split = benchmark_payload.get("reaction_split", {})
                    if (
                        benchmark_payload.get("schema_version")
                        == "scdat.vacuum_arc.benchmark_metrics.v1"
                        and benchmark_payload.get("contract_id")
                        == "vacuum-arc-benchmark-metrics-v1"
                        and math.isfinite(float(benchmark_payload.get("peak_current_a", "nan")))
                        and math.isfinite(
                            float(
                                benchmark_payload.get(
                                    "max_channel_conductivity_s_per_m", "nan"
                                )
                            )
                        )
                        and isinstance(reaction_split, dict)
                    ):
                        benchmark_artifact_status = "PASS"
                    else:
                        failures.append(
                            f"{preset}:{alignment_mode}:invalid_benchmark_metrics_artifact"
                        )
                else:
                    failures.append(
                        f"{preset}:{alignment_mode}:missing_benchmark_metrics_artifact"
                    )
            else:
                failures.append(f"{preset}:{alignment_mode}:missing_metadata_sidecar")

            if (
                metrics_status == "PASS"
                and metadata_status == "PASS"
                and contracts_status == "PASS"
                and benchmark_artifact_status == "PASS"
            ):
                status = "PASS"

        reports.append(
            {
                "preset": preset,
                "alignment_mode": alignment_mode,
                "command_status": command_status,
                "metrics_status": metrics_status,
                "metadata_status": metadata_status,
                "contracts_status": contracts_status,
                "benchmark_artifact_status": benchmark_artifact_status,
                "status": status,
            }
        )

    report = {
        "status": "PASS" if not failures else "FAIL",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "output_root": str(output_root),
        "matrix": str(matrix_path),
        "cases": reports,
        "failures": failures,
    }
    (output_root / "vacuum_arc_benchmark_matrix_gate.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_markdown(report, output_root / "vacuum_arc_benchmark_matrix_gate.md")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
