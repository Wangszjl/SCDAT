#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


PLASMA_DIAGNOSTICS_SCHEMA = "scdat.plasma.physics_diagnostics.v1"
PLASMA_DIAGNOSTICS_CONTRACT = "plasma-physics-diagnostics-v1"
PLASMA_REACTION_CONTRACT = "plasma-reaction-balance-v1"


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


def metadata_sidecar(path: Path) -> Path:
    return Path(str(path) + ".metadata.json")


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Plasma Physics Diagnostics Gate Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- output_root: {report['output_root']}")
    lines.append("")
    lines.append("| preset | command | metadata | diagnostics | thresholds | status |")
    lines.append("|---|---|---|---|---|---|")
    for item in report["cases"]:
        lines.append(
            f"| {item['preset']} | {item['command_status']} | {item['metadata_status']} | "
            f"{item['diagnostics_status']} | {item['thresholds_status']} | {item['status']} |"
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
        description="Run the Plasma Analysis physics diagnostics contract gate."
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--scdat-exe", required=True)
    parser.add_argument(
        "--matrix",
        default="scripts/run/plasma_physics_diagnostics_matrix.json",
    )
    parser.add_argument("--output-root", default="build/plasma_physics_diagnostics_gate")
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
        raise FileNotFoundError(f"missing plasma diagnostics matrix: {matrix_path}")

    matrix = load_json(matrix_path)
    cases = matrix.get("cases", [])
    if not isinstance(cases, list) or not cases:
        raise ValueError(f"matrix must contain non-empty cases list: {matrix_path}")

    output_root.mkdir(parents=True, exist_ok=True)
    case_reports: List[Dict[str, Any]] = []
    failures: List[str] = []

    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            failures.append(f"case[{index}]:invalid_case_type")
            continue

        preset = str(case.get("preset", "")).strip()
        if not preset:
            failures.append(f"case[{index}]:missing_preset")
            continue

        csv_path = output_root / f"{preset}.csv"
        config_path = output_root / f"{preset}.plasma_gate.json"
        payload = {
            "schema_version": "v1",
            "module": "plasma",
            "base_preset": preset,
            "run": {
                "steps": int(case.get("steps", 2)),
                "time_step_s": float(case.get("time_step_s", 2.5e-8)),
                "output_csv": str(csv_path),
            },
            "config": case.get("config", {}),
        }
        config_path.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

        completed = run_command(
            [str(scdat_exe), "plasma-config", str(config_path), str(csv_path)],
            project_root,
        )

        command_status = "PASS" if completed.returncode == 0 else "FAIL"
        metadata_status = "FAIL"
        diagnostics_status = "FAIL"
        thresholds_status = "FAIL"
        status = "FAIL"

        if completed.returncode != 0:
            failures.append(
                f"{preset}:command_failed:code={completed.returncode}:stderr={completed.stderr.strip()}"
            )
        else:
            sidecar_path = metadata_sidecar(csv_path)
            if not sidecar_path.is_file():
                failures.append(f"{preset}:missing_metadata_sidecar:{sidecar_path}")
            else:
                sidecar = load_json(sidecar_path)
                metadata = sidecar.get("metadata", {})
                if (
                    str(metadata.get("environment_model", "")).strip()
                    == str(case.get("expected_environment_model", "")).strip()
                    and str(metadata.get("reaction_registry", "")).strip()
                    == str(case.get("expected_reaction_registry", "")).strip()
                    and str(metadata.get("diagnostic_set", "")).strip()
                    == str(case.get("expected_diagnostic_set", "")).strip()
                    and str(metadata.get("diagnostic_contract_id", "")).strip()
                    == str(case.get("expected_contract_id", PLASMA_DIAGNOSTICS_CONTRACT)).strip()
                    and str(metadata.get("reaction_contract_id", "")).strip()
                    == PLASMA_REACTION_CONTRACT
                ):
                    metadata_status = "PASS"
                else:
                    failures.append(f"{preset}:metadata_contract_or_mode_mismatch")

                diagnostics_path_text = str(
                    metadata.get("plasma_physics_diagnostics_artifact_path", "")
                ).strip()
                diagnostics_path = (
                    csv_path.parent / diagnostics_path_text
                    if diagnostics_path_text
                    else csv_path.with_suffix(".physics_diagnostics.json")
                )
                if not diagnostics_path.is_file():
                    failures.append(
                        f"{preset}:missing_plasma_physics_diagnostics_artifact:{diagnostics_path}"
                    )
                else:
                    diagnostics = load_json(diagnostics_path)
                    boundary = diagnostics.get("boundary_layer_summary", {})
                    reaction = diagnostics.get("reaction_balance_summary", {})
                    closure = diagnostics.get("advanced_closure_summary", {})
                    flags = diagnostics.get("consistency_flags", {})
                    transport = diagnostics.get("transport_summary", {})

                    if (
                        diagnostics.get("schema_version") == PLASMA_DIAGNOSTICS_SCHEMA
                        and diagnostics.get("contract_id") == PLASMA_DIAGNOSTICS_CONTRACT
                        and diagnostics.get("reaction_contract_id") == PLASMA_REACTION_CONTRACT
                        and isinstance(boundary, dict)
                        and isinstance(reaction, dict)
                        and isinstance(closure, dict)
                        and isinstance(flags, dict)
                        and isinstance(transport, dict)
                    ):
                        diagnostics_status = "PASS"
                    else:
                        failures.append(f"{preset}:diagnostics_schema_or_structure_mismatch")

                    if diagnostics_status == "PASS":
                        dense_ok = (
                            case.get("expected_dense_plasma_detected") is None
                            or bool(diagnostics.get("dense_plasma_detected"))
                            == bool(case.get("expected_dense_plasma_detected"))
                        )
                        advanced_ok = (
                            case.get("require_advanced_closure_enabled") is None
                            or bool(closure.get("advanced_closure_enabled"))
                            == bool(case.get("require_advanced_closure_enabled"))
                        )
                        sheath_ok = float(boundary.get("sheath_thickness_m", 0.0)) >= float(
                            case.get("min_sheath_thickness_m", 0.0)
                        )
                        presheath_ok = float(
                            boundary.get("presheath_thickness_m", 0.0)
                        ) >= float(case.get("min_presheath_thickness_m", 0.0))
                        ion_mach_ok = float(
                            boundary.get("ion_mach_at_sheath_edge", 0.0)
                        ) >= float(case.get("min_ion_mach_at_sheath_edge", 0.0))
                        reaction_process_ok = int(
                            reaction.get("reaction_active_processes", 0)
                        ) >= int(case.get("min_reaction_active_processes", 0))
                        collision_ok = float(
                            reaction.get("effective_collision_frequency_hz", 0.0)
                        ) >= float(case.get("min_effective_collision_frequency_hz", 0.0))
                        non_eq_ok = float(closure.get("non_equilibrium_ratio", 0.0)) >= float(
                            case.get("min_non_equilibrium_ratio", 0.0)
                        )
                        turbulence_ok = float(
                            closure.get("turbulence_eddy_diffusivity_m2_per_s", 0.0)
                        ) >= float(case.get("min_turbulence_eddy_diffusivity_m2_per_s", 0.0))
                        closure_terms_ok = int(
                            closure.get("closure_active_terms", 0)
                        ) >= int(case.get("min_closure_active_terms", 0))
                        boundary_flag_ok = (
                            not bool(case.get("require_boundary_layer_resolved", False))
                            or bool(flags.get("boundary_layer_resolved"))
                        )
                        reaction_flag_ok = (
                            not bool(case.get("require_reaction_balance_active", False))
                            or bool(flags.get("reaction_balance_active"))
                        )
                        transport_ok = float(transport.get("average_density_m3", 0.0)) > 0.0

                        if (
                            dense_ok
                            and advanced_ok
                            and sheath_ok
                            and presheath_ok
                            and ion_mach_ok
                            and reaction_process_ok
                            and collision_ok
                            and non_eq_ok
                            and turbulence_ok
                            and closure_terms_ok
                            and boundary_flag_ok
                            and reaction_flag_ok
                            and transport_ok
                        ):
                            thresholds_status = "PASS"
                            status = "PASS"
                        else:
                            failures.append(
                                f"{preset}:threshold_failure:"
                                f"dense_ok={dense_ok}:advanced_ok={advanced_ok}:"
                                f"sheath_ok={sheath_ok}:presheath_ok={presheath_ok}:"
                                f"ion_mach_ok={ion_mach_ok}:reaction_process_ok={reaction_process_ok}:"
                                f"collision_ok={collision_ok}:non_eq_ok={non_eq_ok}:"
                                f"turbulence_ok={turbulence_ok}:closure_terms_ok={closure_terms_ok}:"
                                f"boundary_flag_ok={boundary_flag_ok}:reaction_flag_ok={reaction_flag_ok}:"
                                f"transport_ok={transport_ok}"
                            )

        case_reports.append(
            {
                "preset": preset,
                "command_status": command_status,
                "metadata_status": metadata_status,
                "diagnostics_status": diagnostics_status,
                "thresholds_status": thresholds_status,
                "status": status,
            }
        )

    report = {
        "status": "PASS" if not failures else "FAIL",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "output_root": str(output_root),
        "matrix": str(matrix_path),
        "cases": case_reports,
        "failures": failures,
    }

    (output_root / "plasma_physics_diagnostics_gate.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_markdown(report, output_root / "plasma_physics_diagnostics_gate.md")

    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
