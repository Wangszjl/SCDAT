#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
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


def metadata_sidecar(path: Path) -> Path:
    return Path(str(path) + ".metadata.json")


INTERNAL_DRIVE_PROVENANCE_SCHEMA = "scdat.internal_radiation_drive_provenance.v1"
INTERNAL_DRIVE_PROVENANCE_CONTRACT = "internal-radiation-drive-provenance-v1"
INTERNAL_DRIVE_LAYER_ALIGNMENT_CONTRACT = "internal-radiation-layer-alignment-v1"
INTERNAL_RESPONSE_COUPLED_SOLVE_SCHEMA = "scdat.internal_response_coupled_solve.v1"
INTERNAL_RESPONSE_COUPLED_SOLVE_CONTRACT = "internal-response-coupled-solve-v1"
RADIATION_TRANSPORT_BENCHMARK_SCHEMA = "scdat.radiation_transport_benchmark.v1"
RADIATION_TRANSPORT_BENCHMARK_CONTRACT = "radiation-transport-benchmark-v1"


def load_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def resolve_relative_artifact(path_text: str, project_root: Path, output_root: Path) -> Path:
    candidate = Path(path_text)
    if candidate.is_absolute():
        return candidate
    if candidate.parts[:1] == ("build_codex",):
        return (project_root / candidate).resolve()
    return (output_root / candidate.name).resolve()


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Internal Radiation Coupling Contract Gate Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- output_root: {report['output_root']}")
    lines.append("")
    lines.append("| internal | radiation | command | radiation_artifacts | internal_contract | coupling_history | status |")
    lines.append("|---|---|---|---|---|---|---|")
    for item in report["cases"]:
        lines.append(
            f"| {item['internal_preset']} | {item['radiation_preset']} | {item['command_status']} | "
            f"{item['radiation_artifacts_status']} | {item['internal_contract_status']} | "
            f"{item['coupling_history_status']} | {item['status']} |"
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
        description="Run the Internal-Radiation coupling contract gate."
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--scdat-exe", required=True)
    parser.add_argument(
        "--matrix",
        default="scripts/run/internal_radiation_coupling_matrix.json",
    )
    parser.add_argument("--output-root", default="build/internal_radiation_coupling_contract_gate")
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

        internal_preset = str(case.get("internal_preset", "")).strip()
        radiation_preset = str(case.get("radiation_preset", "")).strip()
        expected_deposition_contract_id = str(
            case.get("expected_deposition_contract_id", "")
        ).strip()
        expected_process_contract_id = str(
            case.get("expected_process_contract_id", "")
        ).strip()
        expected_internal_coupled_solve_contract_id = str(
            case.get(
                "expected_internal_response_coupled_solve_contract_id",
                INTERNAL_RESPONSE_COUPLED_SOLVE_CONTRACT,
            )
        ).strip()
        expected_radiation_transport_benchmark_contract_id = str(
            case.get(
                "expected_radiation_transport_benchmark_contract_id",
                RADIATION_TRANSPORT_BENCHMARK_CONTRACT,
            )
        ).strip()
        if not internal_preset or not radiation_preset:
            failures.append(f"case[{index}]:missing_presets")
            continue

        output_csv = output_root / f"{internal_preset}__{radiation_preset}.csv"
        command = [
            str(scdat_exe),
            "internal-radiation",
            internal_preset,
            radiation_preset,
            str(output_csv),
        ]
        completed = run_command(command, project_root)

        command_status = "PASS" if completed.returncode == 0 else "FAIL"
        radiation_artifacts_status = "FAIL"
        internal_contract_status = "FAIL"
        coupling_history_status = "FAIL"
        status = "FAIL"

        if completed.returncode != 0:
            failures.append(
                f"{internal_preset}:{radiation_preset}:command_failed:code={completed.returncode}:stderr={completed.stderr.strip()}"
            )
        else:
            internal_sidecar_path = metadata_sidecar(output_csv)
            coupling_history_path = output_csv.with_suffix(".coupling_history.csv")
            radiation_drive_csv = output_csv.with_suffix(".radiation_drive.csv")
            radiation_sidecar_path = metadata_sidecar(radiation_drive_csv)

            if coupling_history_path.is_file() and load_csv_rows(coupling_history_path):
                coupling_history_status = "PASS"
            else:
                failures.append(
                    f"{internal_preset}:{radiation_preset}:missing_or_empty_coupling_history:{coupling_history_path}"
                )

            if not internal_sidecar_path.is_file():
                failures.append(
                    f"{internal_preset}:{radiation_preset}:missing_internal_sidecar:{internal_sidecar_path}"
                )
            else:
                internal_sidecar = load_json(internal_sidecar_path)
                internal_metadata = internal_sidecar.get("metadata", {})
                deposition_contract = str(
                    internal_metadata.get("deposition_record_contract_id", "")
                ).strip()
                process_contract = str(
                    internal_metadata.get("process_history_contract_id", "")
                ).strip()
                consumes_history = str(
                    internal_metadata.get("internal_consumes_deposition_history", "")
                ).strip()
                deposition_history_path = resolve_relative_artifact(
                    str(internal_metadata.get("deposition_history_path", "")),
                    project_root,
                    output_root,
                )
                process_history_path = resolve_relative_artifact(
                    str(internal_metadata.get("process_history_path", "")),
                    project_root,
                    output_root,
                )
                provenance_summary_path = resolve_relative_artifact(
                    str(internal_metadata.get("internal_drive_provenance_artifact_path", "")),
                    project_root,
                    output_root,
                )
                layer_alignment_path = resolve_relative_artifact(
                    str(internal_metadata.get("internal_drive_layer_alignment_artifact_path", "")),
                    project_root,
                    output_root,
                )
                provenance_contract = str(
                    internal_metadata.get("internal_drive_provenance_contract_id", "")
                ).strip()
                layer_alignment_contract = str(
                    internal_metadata.get("internal_drive_layer_alignment_contract_id", "")
                ).strip()
                coupled_solve_contract = str(
                    internal_metadata.get("internal_response_coupled_solve_contract_id", "")
                ).strip()
                coupled_solve_path = resolve_relative_artifact(
                    str(
                        internal_metadata.get(
                            "internal_response_coupled_solve_artifact_path", ""
                        )
                    ),
                    project_root,
                    output_root,
                )

                if (
                    deposition_contract == expected_deposition_contract_id
                    and process_contract == expected_process_contract_id
                    and consumes_history == "true"
                    and deposition_history_path.is_file()
                    and process_history_path.is_file()
                    and provenance_contract == INTERNAL_DRIVE_PROVENANCE_CONTRACT
                    and provenance_summary_path.is_file()
                    and layer_alignment_contract == INTERNAL_DRIVE_LAYER_ALIGNMENT_CONTRACT
                    and layer_alignment_path.is_file()
                    and coupled_solve_contract == expected_internal_coupled_solve_contract_id
                    and coupled_solve_path.is_file()
                ):
                    provenance_payload = load_json(provenance_summary_path)
                    alignment_rows = load_csv_rows(layer_alignment_path)
                    coupled_solve_payload = load_json(coupled_solve_path)
                    if (
                        provenance_payload.get("schema") == INTERNAL_DRIVE_PROVENANCE_SCHEMA
                        and provenance_payload.get("contract_id")
                        == INTERNAL_DRIVE_PROVENANCE_CONTRACT
                        and alignment_rows
                        and coupled_solve_payload.get("schema_version")
                        == INTERNAL_RESPONSE_COUPLED_SOLVE_SCHEMA
                        and coupled_solve_payload.get("contract_id")
                        == expected_internal_coupled_solve_contract_id
                    ):
                        internal_contract_status = "PASS"
                    else:
                        failures.append(
                            f"{internal_preset}:{radiation_preset}:internal_provenance_schema_or_contract_mismatch"
                        )
                else:
                    failures.append(
                        f"{internal_preset}:{radiation_preset}:internal_contract_mismatch:"
                        f"deposition_contract={deposition_contract}:process_contract={process_contract}:"
                        f"consumes_history={consumes_history}:deposition_path={deposition_history_path}:"
                        f"process_path={process_history_path}:provenance_contract={provenance_contract}:"
                        f"provenance_path={provenance_summary_path}:layer_alignment_contract={layer_alignment_contract}:"
                        f"layer_alignment_path={layer_alignment_path}:"
                        f"coupled_solve_contract={coupled_solve_contract}:"
                        f"coupled_solve_path={coupled_solve_path}"
                    )

            if not radiation_drive_csv.is_file() or not radiation_sidecar_path.is_file():
                failures.append(
                    f"{internal_preset}:{radiation_preset}:missing_radiation_drive_artifacts:{radiation_drive_csv}"
                )
            else:
                radiation_sidecar = load_json(radiation_sidecar_path)
                radiation_metadata = radiation_sidecar.get("metadata", {})
                deposition_history_path = resolve_relative_artifact(
                    str(radiation_metadata.get("deposition_history_artifact_path", "")),
                    project_root,
                    output_root,
                )
                process_history_path = resolve_relative_artifact(
                    str(radiation_metadata.get("process_history_artifact_path", "")),
                    project_root,
                    output_root,
                )
                radiation_transport_benchmark_contract = str(
                    radiation_metadata.get("radiation_transport_benchmark_contract_id", "")
                ).strip()
                radiation_transport_benchmark_path = resolve_relative_artifact(
                    str(
                        radiation_metadata.get(
                            "radiation_transport_benchmark_artifact_path", ""
                        )
                    ),
                    project_root,
                    output_root,
                )

                artifacts_ok = (
                    str(radiation_metadata.get("deposition_record_contract_id", "")).strip()
                    == expected_deposition_contract_id
                    and str(radiation_metadata.get("process_history_contract_id", "")).strip()
                    == expected_process_contract_id
                    and deposition_history_path.is_file()
                    and process_history_path.is_file()
                    and radiation_transport_benchmark_contract
                    == expected_radiation_transport_benchmark_contract_id
                    and radiation_transport_benchmark_path.is_file()
                )
                if not artifacts_ok:
                    failures.append(
                        f"{internal_preset}:{radiation_preset}:radiation_artifact_contract_mismatch"
                    )
                else:
                    deposition_payload = load_json(deposition_history_path)
                    process_payload = load_json(process_history_path)
                    benchmark_payload = load_json(radiation_transport_benchmark_path)
                    if (
                        deposition_payload.get("schema") == "scdat.radiation.deposition_history.v1"
                        and deposition_payload.get("contract_id") == expected_deposition_contract_id
                        and process_payload.get("schema") == "scdat.radiation.process_history.v1"
                        and process_payload.get("contract_id") == expected_process_contract_id
                        and benchmark_payload.get("schema_version")
                        == RADIATION_TRANSPORT_BENCHMARK_SCHEMA
                        and benchmark_payload.get("contract_id")
                        == expected_radiation_transport_benchmark_contract_id
                    ):
                        radiation_artifacts_status = "PASS"
                    else:
                        failures.append(
                            f"{internal_preset}:{radiation_preset}:radiation_history_schema_or_contract_mismatch"
                        )

            if (
                command_status == "PASS"
                and radiation_artifacts_status == "PASS"
                and internal_contract_status == "PASS"
                and coupling_history_status == "PASS"
            ):
                status = "PASS"

        reports.append(
            {
                "internal_preset": internal_preset,
                "radiation_preset": radiation_preset,
                "command_status": command_status,
                "radiation_artifacts_status": radiation_artifacts_status,
                "internal_contract_status": internal_contract_status,
                "coupling_history_status": coupling_history_status,
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
    (output_root / "internal_radiation_coupling_contract_gate.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_markdown(report, output_root / "internal_radiation_coupling_contract_gate.md")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
