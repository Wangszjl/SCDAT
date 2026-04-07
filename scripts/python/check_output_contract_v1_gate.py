#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple


MODULE_CASES = [
    {
        "module": "plasma",
        "command": "plasma-config",
        "config": "scripts/run/plasma_config_full_example.json",
        "default_dt": 5.0e-10,
    },
    {
        "module": "surface",
        "command": "surface-config",
        "config": "scripts/run/surface_config_full_example.json",
        "default_dt": 2.0e-8,
    },
    {
        "module": "radiation",
        "command": "radiation-config",
        "config": "scripts/run/radiation_config_full_example.json",
        "default_dt": 5.0e-4,
    },
    {
        "module": "arc",
        "command": "arc-config",
        "config": "scripts/run/arc_config_full_example.json",
        "default_dt": 5.0e-10,
    },
]

REQUIRED_TOP_LEVEL_KEYS = {
    "output_contract_schema",
    "output_contract_version",
    "artifact_format",
    "artifact_path",
    "metadata",
}

REQUIRED_METADATA_KEYS = {
    "module",
    "output_contract_schema",
    "output_contract_version",
    "output_axis_name",
    "output_row_count",
    "output_series_count",
    "output_series_names",
    "output_artifact_csv",
    "output_artifact_vtk",
    "output_artifact_h5",
    "output_export_utc",
}

EXPECTED_FORMAT_SUFFIX = {
    "csv": ".csv",
    "vtk": ".vtk",
    "h5": ".h5",
}


def load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected top-level object in {path}")
    return payload


def run_command(command: List[str], cwd: Path) -> Tuple[int, str, str]:
    completed = subprocess.run(
        command,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    return completed.returncode, completed.stdout, completed.stderr


def make_smoke_payload(module_case: Dict[str, Any], payload: Dict[str, Any], output_csv: Path) -> Dict[str, Any]:
    smoke = copy.deepcopy(payload)

    smoke["schema_version"] = "v1"
    smoke["module"] = module_case["module"]

    run_block = smoke.setdefault("run", {})
    if not isinstance(run_block, dict):
        run_block = {}
        smoke["run"] = run_block

    run_block["steps"] = 1
    run_block.setdefault("time_step_s", module_case["default_dt"])
    run_block["output_csv"] = output_csv.as_posix()

    if module_case["module"] == "surface":
        run_block["adaptive_time_stepping"] = False
        run_block["total_duration_s"] = 0.0

    smoke["output_csv"] = output_csv.as_posix()
    return smoke


def metadata_sidecar_path(artifact_path: Path) -> Path:
    return Path(str(artifact_path) + ".metadata.json")


def parse_positive_int(text: str, key: str, errors: List[str]) -> int:
    try:
        value = int(text)
    except (TypeError, ValueError):
        errors.append(f"invalid_integer:{key}={text}")
        return -1

    if value <= 0:
        errors.append(f"non_positive_integer:{key}={value}")
    return value


def normalize_path(text: str) -> Path:
    return Path(text).resolve()


def validate_sidecar(
    artifact_format: str,
    artifact_path: Path,
    expected_artifacts: Dict[str, Path],
) -> Tuple[Dict[str, Any], List[str]]:
    errors: List[str] = []
    sidecar_path = metadata_sidecar_path(artifact_path)

    if not sidecar_path.is_file():
        return {}, [f"missing_sidecar:{sidecar_path}"]

    try:
        payload = load_json(sidecar_path)
    except Exception as exc:  # pylint: disable=broad-except
        return {}, [f"invalid_sidecar_json:{sidecar_path}:{exc}"]

    missing_top_level = sorted(REQUIRED_TOP_LEVEL_KEYS - payload.keys())
    for key in missing_top_level:
        errors.append(f"missing_top_level_key:{sidecar_path}:{key}")

    if payload.get("artifact_format") != artifact_format:
        errors.append(
            f"artifact_format_mismatch:{sidecar_path}:expected={artifact_format}:actual={payload.get('artifact_format')}"
        )

    if payload.get("output_contract_schema") != "scdat.output.contract.v1":
        errors.append(
            f"invalid_contract_schema:{sidecar_path}:{payload.get('output_contract_schema')}"
        )

    if payload.get("output_contract_version") != "v1":
        errors.append(
            f"invalid_contract_version:{sidecar_path}:{payload.get('output_contract_version')}"
        )

    artifact_path_text = payload.get("artifact_path", "")
    try:
        actual_artifact_path = normalize_path(str(artifact_path_text))
    except Exception as exc:  # pylint: disable=broad-except
        errors.append(f"invalid_artifact_path:{sidecar_path}:{exc}")
        actual_artifact_path = Path("<invalid>")

    if actual_artifact_path != artifact_path.resolve():
        errors.append(
            f"artifact_path_mismatch:{sidecar_path}:expected={artifact_path.resolve()}:actual={actual_artifact_path}"
        )

    metadata = payload.get("metadata")
    if not isinstance(metadata, dict):
        errors.append(f"missing_or_invalid_metadata_block:{sidecar_path}")
        return {}, errors

    missing_metadata = sorted(REQUIRED_METADATA_KEYS - metadata.keys())
    for key in missing_metadata:
        errors.append(f"missing_metadata_key:{sidecar_path}:{key}")

    module_label = str(metadata.get("module", "")).strip()
    if not module_label:
        errors.append(f"empty_module_label:{sidecar_path}")

    if metadata.get("output_contract_schema") != "scdat.output.contract.v1":
        errors.append(
            f"metadata_contract_schema_mismatch:{sidecar_path}:{metadata.get('output_contract_schema')}"
        )

    if metadata.get("output_contract_version") != "v1":
        errors.append(
            f"metadata_contract_version_mismatch:{sidecar_path}:{metadata.get('output_contract_version')}"
        )

    axis_name = str(metadata.get("output_axis_name", "")).strip()
    if not axis_name:
        errors.append(f"empty_output_axis_name:{sidecar_path}")

    parse_positive_int(str(metadata.get("output_row_count", "")), "output_row_count", errors)
    parse_positive_int(str(metadata.get("output_series_count", "")), "output_series_count", errors)

    series_names = str(metadata.get("output_series_names", "")).strip()
    if not series_names:
        errors.append(f"empty_output_series_names:{sidecar_path}")

    for key, expected_path in (
        ("output_artifact_csv", expected_artifacts["csv"]),
        ("output_artifact_vtk", expected_artifacts["vtk"]),
        ("output_artifact_h5", expected_artifacts["h5"]),
    ):
        value = metadata.get(key)
        if value is None:
            continue
        try:
            normalized = normalize_path(str(value))
        except Exception as exc:  # pylint: disable=broad-except
            errors.append(f"invalid_metadata_path:{sidecar_path}:{key}:{exc}")
            continue
        if normalized != expected_path.resolve():
            errors.append(
                f"metadata_path_mismatch:{sidecar_path}:{key}:expected={expected_path.resolve()}:actual={normalized}"
            )

    return metadata, errors


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Output Contract v1 Gate Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- output_root: {report['output_root']}")
    lines.append("")
    lines.append("| module | command_status | artifacts_status | sidecar_status | status |")
    lines.append("|---|---|---|---|---|")
    for item in report["modules"]:
        lines.append(
            "| "
            + item["module"]
            + " | "
            + item["command_status"]
            + " | "
            + item["artifacts_status"]
            + " | "
            + item["sidecar_status"]
            + " | "
            + item["status"]
            + " |"
        )

    failures = report.get("failures", [])
    if failures:
        lines.append("")
        lines.append("## Failures")
        lines.append("")
        for failure in failures:
            lines.append(f"- {failure}")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="ENG-004 unified output contract v1 gate: verify csv/vtk/h5 sidecar metadata consistency"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--scdat-exe", required=True)
    parser.add_argument("--output-root", default="build/output_contract_v1_gate")
    parser.add_argument("--report-json", default="")
    parser.add_argument("--report-md", default="")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()

    scdat_exe = Path(args.scdat_exe)
    if not scdat_exe.is_absolute():
        scdat_exe = (project_root / scdat_exe).resolve()

    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = (project_root / output_root).resolve()

    report_json = Path(args.report_json) if args.report_json else output_root / "output_contract_v1_gate.json"
    if not report_json.is_absolute():
        report_json = (project_root / report_json).resolve()

    report_md = Path(args.report_md) if args.report_md else output_root / "output_contract_v1_gate.md"
    if not report_md.is_absolute():
        report_md = (project_root / report_md).resolve()

    try:
        if not scdat_exe.is_file():
            raise FileNotFoundError(f"missing SCDAT executable: {scdat_exe}")

        output_root.mkdir(parents=True, exist_ok=True)
        generated_root = output_root / "generated_configs"
        generated_root.mkdir(parents=True, exist_ok=True)

        module_reports: List[Dict[str, Any]] = []
        failures: List[str] = []

        for module_case in MODULE_CASES:
            module = module_case["module"]
            config_path = (project_root / module_case["config"]).resolve()
            module_failures: List[str] = []

            if not config_path.is_file():
                module_failures.append(f"missing_example_config:{config_path}")
                module_reports.append(
                    {
                        "module": module,
                        "command_status": "SKIP",
                        "artifacts_status": "SKIP",
                        "sidecar_status": "SKIP",
                        "status": "FAIL",
                        "errors": module_failures,
                    }
                )
                failures.extend([f"{module}:{item}" for item in module_failures])
                continue

            payload = load_json(config_path)
            smoke_csv = output_root / f"{module}_output_contract_gate.csv"
            smoke_payload = make_smoke_payload(module_case, payload, smoke_csv)
            smoke_config = generated_root / f"{module}_output_contract_gate.json"
            smoke_config.write_text(json.dumps(smoke_payload, indent=2, ensure_ascii=False), encoding="utf-8")

            command = [
                str(scdat_exe),
                module_case["command"],
                str(smoke_config),
                str(smoke_csv),
            ]
            return_code, stdout_text, stderr_text = run_command(command, project_root)

            command_status = "PASS"
            if return_code != 0:
                command_status = "FAIL"
                module_failures.append(f"command_exit_code={return_code}")
                if stderr_text.strip():
                    module_failures.append(f"command_stderr={stderr_text.strip()}")

            stem = smoke_csv.with_suffix("")
            expected_artifacts = {
                "csv": smoke_csv.resolve(),
                "vtk": Path(str(stem) + EXPECTED_FORMAT_SUFFIX["vtk"]).resolve(),
                "h5": Path(str(stem) + EXPECTED_FORMAT_SUFFIX["h5"]).resolve(),
            }

            artifacts_status = "PASS"
            for artifact in expected_artifacts.values():
                if not artifact.is_file():
                    artifacts_status = "FAIL"
                    module_failures.append(f"missing_artifact:{artifact}")

            sidecar_status = "PASS"
            metadata_by_format: Dict[str, Dict[str, Any]] = {}
            for artifact_format, artifact_path in expected_artifacts.items():
                metadata, sidecar_errors = validate_sidecar(
                    artifact_format=artifact_format,
                    artifact_path=artifact_path,
                    expected_artifacts=expected_artifacts,
                )
                if sidecar_errors:
                    sidecar_status = "FAIL"
                    module_failures.extend(sidecar_errors)
                if metadata:
                    metadata_by_format[artifact_format] = metadata

            if {"csv", "vtk", "h5"}.issubset(metadata_by_format.keys()):
                reference = metadata_by_format["csv"]
                for peer_format in ("vtk", "h5"):
                    peer = metadata_by_format[peer_format]
                    for key in sorted(REQUIRED_METADATA_KEYS):
                        if reference.get(key) != peer.get(key):
                            sidecar_status = "FAIL"
                            module_failures.append(
                                f"cross_format_metadata_mismatch:{module}:{key}:csv={reference.get(key)}:{peer_format}={peer.get(key)}"
                            )

            module_status = "PASS"
            if command_status != "PASS" or artifacts_status != "PASS" or sidecar_status != "PASS":
                module_status = "FAIL"

            module_reports.append(
                {
                    "module": module,
                    "command": command,
                    "example_config": str(config_path),
                    "smoke_config": str(smoke_config),
                    "stdout": stdout_text.strip(),
                    "stderr": stderr_text.strip(),
                    "command_status": command_status,
                    "artifacts_status": artifacts_status,
                    "sidecar_status": sidecar_status,
                    "status": module_status,
                    "errors": module_failures,
                }
            )

            if module_failures:
                failures.extend([f"{module}:{item}" for item in module_failures])

        status = "PASS" if not failures else "FAIL"
        report: Dict[str, Any] = {
            "status": status,
            "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "scdat_exe": str(scdat_exe),
            "output_root": str(output_root),
            "modules": module_reports,
            "failures": failures,
        }

        report_json.parent.mkdir(parents=True, exist_ok=True)
        report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        write_markdown(report, report_md)

        print(f"report_json={report_json}")
        print(f"report_md={report_md}")
        print(f"status={status}")
        return 0 if status == "PASS" else 2

    except Exception as exc:  # pylint: disable=broad-except
        print(f"status=FAIL({exc})")
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
