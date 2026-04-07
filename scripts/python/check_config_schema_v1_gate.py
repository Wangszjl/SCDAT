#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple


MODULE_CASES = [
    {
        "module": "plasma",
        "definition": "plasmaDocument",
        "command": "plasma-config",
        "config": "scripts/run/plasma_config_full_example.json",
        "default_dt": 5.0e-10,
    },
    {
        "module": "surface",
        "definition": "surfaceDocument",
        "command": "surface-config",
        "config": "scripts/run/surface_config_full_example.json",
        "default_dt": 2.0e-8,
    },
    {
        "module": "radiation",
        "definition": "radiationDocument",
        "command": "radiation-config",
        "config": "scripts/run/radiation_config_full_example.json",
        "default_dt": 5.0e-4,
    },
    {
        "module": "arc",
        "definition": "arcDocument",
        "command": "arc-config",
        "config": "scripts/run/arc_config_full_example.json",
        "default_dt": 5.0e-10,
    },
]


def load_json(path: Path) -> Dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid json at {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"expected object json at {path}")
    return payload


def validate_with_jsonschema(schema: Dict[str, Any], definition_name: str, payload: Dict[str, Any]) -> List[str]:
    defs = schema.get("$defs", {})
    if not isinstance(defs, dict) or definition_name not in defs:
        return [f"missing schema definition: {definition_name}"]

    try:
        import jsonschema
    except ImportError:
        errors: List[str] = []
        if "run" not in payload:
            errors.append("missing top-level key: run")
        if "config" not in payload:
            errors.append("missing top-level key: config")
        return errors

    # Validate via root-aware schema so nested refs like #/$defs/runBlock can resolve.
    target_schema: Dict[str, Any] = {
        "$schema": schema.get("$schema", "https://json-schema.org/draft/2020-12/schema"),
        "$defs": defs,
        "$ref": f"#/$defs/{definition_name}",
    }
    validator = jsonschema.Draft202012Validator(target_schema)
    errors = sorted(validator.iter_errors(payload), key=lambda e: list(e.path))
    messages: List[str] = []
    for error in errors:
        path = "$"
        for token in error.path:
            if isinstance(token, int):
                path += f"[{token}]"
            else:
                path += f".{token}"
        messages.append(f"{path}: {error.message}")
    return messages


def make_smoke_payload(module_case: Dict[str, Any], payload: Dict[str, Any], output_csv: Path) -> Dict[str, Any]:
    smoke = copy.deepcopy(payload)

    smoke["schema_version"] = "v1"
    smoke["module"] = module_case["module"]

    run = smoke.setdefault("run", {})
    if not isinstance(run, dict):
        run = {}
        smoke["run"] = run

    run["steps"] = 1
    run.setdefault("time_step_s", module_case["default_dt"])
    run["output_csv"] = output_csv.as_posix()

    if module_case["module"] == "surface":
        run["adaptive_time_stepping"] = False
        run["total_duration_s"] = 0.0

    smoke["output_csv"] = output_csv.as_posix()
    return smoke


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


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Config Schema v1 Gate Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- schema: {report['schema']}")
    lines.append("")
    lines.append("| module | schema_status | command_status | status |")
    lines.append("|---|---|---|---|")
    for item in report["modules"]:
        lines.append(
            f"| {item['module']} | {item['schema_status']} | {item['command_status']} | {item['status']} |"
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
        description="ENG-001 unified config schema v1 gate: schema validation + command compatibility smoke"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument(
        "--schema",
        default="documentation/contracts/config_schema_v1.json",
        help="Schema path relative to project root unless absolute",
    )
    parser.add_argument("--scdat-exe", required=True)
    parser.add_argument("--output-root", default="build/config_schema_v1_gate")
    parser.add_argument("--report-json", default="")
    parser.add_argument("--report-md", default="")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()

    schema_path = Path(args.schema)
    if not schema_path.is_absolute():
        schema_path = (project_root / schema_path).resolve()

    scdat_exe = Path(args.scdat_exe)
    if not scdat_exe.is_absolute():
        scdat_exe = (project_root / scdat_exe).resolve()

    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = (project_root / output_root).resolve()

    report_json = Path(args.report_json) if args.report_json else output_root / "config_schema_v1_gate.json"
    if not report_json.is_absolute():
        report_json = (project_root / report_json).resolve()

    report_md = Path(args.report_md) if args.report_md else output_root / "config_schema_v1_gate.md"
    if not report_md.is_absolute():
        report_md = (project_root / report_md).resolve()

    try:
        if not schema_path.is_file():
            raise FileNotFoundError(f"missing schema file: {schema_path}")
        if not scdat_exe.is_file():
            raise FileNotFoundError(f"missing SCDAT executable: {scdat_exe}")

        schema = load_json(schema_path)
        output_root.mkdir(parents=True, exist_ok=True)
        generated_root = output_root / "generated_configs"
        generated_root.mkdir(parents=True, exist_ok=True)

        module_reports: List[Dict[str, Any]] = []
        failures: List[str] = []

        for module_case in MODULE_CASES:
            module = module_case["module"]
            config_path = project_root / module_case["config"]
            if not config_path.is_file():
                failures.append(f"{module}:missing_example_config:{config_path}")
                module_reports.append(
                    {
                        "module": module,
                        "schema_status": "FAIL",
                        "command_status": "SKIP",
                        "status": "FAIL",
                        "schema_errors": [f"missing example config: {config_path}"],
                        "command_errors": [],
                    }
                )
                continue

            payload = load_json(config_path)
            schema_errors = validate_with_jsonschema(schema, module_case["definition"], payload)
            schema_status = "PASS" if not schema_errors else "FAIL"
            if schema_errors:
                failures.extend([f"{module}:schema:{msg}" for msg in schema_errors])

            smoke_csv = output_root / f"{module}_config_schema_gate.csv"
            smoke_payload = make_smoke_payload(module_case, payload, smoke_csv)
            smoke_config = generated_root / f"{module}_config_schema_gate.json"
            smoke_config.write_text(json.dumps(smoke_payload, indent=2, ensure_ascii=False), encoding="utf-8")

            cmd = [
                str(scdat_exe),
                module_case["command"],
                str(smoke_config),
                str(smoke_csv),
            ]
            code, stdout_text, stderr_text = run_command(cmd, project_root)

            command_errors: List[str] = []
            if code != 0:
                command_errors.append(f"exit_code={code}")
                if stderr_text.strip():
                    command_errors.append(stderr_text.strip())
            if code == 0 and not smoke_csv.is_file():
                command_errors.append(f"missing_output_csv:{smoke_csv}")

            command_status = "PASS" if not command_errors else "FAIL"
            if command_errors:
                failures.extend([f"{module}:command:{msg}" for msg in command_errors])

            module_status = "PASS" if schema_status == "PASS" and command_status == "PASS" else "FAIL"
            module_reports.append(
                {
                    "module": module,
                    "schema_status": schema_status,
                    "command_status": command_status,
                    "status": module_status,
                    "schema_errors": schema_errors,
                    "command_errors": command_errors,
                    "example_config": str(config_path),
                    "smoke_config": str(smoke_config),
                    "smoke_csv": str(smoke_csv),
                    "stdout": stdout_text.strip(),
                    "stderr": stderr_text.strip(),
                }
            )

        status = "PASS" if not failures else "FAIL"
        report: Dict[str, Any] = {
            "status": status,
            "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "schema": str(schema_path),
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
