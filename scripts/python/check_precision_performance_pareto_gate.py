#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

MODULES = ["surface", "vacuum_arc", "radiation", "plasma"]
REQUIRED_TEMPLATE_TOKENS = [
    "{{MODULE_NAME}}",
    "{{GENERATED_UTC}}",
    "{{PHYSICS_GUARDRAIL_STATUS}}",
    "{{FOUR_LAYER_PHYSICS_STATUS}}",
    "{{FOUR_LAYER_STATUS}}",
    "{{DRIFT_STATUS}}",
]


def load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected object json in {path}")
    return payload


def extract_physics_layer(four_layer_report: Dict[str, Any]) -> Dict[str, Any]:
    layers = four_layer_report.get("layers", [])
    if not isinstance(layers, list):
        return {"status": "UNKNOWN", "selected_count": 0}

    for item in layers:
        if isinstance(item, dict) and str(item.get("layer", "")).strip() == "physics":
            return {
                "status": str(item.get("status", "UNKNOWN")).upper(),
                "selected_count": int(item.get("selected_count", 0) or 0),
            }
    return {"status": "UNKNOWN", "selected_count": 0}


def render_template(template_text: str, replacements: Dict[str, str]) -> str:
    output = template_text
    for key, value in replacements.items():
        output = output.replace(key, value)
    return output


def write_markdown_summary(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Precision-Performance Pareto Workflow Gate")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- template: {report['template']}")
    lines.append("")
    lines.append("## Guardrails")
    lines.append("")
    guardrails = report["guardrails"]
    lines.append(f"- four_layer_status: {guardrails['four_layer_status']}")
    lines.append(f"- four_layer_physics_status: {guardrails['four_layer_physics_status']}")
    lines.append(f"- four_layer_physics_selected_count: {guardrails['four_layer_physics_selected_count']}")
    lines.append(f"- drift_dashboard_status: {guardrails['drift_dashboard_status']}")
    lines.append("")
    lines.append("## Module Reports")
    lines.append("")
    lines.append("| module | output_md | status |")
    lines.append("|---|---|---|")
    for item in report["modules"]:
        lines.append(f"| {item['module']} | {item['output_md']} | {item['status']} |")

    if report.get("warnings"):
        lines.append("")
        lines.append("## Warnings")
        lines.append("")
        for warning in report["warnings"]:
            lines.append(f"- {warning}")

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
        description="ENG-008 precision-performance pareto reporting workflow gate"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument(
        "--template",
        default="documentation/templates/precision_performance_pareto_template_v1.md",
    )
    parser.add_argument(
        "--four-layer-report",
        default="build/four_layer_gate_report.json",
    )
    parser.add_argument(
        "--drift-report",
        default="build/drift_monitor_dashboard.json",
    )
    parser.add_argument(
        "--output-root",
        default="build/pareto_reports",
    )
    parser.add_argument(
        "--report-json",
        default="build/pareto_report_workflow_gate.json",
    )
    parser.add_argument(
        "--report-md",
        default="build/pareto_report_workflow_gate.md",
    )
    parser.add_argument(
        "--require-four-layer-overall-pass",
        action="store_true",
        help="Fail when four_layer_report.status is not PASS",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()

    template_path = Path(args.template)
    if not template_path.is_absolute():
        template_path = (project_root / template_path).resolve()

    four_layer_path = Path(args.four_layer_report)
    if not four_layer_path.is_absolute():
        four_layer_path = (project_root / four_layer_path).resolve()

    drift_path = Path(args.drift_report)
    if not drift_path.is_absolute():
        drift_path = (project_root / drift_path).resolve()

    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = (project_root / output_root).resolve()

    report_json = Path(args.report_json)
    if not report_json.is_absolute():
        report_json = (project_root / report_json).resolve()

    report_md = Path(args.report_md)
    if not report_md.is_absolute():
        report_md = (project_root / report_md).resolve()

    failures: List[str] = []
    warnings: List[str] = []
    module_reports: List[Dict[str, str]] = []

    try:
        if not template_path.is_file():
            raise FileNotFoundError(f"missing template: {template_path}")
        if not four_layer_path.is_file():
            raise FileNotFoundError(f"missing four-layer report: {four_layer_path}")
        if not drift_path.is_file():
            raise FileNotFoundError(f"missing drift report: {drift_path}")

        template_text = template_path.read_text(encoding="utf-8")
        for token in REQUIRED_TEMPLATE_TOKENS:
            if token not in template_text:
                failures.append(f"missing_template_token:{token}")

        four_layer_report = load_json(four_layer_path)
        drift_report = load_json(drift_path)

        four_layer_status = str(four_layer_report.get("status", "UNKNOWN")).upper()
        physics_layer = extract_physics_layer(four_layer_report)
        physics_status = str(physics_layer.get("status", "UNKNOWN")).upper()
        physics_selected_count = int(physics_layer.get("selected_count", 0) or 0)
        drift_status = str(drift_report.get("status", "UNKNOWN")).upper()

        if four_layer_status != "PASS":
            message = f"four_layer_status_not_pass:{four_layer_status}"
            if args.require_four_layer_overall_pass:
                failures.append(message)
            else:
                warnings.append(message)
        if physics_status != "PASS":
            failures.append(f"four_layer_physics_status_not_pass:{physics_status}")
        if physics_selected_count <= 0:
            failures.append("four_layer_physics_selected_count_not_positive")
        if drift_status != "PASS":
            failures.append(f"drift_status_not_pass:{drift_status}")

        generated_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        output_root.mkdir(parents=True, exist_ok=True)

        physics_guardrail_status = "PASS" if not failures else "FAIL"
        replacements_common = {
            "{{GENERATED_UTC}}": generated_utc,
            "{{PHYSICS_GUARDRAIL_STATUS}}": physics_guardrail_status,
            "{{FOUR_LAYER_PHYSICS_STATUS}}": physics_status,
            "{{FOUR_LAYER_STATUS}}": four_layer_status,
            "{{DRIFT_STATUS}}": drift_status,
        }

        for module in MODULES:
            replacements = dict(replacements_common)
            replacements["{{MODULE_NAME}}"] = module
            rendered = render_template(template_text, replacements)

            output_file = output_root / f"{module}_pareto_report_v1.md"
            output_file.write_text(rendered, encoding="utf-8")

            module_reports.append(
                {
                    "module": module,
                    "output_md": str(output_file),
                    "status": "PASS" if output_file.is_file() else "FAIL",
                }
            )

        for item in module_reports:
            if item["status"] != "PASS":
                failures.append(f"module_report_not_generated:{item['module']}")

        status = "PASS" if not failures else "FAIL"
        report = {
            "status": status,
            "timestamp_utc": generated_utc,
            "project_root": str(project_root),
            "template": str(template_path),
            "inputs": {
                "four_layer_report": str(four_layer_path),
                "drift_report": str(drift_path),
            },
            "guardrails": {
                "four_layer_status": four_layer_status,
                "four_layer_physics_status": physics_status,
                "four_layer_physics_selected_count": physics_selected_count,
                "drift_dashboard_status": drift_status,
            },
            "modules": module_reports,
            "warnings": warnings,
            "failures": failures,
        }

        report_json.parent.mkdir(parents=True, exist_ok=True)
        report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        write_markdown_summary(report, report_md)

        print(f"report_json={report_json}")
        print(f"report_md={report_md}")
        print(f"status={status}")

        return 0 if status == "PASS" else 2

    except Exception as exc:  # pylint: disable=broad-except
        print(f"status=FAIL({exc})")
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
