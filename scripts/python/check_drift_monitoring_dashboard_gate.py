#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple


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


def load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected json object in {path}")
    return payload


def run_required_gate(test_name: str, project_root: Path, build_dir: Path, config: str) -> None:
    command = [
        "ctest",
        "--test-dir",
        str(build_dir),
        "-C",
        config,
        "-R",
        f"^{re.escape(test_name)}$",
        "--output-on-failure",
    ]
    completed = run_command(command, project_root)
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or "unknown error"
        raise RuntimeError(f"required gate {test_name} failed: {message}")


def ensure_report(
    report_path: Path,
    producer_test_name: str,
    project_root: Path,
    build_dir: Path,
    config: str,
    auto_run_missing: bool,
) -> None:
    if report_path.is_file():
        return
    if not auto_run_missing:
        raise FileNotFoundError(f"missing required report: {report_path}")

    run_required_gate(producer_test_name, project_root, build_dir, config)
    if not report_path.is_file():
        raise FileNotFoundError(
            f"report still missing after running {producer_test_name}: {report_path}"
        )


def extract_schema_metrics(payload: Dict[str, Any]) -> Dict[str, Any]:
    modules = payload.get("modules", [])
    if not isinstance(modules, list):
        modules = []

    total = len(modules)
    passed = 0
    for item in modules:
        if isinstance(item, dict) and str(item.get("status", "")).upper() == "PASS":
            passed += 1

    pass_rate = float(passed) / float(total) if total > 0 else 0.0
    return {
        "status": str(payload.get("status", "UNKNOWN")).upper(),
        "module_total": total,
        "module_passed": passed,
        "module_pass_rate": pass_rate,
    }


def extract_output_contract_metrics(payload: Dict[str, Any]) -> Dict[str, Any]:
    modules = payload.get("modules", [])
    if not isinstance(modules, list):
        modules = []

    total = len(modules)
    passed = 0
    for item in modules:
        if isinstance(item, dict) and str(item.get("status", "")).upper() == "PASS":
            passed += 1

    pass_rate = float(passed) / float(total) if total > 0 else 0.0
    return {
        "status": str(payload.get("status", "UNKNOWN")).upper(),
        "module_total": total,
        "module_passed": passed,
        "module_pass_rate": pass_rate,
    }


def extract_unit_lint_metrics(payload: Dict[str, Any]) -> Dict[str, Any]:
    errors = payload.get("errors", [])
    warnings = payload.get("warnings", [])
    if not isinstance(errors, list):
        errors = []
    if not isinstance(warnings, list):
        warnings = []

    return {
        "status": str(payload.get("status", "UNKNOWN")).upper(),
        "error_count": len(errors),
        "warning_count": len(warnings),
        "checked_keys": int(payload.get("checked_keys", 0) or 0),
    }


def extract_four_layer_metrics(payload: Dict[str, Any]) -> Dict[str, Any]:
    layers_raw = payload.get("layers", [])
    if not isinstance(layers_raw, list):
        layers_raw = []

    metrics: Dict[str, Any] = {
        "status": str(payload.get("status", "UNKNOWN")).upper(),
        "layer_status": {},
        "layer_selected_count": {},
    }

    for item in layers_raw:
        if not isinstance(item, dict):
            continue
        layer = str(item.get("layer", "")).strip()
        if not layer:
            continue
        metrics["layer_status"][layer] = str(item.get("status", "UNKNOWN")).upper()
        metrics["layer_selected_count"][layer] = int(item.get("selected_count", 0) or 0)

    return metrics


def load_history(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        return {"schema_version": "v1", "snapshots": []}

    payload = load_json(path)
    snapshots = payload.get("snapshots", [])
    if not isinstance(snapshots, list):
        snapshots = []

    return {
        "schema_version": str(payload.get("schema_version", "v1")),
        "snapshots": snapshots,
    }


def add_alert(alerts: List[Dict[str, str]], severity: str, code: str, message: str) -> None:
    alerts.append({"severity": severity, "code": code, "message": message})


def evaluate_alerts(
    snapshot: Dict[str, Any],
    previous_snapshot: Dict[str, Any] | None,
    min_layer_counts: Dict[str, int],
    max_layer_count_drop_ratio: float,
) -> List[Dict[str, str]]:
    alerts: List[Dict[str, str]] = []

    schema = snapshot["metrics"]["schema"]
    output_contract = snapshot["metrics"]["output_contract"]
    unit_lint = snapshot["metrics"]["unit_lint"]
    four_layer = snapshot["metrics"]["four_layer"]

    if schema["status"] != "PASS":
        add_alert(alerts, "critical", "schema_status", f"schema gate status={schema['status']}")
    if schema["module_pass_rate"] < 1.0:
        add_alert(
            alerts,
            "critical",
            "schema_module_pass_rate",
            f"schema module pass rate={schema['module_pass_rate']:.3f}",
        )

    if output_contract["status"] != "PASS":
        add_alert(
            alerts,
            "critical",
            "output_contract_status",
            f"output contract gate status={output_contract['status']}",
        )
    if output_contract["module_pass_rate"] < 1.0:
        add_alert(
            alerts,
            "critical",
            "output_contract_module_pass_rate",
            f"output contract module pass rate={output_contract['module_pass_rate']:.3f}",
        )

    if unit_lint["status"] != "PASS":
        add_alert(alerts, "critical", "unit_lint_status", f"unit lint status={unit_lint['status']}")
    if unit_lint["error_count"] > 0:
        add_alert(
            alerts,
            "critical",
            "unit_lint_errors",
            f"unit lint errors={unit_lint['error_count']}",
        )
    if unit_lint["warning_count"] > 0:
        add_alert(
            alerts,
            "warning",
            "unit_lint_warnings",
            f"unit lint warnings={unit_lint['warning_count']}",
        )

    if four_layer["status"] != "PASS":
        add_alert(alerts, "critical", "four_layer_status", f"four-layer status={four_layer['status']}")

    layer_status: Dict[str, str] = four_layer.get("layer_status", {})
    layer_counts: Dict[str, int] = four_layer.get("layer_selected_count", {})

    for layer, min_count in min_layer_counts.items():
        current_count = int(layer_counts.get(layer, 0) or 0)
        current_status = str(layer_status.get(layer, "UNKNOWN")).upper()
        if current_status != "PASS":
            add_alert(alerts, "critical", f"layer_{layer}_status", f"{layer} status={current_status}")
        if current_count < min_count:
            add_alert(
                alerts,
                "critical",
                f"layer_{layer}_count",
                f"{layer} selected_count={current_count} < min={min_count}",
            )

    if previous_snapshot is not None:
        prev_counts = (
            previous_snapshot.get("metrics", {})
            .get("four_layer", {})
            .get("layer_selected_count", {})
        )
        if not isinstance(prev_counts, dict):
            prev_counts = {}

        for layer in min_layer_counts:
            prev = int(prev_counts.get(layer, 0) or 0)
            curr = int(layer_counts.get(layer, 0) or 0)
            if prev <= 0:
                continue
            drop_ratio = float(prev - curr) / float(prev)
            if drop_ratio > max_layer_count_drop_ratio:
                severity = "critical" if drop_ratio > 0.5 else "warning"
                add_alert(
                    alerts,
                    severity,
                    f"layer_{layer}_count_drop",
                    f"{layer} selected_count drop ratio={drop_ratio:.3f} (prev={prev}, curr={curr})",
                )

    return alerts


def summarize_status(alerts: List[Dict[str, str]], fail_on_warning: bool) -> str:
    has_critical = any(alert.get("severity") == "critical" for alert in alerts)
    has_warning = any(alert.get("severity") == "warning" for alert in alerts)

    if has_critical:
        return "FAIL"
    if has_warning and fail_on_warning:
        return "FAIL"
    return "PASS"


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Drift Monitoring Dashboard")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- snapshot_index: {report['snapshot_index']}")
    lines.append(f"- history_size: {report['history_size']}")
    lines.append("")

    metrics = report["metrics"]
    schema = metrics["schema"]
    output_contract = metrics["output_contract"]
    unit_lint = metrics["unit_lint"]
    four_layer = metrics["four_layer"]

    lines.append("## Core Metrics")
    lines.append("")
    lines.append("| metric | value |")
    lines.append("|---|---|")
    lines.append(f"| schema_status | {schema['status']} |")
    lines.append(f"| schema_module_pass_rate | {schema['module_pass_rate']:.3f} |")
    lines.append(f"| output_contract_status | {output_contract['status']} |")
    lines.append(f"| output_contract_module_pass_rate | {output_contract['module_pass_rate']:.3f} |")
    lines.append(f"| unit_lint_status | {unit_lint['status']} |")
    lines.append(f"| unit_lint_error_count | {unit_lint['error_count']} |")
    lines.append(f"| unit_lint_warning_count | {unit_lint['warning_count']} |")
    lines.append(f"| four_layer_status | {four_layer['status']} |")

    layer_counts = four_layer.get("layer_selected_count", {})
    for layer in ["unit", "component", "system", "physics"]:
        lines.append(f"| layer_{layer}_selected_count | {int(layer_counts.get(layer, 0) or 0)} |")

    alerts = report.get("alerts", [])
    if alerts:
        lines.append("")
        lines.append("## Alerts")
        lines.append("")
        for alert in alerts:
            lines.append(
                f"- [{alert.get('severity','unknown').upper()}] {alert.get('code','unknown')}: {alert.get('message','')}"
            )

    recent = report.get("recent_snapshots", [])
    if recent:
        lines.append("")
        lines.append("## Recent Snapshots")
        lines.append("")
        lines.append("| timestamp_utc | status | unit | component | system | physics |")
        lines.append("|---|---|---:|---:|---:|---:|")
        for item in recent:
            layer = item.get("metrics", {}).get("four_layer", {}).get("layer_selected_count", {})
            lines.append(
                f"| {item.get('timestamp_utc','')} | {item.get('status','')} | {int(layer.get('unit',0) or 0)} | {int(layer.get('component',0) or 0)} | {int(layer.get('system',0) or 0)} | {int(layer.get('physics',0) or 0)} |"
            )

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="ENG-006 drift monitoring dashboard gate"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--build-dir", default="build_codex")
    parser.add_argument("--config", default="Debug")
    parser.add_argument("--auto-run-missing", action="store_true")

    parser.add_argument("--schema-report", default="config_schema_v1_gate.json")
    parser.add_argument("--output-contract-report", default="output_contract_v1_gate.json")
    parser.add_argument("--unit-lint-report", default="surface_leo_unit_lint_report.json")
    parser.add_argument("--four-layer-report", default="four_layer_gate_report.json")

    parser.add_argument("--history-json", default="drift_monitor_dashboard_history.json")
    parser.add_argument("--report-json", default="drift_monitor_dashboard.json")
    parser.add_argument("--report-md", default="drift_monitor_dashboard.md")

    parser.add_argument("--max-history-snapshots", type=int, default=200)
    parser.add_argument("--recent-snapshot-window", type=int, default=10)
    parser.add_argument("--max-layer-count-drop-ratio", type=float, default=0.2)

    parser.add_argument("--min-unit-count", type=int, default=20)
    parser.add_argument("--min-component-count", type=int, default=20)
    parser.add_argument("--min-system-count", type=int, default=8)
    parser.add_argument("--min-physics-count", type=int, default=4)
    parser.add_argument("--fail-on-warning", action="store_true")

    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    build_dir = Path(args.build_dir)
    if not build_dir.is_absolute():
        build_dir = (project_root / build_dir).resolve()

    def resolve_from_build(raw_path: str) -> Path:
        path = Path(raw_path)
        if path.is_absolute():
            return path
        return (build_dir / path).resolve()

    schema_report = resolve_from_build(args.schema_report)
    output_contract_report = resolve_from_build(args.output_contract_report)
    unit_lint_report = resolve_from_build(args.unit_lint_report)
    four_layer_report = resolve_from_build(args.four_layer_report)

    history_json = resolve_from_build(args.history_json)
    report_json = resolve_from_build(args.report_json)
    report_md = resolve_from_build(args.report_md)

    min_layer_counts = {
        "unit": args.min_unit_count,
        "component": args.min_component_count,
        "system": args.min_system_count,
        "physics": args.min_physics_count,
    }

    try:
        ensure_report(
            schema_report,
            "UnifiedConfigSchemaV1Gate",
            project_root,
            build_dir,
            args.config,
            args.auto_run_missing,
        )
        ensure_report(
            output_contract_report,
            "UnifiedOutputContractV1Gate",
            project_root,
            build_dir,
            args.config,
            args.auto_run_missing,
        )
        ensure_report(
            unit_lint_report,
            "SurfaceLeoUnitLintGate",
            project_root,
            build_dir,
            args.config,
            args.auto_run_missing,
        )
        ensure_report(
            four_layer_report,
            "FourLayerGateSystemGate",
            project_root,
            build_dir,
            args.config,
            args.auto_run_missing,
        )

        schema_payload = load_json(schema_report)
        output_contract_payload = load_json(output_contract_report)
        unit_lint_payload = load_json(unit_lint_report)
        four_layer_payload = load_json(four_layer_report)

        metrics = {
            "schema": extract_schema_metrics(schema_payload),
            "output_contract": extract_output_contract_metrics(output_contract_payload),
            "unit_lint": extract_unit_lint_metrics(unit_lint_payload),
            "four_layer": extract_four_layer_metrics(four_layer_payload),
        }

        history = load_history(history_json)
        previous_snapshot = history["snapshots"][-1] if history["snapshots"] else None

        snapshot = {
            "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "metrics": metrics,
        }

        alerts = evaluate_alerts(
            snapshot=snapshot,
            previous_snapshot=previous_snapshot,
            min_layer_counts=min_layer_counts,
            max_layer_count_drop_ratio=args.max_layer_count_drop_ratio,
        )
        snapshot["alerts"] = alerts
        snapshot["status"] = summarize_status(alerts, args.fail_on_warning)

        snapshots = history["snapshots"]
        snapshots.append(snapshot)
        if len(snapshots) > max(1, args.max_history_snapshots):
            snapshots[:] = snapshots[-args.max_history_snapshots :]

        history_payload = {
            "schema_version": "v1",
            "snapshots": snapshots,
        }

        history_json.parent.mkdir(parents=True, exist_ok=True)
        history_json.write_text(json.dumps(history_payload, indent=2, ensure_ascii=False), encoding="utf-8")

        status = snapshot["status"]
        report = {
            "status": status,
            "timestamp_utc": snapshot["timestamp_utc"],
            "snapshot_index": len(snapshots),
            "history_size": len(snapshots),
            "metrics": metrics,
            "alerts": alerts,
            "recent_snapshots": snapshots[-max(1, args.recent_snapshot_window) :],
            "history_json": str(history_json),
            "source_reports": {
                "schema_report": str(schema_report),
                "output_contract_report": str(output_contract_report),
                "unit_lint_report": str(unit_lint_report),
                "four_layer_report": str(four_layer_report),
            },
            "thresholds": {
                "max_layer_count_drop_ratio": args.max_layer_count_drop_ratio,
                "min_layer_counts": min_layer_counts,
                "fail_on_warning": args.fail_on_warning,
            },
        }

        report_json.parent.mkdir(parents=True, exist_ok=True)
        report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        write_markdown(report, report_md)

        print(f"report_json={report_json}")
        print(f"report_md={report_md}")
        print(f"history_json={history_json}")
        print(f"status={status}")

        return 0 if status == "PASS" else 2

    except Exception as exc:  # pylint: disable=broad-except
        print(f"status=FAIL({exc})")
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
