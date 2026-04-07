#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


def build_required_paths(project_root: Path) -> List[Path]:
    return [
        project_root / "documentation" / "full_alignment_execution_plan_2026-04-03.md",
        project_root / "documentation" / "contracts" / "surface_leo_acceptance_contract_v1.json",
        project_root / "documentation" / "vacuum_arc_arcpic_mapping_v1.md",
        project_root / "documentation" / "radiation_geant4_process_matrix_v1.md",
        project_root / "documentation" / "plasma_multiscale_architecture_spec_v1.md",
        project_root / "documentation" / "full_alignment_execution_status_board_2026-04-03.md",
        project_root / "results" / "baselines" / "vacuum_arc" / "microgap_flashover_baseline_2026-04-02.csv",
        project_root / "results" / "baselines" / "vacuum_arc" / "triple_junction_flashover_baseline_2026-04-02.csv",
        project_root / "results" / "baselines" / "vacuum_arc" / "restrike_recovery_baseline_2026-04-02.csv",
        project_root / "results" / "baselines" / "vacuum_arc" / "vacuum_arc_baseline_metrics_2026-04-02.json",
    ]


def summarize(required_paths: List[Path], project_root: Path) -> Dict[str, object]:
    present: List[str] = []
    missing: List[str] = []

    for path in required_paths:
        rel = path.relative_to(project_root).as_posix()
        if path.exists():
            present.append(rel)
        else:
            missing.append(rel)

    status = "PASS" if not missing else "FAIL"
    return {
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "required_count": len(required_paths),
        "present_count": len(present),
        "missing_count": len(missing),
        "present": present,
        "missing": missing,
    }


def write_markdown(report: Dict[str, object], output_md: Path) -> None:
    output_md.parent.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    lines.append("# Baseline Artifacts Validation Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- required_count: {report['required_count']}")
    lines.append(f"- present_count: {report['present_count']}")
    lines.append(f"- missing_count: {report['missing_count']}")
    lines.append("")

    lines.append("## Present")
    lines.append("")
    present = report.get("present", [])
    if present:
        for item in present:
            lines.append(f"- {item}")
    else:
        lines.append("- (none)")

    lines.append("")
    lines.append("## Missing")
    lines.append("")
    missing = report.get("missing", [])
    if missing:
        for item in missing:
            lines.append(f"- {item}")
    else:
        lines.append("- (none)")

    output_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate baseline artifacts required by the execution plan.")
    parser.add_argument("--project-root", default=".", help="Project root path.")
    parser.add_argument(
        "--report-json",
        default="results/baseline_artifacts_validation_2026-04-03.json",
        help="Output JSON report path (relative to project root unless absolute).",
    )
    parser.add_argument(
        "--report-md",
        default="results/baseline_artifacts_validation_2026-04-03.md",
        help="Output Markdown report path (relative to project root unless absolute).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return non-zero when any required artifact is missing.",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    report_json = Path(args.report_json)
    report_md = Path(args.report_md)

    if not report_json.is_absolute():
        report_json = project_root / report_json
    if not report_md.is_absolute():
        report_md = project_root / report_md

    required_paths = build_required_paths(project_root)
    report = summarize(required_paths, project_root)

    report_json.parent.mkdir(parents=True, exist_ok=True)
    report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    write_markdown(report, report_md)

    print(f"report_json={report_json}")
    print(f"report_md={report_md}")
    print(f"status={report['status']}")

    if args.strict and report["status"] != "PASS":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
