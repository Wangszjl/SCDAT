#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

TEMPLATE_REQUIREMENTS = {
    "model_description_template_v1.md": {
        "tokens": [
            "{{MODULE_NAME}}",
            "{{MODEL_NAME}}",
            "{{OWNER}}",
            "{{GENERATED_UTC}}",
            "{{OBJECTIVE}}",
            "{{CORE_RULES}}",
        ],
        "headings": [
            "## 1. Model Overview",
            "## 2. Inputs and Outputs",
            "## 3. Assumptions and Preconditions",
            "## 4. Parameters and Configuration",
            "## 5. Scope and Compatibility",
        ],
    },
    "validation_description_template_v1.md": {
        "tokens": [
            "{{MODULE_NAME}}",
            "{{FEATURE_ID}}",
            "{{OWNER}}",
            "{{GENERATED_UTC}}",
            "{{VALIDATION_COMMANDS}}",
            "{{GATE_STATUS}}",
        ],
        "headings": [
            "## 1. Validation Objective",
            "## 2. Baselines and Datasets",
            "## 3. Validation Commands",
            "## 4. Results and Evidence",
            "## 5. Pass/Fail Decision",
        ],
    },
    "limitations_description_template_v1.md": {
        "tokens": [
            "{{MODULE_NAME}}",
            "{{FEATURE_ID}}",
            "{{OWNER}}",
            "{{GENERATED_UTC}}",
            "{{NUMERICAL_LIMITATIONS}}",
            "{{ROLLBACK_STEPS}}",
        ],
        "headings": [
            "## 1. Applicable Scope",
            "## 2. Known Limitations",
            "## 3. Risks and Impact",
            "## 4. Mitigations and Workarounds",
            "## 5. Degradation and Rollback Guidance",
        ],
    },
}

CORE_RUNLOG_SECTION_ALIASES = {
    "任务信息": ["任务信息", "任务背景", "任务概述", "目标"],
    "代码变更": ["代码变更", "代码修改", "实现变更", "变更内容"],
    "构建与验证": [
        "构建与验证",
        "构建与测试",
        "构建与测试记录",
        "验证命令与结果",
        "验证记录",
    ],
    "结论": ["结论", "总结"],
}

STRICT_RUNLOG_SECTION_ALIASES = {
    "报告产物": ["报告产物", "报告与产物", "门禁报告", "关键报告快照", "产物清单"],
}

DEFAULT_POLICY_EFFECTIVE_DATE = "2026-04-07"


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def has_markdown_heading(text: str, keyword: str) -> bool:
    pattern = re.compile(r"^##\s+.*" + re.escape(keyword) + r".*$", re.MULTILINE)
    return pattern.search(text) is not None


def has_any_markdown_heading(text: str, keywords: List[str]) -> bool:
    return any(has_markdown_heading(text, keyword) for keyword in keywords)


def parse_runlog_date(path: Path) -> datetime | None:
    match = re.search(r"(\d{4}-\d{2}-\d{2})\.md$", path.name)
    if not match:
        return None
    try:
        return datetime.strptime(match.group(1), "%Y-%m-%d")
    except ValueError:
        return None


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Documentation Closure Policy Gate")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- template_dir: {report['template_dir']}")
    lines.append(f"- runlog_count: {report['runlog_count']}")
    lines.append("")

    lines.append("## Template Check")
    lines.append("")
    lines.append("| template | status |")
    lines.append("|---|---|")
    for item in report["templates"]:
        lines.append(f"| {item['template']} | {item['status']} |")

    lines.append("")
    lines.append("## Runlog Structure Check")
    lines.append("")
    lines.append("| runlog | status |")
    lines.append("|---|---|")
    for item in report["runlogs"]:
        lines.append(f"| {item['path']} | {item['status']} |")

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
        description="ENG-009 documentation closure policy gate"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--template-dir", default="documentation/templates")
    parser.add_argument("--runlog-glob", default="documentation/*_stage_runlog_*.md")
    parser.add_argument("--min-runlog-count", type=int, default=10)
    parser.add_argument(
        "--policy-effective-date",
        default=DEFAULT_POLICY_EFFECTIVE_DATE,
        help="Runlogs on or after this date must contain strict closure sections",
    )
    parser.add_argument(
        "--report-json",
        default="build/documentation_closure_gate.json",
    )
    parser.add_argument(
        "--report-md",
        default="build/documentation_closure_gate.md",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()

    template_dir = Path(args.template_dir)
    if not template_dir.is_absolute():
        template_dir = (project_root / template_dir).resolve()

    report_json = Path(args.report_json)
    if not report_json.is_absolute():
        report_json = (project_root / report_json).resolve()

    report_md = Path(args.report_md)
    if not report_md.is_absolute():
        report_md = (project_root / report_md).resolve()

    failures: List[str] = []
    warnings: List[str] = []
    template_results: List[Dict[str, str]] = []
    runlog_results: List[Dict[str, str]] = []

    try:
        policy_effective_date = datetime.strptime(args.policy_effective_date, "%Y-%m-%d")

        if not template_dir.is_dir():
            raise FileNotFoundError(f"template dir not found: {template_dir}")

        for template_name, req in TEMPLATE_REQUIREMENTS.items():
            template_path = template_dir / template_name
            if not template_path.is_file():
                failures.append(f"missing_template:{template_path}")
                template_results.append({"template": str(template_path), "status": "FAIL"})
                continue

            text = load_text(template_path)
            status = "PASS"
            for token in req["tokens"]:
                if token not in text:
                    failures.append(f"missing_template_token:{template_name}:{token}")
                    status = "FAIL"
            for heading in req["headings"]:
                if heading not in text:
                    failures.append(f"missing_template_heading:{template_name}:{heading}")
                    status = "FAIL"

            template_results.append({"template": str(template_path), "status": status})

        runlog_paths = sorted(project_root.glob(args.runlog_glob))
        if len(runlog_paths) < max(1, args.min_runlog_count):
            failures.append(
                f"insufficient_runlog_count:{len(runlog_paths)}<min={max(1, args.min_runlog_count)}"
            )

        for runlog in runlog_paths:
            text = load_text(runlog)
            status = "PASS"

            runlog_date = parse_runlog_date(runlog)
            strict_mode = runlog_date is not None and runlog_date >= policy_effective_date

            if runlog_date is None:
                warnings.append(f"runlog_date_not_detected:{runlog}")

            for section, aliases in CORE_RUNLOG_SECTION_ALIASES.items():
                if not has_any_markdown_heading(text, aliases):
                    if strict_mode:
                        failures.append(f"missing_runlog_section:{runlog}:{section}")
                        status = "FAIL"
                    else:
                        warnings.append(f"legacy_runlog_missing_section:{runlog}:{section}")

            for section, aliases in STRICT_RUNLOG_SECTION_ALIASES.items():
                if not has_any_markdown_heading(text, aliases):
                    if strict_mode:
                        failures.append(f"missing_runlog_section:{runlog}:{section}")
                        status = "FAIL"
                    else:
                        warnings.append(f"legacy_runlog_missing_section:{runlog}:{section}")

            runlog_results.append({"path": str(runlog), "status": status})

        status = "PASS" if not failures else "FAIL"
        report = {
            "status": status,
            "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "project_root": str(project_root),
            "template_dir": str(template_dir),
            "runlog_glob": args.runlog_glob,
            "runlog_count": len(runlog_paths),
            "templates": template_results,
            "runlogs": runlog_results,
            "warnings": warnings,
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
