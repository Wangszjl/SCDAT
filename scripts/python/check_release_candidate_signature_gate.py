#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


def resolve_path(base: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path.resolve()
    return (base / path).resolve()


def load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected object json in {path}")
    return payload


def sha256_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            hasher.update(chunk)
    return hasher.hexdigest()


def rel_path(path: Path, project_root: Path) -> str:
    try:
        return path.resolve().relative_to(project_root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def collect_artifact_digests(paths: List[Path], project_root: Path) -> List[Dict[str, Any]]:
    digests: List[Dict[str, Any]] = []
    for path in sorted({p.resolve() for p in paths}, key=lambda item: str(item).lower()):
        if not path.is_file():
            continue
        digests.append(
            {
                "path": rel_path(path, project_root),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return digests


def parse_legacy_warning_inventory(
    warning_entries: List[str], project_root: Path
) -> Dict[str, Any]:
    section_counts: Dict[str, int] = {}
    runlog_map: Dict[str, Dict[str, Any]] = {}
    malformed_entries: List[str] = []

    legacy_prefix = "legacy_runlog_missing_section:"
    date_unknown_prefix = "runlog_date_not_detected:"

    for warning in warning_entries:
        if warning.startswith(legacy_prefix):
            payload = warning[len(legacy_prefix) :]
            parts = payload.split(":")
            if len(parts) < 2:
                malformed_entries.append(warning)
                continue

            section = parts[-1].strip()
            runlog_raw = ":".join(parts[:-1]).strip()
            runlog_path = Path(runlog_raw)
            runlog_key = rel_path(runlog_path, project_root)

            section_counts[section] = section_counts.get(section, 0) + 1
            if runlog_key not in runlog_map:
                runlog_map[runlog_key] = {"path": runlog_key, "count": 0, "sections": set()}
            runlog_map[runlog_key]["count"] += 1
            runlog_map[runlog_key]["sections"].add(section)
            continue

        if warning.startswith(date_unknown_prefix):
            runlog_raw = warning[len(date_unknown_prefix) :].strip()
            runlog_path = Path(runlog_raw)
            runlog_key = rel_path(runlog_path, project_root)

            section = "runlog_date_not_detected"
            section_counts[section] = section_counts.get(section, 0) + 1
            if runlog_key not in runlog_map:
                runlog_map[runlog_key] = {"path": runlog_key, "count": 0, "sections": set()}
            runlog_map[runlog_key]["count"] += 1
            runlog_map[runlog_key]["sections"].add(section)
            continue

        malformed_entries.append(warning)

    runlogs: List[Dict[str, Any]] = []
    for _, item in sorted(
        runlog_map.items(),
        key=lambda pair: (-int(pair[1]["count"]), str(pair[1]["path"]).lower()),
    ):
        runlogs.append(
            {
                "path": str(item["path"]),
                "count": int(item["count"]),
                "sections": sorted(str(value) for value in item["sections"]),
            }
        )

    section_summary = [
        {"section": section, "count": count}
        for section, count in sorted(section_counts.items(), key=lambda pair: (-pair[1], pair[0]))
    ]

    return {
        "legacy_warning_count": len(warning_entries),
        "legacy_warning_parseable_count": sum(section_counts.values()),
        "legacy_warning_section_summary": section_summary,
        "legacy_warning_by_runlog": runlogs,
        "malformed_warning_count": len(malformed_entries),
        "malformed_warning_samples": malformed_entries[:10],
    }


def write_warning_inventory_markdown(report: Dict[str, Any], output_path: Path) -> None:
    lines: List[str] = []
    lines.append("# Documentation Legacy Warning Inventory")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- warning_count: {report['warning_count']}")
    lines.append(f"- warning_budget: {report['warning_budget']}")
    lines.append(f"- within_budget: {report['within_budget']}")
    lines.append("")

    lines.append("## Section Summary")
    lines.append("")
    lines.append("| section | count |")
    lines.append("|---|---:|")
    for item in report.get("section_summary", []):
        lines.append(f"| {item['section']} | {item['count']} |")

    lines.append("")
    lines.append("## Runlog Breakdown")
    lines.append("")
    lines.append("| runlog | count | sections |")
    lines.append("|---|---:|---|")
    for item in report.get("runlogs", []):
        sections = ", ".join(item.get("sections", []))
        lines.append(f"| {item['path']} | {item['count']} | {sections} |")

    malformed_count = int(report.get("malformed_warning_count", 0) or 0)
    if malformed_count > 0:
        lines.append("")
        lines.append("## Malformed Warning Samples")
        lines.append("")
        for entry in report.get("malformed_warning_samples", []):
            lines.append(f"- {entry}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_git_command(project_root: Path, args: List[str]) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=str(project_root),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if completed.returncode != 0:
        return ""
    return completed.stdout.strip()


def collect_git_snapshot(project_root: Path) -> Dict[str, Any]:
    commit = run_git_command(project_root, ["rev-parse", "HEAD"])
    status_output = run_git_command(project_root, ["status", "--porcelain"])
    status_lines = [line for line in status_output.splitlines() if line.strip()]

    return {
        "commit": commit or "UNKNOWN",
        "dirty": len(status_lines) > 0,
        "dirty_file_count": len(status_lines),
        "dirty_samples": status_lines[:20],
    }


def write_gate_markdown(report: Dict[str, Any], output_path: Path) -> None:
    lines: List[str] = []
    lines.append("# Release Candidate Signature Gate")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- candidate_label: {report['candidate_label']}")
    lines.append(f"- signature_sha256: {report['signature_sha256']}")
    lines.append(f"- warning_count: {report['legacy_warning_count']}")
    lines.append(f"- warning_budget: {report['legacy_warning_budget']}")
    lines.append("")

    lines.append("## Gate Snapshot")
    lines.append("")
    lines.append("| gate | status |")
    lines.append("|---|---|")
    for key, value in report.get("gate_snapshot", {}).items():
        lines.append(f"| {key} | {value} |")

    lines.append("")
    lines.append("## Artifact Digests")
    lines.append("")
    lines.append("| path | size_bytes | sha256 |")
    lines.append("|---|---:|---|")
    for item in report.get("artifact_digests", []):
        lines.append(f"| {item['path']} | {item['size_bytes']} | {item['sha256']} |")

    warnings = report.get("warnings", [])
    if warnings:
        lines.append("")
        lines.append("## Warnings")
        lines.append("")
        for warning in warnings:
            lines.append(f"- {warning}")

    failures = report.get("failures", [])
    if failures:
        lines.append("")
        lines.append("## Failures")
        lines.append("")
        for failure in failures:
            lines.append(f"- {failure}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="M6 release candidate signature gate with legacy warning inventory"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--build-dir", default="build_codex")
    parser.add_argument("--candidate-label", default="m6-release-candidate")

    parser.add_argument(
        "--documentation-closure-report",
        default="documentation_closure_gate.json",
    )
    parser.add_argument(
        "--pareto-report",
        default="pareto_report_workflow_gate.json",
    )
    parser.add_argument(
        "--drift-report",
        default="drift_monitor_dashboard.json",
    )
    parser.add_argument(
        "--reproducible-lock-report",
        default="reproducible_build_lock_gate.json",
    )

    parser.add_argument(
        "--release-matrix",
        default="documentation/release_acceptance_matrix_2026-04-07.md",
    )
    parser.add_argument(
        "--execution-plan",
        default="documentation/full_alignment_execution_plan_2026-04-03.md",
    )

    parser.add_argument(
        "--legacy-warning-budget",
        type=int,
        default=-1,
        help="Maximum allowed legacy warning count. Negative means no budget check.",
    )
    parser.add_argument(
        "--strict-legacy-warning-zero",
        action="store_true",
        help="Fail when legacy warning count is not zero.",
    )

    parser.add_argument(
        "--signature-json",
        default="release_candidate_signature.json",
    )
    parser.add_argument(
        "--warning-inventory-json",
        default="documentation_legacy_warning_inventory.json",
    )
    parser.add_argument(
        "--warning-inventory-md",
        default="documentation_legacy_warning_inventory.md",
    )
    parser.add_argument(
        "--report-json",
        default="release_candidate_signature_gate.json",
    )
    parser.add_argument(
        "--report-md",
        default="release_candidate_signature_gate.md",
    )

    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    build_dir = resolve_path(project_root, args.build_dir)

    documentation_closure_report = resolve_path(build_dir, args.documentation_closure_report)
    pareto_report = resolve_path(build_dir, args.pareto_report)
    drift_report = resolve_path(build_dir, args.drift_report)
    reproducible_lock_report = resolve_path(build_dir, args.reproducible_lock_report)

    release_matrix = resolve_path(project_root, args.release_matrix)
    execution_plan = resolve_path(project_root, args.execution_plan)

    signature_json = resolve_path(build_dir, args.signature_json)
    warning_inventory_json = resolve_path(build_dir, args.warning_inventory_json)
    warning_inventory_md = resolve_path(build_dir, args.warning_inventory_md)
    report_json = resolve_path(build_dir, args.report_json)
    report_md = resolve_path(build_dir, args.report_md)

    failures: List[str] = []
    warnings: List[str] = []

    gate_snapshot: Dict[str, str] = {
        "documentation_closure": "UNKNOWN",
        "pareto_workflow": "UNKNOWN",
        "drift_dashboard": "UNKNOWN",
        "reproducible_lock": "UNKNOWN",
    }

    try:
        required_paths = {
            "documentation_closure_report": documentation_closure_report,
            "pareto_report": pareto_report,
            "drift_report": drift_report,
            "reproducible_lock_report": reproducible_lock_report,
            "release_matrix": release_matrix,
            "execution_plan": execution_plan,
        }

        for key, path in required_paths.items():
            if not path.is_file():
                failures.append(f"missing_required_artifact:{key}:{path}")

        if failures:
            status = "FAIL"
            timestamp_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            report = {
                "status": status,
                "timestamp_utc": timestamp_utc,
                "candidate_label": args.candidate_label,
                "gate_snapshot": gate_snapshot,
                "legacy_warning_count": 0,
                "legacy_warning_budget": args.legacy_warning_budget,
                "signature_sha256": "",
                "artifact_digests": [],
                "warnings": warnings,
                "failures": failures,
                "signature_json": str(signature_json),
                "warning_inventory_json": str(warning_inventory_json),
                "warning_inventory_md": str(warning_inventory_md),
            }
            report_json.parent.mkdir(parents=True, exist_ok=True)
            report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
            write_gate_markdown(report, report_md)
            print(f"report_json={report_json}")
            print(f"report_md={report_md}")
            print("status=FAIL")
            return 2

        documentation_closure_payload = load_json(documentation_closure_report)
        pareto_payload = load_json(pareto_report)
        drift_payload = load_json(drift_report)
        reproducible_payload = load_json(reproducible_lock_report)

        gate_snapshot = {
            "documentation_closure": str(
                documentation_closure_payload.get("status", "UNKNOWN")
            ).upper(),
            "pareto_workflow": str(pareto_payload.get("status", "UNKNOWN")).upper(),
            "drift_dashboard": str(drift_payload.get("status", "UNKNOWN")).upper(),
            "reproducible_lock": str(reproducible_payload.get("status", "UNKNOWN")).upper(),
        }

        for gate_name, gate_status in gate_snapshot.items():
            if gate_status != "PASS":
                failures.append(f"gate_not_pass:{gate_name}:{gate_status}")

        raw_warnings = documentation_closure_payload.get("warnings", [])
        if not isinstance(raw_warnings, list):
            failures.append("documentation_closure_warnings_not_list")
            raw_warnings = []

        warning_entries = [str(item) for item in raw_warnings]
        warning_inventory = parse_legacy_warning_inventory(warning_entries, project_root)
        legacy_warning_count = int(warning_inventory["legacy_warning_count"])

        if args.strict_legacy_warning_zero and legacy_warning_count > 0:
            failures.append(
                f"legacy_warning_nonzero_under_strict_zero:{legacy_warning_count}"
            )

        if args.legacy_warning_budget >= 0 and legacy_warning_count > args.legacy_warning_budget:
            failures.append(
                f"legacy_warning_budget_exceeded:{legacy_warning_count}>budget={args.legacy_warning_budget}"
            )

        if args.legacy_warning_budget >= 0 and legacy_warning_count < args.legacy_warning_budget:
            warnings.append(
                f"legacy_warning_below_budget:{legacy_warning_count}<budget={args.legacy_warning_budget}"
            )

        if int(warning_inventory.get("malformed_warning_count", 0) or 0) > 0:
            warnings.append(
                f"legacy_warning_parse_partial:malformed_count={warning_inventory['malformed_warning_count']}"
            )

        git_snapshot = collect_git_snapshot(project_root)
        if str(git_snapshot.get("commit", "UNKNOWN")).upper() == "UNKNOWN":
            warnings.append("git_commit_unavailable")
        if bool(git_snapshot.get("dirty", False)):
            warnings.append(
                f"git_worktree_dirty:file_count={int(git_snapshot.get('dirty_file_count', 0) or 0)}"
            )

        artifact_paths: List[Path] = [
            documentation_closure_report,
            pareto_report,
            drift_report,
            reproducible_lock_report,
            release_matrix,
            execution_plan,
        ]

        for candidate in [
            documentation_closure_report.with_suffix(".md"),
            pareto_report.with_suffix(".md"),
            drift_report.with_suffix(".md"),
            reproducible_lock_report.with_suffix(".md"),
        ]:
            if candidate.is_file():
                artifact_paths.append(candidate)

        artifact_digests = collect_artifact_digests(artifact_paths, project_root)

        timestamp_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        signature_payload: Dict[str, Any] = {
            "schema_version": "release_candidate_signature_v1",
            "timestamp_utc": timestamp_utc,
            "candidate_label": args.candidate_label,
            "project_root": str(project_root),
            "build_dir": str(build_dir),
            "gate_snapshot": gate_snapshot,
            "legacy_warning_summary": {
                "warning_count": legacy_warning_count,
                "warning_budget": args.legacy_warning_budget,
                "strict_zero": bool(args.strict_legacy_warning_zero),
                "section_summary": warning_inventory["legacy_warning_section_summary"],
            },
            "artifact_digests": artifact_digests,
            "git": git_snapshot,
        }

        signature_digest_source = json.dumps(
            signature_payload, ensure_ascii=False, sort_keys=True
        ).encode("utf-8")
        signature_sha256 = hashlib.sha256(signature_digest_source).hexdigest()
        signature_payload["signature_sha256"] = signature_sha256

        signature_json.parent.mkdir(parents=True, exist_ok=True)
        signature_json.write_text(
            json.dumps(signature_payload, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

        warning_inventory_report = {
            "status": "PASS",
            "timestamp_utc": timestamp_utc,
            "warning_count": legacy_warning_count,
            "warning_budget": args.legacy_warning_budget,
            "within_budget": args.legacy_warning_budget < 0
            or legacy_warning_count <= args.legacy_warning_budget,
            "section_summary": warning_inventory["legacy_warning_section_summary"],
            "runlogs": warning_inventory["legacy_warning_by_runlog"],
            "malformed_warning_count": warning_inventory["malformed_warning_count"],
            "malformed_warning_samples": warning_inventory["malformed_warning_samples"],
        }

        warning_inventory_json.parent.mkdir(parents=True, exist_ok=True)
        warning_inventory_json.write_text(
            json.dumps(warning_inventory_report, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        write_warning_inventory_markdown(warning_inventory_report, warning_inventory_md)

        status = "PASS" if not failures else "FAIL"
        report = {
            "status": status,
            "timestamp_utc": timestamp_utc,
            "candidate_label": args.candidate_label,
            "gate_snapshot": gate_snapshot,
            "legacy_warning_count": legacy_warning_count,
            "legacy_warning_budget": args.legacy_warning_budget,
            "strict_legacy_warning_zero": bool(args.strict_legacy_warning_zero),
            "signature_json": str(signature_json),
            "signature_sha256": signature_sha256,
            "warning_inventory_json": str(warning_inventory_json),
            "warning_inventory_md": str(warning_inventory_md),
            "artifact_digests": artifact_digests,
            "warnings": warnings,
            "failures": failures,
        }

        report_json.parent.mkdir(parents=True, exist_ok=True)
        report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        write_gate_markdown(report, report_md)

        print(f"report_json={report_json}")
        print(f"report_md={report_md}")
        print(f"signature_json={signature_json}")
        print(f"warning_inventory_json={warning_inventory_json}")
        print(f"status={status}")
        return 0 if status == "PASS" else 2

    except Exception as exc:  # pylint: disable=broad-except
        print(f"status=FAIL({exc})")
        return 3


if __name__ == "__main__":
    raise SystemExit(main())