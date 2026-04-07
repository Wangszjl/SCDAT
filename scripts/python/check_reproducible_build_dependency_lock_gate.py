#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple


def load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected object json in {path}")
    return payload


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_version_command(command: List[str]) -> str:
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if completed.returncode != 0:
        return "not_found"

    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if lines:
        return lines[0]

    err_lines = [line.strip() for line in completed.stderr.splitlines() if line.strip()]
    if err_lines:
        return err_lines[0]

    return "unknown"


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Reproducible Build Dependency Lock Gate")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- lock_json: {report['lock_json']}")
    lines.append("")

    lines.append("## Locked Artifacts")
    lines.append("")
    lines.append("| artifact | path | expected_sha256 | actual_sha256 | status |")
    lines.append("|---|---|---|---|---|")
    for item in report["artifacts"]:
        lines.append(
            f"| {item['artifact']} | {item['path']} | {item['expected_sha256']} | {item['actual_sha256']} | {item['status']} |"
        )

    lines.append("")
    lines.append("## Host Toolchain Snapshot")
    lines.append("")
    lines.append("| tool | version |")
    lines.append("|---|---|")
    for key, value in report["host_toolchain"].items():
        lines.append(f"| {key} | {value} |")

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
        description="ENG-007 reproducible build/dependency lock gate"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument(
        "--lock-json",
        default="scripts/reproducible/toolchain_lock_v1.json",
    )
    parser.add_argument(
        "--report-json",
        default="build/reproducible_build_lock_gate.json",
    )
    parser.add_argument(
        "--report-md",
        default="build/reproducible_build_lock_gate.md",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()

    lock_json = Path(args.lock_json)
    if not lock_json.is_absolute():
        lock_json = (project_root / lock_json).resolve()

    report_json = Path(args.report_json)
    if not report_json.is_absolute():
        report_json = (project_root / report_json).resolve()

    report_md = Path(args.report_md)
    if not report_md.is_absolute():
        report_md = (project_root / report_md).resolve()

    failures: List[str] = []
    artifact_reports: List[Dict[str, str]] = []

    try:
        if not lock_json.is_file():
            raise FileNotFoundError(f"missing lock json: {lock_json}")

        lock_payload = load_json(lock_json)
        if str(lock_payload.get("schema_version", "")).strip() != "v1":
            failures.append("invalid_schema_version")

        paths_raw = lock_payload.get("paths", {})
        hashes_raw = lock_payload.get("sha256", {})
        if not isinstance(paths_raw, dict):
            failures.append("paths_block_missing_or_invalid")
            paths_raw = {}
        if not isinstance(hashes_raw, dict):
            failures.append("sha256_block_missing_or_invalid")
            hashes_raw = {}

        required_artifacts = [
            "python_requirements_lock",
            "dockerfile",
            "run_script_ps1",
            "run_script_sh",
        ]

        for artifact in required_artifacts:
            path_text = str(paths_raw.get(artifact, "")).strip()
            expected_sha = str(hashes_raw.get(artifact, "")).strip().lower()

            if not path_text:
                failures.append(f"missing_path:{artifact}")
                continue
            if not expected_sha or expected_sha == "to_be_filled":
                failures.append(f"missing_expected_sha256:{artifact}")

            artifact_path = Path(path_text)
            if not artifact_path.is_absolute():
                artifact_path = (project_root / artifact_path).resolve()

            if not artifact_path.is_file():
                failures.append(f"missing_artifact_file:{artifact}:{artifact_path}")
                artifact_reports.append(
                    {
                        "artifact": artifact,
                        "path": str(artifact_path),
                        "expected_sha256": expected_sha or "",
                        "actual_sha256": "missing",
                        "status": "FAIL",
                    }
                )
                continue

            actual_sha = sha256_file(artifact_path)
            status = "PASS"
            if expected_sha and expected_sha != "to_be_filled" and actual_sha != expected_sha:
                status = "FAIL"
                failures.append(
                    f"sha256_mismatch:{artifact}:expected={expected_sha}:actual={actual_sha}"
                )

            artifact_reports.append(
                {
                    "artifact": artifact,
                    "path": str(artifact_path),
                    "expected_sha256": expected_sha,
                    "actual_sha256": actual_sha,
                    "status": status,
                }
            )

        host_toolchain = {
            "python": sys.version.split()[0],
            "cmake": run_version_command(["cmake", "--version"]),
            "ninja": run_version_command(["ninja", "--version"]),
            "g++": run_version_command(["g++", "--version"]),
            "clang++": run_version_command(["clang++", "--version"]),
        }

        status = "PASS" if not failures else "FAIL"
        report = {
            "status": status,
            "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "project_root": str(project_root),
            "lock_json": str(lock_json),
            "artifacts": artifact_reports,
            "host_toolchain": host_toolchain,
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
