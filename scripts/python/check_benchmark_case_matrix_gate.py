#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


REQUIRED_TRACKS = ("surface", "vacuum_arc", "internal_radiation")
MIN_CASES_PER_TRACK = 3


def load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected object json in {path}")
    return payload


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Benchmark Case Matrix Gate Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- matrix: {report['matrix_path']}")
    lines.append("")
    lines.append("| track | status | case_count |")
    lines.append("|---|---|---:|")
    for row in report["tracks"]:
        lines.append(f"| {row['track']} | {row['status']} | {row['case_count']} |")

    failures = report.get("failures", [])
    if failures:
        lines.append("")
        lines.append("## Failures")
        lines.append("")
        for failure in failures:
            lines.append(f"- {failure}")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_case(track: str, case: Any, failures: List[str]) -> None:
    if not isinstance(case, dict):
        failures.append(f"{track}:invalid_case_record")
        return
    case_id = str(case.get("id", "")).strip()
    if not case_id:
        failures.append(f"{track}:missing_case_id")
    if not str(case.get("reference_source", "")).strip():
        failures.append(f"{track}:missing_reference_source:{case_id or '<unknown>'}")
    metrics = case.get("metrics")
    if not isinstance(metrics, list) or not metrics:
        failures.append(f"{track}:missing_metrics:{case_id or '<unknown>'}")
    tolerance = case.get("tolerance")
    if not isinstance(tolerance, dict):
        failures.append(f"{track}:missing_tolerance:{case_id or '<unknown>'}")
    else:
        if not any(key.endswith("_max") for key in tolerance.keys()):
            failures.append(f"{track}:invalid_tolerance_keys:{case_id or '<unknown>'}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate benchmark case matrix coverage and schema completeness."
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument(
        "--matrix",
        default="documentation/contracts/benchmark_case_matrix_v1.json",
    )
    parser.add_argument(
        "--report-json",
        default="build/benchmark_case_matrix_gate.json",
    )
    parser.add_argument(
        "--report-md",
        default="build/benchmark_case_matrix_gate.md",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    matrix_path = Path(args.matrix)
    if not matrix_path.is_absolute():
        matrix_path = (project_root / matrix_path).resolve()

    report_json = Path(args.report_json)
    if not report_json.is_absolute():
        report_json = (project_root / report_json).resolve()
    report_md = Path(args.report_md)
    if not report_md.is_absolute():
        report_md = (project_root / report_md).resolve()

    failures: List[str] = []
    track_rows: List[Dict[str, Any]] = []

    try:
        if not matrix_path.is_file():
            failures.append(f"missing_matrix:{matrix_path}")
            payload: Dict[str, Any] = {}
        else:
            payload = load_json(matrix_path)

        if payload:
            if payload.get("schema_version") != "scdat.benchmark_case_matrix.v1":
                failures.append(
                    f"schema_version_mismatch:{payload.get('schema_version')}"
                )
            if payload.get("matrix_version") != "v1":
                failures.append(f"matrix_version_mismatch:{payload.get('matrix_version')}")

            tracks = payload.get("tracks")
            if not isinstance(tracks, list):
                failures.append("missing_tracks_array")
                tracks = []

            track_map = {}
            for entry in tracks:
                if not isinstance(entry, dict):
                    failures.append("invalid_track_record")
                    continue
                track = str(entry.get("track", "")).strip()
                if not track:
                    failures.append("missing_track_name")
                    continue
                cases = entry.get("cases")
                if not isinstance(cases, list):
                    failures.append(f"{track}:missing_cases_array")
                    cases = []
                for case in cases:
                    validate_case(track, case, failures)
                track_map[track] = cases
                track_rows.append(
                    {
                        "track": track,
                        "status": "PASS",
                        "case_count": len(cases),
                    }
                )

            for track in REQUIRED_TRACKS:
                if track not in track_map:
                    failures.append(f"missing_required_track:{track}")
                    track_rows.append({"track": track, "status": "FAIL", "case_count": 0})
                    continue
                if len(track_map[track]) < MIN_CASES_PER_TRACK:
                    failures.append(
                        f"{track}:insufficient_cases:{len(track_map[track])}"
                    )
                    for row in track_rows:
                        if row["track"] == track:
                            row["status"] = "FAIL"
                            break

            for row in track_rows:
                if any(failure.startswith(row["track"] + ":") for failure in failures):
                    row["status"] = "FAIL"

    except Exception as exc:  # pylint: disable=broad-except
        failures.append(f"unexpected_error:{exc}")

    status = "PASS" if not failures else "FAIL"
    report = {
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "matrix_path": str(matrix_path),
        "tracks": track_rows,
        "failures": failures,
    }
    report_json.parent.mkdir(parents=True, exist_ok=True)
    report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    write_markdown(report, report_md)

    print(f"report_json={report_json}")
    print(f"report_md={report_md}")
    print(f"status={status}")
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())

