#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


def load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected top-level object in {path}")
    return payload


def rel_path(path: Path, project_root: Path) -> str:
    try:
        return str(path.relative_to(project_root)).replace("\\", "/")
    except ValueError:
        return str(path).replace("\\", "/")


def load_text(path: Path) -> str:
    raw = path.read_bytes()
    if raw.startswith(b"\xff\xfe") or raw.startswith(b"\xfe\xff"):
        return raw.decode("utf-16", errors="replace")
    text = raw.decode("utf-8", errors="replace")
    if "\x00" in text:
        try:
            return raw.decode("utf-16", errors="replace")
        except UnicodeDecodeError:
            return text
    return text


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Surface Family Coverage Gate Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- matrix: {report['matrix']}")
    lines.append(f"- tracked_family_count: {report['tracked_family_count']}")
    lines.append("")
    lines.append("| family_id | domain | code | test | config_route | status |")
    lines.append("|---|---|---|---|---|---|")
    for item in report["families"]:
        lines.append(
            f"| {item['family_id']} | {item['domain']} | {item['code_status']} | "
            f"{item['test_status']} | {item['config_route_status']} | {item['status']} |"
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
        description="Validate the Surface/SPIS family coverage matrix artifact."
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument(
        "--matrix", default="scripts/run/surface_spis_family_coverage_matrix.json"
    )
    parser.add_argument(
        "--reference-matrix", default="scripts/run/surface_reference_matrix.json"
    )
    parser.add_argument("--output-root", default="build/surface_family_coverage_gate")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    matrix_path = (project_root / args.matrix).resolve()
    reference_matrix_path = (project_root / args.reference_matrix).resolve()
    output_root = (project_root / args.output_root).resolve()
    report_json_path = output_root / "surface_family_coverage_gate.json"
    report_md_path = output_root / "surface_family_coverage_gate.md"

    if not matrix_path.is_file():
        raise FileNotFoundError(f"missing matrix artifact: {matrix_path}")
    if not reference_matrix_path.is_file():
        raise FileNotFoundError(f"missing reference matrix: {reference_matrix_path}")

    matrix = load_json(matrix_path)
    reference_matrix = load_json(reference_matrix_path)

    target_families = matrix.get("target_families", [])
    family_to_code = matrix.get("family_to_code", [])
    family_to_test = matrix.get("family_to_test", [])
    family_to_config_route = matrix.get("family_to_config_route", [])
    if not isinstance(target_families, list):
        raise ValueError("target_families must be a list")
    if not isinstance(family_to_code, list):
        raise ValueError("family_to_code must be a list")
    if not isinstance(family_to_test, list):
        raise ValueError("family_to_test must be a list")
    if not isinstance(family_to_config_route, list):
        raise ValueError("family_to_config_route must be a list")

    reference_cases = reference_matrix.get("cases", [])
    if not isinstance(reference_cases, list):
        raise ValueError("reference matrix cases must be a list")

    cases_by_label: Dict[str, Dict[str, Any]] = {}
    for case in reference_cases:
        if not isinstance(case, dict):
            continue
        label = str(case.get("case_label", case.get("preset", ""))).strip()
        preset = str(case.get("preset", "")).strip()
        if label:
            cases_by_label[label] = case
        if preset:
            cases_by_label.setdefault(preset, case)

    code_by_family = {
        str(entry.get("family_id", "")).strip(): entry
        for entry in family_to_code
        if isinstance(entry, dict)
    }
    test_by_family = {
        str(entry.get("family_id", "")).strip(): entry
        for entry in family_to_test
        if isinstance(entry, dict)
    }
    route_by_family = {
        str(entry.get("family_id", "")).strip(): entry
        for entry in family_to_config_route
        if isinstance(entry, dict)
    }

    failures: List[str] = []
    family_reports: List[Dict[str, Any]] = []

    for family in target_families:
        if not isinstance(family, dict):
            failures.append("target_families:invalid_entry_type")
            continue

        family_id = str(family.get("id", "")).strip()
        domain = str(family.get("domain", "")).strip()
        if not family_id:
            failures.append("target_families:missing_family_id")
            continue

        code_ok = True
        test_ok = True
        route_ok = True

        code_entry = code_by_family.get(family_id)
        if code_entry is None:
            failures.append(f"{family_id}:missing_family_to_code_entry")
            code_ok = False
        else:
            spis_sources = code_entry.get("spis_sources", [])
            code_paths = code_entry.get("code_paths", [])
            documentation_paths = code_entry.get("documentation_paths", [])

            if not isinstance(spis_sources, list) or not spis_sources:
                failures.append(f"{family_id}:missing_spis_sources")
                code_ok = False
            else:
                for pattern in spis_sources:
                    pattern_text = str(pattern).strip()
                    if not pattern_text:
                        failures.append(f"{family_id}:invalid_spis_source_pattern")
                        code_ok = False
                        continue
                    if not list(project_root.glob(pattern_text)):
                        failures.append(f"{family_id}:spis_source_missing:{pattern_text}")
                        code_ok = False

            if not isinstance(code_paths, list) or not code_paths:
                failures.append(f"{family_id}:missing_code_paths")
                code_ok = False
            else:
                for entry in code_paths:
                    resolved = (project_root / str(entry)).resolve()
                    if not resolved.is_file():
                        failures.append(f"{family_id}:missing_code_path:{entry}")
                        code_ok = False

            if not isinstance(documentation_paths, list) or not documentation_paths:
                failures.append(f"{family_id}:missing_documentation_paths")
                code_ok = False
            else:
                for entry in documentation_paths:
                    resolved = (project_root / str(entry)).resolve()
                    if not resolved.is_file():
                        failures.append(f"{family_id}:missing_documentation_path:{entry}")
                        code_ok = False

        test_entry = test_by_family.get(family_id)
        if test_entry is None:
            failures.append(f"{family_id}:missing_family_to_test_entry")
            test_ok = False
        else:
            test_evidence = test_entry.get("test_evidence", [])
            gate_artifacts = test_entry.get("gate_artifacts", [])
            if not isinstance(test_evidence, list) or not test_evidence:
                failures.append(f"{family_id}:missing_test_evidence")
                test_ok = False
            else:
                for evidence in test_evidence:
                    if not isinstance(evidence, dict):
                        failures.append(f"{family_id}:invalid_test_evidence_entry")
                        test_ok = False
                        continue
                    path_text = str(evidence.get("path", "")).strip()
                    if not path_text:
                        failures.append(f"{family_id}:test_evidence_missing_path")
                        test_ok = False
                        continue
                    resolved = (project_root / path_text).resolve()
                    if not resolved.is_file():
                        failures.append(f"{family_id}:missing_test_file:{path_text}")
                        test_ok = False
                        continue
                    content = load_text(resolved)
                    selectors = evidence.get("selectors", [])
                    if isinstance(selectors, list):
                        for selector in selectors:
                            selector_text = str(selector).strip()
                            if selector_text and selector_text not in content:
                                failures.append(
                                    f"{family_id}:missing_test_selector:{path_text}:{selector_text}"
                                )
                                test_ok = False
            if not isinstance(gate_artifacts, list) or not gate_artifacts:
                failures.append(f"{family_id}:missing_gate_artifacts")
                test_ok = False
            else:
                for artifact in gate_artifacts:
                    resolved = (project_root / str(artifact)).resolve()
                    if not resolved.exists():
                        failures.append(f"{family_id}:missing_gate_artifact:{artifact}")
                        test_ok = False

        route_entry = route_by_family.get(family_id)
        if route_entry is None:
            failures.append(f"{family_id}:missing_family_to_config_route_entry")
            route_ok = False
        else:
            metadata_keys = route_entry.get("metadata_keys", [])
            metadata_source_files = route_entry.get("metadata_source_files", [])
            config_keys = route_entry.get("config_keys", [])
            config_source_files = route_entry.get("config_source_files", [])
            preset_names = route_entry.get("preset_names", [])
            reference_matrix_cases = route_entry.get("reference_matrix_cases", [])

            if not isinstance(metadata_keys, list) or not metadata_keys:
                failures.append(f"{family_id}:missing_metadata_keys")
                route_ok = False
            if not isinstance(metadata_source_files, list) or not metadata_source_files:
                failures.append(f"{family_id}:missing_metadata_source_files")
                route_ok = False
            else:
                metadata_texts: List[str] = []
                for entry in metadata_source_files:
                    resolved = (project_root / str(entry)).resolve()
                    if not resolved.is_file():
                        failures.append(f"{family_id}:missing_metadata_source_file:{entry}")
                        route_ok = False
                    else:
                        metadata_texts.append(load_text(resolved))
                combined = "\n".join(metadata_texts)
                for key in metadata_keys:
                    key_text = str(key).strip()
                    if key_text and key_text not in combined:
                        failures.append(f"{family_id}:metadata_key_not_found:{key_text}")
                        route_ok = False

            if not isinstance(config_keys, list) or not config_keys:
                failures.append(f"{family_id}:missing_config_keys")
                route_ok = False
            if not isinstance(config_source_files, list) or not config_source_files:
                failures.append(f"{family_id}:missing_config_source_files")
                route_ok = False
            else:
                config_texts: List[str] = []
                for entry in config_source_files:
                    resolved = (project_root / str(entry)).resolve()
                    if not resolved.is_file():
                        failures.append(f"{family_id}:missing_config_source_file:{entry}")
                        route_ok = False
                    else:
                        config_texts.append(load_text(resolved))
                combined = "\n".join(config_texts)
                for key in config_keys:
                    key_text = str(key).strip()
                    if key_text and key_text not in combined:
                        failures.append(f"{family_id}:config_key_not_found:{key_text}")
                        route_ok = False

            if not isinstance(preset_names, list) or not preset_names:
                failures.append(f"{family_id}:missing_preset_names")
                route_ok = False
            if not isinstance(reference_matrix_cases, list) or not reference_matrix_cases:
                failures.append(f"{family_id}:missing_reference_matrix_cases")
                route_ok = False
            else:
                for case_label in reference_matrix_cases:
                    label = str(case_label).strip()
                    if label and label not in cases_by_label:
                        failures.append(f"{family_id}:missing_reference_matrix_case:{label}")
                        route_ok = False

        family_reports.append(
            {
                "family_id": family_id,
                "domain": domain,
                "code_status": "PASS" if code_ok else "FAIL",
                "test_status": "PASS" if test_ok else "FAIL",
                "config_route_status": "PASS" if route_ok else "FAIL",
                "status": "PASS" if code_ok and test_ok and route_ok else "FAIL",
            }
        )

    report = {
        "status": "PASS" if not failures else "FAIL",
        "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "matrix": rel_path(matrix_path, project_root),
        "reference_matrix": rel_path(reference_matrix_path, project_root),
        "tracked_family_count": len(family_reports),
        "families": family_reports,
        "failures": failures,
    }

    output_root.mkdir(parents=True, exist_ok=True)
    report_json_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_markdown(report, report_md_path)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
