#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
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


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Surface Round 26 Sealoff Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- family_coverage_gate_status: {report['family_coverage_gate_status']}")
    lines.append(f"- reference_matrix_gate_status: {report['reference_matrix_gate_status']}")
    lines.append(f"- final_alignment_verdict: {report['final_alignment_verdict']}")
    lines.append("")
    lines.append("| alignment_status | count |")
    lines.append("|---|---|")
    for key, value in sorted(report["alignment_status_counts"].items()):
        lines.append(f"| {key} | {value} |")
    lines.append("")
    lines.append("## Unresolved Families")
    lines.append("")
    if report["unresolved_families"]:
        for item in report["unresolved_families"]:
            lines.append(f"- {item['family_id']}: {item['alignment_status']}")
    else:
        lines.append("- none")
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
        description="Build the final Round 26 sealoff report from coverage and reference gates."
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument(
        "--coverage-matrix", default="scripts/run/surface_spis_family_coverage_matrix.json"
    )
    parser.add_argument(
        "--coverage-gate",
        default="build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json",
    )
    parser.add_argument(
        "--reference-gate",
        default="build_codex/surface_reference_matrix_gate_round26_full/surface_reference_matrix_gate.json",
    )
    parser.add_argument(
        "--contract-catalog",
        default="documentation/contracts/kernel_contract_catalog_v1.json",
    )
    parser.add_argument("--output-root", default="build_codex/surface_round26_sealoff")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    coverage_matrix_path = (project_root / args.coverage_matrix).resolve()
    coverage_gate_path = (project_root / args.coverage_gate).resolve()
    reference_gate_path = (project_root / args.reference_gate).resolve()
    contract_catalog_path = (project_root / args.contract_catalog).resolve()
    output_root = (project_root / args.output_root).resolve()
    report_json_path = output_root / "surface_round26_sealoff.json"
    report_md_path = output_root / "surface_round26_sealoff.md"

    for path in [
        coverage_matrix_path,
        coverage_gate_path,
        reference_gate_path,
        contract_catalog_path,
    ]:
        if not path.is_file():
            raise FileNotFoundError(f"missing required input: {path}")

    coverage_matrix = load_json(coverage_matrix_path)
    coverage_gate = load_json(coverage_gate_path)
    reference_gate = load_json(reference_gate_path)
    contract_catalog = load_json(contract_catalog_path)

    target_families = coverage_matrix.get("target_families", [])
    if not isinstance(target_families, list):
        raise ValueError("coverage matrix target_families must be a list")

    family_counter: Counter[str] = Counter()
    unresolved_families: List[Dict[str, str]] = []
    for item in target_families:
        if not isinstance(item, dict):
            continue
        family_id = str(item.get("id", "")).strip()
        alignment_status = str(item.get("alignment_status", "")).strip() or "unknown"
        family_counter[alignment_status] += 1
        if alignment_status != "completed":
            unresolved_families.append(
                {"family_id": family_id, "alignment_status": alignment_status}
            )

    failures: List[str] = []
    if coverage_gate.get("status") != "PASS":
        failures.append(f"family_coverage_gate_not_pass:{coverage_gate.get('status')}")
    if reference_gate.get("status") != "PASS":
        failures.append(f"reference_matrix_gate_not_pass:{reference_gate.get('status')}")

    catalog_modules = contract_catalog.get("modules", [])
    catalog_has_family_coverage_contract = False
    if isinstance(catalog_modules, list):
        for module in catalog_modules:
            if not isinstance(module, dict):
                continue
            if str(module.get("module", "")).strip() != "surface":
                continue
            contracts = module.get("contracts", [])
            if not isinstance(contracts, list):
                continue
            for contract in contracts:
                if not isinstance(contract, dict):
                    continue
                if str(contract.get("id", "")).strip() == "surface-spis-family-coverage-matrix-v1":
                    catalog_has_family_coverage_contract = True
                    break
    if not catalog_has_family_coverage_contract:
        failures.append("kernel_contract_catalog_missing_surface_family_coverage_contract")

    final_alignment_verdict = (
        "achieved" if not unresolved_families else "not_yet_achieved"
    )
    status = "PASS" if not failures else "FAIL"

    report = {
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "coverage_matrix": rel_path(coverage_matrix_path, project_root),
        "coverage_gate": rel_path(coverage_gate_path, project_root),
        "reference_gate": rel_path(reference_gate_path, project_root),
        "contract_catalog": rel_path(contract_catalog_path, project_root),
        "family_coverage_gate_status": str(coverage_gate.get("status", "unknown")),
        "reference_matrix_gate_status": str(reference_gate.get("status", "unknown")),
        "catalog_has_family_coverage_contract": catalog_has_family_coverage_contract,
        "tracked_family_count": len(target_families),
        "alignment_status_counts": dict(family_counter),
        "unresolved_families": unresolved_families,
        "final_alignment_verdict": final_alignment_verdict,
        "failures": failures,
    }

    output_root.mkdir(parents=True, exist_ok=True)
    report_json_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_markdown(report, report_md_path)
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
