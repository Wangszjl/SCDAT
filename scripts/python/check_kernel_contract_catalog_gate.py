#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple


REQUIRED_MODULES = ("surface", "vacuum_arc", "radiation", "internal", "plasma")
REQUIRED_CONTRACT_IDS = {
    "surface": {
        "surface-pic-shared-particle-transport-domain-v1",
        "surface-pic-global-particle-domain-v1",
        "surface-pic-global-particle-repository-v1",
        "surface-pic-global-sheath-field-solve-v1",
    },
    "vacuum_arc": {
        "arcpic-6-stage-v1",
        "vacuum-arc-collision-emission-channel-v1",
        "vacuum-arc-benchmark-metrics-v1",
    },
    "radiation": {
        "geant4-aligned-deposition-record-v1",
        "geant4-aligned-process-history-v1",
        "radiation-material-governance-v1",
        "radiation-transport-benchmark-v1",
    },
    "internal": {
        "internal-radiation-drive-provenance-v1",
        "internal-radiation-layer-alignment-v1",
        "internal-response-coupled-solve-v1",
    },
    "plasma": {
        "plasma-physics-diagnostics-v1",
        "plasma-reaction-balance-v1",
    },
}


def load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected object json in {path}")
    return payload


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Kernel Contract Catalog Gate Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- catalog: {report['catalog_path']}")
    lines.append("")
    lines.append("| module | status |")
    lines.append("|---|---|")
    for item in report["modules"]:
        lines.append(f"| {item['module']} | {item['status']} |")

    failures = report.get("failures", [])
    if failures:
        lines.append("")
        lines.append("## Failures")
        lines.append("")
        for failure in failures:
            lines.append(f"- {failure}")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_contract_record(record: Any, module: str, failures: List[str], seen_ids: Set[str]) -> None:
    if not isinstance(record, dict):
        failures.append(f"{module}:invalid_contract_record_type")
        return
    contract_id = str(record.get("id", "")).strip()
    if not contract_id:
        failures.append(f"{module}:missing_contract_id")
        return
    if contract_id in seen_ids:
        failures.append(f"{module}:duplicate_contract_id:{contract_id}")
    seen_ids.add(contract_id)

    scope = str(record.get("scope", "")).strip()
    if not scope:
        failures.append(f"{module}:missing_scope:{contract_id}")
    artifact_kind = str(record.get("artifact_kind", "")).strip()
    if artifact_kind not in {"json", "csv", "metadata"}:
        failures.append(f"{module}:invalid_artifact_kind:{contract_id}:{artifact_kind}")
    if artifact_kind in {"json", "csv"}:
        suffix = str(record.get("artifact_suffix", "")).strip()
        if not suffix.startswith("."):
            failures.append(f"{module}:invalid_artifact_suffix:{contract_id}:{suffix}")


def collect_module_contracts(payload: Dict[str, Any], failures: List[str]) -> Tuple[Dict[str, Set[str]], List[Dict[str, str]]]:
    modules = payload.get("modules")
    module_status: List[Dict[str, str]] = []
    contract_map: Dict[str, Set[str]] = {}

    if not isinstance(modules, list):
        failures.append("missing_or_invalid_modules_array")
        return contract_map, module_status

    for module_record in modules:
        if not isinstance(module_record, dict):
            failures.append("invalid_module_record_type")
            continue
        module_name = str(module_record.get("module", "")).strip()
        if not module_name:
            failures.append("missing_module_name")
            continue
        contracts = module_record.get("contracts")
        if not isinstance(contracts, list):
            failures.append(f"{module_name}:missing_contracts_array")
            module_status.append({"module": module_name, "status": "FAIL"})
            continue

        seen_ids: Set[str] = set()
        for contract in contracts:
            validate_contract_record(contract, module_name, failures, seen_ids)
        contract_map[module_name] = seen_ids
        module_status.append({"module": module_name, "status": "PASS"})

    return contract_map, module_status


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate kernel contract catalog integrity and required contract coverage."
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument(
        "--catalog",
        default="documentation/contracts/kernel_contract_catalog_v1.json",
    )
    parser.add_argument(
        "--report-json",
        default="build/kernel_contract_catalog_gate.json",
    )
    parser.add_argument(
        "--report-md",
        default="build/kernel_contract_catalog_gate.md",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    catalog_path = Path(args.catalog)
    if not catalog_path.is_absolute():
        catalog_path = (project_root / catalog_path).resolve()
    report_json = Path(args.report_json)
    if not report_json.is_absolute():
        report_json = (project_root / report_json).resolve()
    report_md = Path(args.report_md)
    if not report_md.is_absolute():
        report_md = (project_root / report_md).resolve()

    failures: List[str] = []
    module_status: List[Dict[str, str]] = []

    try:
        if not catalog_path.is_file():
            failures.append(f"missing_catalog:{catalog_path}")
            payload: Dict[str, Any] = {}
        else:
            payload = load_json(catalog_path)

        if payload:
            if payload.get("schema_version") != "scdat.kernel_contract_catalog.v1":
                failures.append(
                    "catalog_schema_version_mismatch:"
                    f"{payload.get('schema_version')}"
                )
            if payload.get("catalog_version") != "v1":
                failures.append(
                    f"catalog_version_mismatch:{payload.get('catalog_version')}"
                )
            contract_map, module_status = collect_module_contracts(payload, failures)
            module_names = set(contract_map.keys())

            for module in REQUIRED_MODULES:
                if module not in module_names:
                    failures.append(f"missing_required_module:{module}")
                    continue
                missing_ids = REQUIRED_CONTRACT_IDS[module] - contract_map[module]
                if missing_ids:
                    failures.append(
                        f"{module}:missing_required_contract_ids:{','.join(sorted(missing_ids))}"
                    )

            if module_status:
                failed_modules = {
                    failure.split(":", 1)[0]
                    for failure in failures
                    if ":" in failure
                }
                for row in module_status:
                    if row["module"] in failed_modules:
                        row["status"] = "FAIL"

    except Exception as exc:  # pylint: disable=broad-except
        failures.append(f"unexpected_error:{exc}")

    if not module_status:
        module_status = [{"module": module, "status": "FAIL"} for module in REQUIRED_MODULES]

    status = "PASS" if not failures else "FAIL"
    report = {
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "catalog_path": str(catalog_path),
        "modules": module_status,
        "failures": failures,
    }

    report_json.parent.mkdir(parents=True, exist_ok=True)
    report_json.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    write_markdown(report, report_md)

    print(f"report_json={report_json}")
    print(f"report_md={report_md}")
    print(f"status={status}")
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
