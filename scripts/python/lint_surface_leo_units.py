#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple


MATERIAL_ALIASES: Dict[str, str] = {
    "vac": "vacuum",
    "space": "vacuum",
    "al": "aluminum",
    "al6061": "aluminum",
    "aluminum6061": "aluminum",
    "aluminium": "aluminum",
    "polyimide": "kapton",
    "kaptonhn": "kapton",
    "kapton_hn": "kapton",
    "kapton-hn": "kapton",
    "teflon": "ptfe",
    "polytetrafluoroethylene": "ptfe",
}

CANONICAL_MATERIALS = {"vacuum", "aluminum", "kapton", "ptfe"}

UNITLESS_KEY_CONFLICTS: Dict[str, str] = {
    "bulk_flow_velocity": "bulk_flow_velocity_m_per_s",
    "electron_density": "electron_density_m3",
    "ion_density": "ion_density_m3",
    "electron_temperature": "electron_temperature_ev",
    "ion_temperature": "ion_temperature_ev",
    "dielectric_thickness": "dielectric_thickness_m",
    "surface_temperature": "surface_temperature_k",
    "photon_flux": "photon_flux_m2_s",
}


def normalize_token(text: str) -> str:
    return "".join(ch.lower() for ch in text if ch.isalnum())


def load_json(path: Path) -> Dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid json: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError("top-level json must be an object")
    return payload


def walk_keys(node: Any, prefix: str = "") -> List[Tuple[str, Any]]:
    items: List[Tuple[str, Any]] = []
    if isinstance(node, dict):
        for key, value in node.items():
            path = f"{prefix}.{key}" if prefix else key
            items.append((path, value))
            items.extend(walk_keys(value, path))
    elif isinstance(node, list):
        for index, value in enumerate(node):
            path = f"{prefix}[{index}]"
            items.extend(walk_keys(value, path))
    return items


def get_nested(node: Dict[str, Any], path: str) -> Any:
    current: Any = node
    for key in path.split("."):
        if not isinstance(current, dict) or key not in current:
            return None
        current = current[key]
    return current


def has_explicit_material_definition(material_obj: Dict[str, Any]) -> bool:
    explicit_keys = {
        "id",
        "type",
        "permittivity",
        "conductivity",
        "permeability",
        "work_function_ev",
        "breakdown_field_v_per_m",
        "secondary_electron_yield",
        "mass_density_kg_per_m3",
        "scalar_properties",
    }
    return any(key in material_obj for key in explicit_keys)


def resolve_material_name(raw_name: str) -> str | None:
    token = normalize_token(raw_name)
    if token in CANONICAL_MATERIALS:
        return token
    return MATERIAL_ALIASES.get(token)


def lint_surface_config(payload: Dict[str, Any]) -> Dict[str, Any]:
    errors: List[str] = []
    warnings: List[str] = []

    config = payload.get("config", payload)
    if not isinstance(config, dict):
        errors.append("config must be an object")
        return {"status": "FAIL", "errors": errors, "warnings": warnings}

    all_key_paths = walk_keys(config)
    key_set = {path.split(".")[-1] for path, _ in all_key_paths}

    for unitless_key, canonical_key in UNITLESS_KEY_CONFLICTS.items():
        if unitless_key in key_set:
            errors.append(
                f"unitless key '{unitless_key}' is not allowed; use '{canonical_key}'"
            )

    bulk_flow = config.get("bulk_flow_velocity_m_per_s")
    if isinstance(bulk_flow, (int, float)):
        if not math.isfinite(float(bulk_flow)):
            errors.append("bulk_flow_velocity_m_per_s must be finite")
        elif not (1.0e3 <= float(bulk_flow) <= 1.2e4):
            warnings.append(
                "bulk_flow_velocity_m_per_s outside typical LEO range [1e3, 1.2e4]"
            )

    dielectric_thickness = config.get("dielectric_thickness_m")
    if dielectric_thickness is not None:
        if not isinstance(dielectric_thickness, (int, float)) or float(dielectric_thickness) <= 0.0:
            errors.append("dielectric_thickness_m must be > 0")

    plasma = config.get("plasma_model", config.get("plasma", {}))
    if isinstance(plasma, dict):
        for key in ("electron_density_m3", "ion_density_m3"):
            if key in plasma:
                value = plasma[key]
                if not isinstance(value, (int, float)) or float(value) <= 0.0:
                    errors.append(f"plasma.{key} must be > 0")

        for key in ("electron_temperature_ev", "ion_temperature_ev"):
            if key in plasma:
                value = plasma[key]
                if not isinstance(value, (int, float)) or float(value) <= 0.0:
                    errors.append(f"plasma.{key} must be > 0")

    # Material alias governance checks.
    material_name = config.get("material_name")
    if isinstance(material_name, str):
        resolved = resolve_material_name(material_name)
        if resolved is None:
            errors.append(f"unknown material_name '{material_name}'")

    material_alias = config.get("material_alias")
    if isinstance(material_alias, str):
        resolved = resolve_material_name(material_alias)
        if resolved is None:
            errors.append(f"unknown material_alias '{material_alias}'")

    for material_key in ("surface_material", "material"):
        material_obj = config.get(material_key)
        if not isinstance(material_obj, dict):
            continue

        alias = material_obj.get("alias")
        if isinstance(alias, str) and resolve_material_name(alias) is None:
            errors.append(f"unknown {material_key}.alias '{alias}'")

        name = material_obj.get("name")
        if isinstance(name, str):
            resolved = resolve_material_name(name)
            if resolved is None and not has_explicit_material_definition(material_obj):
                errors.append(
                    f"unknown {material_key}.name '{name}' without explicit material properties"
                )
            elif resolved is None:
                warnings.append(
                    f"{material_key}.name '{name}' is treated as custom material with explicit properties"
                )

    status = "PASS" if not errors else "FAIL"
    return {
        "status": status,
        "errors": errors,
        "warnings": warnings,
        "checked_keys": len(all_key_paths),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Lint surface LEO configuration units and material alias mapping governance."
    )
    parser.add_argument("--config", required=True, help="surface config json path")
    parser.add_argument("--report-json", default="", help="optional report json path")
    parser.add_argument("--strict", action="store_true", help="fail on lint errors")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.is_file():
        print("status=FAIL")
        print(f"error=missing_config:{config_path}")
        return 2

    try:
        payload = load_json(config_path)
        report = lint_surface_config(payload)
    except Exception as exc:  # pylint: disable=broad-except
        print("status=FAIL")
        print(f"error={exc}")
        return 3

    report["config"] = str(config_path)

    if args.report_json:
        report_path = Path(args.report_json)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"report_json={report_path}")

    for warning in report["warnings"]:
        print(f"warning={warning}")
    for error in report["errors"]:
        print(f"error={error}")

    print(f"status={report['status']}")
    if args.strict and report["status"] != "PASS":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
