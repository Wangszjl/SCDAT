#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import re
from dataclasses import dataclass
from typing import Any


def _load_netcdf_module():
    try:
        import netCDF4  # type: ignore
    except ModuleNotFoundError as exc:  # pragma: no cover - environment dependent
        raise SystemExit(
            "netCDF4 is required to read SPIS Java outputs. Install it in the active Python environment."
        ) from exc
    return netCDF4


@dataclass
class SeriesMatch:
    status: str
    series_name: str | None = None
    values: list[float] | None = None
    source_kind: str = "missing"
    scale: float = 1.0
    note: str = ""


def load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_scdat_csv(path: pathlib.Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    return rows, fieldnames


def column_as_floats(rows: list[dict[str, str]], field: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        text = (row.get(field) or "").strip()
        if not text:
            values.append(float("nan"))
            continue
        try:
            values.append(float(text))
        except ValueError:
            values.append(float("nan"))
    return values


def interpolate_series(x_values: list[float], y_values: list[float], query: float) -> float:
    if not x_values or not y_values:
        return float("nan")
    if query <= x_values[0]:
        return y_values[0]
    if query >= x_values[-1]:
        return y_values[-1]
    for index in range(1, len(x_values)):
        x0 = x_values[index - 1]
        x1 = x_values[index]
        if query > x1:
            continue
        y0 = y_values[index - 1]
        y1 = y_values[index]
        if abs(x1 - x0) <= 1.0e-30:
            return y1
        ratio = (query - x0) / (x1 - x0)
        return y0 + (y1 - y0) * ratio
    return y_values[-1]


def discover_reference_fields(reference_root: pathlib.Path) -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []
    for family in ("DataFieldMonitored", "DataFieldExtracted"):
        family_root = reference_root / family
        if not family_root.exists():
            continue
        for candidate in sorted(family_root.glob("*.nc")):
            if "-VERTEXMask-" in candidate.name or "--Compared" in candidate.name or "--Reference" in candidate.name:
                continue
            items.append({"family": family, "file": candidate.name, "path": candidate.as_posix()})
    return items


def load_java_series(path: pathlib.Path) -> tuple[list[float], str]:
    netcdf4 = _load_netcdf_module()
    with netcdf4.Dataset(path) as dataset:
        values = dataset.variables["dataArray"][:]
        flat = [float(item) for item in values.reshape(-1)]
        unit = str(getattr(dataset, "unit", "")) or str(getattr(dataset.variables["dataArray"], "unit", ""))
    return flat, unit


def resolve_target_file(reference_root: pathlib.Path, family: str, pattern: str) -> pathlib.Path | None:
    family_root = reference_root / family
    if not family_root.exists():
        return None
    candidates = [
        candidate
        for candidate in sorted(family_root.glob(pattern))
        if "-VERTEXMask-" not in candidate.name and "--Compared" not in candidate.name and "--Reference" not in candidate.name
    ]
    return candidates[0] if candidates else None


def load_optional_metadata(csv_path: pathlib.Path) -> dict[str, str]:
    metadata_path = csv_path.with_name(csv_path.name + ".metadata.json")
    if not metadata_path.exists():
        return {}
    document = load_json(metadata_path)
    metadata = document.get("metadata", {})
    if not isinstance(metadata, dict):
        return {}
    return {str(key): str(value) for key, value in metadata.items()}


def parse_boolish(value: str | None) -> bool | None:
    if value is None:
        return None
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes"}:
        return True
    if normalized in {"0", "false", "no"}:
        return False
    return None


def expected_source_resolved_series(config: dict[str, Any]) -> list[str]:
    config_block = config.get("config", {})
    targets = config_block.get("comparison_targets", [])
    expected: list[str] = []
    for target in targets:
        source_key = (target.get("source_key") or "").strip()
        if not source_key:
            continue
        hinted = (target.get("scdat_series_hint") or "").strip()
        if hinted and hinted not in expected:
            expected.append(hinted)
    return expected


def build_source_resolved_checks(
    config: dict[str, Any],
    fieldnames: list[str],
    metadata: dict[str, str],
) -> dict[str, Any]:
    config_block = config.get("config", {})
    spis_import = config_block.get("spis_import", {})
    source_keys = [str(item).strip() for item in spis_import.get("source_keys", []) if str(item).strip()]
    comparison_target_source_keys = [
        str(item).strip()
        for item in spis_import.get("comparison_target_source_keys", [])
        if str(item).strip()
    ]
    numkernel_only_source_keys = [
        str(item).strip()
        for item in spis_import.get("numkernel_only_source_keys", [])
        if str(item).strip()
    ]
    expected_series = expected_source_resolved_series(config)
    missing_series = [name for name in expected_series if name not in fieldnames]
    return {
        "source_keys": source_keys,
        "source_key_count": len(source_keys),
        "reference_source_keys": comparison_target_source_keys,
        "reference_source_key_count": len(comparison_target_source_keys),
        "numkernel_only_source_keys": numkernel_only_source_keys,
        "numkernel_only_source_key_count": len(numkernel_only_source_keys),
        "expected_series_count": len(expected_series),
        "missing_required_series": missing_series,
        "missing_required_series_count": len(missing_series),
        "bookkeeping_mode": metadata.get("surface_source_resolved_bookkeeping_mode", ""),
        "slot_count_complete": parse_boolish(metadata.get("surface_source_resolved_slot_count_complete")),
        "strict_attribution_available": parse_boolish(
            metadata.get("surface_source_resolved_strict_attribution_available")
        ),
    }


def resolve_scdat_series(
    target: dict[str, Any],
    fieldnames: list[str],
    rows: list[dict[str, str]],
    config: dict[str, Any],
) -> SeriesMatch:
    target_name = target["name"]
    field_set = set(fieldnames)
    config_block = config.get("config", {})
    surface_area_m2 = float(config_block.get("surface_area_m2", 1.0))
    node_areas = [
        float(node.get("area_m2", surface_area_m2))
        for node in config_block.get("surface_nodes", [])
    ]

    def exact(name: str) -> SeriesMatch:
        if name in field_set:
            return SeriesMatch(
                status="matched",
                series_name=name,
                values=column_as_floats(rows, name),
                source_kind="exact",
            )
        return SeriesMatch(status="unsupported", note=f"Missing SCDAT series: {name}")

    def density(name: str, area: float) -> SeriesMatch:
        if name in field_set:
            base = column_as_floats(rows, name)
            return SeriesMatch(
                status="matched",
                series_name=name,
                values=[value * area for value in base],
                source_kind="density_scaled",
                scale=area,
                note=f"Scaled by area_m2={area}",
            )
        return SeriesMatch(status="unsupported", note=f"Missing SCDAT density series: {name}")

    scdat_series_hint = (target.get("scdat_series_hint") or "").strip()
    if scdat_series_hint:
        hinted = exact(scdat_series_hint)
        if hinted.status == "matched":
            return hinted

    source_key = (target.get("source_key") or "").strip()
    node_index = target.get("node_index")
    if source_key:
        if re.fullmatch(rf"Individual_current_on_spacecraft_Collected_{re.escape(source_key)}", target_name):
            return exact(f"surface_source_{source_key}_spacecraft_collected_current_a")
        if re.fullmatch(rf"Individual_current_on_spacecraft_Interactor_{re.escape(source_key)}", target_name):
            return exact(f"surface_source_{source_key}_spacecraft_interactor_current_a")
        if re.fullmatch(rf"Number_of_{re.escape(source_key)}", target_name):
            return exact(f"surface_source_{source_key}_superparticle_count")
        if node_index is not None:
            try:
                node_id = int(node_index)
            except (TypeError, ValueError):
                node_id = None
            if node_id is not None:
                if re.fullmatch(rf"Collected_current_{re.escape(source_key)}_node_{node_id}", target_name):
                    return exact(f"surface_node_{node_id}_source_{source_key}_collected_current_a")
                if re.fullmatch(rf"Interactor_current_{re.escape(source_key)}_node_{node_id}", target_name):
                    return exact(f"surface_node_{node_id}_source_{source_key}_interactor_current_a")

    mapping_table = {
        "Average_surface_potential_of_node_0": lambda: exact("surface_node_0_potential_v"),
        "Average_surface_potential_of_node_1": lambda: exact("surface_node_1_potential_v"),
        "Total_current_on_node_0": lambda: density(
            "surface_node_0_total_current_density_a_per_m2",
            node_areas[0] if len(node_areas) > 0 else surface_area_m2,
        ),
        "Total_current_on_node_1": lambda: density(
            "surface_node_1_total_current_density_a_per_m2",
            node_areas[1] if len(node_areas) > 1 else surface_area_m2,
        ),
        "Total_current_on_spacecraft": lambda: exact("surface_interface_0_current_a")
        if "surface_interface_0_current_a" in field_set
        else exact("circuit_branch_current_a")
        if "circuit_branch_current_a" in field_set
        else density(
            "surface_node_0_total_current_density_a_per_m2",
            node_areas[0] if len(node_areas) > 0 else surface_area_m2,
        )
        if len(node_areas) == 1 and "surface_node_0_total_current_density_a_per_m2" in field_set
        else density("total_current_density_a_per_m2", surface_area_m2),
    }
    resolver = mapping_table.get(target_name)
    if resolver is None:
        return SeriesMatch(
            status="unsupported",
            note="No SCDAT mapping rule has been implemented for this SPIS field.",
        )
    return resolver()


def finite_pairs(left: list[float], right: list[float]) -> list[tuple[float, float]]:
    pairs: list[tuple[float, float]] = []
    for left_value, right_value in zip(left, right):
        if math.isfinite(left_value) and math.isfinite(right_value):
            pairs.append((left_value, right_value))
    return pairs


def build_error_stats(reference: list[float], candidate: list[float]) -> dict[str, Any]:
    pairs = finite_pairs(reference, candidate)
    if not pairs:
        return {"sample_count": 0, "mae": None, "max_abs_error": None, "rmse": None}
    abs_errors = [abs(left - right) for left, right in pairs]
    sq_errors = [(left - right) ** 2 for left, right in pairs]
    return {
        "sample_count": len(pairs),
        "mae": sum(abs_errors) / len(abs_errors),
        "max_abs_error": max(abs_errors),
        "rmse": math.sqrt(sum(sq_errors) / len(sq_errors)),
    }


def summarize_discrete_case(
    manifest_entry: dict[str, Any],
    field_name: str,
) -> tuple[float | None, float | None]:
    csv_path = pathlib.Path(manifest_entry["output_csv"])
    if not csv_path.exists():
        return None, None
    rows, _ = load_scdat_csv(csv_path)
    if not rows:
        return None, None
    last = rows[-1]
    try:
        time_value = float(last.get("time_s", "nan"))
    except ValueError:
        time_value = float("nan")
    try:
        field_value = float(last.get(field_name, "nan"))
    except ValueError:
        field_value = float("nan")
    return (time_value if math.isfinite(time_value) else None, field_value if math.isfinite(field_value) else None)


def compare(config_path: pathlib.Path, scdat_csv_path: pathlib.Path, discrete_manifest_path: pathlib.Path | None) -> dict[str, Any]:
    config = load_json(config_path)
    config_block = config.get("config", {})
    reference_cfg = config_block.get("spis_reference_output", {})
    reference_root = pathlib.Path(reference_cfg["output_root"])
    comparison_csv_path = pathlib.Path(reference_cfg["comparison_csv_path"])
    comparison_summary_path = pathlib.Path(reference_cfg["comparison_summary_path"])
    comparison_json_path = pathlib.Path(reference_cfg["comparison_json_path"])
    comparison_csv_path.parent.mkdir(parents=True, exist_ok=True)

    rows, fieldnames = load_scdat_csv(scdat_csv_path)
    metadata = load_optional_metadata(scdat_csv_path)
    time_series = column_as_floats(rows, "time_s")
    reference_fields = discover_reference_fields(reference_root)
    discrete_manifest = load_json(discrete_manifest_path) if discrete_manifest_path and discrete_manifest_path.exists() else {"cases": []}
    source_resolved_checks = build_source_resolved_checks(config, fieldnames, metadata)

    unified_rows: list[dict[str, Any]] = []
    field_summaries: list[dict[str, Any]] = []
    aligned_count = 0
    adapted_count = 0
    unsupported_count = 0

    for target in config_block.get("comparison_targets", []):
        target_name = target["name"]
        family = target["source_family"]
        pattern = target["source_pattern"]
        java_file = resolve_target_file(reference_root, family, pattern)
        if java_file is None:
            unsupported_count += 1
            field_summaries.append(
                {
                    "name": target_name,
                    "status": "unsupported",
                    "reason": f"No Java reference file matched pattern {pattern}",
                }
            )
            continue

        java_values, java_unit = load_java_series(java_file)
        sample_count = min(len(time_series), len(java_values))
        java_values = java_values[:sample_count]
        aligned_time = time_series[:sample_count]
        scdat_match = resolve_scdat_series(target, fieldnames, rows[:sample_count], config)
        if scdat_match.status != "matched" or scdat_match.values is None:
            unsupported_count += 1
            field_summaries.append(
                {
                    "name": target_name,
                    "status": "unsupported",
                    "reason": scdat_match.note,
                    "java_file": java_file.as_posix(),
                    "unit": target.get("unit") or java_unit,
                }
            )
            continue

        scdat_values = scdat_match.values[:sample_count]
        stats = build_error_stats(java_values, scdat_values)
        if scdat_match.source_kind == "exact":
            aligned_count += 1
            status = "aligned"
        else:
            adapted_count += 1
            status = "adapted"

        for time_value, java_value, scdat_value in zip(aligned_time, java_values, scdat_values):
            if not (math.isfinite(time_value) and math.isfinite(java_value) and math.isfinite(scdat_value)):
                continue
            abs_error = abs(scdat_value - java_value)
            rel_error = abs_error / max(abs(java_value), 1.0e-30)
            unified_rows.append(
                {
                    "comparison_kind": "java_vs_scdat",
                    "target_name": target_name,
                    "status": status,
                    "alignment_axis": "time_s",
                    "axis_value": time_value,
                    "java_value": java_value,
                    "scdat_value": scdat_value,
                    "abs_error": abs_error,
                    "rel_error": rel_error,
                    "unit": target.get("unit") or java_unit,
                    "java_file": java_file.name,
                    "scdat_series": scdat_match.series_name or "",
                    "mapping_kind": scdat_match.source_kind,
                }
            )

        field_summary = {
            "name": target_name,
            "status": status,
            "unit": target.get("unit") or java_unit,
            "java_file": java_file.as_posix(),
            "scdat_series": scdat_match.series_name,
            "mapping_kind": scdat_match.source_kind,
            "mapping_note": scdat_match.note,
            "error_stats": stats,
        }

        discrete_rows: list[dict[str, Any]] = []
        if discrete_manifest.get("cases") and scdat_match.series_name:
            for entry in discrete_manifest["cases"]:
                discrete_time, discrete_value = summarize_discrete_case(entry, scdat_match.series_name)
                if discrete_time is None or discrete_value is None:
                    continue
                pwl_value = interpolate_series(aligned_time, scdat_values, float(entry["time_s"]))
                abs_error = abs(discrete_value - pwl_value)
                rel_error = abs_error / max(abs(pwl_value), 1.0e-30)
                record = {
                    "comparison_kind": "pwl_vs_discrete",
                    "target_name": target_name,
                    "status": "crosscheck",
                    "alignment_axis": "time_s",
                    "axis_value": float(entry["time_s"]),
                    "java_value": "",
                    "scdat_value": pwl_value,
                    "abs_error": abs_error,
                    "rel_error": rel_error,
                    "unit": target.get("unit") or java_unit,
                    "java_file": "",
                    "scdat_series": scdat_match.series_name,
                    "mapping_kind": "crosscheck_discrete",
                    "discrete_value": discrete_value,
                    "discrete_output_csv": entry["output_csv"],
                }
                unified_rows.append(record)
                discrete_rows.append(record)
            if discrete_rows:
                field_summary["discrete_crosscheck"] = build_error_stats(
                    [float(row["scdat_value"]) for row in discrete_rows],
                    [float(row["discrete_value"]) for row in discrete_rows],
                )

        field_summaries.append(field_summary)

    all_target_names = {item["name"] for item in config_block.get("comparison_targets", [])}
    extra_reference_fields = []
    for item in reference_fields:
        simple_name = pathlib.Path(item["file"]).stem
        if not any(simple_name.startswith(prefix) for prefix in all_target_names):
            extra_reference_fields.append(item)

    with comparison_csv_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames_out = [
            "comparison_kind",
            "target_name",
            "status",
            "alignment_axis",
            "axis_value",
            "java_value",
            "scdat_value",
            "discrete_value",
            "abs_error",
            "rel_error",
            "unit",
            "java_file",
            "scdat_series",
            "mapping_kind",
            "discrete_output_csv",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames_out)
        writer.writeheader()
        for row in unified_rows:
            writer.writerow(row)

    summary_payload = {
        "config_path": config_path.as_posix(),
        "scdat_csv_path": scdat_csv_path.as_posix(),
        "reference_root": reference_root.as_posix(),
        "field_coverage": {
            "aligned": aligned_count,
            "adapted": adapted_count,
            "unsupported": unsupported_count,
            "total_targets": len(config_block.get("comparison_targets", [])),
        },
        "fields": field_summaries,
        "source_resolved_checks": source_resolved_checks,
        "discovered_reference_fields": reference_fields,
        "extra_reference_fields": extra_reference_fields,
    }
    comparison_json_path.write_text(json.dumps(summary_payload, indent=2), encoding="utf-8")

    case_name = str(config.get("name") or config.get("module") or "SPIS Case")
    markdown_lines = [
        f"# {case_name} SPIS Comparison Summary",
        "",
        f"- Config: `{config_path.as_posix()}`",
        f"- SCDAT CSV: `{scdat_csv_path.as_posix()}`",
        f"- Java Output: `{reference_root.as_posix()}`",
        "",
        "## Coverage",
        "",
        f"- Aligned fields: {aligned_count}",
        f"- Adapted fields: {adapted_count}",
        f"- Unsupported fields: {unsupported_count}",
        f"- Total requested targets: {len(config_block.get('comparison_targets', []))}",
        "",
        "## Source-Resolved Checks",
        "",
        f"- Source keys discovered in config: {source_resolved_checks['source_key_count']}",
        f"- Java-reference-backed source keys: {source_resolved_checks['reference_source_key_count']}",
        f"- NumKernel-only source keys: {source_resolved_checks['numkernel_only_source_key_count']}",
        f"- Expected source-resolved series: {source_resolved_checks['expected_series_count']}",
        f"- Missing required source-resolved series: {source_resolved_checks['missing_required_series_count']}",
        f"- Bookkeeping mode: `{source_resolved_checks['bookkeeping_mode'] or 'unknown'}`",
        f"- Slot count complete: {source_resolved_checks['slot_count_complete']}",
        f"- Strict attribution available: {source_resolved_checks['strict_attribution_available']}",
        "",
        "## Fields",
        "",
    ]
    if source_resolved_checks["missing_required_series"]:
        for series_name in source_resolved_checks["missing_required_series"]:
            markdown_lines.append(f"- Missing source-resolved series: `{series_name}`")
        markdown_lines.append("")
    for field in field_summaries:
        markdown_lines.append(
            f"- `{field['name']}`: {field['status']}"
            + (f", series `{field.get('scdat_series', '')}`" if field.get("scdat_series") else "")
            + (f", reason: {field['reason']}" if field.get("reason") else "")
        )
        stats = field.get("error_stats")
        if stats and stats.get("sample_count"):
            markdown_lines.append(
                f"  sample_count={stats['sample_count']}, mae={stats['mae']:.6e}, max_abs_error={stats['max_abs_error']:.6e}, rmse={stats['rmse']:.6e}"
            )
        crosscheck = field.get("discrete_crosscheck")
        if crosscheck and crosscheck.get("sample_count"):
            markdown_lines.append(
                f"  pwl_vs_discrete sample_count={crosscheck['sample_count']}, mae={crosscheck['mae']:.6e}, max_abs_error={crosscheck['max_abs_error']:.6e}, rmse={crosscheck['rmse']:.6e}"
            )
    markdown_lines.extend(
        [
            "",
            "## Extra Java Fields",
            "",
            f"- Readable reference fields discovered: {len(reference_fields)}",
            f"- Not explicitly targeted by current mapping: {len(extra_reference_fields)}",
        ]
    )
    comparison_summary_path.write_text("\n".join(markdown_lines) + "\n", encoding="utf-8")
    return summary_payload


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare SCDAT surface-config results with SPIS Java netCDF outputs."
    )
    parser.add_argument("config_json", type=pathlib.Path, help="Generated surface-config JSON path.")
    parser.add_argument("scdat_csv", type=pathlib.Path, help="SCDAT output CSV produced by surface-config.")
    parser.add_argument(
        "--discrete-manifest",
        type=pathlib.Path,
        default=None,
        help="Optional discrete-point manifest JSON emitted by the importer.",
    )
    args = parser.parse_args()
    summary = compare(args.config_json, args.scdat_csv, args.discrete_manifest)
    print(json.dumps(summary["field_coverage"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
