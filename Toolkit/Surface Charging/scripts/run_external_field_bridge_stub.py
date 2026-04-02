#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import pathlib
from collections import defaultdict
from typing import Any, Dict, List


def clamp(value: float, minimum: float, maximum: float) -> float:
    return max(minimum, min(maximum, value))


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def derive_node_result(node: Dict[str, Any]) -> Dict[str, float | str]:
    potential_v = safe_float(node.get("potential_v"))
    total_current_density = safe_float(node.get("total_current_density_a_per_m2"))
    normal_field_v_per_m = safe_float(node.get("normal_field_v_per_m"))
    local_charge_density_c_per_m3 = safe_float(node.get("local_charge_density_c_per_m3"))
    diagonal_f = max(abs(safe_float(node.get("graph_capacitance_diagonal_f"), 1.0)), 1.0e-30)
    row_sum_f = max(abs(safe_float(node.get("graph_capacitance_row_sum_f"), diagonal_f)), diagonal_f)
    graph_ratio = row_sum_f / diagonal_f
    graph_coupling = max(0.0, graph_ratio - 1.0)
    incoming_scale = safe_float(node.get("field_solver_capacitance_scale"), 1.0)

    reference_shift_v = clamp(total_current_density, -25.0, 25.0) * (0.02 + 0.03 * min(graph_coupling, 5.0))
    reference_potential_v = potential_v + reference_shift_v
    capacitance_scale = clamp(incoming_scale * (1.0 + 0.10 * min(graph_coupling, 2.0)), 0.5, 2.0)
    adjusted_field_v_per_m = normal_field_v_per_m * (1.0 + 0.05 * min(graph_coupling, 3.0))
    adjusted_charge_density_c_per_m3 = local_charge_density_c_per_m3 * (1.0 + 0.05 * min(graph_coupling, 3.0))

    return {
        "node_id": str(node.get("node_id", "")),
        "reference_potential_v": reference_potential_v,
        "normal_field_v_per_m": adjusted_field_v_per_m,
        "local_charge_density_c_per_m3": adjusted_charge_density_c_per_m3,
        "capacitance_scale": capacitance_scale,
    }


def aggregate_boundary_group_results(request: Dict[str, Any], node_results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    mapping_groups: Dict[str, List[str]] = defaultdict(list)
    for mapping in request.get("boundary_mappings", []):
        node_id = str(mapping.get("node_id", ""))
        boundary_group_id = str(mapping.get("boundary_group_id", ""))
        if node_id and boundary_group_id:
            mapping_groups[boundary_group_id].append(node_id)

    node_result_by_id = {str(entry.get("node_id", "")): entry for entry in node_results}
    aggregated: List[Dict[str, Any]] = []
    for boundary_group_id, node_ids in mapping_groups.items():
        matched = [node_result_by_id[node_id] for node_id in node_ids if node_id in node_result_by_id]
        if not matched:
            continue
        count = float(len(matched))
        aggregated.append(
            {
                "node_id": boundary_group_id,
                "reference_potential_v": sum(safe_float(item.get("reference_potential_v")) for item in matched) / count,
                "normal_field_v_per_m": sum(safe_float(item.get("normal_field_v_per_m")) for item in matched) / count,
                "local_charge_density_c_per_m3": sum(safe_float(item.get("local_charge_density_c_per_m3")) for item in matched) / count,
                "capacitance_scale": sum(safe_float(item.get("capacitance_scale"), 1.0) for item in matched) / count,
            }
        )
    return aggregated


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate a stub external field-solver result from a SCDAT field-request JSON."
    )
    parser.add_argument("request_json")
    parser.add_argument("--output", default="", help="Optional explicit result JSON path.")
    parser.add_argument(
        "--include-boundary-group-results",
        action="store_true",
        help="Also emit aggregated results keyed by boundary_group_id for mapped nodes.",
    )
    args = parser.parse_args()

    request_path = pathlib.Path(args.request_json)
    request = json.loads(request_path.read_text(encoding="utf-8"))
    result_path = pathlib.Path(args.output) if args.output else request_path.with_name(
        request_path.name.replace(".field_request.json", ".field_result.json")
    )

    node_results = [derive_node_result(node) for node in request.get("nodes", [])]
    if args.include_boundary_group_results:
        node_results.extend(aggregate_boundary_group_results(request, node_results))

    result = {
        "schema_version": "scdat.external_field_result.v1",
        "generator": "run_external_field_bridge_stub.py",
        "source_request": str(request_path),
        "nodes": node_results,
    }
    result_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print(f"request_json={request_path}")
    print(f"result_json={result_path}")
    print(f"node_result_count={len(node_results)}")
    if node_results:
        print(f"first_node_id={node_results[0]['node_id']}")
        print(f"first_capacitance_scale={safe_float(node_results[0]['capacitance_scale']):.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
