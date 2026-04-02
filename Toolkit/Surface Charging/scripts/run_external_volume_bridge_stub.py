#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pathlib
from typing import Any, Dict, List


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def load_json(path: pathlib.Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate a stub external volume-solver result from a SCDAT volume-request JSON."
    )
    parser.add_argument("volume_request_json")
    parser.add_argument("--output", default="", help="Optional explicit result JSON path.")
    args = parser.parse_args()

    request_path = pathlib.Path(args.volume_request_json)
    request = load_json(request_path)
    artifacts = request.get("artifacts", {})
    volume_stub_path = pathlib.Path(artifacts["volume_stub_json"])
    volume_stub = load_json(volume_stub_path)
    volume_mesh_stub_path = pathlib.Path(artifacts.get("volume_mesh_stub_json", ""))
    volume_mesh_stub = load_json(volume_mesh_stub_path) if volume_mesh_stub_path.exists() else {}
    projection_path = pathlib.Path(artifacts.get("surface_volume_projection_json", ""))
    projection = load_json(projection_path) if projection_path.exists() else {}

    projection_weights = {
        str(entry.get("cell_id", "")): safe_float(entry.get("surface_to_volume_weight"), 1.0)
        for entry in projection.get("projection_weights", [])
        if str(entry.get("cell_id", ""))
    }
    cell_face_counts: Dict[str, int] = {}
    for entry in volume_mesh_stub.get("cell_faces", []):
        cell_id = str(entry.get("cell_id", ""))
        if cell_id:
            cell_face_counts[cell_id] = cell_face_counts.get(cell_id, 0) + 1
    neighbor_counts: Dict[str, int] = {}
    neighbor_geometry: Dict[str, float] = {}
    for entry in volume_mesh_stub.get("cell_neighbors", []):
        source_cell = str(entry.get("source_cell_id", ""))
        target_cell = str(entry.get("target_cell_id", ""))
        shared_face_area = safe_float(entry.get("shared_face_area_m2"), 0.0)
        face_distance = max(1.0e-9, safe_float(entry.get("face_distance_m"), 0.0))
        geometry_gain = shared_face_area / face_distance if shared_face_area > 0.0 else 0.0
        if source_cell:
            neighbor_counts[source_cell] = neighbor_counts.get(source_cell, 0) + 1
            neighbor_geometry[source_cell] = neighbor_geometry.get(source_cell, 0.0) + geometry_gain
        if target_cell:
            neighbor_counts[target_cell] = neighbor_counts.get(target_cell, 0) + 1
            neighbor_geometry[target_cell] = neighbor_geometry.get(target_cell, 0.0) + geometry_gain
    mesh_cell_by_id = {
        str(entry.get("cell_id", "")): entry
        for entry in volume_mesh_stub.get("cells", [])
        if str(entry.get("cell_id", ""))
    }

    result_path = pathlib.Path(args.output) if args.output else request_path.with_name(
        request_path.name.replace(".volume_request.json", ".volume_result.json")
    )

    cells: List[Dict[str, Any]] = []
    for cell in volume_stub.get("cells", []):
        potential_v = safe_float(cell.get("potential_v"))
        reference_potential_v = safe_float(cell.get("reference_potential_v"))
        normal_field_v_per_m = safe_float(cell.get("normal_field_v_per_m"))
        local_charge_density_c_per_m3 = safe_float(cell.get("local_charge_density_c_per_m3"))
        field_solver_capacitance_scale = safe_float(cell.get("field_solver_capacitance_scale"), 1.0)
        cell_id = str(cell.get("cell_id", ""))
        mesh_cell = mesh_cell_by_id.get(cell_id, {})
        linked_face_count = safe_float(mesh_cell.get("linked_face_count"), 0.0)
        linked_face_count = max(linked_face_count, float(cell_face_counts.get(cell_id, 0)))
        characteristic_length_m = safe_float(mesh_cell.get("characteristic_length_m"), 0.0)
        center_x_m = safe_float(mesh_cell.get("center_x_m"), 0.0)
        center_y_m = safe_float(mesh_cell.get("center_y_m"), 0.0)
        center_z_m = safe_float(mesh_cell.get("center_z_m"), 0.0)
        neighbor_count = float(neighbor_counts.get(cell_id, 0))
        geometry_gain = neighbor_geometry.get(cell_id, 0.0)
        projection_weight = projection_weights.get(cell_id, 1.0)
        center_radius_m = (center_x_m ** 2 + center_y_m ** 2 + center_z_m ** 2) ** 0.5

        coupling_gain = max(
            0.0,
            min(1.0, 0.12 * neighbor_count + 0.06 * linked_face_count + 0.02 * geometry_gain),
        )
        blended_potential_v = (
            0.60 * potential_v
            + 0.40 * reference_potential_v
            - 0.02 * coupling_gain * potential_v
            + 0.01 * center_radius_m
        )
        blended_reference_potential_v = (
            0.35 * potential_v + 0.65 * reference_potential_v + 0.005 * geometry_gain
        )
        adjusted_field_v_per_m = normal_field_v_per_m * (1.03 + 0.04 * coupling_gain)
        adjusted_charge_density = local_charge_density_c_per_m3 * (1.02 + 0.05 * coupling_gain)
        adjusted_capacitance_scale = max(
            0.5,
            min(2.5, field_solver_capacitance_scale * (1.06 + 0.08 * coupling_gain + 0.01 * geometry_gain)),
        )
        projection_weight_scale = max(1.0, min(3.0, projection_weight * (1.0 + 0.05 * neighbor_count)))
        sheath_length_scale = max(
            0.75,
            min(1.35, 1.0 + 0.10 * coupling_gain - 0.02 * max(0.0, characteristic_length_m / 0.1 - 1.0)),
        )
        cells.append(
            {
                "cell_id": cell_id,
                "potential_v": blended_potential_v,
                "reference_potential_v": blended_reference_potential_v,
                "normal_field_v_per_m": adjusted_field_v_per_m,
                "local_charge_density_c_per_m3": adjusted_charge_density,
                "capacitance_scale": adjusted_capacitance_scale,
                "coupling_gain": coupling_gain,
                "projection_weight_scale": projection_weight_scale,
                "sheath_length_scale": sheath_length_scale,
            }
        )

    result = {
        "schema_version": "scdat.external_volume_result.v1",
        "generator": "run_external_volume_bridge_stub.py",
        "source_request": str(request_path),
        "cells": cells,
    }
    result_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print(f"volume_request_json={request_path}")
    print(f"volume_result_json={result_path}")
    print(f"cell_result_count={len(cells)}")
    if cells:
        print(f"first_cell_id={cells[0]['cell_id']}")
        print(f"first_potential_v={safe_float(cells[0]['potential_v']):.12g}")
        print(f"first_capacitance_scale={safe_float(cells[0]['capacitance_scale'], 1.0):.12g}")
        print(f"first_coupling_gain={safe_float(cells[0]['coupling_gain'], 0.0):.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
