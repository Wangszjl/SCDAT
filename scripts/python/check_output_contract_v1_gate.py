#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple


MODULE_CASES = [
    {
        "module": "plasma",
        "command": "plasma-config",
        "config": "scripts/run/plasma_config_full_example.json",
        "default_dt": 5.0e-10,
    },
    {
        "module": "surface",
        "command": "surface-config",
        "config": "scripts/run/surface_config_full_example.json",
        "default_dt": 2.0e-8,
    },
    {
        "module": "internal",
        "command": "internal-config",
        "config": "scripts/run/internal_config_full_example.json",
        "default_dt": 1.0e-3,
    },
    {
        "module": "radiation",
        "command": "radiation-config",
        "config": "scripts/run/radiation_config_full_example.json",
        "default_dt": 5.0e-4,
    },
    {
        "module": "arc",
        "command": "arc-config",
        "config": "scripts/run/arc_config_full_example.json",
        "default_dt": 5.0e-10,
    },
]

REQUIRED_TOP_LEVEL_KEYS = {
    "output_contract_schema",
    "output_contract_version",
    "artifact_format",
    "artifact_path",
    "metadata",
}

REQUIRED_METADATA_KEYS = {
    "module",
    "output_contract_schema",
    "output_contract_version",
    "output_axis_name",
    "output_row_count",
    "output_series_count",
    "output_series_names",
    "output_artifact_csv",
    "output_artifact_vtk",
    "output_artifact_h5",
    "output_export_utc",
}

EXPECTED_FORMAT_SUFFIX = {
    "csv": ".csv",
    "vtk": ".vtk",
    "h5": ".h5",
}

SURFACE_BOUNDARY_OBSERVER_SCHEMA = "scdat.surface_shared_runtime_observer.v1"
SURFACE_BOUNDARY_OBSERVER_CONTRACT = "surface-pic-boundary-observer-v1"
SURFACE_RUNTIME_CONSISTENCY_SCHEMA = "scdat.surface_shared_runtime_consistency.v1"
SURFACE_RUNTIME_CONSISTENCY_CONTRACT = "surface-pic-shared-runtime-consistency-v1"
SURFACE_PARTICLE_TRANSPORT_DOMAIN_SCHEMA = "scdat.surface_shared_particle_transport_domain.v1"
SURFACE_PARTICLE_TRANSPORT_DOMAIN_CONTRACT = (
    "surface-pic-shared-particle-transport-domain-v1"
)
SURFACE_GLOBAL_PARTICLE_DOMAIN_SCHEMA = "scdat.surface_global_particle_domain.v1"
SURFACE_GLOBAL_PARTICLE_DOMAIN_CONTRACT = "surface-pic-global-particle-domain-v1"
SURFACE_GLOBAL_PARTICLE_REPOSITORY_SCHEMA = (
    "scdat.surface_global_particle_repository.v1"
)
SURFACE_GLOBAL_PARTICLE_REPOSITORY_CONTRACT = (
    "surface-pic-global-particle-repository-v1"
)
SURFACE_GLOBAL_SHEATH_FIELD_SOLVE_SCHEMA = "scdat.surface_global_sheath_field_solve.v1"
SURFACE_GLOBAL_SHEATH_FIELD_SOLVE_CONTRACT = (
    "surface-pic-global-sheath-field-solve-v1"
)
RADIATION_DEPOSITION_HISTORY_SCHEMA = "scdat.radiation.deposition_history.v1"
RADIATION_PROCESS_HISTORY_SCHEMA = "scdat.radiation.process_history.v1"
RADIATION_DEPOSITION_HISTORY_CONTRACT = "geant4-aligned-deposition-record-v1"
RADIATION_PROCESS_HISTORY_CONTRACT = "geant4-aligned-process-history-v1"
RADIATION_MATERIAL_GOVERNANCE_CONTRACT = "radiation-material-governance-v1"
RADIATION_TRANSPORT_BENCHMARK_SCHEMA = "scdat.radiation_transport_benchmark.v1"
RADIATION_TRANSPORT_BENCHMARK_CONTRACT = "radiation-transport-benchmark-v1"
VACUUM_ARC_BENCHMARK_METRICS_SCHEMA = "scdat.vacuum_arc.benchmark_metrics.v1"
VACUUM_ARC_BENCHMARK_METRICS_CONTRACT = "vacuum-arc-benchmark-metrics-v1"
VACUUM_ARC_PIPELINE_CONTRACT = "arcpic-6-stage-v1"
VACUUM_ARC_COLLISION_DIAGNOSTIC_CONTRACT = (
    "vacuum-arc-collision-emission-channel-v1"
)
PLASMA_PHYSICS_DIAGNOSTICS_SCHEMA = "scdat.plasma.physics_diagnostics.v1"
PLASMA_PHYSICS_DIAGNOSTICS_CONTRACT = "plasma-physics-diagnostics-v1"
PLASMA_REACTION_BALANCE_CONTRACT = "plasma-reaction-balance-v1"
INTERNAL_DRIVE_PROVENANCE_SCHEMA = "scdat.internal_radiation_drive_provenance.v1"
INTERNAL_DRIVE_PROVENANCE_CONTRACT = "internal-radiation-drive-provenance-v1"
INTERNAL_DRIVE_LAYER_ALIGNMENT_CONTRACT = "internal-radiation-layer-alignment-v1"
INTERNAL_RESPONSE_COUPLED_SOLVE_SCHEMA = "scdat.internal_response_coupled_solve.v1"
INTERNAL_RESPONSE_COUPLED_SOLVE_CONTRACT = "internal-response-coupled-solve-v1"


def load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected top-level object in {path}")
    return payload


def run_command(command: List[str], cwd: Path) -> Tuple[int, str, str]:
    completed = subprocess.run(
        command,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    return completed.returncode, completed.stdout, completed.stderr


def make_smoke_payload(module_case: Dict[str, Any], payload: Dict[str, Any], output_csv: Path) -> Dict[str, Any]:
    smoke = copy.deepcopy(payload)

    smoke["schema_version"] = "v1"
    smoke["module"] = module_case["module"]

    run_block = smoke.setdefault("run", {})
    if not isinstance(run_block, dict):
        run_block = {}
        smoke["run"] = run_block

    run_block["steps"] = 1
    run_block.setdefault("time_step_s", module_case["default_dt"])
    run_block["output_csv"] = output_csv.as_posix()

    if module_case["module"] == "surface":
        run_block["adaptive_time_stepping"] = False
        run_block["total_duration_s"] = 0.0

    smoke["output_csv"] = output_csv.as_posix()
    return smoke


def metadata_sidecar_path(artifact_path: Path) -> Path:
    return Path(str(artifact_path) + ".metadata.json")


def surface_boundary_observer_path(metadata: Dict[str, Any], artifact_path: Path) -> Path:
    artifact_name = str(metadata.get("surface_pic_runtime_boundary_observer_artifact", "")).strip()
    if artifact_name:
        return artifact_path.parent / artifact_name
    return artifact_path.with_name(artifact_path.stem + ".shared_runtime_observer.json")


def parse_positive_int(text: str, key: str, errors: List[str]) -> int:
    try:
        value = int(text)
    except (TypeError, ValueError):
        errors.append(f"invalid_integer:{key}={text}")
        return -1

    if value <= 0:
        errors.append(f"non_positive_integer:{key}={value}")
    return value


def normalize_path(text: str) -> Path:
    return Path(text).resolve()


def validate_sidecar(
    artifact_format: str,
    artifact_path: Path,
    expected_artifacts: Dict[str, Path],
) -> Tuple[Dict[str, Any], List[str]]:
    errors: List[str] = []
    sidecar_path = metadata_sidecar_path(artifact_path)

    if not sidecar_path.is_file():
        return {}, [f"missing_sidecar:{sidecar_path}"]

    try:
        payload = load_json(sidecar_path)
    except Exception as exc:  # pylint: disable=broad-except
        return {}, [f"invalid_sidecar_json:{sidecar_path}:{exc}"]

    missing_top_level = sorted(REQUIRED_TOP_LEVEL_KEYS - payload.keys())
    for key in missing_top_level:
        errors.append(f"missing_top_level_key:{sidecar_path}:{key}")

    if payload.get("artifact_format") != artifact_format:
        errors.append(
            f"artifact_format_mismatch:{sidecar_path}:expected={artifact_format}:actual={payload.get('artifact_format')}"
        )

    if payload.get("output_contract_schema") != "scdat.output.contract.v1":
        errors.append(
            f"invalid_contract_schema:{sidecar_path}:{payload.get('output_contract_schema')}"
        )

    if payload.get("output_contract_version") != "v1":
        errors.append(
            f"invalid_contract_version:{sidecar_path}:{payload.get('output_contract_version')}"
        )

    artifact_path_text = payload.get("artifact_path", "")
    try:
        actual_artifact_path = normalize_path(str(artifact_path_text))
    except Exception as exc:  # pylint: disable=broad-except
        errors.append(f"invalid_artifact_path:{sidecar_path}:{exc}")
        actual_artifact_path = Path("<invalid>")

    if actual_artifact_path != artifact_path.resolve():
        errors.append(
            f"artifact_path_mismatch:{sidecar_path}:expected={artifact_path.resolve()}:actual={actual_artifact_path}"
        )

    metadata = payload.get("metadata")
    if not isinstance(metadata, dict):
        errors.append(f"missing_or_invalid_metadata_block:{sidecar_path}")
        return {}, errors

    missing_metadata = sorted(REQUIRED_METADATA_KEYS - metadata.keys())
    for key in missing_metadata:
        errors.append(f"missing_metadata_key:{sidecar_path}:{key}")

    module_label = str(metadata.get("module", "")).strip()
    if not module_label:
        errors.append(f"empty_module_label:{sidecar_path}")

    if metadata.get("output_contract_schema") != "scdat.output.contract.v1":
        errors.append(
            f"metadata_contract_schema_mismatch:{sidecar_path}:{metadata.get('output_contract_schema')}"
        )

    if metadata.get("output_contract_version") != "v1":
        errors.append(
            f"metadata_contract_version_mismatch:{sidecar_path}:{metadata.get('output_contract_version')}"
        )

    axis_name = str(metadata.get("output_axis_name", "")).strip()
    if not axis_name:
        errors.append(f"empty_output_axis_name:{sidecar_path}")

    parse_positive_int(str(metadata.get("output_row_count", "")), "output_row_count", errors)
    parse_positive_int(str(metadata.get("output_series_count", "")), "output_series_count", errors)

    series_names = str(metadata.get("output_series_names", "")).strip()
    if not series_names:
        errors.append(f"empty_output_series_names:{sidecar_path}")

    for key, expected_path in (
        ("output_artifact_csv", expected_artifacts["csv"]),
        ("output_artifact_vtk", expected_artifacts["vtk"]),
        ("output_artifact_h5", expected_artifacts["h5"]),
    ):
        value = metadata.get(key)
        if value is None:
            continue
        try:
            normalized = normalize_path(str(value))
        except Exception as exc:  # pylint: disable=broad-except
            errors.append(f"invalid_metadata_path:{sidecar_path}:{key}:{exc}")
            continue
        if normalized != expected_path.resolve():
            errors.append(
                f"metadata_path_mismatch:{sidecar_path}:{key}:expected={expected_path.resolve()}:actual={normalized}"
            )

    return metadata, errors


def validate_surface_boundary_observer(metadata: Dict[str, Any], artifact_path: Path) -> List[str]:
    errors: List[str] = []
    contract_id = str(metadata.get("surface_pic_runtime_boundary_observer_contract_id", "")).strip()
    if not contract_id:
        return errors

    if contract_id != SURFACE_BOUNDARY_OBSERVER_CONTRACT:
        errors.append(
            f"invalid_surface_boundary_observer_contract:{artifact_path}:{contract_id}"
        )
        return errors

    observer_path = surface_boundary_observer_path(metadata, artifact_path)
    if not observer_path.is_file():
        errors.append(f"missing_surface_boundary_observer:{observer_path}")
        return errors

    try:
        observer_payload = load_json(observer_path)
    except Exception as exc:  # pylint: disable=broad-except
        errors.append(f"invalid_surface_boundary_observer_json:{observer_path}:{exc}")
        return errors

    actual_schema = str(observer_payload.get("schema_version", "")).strip()
    if actual_schema != SURFACE_BOUNDARY_OBSERVER_SCHEMA:
        errors.append(
            f"surface_boundary_observer_schema_mismatch:{observer_path}:expected={SURFACE_BOUNDARY_OBSERVER_SCHEMA}:actual={actual_schema}"
        )

    actual_contract = str(observer_payload.get("contract_id", "")).strip()
    if actual_contract != SURFACE_BOUNDARY_OBSERVER_CONTRACT:
        errors.append(
            f"surface_boundary_observer_payload_contract_mismatch:{observer_path}:expected={SURFACE_BOUNDARY_OBSERVER_CONTRACT}:actual={actual_contract}"
        )

    module_runtime = str(metadata.get("surface_pic_runtime", "")).strip()
    observer_runtime = str(observer_payload.get("surface_pic_runtime", "")).strip()
    if module_runtime and observer_runtime != module_runtime:
        errors.append(
            f"surface_boundary_observer_runtime_mismatch:{observer_path}:expected={module_runtime}:actual={observer_runtime}"
        )

    if not isinstance(observer_payload.get("shared_runtime_enabled"), bool):
        errors.append(f"surface_boundary_observer_missing_bool_enabled:{observer_path}")

    nodes = observer_payload.get("nodes")
    if not isinstance(nodes, list) or not nodes:
        errors.append(f"surface_boundary_observer_missing_nodes:{observer_path}")

    return errors


def validate_surface_runtime_consistency(metadata: Dict[str, Any], artifact_path: Path) -> List[str]:
    errors: List[str] = []
    contract_id = str(metadata.get("surface_pic_runtime_consistency_contract_id", "")).strip()
    if not contract_id:
        return errors

    if contract_id != SURFACE_RUNTIME_CONSISTENCY_CONTRACT:
        errors.append(
            f"invalid_surface_runtime_consistency_contract:{artifact_path}:{contract_id}"
        )
        return errors

    artifact_name = str(metadata.get("surface_pic_runtime_consistency_artifact", "")).strip()
    consistency_path = (
        artifact_path.parent / artifact_name
        if artifact_name
        else artifact_path.with_name(artifact_path.stem + ".shared_runtime_consistency.json")
    )
    if not consistency_path.is_file():
        errors.append(f"missing_surface_runtime_consistency:{consistency_path}")
        return errors

    payload = load_json(consistency_path)
    if payload.get("schema_version") != SURFACE_RUNTIME_CONSISTENCY_SCHEMA:
        errors.append(
            f"surface_runtime_consistency_schema_mismatch:{consistency_path}:actual={payload.get('schema_version')}"
        )
    if payload.get("contract_id") != SURFACE_RUNTIME_CONSISTENCY_CONTRACT:
        errors.append(
            f"surface_runtime_consistency_payload_contract_mismatch:{consistency_path}:actual={payload.get('contract_id')}"
        )
    if payload.get("consistency_pass") is not True:
        errors.append(f"surface_runtime_consistency_not_passed:{consistency_path}")
    if payload.get("sheath_consistency_pass") is not True:
        errors.append(f"surface_runtime_sheath_consistency_not_passed:{consistency_path}")
    if payload.get("shared_particle_pool_consistency_pass") is not True:
        errors.append(
            f"surface_runtime_shared_particle_pool_consistency_not_passed:{consistency_path}"
        )
    shared_patch_node_count = payload.get("shared_patch_node_count")
    if not isinstance(shared_patch_node_count, int) or shared_patch_node_count <= 0:
        errors.append(
            f"surface_runtime_consistency_invalid_shared_patch_node_count:{consistency_path}"
        )
    patch_potential_spread_v = payload.get("patch_potential_spread_v")
    if not isinstance(patch_potential_spread_v, (int, float)):
        errors.append(
            f"surface_runtime_consistency_invalid_patch_potential_spread:{consistency_path}"
        )
    pre_patch_spread_v = payload.get("pre_global_sheath_proxy_patch_potential_spread_v")
    if not isinstance(pre_patch_spread_v, (int, float)):
        errors.append(
            f"surface_runtime_consistency_invalid_pre_global_sheath_proxy_patch_spread:{consistency_path}"
        )
    patch_spread_reduction_v = payload.get("patch_potential_spread_reduction_v")
    if not isinstance(patch_spread_reduction_v, (int, float)):
        errors.append(
            f"surface_runtime_consistency_invalid_patch_spread_reduction_v:{consistency_path}"
        )
    patch_spread_reduction_ratio = payload.get("patch_potential_spread_reduction_ratio")
    if not isinstance(patch_spread_reduction_ratio, (int, float)):
        errors.append(
            f"surface_runtime_consistency_invalid_patch_spread_reduction_ratio:{consistency_path}"
        )
    elif float(patch_spread_reduction_ratio) < 0.0 or float(patch_spread_reduction_ratio) > 1.0:
        errors.append(
            f"surface_runtime_consistency_patch_spread_reduction_ratio_out_of_range:{consistency_path}:actual={patch_spread_reduction_ratio}"
        )
    global_solve_weight = payload.get("global_sheath_proxy_solve_weight")
    if not isinstance(global_solve_weight, (int, float)):
        errors.append(
            f"surface_runtime_consistency_invalid_global_sheath_proxy_solve_weight:{consistency_path}"
        )
    global_solve_active = payload.get("global_sheath_proxy_solve_active")
    if not isinstance(global_solve_active, bool):
        errors.append(
            f"surface_runtime_consistency_invalid_global_sheath_proxy_solve_active:{consistency_path}"
        )
    shared_current_matrix_coupling_active = payload.get(
        "shared_current_matrix_coupling_active"
    )
    if not isinstance(shared_current_matrix_coupling_active, bool):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_current_matrix_coupling_active:{consistency_path}"
        )
    shared_current_matrix_coupling_offdiag_entry_count = payload.get(
        "shared_current_matrix_coupling_offdiag_entry_count"
    )
    if not isinstance(shared_current_matrix_coupling_offdiag_entry_count, int) or (
        shared_current_matrix_coupling_offdiag_entry_count < 0
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_current_matrix_coupling_offdiag_entry_count:{consistency_path}"
        )
    shared_global_coupled_solve_active = payload.get(
        "shared_global_coupled_solve_active"
    )
    if not isinstance(shared_global_coupled_solve_active, bool):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_global_coupled_solve_active:{consistency_path}"
        )
    shared_global_coupled_solve_iterations = payload.get(
        "shared_global_coupled_solve_iterations"
    )
    if not isinstance(shared_global_coupled_solve_iterations, int) or (
        shared_global_coupled_solve_iterations < 0
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_global_coupled_solve_iterations:{consistency_path}"
        )
    shared_global_coupled_solve_converged = payload.get(
        "shared_global_coupled_solve_converged"
    )
    if not isinstance(shared_global_coupled_solve_converged, bool):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_global_coupled_solve_converged:{consistency_path}"
        )
    shared_global_coupled_solve_max_delta_v = payload.get(
        "shared_global_coupled_solve_max_delta_v"
    )
    if not isinstance(shared_global_coupled_solve_max_delta_v, (int, float)) or (
        float(shared_global_coupled_solve_max_delta_v) < 0.0
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_global_coupled_solve_max_delta_v:{consistency_path}"
        )
    normal_electric_field_spread_v_per_m = payload.get(
        "normal_electric_field_spread_v_per_m"
    )
    if not isinstance(normal_electric_field_spread_v_per_m, (int, float)) or (
        float(normal_electric_field_spread_v_per_m) < 0.0
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_normal_electric_field_spread:{consistency_path}"
        )
    local_charge_density_spread_c_per_m3 = payload.get(
        "local_charge_density_spread_c_per_m3"
    )
    if not isinstance(local_charge_density_spread_c_per_m3, (int, float)) or (
        float(local_charge_density_spread_c_per_m3) < 0.0
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_local_charge_density_spread:{consistency_path}"
        )
    shared_live_pic_coupled_refresh_active = payload.get(
        "shared_live_pic_coupled_refresh_active"
    )
    if not isinstance(shared_live_pic_coupled_refresh_active, bool):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_live_pic_coupled_refresh_active:{consistency_path}"
        )
    shared_live_pic_coupled_refresh_count = payload.get(
        "shared_live_pic_coupled_refresh_count"
    )
    if not isinstance(shared_live_pic_coupled_refresh_count, int) or (
        shared_live_pic_coupled_refresh_count < 0
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_live_pic_coupled_refresh_count:{consistency_path}"
        )
    shared_particle_transport_coupling_active = payload.get(
        "shared_particle_transport_coupling_active"
    )
    if not isinstance(shared_particle_transport_coupling_active, bool):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_coupling_active:{consistency_path}"
        )
    shared_particle_transport_offdiag_entry_count = payload.get(
        "shared_particle_transport_offdiag_entry_count"
    )
    if not isinstance(shared_particle_transport_offdiag_entry_count, int) or (
        shared_particle_transport_offdiag_entry_count < 0
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_offdiag_entry_count:{consistency_path}"
        )
    shared_particle_transport_total_conductance_s = payload.get(
        "shared_particle_transport_total_conductance_s"
    )
    if not isinstance(shared_particle_transport_total_conductance_s, (int, float)) or (
        float(shared_particle_transport_total_conductance_s) < 0.0
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_total_conductance_s:{consistency_path}"
        )
    shared_particle_transport_conservation_error_a_per_v = payload.get(
        "shared_particle_transport_conservation_error_a_per_v"
    )
    if not isinstance(
        shared_particle_transport_conservation_error_a_per_v, (int, float)
    ) or (float(shared_particle_transport_conservation_error_a_per_v) < 0.0):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_conservation_error_a_per_v:{consistency_path}"
        )
    shared_particle_transport_charge_c = payload.get(
        "shared_particle_transport_charge_c"
    )
    if not isinstance(shared_particle_transport_charge_c, (int, float)):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_charge_c:{consistency_path}"
        )
    shared_particle_transport_reference_shift_v = payload.get(
        "shared_particle_transport_reference_shift_v"
    )
    if not isinstance(shared_particle_transport_reference_shift_v, (int, float)):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_reference_shift_v:{consistency_path}"
        )
    shared_particle_transport_distribution_active = payload.get(
        "shared_particle_transport_distribution_active"
    )
    if not isinstance(shared_particle_transport_distribution_active, bool):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_distribution_active:{consistency_path}"
        )
    shared_particle_transport_distribution_charge_spread_c = payload.get(
        "shared_particle_transport_distribution_charge_spread_c"
    )
    if not isinstance(
        shared_particle_transport_distribution_charge_spread_c, (int, float)
    ) or (float(shared_particle_transport_distribution_charge_spread_c) < 0.0):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_distribution_charge_spread_c:{consistency_path}"
        )
    shared_particle_transport_distribution_reference_shift_spread_v = payload.get(
        "shared_particle_transport_distribution_reference_shift_spread_v"
    )
    if not isinstance(
        shared_particle_transport_distribution_reference_shift_spread_v, (int, float)
    ) or (float(shared_particle_transport_distribution_reference_shift_spread_v) < 0.0):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_distribution_reference_shift_spread_v:{consistency_path}"
        )
    shared_particle_transport_distribution_conservation_error_c = payload.get(
        "shared_particle_transport_distribution_conservation_error_c"
    )
    if not isinstance(
        shared_particle_transport_distribution_conservation_error_c, (int, float)
    ) or (float(shared_particle_transport_distribution_conservation_error_c) < 0.0):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_distribution_conservation_error_c:{consistency_path}"
        )
    shared_particle_transport_exchange_active = payload.get(
        "shared_particle_transport_exchange_active"
    )
    if not isinstance(shared_particle_transport_exchange_active, bool):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_exchange_active:{consistency_path}"
        )
    shared_particle_transport_exchange_flux_spread_a = payload.get(
        "shared_particle_transport_exchange_flux_spread_a"
    )
    if not isinstance(
        shared_particle_transport_exchange_flux_spread_a, (int, float)
    ) or (float(shared_particle_transport_exchange_flux_spread_a) < 0.0):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_exchange_flux_spread_a:{consistency_path}"
        )
    shared_particle_transport_exchange_flux_conservation_error_a = payload.get(
        "shared_particle_transport_exchange_flux_conservation_error_a"
    )
    if not isinstance(
        shared_particle_transport_exchange_flux_conservation_error_a, (int, float)
    ) or (float(shared_particle_transport_exchange_flux_conservation_error_a) < 0.0):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_exchange_flux_conservation_error_a:{consistency_path}"
        )
    shared_particle_transport_edge_domain_active = payload.get(
        "shared_particle_transport_edge_domain_active"
    )
    if not isinstance(shared_particle_transport_edge_domain_active, bool):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_domain_active:{consistency_path}"
        )
    shared_particle_transport_edge_charge_total_abs_c = payload.get(
        "shared_particle_transport_edge_charge_total_abs_c"
    )
    if not isinstance(shared_particle_transport_edge_charge_total_abs_c, (int, float)) or (
        float(shared_particle_transport_edge_charge_total_abs_c) < 0.0
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_charge_total_abs_c:{consistency_path}"
        )
    shared_particle_transport_edge_charge_spread_c = payload.get(
        "shared_particle_transport_edge_charge_spread_c"
    )
    if not isinstance(shared_particle_transport_edge_charge_spread_c, (int, float)) or (
        float(shared_particle_transport_edge_charge_spread_c) < 0.0
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_charge_spread_c:{consistency_path}"
        )
    shared_particle_transport_edge_charge_conservation_error_c = payload.get(
        "shared_particle_transport_edge_charge_conservation_error_c"
    )
    if not isinstance(
        shared_particle_transport_edge_charge_conservation_error_c, (int, float)
    ) or (float(shared_particle_transport_edge_charge_conservation_error_c) < 0.0):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_charge_conservation_error_c:{consistency_path}"
        )
    shared_particle_transport_edge_operator_active = payload.get(
        "shared_particle_transport_edge_operator_active"
    )
    if not isinstance(shared_particle_transport_edge_operator_active, bool):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_operator_active:{consistency_path}"
        )
    shared_particle_transport_edge_target_charge_total_abs_c = payload.get(
        "shared_particle_transport_edge_target_charge_total_abs_c"
    )
    if not isinstance(
        shared_particle_transport_edge_target_charge_total_abs_c, (int, float)
    ) or (float(shared_particle_transport_edge_target_charge_total_abs_c) < 0.0):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_target_charge_total_abs_c:{consistency_path}"
        )
    shared_particle_transport_edge_operator_total_abs_drive_charge_c = payload.get(
        "shared_particle_transport_edge_operator_total_abs_drive_charge_c"
    )
    if not isinstance(
        shared_particle_transport_edge_operator_total_abs_drive_charge_c, (int, float)
    ) or (float(shared_particle_transport_edge_operator_total_abs_drive_charge_c) < 0.0):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_operator_total_abs_drive_charge_c:{consistency_path}"
        )
    shared_particle_transport_edge_operator_drive_conservation_error_c = payload.get(
        "shared_particle_transport_edge_operator_drive_conservation_error_c"
    )
    if not isinstance(
        shared_particle_transport_edge_operator_drive_conservation_error_c, (int, float)
    ) or (float(shared_particle_transport_edge_operator_drive_conservation_error_c) < 0.0):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_operator_drive_conservation_error_c:{consistency_path}"
        )
    shared_particle_transport_edge_graph_operator_iterations = payload.get(
        "shared_particle_transport_edge_graph_operator_iterations"
    )
    if not isinstance(shared_particle_transport_edge_graph_operator_iterations, (int, float)) or (
        float(shared_particle_transport_edge_graph_operator_iterations) < 0.0
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_graph_operator_iterations:{consistency_path}"
        )
    shared_particle_transport_edge_graph_operator_converged = payload.get(
        "shared_particle_transport_edge_graph_operator_converged"
    )
    if not isinstance(shared_particle_transport_edge_graph_operator_converged, bool):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_graph_operator_converged:{consistency_path}"
        )
    shared_particle_transport_edge_graph_operator_max_balance_residual_c = payload.get(
        "shared_particle_transport_edge_graph_operator_max_balance_residual_c"
    )
    if not isinstance(
        shared_particle_transport_edge_graph_operator_max_balance_residual_c, (int, float)
    ) or (float(shared_particle_transport_edge_graph_operator_max_balance_residual_c) < 0.0):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_graph_operator_max_balance_residual_c:{consistency_path}"
        )
    shared_particle_transport_edge_graph_operator_branch_graph_active = payload.get(
        "shared_particle_transport_edge_graph_operator_branch_graph_active"
    )
    if not isinstance(
        shared_particle_transport_edge_graph_operator_branch_graph_active, bool
    ):
        errors.append(
            f"surface_runtime_consistency_invalid_shared_particle_transport_edge_graph_operator_branch_graph_active:{consistency_path}"
        )
    for key in (
        "shared_particle_transport_edge_graph_operator_branch_graph_edge_count",
        "shared_particle_transport_edge_graph_operator_branch_graph_pair_count",
        "shared_particle_transport_edge_graph_operator_effective_pair_count",
        "shared_particle_transport_edge_graph_operator_total_pair_weight_f",
        "shared_particle_transport_edge_graph_operator_total_conductance_weight_f",
        "shared_particle_transport_edge_graph_operator_min_node_preconditioner",
        "shared_particle_transport_edge_graph_operator_max_node_preconditioner",
    ):
        value = payload.get(key)
        if not isinstance(value, (int, float)) or float(value) < 0.0:
            errors.append(
                f"surface_runtime_consistency_invalid_{key}:{consistency_path}"
            )
    return errors


def validate_surface_particle_transport_domain(
    metadata: Dict[str, Any], artifact_path: Path
) -> List[str]:
    errors: List[str] = []
    contract_id = str(
        metadata.get("surface_pic_runtime_shared_particle_transport_domain_contract_id", "")
    ).strip()
    if not contract_id:
        return errors

    if contract_id != SURFACE_PARTICLE_TRANSPORT_DOMAIN_CONTRACT:
        errors.append(
            f"invalid_surface_particle_transport_domain_contract:{artifact_path}:{contract_id}"
        )
        return errors

    artifact_name = str(
        metadata.get("surface_pic_runtime_shared_particle_transport_domain_artifact", "")
    ).strip()
    domain_path = (
        artifact_path.parent / artifact_name
        if artifact_name
        else artifact_path.with_name(
            artifact_path.stem + ".shared_particle_transport_domain.json"
        )
    )
    if not domain_path.is_file():
        errors.append(f"missing_surface_particle_transport_domain:{domain_path}")
        return errors

    payload = load_json(domain_path)
    if payload.get("schema_version") != SURFACE_PARTICLE_TRANSPORT_DOMAIN_SCHEMA:
        errors.append(
            f"surface_particle_transport_domain_schema_mismatch:{domain_path}:actual={payload.get('schema_version')}"
        )
    if payload.get("contract_id") != SURFACE_PARTICLE_TRANSPORT_DOMAIN_CONTRACT:
        errors.append(
            f"surface_particle_transport_domain_payload_contract_mismatch:{domain_path}:actual={payload.get('contract_id')}"
        )

    domain_active = payload.get("shared_particle_transport_domain_active")
    if not isinstance(domain_active, bool):
        errors.append(
            f"surface_particle_transport_domain_invalid_active:{domain_path}"
        )
    shared_patch_node_count = payload.get("shared_patch_node_count")
    if not isinstance(shared_patch_node_count, int) or shared_patch_node_count <= 0:
        errors.append(
            f"surface_particle_transport_domain_invalid_shared_patch_node_count:{domain_path}"
        )
    exchange_edge_count = payload.get("exchange_edge_count")
    if not isinstance(exchange_edge_count, int) or exchange_edge_count < 0:
        errors.append(
            f"surface_particle_transport_domain_invalid_exchange_edge_count:{domain_path}"
        )
    shared_particle_transport_charge_c = payload.get("shared_particle_transport_charge_c")
    if not isinstance(shared_particle_transport_charge_c, (int, float)):
        errors.append(
            f"surface_particle_transport_domain_invalid_shared_particle_transport_charge_c:{domain_path}"
        )
    shared_particle_transport_reference_shift_v = payload.get(
        "shared_particle_transport_reference_shift_v"
    )
    if not isinstance(shared_particle_transport_reference_shift_v, (int, float)):
        errors.append(
            f"surface_particle_transport_domain_invalid_shared_particle_transport_reference_shift_v:{domain_path}"
        )
    distributed_transport_charge_total_c = payload.get(
        "distributed_particle_transport_charge_total_c"
    )
    if not isinstance(distributed_transport_charge_total_c, (int, float)):
        errors.append(
            f"surface_particle_transport_domain_invalid_distributed_particle_transport_charge_total_c:{domain_path}"
        )
    distributed_transport_charge_conservation_error_c = payload.get(
        "distributed_particle_transport_charge_conservation_error_c"
    )
    if not isinstance(
        distributed_transport_charge_conservation_error_c, (int, float)
    ) or (float(distributed_transport_charge_conservation_error_c) < 0.0):
        errors.append(
            f"surface_particle_transport_domain_invalid_distributed_particle_transport_charge_conservation_error_c:{domain_path}"
        )
    exchange_flux_conservation_error_a = payload.get(
        "exchange_flux_conservation_error_a"
    )
    if not isinstance(exchange_flux_conservation_error_a, (int, float)) or (
        float(exchange_flux_conservation_error_a) < 0.0
    ):
        errors.append(
            f"surface_particle_transport_domain_invalid_exchange_flux_conservation_error_a:{domain_path}"
        )
    edge_charge_total_abs_c = payload.get("edge_charge_total_abs_c")
    if not isinstance(edge_charge_total_abs_c, (int, float)) or (
        float(edge_charge_total_abs_c) < 0.0
    ):
        errors.append(
            f"surface_particle_transport_domain_invalid_edge_charge_total_abs_c:{domain_path}"
        )
    edge_charge_conservation_error_c = payload.get("edge_charge_conservation_error_c")
    if not isinstance(edge_charge_conservation_error_c, (int, float)) or (
        float(edge_charge_conservation_error_c) < 0.0
    ):
        errors.append(
            f"surface_particle_transport_domain_invalid_edge_charge_conservation_error_c:{domain_path}"
        )
    edge_target_charge_total_abs_c = payload.get("edge_target_charge_total_abs_c")
    if not isinstance(edge_target_charge_total_abs_c, (int, float)) or (
        float(edge_target_charge_total_abs_c) < 0.0
    ):
        errors.append(
            f"surface_particle_transport_domain_invalid_edge_target_charge_total_abs_c:{domain_path}"
        )
    edge_operator_drive_total_abs_c = payload.get("edge_operator_drive_total_abs_c")
    if not isinstance(edge_operator_drive_total_abs_c, (int, float)) or (
        float(edge_operator_drive_total_abs_c) < 0.0
    ):
        errors.append(
            f"surface_particle_transport_domain_invalid_edge_operator_drive_total_abs_c:{domain_path}"
        )
    edge_operator_drive_conservation_error_c = payload.get(
        "edge_operator_drive_conservation_error_c"
    )
    if not isinstance(edge_operator_drive_conservation_error_c, (int, float)) or (
        float(edge_operator_drive_conservation_error_c) < 0.0
    ):
        errors.append(
            f"surface_particle_transport_domain_invalid_edge_operator_drive_conservation_error_c:{domain_path}"
        )
    edge_graph_operator_iterations = payload.get("edge_graph_operator_iterations")
    if not isinstance(edge_graph_operator_iterations, (int, float)) or (
        float(edge_graph_operator_iterations) < 0.0
    ):
        errors.append(
            f"surface_particle_transport_domain_invalid_edge_graph_operator_iterations:{domain_path}"
        )
    edge_graph_operator_converged = payload.get("edge_graph_operator_converged")
    if not isinstance(edge_graph_operator_converged, bool):
        errors.append(
            f"surface_particle_transport_domain_invalid_edge_graph_operator_converged:{domain_path}"
        )
    edge_graph_operator_max_balance_residual_c = payload.get(
        "edge_graph_operator_max_balance_residual_c"
    )
    if not isinstance(edge_graph_operator_max_balance_residual_c, (int, float)) or (
        float(edge_graph_operator_max_balance_residual_c) < 0.0
    ):
        errors.append(
            f"surface_particle_transport_domain_invalid_edge_graph_operator_max_balance_residual_c:{domain_path}"
        )
    edge_graph_operator_branch_graph_active = payload.get(
        "edge_graph_operator_branch_graph_active"
    )
    if not isinstance(edge_graph_operator_branch_graph_active, bool):
        errors.append(
            f"surface_particle_transport_domain_invalid_edge_graph_operator_branch_graph_active:{domain_path}"
        )
    for key in (
        "edge_graph_operator_branch_graph_edge_count",
        "edge_graph_operator_branch_graph_pair_count",
        "edge_graph_operator_effective_pair_count",
        "edge_graph_operator_total_pair_weight_f",
        "edge_graph_operator_total_conductance_weight_f",
        "edge_graph_operator_min_node_preconditioner",
        "edge_graph_operator_max_node_preconditioner",
    ):
        value = payload.get(key)
        if not isinstance(value, (int, float)) or float(value) < 0.0:
            errors.append(
                f"surface_particle_transport_domain_invalid_{key}:{domain_path}"
            )

    nodes = payload.get("nodes")
    if not isinstance(nodes, list) or not nodes:
        errors.append(f"surface_particle_transport_domain_missing_nodes:{domain_path}")
    elif isinstance(shared_patch_node_count, int) and len(nodes) != shared_patch_node_count:
        errors.append(
            f"surface_particle_transport_domain_node_count_mismatch:{domain_path}:expected={shared_patch_node_count}:actual={len(nodes)}"
        )
    else:
        for node in nodes:
            if not isinstance(node, dict):
                errors.append(
                    f"surface_particle_transport_domain_invalid_node_entry:{domain_path}"
                )
                break
            for key in (
                "index",
                "name",
                "latest_patch_potential_v",
                "latest_shared_reference_potential_v",
                "latest_distributed_particle_transport_charge_c",
                "latest_distributed_particle_transport_reference_shift_v",
                "latest_distributed_particle_transport_net_flux_a",
                "net_edge_stored_charge_c",
                "net_edge_target_charge_c",
                "net_edge_operator_drive_charge_c",
            ):
                if key not in node:
                    errors.append(
                        f"surface_particle_transport_domain_missing_node_field:{domain_path}:{key}"
                    )
                    break

    edges = payload.get("exchange_edges")
    if not isinstance(edges, list):
        errors.append(
            f"surface_particle_transport_domain_missing_exchange_edges:{domain_path}"
        )
    elif isinstance(exchange_edge_count, int) and len(edges) != exchange_edge_count:
        errors.append(
            f"surface_particle_transport_domain_exchange_edge_count_mismatch:{domain_path}:expected={exchange_edge_count}:actual={len(edges)}"
        )
    else:
        for edge in edges:
            if not isinstance(edge, dict):
                errors.append(
                    f"surface_particle_transport_domain_invalid_exchange_edge:{domain_path}"
                )
                break
            for key in (
                "from_index",
                "from_name",
                "to_index",
                "to_name",
                "net_flux_a",
                "abs_flux_a",
                "stored_charge_c",
                "abs_stored_charge_c",
                "target_charge_c",
                "abs_target_charge_c",
                "operator_drive_charge_c",
                "abs_operator_drive_charge_c",
            ):
                if key not in edge:
                    errors.append(
                        f"surface_particle_transport_domain_missing_exchange_edge_field:{domain_path}:{key}"
                    )
                    break
    return errors


def validate_surface_global_particle_domain(
    metadata: Dict[str, Any], artifact_path: Path
) -> List[str]:
    errors: List[str] = []
    contract_id = str(
        metadata.get("surface_pic_runtime_global_particle_domain_contract_id", "")
    ).strip()
    if not contract_id:
        return errors
    if contract_id != SURFACE_GLOBAL_PARTICLE_DOMAIN_CONTRACT:
        errors.append(
            f"invalid_surface_global_particle_domain_contract:{artifact_path}:{contract_id}"
        )
        return errors

    artifact_name = str(
        metadata.get("surface_pic_runtime_global_particle_domain_artifact", "")
    ).strip()
    domain_path = (
        artifact_path.parent / artifact_name
        if artifact_name
        else artifact_path.with_name(artifact_path.stem + ".global_particle_domain.json")
    )
    if not domain_path.is_file():
        errors.append(f"missing_surface_global_particle_domain:{domain_path}")
        return errors

    payload = load_json(domain_path)
    if payload.get("schema_version") != SURFACE_GLOBAL_PARTICLE_DOMAIN_SCHEMA:
        errors.append(
            f"surface_global_particle_domain_schema_mismatch:{domain_path}:actual={payload.get('schema_version')}"
        )
    if payload.get("contract_id") != SURFACE_GLOBAL_PARTICLE_DOMAIN_CONTRACT:
        errors.append(
            f"surface_global_particle_domain_payload_contract_mismatch:{domain_path}:actual={payload.get('contract_id')}"
        )
    if not isinstance(payload.get("global_particle_domain_active"), bool):
        errors.append(
            f"surface_global_particle_domain_invalid_active:{domain_path}"
        )
    for key in ("shared_patch_node_count", "domain_edge_count"):
        value = payload.get(key)
        if not isinstance(value, int) or value < 0:
            errors.append(f"surface_global_particle_domain_invalid_{key}:{domain_path}")
    for key in (
        "global_particle_charge_c",
        "global_particle_reference_shift_v",
        "global_particle_charge_conservation_error_c",
        "global_particle_flux_conservation_error_a",
        "global_particle_edge_charge_total_abs_c",
        "global_particle_edge_target_charge_total_abs_c",
        "global_particle_edge_operator_drive_total_abs_c",
        "edge_graph_operator_iterations",
        "edge_graph_operator_max_balance_residual_c",
        "edge_graph_operator_branch_graph_edge_count",
        "edge_graph_operator_branch_graph_pair_count",
        "edge_graph_operator_effective_pair_count",
        "edge_graph_operator_total_pair_weight_f",
        "edge_graph_operator_total_conductance_weight_f",
        "edge_graph_operator_min_node_preconditioner",
        "edge_graph_operator_max_node_preconditioner",
    ):
        value = payload.get(key)
        if not isinstance(value, (int, float)):
            errors.append(f"surface_global_particle_domain_invalid_{key}:{domain_path}")
        elif (
            "error" in key or "total_abs" in key or key.endswith("_count")
        ) and float(value) < 0.0:
            errors.append(f"surface_global_particle_domain_negative_{key}:{domain_path}")
    if not isinstance(payload.get("edge_graph_operator_converged"), bool):
        errors.append(
            f"surface_global_particle_domain_invalid_edge_graph_operator_converged:{domain_path}"
        )
    if not isinstance(payload.get("edge_graph_operator_branch_graph_active"), bool):
        errors.append(
            f"surface_global_particle_domain_invalid_edge_graph_operator_branch_graph_active:{domain_path}"
        )
    nodes = payload.get("nodes")
    edges = payload.get("domain_edges")
    if not isinstance(nodes, list) or not nodes:
        errors.append(f"surface_global_particle_domain_missing_nodes:{domain_path}")
    if not isinstance(edges, list):
        errors.append(f"surface_global_particle_domain_missing_domain_edges:{domain_path}")
    return errors


def validate_surface_global_particle_repository(
    metadata: Dict[str, Any], artifact_path: Path
) -> List[str]:
    errors: List[str] = []
    contract_id = str(
        metadata.get("surface_pic_runtime_global_particle_repository_contract_id", "")
    ).strip()
    if not contract_id:
        return errors
    if contract_id != SURFACE_GLOBAL_PARTICLE_REPOSITORY_CONTRACT:
        errors.append(
            f"invalid_surface_global_particle_repository_contract:{artifact_path}:{contract_id}"
        )
        return errors

    artifact_name = str(
        metadata.get("surface_pic_runtime_global_particle_repository_artifact", "")
    ).strip()
    repository_path = (
        artifact_path.parent / artifact_name
        if artifact_name
        else artifact_path.with_name(
            artifact_path.stem + ".global_particle_repository.json"
        )
    )
    if not repository_path.is_file():
        errors.append(f"missing_surface_global_particle_repository:{repository_path}")
        return errors

    payload = load_json(repository_path)
    if payload.get("schema_version") != SURFACE_GLOBAL_PARTICLE_REPOSITORY_SCHEMA:
        errors.append(
            "surface_global_particle_repository_schema_mismatch:"
            f"{repository_path}:actual={payload.get('schema_version')}"
        )
    if payload.get("contract_id") != SURFACE_GLOBAL_PARTICLE_REPOSITORY_CONTRACT:
        errors.append(
            "surface_global_particle_repository_payload_contract_mismatch:"
            f"{repository_path}:actual={payload.get('contract_id')}"
        )

    if not isinstance(payload.get("runtime_state_backed"), bool):
        errors.append(
            f"surface_global_particle_repository_invalid_runtime_state_backed:{repository_path}"
        )
    if not isinstance(payload.get("global_particle_repository_active"), bool):
        errors.append(
            f"surface_global_particle_repository_invalid_active:{repository_path}"
        )

    bookkeeping_mode = payload.get("global_particle_repository_bookkeeping_mode")
    if not isinstance(bookkeeping_mode, str) or not bookkeeping_mode.strip():
        errors.append(
            f"surface_global_particle_repository_invalid_bookkeeping_mode:{repository_path}"
        )
    lifecycle_mode = payload.get("global_particle_repository_lifecycle_mode")
    if not isinstance(lifecycle_mode, str) or not lifecycle_mode.strip():
        errors.append(
            f"surface_global_particle_repository_invalid_lifecycle_mode:{repository_path}"
        )

    for key in ("shared_patch_node_count", "repository_edge_count"):
        value = payload.get(key)
        if not isinstance(value, int) or value < 0:
            errors.append(
                f"surface_global_particle_repository_invalid_{key}:{repository_path}"
            )

    signed_numeric_keys = (
        "total_reservoir_charge_c",
        "total_target_reservoir_charge_c",
    )
    nonnegative_numeric_keys = (
        "total_migration_delta_abs_charge_c",
        "total_edge_feedback_abs_charge_c",
        "total_conservation_correction_abs_charge_c",
        "total_migration_edge_abs_charge_c",
        "charge_conservation_error_c",
        "migration_charge_conservation_error_c",
    )
    for key in signed_numeric_keys + nonnegative_numeric_keys:
        value = payload.get(key)
        if not isinstance(value, (int, float)):
            errors.append(
                f"surface_global_particle_repository_invalid_{key}:{repository_path}"
            )
        elif key in nonnegative_numeric_keys and float(value) < 0.0:
            errors.append(
                f"surface_global_particle_repository_negative_{key}:{repository_path}"
            )

    nodes = payload.get("nodes")
    if not isinstance(nodes, list):
        errors.append(
            f"surface_global_particle_repository_missing_nodes:{repository_path}"
        )
    else:
        for node in nodes:
            if not isinstance(node, dict):
                errors.append(
                    f"surface_global_particle_repository_invalid_node:{repository_path}"
                )
                break
            for key in (
                "index",
                "name",
                "area_m2",
                "latest_patch_potential_v",
                "latest_shared_reference_potential_v",
                "reference_shift_v",
                "reservoir_charge_c",
                "target_reservoir_charge_c",
                "migration_delta_charge_c",
                "edge_feedback_charge_c",
                "conservation_correction_charge_c",
                "net_flux_a",
            ):
                if key not in node:
                    errors.append(
                        "surface_global_particle_repository_missing_node_field:"
                        f"{repository_path}:{key}"
                    )
                    break

    edges = payload.get("edges")
    if not isinstance(edges, list):
        errors.append(f"surface_global_particle_repository_missing_edges:{repository_path}")
    else:
        for edge in edges:
            if not isinstance(edge, dict):
                errors.append(
                    f"surface_global_particle_repository_invalid_edge:{repository_path}"
                )
                break
            for key in (
                "from_index",
                "from_name",
                "to_index",
                "to_name",
                "conductance_s",
                "migration_flux_a",
                "migration_charge_c",
                "stored_charge_c",
                "target_charge_c",
                "operator_drive_charge_c",
            ):
                if key not in edge:
                    errors.append(
                        "surface_global_particle_repository_missing_edge_field:"
                        f"{repository_path}:{key}"
                    )
                    break
    return errors


def validate_surface_global_sheath_field_solve(
    metadata: Dict[str, Any], artifact_path: Path
) -> List[str]:
    errors: List[str] = []
    contract_id = str(
        metadata.get("surface_pic_runtime_global_sheath_field_solve_contract_id", "")
    ).strip()
    if not contract_id:
        return errors
    if contract_id != SURFACE_GLOBAL_SHEATH_FIELD_SOLVE_CONTRACT:
        errors.append(
            f"invalid_surface_global_sheath_field_solve_contract:{artifact_path}:{contract_id}"
        )
        return errors

    artifact_name = str(
        metadata.get("surface_pic_runtime_global_sheath_field_solve_artifact", "")
    ).strip()
    solve_path = (
        artifact_path.parent / artifact_name
        if artifact_name
        else artifact_path.with_name(
            artifact_path.stem + ".global_sheath_field_solve.json"
        )
    )
    if not solve_path.is_file():
        errors.append(f"missing_surface_global_sheath_field_solve:{solve_path}")
        return errors

    payload = load_json(solve_path)
    if payload.get("schema_version") != SURFACE_GLOBAL_SHEATH_FIELD_SOLVE_SCHEMA:
        errors.append(
            f"surface_global_sheath_field_solve_schema_mismatch:{solve_path}:actual={payload.get('schema_version')}"
        )
    if payload.get("contract_id") != SURFACE_GLOBAL_SHEATH_FIELD_SOLVE_CONTRACT:
        errors.append(
            f"surface_global_sheath_field_solve_payload_contract_mismatch:{solve_path}:actual={payload.get('contract_id')}"
        )
    if not isinstance(payload.get("global_sheath_field_solve_active"), bool):
        errors.append(
            f"surface_global_sheath_field_solve_invalid_active:{solve_path}"
        )
    if not isinstance(payload.get("shared_global_coupled_solve_converged"), bool):
        errors.append(
            f"surface_global_sheath_field_solve_invalid_converged:{solve_path}"
        )
    for key in (
        "shared_patch_node_count",
        "global_reference_potential_v",
        "global_patch_potential_v",
        "global_effective_sheath_length_m",
        "global_normal_electric_field_v_per_m",
        "global_local_charge_density_c_per_m3",
        "shared_global_coupled_solve_iterations",
        "shared_global_coupled_solve_max_delta_v",
        "global_reference_spread_v",
        "global_patch_potential_spread_v",
        "global_normal_electric_field_spread_v_per_m",
        "global_local_charge_density_spread_c_per_m3",
        "field_residual_v_per_m",
        "particle_field_coupled_residual_v",
        "multi_step_stability_metric_v",
        "global_particle_charge_conservation_error_c",
        "global_particle_flux_conservation_error_a",
    ):
        value = payload.get(key)
        if not isinstance(value, (int, float)):
            errors.append(
                f"surface_global_sheath_field_solve_invalid_{key}:{solve_path}"
            )
        elif (
            "spread" in key
            or "residual" in key
            or "metric" in key
            or "error" in key
            or key.endswith("_iterations")
            or key.endswith("_length_m")
        ) and float(value) < 0.0:
            errors.append(
                f"surface_global_sheath_field_solve_negative_{key}:{solve_path}"
            )
    if not isinstance(payload.get("nodes"), list):
        errors.append(f"surface_global_sheath_field_solve_missing_nodes:{solve_path}")
    return errors


def resolve_metadata_artifact_path(metadata: Dict[str, Any], key: str, artifact_path: Path) -> Path | None:
    value = str(metadata.get(key, "")).strip()
    if not value:
        return None
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    return (artifact_path.parent.parent / candidate).resolve() if candidate.parts[:1] == ("build_codex",) else (artifact_path.parent / candidate).resolve()


def validate_radiation_history_artifacts(metadata: Dict[str, Any], artifact_path: Path) -> List[str]:
    errors: List[str] = []
    deposition_contract = str(metadata.get("deposition_record_contract_id", "")).strip()
    process_contract = str(metadata.get("process_history_contract_id", "")).strip()
    material_governance_contract = str(
        metadata.get("material_governance_contract_id", "")
    ).strip()
    benchmark_contract = str(
        metadata.get("radiation_transport_benchmark_contract_id", "")
    ).strip()
    deposition_path = resolve_metadata_artifact_path(
        metadata, "deposition_history_artifact_path", artifact_path
    )
    process_path = resolve_metadata_artifact_path(
        metadata, "process_history_artifact_path", artifact_path
    )
    benchmark_path = resolve_metadata_artifact_path(
        metadata, "radiation_transport_benchmark_artifact_path", artifact_path
    )

    if deposition_contract == RADIATION_DEPOSITION_HISTORY_CONTRACT:
        if deposition_path is None or not deposition_path.is_file():
            errors.append(f"missing_radiation_deposition_history:{deposition_path}")
        else:
            payload = load_json(deposition_path)
            if payload.get("schema") != RADIATION_DEPOSITION_HISTORY_SCHEMA:
                errors.append(
                    f"radiation_deposition_history_schema_mismatch:{deposition_path}:actual={payload.get('schema')}"
                )
            if payload.get("contract_id") != RADIATION_DEPOSITION_HISTORY_CONTRACT:
                errors.append(
                    f"radiation_deposition_history_contract_mismatch:{deposition_path}:actual={payload.get('contract_id')}"
                )

    if process_contract == RADIATION_PROCESS_HISTORY_CONTRACT:
        if process_path is None or not process_path.is_file():
            errors.append(f"missing_radiation_process_history:{process_path}")
        else:
            payload = load_json(process_path)
            if payload.get("schema") != RADIATION_PROCESS_HISTORY_SCHEMA:
                errors.append(
                    f"radiation_process_history_schema_mismatch:{process_path}:actual={payload.get('schema')}"
                )
            if payload.get("contract_id") != RADIATION_PROCESS_HISTORY_CONTRACT:
                errors.append(
                    f"radiation_process_history_contract_mismatch:{process_path}:actual={payload.get('contract_id')}"
                )
    if (
        material_governance_contract
        and material_governance_contract != RADIATION_MATERIAL_GOVERNANCE_CONTRACT
    ):
        errors.append(
            "radiation_material_governance_contract_mismatch:"
            f"{artifact_path}:actual={material_governance_contract}"
        )
    if benchmark_contract and benchmark_contract != RADIATION_TRANSPORT_BENCHMARK_CONTRACT:
        errors.append(
            f"radiation_transport_benchmark_contract_mismatch:{artifact_path}:actual={benchmark_contract}"
        )
    if benchmark_contract == RADIATION_TRANSPORT_BENCHMARK_CONTRACT:
        if benchmark_path is None or not benchmark_path.is_file():
            errors.append(f"missing_radiation_transport_benchmark:{benchmark_path}")
        else:
            payload = load_json(benchmark_path)
            if payload.get("schema_version") != RADIATION_TRANSPORT_BENCHMARK_SCHEMA:
                errors.append(
                    "radiation_transport_benchmark_schema_mismatch:"
                    f"{benchmark_path}:actual={payload.get('schema_version')}"
                )
            if payload.get("contract_id") != RADIATION_TRANSPORT_BENCHMARK_CONTRACT:
                errors.append(
                    "radiation_transport_benchmark_payload_contract_mismatch:"
                    f"{benchmark_path}:actual={payload.get('contract_id')}"
                )
            for key in (
                "transport_backend_mode",
                "energy_conservation_error",
                "layer_count",
                "process_coverage",
            ):
                if key not in payload:
                    errors.append(
                        f"radiation_transport_benchmark_missing_field:{benchmark_path}:{key}"
                    )
    return errors


def validate_vacuum_arc_benchmark_artifact(
    metadata: Dict[str, Any], artifact_path: Path
) -> List[str]:
    errors: List[str] = []
    pipeline_contract = str(metadata.get("pipeline_contract_id", "")).strip()
    if pipeline_contract and pipeline_contract != VACUUM_ARC_PIPELINE_CONTRACT:
        errors.append(
            f"vacuum_arc_pipeline_contract_mismatch:{artifact_path}:actual={pipeline_contract}"
        )
    collision_contract = str(
        metadata.get("collision_diagnostic_contract_id", "")
    ).strip()
    if (
        collision_contract
        and collision_contract != VACUUM_ARC_COLLISION_DIAGNOSTIC_CONTRACT
    ):
        errors.append(
            "vacuum_arc_collision_diagnostic_contract_mismatch:"
            f"{artifact_path}:actual={collision_contract}"
        )
    contract_id = str(metadata.get("benchmark_metrics_contract_id", "")).strip()
    if not contract_id:
        return errors
    if contract_id != VACUUM_ARC_BENCHMARK_METRICS_CONTRACT:
        errors.append(
            f"vacuum_arc_benchmark_contract_mismatch:{artifact_path}:actual={contract_id}"
        )
        return errors

    artifact_name = str(metadata.get("benchmark_metrics_artifact_path", "")).strip()
    benchmark_path = (
        artifact_path.parent / artifact_name
        if artifact_name
        else artifact_path.with_name(artifact_path.stem + ".benchmark_metrics.json")
    )
    if not benchmark_path.is_file():
        errors.append(f"missing_vacuum_arc_benchmark_artifact:{benchmark_path}")
        return errors

    payload = load_json(benchmark_path)
    if payload.get("schema_version") != VACUUM_ARC_BENCHMARK_METRICS_SCHEMA:
        errors.append(
            f"vacuum_arc_benchmark_schema_mismatch:{benchmark_path}:actual={payload.get('schema_version')}"
        )
    if payload.get("contract_id") != VACUUM_ARC_BENCHMARK_METRICS_CONTRACT:
        errors.append(
            f"vacuum_arc_benchmark_payload_contract_mismatch:{benchmark_path}:actual={payload.get('contract_id')}"
        )
    return errors


def validate_plasma_diagnostics_artifact(
    metadata: Dict[str, Any], artifact_path: Path
) -> List[str]:
    errors: List[str] = []
    reaction_contract = str(metadata.get("reaction_contract_id", "")).strip()
    if reaction_contract and reaction_contract != PLASMA_REACTION_BALANCE_CONTRACT:
        errors.append(
            f"plasma_reaction_contract_mismatch:{artifact_path}:actual={reaction_contract}"
        )
    contract_id = str(metadata.get("diagnostic_contract_id", "")).strip()
    if not contract_id:
        return errors
    if contract_id != PLASMA_PHYSICS_DIAGNOSTICS_CONTRACT:
        errors.append(
            f"plasma_diagnostics_contract_mismatch:{artifact_path}:actual={contract_id}"
        )
        return errors

    diagnostics_path = resolve_metadata_artifact_path(
        metadata, "plasma_physics_diagnostics_artifact_path", artifact_path
    )
    if diagnostics_path is None or not diagnostics_path.is_file():
        errors.append(f"missing_plasma_diagnostics_artifact:{diagnostics_path}")
        return errors

    payload = load_json(diagnostics_path)
    if payload.get("schema_version") != PLASMA_PHYSICS_DIAGNOSTICS_SCHEMA:
        errors.append(
            f"plasma_diagnostics_schema_mismatch:{diagnostics_path}:actual={payload.get('schema_version')}"
        )
    if payload.get("contract_id") != PLASMA_PHYSICS_DIAGNOSTICS_CONTRACT:
        errors.append(
            f"plasma_diagnostics_payload_contract_mismatch:{diagnostics_path}:actual={payload.get('contract_id')}"
        )
    if payload.get("reaction_contract_id") != str(metadata.get("reaction_contract_id", "")).strip():
        errors.append(
            f"plasma_diagnostics_reaction_contract_mismatch:{diagnostics_path}:actual={payload.get('reaction_contract_id')}"
        )
    if payload.get("diagnostic_set") != str(metadata.get("diagnostic_set", "")).strip():
        errors.append(
            f"plasma_diagnostics_set_mismatch:{diagnostics_path}:actual={payload.get('diagnostic_set')}"
        )
    if payload.get("environment_model") != str(metadata.get("environment_model", "")).strip():
        errors.append(
            f"plasma_diagnostics_environment_mismatch:{diagnostics_path}:actual={payload.get('environment_model')}"
        )
    if not isinstance(payload.get("boundary_layer_summary"), dict):
        errors.append(f"plasma_diagnostics_missing_boundary_layer_summary:{diagnostics_path}")
    if not isinstance(payload.get("reaction_balance_summary"), dict):
        errors.append(
            f"plasma_diagnostics_missing_reaction_balance_summary:{diagnostics_path}"
        )
    if not isinstance(payload.get("advanced_closure_summary"), dict):
        errors.append(
            f"plasma_diagnostics_missing_advanced_closure_summary:{diagnostics_path}"
        )
    if not isinstance(payload.get("consistency_flags"), dict):
        errors.append(f"plasma_diagnostics_missing_consistency_flags:{diagnostics_path}")
    return errors


def validate_internal_drive_artifacts(metadata: Dict[str, Any], artifact_path: Path) -> List[str]:
    errors: List[str] = []
    consumes_history = str(metadata.get("internal_consumes_deposition_history", "")).strip().lower()
    if consumes_history != "true":
        return errors

    deposition_path = resolve_metadata_artifact_path(metadata, "deposition_history_path", artifact_path)
    process_path = resolve_metadata_artifact_path(metadata, "process_history_path", artifact_path)
    provenance_contract = str(metadata.get("internal_drive_provenance_contract_id", "")).strip()
    provenance_path = resolve_metadata_artifact_path(
        metadata, "internal_drive_provenance_artifact_path", artifact_path
    )
    layer_alignment_contract = str(
        metadata.get("internal_drive_layer_alignment_contract_id", "")
    ).strip()
    layer_alignment_path = resolve_metadata_artifact_path(
        metadata, "internal_drive_layer_alignment_artifact_path", artifact_path
    )
    coupled_solve_contract = str(
        metadata.get("internal_response_coupled_solve_contract_id", "")
    ).strip()
    coupled_solve_path = resolve_metadata_artifact_path(
        metadata, "internal_response_coupled_solve_artifact_path", artifact_path
    )

    if deposition_path is None or not deposition_path.is_file():
        errors.append(f"missing_internal_deposition_history:{deposition_path}")
    if process_path is None or not process_path.is_file():
        errors.append(f"missing_internal_process_history:{process_path}")
    if provenance_contract and provenance_contract != INTERNAL_DRIVE_PROVENANCE_CONTRACT:
        errors.append(
            f"internal_drive_provenance_contract_mismatch:{artifact_path}:actual={provenance_contract}"
        )
    if provenance_contract == INTERNAL_DRIVE_PROVENANCE_CONTRACT:
        if provenance_path is None or not provenance_path.is_file():
            errors.append(f"missing_internal_drive_provenance_summary:{provenance_path}")
        else:
            payload = load_json(provenance_path)
            if payload.get("schema") != INTERNAL_DRIVE_PROVENANCE_SCHEMA:
                errors.append(
                    f"internal_drive_provenance_schema_mismatch:{provenance_path}:actual={payload.get('schema')}"
                )
            if payload.get("contract_id") != INTERNAL_DRIVE_PROVENANCE_CONTRACT:
                errors.append(
                    f"internal_drive_provenance_payload_contract_mismatch:{provenance_path}:actual={payload.get('contract_id')}"
                )
    if layer_alignment_contract and layer_alignment_contract != INTERNAL_DRIVE_LAYER_ALIGNMENT_CONTRACT:
        errors.append(
            f"internal_drive_layer_alignment_contract_mismatch:{artifact_path}:actual={layer_alignment_contract}"
        )
    if layer_alignment_contract == INTERNAL_DRIVE_LAYER_ALIGNMENT_CONTRACT:
        if layer_alignment_path is None or not layer_alignment_path.is_file():
            errors.append(f"missing_internal_drive_layer_alignment:{layer_alignment_path}")
        else:
            rows = layer_alignment_path.read_text(encoding="utf-8").splitlines()
            if len(rows) < 2:
                errors.append(
                    f"internal_drive_layer_alignment_empty:{layer_alignment_path}"
                )
            header = rows[0] if rows else ""
            for required_column in (
                "internal_layer_index",
                "radiation_layer_index",
                "normalized_depth_offset",
            ):
                if required_column not in header:
                    errors.append(
                        f"internal_drive_layer_alignment_missing_column:{layer_alignment_path}:{required_column}"
                    )
    if coupled_solve_contract and coupled_solve_contract != INTERNAL_RESPONSE_COUPLED_SOLVE_CONTRACT:
        errors.append(
            "internal_response_coupled_solve_contract_mismatch:"
            f"{artifact_path}:actual={coupled_solve_contract}"
        )
    if coupled_solve_contract == INTERNAL_RESPONSE_COUPLED_SOLVE_CONTRACT:
        if coupled_solve_path is None or not coupled_solve_path.is_file():
            errors.append(
                f"missing_internal_response_coupled_solve:{coupled_solve_path}"
            )
        else:
            payload = load_json(coupled_solve_path)
            if payload.get("schema_version") != INTERNAL_RESPONSE_COUPLED_SOLVE_SCHEMA:
                errors.append(
                    "internal_response_coupled_solve_schema_mismatch:"
                    f"{coupled_solve_path}:actual={payload.get('schema_version')}"
                )
            if payload.get("contract_id") != INTERNAL_RESPONSE_COUPLED_SOLVE_CONTRACT:
                errors.append(
                    "internal_response_coupled_solve_payload_contract_mismatch:"
                    f"{coupled_solve_path}:actual={payload.get('contract_id')}"
                )
            for key in (
                "coupled_solve_converged",
                "coupled_residual_metric",
                "alignment_max_normalized_depth_offset",
                "internal_response",
            ):
                if key not in payload:
                    errors.append(
                        f"internal_response_coupled_solve_missing_field:{coupled_solve_path}:{key}"
                    )
    return errors


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Output Contract v1 Gate Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- output_root: {report['output_root']}")
    lines.append("")
    lines.append("| module | command_status | artifacts_status | sidecar_status | contract_artifact_status | status |")
    lines.append("|---|---|---|---|---|---|")
    for item in report["modules"]:
        lines.append(
            "| "
            + item["module"]
            + " | "
            + item["command_status"]
            + " | "
            + item["artifacts_status"]
            + " | "
            + item["sidecar_status"]
            + " | "
            + item["contract_artifact_status"]
            + " | "
            + item["status"]
            + " |"
        )

    failures = report.get("failures", [])
    if failures:
        lines.append("")
        lines.append("## Failures")
        lines.append("")
        for failure in failures:
            lines.append(f"- {failure}")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="ENG-004 unified output contract v1 gate: verify csv/vtk/h5 sidecar metadata consistency"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--scdat-exe", required=True)
    parser.add_argument("--output-root", default="build/output_contract_v1_gate")
    parser.add_argument("--report-json", default="")
    parser.add_argument("--report-md", default="")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()

    scdat_exe = Path(args.scdat_exe)
    if not scdat_exe.is_absolute():
        scdat_exe = (project_root / scdat_exe).resolve()

    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = (project_root / output_root).resolve()

    report_json = Path(args.report_json) if args.report_json else output_root / "output_contract_v1_gate.json"
    if not report_json.is_absolute():
        report_json = (project_root / report_json).resolve()

    report_md = Path(args.report_md) if args.report_md else output_root / "output_contract_v1_gate.md"
    if not report_md.is_absolute():
        report_md = (project_root / report_md).resolve()

    try:
        if not scdat_exe.is_file():
            raise FileNotFoundError(f"missing SCDAT executable: {scdat_exe}")

        output_root.mkdir(parents=True, exist_ok=True)
        generated_root = output_root / "generated_configs"
        generated_root.mkdir(parents=True, exist_ok=True)

        module_reports: List[Dict[str, Any]] = []
        failures: List[str] = []

        for module_case in MODULE_CASES:
            module = module_case["module"]
            config_path = (project_root / module_case["config"]).resolve()
            module_failures: List[str] = []

            if not config_path.is_file():
                module_failures.append(f"missing_example_config:{config_path}")
                module_reports.append(
                    {
                        "module": module,
                        "command_status": "SKIP",
                        "artifacts_status": "SKIP",
                        "sidecar_status": "SKIP",
                        "contract_artifact_status": "SKIP",
                        "status": "FAIL",
                        "errors": module_failures,
                    }
                )
                failures.extend([f"{module}:{item}" for item in module_failures])
                continue

            payload = load_json(config_path)
            smoke_csv = output_root / f"{module}_output_contract_gate.csv"
            smoke_payload = make_smoke_payload(module_case, payload, smoke_csv)
            smoke_config = generated_root / f"{module}_output_contract_gate.json"
            smoke_config.write_text(json.dumps(smoke_payload, indent=2, ensure_ascii=False), encoding="utf-8")

            command = [
                str(scdat_exe),
                module_case["command"],
                str(smoke_config),
                str(smoke_csv),
            ]
            return_code, stdout_text, stderr_text = run_command(command, project_root)

            command_status = "PASS"
            if return_code != 0:
                command_status = "FAIL"
                module_failures.append(f"command_exit_code={return_code}")
                if stderr_text.strip():
                    module_failures.append(f"command_stderr={stderr_text.strip()}")

            stem = smoke_csv.with_suffix("")
            expected_artifacts = {
                "csv": smoke_csv.resolve(),
                "vtk": Path(str(stem) + EXPECTED_FORMAT_SUFFIX["vtk"]).resolve(),
                "h5": Path(str(stem) + EXPECTED_FORMAT_SUFFIX["h5"]).resolve(),
            }

            artifacts_status = "PASS"
            for artifact in expected_artifacts.values():
                if not artifact.is_file():
                    artifacts_status = "FAIL"
                    module_failures.append(f"missing_artifact:{artifact}")

            sidecar_status = "PASS"
            contract_artifact_status = "PASS"
            metadata_by_format: Dict[str, Dict[str, Any]] = {}
            for artifact_format, artifact_path in expected_artifacts.items():
                metadata, sidecar_errors = validate_sidecar(
                    artifact_format=artifact_format,
                    artifact_path=artifact_path,
                    expected_artifacts=expected_artifacts,
                )
                if sidecar_errors:
                    sidecar_status = "FAIL"
                    module_failures.extend(sidecar_errors)
                if metadata:
                    metadata_by_format[artifact_format] = metadata
                    if artifact_format == "csv":
                        extra_errors: List[str] = []
                        extra_errors.extend(validate_surface_boundary_observer(metadata, artifact_path))
                        extra_errors.extend(validate_surface_runtime_consistency(metadata, artifact_path))
                        extra_errors.extend(
                            validate_surface_particle_transport_domain(metadata, artifact_path)
                        )
                        extra_errors.extend(
                            validate_surface_global_particle_domain(metadata, artifact_path)
                        )
                        extra_errors.extend(
                            validate_surface_global_particle_repository(
                                metadata, artifact_path
                            )
                        )
                        extra_errors.extend(
                            validate_surface_global_sheath_field_solve(metadata, artifact_path)
                        )
                        extra_errors.extend(validate_vacuum_arc_benchmark_artifact(metadata, artifact_path))
                        extra_errors.extend(validate_plasma_diagnostics_artifact(metadata, artifact_path))
                        extra_errors.extend(validate_radiation_history_artifacts(metadata, artifact_path))
                        extra_errors.extend(validate_internal_drive_artifacts(metadata, artifact_path))
                        if extra_errors:
                            contract_artifact_status = "FAIL"
                            module_failures.extend(extra_errors)

            if {"csv", "vtk", "h5"}.issubset(metadata_by_format.keys()):
                reference = metadata_by_format["csv"]
                for peer_format in ("vtk", "h5"):
                    peer = metadata_by_format[peer_format]
                    for key in sorted(REQUIRED_METADATA_KEYS):
                        if reference.get(key) != peer.get(key):
                            sidecar_status = "FAIL"
                            module_failures.append(
                                f"cross_format_metadata_mismatch:{module}:{key}:csv={reference.get(key)}:{peer_format}={peer.get(key)}"
                            )

            module_status = "PASS"
            if (
                command_status != "PASS"
                or artifacts_status != "PASS"
                or sidecar_status != "PASS"
                or contract_artifact_status != "PASS"
            ):
                module_status = "FAIL"

            module_reports.append(
                {
                    "module": module,
                    "command": command,
                    "example_config": str(config_path),
                    "smoke_config": str(smoke_config),
                    "stdout": stdout_text.strip(),
                    "stderr": stderr_text.strip(),
                    "command_status": command_status,
                    "artifacts_status": artifacts_status,
                    "sidecar_status": sidecar_status,
                    "contract_artifact_status": contract_artifact_status,
                    "status": module_status,
                    "errors": module_failures,
                }
            )

            if module_failures:
                failures.extend([f"{module}:{item}" for item in module_failures])

        status = "PASS" if not failures else "FAIL"
        report: Dict[str, Any] = {
            "status": status,
            "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "scdat_exe": str(scdat_exe),
            "output_root": str(output_root),
            "modules": module_reports,
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
