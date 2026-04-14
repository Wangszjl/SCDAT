#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


def load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected top-level object in {path}")
    return payload


def deep_merge_dict(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    merged: Dict[str, Any] = dict(base)
    for key, value in override.items():
        existing = merged.get(key)
        if isinstance(existing, dict) and isinstance(value, dict):
            merged[key] = deep_merge_dict(existing, value)
        else:
            merged[key] = value
    return merged


def run_command(command: List[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )


def metadata_sidecar(csv_path: Path) -> Path:
    return Path(str(csv_path) + ".metadata.json")


def benchmark_case_sidecar(output_root: Path, metadata: Dict[str, Any], csv_path: Path) -> Path:
    artifact_name = str(metadata.get("benchmark_case_path", "")).strip()
    if artifact_name:
        return output_root / artifact_name
    return csv_path.with_suffix(".benchmark_case.json")


def simulation_artifact_sidecar(
    output_root: Path, metadata: Dict[str, Any], csv_path: Path
) -> Path:
    artifact_name = str(metadata.get("simulation_artifact_path", "")).strip()
    if artifact_name:
        return output_root / artifact_name
    return csv_path.with_suffix(".simulation_artifact.json")


def observer_sidecar(output_root: Path, metadata: Dict[str, Any], csv_path: Path) -> Path:
    artifact_name = str(metadata.get("surface_pic_runtime_boundary_observer_artifact", "")).strip()
    if artifact_name:
        return output_root / artifact_name
    return csv_path.with_suffix(".shared_runtime_observer.json")


def consistency_sidecar(output_root: Path, metadata: Dict[str, Any], csv_path: Path) -> Path:
    artifact_name = str(metadata.get("surface_pic_runtime_consistency_artifact", "")).strip()
    if artifact_name:
        return output_root / artifact_name
    return csv_path.with_suffix(".shared_runtime_consistency.json")


def transport_domain_sidecar(output_root: Path, metadata: Dict[str, Any], csv_path: Path) -> Path:
    artifact_name = str(
        metadata.get("surface_pic_runtime_shared_particle_transport_domain_artifact", "")
    ).strip()
    if artifact_name:
        return output_root / artifact_name
    return csv_path.with_suffix(".shared_particle_transport_domain.json")


def global_particle_domain_sidecar(
    output_root: Path, metadata: Dict[str, Any], csv_path: Path
) -> Path:
    artifact_name = str(
        metadata.get("surface_pic_runtime_global_particle_domain_artifact", "")
    ).strip()
    if artifact_name:
        return output_root / artifact_name
    return csv_path.with_suffix(".global_particle_domain.json")


def global_particle_repository_sidecar(
    output_root: Path, metadata: Dict[str, Any], csv_path: Path
) -> Path:
    artifact_name = str(
        metadata.get("surface_pic_runtime_global_particle_repository_artifact", "")
    ).strip()
    if artifact_name:
        return output_root / artifact_name
    return csv_path.with_suffix(".global_particle_repository.json")


def global_sheath_field_solve_sidecar(
    output_root: Path, metadata: Dict[str, Any], csv_path: Path
) -> Path:
    artifact_name = str(
        metadata.get("surface_pic_runtime_global_sheath_field_solve_artifact", "")
    ).strip()
    if artifact_name:
        return output_root / artifact_name
    return csv_path.with_suffix(".global_sheath_field_solve.json")


def validation_metric_values(validation_metrics: Any) -> Dict[str, float]:
    values: Dict[str, float] = {}
    if not isinstance(validation_metrics, list):
        return values

    for entry in validation_metrics:
        if not isinstance(entry, dict):
            continue
        metric_id = str(entry.get("id", "")).strip()
        if not metric_id:
            continue
        try:
            values[metric_id] = float(entry.get("value"))
        except (TypeError, ValueError):
            continue
    return values


def lookup_nested_value(payload: Dict[str, Any], dotted_path: str) -> Any:
    current: Any = payload
    for part in dotted_path.split("."):
        key = str(part).strip()
        if not key or not isinstance(current, dict) or key not in current:
            return None
        current = current[key]
    return current


def bool_to_yes_no(value: bool) -> str:
    return "yes" if value else "no"


def optional_bool_to_yes_no(value: Any) -> str:
    if value is True:
        return "yes"
    if value is False:
        return "no"
    return "n/a"


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    lines: List[str] = []
    lines.append("# Surface Reference Matrix Gate Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- output_root: {report['output_root']}")
    lines.append(f"- shim_removed_cases: {report.get('shim_removed_cases', 0)}")
    lines.append(f"- parity_passed_cases: {report.get('parity_passed_cases', 0)}")
    lines.append(
        f"- parity_applicable_cases: {report.get('parity_applicable_cases', 0)}"
    )
    lines.append("")
    lines.append("| preset | command | metadata | benchmark_case | observer | consistency | transport_domain | global_particle_domain | global_particle_repository | global_sheath_field | shim_removed | parity_passed | status |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for item in report["cases"]:
        shim_removed_label = bool_to_yes_no(bool(item.get("shim_removed", False)))
        parity_passed_label = optional_bool_to_yes_no(item.get("parity_passed"))
        lines.append(
            f"| {item['preset']} | {item['command_status']} | {item['metadata_status']} | {item['benchmark_case_status']} | {item['observer_status']} | {item['consistency_status']} | {item['transport_domain_status']} | {item['global_particle_domain_status']} | {item['global_particle_repository_status']} | {item['global_sheath_field_status']} | {shim_removed_label} | {parity_passed_label} | {item['status']} |"
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
        description="Run the Surface Charging reference family matrix gate."
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--scdat-exe", required=True)
    parser.add_argument(
        "--matrix",
        default="scripts/run/surface_reference_matrix.json",
        help="Matrix json path relative to project root unless absolute",
    )
    parser.add_argument("--output-root", default="build/surface_reference_matrix_gate")
    parser.add_argument(
        "--report-json",
        default=None,
        help="Optional report json path; defaults to <output-root>/surface_reference_matrix_gate.json",
    )
    parser.add_argument(
        "--report-md",
        default=None,
        help="Optional report markdown path; defaults to <output-root>/surface_reference_matrix_gate.md",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    scdat_exe = Path(args.scdat_exe)
    if not scdat_exe.is_absolute():
        scdat_exe = (project_root / scdat_exe).resolve()

    matrix_path = Path(args.matrix)
    if not matrix_path.is_absolute():
        matrix_path = (project_root / matrix_path).resolve()

    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = (project_root / output_root).resolve()

    report_json_path = (
        Path(args.report_json)
        if args.report_json
        else output_root / "surface_reference_matrix_gate.json"
    )
    if not report_json_path.is_absolute():
        report_json_path = (project_root / report_json_path).resolve()

    report_md_path = (
        Path(args.report_md)
        if args.report_md
        else output_root / "surface_reference_matrix_gate.md"
    )
    if not report_md_path.is_absolute():
        report_md_path = (project_root / report_md_path).resolve()

    if not scdat_exe.is_file():
        raise FileNotFoundError(f"missing SCDAT executable: {scdat_exe}")
    if not matrix_path.is_file():
        raise FileNotFoundError(f"missing surface reference matrix: {matrix_path}")

    matrix = load_json(matrix_path)
    cases = matrix.get("cases", [])
    if not isinstance(cases, list) or not cases:
        raise ValueError(f"matrix must contain non-empty cases list: {matrix_path}")

    output_root.mkdir(parents=True, exist_ok=True)
    case_reports: List[Dict[str, Any]] = []
    failures: List[str] = []

    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            failures.append(f"case[{index}]:invalid_case_type")
            continue

        preset = str(case.get("preset", "")).strip()
        case_label = str(case.get("case_label", preset)).strip() or preset
        expected_family = str(case.get("expected_reference_family", "")).strip()
        expected_matrix_case = str(case.get("expected_reference_matrix_case_id", "")).strip()
        expected_boundary_observer_contract = str(
            case.get("expected_boundary_observer_contract_id", "")
        ).strip()
        expected_metadata_values_raw = case.get("expected_metadata_values", {})
        expected_metadata_values = (
            expected_metadata_values_raw
            if isinstance(expected_metadata_values_raw, dict)
            else {}
        )
        expected_shim_removed_raw = case.get("expected_shim_removed")
        expected_shim_removed = (
            bool(expected_shim_removed_raw)
            if isinstance(expected_shim_removed_raw, bool)
            else False
        )
        expected_surface_pic_runtime = str(case.get("expected_surface_pic_runtime", "")).strip()
        expected_shared_runtime_enabled = case.get("expected_shared_runtime_enabled")
        expected_consistency_contract = str(
            case.get("expected_consistency_contract_id", "")
        ).strip()
        require_surface_runtime_consistency_pass = case.get(
            "require_surface_runtime_consistency_pass", True
        )
        max_reference_spread_v = case.get("max_reference_spread_v")
        max_sheath_charge_spread_c = case.get("max_sheath_charge_spread_c")
        max_patch_potential_spread_v = case.get("max_patch_potential_spread_v")
        min_patch_potential_spread_reduction_ratio = case.get(
            "min_patch_potential_spread_reduction_ratio"
        )
        max_live_pic_net_current_spread_a_per_m2 = case.get(
            "max_live_pic_net_current_spread_a_per_m2"
        )
        max_electron_calibration_factor_spread = case.get(
            "max_electron_calibration_factor_spread"
        )
        require_global_sheath_proxy_solve_active = case.get(
            "require_global_sheath_proxy_solve_active"
        )
        require_shared_current_matrix_coupling_active = case.get(
            "require_shared_current_matrix_coupling_active"
        )
        min_shared_current_matrix_coupling_offdiag_entry_count = case.get(
            "min_shared_current_matrix_coupling_offdiag_entry_count"
        )
        require_shared_global_coupled_solve_active = case.get(
            "require_shared_global_coupled_solve_active"
        )
        require_shared_global_coupled_solve_converged = case.get(
            "require_shared_global_coupled_solve_converged"
        )
        min_shared_global_coupled_solve_iterations = case.get(
            "min_shared_global_coupled_solve_iterations"
        )
        max_shared_global_coupled_solve_max_delta_v = case.get(
            "max_shared_global_coupled_solve_max_delta_v"
        )
        max_normal_electric_field_spread_v_per_m = case.get(
            "max_normal_electric_field_spread_v_per_m"
        )
        max_local_charge_density_spread_c_per_m3 = case.get(
            "max_local_charge_density_spread_c_per_m3"
        )
        require_shared_live_pic_coupled_refresh_active = case.get(
            "require_shared_live_pic_coupled_refresh_active"
        )
        min_shared_live_pic_coupled_refresh_count = case.get(
            "min_shared_live_pic_coupled_refresh_count"
        )
        require_shared_particle_transport_coupling_active = case.get(
            "require_shared_particle_transport_coupling_active"
        )
        min_shared_particle_transport_offdiag_entry_count = case.get(
            "min_shared_particle_transport_offdiag_entry_count"
        )
        max_shared_particle_transport_conservation_error_a_per_v = case.get(
            "max_shared_particle_transport_conservation_error_a_per_v"
        )
        min_abs_shared_particle_transport_reference_shift_v = case.get(
            "min_abs_shared_particle_transport_reference_shift_v"
        )
        min_abs_shared_particle_transport_charge_c = case.get(
            "min_abs_shared_particle_transport_charge_c"
        )
        require_shared_particle_transport_distribution_active = case.get(
            "require_shared_particle_transport_distribution_active"
        )
        min_shared_particle_transport_distribution_charge_spread_c = case.get(
            "min_shared_particle_transport_distribution_charge_spread_c"
        )
        min_shared_particle_transport_distribution_reference_shift_spread_v = case.get(
            "min_shared_particle_transport_distribution_reference_shift_spread_v"
        )
        max_shared_particle_transport_distribution_conservation_error_c = case.get(
            "max_shared_particle_transport_distribution_conservation_error_c"
        )
        require_shared_particle_transport_exchange_active = case.get(
            "require_shared_particle_transport_exchange_active"
        )
        min_shared_particle_transport_exchange_flux_spread_a = case.get(
            "min_shared_particle_transport_exchange_flux_spread_a"
        )
        max_shared_particle_transport_exchange_flux_conservation_error_a = case.get(
            "max_shared_particle_transport_exchange_flux_conservation_error_a"
        )
        expected_transport_domain_contract = str(
            case.get("expected_transport_domain_contract_id", "")
        ).strip()
        expected_shared_particle_transport_bookkeeping_mode = str(
            case.get("expected_shared_particle_transport_bookkeeping_mode", "")
        ).strip()
        require_shared_particle_transport_domain_runtime_state_backed = case.get(
            "require_shared_particle_transport_domain_runtime_state_backed"
        )
        require_shared_particle_transport_domain_active = case.get(
            "require_shared_particle_transport_domain_active"
        )
        min_shared_particle_transport_domain_node_count = case.get(
            "min_shared_particle_transport_domain_node_count"
        )
        min_shared_particle_transport_domain_edge_count = case.get(
            "min_shared_particle_transport_domain_edge_count"
        )
        max_shared_particle_transport_domain_exchange_flux_conservation_error_a = case.get(
            "max_shared_particle_transport_domain_exchange_flux_conservation_error_a"
        )
        min_shared_particle_transport_domain_edge_charge_total_abs_c = case.get(
            "min_shared_particle_transport_domain_edge_charge_total_abs_c"
        )
        max_shared_particle_transport_domain_edge_charge_conservation_error_c = case.get(
            "max_shared_particle_transport_domain_edge_charge_conservation_error_c"
        )
        require_shared_particle_transport_edge_operator_active = case.get(
            "require_shared_particle_transport_edge_operator_active"
        )
        min_shared_particle_transport_edge_operator_total_abs_drive_charge_c = case.get(
            "min_shared_particle_transport_edge_operator_total_abs_drive_charge_c"
        )
        min_shared_particle_transport_domain_edge_operator_drive_total_abs_c = case.get(
            "min_shared_particle_transport_domain_edge_operator_drive_total_abs_c"
        )
        max_shared_particle_transport_domain_edge_operator_drive_conservation_error_c = case.get(
            "max_shared_particle_transport_domain_edge_operator_drive_conservation_error_c"
        )
        require_shared_particle_transport_edge_graph_operator_converged = case.get(
            "require_shared_particle_transport_edge_graph_operator_converged"
        )
        min_shared_particle_transport_edge_graph_operator_iterations = case.get(
            "min_shared_particle_transport_edge_graph_operator_iterations"
        )
        max_shared_particle_transport_edge_graph_operator_max_balance_residual_c = case.get(
            "max_shared_particle_transport_edge_graph_operator_max_balance_residual_c"
        )
        require_shared_particle_transport_domain_edge_graph_operator_converged = case.get(
            "require_shared_particle_transport_domain_edge_graph_operator_converged"
        )
        min_shared_particle_transport_domain_edge_graph_operator_iterations = case.get(
            "min_shared_particle_transport_domain_edge_graph_operator_iterations"
        )
        max_shared_particle_transport_domain_edge_graph_operator_max_balance_residual_c = case.get(
            "max_shared_particle_transport_domain_edge_graph_operator_max_balance_residual_c"
        )
        require_shared_particle_transport_edge_graph_operator_branch_graph_active = case.get(
            "require_shared_particle_transport_edge_graph_operator_branch_graph_active"
        )
        min_shared_particle_transport_edge_graph_operator_branch_graph_edge_count = case.get(
            "min_shared_particle_transport_edge_graph_operator_branch_graph_edge_count"
        )
        min_shared_particle_transport_edge_graph_operator_branch_graph_pair_count = case.get(
            "min_shared_particle_transport_edge_graph_operator_branch_graph_pair_count"
        )
        min_shared_particle_transport_edge_graph_operator_total_conductance_weight_f = case.get(
            "min_shared_particle_transport_edge_graph_operator_total_conductance_weight_f"
        )
        min_shared_particle_transport_edge_graph_operator_max_node_preconditioner = case.get(
            "min_shared_particle_transport_edge_graph_operator_max_node_preconditioner"
        )
        require_shared_particle_transport_domain_edge_graph_operator_branch_graph_active = case.get(
            "require_shared_particle_transport_domain_edge_graph_operator_branch_graph_active"
        )
        min_shared_particle_transport_domain_edge_graph_operator_branch_graph_edge_count = case.get(
            "min_shared_particle_transport_domain_edge_graph_operator_branch_graph_edge_count"
        )
        min_shared_particle_transport_domain_edge_graph_operator_branch_graph_pair_count = case.get(
            "min_shared_particle_transport_domain_edge_graph_operator_branch_graph_pair_count"
        )
        min_shared_particle_transport_domain_edge_graph_operator_total_conductance_weight_f = case.get(
            "min_shared_particle_transport_domain_edge_graph_operator_total_conductance_weight_f"
        )
        min_shared_particle_transport_domain_edge_graph_operator_max_node_preconditioner = case.get(
            "min_shared_particle_transport_domain_edge_graph_operator_max_node_preconditioner"
        )
        expected_global_particle_domain_contract = str(
            case.get("expected_global_particle_domain_contract_id", "")
        ).strip()
        expected_global_particle_bookkeeping_mode = str(
            case.get("expected_global_particle_bookkeeping_mode", "")
        ).strip()
        require_global_particle_runtime_state_backed = case.get(
            "require_global_particle_runtime_state_backed"
        )
        require_global_particle_domain_active = case.get(
            "require_global_particle_domain_active"
        )
        min_global_particle_domain_node_count = case.get(
            "min_global_particle_domain_node_count"
        )
        min_global_particle_domain_edge_count = case.get(
            "min_global_particle_domain_edge_count"
        )
        max_global_particle_domain_charge_conservation_error_c = case.get(
            "max_global_particle_domain_charge_conservation_error_c"
        )
        max_global_particle_domain_flux_conservation_error_a = case.get(
            "max_global_particle_domain_flux_conservation_error_a"
        )
        min_global_particle_domain_edge_charge_total_abs_c = case.get(
            "min_global_particle_domain_edge_charge_total_abs_c"
        )
        min_global_particle_domain_edge_target_charge_total_abs_c = case.get(
            "min_global_particle_domain_edge_target_charge_total_abs_c"
        )
        min_global_particle_domain_edge_operator_drive_total_abs_c = case.get(
            "min_global_particle_domain_edge_operator_drive_total_abs_c"
        )
        min_global_particle_domain_edge_conductance_total_s = case.get(
            "min_global_particle_domain_edge_conductance_total_s"
        )
        require_global_particle_domain_edge_graph_operator_converged = case.get(
            "require_global_particle_domain_edge_graph_operator_converged"
        )
        expected_global_particle_repository_contract = str(
            case.get("expected_global_particle_repository_contract_id", "")
        ).strip()
        expected_global_particle_repository_bookkeeping_mode = str(
            case.get("expected_global_particle_repository_bookkeeping_mode", "")
        ).strip()
        expected_global_particle_repository_lifecycle_mode = str(
            case.get("expected_global_particle_repository_lifecycle_mode", "")
        ).strip()
        require_global_particle_repository_runtime_state_backed = case.get(
            "require_global_particle_repository_runtime_state_backed"
        )
        require_global_particle_repository_active = case.get(
            "require_global_particle_repository_active"
        )
        min_global_particle_repository_node_count = case.get(
            "min_global_particle_repository_node_count"
        )
        min_global_particle_repository_edge_count = case.get(
            "min_global_particle_repository_edge_count"
        )
        max_global_particle_repository_charge_conservation_error_c = case.get(
            "max_global_particle_repository_charge_conservation_error_c"
        )
        max_global_particle_repository_migration_charge_conservation_error_c = case.get(
            "max_global_particle_repository_migration_charge_conservation_error_c"
        )
        min_global_particle_repository_total_migration_delta_abs_charge_c = case.get(
            "min_global_particle_repository_total_migration_delta_abs_charge_c"
        )
        min_global_particle_repository_total_edge_feedback_abs_charge_c = case.get(
            "min_global_particle_repository_total_edge_feedback_abs_charge_c"
        )
        min_global_particle_repository_total_conservation_correction_abs_charge_c = case.get(
            "min_global_particle_repository_total_conservation_correction_abs_charge_c"
        )
        min_global_particle_repository_total_migration_edge_abs_charge_c = case.get(
            "min_global_particle_repository_total_migration_edge_abs_charge_c"
        )
        expected_global_sheath_field_solve_contract = str(
            case.get("expected_global_sheath_field_solve_contract_id", "")
        ).strip()
        expected_global_sheath_field_solve_mode = str(
            case.get("expected_global_sheath_field_solve_mode", "")
        ).strip()
        require_global_sheath_field_runtime_state_backed = case.get(
            "require_global_sheath_field_runtime_state_backed"
        )
        require_global_sheath_field_solve_active = case.get(
            "require_global_sheath_field_solve_active"
        )
        require_global_sheath_field_solve_converged = case.get(
            "require_global_sheath_field_solve_converged"
        )
        min_global_sheath_field_solve_node_count = case.get(
            "min_global_sheath_field_solve_node_count"
        )
        max_global_sheath_field_residual_v_per_m = case.get(
            "max_global_sheath_field_residual_v_per_m"
        )
        max_global_particle_field_coupled_residual_v = case.get(
            "max_global_particle_field_coupled_residual_v"
        )
        max_global_sheath_field_multi_step_stability_metric_v = case.get(
            "max_global_sheath_field_multi_step_stability_metric_v"
        )
        min_global_sheath_field_matrix_row_count = case.get(
            "min_global_sheath_field_matrix_row_count"
        )
        min_global_sheath_field_matrix_nonzeros = case.get(
            "min_global_sheath_field_matrix_nonzeros"
        )
        max_global_sheath_field_linear_residual_norm_v = case.get(
            "max_global_sheath_field_linear_residual_norm_v"
        )
        min_shared_patch_node_count = case.get("min_shared_patch_node_count")
        max_benchmark_case_charge_conservation_drift_ratio = case.get(
            "max_benchmark_case_charge_conservation_drift_ratio"
        )
        max_benchmark_case_convergence_failure_rate = case.get(
            "max_benchmark_case_convergence_failure_rate"
        )
        max_benchmark_case_curve_relative_rmse = case.get(
            "max_benchmark_case_curve_relative_rmse"
        )
        max_benchmark_case_drift_tolerance = case.get(
            "max_benchmark_case_drift_tolerance"
        )
        expected_benchmark_reference_source_resolved = str(
            case.get("expected_benchmark_reference_source_resolved", "")
        ).strip()
        min_benchmark_case_reference_dataset_count = case.get(
            "min_benchmark_case_reference_dataset_count"
        )
        require_benchmark_case_curve_relative_rmse = case.get(
            "require_benchmark_case_curve_relative_rmse"
        )
        require_benchmark_case_acceptance_gate_metric = case.get(
            "require_benchmark_case_acceptance_gate_metric"
        )
        required_benchmark_case_metric_ids = case.get(
            "required_benchmark_case_metric_ids", []
        )
        max_benchmark_case_metric_values = case.get(
            "max_benchmark_case_metric_values", {}
        )
        required_simulation_artifact_metric_paths = case.get(
            "required_simulation_artifact_metric_paths", []
        )
        max_simulation_artifact_metric_values = case.get(
            "max_simulation_artifact_metric_values", {}
        )
        config_overrides = case.get("config_overrides", {})
        run_overrides = case.get("run_overrides", {})
        if not preset or not expected_family or not expected_matrix_case:
            failures.append(f"case[{index}]:missing_required_fields")
            continue
        if not isinstance(config_overrides, dict):
            failures.append(f"case[{index}]:config_overrides_must_be_object")
            continue
        if not isinstance(run_overrides, dict):
            failures.append(f"case[{index}]:run_overrides_must_be_object")
            continue
        if not isinstance(required_benchmark_case_metric_ids, list):
            failures.append(f"case[{index}]:required_benchmark_case_metric_ids_must_be_array")
            continue
        if not isinstance(max_benchmark_case_metric_values, dict):
            failures.append(f"case[{index}]:max_benchmark_case_metric_values_must_be_object")
            continue
        if not isinstance(required_simulation_artifact_metric_paths, list):
            failures.append(
                f"case[{index}]:required_simulation_artifact_metric_paths_must_be_array"
            )
            continue
        if not isinstance(max_simulation_artifact_metric_values, dict):
            failures.append(
                f"case[{index}]:max_simulation_artifact_metric_values_must_be_object"
            )
            continue

        csv_path = output_root / f"{case_label}.csv"
        config_path = output_root / f"{case_label}.surface_matrix_gate.json"
        config_payload: Dict[str, Any] = {
            "schema_version": "v1",
            "module": "surface",
            "base_preset": preset,
            "run": {
                "steps": 1,
                "adaptive_time_stepping": False,
                "total_duration_s": 0.0,
                "output_csv": str(csv_path),
            },
            "config": {
                "reference_family": expected_family,
                "reference_matrix_case_id": expected_matrix_case,
            },
        }
        if run_overrides:
            config_payload["run"] = deep_merge_dict(config_payload["run"], run_overrides)
        if config_overrides:
            config_payload["config"] = deep_merge_dict(
                config_payload["config"], config_overrides
            )
        config_path.write_text(
            json.dumps(config_payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

        command = [str(scdat_exe), "surface-config", str(config_path), str(csv_path)]
        completed = run_command(command, project_root)
        case_failure_start = len(failures)

        command_status = "PASS" if completed.returncode == 0 else "FAIL"
        metadata_status = "SKIP"
        benchmark_case_status = "SKIP"
        observer_status = "SKIP"
        consistency_status = "SKIP"
        transport_domain_status = "SKIP"
        global_particle_domain_status = "SKIP"
        global_particle_repository_status = "SKIP"
        global_sheath_field_status = "SKIP"
        status = "FAIL"
        metadata: Dict[str, Any] = {}
        shim_removed = expected_shim_removed
        parity_passed: Any = None

        if completed.returncode != 0:
            failures.append(
                f"{preset}:command_failed:code={completed.returncode}:stderr={completed.stderr.strip()}"
            )
        else:
            sidecar_path = metadata_sidecar(csv_path)
            if not sidecar_path.is_file():
                failures.append(f"{preset}:missing_metadata_sidecar:{sidecar_path}")
            else:
                sidecar = load_json(sidecar_path)
                metadata = sidecar.get("metadata", {})
                surface_legacy_input_adapter = str(
                    metadata.get("surface_legacy_input_adapter", "")
                ).strip()
                if surface_legacy_input_adapter:
                    shim_removed = surface_legacy_input_adapter.lower() == "none"
                actual_family = str(metadata.get("surface_reference_family", ""))
                actual_matrix_case = str(metadata.get("surface_reference_matrix_case_id", ""))
                if actual_family != expected_family:
                    failures.append(
                        f"{preset}:reference_family_mismatch:expected={expected_family}:actual={actual_family}"
                    )
                elif actual_matrix_case != expected_matrix_case:
                    failures.append(
                        f"{preset}:reference_matrix_case_id_mismatch:expected={expected_matrix_case}:actual={actual_matrix_case}"
                    )
                elif (
                    expected_benchmark_reference_source_resolved
                    and str(
                        metadata.get("benchmark_reference_source_resolved", "")
                    ).strip()
                    != expected_benchmark_reference_source_resolved
                ):
                    failures.append(
                        f"{preset}:benchmark_reference_source_resolved_mismatch:"
                        f"expected={expected_benchmark_reference_source_resolved}:"
                        f"actual={str(metadata.get('benchmark_reference_source_resolved', '')).strip()}"
                    )
                elif expected_metadata_values and any(
                    str(metadata.get(str(key), "")).strip() != str(value).strip()
                    for key, value in expected_metadata_values.items()
                ):
                    for key, value in expected_metadata_values.items():
                        actual_value = str(metadata.get(str(key), "")).strip()
                        expected_value = str(value).strip()
                        if actual_value != expected_value:
                            failures.append(
                                f"{preset}:metadata_value_mismatch:key={key}:"
                                f"expected={expected_value}:actual={actual_value}"
                            )
                else:
                    metadata_status = "PASS"
                    benchmark_case_contract_id = str(
                        metadata.get("benchmark_case_contract_id", "")
                    ).strip()
                    if benchmark_case_contract_id != "benchmark-case-v1":
                        failures.append(
                            f"{preset}:benchmark_case_contract_mismatch:"
                            f"expected=benchmark-case-v1:actual={benchmark_case_contract_id}"
                        )
                    else:
                        benchmark_case_path = benchmark_case_sidecar(
                            output_root, metadata, csv_path
                        )
                        if not benchmark_case_path.is_file():
                            failures.append(
                                f"{preset}:missing_benchmark_case_sidecar:{benchmark_case_path}"
                            )
                        else:
                            benchmark_case_payload = load_json(benchmark_case_path)
                            actual_schema = str(
                                benchmark_case_payload.get("schema_version", "")
                            )
                            actual_contract_id = str(
                                benchmark_case_payload.get("contract_id", "")
                            )
                            actual_case_id = str(benchmark_case_payload.get("id", ""))
                            actual_module = str(benchmark_case_payload.get("module", ""))
                            actual_reference_family = str(
                                benchmark_case_payload.get("reference_family", "")
                            )
                            validation_metrics = benchmark_case_payload.get(
                                "validation_metrics"
                            )
                            tolerance_profile = benchmark_case_payload.get(
                                "tolerance_profile"
                            )
                            reference_datasets = benchmark_case_payload.get(
                                "reference_datasets"
                            )
                            metric_values = validation_metric_values(validation_metrics)

                            if actual_schema != "scdat.benchmark_case.v1":
                                failures.append(
                                    f"{preset}:benchmark_case_schema_mismatch:"
                                    f"actual={actual_schema}"
                                )
                            elif actual_contract_id != "benchmark-case-v1":
                                failures.append(
                                    f"{preset}:benchmark_case_payload_contract_mismatch:"
                                    f"expected=benchmark-case-v1:actual={actual_contract_id}"
                                )
                            elif actual_case_id != expected_matrix_case:
                                failures.append(
                                    f"{preset}:benchmark_case_id_mismatch:"
                                    f"expected={expected_matrix_case}:actual={actual_case_id}"
                                )
                            elif actual_module != "surface":
                                failures.append(
                                    f"{preset}:benchmark_case_module_mismatch:"
                                    f"expected=surface:actual={actual_module}"
                                )
                            elif actual_reference_family != expected_family:
                                failures.append(
                                    f"{preset}:benchmark_case_reference_family_mismatch:"
                                    f"expected={expected_family}:actual={actual_reference_family}"
                                )
                            elif not isinstance(validation_metrics, list) or not validation_metrics:
                                failures.append(
                                    f"{preset}:benchmark_case_missing_validation_metrics"
                                )
                            elif not isinstance(tolerance_profile, dict):
                                failures.append(
                                    f"{preset}:benchmark_case_missing_tolerance_profile"
                                )
                            elif "rmse_tolerance" not in tolerance_profile:
                                failures.append(
                                    f"{preset}:benchmark_case_missing_rmse_tolerance"
                                )
                            elif (
                                max_benchmark_case_drift_tolerance is not None
                                and "drift_tolerance" not in tolerance_profile
                            ):
                                failures.append(
                                    f"{preset}:benchmark_case_missing_drift_tolerance"
                                )
                            elif (
                                max_benchmark_case_drift_tolerance is not None
                                and float(tolerance_profile.get("drift_tolerance", 0.0))
                                > float(max_benchmark_case_drift_tolerance)
                            ):
                                failures.append(
                                    f"{case_label}:benchmark_case_drift_tolerance_exceeded:"
                                    f"limit={float(max_benchmark_case_drift_tolerance)}:"
                                    f"actual={float(tolerance_profile.get('drift_tolerance', 0.0))}"
                                )
                            elif (
                                max_benchmark_case_charge_conservation_drift_ratio is not None
                                and "charge_conservation_drift_ratio" not in metric_values
                            ):
                                failures.append(
                                    f"{preset}:benchmark_case_missing_charge_conservation_drift_ratio"
                                )
                            elif (
                                max_benchmark_case_charge_conservation_drift_ratio is not None
                                and metric_values["charge_conservation_drift_ratio"]
                                > float(max_benchmark_case_charge_conservation_drift_ratio)
                            ):
                                failures.append(
                                    f"{case_label}:benchmark_case_charge_conservation_drift_ratio_exceeded:"
                                    f"limit={float(max_benchmark_case_charge_conservation_drift_ratio)}:"
                                    f"actual={metric_values['charge_conservation_drift_ratio']}"
                                )
                            elif (
                                max_benchmark_case_convergence_failure_rate is not None
                                and "shared_global_coupled_convergence_failure_rate"
                                not in metric_values
                            ):
                                failures.append(
                                    f"{preset}:benchmark_case_missing_convergence_failure_rate"
                                )
                            elif (
                                max_benchmark_case_convergence_failure_rate is not None
                                and metric_values[
                                    "shared_global_coupled_convergence_failure_rate"
                                ]
                                > float(max_benchmark_case_convergence_failure_rate)
                            ):
                                failures.append(
                                    f"{case_label}:benchmark_case_convergence_failure_rate_exceeded:"
                                    f"limit={float(max_benchmark_case_convergence_failure_rate)}:"
                                    f"actual={metric_values['shared_global_coupled_convergence_failure_rate']}"
                                )
                            elif (
                                max_benchmark_case_curve_relative_rmse is not None
                                and "benchmark_curve_relative_rmse_max" not in metric_values
                            ):
                                failures.append(
                                    f"{preset}:benchmark_case_missing_curve_relative_rmse_max"
                                )
                            elif (
                                max_benchmark_case_curve_relative_rmse is not None
                                and metric_values["benchmark_curve_relative_rmse_max"]
                                > float(max_benchmark_case_curve_relative_rmse)
                            ):
                                failures.append(
                                    f"{case_label}:benchmark_case_curve_relative_rmse_exceeded:"
                                    f"limit={float(max_benchmark_case_curve_relative_rmse)}:"
                                    f"actual={metric_values['benchmark_curve_relative_rmse_max']}"
                                )
                            elif not isinstance(reference_datasets, list):
                                failures.append(
                                    f"{preset}:benchmark_case_invalid_reference_datasets"
                                )
                            elif (
                                min_benchmark_case_reference_dataset_count is not None
                                and len(reference_datasets)
                                < int(min_benchmark_case_reference_dataset_count)
                            ):
                                failures.append(
                                    f"{case_label}:benchmark_case_reference_dataset_count_too_small:"
                                    f"minimum={int(min_benchmark_case_reference_dataset_count)}:"
                                    f"actual={len(reference_datasets)}"
                                )
                            else:
                                if (
                                    require_benchmark_case_curve_relative_rmse
                                    and "benchmark_curve_relative_rmse_max"
                                    not in metric_values
                                ):
                                    failures.append(
                                        f"{preset}:benchmark_case_missing_curve_relative_rmse_metric"
                                    )
                                    continue
                                if (
                                    require_benchmark_case_acceptance_gate_metric
                                    and "benchmark_acceptance_gate_pass"
                                    not in metric_values
                                ):
                                    failures.append(
                                        f"{preset}:benchmark_case_missing_acceptance_gate_metric"
                                    )
                                    continue
                                missing_metric_ids = [
                                    str(metric_id)
                                    for metric_id in required_benchmark_case_metric_ids
                                    if str(metric_id) not in metric_values
                                ]
                                if missing_metric_ids:
                                    failures.append(
                                        f"{preset}:benchmark_case_missing_required_metrics:"
                                        + ",".join(sorted(missing_metric_ids))
                                    )
                                    continue
                                metric_limit_exceeded = False
                                for metric_id, raw_limit in max_benchmark_case_metric_values.items():
                                    metric_key = str(metric_id)
                                    if metric_key not in metric_values:
                                        failures.append(
                                            f"{preset}:benchmark_case_missing_threshold_metric:{metric_key}"
                                        )
                                        metric_limit_exceeded = True
                                        break
                                    try:
                                        limit_value = float(raw_limit)
                                    except (TypeError, ValueError):
                                        failures.append(
                                            f"{preset}:benchmark_case_invalid_threshold_value:{metric_key}"
                                        )
                                        metric_limit_exceeded = True
                                        break
                                    if metric_values[metric_key] > limit_value:
                                        failures.append(
                                            f"{case_label}:benchmark_case_metric_exceeded:"
                                            f"id={metric_key}:limit={limit_value}:actual={metric_values[metric_key]}"
                                        )
                                        metric_limit_exceeded = True
                                        break
                                if metric_limit_exceeded:
                                    continue
                                simulation_artifact_contract_id = str(
                                    metadata.get("simulation_artifact_contract_id", "")
                                ).strip()
                                if simulation_artifact_contract_id != "simulation-artifact-v1":
                                    failures.append(
                                        f"{preset}:simulation_artifact_contract_mismatch:"
                                        f"expected=simulation-artifact-v1:actual={simulation_artifact_contract_id}"
                                    )
                                    continue
                                simulation_artifact_path = simulation_artifact_sidecar(
                                    output_root, metadata, csv_path
                                )
                                if not simulation_artifact_path.is_file():
                                    failures.append(
                                        f"{preset}:missing_simulation_artifact_sidecar:{simulation_artifact_path}"
                                    )
                                    continue
                                simulation_artifact_payload = load_json(
                                    simulation_artifact_path
                                )
                                if (
                                    str(
                                        simulation_artifact_payload.get(
                                            "schema_version", ""
                                        )
                                    )
                                    != "scdat.simulation_artifact.v1"
                                ):
                                    failures.append(
                                        f"{preset}:simulation_artifact_schema_mismatch:"
                                        f"actual={simulation_artifact_payload.get('schema_version')}"
                                    )
                                    continue
                                if (
                                    str(
                                        simulation_artifact_payload.get("contract_id", "")
                                    )
                                    != "simulation-artifact-v1"
                                ):
                                    failures.append(
                                        f"{preset}:simulation_artifact_payload_contract_mismatch:"
                                        f"expected=simulation-artifact-v1:"
                                        f"actual={simulation_artifact_payload.get('contract_id')}"
                                    )
                                    continue

                                missing_simulation_metric_paths = [
                                    str(metric_path)
                                    for metric_path in required_simulation_artifact_metric_paths
                                    if lookup_nested_value(
                                        simulation_artifact_payload,
                                        str(metric_path),
                                    )
                                    is None
                                ]
                                if missing_simulation_metric_paths:
                                    failures.append(
                                        f"{preset}:simulation_artifact_missing_required_metrics:"
                                        + ",".join(
                                            sorted(missing_simulation_metric_paths)
                                        )
                                    )
                                    continue

                                simulation_metric_limit_exceeded = False
                                for metric_path, raw_limit in (
                                    max_simulation_artifact_metric_values.items()
                                ):
                                    metric_key = str(metric_path)
                                    actual_value = lookup_nested_value(
                                        simulation_artifact_payload, metric_key
                                    )
                                    if actual_value is None:
                                        failures.append(
                                            f"{preset}:simulation_artifact_missing_threshold_metric:{metric_key}"
                                        )
                                        simulation_metric_limit_exceeded = True
                                        break
                                    try:
                                        actual_numeric = float(actual_value)
                                        limit_value = float(raw_limit)
                                    except (TypeError, ValueError):
                                        failures.append(
                                            f"{preset}:simulation_artifact_invalid_threshold_value:{metric_key}"
                                        )
                                        simulation_metric_limit_exceeded = True
                                        break
                                    if actual_numeric > limit_value:
                                        failures.append(
                                            f"{case_label}:simulation_artifact_metric_exceeded:"
                                            f"id={metric_key}:limit={limit_value}:actual={actual_numeric}"
                                        )
                                        simulation_metric_limit_exceeded = True
                                        break
                                if simulation_metric_limit_exceeded:
                                    continue
                                benchmark_case_status = "PASS"
                                observer_contract_id = str(
                                    metadata.get(
                                        "surface_pic_runtime_boundary_observer_contract_id", ""
                                    )
                                )
                                if (
                                    expected_boundary_observer_contract
                                    and observer_contract_id
                                    != expected_boundary_observer_contract
                                ):
                                    failures.append(
                                        f"{preset}:boundary_observer_contract_mismatch:"
                                        f"expected={expected_boundary_observer_contract}:actual={observer_contract_id}"
                                    )
                                else:
                                    observer_path = observer_sidecar(
                                        output_root, metadata, csv_path
                                    )
                                    if not observer_path.is_file():
                                        failures.append(
                                            f"{preset}:missing_boundary_observer_sidecar:{observer_path}"
                                        )
                                    else:
                                        observer_payload = load_json(observer_path)
                                        actual_schema = str(
                                            observer_payload.get("schema_version", "")
                                        )
                                        actual_contract_id = str(
                                            observer_payload.get("contract_id", "")
                                        )
                                        actual_runtime = str(
                                            observer_payload.get("surface_pic_runtime", "")
                                        )
                                        actual_shared_runtime_enabled = observer_payload.get(
                                            "shared_runtime_enabled"
                                        )

                                        if (
                                            actual_schema
                                            != "scdat.surface_shared_runtime_observer.v1"
                                        ):
                                            failures.append(
                                                f"{preset}:boundary_observer_schema_mismatch:"
                                                f"actual={actual_schema}"
                                            )
                                        elif (
                                            expected_boundary_observer_contract
                                            and actual_contract_id
                                            != expected_boundary_observer_contract
                                        ):
                                            failures.append(
                                                f"{preset}:boundary_observer_payload_contract_mismatch:"
                                                f"expected={expected_boundary_observer_contract}:actual={actual_contract_id}"
                                            )
                                        elif (
                                            expected_surface_pic_runtime
                                            and actual_runtime != expected_surface_pic_runtime
                                        ):
                                            failures.append(
                                                f"{preset}:surface_pic_runtime_mismatch:"
                                                f"expected={expected_surface_pic_runtime}:actual={actual_runtime}"
                                            )
                                        elif (
                                            expected_shared_runtime_enabled is not None
                                            and bool(actual_shared_runtime_enabled)
                                            != bool(expected_shared_runtime_enabled)
                                        ):
                                            failures.append(
                                                f"{preset}:shared_runtime_enabled_mismatch:"
                                                f"expected={bool(expected_shared_runtime_enabled)}:"
                                                f"actual={bool(actual_shared_runtime_enabled)}"
                                            )
                                        else:
                                            observer_status = "PASS"
                                            if not expected_consistency_contract:
                                                consistency_status = "PASS"
                                                status = "PASS"
                                            else:
                                                actual_consistency_contract = str(
                                                    metadata.get(
                                                        "surface_pic_runtime_consistency_contract_id",
                                                        "",
                                                    )
                                                )
                                                if (
                                                    actual_consistency_contract
                                                    != expected_consistency_contract
                                                ):
                                                    failures.append(
                                                        f"{preset}:surface_runtime_consistency_contract_mismatch:"
                                                        f"expected={expected_consistency_contract}:actual={actual_consistency_contract}"
                                                    )
                                                else:
                                                    consistency_path = consistency_sidecar(
                                                        output_root, metadata, csv_path
                                                    )
                                                    if not consistency_path.is_file():
                                                        failures.append(
                                                            f"{preset}:missing_surface_runtime_consistency_sidecar:{consistency_path}"
                                                        )
                                                    else:
                                                        consistency_payload = load_json(
                                                            consistency_path
                                                        )
                                                        actual_consistency_schema = str(
                                                            consistency_payload.get(
                                                                "schema_version", ""
                                                            )
                                                        )
                                                        actual_consistency_contract = str(
                                                            consistency_payload.get(
                                                                "contract_id", ""
                                                            )
                                                        )
                                                        shared_patch_node_count = int(
                                                            consistency_payload.get(
                                                                "shared_patch_node_count", 0
                                                            )
                                                        )
                                                        reference_spread_v = float(
                                                            consistency_payload.get(
                                                                "reference_potential_spread_v",
                                                                0.0,
                                                            )
                                                        )
                                                        patch_potential_spread_v = float(
                                                            consistency_payload.get(
                                                                "patch_potential_spread_v", 0.0
                                                            )
                                                        )
                                                        patch_potential_spread_reduction_ratio = float(
                                                            consistency_payload.get(
                                                                "patch_potential_spread_reduction_ratio", 0.0
                                                            )
                                                        )
                                                        sheath_charge_spread_c = float(
                                                            consistency_payload.get(
                                                                "sheath_charge_spread_c", 0.0
                                                            )
                                                        )
                                                        live_pic_net_current_spread_a_per_m2 = float(
                                                            consistency_payload.get(
                                                                "live_pic_net_current_spread_a_per_m2", 0.0
                                                            )
                                                        )
                                                        electron_calibration_factor_spread = float(
                                                            consistency_payload.get(
                                                                "electron_calibration_factor_spread", 0.0
                                                            )
                                                        )
                                                        consistency_pass = bool(
                                                            consistency_payload.get("consistency_pass")
                                                        )
                                                        shared_particle_pool_consistency_pass = bool(
                                                            consistency_payload.get(
                                                                "shared_particle_pool_consistency_pass", False
                                                            )
                                                        )
                                                        global_sheath_proxy_solve_active = bool(
                                                            consistency_payload.get(
                                                                "global_sheath_proxy_solve_active", False
                                                            )
                                                        )
                                                        shared_current_matrix_coupling_active = bool(
                                                            consistency_payload.get(
                                                                "shared_current_matrix_coupling_active", False
                                                            )
                                                        )
                                                        shared_current_matrix_coupling_offdiag_entry_count = int(
                                                            consistency_payload.get(
                                                                "shared_current_matrix_coupling_offdiag_entry_count",
                                                                0,
                                                            )
                                                        )
                                                        shared_global_coupled_solve_active = bool(
                                                            consistency_payload.get(
                                                                "shared_global_coupled_solve_active", False
                                                            )
                                                        )
                                                        shared_global_coupled_solve_iterations = int(
                                                            consistency_payload.get(
                                                                "shared_global_coupled_solve_iterations", 0
                                                            )
                                                        )
                                                        shared_global_coupled_solve_converged = bool(
                                                            consistency_payload.get(
                                                                "shared_global_coupled_solve_converged", False
                                                            )
                                                        )
                                                        shared_global_coupled_solve_max_delta_v = float(
                                                            consistency_payload.get(
                                                                "shared_global_coupled_solve_max_delta_v", 0.0
                                                            )
                                                        )
                                                        normal_electric_field_spread_v_per_m = float(
                                                            consistency_payload.get(
                                                                "normal_electric_field_spread_v_per_m", 0.0
                                                            )
                                                        )
                                                        local_charge_density_spread_c_per_m3 = float(
                                                            consistency_payload.get(
                                                                "local_charge_density_spread_c_per_m3", 0.0
                                                            )
                                                        )
                                                        shared_live_pic_coupled_refresh_active = bool(
                                                            consistency_payload.get(
                                                                "shared_live_pic_coupled_refresh_active", False
                                                            )
                                                        )
                                                        shared_live_pic_coupled_refresh_count = int(
                                                            consistency_payload.get(
                                                                "shared_live_pic_coupled_refresh_count", 0
                                                            )
                                                        )
                                                        shared_particle_transport_coupling_active = bool(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_coupling_active",
                                                                False,
                                                            )
                                                        )
                                                        shared_particle_transport_offdiag_entry_count = int(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_offdiag_entry_count",
                                                                0,
                                                            )
                                                        )
                                                        shared_particle_transport_conservation_error_a_per_v = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_conservation_error_a_per_v",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_charge_c = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_charge_c", 0.0
                                                            )
                                                        )
                                                        shared_particle_transport_reference_shift_v = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_reference_shift_v",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_distribution_active = bool(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_distribution_active",
                                                                False,
                                                            )
                                                        )
                                                        shared_particle_transport_distribution_charge_spread_c = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_distribution_charge_spread_c",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_distribution_reference_shift_spread_v = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_distribution_reference_shift_spread_v",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_distribution_conservation_error_c = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_distribution_conservation_error_c",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_exchange_active = bool(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_exchange_active",
                                                                False,
                                                            )
                                                        )
                                                        shared_particle_transport_exchange_flux_spread_a = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_exchange_flux_spread_a",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_exchange_flux_conservation_error_a = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_exchange_flux_conservation_error_a",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_edge_domain_active = bool(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_edge_domain_active",
                                                                False,
                                                            )
                                                        )
                                                        shared_particle_transport_edge_operator_active = bool(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_edge_operator_active",
                                                                False,
                                                            )
                                                        )
                                                        shared_particle_transport_edge_operator_total_abs_drive_charge_c = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_edge_operator_total_abs_drive_charge_c",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_edge_graph_operator_iterations = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_edge_graph_operator_iterations",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_edge_graph_operator_converged = bool(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_edge_graph_operator_converged",
                                                                False,
                                                            )
                                                        )
                                                        shared_particle_transport_edge_graph_operator_max_balance_residual_c = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_edge_graph_operator_max_balance_residual_c",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_edge_graph_operator_branch_graph_active = bool(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_edge_graph_operator_branch_graph_active",
                                                                False,
                                                            )
                                                        )
                                                        shared_particle_transport_edge_graph_operator_branch_graph_edge_count = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_edge_graph_operator_branch_graph_edge_count",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_edge_graph_operator_branch_graph_pair_count = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_edge_graph_operator_branch_graph_pair_count",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_edge_graph_operator_total_conductance_weight_f = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_edge_graph_operator_total_conductance_weight_f",
                                                                0.0,
                                                            )
                                                        )
                                                        shared_particle_transport_edge_graph_operator_max_node_preconditioner = float(
                                                            consistency_payload.get(
                                                                "shared_particle_transport_edge_graph_operator_max_node_preconditioner",
                                                                0.0,
                                                            )
                                                        )
                                                        if (
                                                            actual_consistency_schema
                                                            != "scdat.surface_shared_runtime_consistency.v1"
                                                        ):
                                                            failures.append(
                                                                f"{preset}:surface_runtime_consistency_schema_mismatch:"
                                                                f"actual={actual_consistency_schema}"
                                                            )
                                                        elif (
                                                            actual_consistency_contract
                                                            != expected_consistency_contract
                                                        ):
                                                            failures.append(
                                                                f"{preset}:surface_runtime_consistency_payload_contract_mismatch:"
                                                                f"expected={expected_consistency_contract}:actual={actual_consistency_contract}"
                                                            )
                                                        elif (
                                                            require_surface_runtime_consistency_pass
                                                            and not consistency_pass
                                                        ):
                                                            failures.append(
                                                                f"{preset}:surface_runtime_consistency_failed"
                                                            )
                                                        elif not shared_particle_pool_consistency_pass:
                                                            failures.append(
                                                                f"{preset}:shared_particle_pool_consistency_failed"
                                                            )
                                                        elif (
                                                            require_global_sheath_proxy_solve_active
                                                            and not global_sheath_proxy_solve_active
                                                        ):
                                                            failures.append(
                                                                f"{preset}:global_sheath_proxy_solve_not_active"
                                                            )
                                                        elif (
                                                            require_shared_current_matrix_coupling_active
                                                            and not shared_current_matrix_coupling_active
                                                        ):
                                                            failures.append(
                                                                f"{preset}:shared_current_matrix_coupling_not_active"
                                                            )
                                                        elif (
                                                            min_shared_current_matrix_coupling_offdiag_entry_count
                                                            is not None
                                                            and shared_current_matrix_coupling_offdiag_entry_count
                                                            < int(
                                                                min_shared_current_matrix_coupling_offdiag_entry_count
                                                            )
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_current_matrix_coupling_offdiag_entry_count_too_small:"
                                                                f"minimum={int(min_shared_current_matrix_coupling_offdiag_entry_count)}:"
                                                                f"actual={shared_current_matrix_coupling_offdiag_entry_count}"
                                                            )
                                                        elif (
                                                            require_shared_global_coupled_solve_active
                                                            and not shared_global_coupled_solve_active
                                                        ):
                                                            failures.append(
                                                                f"{preset}:shared_global_coupled_solve_not_active"
                                                            )
                                                        elif (
                                                            require_shared_global_coupled_solve_converged
                                                            and not shared_global_coupled_solve_converged
                                                        ):
                                                            failures.append(
                                                                f"{preset}:shared_global_coupled_solve_not_converged"
                                                            )
                                                        elif (
                                                            min_shared_global_coupled_solve_iterations
                                                            is not None
                                                            and shared_global_coupled_solve_iterations
                                                            < int(min_shared_global_coupled_solve_iterations)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_global_coupled_solve_iterations_too_small:"
                                                                f"minimum={int(min_shared_global_coupled_solve_iterations)}:"
                                                                f"actual={shared_global_coupled_solve_iterations}"
                                                            )
                                                        elif (
                                                            max_shared_global_coupled_solve_max_delta_v
                                                            is not None
                                                            and shared_global_coupled_solve_max_delta_v
                                                            > float(max_shared_global_coupled_solve_max_delta_v)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_global_coupled_solve_max_delta_exceeded:"
                                                                f"limit={float(max_shared_global_coupled_solve_max_delta_v)}:"
                                                                f"actual={shared_global_coupled_solve_max_delta_v}"
                                                            )
                                                        elif (
                                                            max_normal_electric_field_spread_v_per_m is not None
                                                            and normal_electric_field_spread_v_per_m
                                                            > float(max_normal_electric_field_spread_v_per_m)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:normal_electric_field_spread_exceeded:"
                                                                f"limit={float(max_normal_electric_field_spread_v_per_m)}:"
                                                                f"actual={normal_electric_field_spread_v_per_m}"
                                                            )
                                                        elif (
                                                            max_local_charge_density_spread_c_per_m3 is not None
                                                            and local_charge_density_spread_c_per_m3
                                                            > float(max_local_charge_density_spread_c_per_m3)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:local_charge_density_spread_exceeded:"
                                                                f"limit={float(max_local_charge_density_spread_c_per_m3)}:"
                                                                f"actual={local_charge_density_spread_c_per_m3}"
                                                            )
                                                        elif (
                                                            require_shared_live_pic_coupled_refresh_active
                                                            and not shared_live_pic_coupled_refresh_active
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_live_pic_coupled_refresh_not_active"
                                                            )
                                                        elif (
                                                            min_shared_live_pic_coupled_refresh_count
                                                            is not None
                                                            and shared_live_pic_coupled_refresh_count
                                                            < int(min_shared_live_pic_coupled_refresh_count)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_live_pic_coupled_refresh_count_too_small:"
                                                                f"minimum={int(min_shared_live_pic_coupled_refresh_count)}:"
                                                                f"actual={shared_live_pic_coupled_refresh_count}"
                                                            )
                                                        elif (
                                                            require_shared_particle_transport_coupling_active
                                                            and not shared_particle_transport_coupling_active
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_coupling_not_active"
                                                            )
                                                        elif (
                                                            min_shared_particle_transport_offdiag_entry_count
                                                            is not None
                                                            and shared_particle_transport_offdiag_entry_count
                                                            < int(min_shared_particle_transport_offdiag_entry_count)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_offdiag_entry_count_too_small:"
                                                                f"minimum={int(min_shared_particle_transport_offdiag_entry_count)}:"
                                                                f"actual={shared_particle_transport_offdiag_entry_count}"
                                                            )
                                                        elif (
                                                            max_shared_particle_transport_conservation_error_a_per_v
                                                            is not None
                                                            and shared_particle_transport_conservation_error_a_per_v
                                                            > float(max_shared_particle_transport_conservation_error_a_per_v)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_conservation_error_exceeded:"
                                                                f"limit={float(max_shared_particle_transport_conservation_error_a_per_v)}:"
                                                                f"actual={shared_particle_transport_conservation_error_a_per_v}"
                                                            )
                                                        elif (
                                                            min_abs_shared_particle_transport_reference_shift_v
                                                            is not None
                                                            and abs(shared_particle_transport_reference_shift_v)
                                                            < float(min_abs_shared_particle_transport_reference_shift_v)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_reference_shift_too_small:"
                                                                f"minimum_abs={float(min_abs_shared_particle_transport_reference_shift_v)}:"
                                                                f"actual={shared_particle_transport_reference_shift_v}"
                                                            )
                                                        elif (
                                                            min_abs_shared_particle_transport_charge_c is not None
                                                            and abs(shared_particle_transport_charge_c)
                                                            < float(min_abs_shared_particle_transport_charge_c)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_charge_too_small:"
                                                                f"minimum_abs={float(min_abs_shared_particle_transport_charge_c)}:"
                                                                f"actual={shared_particle_transport_charge_c}"
                                                            )
                                                        elif (
                                                            require_shared_particle_transport_distribution_active
                                                            and not shared_particle_transport_distribution_active
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_distribution_not_active"
                                                            )
                                                        elif (
                                                            min_shared_particle_transport_distribution_charge_spread_c
                                                            is not None
                                                            and shared_particle_transport_distribution_charge_spread_c
                                                            < float(min_shared_particle_transport_distribution_charge_spread_c)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_distribution_charge_spread_too_small:"
                                                                f"minimum={float(min_shared_particle_transport_distribution_charge_spread_c)}:"
                                                                f"actual={shared_particle_transport_distribution_charge_spread_c}"
                                                            )
                                                        elif (
                                                            min_shared_particle_transport_distribution_reference_shift_spread_v
                                                            is not None
                                                            and shared_particle_transport_distribution_reference_shift_spread_v
                                                            < float(min_shared_particle_transport_distribution_reference_shift_spread_v)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_distribution_reference_shift_spread_too_small:"
                                                                f"minimum={float(min_shared_particle_transport_distribution_reference_shift_spread_v)}:"
                                                                f"actual={shared_particle_transport_distribution_reference_shift_spread_v}"
                                                            )
                                                        elif (
                                                            max_shared_particle_transport_distribution_conservation_error_c
                                                            is not None
                                                            and shared_particle_transport_distribution_conservation_error_c
                                                            > float(max_shared_particle_transport_distribution_conservation_error_c)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_distribution_conservation_error_exceeded:"
                                                                f"limit={float(max_shared_particle_transport_distribution_conservation_error_c)}:"
                                                                f"actual={shared_particle_transport_distribution_conservation_error_c}"
                                                            )
                                                        elif (
                                                            require_shared_particle_transport_exchange_active
                                                            and not shared_particle_transport_exchange_active
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_exchange_not_active"
                                                            )
                                                        elif (
                                                            min_shared_particle_transport_exchange_flux_spread_a
                                                            is not None
                                                            and shared_particle_transport_exchange_flux_spread_a
                                                            < float(min_shared_particle_transport_exchange_flux_spread_a)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_exchange_flux_spread_too_small:"
                                                                f"minimum={float(min_shared_particle_transport_exchange_flux_spread_a)}:"
                                                                f"actual={shared_particle_transport_exchange_flux_spread_a}"
                                                            )
                                                        elif (
                                                            max_shared_particle_transport_exchange_flux_conservation_error_a
                                                            is not None
                                                            and shared_particle_transport_exchange_flux_conservation_error_a
                                                            > float(max_shared_particle_transport_exchange_flux_conservation_error_a)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_exchange_flux_conservation_error_exceeded:"
                                                                f"limit={float(max_shared_particle_transport_exchange_flux_conservation_error_a)}:"
                                                                f"actual={shared_particle_transport_exchange_flux_conservation_error_a}"
                                                            )
                                                        elif (
                                                            require_shared_particle_transport_edge_operator_active
                                                            and not shared_particle_transport_edge_operator_active
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_edge_operator_not_active"
                                                            )
                                                        elif (
                                                            min_shared_particle_transport_edge_operator_total_abs_drive_charge_c
                                                            is not None
                                                            and shared_particle_transport_edge_operator_total_abs_drive_charge_c
                                                            < float(
                                                                min_shared_particle_transport_edge_operator_total_abs_drive_charge_c
                                                            )
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_edge_operator_total_abs_drive_charge_too_small:"
                                                                f"minimum={float(min_shared_particle_transport_edge_operator_total_abs_drive_charge_c)}:"
                                                                f"actual={shared_particle_transport_edge_operator_total_abs_drive_charge_c}"
                                                            )
                                                        elif (
                                                            require_shared_particle_transport_edge_graph_operator_converged
                                                            and not shared_particle_transport_edge_graph_operator_converged
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_edge_graph_operator_not_converged"
                                                            )
                                                        elif (
                                                            min_shared_particle_transport_edge_graph_operator_iterations
                                                            is not None
                                                            and shared_particle_transport_edge_graph_operator_iterations
                                                            < float(
                                                                min_shared_particle_transport_edge_graph_operator_iterations
                                                            )
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_edge_graph_operator_iterations_too_small:"
                                                                f"minimum={float(min_shared_particle_transport_edge_graph_operator_iterations)}:"
                                                                f"actual={shared_particle_transport_edge_graph_operator_iterations}"
                                                            )
                                                        elif (
                                                            max_shared_particle_transport_edge_graph_operator_max_balance_residual_c
                                                            is not None
                                                            and shared_particle_transport_edge_graph_operator_max_balance_residual_c
                                                            > float(
                                                                max_shared_particle_transport_edge_graph_operator_max_balance_residual_c
                                                            )
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_edge_graph_operator_max_balance_residual_exceeded:"
                                                                f"limit={float(max_shared_particle_transport_edge_graph_operator_max_balance_residual_c)}:"
                                                                f"actual={shared_particle_transport_edge_graph_operator_max_balance_residual_c}"
                                                            )
                                                        elif (
                                                            require_shared_particle_transport_edge_graph_operator_branch_graph_active
                                                            and not shared_particle_transport_edge_graph_operator_branch_graph_active
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_edge_graph_operator_branch_graph_not_active"
                                                            )
                                                        elif (
                                                            min_shared_particle_transport_edge_graph_operator_branch_graph_edge_count
                                                            is not None
                                                            and shared_particle_transport_edge_graph_operator_branch_graph_edge_count
                                                            < float(
                                                                min_shared_particle_transport_edge_graph_operator_branch_graph_edge_count
                                                            )
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_edge_graph_operator_branch_graph_edge_count_too_small:"
                                                                f"minimum={float(min_shared_particle_transport_edge_graph_operator_branch_graph_edge_count)}:"
                                                                f"actual={shared_particle_transport_edge_graph_operator_branch_graph_edge_count}"
                                                            )
                                                        elif (
                                                            min_shared_particle_transport_edge_graph_operator_branch_graph_pair_count
                                                            is not None
                                                            and shared_particle_transport_edge_graph_operator_branch_graph_pair_count
                                                            < float(
                                                                min_shared_particle_transport_edge_graph_operator_branch_graph_pair_count
                                                            )
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_edge_graph_operator_branch_graph_pair_count_too_small:"
                                                                f"minimum={float(min_shared_particle_transport_edge_graph_operator_branch_graph_pair_count)}:"
                                                                f"actual={shared_particle_transport_edge_graph_operator_branch_graph_pair_count}"
                                                            )
                                                        elif (
                                                            min_shared_particle_transport_edge_graph_operator_total_conductance_weight_f
                                                            is not None
                                                            and shared_particle_transport_edge_graph_operator_total_conductance_weight_f
                                                            < float(
                                                                min_shared_particle_transport_edge_graph_operator_total_conductance_weight_f
                                                            )
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_edge_graph_operator_total_conductance_weight_too_small:"
                                                                f"minimum={float(min_shared_particle_transport_edge_graph_operator_total_conductance_weight_f)}:"
                                                                f"actual={shared_particle_transport_edge_graph_operator_total_conductance_weight_f}"
                                                            )
                                                        elif (
                                                            min_shared_particle_transport_edge_graph_operator_max_node_preconditioner
                                                            is not None
                                                            and shared_particle_transport_edge_graph_operator_max_node_preconditioner
                                                            < float(
                                                                min_shared_particle_transport_edge_graph_operator_max_node_preconditioner
                                                            )
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_particle_transport_edge_graph_operator_max_node_preconditioner_too_small:"
                                                                f"minimum={float(min_shared_particle_transport_edge_graph_operator_max_node_preconditioner)}:"
                                                                f"actual={shared_particle_transport_edge_graph_operator_max_node_preconditioner}"
                                                            )
                                                        elif (
                                                            min_shared_patch_node_count is not None
                                                            and shared_patch_node_count
                                                            < int(min_shared_patch_node_count)
                                                        ):
                                                            failures.append(
                                                                f"{case_label}:shared_patch_node_count_too_small:"
                                                                f"minimum={int(min_shared_patch_node_count)}:"
                                                                f"actual={shared_patch_node_count}"
                                                            )
                                                        elif (
                                                            max_reference_spread_v is not None
                                                            and reference_spread_v
                                                            > float(max_reference_spread_v)
                                                        ):
                                                            failures.append(
                                                                f"{preset}:reference_spread_exceeded:"
                                                                f"limit={float(max_reference_spread_v)}:actual={reference_spread_v}"
                                                            )
                                                        elif (
                                                            min_patch_potential_spread_reduction_ratio
                                                            is not None
                                                            and patch_potential_spread_reduction_ratio
                                                            < float(min_patch_potential_spread_reduction_ratio)
                                                        ):
                                                            failures.append(
                                                                f"{preset}:patch_potential_spread_reduction_ratio_too_small:"
                                                                f"minimum={float(min_patch_potential_spread_reduction_ratio)}:"
                                                                f"actual={patch_potential_spread_reduction_ratio}"
                                                            )
                                                        elif (
                                                            max_patch_potential_spread_v is not None
                                                            and patch_potential_spread_v
                                                            > float(max_patch_potential_spread_v)
                                                        ):
                                                            failures.append(
                                                                f"{preset}:patch_potential_spread_exceeded:"
                                                                f"limit={float(max_patch_potential_spread_v)}:"
                                                                f"actual={patch_potential_spread_v}"
                                                            )
                                                        elif (
                                                            max_sheath_charge_spread_c is not None
                                                            and sheath_charge_spread_c
                                                            > float(max_sheath_charge_spread_c)
                                                        ):
                                                            failures.append(
                                                                f"{preset}:sheath_charge_spread_exceeded:"
                                                                f"limit={float(max_sheath_charge_spread_c)}:actual={sheath_charge_spread_c}"
                                                            )
                                                        elif (
                                                            max_live_pic_net_current_spread_a_per_m2 is not None
                                                            and live_pic_net_current_spread_a_per_m2
                                                            > float(max_live_pic_net_current_spread_a_per_m2)
                                                        ):
                                                            failures.append(
                                                                f"{preset}:live_pic_net_current_spread_exceeded:"
                                                                f"limit={float(max_live_pic_net_current_spread_a_per_m2)}:"
                                                                f"actual={live_pic_net_current_spread_a_per_m2}"
                                                            )
                                                        elif (
                                                            max_electron_calibration_factor_spread is not None
                                                            and electron_calibration_factor_spread
                                                            > float(max_electron_calibration_factor_spread)
                                                        ):
                                                            failures.append(
                                                                f"{preset}:electron_calibration_factor_spread_exceeded:"
                                                                f"limit={float(max_electron_calibration_factor_spread)}:"
                                                                f"actual={electron_calibration_factor_spread}"
                                                            )
                                                        else:
                                                            consistency_status = "PASS"
                                                            status = "PASS"
        if status == "PASS" and expected_transport_domain_contract:
            actual_transport_domain_contract = str(
                metadata.get(
                    "surface_pic_runtime_shared_particle_transport_domain_contract_id", ""
                )
            )
            if actual_transport_domain_contract != expected_transport_domain_contract:
                failures.append(
                    f"{case_label}:shared_particle_transport_domain_contract_mismatch:"
                    f"expected={expected_transport_domain_contract}:actual={actual_transport_domain_contract}"
                )
                status = "FAIL"
            else:
                domain_path = transport_domain_sidecar(output_root, metadata, csv_path)
                if not domain_path.is_file():
                    failures.append(
                        f"{case_label}:missing_shared_particle_transport_domain_sidecar:{domain_path}"
                    )
                    status = "FAIL"
                else:
                    domain_payload = load_json(domain_path)
                    actual_domain_schema = str(domain_payload.get("schema_version", ""))
                    actual_domain_contract = str(domain_payload.get("contract_id", ""))
                    shared_particle_transport_domain_active = bool(
                        domain_payload.get("shared_particle_transport_domain_active", False)
                    )
                    transport_domain_runtime_state_backed = domain_payload.get(
                        "runtime_state_backed", False
                    )
                    shared_particle_transport_bookkeeping_mode = str(
                        domain_payload.get(
                            "shared_particle_transport_bookkeeping_mode", ""
                        )
                    )
                    transport_domain_node_count = int(
                        domain_payload.get("shared_patch_node_count", 0)
                    )
                    transport_domain_edge_count = int(
                        domain_payload.get("exchange_edge_count", 0)
                    )
                    transport_domain_exchange_flux_conservation_error_a = float(
                        domain_payload.get("exchange_flux_conservation_error_a", 0.0)
                    )
                    transport_domain_edge_charge_total_abs_c = float(
                        domain_payload.get("edge_charge_total_abs_c", 0.0)
                    )
                    transport_domain_edge_target_charge_total_abs_c = float(
                        domain_payload.get("edge_target_charge_total_abs_c", 0.0)
                    )
                    transport_domain_edge_charge_conservation_error_c = float(
                        domain_payload.get("edge_charge_conservation_error_c", 0.0)
                    )
                    transport_domain_edge_operator_drive_total_abs_c = float(
                        domain_payload.get("edge_operator_drive_total_abs_c", 0.0)
                    )
                    transport_domain_edge_operator_drive_conservation_error_c = float(
                        domain_payload.get(
                            "edge_operator_drive_conservation_error_c", 0.0
                        )
                    )
                    transport_domain_edge_graph_operator_iterations = float(
                        domain_payload.get("edge_graph_operator_iterations", 0.0)
                    )
                    transport_domain_edge_graph_operator_converged = bool(
                        domain_payload.get("edge_graph_operator_converged", False)
                    )
                    transport_domain_edge_graph_operator_max_balance_residual_c = float(
                        domain_payload.get(
                            "edge_graph_operator_max_balance_residual_c", 0.0
                        )
                    )
                    transport_domain_edge_graph_operator_branch_graph_active = bool(
                        domain_payload.get("edge_graph_operator_branch_graph_active", False)
                    )
                    transport_domain_edge_graph_operator_branch_graph_edge_count = float(
                        domain_payload.get(
                            "edge_graph_operator_branch_graph_edge_count", 0.0
                        )
                    )
                    transport_domain_edge_graph_operator_branch_graph_pair_count = float(
                        domain_payload.get(
                            "edge_graph_operator_branch_graph_pair_count", 0.0
                        )
                    )
                    transport_domain_edge_graph_operator_total_conductance_weight_f = float(
                        domain_payload.get(
                            "edge_graph_operator_total_conductance_weight_f", 0.0
                        )
                    )
                    transport_domain_edge_graph_operator_max_node_preconditioner = float(
                        domain_payload.get(
                            "edge_graph_operator_max_node_preconditioner", 0.0
                        )
                    )
                    domain_nodes = domain_payload.get("nodes", [])
                    domain_edges = domain_payload.get("exchange_edges", [])

                    if (
                        actual_domain_schema
                        != "scdat.surface_shared_particle_transport_domain.v1"
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_schema_mismatch:"
                            f"actual={actual_domain_schema}"
                        )
                        status = "FAIL"
                    elif actual_domain_contract != expected_transport_domain_contract:
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_payload_contract_mismatch:"
                            f"expected={expected_transport_domain_contract}:actual={actual_domain_contract}"
                        )
                        status = "FAIL"
                    elif (
                        require_shared_particle_transport_domain_runtime_state_backed is True
                        and transport_domain_runtime_state_backed is not True
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_not_runtime_state_backed"
                        )
                        status = "FAIL"
                    elif (
                        expected_shared_particle_transport_bookkeeping_mode
                        and shared_particle_transport_bookkeeping_mode
                        != expected_shared_particle_transport_bookkeeping_mode
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_bookkeeping_mode_mismatch:"
                            f"expected={expected_shared_particle_transport_bookkeeping_mode}:"
                            f"actual={shared_particle_transport_bookkeeping_mode}"
                        )
                        status = "FAIL"
                    elif (
                        require_shared_particle_transport_domain_active
                        and not shared_particle_transport_domain_active
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_not_active"
                        )
                        status = "FAIL"
                    elif (
                        min_shared_particle_transport_domain_node_count is not None
                        and transport_domain_node_count
                        < int(min_shared_particle_transport_domain_node_count)
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_node_count_too_small:"
                            f"minimum={int(min_shared_particle_transport_domain_node_count)}:"
                            f"actual={transport_domain_node_count}"
                        )
                        status = "FAIL"
                    elif not isinstance(domain_nodes, list) or len(domain_nodes) != transport_domain_node_count:
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_nodes_mismatch:"
                            f"expected={transport_domain_node_count}:actual={len(domain_nodes) if isinstance(domain_nodes, list) else 'invalid'}"
                        )
                        status = "FAIL"
                    elif (
                        min_shared_particle_transport_domain_edge_count is not None
                        and transport_domain_edge_count
                        < int(min_shared_particle_transport_domain_edge_count)
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_count_too_small:"
                            f"minimum={int(min_shared_particle_transport_domain_edge_count)}:"
                            f"actual={transport_domain_edge_count}"
                        )
                        status = "FAIL"
                    elif not isinstance(domain_edges, list) or len(domain_edges) != transport_domain_edge_count:
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edges_mismatch:"
                            f"expected={transport_domain_edge_count}:actual={len(domain_edges) if isinstance(domain_edges, list) else 'invalid'}"
                        )
                        status = "FAIL"
                    elif (
                        max_shared_particle_transport_domain_exchange_flux_conservation_error_a
                        is not None
                        and transport_domain_exchange_flux_conservation_error_a
                        > float(
                            max_shared_particle_transport_domain_exchange_flux_conservation_error_a
                        )
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_exchange_flux_conservation_error_exceeded:"
                            f"limit={float(max_shared_particle_transport_domain_exchange_flux_conservation_error_a)}:"
                            f"actual={transport_domain_exchange_flux_conservation_error_a}"
                        )
                        status = "FAIL"
                    elif (
                        min_shared_particle_transport_domain_edge_charge_total_abs_c
                        is not None
                        and transport_domain_edge_charge_total_abs_c
                        < float(min_shared_particle_transport_domain_edge_charge_total_abs_c)
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_charge_total_abs_too_small:"
                            f"minimum={float(min_shared_particle_transport_domain_edge_charge_total_abs_c)}:"
                            f"actual={transport_domain_edge_charge_total_abs_c}"
                        )
                        status = "FAIL"
                    elif (
                        max_shared_particle_transport_domain_edge_charge_conservation_error_c
                        is not None
                        and transport_domain_edge_charge_conservation_error_c
                        > float(
                            max_shared_particle_transport_domain_edge_charge_conservation_error_c
                        )
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_charge_conservation_error_exceeded:"
                            f"limit={float(max_shared_particle_transport_domain_edge_charge_conservation_error_c)}:"
                            f"actual={transport_domain_edge_charge_conservation_error_c}"
                        )
                        status = "FAIL"
                    elif (
                        min_shared_particle_transport_domain_edge_operator_drive_total_abs_c
                        is not None
                        and transport_domain_edge_operator_drive_total_abs_c
                        < float(
                            min_shared_particle_transport_domain_edge_operator_drive_total_abs_c
                        )
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_operator_drive_total_abs_too_small:"
                            f"minimum={float(min_shared_particle_transport_domain_edge_operator_drive_total_abs_c)}:"
                            f"actual={transport_domain_edge_operator_drive_total_abs_c}"
                        )
                        status = "FAIL"
                    elif (
                        max_shared_particle_transport_domain_edge_operator_drive_conservation_error_c
                        is not None
                        and transport_domain_edge_operator_drive_conservation_error_c
                        > float(
                            max_shared_particle_transport_domain_edge_operator_drive_conservation_error_c
                        )
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_operator_drive_conservation_error_exceeded:"
                            f"limit={float(max_shared_particle_transport_domain_edge_operator_drive_conservation_error_c)}:"
                            f"actual={transport_domain_edge_operator_drive_conservation_error_c}"
                        )
                        status = "FAIL"
                    elif (
                        require_shared_particle_transport_domain_edge_graph_operator_converged
                        and not transport_domain_edge_graph_operator_converged
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_graph_operator_not_converged"
                        )
                        status = "FAIL"
                    elif (
                        min_shared_particle_transport_domain_edge_graph_operator_iterations
                        is not None
                        and transport_domain_edge_graph_operator_iterations
                        < float(
                            min_shared_particle_transport_domain_edge_graph_operator_iterations
                        )
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_graph_operator_iterations_too_small:"
                            f"minimum={float(min_shared_particle_transport_domain_edge_graph_operator_iterations)}:"
                            f"actual={transport_domain_edge_graph_operator_iterations}"
                        )
                        status = "FAIL"
                    elif (
                        max_shared_particle_transport_domain_edge_graph_operator_max_balance_residual_c
                        is not None
                        and transport_domain_edge_graph_operator_max_balance_residual_c
                        > float(
                            max_shared_particle_transport_domain_edge_graph_operator_max_balance_residual_c
                        )
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_graph_operator_max_balance_residual_exceeded:"
                            f"limit={float(max_shared_particle_transport_domain_edge_graph_operator_max_balance_residual_c)}:"
                            f"actual={transport_domain_edge_graph_operator_max_balance_residual_c}"
                        )
                        status = "FAIL"
                    elif (
                        require_shared_particle_transport_domain_edge_graph_operator_branch_graph_active
                        and not transport_domain_edge_graph_operator_branch_graph_active
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_graph_operator_branch_graph_not_active"
                        )
                        status = "FAIL"
                    elif (
                        min_shared_particle_transport_domain_edge_graph_operator_branch_graph_edge_count
                        is not None
                        and transport_domain_edge_graph_operator_branch_graph_edge_count
                        < float(
                            min_shared_particle_transport_domain_edge_graph_operator_branch_graph_edge_count
                        )
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_graph_operator_branch_graph_edge_count_too_small:"
                            f"minimum={float(min_shared_particle_transport_domain_edge_graph_operator_branch_graph_edge_count)}:"
                            f"actual={transport_domain_edge_graph_operator_branch_graph_edge_count}"
                        )
                        status = "FAIL"
                    elif (
                        min_shared_particle_transport_domain_edge_graph_operator_branch_graph_pair_count
                        is not None
                        and transport_domain_edge_graph_operator_branch_graph_pair_count
                        < float(
                            min_shared_particle_transport_domain_edge_graph_operator_branch_graph_pair_count
                        )
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_graph_operator_branch_graph_pair_count_too_small:"
                            f"minimum={float(min_shared_particle_transport_domain_edge_graph_operator_branch_graph_pair_count)}:"
                            f"actual={transport_domain_edge_graph_operator_branch_graph_pair_count}"
                        )
                        status = "FAIL"
                    elif (
                        min_shared_particle_transport_domain_edge_graph_operator_total_conductance_weight_f
                        is not None
                        and transport_domain_edge_graph_operator_total_conductance_weight_f
                        < float(
                            min_shared_particle_transport_domain_edge_graph_operator_total_conductance_weight_f
                        )
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_graph_operator_total_conductance_weight_too_small:"
                            f"minimum={float(min_shared_particle_transport_domain_edge_graph_operator_total_conductance_weight_f)}:"
                            f"actual={transport_domain_edge_graph_operator_total_conductance_weight_f}"
                        )
                        status = "FAIL"
                    elif (
                        min_shared_particle_transport_domain_edge_graph_operator_max_node_preconditioner
                        is not None
                        and transport_domain_edge_graph_operator_max_node_preconditioner
                        < float(
                            min_shared_particle_transport_domain_edge_graph_operator_max_node_preconditioner
                        )
                    ):
                        failures.append(
                            f"{case_label}:shared_particle_transport_domain_edge_graph_operator_max_node_preconditioner_too_small:"
                            f"minimum={float(min_shared_particle_transport_domain_edge_graph_operator_max_node_preconditioner)}:"
                            f"actual={transport_domain_edge_graph_operator_max_node_preconditioner}"
                        )
                        status = "FAIL"
                    else:
                        transport_domain_status = "PASS"

        if status == "PASS" and expected_global_particle_domain_contract:
            actual_global_particle_domain_contract = str(
                metadata.get("surface_pic_runtime_global_particle_domain_contract_id", "")
            )
            if actual_global_particle_domain_contract != expected_global_particle_domain_contract:
                failures.append(
                    f"{case_label}:global_particle_domain_contract_mismatch:"
                    f"expected={expected_global_particle_domain_contract}:actual={actual_global_particle_domain_contract}"
                )
                status = "FAIL"
            else:
                global_particle_path = global_particle_domain_sidecar(
                    output_root, metadata, csv_path
                )
                if not global_particle_path.is_file():
                    failures.append(
                        f"{case_label}:missing_global_particle_domain_sidecar:{global_particle_path}"
                    )
                    status = "FAIL"
                else:
                    global_particle_payload = load_json(global_particle_path)
                    global_particle_domain_active = bool(
                        global_particle_payload.get("global_particle_domain_active", False)
                    )
                    global_particle_bookkeeping_mode = str(
                        global_particle_payload.get("global_particle_bookkeeping_mode", "")
                    )
                    global_particle_runtime_state_backed = global_particle_payload.get(
                        "runtime_state_backed", False
                    )
                    global_particle_domain_node_count = int(
                        global_particle_payload.get("shared_patch_node_count", 0)
                    )
                    global_particle_domain_edge_count = int(
                        global_particle_payload.get("domain_edge_count", 0)
                    )
                    global_particle_charge_conservation_error_c = float(
                        global_particle_payload.get(
                            "global_particle_charge_conservation_error_c", 0.0
                        )
                    )
                    global_particle_flux_conservation_error_a = float(
                        global_particle_payload.get(
                            "global_particle_flux_conservation_error_a", 0.0
                        )
                    )
                    global_particle_edge_charge_total_abs_c = float(
                        global_particle_payload.get(
                            "global_particle_edge_charge_total_abs_c", 0.0
                        )
                    )
                    global_particle_edge_target_charge_total_abs_c = float(
                        global_particle_payload.get(
                            "global_particle_edge_target_charge_total_abs_c", 0.0
                        )
                    )
                    global_particle_edge_operator_drive_total_abs_c = float(
                        global_particle_payload.get(
                            "global_particle_edge_operator_drive_total_abs_c", 0.0
                        )
                    )
                    global_particle_edge_conductance_total_s = float(
                        global_particle_payload.get(
                            "global_particle_edge_conductance_total_s", 0.0
                        )
                    )
                    global_particle_edge_graph_operator_converged = bool(
                        global_particle_payload.get("edge_graph_operator_converged", False)
                    )
                    global_particle_nodes = global_particle_payload.get("nodes", [])
                    global_particle_edges = global_particle_payload.get("domain_edges", [])
                    if (
                        str(global_particle_payload.get("schema_version", ""))
                        != "scdat.surface_global_particle_domain.v1"
                    ):
                        failures.append(
                            f"{case_label}:global_particle_domain_schema_mismatch:"
                            f"actual={global_particle_payload.get('schema_version')}"
                        )
                        status = "FAIL"
                    elif str(global_particle_payload.get("contract_id", "")) != expected_global_particle_domain_contract:
                        failures.append(
                            f"{case_label}:global_particle_domain_payload_contract_mismatch:"
                            f"expected={expected_global_particle_domain_contract}:actual={global_particle_payload.get('contract_id')}"
                        )
                        status = "FAIL"
                    elif require_global_particle_domain_active and not global_particle_domain_active:
                        failures.append(f"{case_label}:global_particle_domain_not_active")
                        status = "FAIL"
                    elif (
                        expected_global_particle_bookkeeping_mode
                        and global_particle_bookkeeping_mode
                        != expected_global_particle_bookkeeping_mode
                    ):
                        failures.append(
                            f"{case_label}:global_particle_bookkeeping_mode_mismatch:"
                            f"expected={expected_global_particle_bookkeeping_mode}:"
                            f"actual={global_particle_bookkeeping_mode}"
                        )
                        status = "FAIL"
                    elif (
                        require_global_particle_runtime_state_backed is True
                        and global_particle_runtime_state_backed is not True
                    ):
                        failures.append(
                            f"{case_label}:global_particle_domain_not_runtime_state_backed"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_domain_node_count is not None
                        and global_particle_domain_node_count
                        < int(min_global_particle_domain_node_count)
                    ):
                        failures.append(
                            f"{case_label}:global_particle_domain_node_count_too_small:"
                            f"minimum={int(min_global_particle_domain_node_count)}:"
                            f"actual={global_particle_domain_node_count}"
                        )
                        status = "FAIL"
                    elif not isinstance(global_particle_nodes, list) or len(global_particle_nodes) != global_particle_domain_node_count:
                        failures.append(
                            f"{case_label}:global_particle_domain_nodes_mismatch:"
                            f"expected={global_particle_domain_node_count}:actual={len(global_particle_nodes) if isinstance(global_particle_nodes, list) else 'invalid'}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_domain_edge_count is not None
                        and global_particle_domain_edge_count
                        < int(min_global_particle_domain_edge_count)
                    ):
                        failures.append(
                            f"{case_label}:global_particle_domain_edge_count_too_small:"
                            f"minimum={int(min_global_particle_domain_edge_count)}:"
                            f"actual={global_particle_domain_edge_count}"
                        )
                        status = "FAIL"
                    elif not isinstance(global_particle_edges, list) or len(global_particle_edges) != global_particle_domain_edge_count:
                        failures.append(
                            f"{case_label}:global_particle_domain_edges_mismatch:"
                            f"expected={global_particle_domain_edge_count}:actual={len(global_particle_edges) if isinstance(global_particle_edges, list) else 'invalid'}"
                        )
                        status = "FAIL"
                    elif (
                        max_global_particle_domain_charge_conservation_error_c is not None
                        and global_particle_charge_conservation_error_c
                        > float(max_global_particle_domain_charge_conservation_error_c)
                    ):
                        failures.append(
                            f"{case_label}:global_particle_domain_charge_conservation_error_exceeded:"
                            f"limit={float(max_global_particle_domain_charge_conservation_error_c)}:"
                            f"actual={global_particle_charge_conservation_error_c}"
                        )
                        status = "FAIL"
                    elif (
                        max_global_particle_domain_flux_conservation_error_a is not None
                        and global_particle_flux_conservation_error_a
                        > float(max_global_particle_domain_flux_conservation_error_a)
                    ):
                        failures.append(
                            f"{case_label}:global_particle_domain_flux_conservation_error_exceeded:"
                            f"limit={float(max_global_particle_domain_flux_conservation_error_a)}:"
                            f"actual={global_particle_flux_conservation_error_a}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_domain_edge_charge_total_abs_c is not None
                        and global_particle_edge_charge_total_abs_c
                        < float(min_global_particle_domain_edge_charge_total_abs_c)
                    ):
                        failures.append(
                            f"{case_label}:global_particle_domain_edge_charge_total_abs_too_small:"
                            f"minimum={float(min_global_particle_domain_edge_charge_total_abs_c)}:"
                            f"actual={global_particle_edge_charge_total_abs_c}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_domain_edge_target_charge_total_abs_c is not None
                        and global_particle_edge_target_charge_total_abs_c
                        < float(min_global_particle_domain_edge_target_charge_total_abs_c)
                    ):
                        failures.append(
                            f"{case_label}:global_particle_domain_edge_target_charge_total_abs_too_small:"
                            f"minimum={float(min_global_particle_domain_edge_target_charge_total_abs_c)}:"
                            f"actual={global_particle_edge_target_charge_total_abs_c}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_domain_edge_operator_drive_total_abs_c is not None
                        and global_particle_edge_operator_drive_total_abs_c
                        < float(min_global_particle_domain_edge_operator_drive_total_abs_c)
                    ):
                        failures.append(
                            f"{case_label}:global_particle_domain_edge_operator_drive_total_abs_too_small:"
                            f"minimum={float(min_global_particle_domain_edge_operator_drive_total_abs_c)}:"
                            f"actual={global_particle_edge_operator_drive_total_abs_c}"
                        )
                        status = "FAIL"
                    elif (
                        expected_transport_domain_contract
                        and transport_domain_status == "PASS"
                        and abs(
                            global_particle_edge_charge_total_abs_c
                            - transport_domain_edge_charge_total_abs_c
                        )
                        > 1.0e-15
                    ):
                        failures.append(
                            f"{case_label}:global_particle_vs_transport_edge_charge_total_mismatch:"
                            f"global={global_particle_edge_charge_total_abs_c}:"
                            f"transport={transport_domain_edge_charge_total_abs_c}"
                        )
                        status = "FAIL"
                    elif (
                        expected_transport_domain_contract
                        and transport_domain_status == "PASS"
                        and abs(
                            global_particle_edge_target_charge_total_abs_c
                            - transport_domain_edge_target_charge_total_abs_c
                        )
                        > 1.0e-15
                    ):
                        failures.append(
                            f"{case_label}:global_particle_vs_transport_edge_target_charge_total_mismatch:"
                            f"global={global_particle_edge_target_charge_total_abs_c}:"
                            f"transport={transport_domain_edge_target_charge_total_abs_c}"
                        )
                        status = "FAIL"
                    elif (
                        expected_transport_domain_contract
                        and transport_domain_status == "PASS"
                        and abs(
                            global_particle_edge_operator_drive_total_abs_c
                            - transport_domain_edge_operator_drive_total_abs_c
                        )
                        > 1.0e-15
                    ):
                        failures.append(
                            f"{case_label}:global_particle_vs_transport_edge_operator_drive_total_mismatch:"
                            f"global={global_particle_edge_operator_drive_total_abs_c}:"
                            f"transport={transport_domain_edge_operator_drive_total_abs_c}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_domain_edge_conductance_total_s is not None
                        and global_particle_edge_conductance_total_s
                        < float(min_global_particle_domain_edge_conductance_total_s)
                    ):
                        failures.append(
                            f"{case_label}:global_particle_domain_edge_conductance_total_too_small:"
                            f"minimum={float(min_global_particle_domain_edge_conductance_total_s)}:"
                            f"actual={global_particle_edge_conductance_total_s}"
                        )
                        status = "FAIL"
                    elif (
                        require_global_particle_domain_edge_graph_operator_converged
                        and not global_particle_edge_graph_operator_converged
                    ):
                        failures.append(
                            f"{case_label}:global_particle_domain_edge_graph_operator_not_converged"
                        )
                        status = "FAIL"
                    else:
                        global_particle_domain_status = "PASS"
        elif status == "PASS":
            global_particle_domain_status = "PASS"

        if status == "PASS" and expected_global_particle_repository_contract:
            actual_global_particle_repository_contract = str(
                metadata.get("surface_pic_runtime_global_particle_repository_contract_id", "")
            )
            if (
                actual_global_particle_repository_contract
                != expected_global_particle_repository_contract
            ):
                failures.append(
                    f"{case_label}:global_particle_repository_contract_mismatch:"
                    f"expected={expected_global_particle_repository_contract}:actual={actual_global_particle_repository_contract}"
                )
                status = "FAIL"
            else:
                global_particle_repository_path = global_particle_repository_sidecar(
                    output_root, metadata, csv_path
                )
                if not global_particle_repository_path.is_file():
                    failures.append(
                        f"{case_label}:missing_global_particle_repository_sidecar:{global_particle_repository_path}"
                    )
                    status = "FAIL"
                else:
                    global_particle_repository_payload = load_json(
                        global_particle_repository_path
                    )
                    global_particle_repository_active = bool(
                        global_particle_repository_payload.get(
                            "global_particle_repository_active", False
                        )
                    )
                    global_particle_repository_bookkeeping_mode = str(
                        global_particle_repository_payload.get(
                            "global_particle_repository_bookkeeping_mode", ""
                        )
                    )
                    global_particle_repository_lifecycle_mode = str(
                        global_particle_repository_payload.get(
                            "global_particle_repository_lifecycle_mode", ""
                        )
                    )
                    global_particle_repository_runtime_state_backed = (
                        global_particle_repository_payload.get(
                            "runtime_state_backed", False
                        )
                    )
                    global_particle_repository_node_count = int(
                        global_particle_repository_payload.get("shared_patch_node_count", 0)
                    )
                    global_particle_repository_edge_count = int(
                        global_particle_repository_payload.get("repository_edge_count", 0)
                    )
                    global_particle_repository_charge_conservation_error_c = float(
                        global_particle_repository_payload.get(
                            "charge_conservation_error_c", 0.0
                        )
                    )
                    global_particle_repository_migration_charge_conservation_error_c = float(
                        global_particle_repository_payload.get(
                            "migration_charge_conservation_error_c", 0.0
                        )
                    )
                    global_particle_repository_total_migration_delta_abs_charge_c = float(
                        global_particle_repository_payload.get(
                            "total_migration_delta_abs_charge_c", 0.0
                        )
                    )
                    global_particle_repository_total_edge_feedback_abs_charge_c = float(
                        global_particle_repository_payload.get(
                            "total_edge_feedback_abs_charge_c", 0.0
                        )
                    )
                    global_particle_repository_total_conservation_correction_abs_charge_c = float(
                        global_particle_repository_payload.get(
                            "total_conservation_correction_abs_charge_c", 0.0
                        )
                    )
                    global_particle_repository_total_migration_edge_abs_charge_c = float(
                        global_particle_repository_payload.get(
                            "total_migration_edge_abs_charge_c", 0.0
                        )
                    )
                    global_particle_repository_nodes = (
                        global_particle_repository_payload.get("nodes", [])
                    )
                    global_particle_repository_edges = (
                        global_particle_repository_payload.get("edges", [])
                    )

                    if (
                        str(global_particle_repository_payload.get("schema_version", ""))
                        != "scdat.surface_global_particle_repository.v1"
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_schema_mismatch:"
                            f"actual={global_particle_repository_payload.get('schema_version')}"
                        )
                        status = "FAIL"
                    elif (
                        str(global_particle_repository_payload.get("contract_id", ""))
                        != expected_global_particle_repository_contract
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_payload_contract_mismatch:"
                            f"expected={expected_global_particle_repository_contract}:actual={global_particle_repository_payload.get('contract_id')}"
                        )
                        status = "FAIL"
                    elif (
                        require_global_particle_repository_runtime_state_backed is True
                        and global_particle_repository_runtime_state_backed is not True
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_not_runtime_state_backed"
                        )
                        status = "FAIL"
                    elif (
                        expected_global_particle_repository_bookkeeping_mode
                        and global_particle_repository_bookkeeping_mode
                        != expected_global_particle_repository_bookkeeping_mode
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_bookkeeping_mode_mismatch:"
                            f"expected={expected_global_particle_repository_bookkeeping_mode}:"
                            f"actual={global_particle_repository_bookkeeping_mode}"
                        )
                        status = "FAIL"
                    elif (
                        expected_global_particle_repository_lifecycle_mode
                        and global_particle_repository_lifecycle_mode
                        != expected_global_particle_repository_lifecycle_mode
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_lifecycle_mode_mismatch:"
                            f"expected={expected_global_particle_repository_lifecycle_mode}:"
                            f"actual={global_particle_repository_lifecycle_mode}"
                        )
                        status = "FAIL"
                    elif (
                        require_global_particle_repository_active
                        and not global_particle_repository_active
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_not_active"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_repository_node_count is not None
                        and global_particle_repository_node_count
                        < int(min_global_particle_repository_node_count)
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_node_count_too_small:"
                            f"minimum={int(min_global_particle_repository_node_count)}:"
                            f"actual={global_particle_repository_node_count}"
                        )
                        status = "FAIL"
                    elif not isinstance(global_particle_repository_nodes, list) or len(
                        global_particle_repository_nodes
                    ) != global_particle_repository_node_count:
                        failures.append(
                            f"{case_label}:global_particle_repository_nodes_mismatch:"
                            f"expected={global_particle_repository_node_count}:actual={len(global_particle_repository_nodes) if isinstance(global_particle_repository_nodes, list) else 'invalid'}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_repository_edge_count is not None
                        and global_particle_repository_edge_count
                        < int(min_global_particle_repository_edge_count)
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_edge_count_too_small:"
                            f"minimum={int(min_global_particle_repository_edge_count)}:"
                            f"actual={global_particle_repository_edge_count}"
                        )
                        status = "FAIL"
                    elif not isinstance(global_particle_repository_edges, list) or len(
                        global_particle_repository_edges
                    ) != global_particle_repository_edge_count:
                        failures.append(
                            f"{case_label}:global_particle_repository_edges_mismatch:"
                            f"expected={global_particle_repository_edge_count}:actual={len(global_particle_repository_edges) if isinstance(global_particle_repository_edges, list) else 'invalid'}"
                        )
                        status = "FAIL"
                    elif (
                        max_global_particle_repository_charge_conservation_error_c is not None
                        and global_particle_repository_charge_conservation_error_c
                        > float(
                            max_global_particle_repository_charge_conservation_error_c
                        )
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_charge_conservation_error_exceeded:"
                            f"limit={float(max_global_particle_repository_charge_conservation_error_c)}:"
                            f"actual={global_particle_repository_charge_conservation_error_c}"
                        )
                        status = "FAIL"
                    elif (
                        max_global_particle_repository_migration_charge_conservation_error_c
                        is not None
                        and global_particle_repository_migration_charge_conservation_error_c
                        > float(
                            max_global_particle_repository_migration_charge_conservation_error_c
                        )
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_migration_charge_conservation_error_exceeded:"
                            f"limit={float(max_global_particle_repository_migration_charge_conservation_error_c)}:"
                            f"actual={global_particle_repository_migration_charge_conservation_error_c}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_repository_total_migration_delta_abs_charge_c
                        is not None
                        and global_particle_repository_total_migration_delta_abs_charge_c
                        < float(
                            min_global_particle_repository_total_migration_delta_abs_charge_c
                        )
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_total_migration_delta_abs_charge_too_small:"
                            f"minimum={float(min_global_particle_repository_total_migration_delta_abs_charge_c)}:"
                            f"actual={global_particle_repository_total_migration_delta_abs_charge_c}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_repository_total_edge_feedback_abs_charge_c
                        is not None
                        and global_particle_repository_total_edge_feedback_abs_charge_c
                        < float(
                            min_global_particle_repository_total_edge_feedback_abs_charge_c
                        )
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_total_edge_feedback_abs_charge_too_small:"
                            f"minimum={float(min_global_particle_repository_total_edge_feedback_abs_charge_c)}:"
                            f"actual={global_particle_repository_total_edge_feedback_abs_charge_c}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_repository_total_conservation_correction_abs_charge_c
                        is not None
                        and global_particle_repository_total_conservation_correction_abs_charge_c
                        < float(
                            min_global_particle_repository_total_conservation_correction_abs_charge_c
                        )
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_total_conservation_correction_abs_charge_too_small:"
                            f"minimum={float(min_global_particle_repository_total_conservation_correction_abs_charge_c)}:"
                            f"actual={global_particle_repository_total_conservation_correction_abs_charge_c}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_particle_repository_total_migration_edge_abs_charge_c
                        is not None
                        and global_particle_repository_total_migration_edge_abs_charge_c
                        < float(
                            min_global_particle_repository_total_migration_edge_abs_charge_c
                        )
                    ):
                        failures.append(
                            f"{case_label}:global_particle_repository_total_migration_edge_abs_charge_too_small:"
                            f"minimum={float(min_global_particle_repository_total_migration_edge_abs_charge_c)}:"
                            f"actual={global_particle_repository_total_migration_edge_abs_charge_c}"
                        )
                        status = "FAIL"
                    else:
                        global_particle_repository_status = "PASS"
        elif status == "PASS":
            global_particle_repository_status = "PASS"

        if status == "PASS" and expected_global_sheath_field_solve_contract:
            actual_global_sheath_field_solve_contract = str(
                metadata.get("surface_pic_runtime_global_sheath_field_solve_contract_id", "")
            )
            if actual_global_sheath_field_solve_contract != expected_global_sheath_field_solve_contract:
                failures.append(
                    f"{case_label}:global_sheath_field_solve_contract_mismatch:"
                    f"expected={expected_global_sheath_field_solve_contract}:actual={actual_global_sheath_field_solve_contract}"
                )
                status = "FAIL"
            else:
                global_sheath_field_path = global_sheath_field_solve_sidecar(
                    output_root, metadata, csv_path
                )
                if not global_sheath_field_path.is_file():
                    failures.append(
                        f"{case_label}:missing_global_sheath_field_solve_sidecar:{global_sheath_field_path}"
                    )
                    status = "FAIL"
                else:
                    global_sheath_field_payload = load_json(global_sheath_field_path)
                    global_sheath_field_active = bool(
                        global_sheath_field_payload.get(
                            "global_sheath_field_solve_active", False
                        )
                    )
                    global_sheath_field_mode = str(
                        global_sheath_field_payload.get(
                            "global_sheath_field_solve_mode", ""
                        )
                    )
                    global_sheath_field_runtime_state_backed = global_sheath_field_payload.get(
                        "runtime_state_backed", False
                    )
                    global_sheath_field_converged = bool(
                        global_sheath_field_payload.get(
                            "shared_global_coupled_solve_converged", False
                        )
                    )
                    global_sheath_field_node_count = int(
                        global_sheath_field_payload.get("shared_patch_node_count", 0)
                    )
                    global_sheath_field_residual_v_per_m = float(
                        global_sheath_field_payload.get("field_residual_v_per_m", 0.0)
                    )
                    global_particle_field_coupled_residual_v = float(
                        global_sheath_field_payload.get(
                            "particle_field_coupled_residual_v", 0.0
                        )
                    )
                    global_sheath_field_multi_step_stability_metric_v = float(
                        global_sheath_field_payload.get(
                            "multi_step_stability_metric_v", 0.0
                        )
                    )
                    global_sheath_field_linear_residual_norm_v = float(
                        global_sheath_field_payload.get("linear_residual_norm_v", 0.0)
                    )
                    global_sheath_field_matrix_row_count = float(
                        global_sheath_field_payload.get("matrix_row_count", 0.0)
                    )
                    global_sheath_field_matrix_nonzeros = float(
                        global_sheath_field_payload.get("matrix_nonzeros", 0.0)
                    )
                    global_sheath_field_nodes = global_sheath_field_payload.get("nodes", [])
                    if (
                        str(global_sheath_field_payload.get("schema_version", ""))
                        != "scdat.surface_global_sheath_field_solve.v1"
                    ):
                        failures.append(
                            f"{case_label}:global_sheath_field_solve_schema_mismatch:"
                            f"actual={global_sheath_field_payload.get('schema_version')}"
                        )
                        status = "FAIL"
                    elif str(global_sheath_field_payload.get("contract_id", "")) != expected_global_sheath_field_solve_contract:
                        failures.append(
                            f"{case_label}:global_sheath_field_solve_payload_contract_mismatch:"
                            f"expected={expected_global_sheath_field_solve_contract}:actual={global_sheath_field_payload.get('contract_id')}"
                        )
                        status = "FAIL"
                    elif require_global_sheath_field_solve_active and not global_sheath_field_active:
                        failures.append(
                            f"{case_label}:global_sheath_field_solve_not_active"
                        )
                        status = "FAIL"
                    elif (
                        expected_global_sheath_field_solve_mode
                        and global_sheath_field_mode
                        != expected_global_sheath_field_solve_mode
                    ):
                        failures.append(
                            f"{case_label}:global_sheath_field_solve_mode_mismatch:"
                            f"expected={expected_global_sheath_field_solve_mode}:"
                            f"actual={global_sheath_field_mode}"
                        )
                        status = "FAIL"
                    elif (
                        require_global_sheath_field_runtime_state_backed is True
                        and global_sheath_field_runtime_state_backed is not True
                    ):
                        failures.append(
                            f"{case_label}:global_sheath_field_not_runtime_state_backed"
                        )
                        status = "FAIL"
                    elif require_global_sheath_field_solve_converged and not global_sheath_field_converged:
                        failures.append(
                            f"{case_label}:global_sheath_field_solve_not_converged"
                        )
                        status = "FAIL"
                    elif (
                        min_global_sheath_field_solve_node_count is not None
                        and global_sheath_field_node_count
                        < int(min_global_sheath_field_solve_node_count)
                    ):
                        failures.append(
                            f"{case_label}:global_sheath_field_solve_node_count_too_small:"
                            f"minimum={int(min_global_sheath_field_solve_node_count)}:"
                            f"actual={global_sheath_field_node_count}"
                        )
                        status = "FAIL"
                    elif not isinstance(global_sheath_field_nodes, list) or len(global_sheath_field_nodes) != global_sheath_field_node_count:
                        failures.append(
                            f"{case_label}:global_sheath_field_solve_nodes_mismatch:"
                            f"expected={global_sheath_field_node_count}:actual={len(global_sheath_field_nodes) if isinstance(global_sheath_field_nodes, list) else 'invalid'}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_sheath_field_matrix_row_count is not None
                        and global_sheath_field_matrix_row_count
                        < float(min_global_sheath_field_matrix_row_count)
                    ):
                        failures.append(
                            f"{case_label}:global_sheath_field_matrix_row_count_too_small:"
                            f"minimum={float(min_global_sheath_field_matrix_row_count)}:"
                            f"actual={global_sheath_field_matrix_row_count}"
                        )
                        status = "FAIL"
                    elif (
                        min_global_sheath_field_matrix_nonzeros is not None
                        and global_sheath_field_matrix_nonzeros
                        < float(min_global_sheath_field_matrix_nonzeros)
                    ):
                        failures.append(
                            f"{case_label}:global_sheath_field_matrix_nonzeros_too_small:"
                            f"minimum={float(min_global_sheath_field_matrix_nonzeros)}:"
                            f"actual={global_sheath_field_matrix_nonzeros}"
                        )
                        status = "FAIL"
                    elif (
                        max_global_sheath_field_linear_residual_norm_v is not None
                        and global_sheath_field_linear_residual_norm_v
                        > float(max_global_sheath_field_linear_residual_norm_v)
                    ):
                        failures.append(
                            f"{case_label}:global_sheath_field_linear_residual_exceeded:"
                            f"limit={float(max_global_sheath_field_linear_residual_norm_v)}:"
                            f"actual={global_sheath_field_linear_residual_norm_v}"
                        )
                        status = "FAIL"
                    elif (
                        max_global_sheath_field_residual_v_per_m is not None
                        and global_sheath_field_residual_v_per_m
                        > float(max_global_sheath_field_residual_v_per_m)
                    ):
                        failures.append(
                            f"{case_label}:global_sheath_field_residual_exceeded:"
                            f"limit={float(max_global_sheath_field_residual_v_per_m)}:"
                            f"actual={global_sheath_field_residual_v_per_m}"
                        )
                        status = "FAIL"
                    elif (
                        max_global_particle_field_coupled_residual_v is not None
                        and global_particle_field_coupled_residual_v
                        > float(max_global_particle_field_coupled_residual_v)
                    ):
                        failures.append(
                            f"{case_label}:global_particle_field_coupled_residual_exceeded:"
                            f"limit={float(max_global_particle_field_coupled_residual_v)}:"
                            f"actual={global_particle_field_coupled_residual_v}"
                        )
                        status = "FAIL"
                    elif (
                        max_global_sheath_field_multi_step_stability_metric_v is not None
                        and global_sheath_field_multi_step_stability_metric_v
                        > float(max_global_sheath_field_multi_step_stability_metric_v)
                    ):
                        failures.append(
                            f"{case_label}:global_sheath_field_multi_step_stability_exceeded:"
                            f"limit={float(max_global_sheath_field_multi_step_stability_metric_v)}:"
                            f"actual={global_sheath_field_multi_step_stability_metric_v}"
                        )
                        status = "FAIL"
                    else:
                        global_sheath_field_status = "PASS"
        elif status == "PASS":
            transport_domain_status = "PASS"
            global_particle_domain_status = "PASS"
            global_particle_repository_status = "PASS"
            global_sheath_field_status = "PASS"

        if len(failures) == case_failure_start:
            status = "PASS"

        if shim_removed:
            parity_passed = status == "PASS"

        case_reports.append(
            {
                "preset": case_label,
                "command_status": command_status,
                "metadata_status": metadata_status,
                "benchmark_case_status": benchmark_case_status,
                "observer_status": observer_status,
                "consistency_status": consistency_status,
                "transport_domain_status": transport_domain_status,
                "global_particle_domain_status": global_particle_domain_status,
                "global_particle_repository_status": global_particle_repository_status,
                "global_sheath_field_status": global_sheath_field_status,
                "shim_removed": shim_removed,
                "parity_passed": parity_passed,
                "status": status,
            }
        )

    shim_removed_cases = sum(
        1 for item in case_reports if bool(item.get("shim_removed", False))
    )
    parity_applicable_cases = sum(
        1 for item in case_reports if item.get("parity_passed") is not None
    )
    parity_passed_cases = sum(1 for item in case_reports if item.get("parity_passed") is True)

    report = {
        "status": "PASS" if not failures else "FAIL",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "output_root": str(output_root),
        "matrix": str(matrix_path),
        "shim_removed_cases": shim_removed_cases,
        "parity_passed_cases": parity_passed_cases,
        "parity_applicable_cases": parity_applicable_cases,
        "cases": case_reports,
        "failures": failures,
    }

    report_json_path.parent.mkdir(parents=True, exist_ok=True)
    report_json_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_markdown(report, report_md_path)

    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
