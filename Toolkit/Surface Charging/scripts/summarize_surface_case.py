#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import pathlib
from typing import Dict


def parse_key_value_report(path: pathlib.Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    if not path.exists():
        return values
    for line in path.read_text(encoding='utf-8').splitlines():
        if '=' not in line:
            continue
        key, value = line.split('=', 1)
        values[key.strip()] = value.strip()
    return values


def main() -> int:
    parser = argparse.ArgumentParser(
        description='Summarize a SCDAT surface charging case from CSV plus optional sidecar reports.'
    )
    parser.add_argument('csv_path')
    parser.add_argument('--benchmark-report', default='')
    parser.add_argument('--graph-report', default='')
    parser.add_argument('--monitor-json', default='')
    parser.add_argument('--field-adapter-report', default='')
    parser.add_argument('--field-adapter-json', default='')
    parser.add_argument('--graph-matrix-json', default='')
    parser.add_argument('--boundary-mapping-json', default='')
    parser.add_argument('--field-request-json', default='')
    parser.add_argument('--field-result-json', default='')
    parser.add_argument('--field-result-template-json', default='')
    parser.add_argument('--volume-stub-json', default='')
    parser.add_argument('--volume-mesh-stub-json', default='')
    parser.add_argument('--field-bridge-manifest-json', default='')
    parser.add_argument('--surface-volume-projection-json', default='')
    parser.add_argument('--volumetric-adapter-json', default='')
    parser.add_argument('--volume-request-json', default='')
    parser.add_argument('--volume-result-json', default='')
    parser.add_argument('--volume-result-template-json', default='')
    parser.add_argument('--volume-history-json', default='')
    args = parser.parse_args()

    csv_path = pathlib.Path(args.csv_path)
    benchmark_report = pathlib.Path(args.benchmark_report) if args.benchmark_report else csv_path.with_suffix('.benchmark.txt')
    graph_report = pathlib.Path(args.graph_report) if args.graph_report else csv_path.with_suffix('.graph.txt')
    monitor_json = pathlib.Path(args.monitor_json) if args.monitor_json else csv_path.with_suffix('.monitor.json')
    field_adapter_report = pathlib.Path(args.field_adapter_report) if args.field_adapter_report else csv_path.with_suffix('.field_adapter.txt')
    field_adapter_json = pathlib.Path(args.field_adapter_json) if args.field_adapter_json else csv_path.with_suffix('.field_adapter.json')
    graph_matrix_json = pathlib.Path(args.graph_matrix_json) if args.graph_matrix_json else csv_path.with_suffix('.graph_matrix.json')
    boundary_mapping_json = pathlib.Path(args.boundary_mapping_json) if args.boundary_mapping_json else csv_path.with_suffix('.boundary_mapping.json')
    field_request_json = pathlib.Path(args.field_request_json) if args.field_request_json else csv_path.with_suffix('.field_request.json')
    field_result_json = pathlib.Path(args.field_result_json) if args.field_result_json else csv_path.with_suffix('.field_result.json')
    field_result_template_json = pathlib.Path(args.field_result_template_json) if args.field_result_template_json else csv_path.with_suffix('.field_result_template.json')
    volume_stub_json = pathlib.Path(args.volume_stub_json) if args.volume_stub_json else csv_path.with_suffix('.volume_stub.json')
    volume_mesh_stub_json = pathlib.Path(args.volume_mesh_stub_json) if args.volume_mesh_stub_json else csv_path.with_suffix('.volume_mesh_stub.json')
    field_bridge_manifest_json = pathlib.Path(args.field_bridge_manifest_json) if args.field_bridge_manifest_json else csv_path.with_suffix('.field_bridge_manifest.json')
    surface_volume_projection_json = pathlib.Path(args.surface_volume_projection_json) if args.surface_volume_projection_json else csv_path.with_suffix('.surface_volume_projection.json')
    volumetric_adapter_json = pathlib.Path(args.volumetric_adapter_json) if args.volumetric_adapter_json else csv_path.with_suffix('.volumetric_adapter.json')
    volume_request_json = pathlib.Path(args.volume_request_json) if args.volume_request_json else csv_path.with_suffix('.volume_request.json')
    volume_result_json = pathlib.Path(args.volume_result_json) if args.volume_result_json else csv_path.with_suffix('.volume_result.json')
    volume_result_template_json = pathlib.Path(args.volume_result_template_json) if args.volume_result_template_json else csv_path.with_suffix('.volume_result_template.json')
    volume_history_json = pathlib.Path(args.volume_history_json) if args.volume_history_json else csv_path.with_suffix('.volume_history.json')

    with csv_path.open('r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    print(f'csv_path={csv_path}')
    print(f'row_count={len(rows)}')
    node_potential_columns = [name for name in fieldnames if name.startswith('surface_node_') and name.endswith('_potential_v')]
    interface_current_columns = [name for name in fieldnames if name.startswith('surface_interface_') and name.endswith('_current_a')]
    print(f'node_potential_column_count={len(node_potential_columns)}')
    print(f'interface_current_column_count={len(interface_current_columns)}')

    if rows:
        last = rows[-1]
        for key in [
            'surface_potential_v',
            'body_potential_v',
            'patch_potential_v',
            'floating_equilibrium_potential_v',
            'equilibrium_error_v',
        ]:
            if key in last and last[key] != '':
                print(f'last_{key}={last[key]}')
        for key in [
            'surface_node_2_pseudo_volume_m3',
            'surface_node_2_volume_projection_weight_sum',
            'surface_node_2_volume_mesh_coupling_gain',
            'surface_node_2_volume_potential_v',
            'surface_node_2_deposited_charge_c',
            'surface_node_2_poisson_residual_v_m',
        ]:
            if key in last and last[key] != '':
                print(f'last_{key}={last[key]}')

    benchmark_values = parse_key_value_report(benchmark_report)
    if benchmark_values:
        print(f'benchmark_report={benchmark_report}')
        for key in [
            'benchmark_source',
            'execution_mode',
            'consistency_status',
            'consistency_authority',
            'patch_rmse_v',
            'body_rmse_v',
            'body_je_rmse_a_per_m2',
        ]:
            if key in benchmark_values:
                print(f'{key}={benchmark_values[key]}')

    graph_values = parse_key_value_report(graph_report)
    if graph_values:
        print(f'graph_report={graph_report}')
        for key in [
            'surface_circuit_model',
            'surface_reference_graph_propagator',
            'surface_graph_capacitance_matrix_provider',
            'surface_field_solver_adapter',
            'surface_graph_capacitance_matrix_family',
            'node_count',
            'body_count',
            'patch_count',
            'branch_count',
            'body_body_interface_count',
            'body_patch_interface_count',
            'patch_patch_interface_count',
            'max_node_degree',
            'total_graph_capacitance_diagonal_f',
            'total_interface_mutual_capacitance_f',
            'max_abs_neighbor_potential_delta_v',
            'max_abs_neighbor_field_contrast_v_per_m',
            'max_abs_node_reference_offset_v',
            'max_abs_node_field_solver_reference_offset_v',
            'max_field_solver_coupling_gain',
            'strongest_field_node_name',
            'strongest_field_node_abs_field_v_per_m',
            'strongest_coupling_interface_name',
            'strongest_coupling_interface_metric',
            'graph_coupling_metric',
            'max_abs_interface_current_a',
            'max_abs_interface_power_w',
        ]:
            if key in graph_values:
                print(f'{key}={graph_values[key]}')

    if monitor_json.exists():
        print(f'monitor_json={monitor_json}')
        payload = json.loads(monitor_json.read_text(encoding='utf-8'))
        graph_payload = payload.get('graph', {})
        field_payload = payload.get('field', {})
        benchmark_payload = payload.get('benchmark', {})
        for key in [
            'graph_coupling_metric',
            'strongest_field_node_name',
            'strongest_coupling_interface_name',
        ]:
            if key in graph_payload:
                print(f'monitor_graph_{key}={graph_payload[key]}')
        for key in [
            'reference_offset_envelope_v',
            'max_field_solver_coupling_gain',
        ]:
            if key in field_payload:
                print(f'monitor_field_{key}={field_payload[key]}')
        for key in [
            'benchmark_mode',
            'benchmark_source',
            'consistency_status',
            'patch_rmse_v',
            'body_rmse_v',
        ]:
            if key in benchmark_payload:
                print(f'monitor_benchmark_{key}={benchmark_payload[key]}')

    field_adapter_values = parse_key_value_report(field_adapter_report)
    if field_adapter_values:
        print(f'field_adapter_report={field_adapter_report}')
        for key in [
            'surface_field_solver_adapter',
            'surface_graph_capacitance_matrix_provider',
            'surface_graph_capacitance_matrix_family',
            'surface_graph_capacitance_solver_adapter_hint',
            'expects_graph_capacitance_inputs',
            'supports_mutual_matrix',
            'supports_diagonal_matrix',
            'node_count',
            'branch_count',
        ]:
            if key in field_adapter_values:
                print(f'{key}={field_adapter_values[key]}')

    print(f'boundary_mapping_json_exists={int(boundary_mapping_json.exists())}')
    print(f'graph_matrix_json_exists={int(graph_matrix_json.exists())}')
    print(f'monitor_json_exists={int(monitor_json.exists())}')
    print(f'field_adapter_json_exists={int(field_adapter_json.exists())}')
    print(f'field_request_json_exists={int(field_request_json.exists())}')
    print(f'field_result_json_exists={int(field_result_json.exists())}')
    print(f'field_result_template_json_exists={int(field_result_template_json.exists())}')
    print(f'volume_stub_json_exists={int(volume_stub_json.exists())}')
    print(f'volume_mesh_stub_json_exists={int(volume_mesh_stub_json.exists())}')
    print(f'field_bridge_manifest_json_exists={int(field_bridge_manifest_json.exists())}')
    print(f'surface_volume_projection_json_exists={int(surface_volume_projection_json.exists())}')
    print(f'volumetric_adapter_json_exists={int(volumetric_adapter_json.exists())}')
    print(f'volume_request_json_exists={int(volume_request_json.exists())}')
    print(f'volume_result_json_exists={int(volume_result_json.exists())}')
    print(f'volume_result_template_json_exists={int(volume_result_template_json.exists())}')
    print(f'volume_history_json_exists={int(volume_history_json.exists())}')

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
