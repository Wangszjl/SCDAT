#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import pathlib
import re


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize surface node/interface columns from a SCDAT CSV export.")
    parser.add_argument("csv_path")
    args = parser.parse_args()

    path = pathlib.Path(args.csv_path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])

    node_pattern = re.compile(r"surface_node_(\d+)_")
    interface_pattern = re.compile(r"surface_interface_(\d+)_")
    node_indices = set()
    interface_indices = set()
    propagated_reference_columns = 0
    mutual_capacitance_columns = 0
    field_solver_reference_columns = 0
    graph_capacitance_row_sum_columns = 0
    field_solver_coupling_gain_columns = 0
    for column in header:
        node_match = node_pattern.match(column)
        if node_match:
            node_indices.add(int(node_match.group(1)))
        interface_match = interface_pattern.match(column)
        if interface_match:
            interface_indices.add(int(interface_match.group(1)))
        if column.endswith("_propagated_reference_potential_v"):
            propagated_reference_columns += 1
        if column.endswith("_mutual_capacitance_f"):
            mutual_capacitance_columns += 1
        if column.endswith("_field_solver_reference_potential_v"):
            field_solver_reference_columns += 1
        if column.endswith("_graph_capacitance_row_sum_f"):
            graph_capacitance_row_sum_columns += 1
        if column.endswith("_field_solver_coupling_gain"):
            field_solver_coupling_gain_columns += 1

    print(f"path={path}")
    print(f"surface_node_count={len(node_indices)}")
    print(f"surface_interface_count={len(interface_indices)}")
    print(f"surface_node_propagated_reference_columns={propagated_reference_columns}")
    print(f"surface_node_field_solver_reference_columns={field_solver_reference_columns}")
    print(f"surface_node_graph_capacitance_row_sum_columns={graph_capacitance_row_sum_columns}")
    print(f"surface_node_field_solver_coupling_gain_columns={field_solver_coupling_gain_columns}")
    print(f"surface_interface_mutual_capacitance_columns={mutual_capacitance_columns}")
    print(f"node_indices={sorted(node_indices)}")
    print(f"interface_indices={sorted(interface_indices)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
