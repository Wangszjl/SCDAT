#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import pathlib
from typing import Any


def load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_csv(path: pathlib.Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    return rows, fieldnames


def save_csv(path: pathlib.Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def as_float(text: str | None) -> float:
    if text is None:
        return float("nan")
    try:
        return float(text.strip())
    except (AttributeError, ValueError):
        return float("nan")


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


def parse_header_tokens(line: str) -> list[str]:
    return [token.strip() for token in line.lstrip("#").split(",")]


def normalize_source_alias(source_key: str) -> set[str]:
    lowered = source_key.lower()
    aliases = {lowered}
    if lowered.startswith("source"):
        suffix = lowered[len("source") :]
        aliases.add(f"0-{lowered}")
        if suffix:
            aliases.add(f"0-source{suffix}")
    return aliases


def parse_current_matrix(
    path: pathlib.Path,
    aliases: set[str],
) -> dict[str, tuple[list[float], list[float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 3:
        raise ValueError(f"Expected at least three lines in {path}")

    groups = parse_header_tokens(lines[0])
    nodes = parse_header_tokens(lines[1])
    width = min(len(groups), len(nodes))
    groups = groups[:width]
    nodes = nodes[:width]
    if width < 2:
        raise ValueError(f"Malformed matrix header in {path}")

    group_match = [
        index
        for index, token in enumerate(groups[1:], start=1)
        if token.lower() in aliases
    ]
    if not group_match:
        raise ValueError(f"Unable to find source columns for {sorted(aliases)} in {path}")

    selected: dict[str, int] = {}
    for index in group_match:
        node_name = nodes[index].lower()
        selected[node_name] = index

    time_values: list[float] = []
    series: dict[str, list[float]] = {name: [] for name in selected}
    for line in lines[2:]:
        stripped = line.strip()
        if not stripped:
            continue
        tokens = [token.strip() for token in stripped.split(",") if token.strip()]
        if len(tokens) <= max(selected.values(), default=0):
            continue
        time_values.append(float(tokens[0]))
        for node_name, index in selected.items():
            series[node_name].append(float(tokens[index]))

    return {name: (time_values, values) for name, values in series.items()}


def parse_count_series(path: pathlib.Path) -> tuple[list[float], list[float]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    time_values: list[float] = []
    counts: list[float] = []
    for line in lines[1:]:
        stripped = line.strip()
        if not stripped:
            continue
        tokens = [token.strip() for token in stripped.split(",") if token.strip()]
        if len(tokens) < 2:
            continue
        time_values.append(float(tokens[0]))
        counts.append(float(tokens[1]))
    return time_values, counts


def load_count_series_for_source(
    numkernel_root: pathlib.Path,
    source_key: str,
) -> tuple[list[float], list[float]]:
    matches = sorted(numkernel_root.glob(f"Number_of_{source_key}_*.txt"))
    if not matches:
        raise FileNotFoundError(f"Unable to find count series for source {source_key} in {numkernel_root}")

    combined_time: list[float] | None = None
    combined_values: list[float] | None = None
    for path in matches:
        time_values, count_values = parse_count_series(path)
        if combined_time is None:
            combined_time = time_values
            combined_values = list(count_values)
            continue
        if len(time_values) != len(combined_time) or any(
            abs(left - right) > 1.0e-18 for left, right in zip(time_values, combined_time)
        ):
            raise ValueError(f"Mismatched count time axes for source {source_key}: {path}")
        combined_values = [
            left + right for left, right in zip(combined_values or [], count_values)
        ]
    return combined_time or [], combined_values or []


def build_augmented_series(config: dict[str, Any]) -> dict[str, tuple[list[float], list[float]]]:
    config_block = config.get("config", {})
    spis_import = config_block.get("spis_import", {})
    numkernel_root = pathlib.Path(spis_import.get("numkernel_output_root", ""))
    source_keys = list(spis_import.get("source_keys", []))
    if not numkernel_root.exists():
        raise FileNotFoundError(f"NumKernel output root does not exist: {numkernel_root}")

    series_map: dict[str, tuple[list[float], list[float]]] = {}
    for source_key in source_keys:
        aliases = normalize_source_alias(source_key)
        collected: dict[str, tuple[list[float], list[float]]] = {}
        emitted: dict[str, tuple[list[float], list[float]]] = {}
        try:
            collected = parse_current_matrix(numkernel_root / "collectedCurrents.txt", aliases)
        except ValueError:
            collected = {}
        try:
            emitted = parse_current_matrix(numkernel_root / "emittedCurrents.txt", aliases)
        except ValueError:
            emitted = {}
        count_time, count_values = load_count_series_for_source(numkernel_root, source_key)

        if "allnodes" in collected:
            series_map[f"surface_source_{source_key}_spacecraft_collected_current_a"] = collected["allnodes"]
        if "node0" in collected:
            series_map[f"surface_node_0_source_{source_key}_collected_current_a"] = collected["node0"]
        if "node1" in collected:
            series_map[f"surface_node_1_source_{source_key}_collected_current_a"] = collected["node1"]

        if "allnodes" in emitted:
            series_map[f"surface_source_{source_key}_spacecraft_interactor_current_a"] = emitted["allnodes"]
        if "node0" in emitted:
            series_map[f"surface_node_0_source_{source_key}_interactor_current_a"] = emitted["node0"]
        if "node1" in emitted:
            series_map[f"surface_node_1_source_{source_key}_interactor_current_a"] = emitted["node1"]

        series_map[f"surface_source_{source_key}_superparticle_count"] = (count_time, count_values)
    return series_map


def expected_source_resolved_series(config: dict[str, Any]) -> list[str]:
    config_block = config.get("config", {})
    expected: list[str] = []
    for target in config_block.get("comparison_targets", []):
        hinted = (target.get("scdat_series_hint") or "").strip()
        source_key = (target.get("source_key") or "").strip()
        if source_key and hinted and hinted not in expected:
            expected.append(hinted)
    return expected


def augment(config_path: pathlib.Path, csv_path: pathlib.Path) -> dict[str, Any]:
    config = load_json(config_path)
    rows, fieldnames = load_csv(csv_path)
    if not rows:
        return {"csv_path": csv_path.as_posix(), "augmented_columns": []}

    time_values = [as_float(row.get("time_s")) for row in rows]
    series_map = build_augmented_series(config)
    augmented_columns: list[str] = []

    for column_name, (series_time, series_values) in series_map.items():
        if column_name not in fieldnames:
            fieldnames.append(column_name)
            augmented_columns.append(column_name)
        for row, query_time in zip(rows, time_values):
            value = interpolate_series(series_time, series_values, query_time)
            row[column_name] = format(value, ".17g")

    save_csv(csv_path, rows, fieldnames)

    metadata_path = csv_path.with_name(csv_path.name + ".metadata.json")
    if metadata_path.exists():
        metadata_document = load_json(metadata_path)
        metadata_block = metadata_document.setdefault("metadata", {})
        config_block = config.get("config", {})
        spis_import = config_block.get("spis_import", {})
        source_keys = list(spis_import.get("source_keys", []))
        comparison_target_source_keys = list(spis_import.get("comparison_target_source_keys", []))
        numkernel_only_source_keys = list(spis_import.get("numkernel_only_source_keys", []))
        expected_series = expected_source_resolved_series(config)
        missing_required_series = [name for name in expected_series if name not in fieldnames]
        metadata_block["surface_source_resolved_series_count"] = str(len(series_map))
        metadata_block["surface_source_resolved_augmentation_origin"] = "spis_numkernel_ascii"
        metadata_block["surface_source_resolved_bookkeeping_mode"] = "spis_numkernel_ascii_source_ledger_v1"
        metadata_block["surface_source_resolved_slot_count_complete"] = str(len(missing_required_series) == 0).lower()
        metadata_block["surface_source_resolved_strict_attribution_available"] = "true"
        metadata_block["surface_source_resolved_numkernel_output_root"] = str(
            spis_import.get("numkernel_output_root", "")
        )
        metadata_block["surface_source_resolved_source_key_count"] = str(len(source_keys))
        for index, source_key in enumerate(source_keys):
            metadata_block[f"surface_source_resolved_source_key_{index}"] = source_key
        metadata_block["surface_source_resolved_reference_source_key_count"] = str(
            len(comparison_target_source_keys)
        )
        for index, source_key in enumerate(comparison_target_source_keys):
            metadata_block[f"surface_source_resolved_reference_source_key_{index}"] = source_key
        metadata_block["surface_source_resolved_numkernel_only_source_key_count"] = str(
            len(numkernel_only_source_keys)
        )
        for index, source_key in enumerate(numkernel_only_source_keys):
            metadata_block[f"surface_source_resolved_numkernel_only_source_key_{index}"] = source_key
        metadata_block["surface_source_resolved_missing_required_series_count"] = str(len(missing_required_series))
        for index, column_name in enumerate(missing_required_series):
            metadata_block[f"surface_source_resolved_missing_required_series_{index}"] = column_name
        for index, column_name in enumerate(series_map):
            metadata_block[f"surface_source_resolved_series_{index}"] = column_name
        metadata_path.write_text(json.dumps(metadata_document, indent=2), encoding="utf-8")

    return {
        "csv_path": csv_path.as_posix(),
        "augmented_columns": augmented_columns,
        "series_count": len(series_map),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Augment SCDAT CSV outputs with source-resolved SPIS NumKernel traces."
    )
    parser.add_argument("config_json", type=pathlib.Path, help="Surface-config JSON path.")
    parser.add_argument("scdat_csv", type=pathlib.Path, help="SCDAT CSV to augment in place.")
    args = parser.parse_args()
    summary = augment(args.config_json, args.scdat_csv)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
