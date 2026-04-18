#!/usr/bin/env python3
from __future__ import annotations

import json
import pathlib
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Any


@dataclass
class SpisCaseArtifacts:
    case_root: pathlib.Path
    output_root: pathlib.Path
    config_pwl_path: pathlib.Path
    discrete_manifest_path: pathlib.Path
    import_manifest_path: pathlib.Path
    scan_plan_path: pathlib.Path
    spis_reference_output_root: pathlib.Path
    discrete_config_paths: list[pathlib.Path]


def _text_value(parameter: ET.Element) -> Any:
    string_value = (parameter.findtext("valueAsString") or "").strip()
    if string_value:
        return string_value

    bool_value = (parameter.findtext("valueAsBoolean") or "").strip().lower()
    if bool_value in {"true", "false"}:
        return bool_value == "true"

    for tag in ("valueAsDouble", "valueAsFloat", "valueAsInt", "valueAsLong"):
        text = (parameter.findtext(tag) or "").strip()
        if not text:
            continue
        try:
            numeric = float(text)
        except ValueError:
            continue
        if abs(numeric) > 0.0 or text in {"0", "0.0", "0.0f"}:
            return numeric
    return ""


def parse_global_parameters(path: pathlib.Path) -> dict[str, Any]:
    root = ET.parse(path).getroot()
    parameters: dict[str, Any] = {}
    for parameter in root.findall(".//GlobalParameter"):
        key = (parameter.findtext("keyName") or "").strip()
        if not key:
            continue
        parameters[key] = _text_value(parameter)
    return parameters


def parse_model_paths(case_root: pathlib.Path) -> dict[str, pathlib.Path]:
    model_path = case_root / "model.xml"
    root = ET.parse(model_path).getroot()

    def find_first(tag_name: str) -> pathlib.Path | None:
        element = root.find(f".//{tag_name}")
        if element is None or not (element.text or "").strip():
            return None
        return case_root / element.text.strip()

    global_parameters = find_first("GlobalParameters")
    groups = find_first("Groups")
    circuit = find_first("ElectricalCircuit")
    run_root = global_parameters.parent.parent if global_parameters else case_root / "DefaultStudy/Simulations/Run1"
    output_root = run_root / "OutputFolder"
    return {
        "model_xml_path": model_path,
        "global_parameters_path": global_parameters or run_root / "GlobalParameters/globalParameters.xml",
        "groups_path": groups or case_root / "DefaultStudy/Preprocessing/Groups/groups.xml",
        "circuit_path": circuit or case_root / "DefaultStudy/Preprocessing/ElectricalCircuit/circuit.txt",
        "run_root": run_root,
        "output_root": output_root,
    }


def parse_groups(groups_path: pathlib.Path) -> dict[str, Any]:
    root = ET.parse(groups_path).getroot()
    electrical_nodes: list[dict[str, Any]] = []
    materials: list[dict[str, Any]] = []

    for prop in root.findall(".//Property"):
        characteristics: dict[str, Any] = {}
        for entry in prop.findall("./characteristicList/entry"):
            key = (entry.findtext("string") or "").strip()
            if not key:
                continue
            value_element = next((child for child in entry if child.tag != "string"), None)
            if value_element is None:
                continue
            value_text = (value_element.findtext("value") or "").strip()
            if not value_text:
                continue
            try:
                characteristics[key] = float(value_text)
            except ValueError:
                characteristics[key] = value_text

        if "ElecNodeId" in characteristics:
            electrical_nodes.append(
                {
                    "id": int(characteristics["ElecNodeId"]),
                    "name": (prop.findtext("name") or "").strip(),
                    "type": (prop.findtext("type") or "").strip(),
                    "characteristics": characteristics,
                }
            )

        if "MatTypeId" in characteristics or "MatModelId" in characteristics:
            material_entry = {
                "name": (prop.findtext("name") or "").strip(),
                "characteristics": characteristics,
                "material_scalars": {},
            }
            for child in prop.findall("./subPropertiesList/entry/Property/characteristicList/entry"):
                key = (child.findtext("string") or "").strip()
                value_element = next((grand for grand in child if grand.tag != "string"), None)
                value_text = (value_element.findtext("value") or "").strip() if value_element is not None else ""
                if not key or not value_text:
                    continue
                try:
                    material_entry["material_scalars"][key] = float(value_text)
                except ValueError:
                    material_entry["material_scalars"][key] = value_text
            materials.append(material_entry)

    electrical_nodes.sort(key=lambda item: item["id"])
    return {"electrical_nodes": electrical_nodes, "materials": materials}


def parse_pwl(circuit_path: pathlib.Path) -> dict[str, Any]:
    text = circuit_path.read_text(encoding="utf-8")
    match = re.search(r"V\s+(\d+)\s+(\d+)\s+PWL\(([^)]*)\)", text, re.IGNORECASE | re.MULTILINE)
    if not match:
        raise ValueError(f"Unable to find PWL voltage source in {circuit_path}")

    from_node = int(match.group(1))
    to_node = int(match.group(2))
    payload = re.sub(r"\b[A-Za-z_]+\s*=\s*[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", "", match.group(3))
    numbers = [float(token) for token in re.findall(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", payload)]
    if len(numbers) < 2 or len(numbers) % 2 != 0:
        raise ValueError(f"Malformed PWL payload in {circuit_path}")

    points = [
        {"time_s": numbers[index], "value": numbers[index + 1], "label": f"point_{index // 2:03d}"}
        for index in range(0, len(numbers), 2)
    ]
    return {"from_node": from_node, "to_node": to_node, "points": points, "source_text": text.strip()}


def _material_from_groups(groups: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    unsupported: list[str] = []
    material = {
        "id": 85,
        "type": "dielectric",
        "name": "child_langmuir_ito",
        "permittivity": 3.0,
        "conductivity": 1.0e-12,
        "work_function_ev": 4.7,
        "secondary_electron_yield": 0.49,
        "scalar_properties": {},
    }
    if not groups["materials"]:
        unsupported.append("groups.materials")
        return material, unsupported

    scalars = groups["materials"][0].get("material_scalars", {})
    if "RDC" in scalars:
        material["permittivity"] = float(scalars["RDC"])
    else:
        unsupported.append("material.RDC")
    if "BULKDCONDUCT" in scalars:
        material["conductivity"] = max(0.0, float(scalars["BULKDCONDUCT"]))
    elif "BUC" in scalars:
        conductivity = float(scalars["BUC"])
        material["conductivity"] = 1.0e12 if conductivity < 0.0 else max(0.0, conductivity)
    else:
        unsupported.append("material.BULKDCONDUCT")
    if "WFE" in scalars:
        material["work_function_ev"] = float(scalars["WFE"])
    else:
        unsupported.append("material.WFE")
    if "SEY" in scalars:
        material["secondary_electron_yield"] = max(0.0, float(scalars["SEY"]))
    else:
        unsupported.append("material.SEY")
    if "PEY" in scalars:
        material["scalar_properties"]["photoelectron_yield"] = max(0.0, float(scalars["PEY"]))
    return material, unsupported


def _natural_sort_key(value: str) -> list[Any]:
    return [int(token) if token.isdigit() else token.lower() for token in re.split(r"(\d+)", value)]


def _append_unique(items: list[str], value: str) -> None:
    candidate = value.strip()
    if candidate and candidate not in items:
        items.append(candidate)


def _is_aggregate_source_key(value: str) -> bool:
    normalized = value.strip().lower()
    return normalized.startswith("all_") or normalized in {"allnodes", "all_nodes", "all"}


def discover_source_keys(reference_output_root: pathlib.Path, numkernel_output_root: pathlib.Path) -> list[str]:
    discovered: list[str] = []
    monitored_root = reference_output_root / "DataFieldMonitored"
    extracted_root = reference_output_root / "DataFieldExtracted"
    patterns = [
        (monitored_root, "Individual_current_on_spacecraft_-_Collected_*.nc", r"Collected_(.+?)(?:_\(|--)"),
        (monitored_root, "Individual_current_on_spacecraft_-_Interactor_*.nc", r"Interactor_(.+?)(?:_\(|--)"),
        (monitored_root, "Number_of_*.nc", r"Number_of_(.+?)(?:_\(|--)"),
        (extracted_root, "Collected_current_A_versus_time_s_*.nc", r"Collected_current_A_versus_time_s_(.+?),_node_"),
        (extracted_root, "Emitted_current_A_versus_time_s_*.nc", r"Emitted_current_A_versus_time_s_(.+?),_node_"),
    ]
    for root, glob_pattern, regex in patterns:
        if not root.exists():
            continue
        for candidate in sorted(root.glob(glob_pattern)):
            if "-VERTEXMask-" in candidate.name or "--Compared" in candidate.name or "--Reference" in candidate.name:
                continue
            match = re.search(regex, candidate.name)
            if match and not _is_aggregate_source_key(match.group(1)):
                _append_unique(discovered, match.group(1))

    if numkernel_output_root.exists():
        for candidate in sorted(numkernel_output_root.glob("Number_of_*_*.txt")):
            match = re.match(r"Number_of_(.+?)_.+\.txt$", candidate.name)
            if match and not _is_aggregate_source_key(match.group(1)):
                _append_unique(discovered, match.group(1))

    if not discovered:
        discovered.append("source1")
    return sorted(discovered, key=_natural_sort_key)


def _has_reference_file(reference_output_root: pathlib.Path, family: str, pattern: str) -> bool:
    family_root = reference_output_root / family
    if not family_root.exists():
        return False
    for candidate in sorted(family_root.glob(pattern)):
        if "-VERTEXMask-" in candidate.name or "--Compared" in candidate.name or "--Reference" in candidate.name:
            continue
        return True
    return False


def build_source_resolved_comparison_targets(
    source_keys: list[str],
    node_indices: list[int],
    reference_output_root: pathlib.Path,
) -> list[dict[str, Any]]:
    targets: list[dict[str, Any]] = []
    for source_key in source_keys:
        spacecraft_collected_pattern = f"Individual_current_on_spacecraft_-_Collected_{source_key}*.nc"
        if _has_reference_file(reference_output_root, "DataFieldMonitored", spacecraft_collected_pattern):
            targets.append(
                {
                    "name": f"Individual_current_on_spacecraft_Collected_{source_key}",
                    "source_family": "DataFieldMonitored",
                    "source_pattern": spacecraft_collected_pattern,
                    "unit": "A",
                    "source_key": source_key,
                    "scdat_series_hint": f"surface_source_{source_key}_spacecraft_collected_current_a",
                }
            )
        spacecraft_interactor_pattern = f"Individual_current_on_spacecraft_-_Interactor_{source_key}*.nc"
        if _has_reference_file(reference_output_root, "DataFieldMonitored", spacecraft_interactor_pattern):
            targets.append(
                {
                    "name": f"Individual_current_on_spacecraft_Interactor_{source_key}",
                    "source_family": "DataFieldMonitored",
                    "source_pattern": spacecraft_interactor_pattern,
                    "unit": "A",
                    "source_key": source_key,
                    "scdat_series_hint": f"surface_source_{source_key}_spacecraft_interactor_current_a",
                }
            )
        for node_index in node_indices:
            collected_pattern = f"Collected_current_A_versus_time_s_{source_key},_node_{node_index}.nc"
            if _has_reference_file(reference_output_root, "DataFieldExtracted", collected_pattern):
                targets.append(
                    {
                        "name": f"Collected_current_{source_key}_node_{node_index}",
                        "source_family": "DataFieldExtracted",
                        "source_pattern": collected_pattern,
                        "unit": "A",
                        "source_key": source_key,
                        "node_index": node_index,
                        "scdat_series_hint": f"surface_node_{node_index}_source_{source_key}_collected_current_a",
                    }
                )
            interactor_pattern = f"Emitted_current_A_versus_time_s_{source_key},_node_{node_index}.nc"
            if _has_reference_file(reference_output_root, "DataFieldExtracted", interactor_pattern):
                targets.append(
                    {
                        "name": f"Interactor_current_{source_key}_node_{node_index}",
                        "source_family": "DataFieldExtracted",
                        "source_pattern": interactor_pattern,
                        "unit": "A",
                        "source_key": source_key,
                        "node_index": node_index,
                        "scdat_series_hint": f"surface_node_{node_index}_source_{source_key}_interactor_current_a",
                    }
                )
        count_pattern = f"Number_of_{source_key}*.nc"
        if _has_reference_file(reference_output_root, "DataFieldMonitored", count_pattern):
            targets.append(
                {
                    "name": f"Number_of_{source_key}",
                    "source_family": "DataFieldMonitored",
                    "source_pattern": count_pattern,
                    "unit": "count",
                    "source_key": source_key,
                    "scdat_series_hint": f"surface_source_{source_key}_superparticle_count",
                }
            )
    return targets


def build_child_langmuir_documents(case_root: pathlib.Path, output_root: pathlib.Path) -> SpisCaseArtifacts:
    case_root = case_root.resolve()
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    generated_root = output_root / "generated"
    generated_root.mkdir(parents=True, exist_ok=True)
    discrete_root = generated_root / "discrete_points"
    discrete_root.mkdir(parents=True, exist_ok=True)

    model_paths = parse_model_paths(case_root)
    global_parameters = parse_global_parameters(model_paths["global_parameters_path"])
    groups = parse_groups(model_paths["groups_path"])
    pwl = parse_pwl(model_paths["circuit_path"])
    material, unsupported_fields = _material_from_groups(groups)

    sample_period_s = float(global_parameters.get("instrumentSamplePeriod", 4.0e-8) or 4.0e-8)
    duration_s = float(global_parameters.get("duration", pwl["points"][-1]["time_s"]) or pwl["points"][-1]["time_s"])
    step_count = max(1, int(round(duration_s / max(sample_period_s, 1.0e-12))))

    electrical_nodes = groups["electrical_nodes"]
    node0 = electrical_nodes[0] if electrical_nodes else {"id": 0, "name": "node_0"}
    node1 = electrical_nodes[1] if len(electrical_nodes) > 1 else {"id": 1, "name": "node_1"}

    reference_output_root = model_paths["output_root"]
    numkernel_output_root = model_paths["run_root"] / "NumKernel/Output"
    source_keys = discover_source_keys(reference_output_root, numkernel_output_root)
    node_indices = list(range(max(2, len(electrical_nodes) if electrical_nodes else 0)))
    comparison_targets = [
        {
            "name": "Average_surface_potential_of_node_0",
            "source_family": "DataFieldMonitored",
            "source_pattern": "Average_surface_potential_of_node_0*.nc",
            "unit": "V",
        },
        {
            "name": "Average_surface_potential_of_node_1",
            "source_family": "DataFieldMonitored",
            "source_pattern": "Average_surface_potential_of_node_1*.nc",
            "unit": "V",
        },
        {
            "name": "Total_current_on_node_0",
            "source_family": "DataFieldMonitored",
            "source_pattern": "Total_current_on_node_0*.nc",
            "unit": "A",
        },
        {
            "name": "Total_current_on_node_1",
            "source_family": "DataFieldMonitored",
            "source_pattern": "Total_current_on_node_1*.nc",
            "unit": "A",
        },
        {
            "name": "Total_current_on_spacecraft",
            "source_family": "DataFieldMonitored",
            "source_pattern": "Total_current_on_spacecraft*.nc",
            "unit": "A",
        },
    ]
    comparison_targets.extend(
        build_source_resolved_comparison_targets(source_keys, node_indices, reference_output_root)
    )
    comparison_target_source_keys = sorted(
        {
            str(item.get("source_key", "")).strip()
            for item in comparison_targets
            if str(item.get("source_key", "")).strip()
        },
        key=_natural_sort_key,
    )
    numkernel_only_source_keys = [
        source_key for source_key in source_keys if source_key not in comparison_target_source_keys
    ]

    import_manifest = {
        "case_name": "Child_Langmuir",
        "case_id": "child_langmuir",
        "case_root": case_root.as_posix(),
        "source_files": {key: value.as_posix() for key, value in model_paths.items()},
        "source_keys": source_keys,
        "comparison_target_source_keys": comparison_target_source_keys,
        "numkernel_only_source_keys": numkernel_only_source_keys,
        "numkernel_output_root": numkernel_output_root.as_posix(),
        "global_parameters_subset": {
            key: global_parameters.get(key)
            for key in (
                "environmentType",
                "electronDistrib2",
                "ionDistrib2",
                "electronTemperature",
                "electronTemperature2",
                "ionTemperature",
                "ionTemperature2",
                "duration",
                "instrumentSamplePeriod",
            )
        },
        "electrical_nodes": electrical_nodes,
        "unsupported_fields": unsupported_fields,
        "mapping_notes": [
            "SPIS electrical nodes are mapped onto surface_nodes/body-patch circuit nodes.",
            "PWL source is mapped onto circuit_waveform targeting branch_bias.",
            "Discrete sweep configs are emitted as per-point constant-bias cases for cross-validation.",
        ],
    }

    scan_plan = {
        "mode": "both",
        "timeline_points": pwl["points"],
        "discrete_points": [
            {"index": index, "time_s": point["time_s"], "bias_v": point["value"], "label": point["label"]}
            for index, point in enumerate(pwl["points"])
        ],
    }

    common_config = {
        "schema_version": "v1",
        "module": "surface",
        "name": "child_langmuir_spis_import",
        "base_preset": "geo_ecss_kapton_pic_circuit",
        "run": {
            "time_step_s": sample_period_s,
            "steps": step_count,
            "adaptive_time_stepping": False,
        },
        "config": {
            "runtime_route": "unified",
            "surface_pic_runtime": "graph_coupled_shared_surface",
            "surface_instrument_set": "surface_pic_observer_set",
            "current_algorithm_mode": "unified",
            "benchmark_mode": "unified",
            "enable_live_pic_window": True,
            "enable_live_pic_mcc": False,
            "enable_pic_calibration": False,
            "enable_body_patch_circuit": True,
            "body_floating": False,
            "body_initial_potential_v": 0.0,
            "body_capacitance_f": 1.0e-10,
            "surface_area_m2": 1.0e-4,
            "internal_substeps": 1,
            "surface_material": material,
            "surface_nodes": [
                {
                    "name": f"node_{node0['id']}_ground",
                    "area_m2": 1.0e-4,
                    "is_patch": False,
                    "initial_potential_v": 0.0,
                    "capacitance_f": 1.0e-10,
                    "fixed_potential": True,
                    "fixed_value_v": 0.0,
                },
                {
                    "name": f"node_{node1['id']}_probe",
                    "area_m2": 1.0e-4,
                    "is_patch": True,
                    "initial_potential_v": 0.0,
                    "capacitance_f": 1.0e-10,
                    "fixed_potential": False,
                    "fixed_value_v": 0.0,
                },
            ],
            "surface_branches": [
                {
                    "from_node": 1,
                    "to_node": 0,
                    "conductance_s": 1.0e-9,
                    "bias_v": 0.0,
                }
            ],
            "spis_import": {
                "case_name": "Child_Langmuir",
                "case_id": "child_langmuir",
                "case_root": case_root.as_posix(),
                "model_xml_path": model_paths["model_xml_path"].as_posix(),
                "study_root": (case_root / "DefaultStudy").as_posix(),
                "run_root": model_paths["run_root"].as_posix(),
                "global_parameters_path": model_paths["global_parameters_path"].as_posix(),
                "groups_path": model_paths["groups_path"].as_posix(),
                "circuit_path": model_paths["circuit_path"].as_posix(),
                "import_manifest_path": (generated_root / "child_langmuir.import_manifest.json").as_posix(),
                "scan_plan_path": (generated_root / "child_langmuir.scan_plan.json").as_posix(),
                "numkernel_output_root": numkernel_output_root.as_posix(),
                "source_keys": source_keys,
                "comparison_target_source_keys": comparison_target_source_keys,
                "numkernel_only_source_keys": numkernel_only_source_keys,
            },
            "comparison_targets": comparison_targets,
            "spis_reference_output": {
                "output_root": reference_output_root.as_posix(),
                "monitored_root": (reference_output_root / "DataFieldMonitored").as_posix(),
                "extracted_root": (reference_output_root / "DataFieldExtracted").as_posix(),
                "comparison_csv_path": (output_root / "child_langmuir.comparison.csv").as_posix(),
                "comparison_summary_path": (output_root / "child_langmuir.comparison.md").as_posix(),
                "comparison_json_path": (output_root / "child_langmuir.comparison.json").as_posix(),
            },
        },
    }

    pwl_config = json.loads(json.dumps(common_config))
    pwl_config["output_csv"] = (output_root / "child_langmuir_pwl.csv").as_posix()
    pwl_config["run"]["output_csv"] = pwl_config["output_csv"]
    pwl_config["config"]["potential_sweep"] = {
        "mode": "pwl_timeline",
        "active_point_index": 0,
        "points": pwl["points"],
    }
    pwl_config["config"]["circuit_waveform"] = [
        {
            "name": "child_langmuir_probe_bias",
            "target_kind": "branch_bias",
            "target_index": 0,
            "waveform": "pwl",
            "pwl_points": [{"time_s": point["time_s"], "value": point["value"]} for point in pwl["points"]],
        }
    ]

    discrete_paths: list[pathlib.Path] = []
    discrete_entries: list[dict[str, Any]] = []
    for entry in scan_plan["discrete_points"]:
        config_document = json.loads(json.dumps(common_config))
        config_document["name"] = f"child_langmuir_bias_{entry['index']:03d}"
        csv_path = output_root / "discrete_points" / f"child_langmuir_bias_{entry['index']:03d}.csv"
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        config_document["output_csv"] = csv_path.as_posix()
        config_document["run"]["output_csv"] = csv_path.as_posix()
        config_document["config"]["potential_sweep"] = {
            "mode": "discrete_points",
            "active_point_index": entry["index"],
            "points": [
                {
                    "time_s": entry["time_s"],
                    "value": entry["bias_v"],
                    "label": entry["label"],
                }
            ],
        }
        config_document["config"]["surface_branches"][0]["bias_v"] = entry["bias_v"]
        config_document["config"]["circuit_waveform"] = []
        path = discrete_root / f"child_langmuir_bias_{entry['index']:03d}.surface.json"
        path.write_text(json.dumps(config_document, indent=2), encoding="utf-8")
        discrete_paths.append(path)
        discrete_entries.append(
            {
                "index": entry["index"],
                "bias_v": entry["bias_v"],
                "time_s": entry["time_s"],
                "config_path": path.as_posix(),
                "output_csv": csv_path.as_posix(),
            }
        )

    import_manifest_path = generated_root / "child_langmuir.import_manifest.json"
    scan_plan_path = generated_root / "child_langmuir.scan_plan.json"
    config_pwl_path = generated_root / "child_langmuir.pwl.surface.json"
    discrete_manifest_path = generated_root / "child_langmuir.discrete_points.json"

    import_manifest_path.write_text(json.dumps(import_manifest, indent=2), encoding="utf-8")
    scan_plan_path.write_text(json.dumps(scan_plan, indent=2), encoding="utf-8")
    config_pwl_path.write_text(json.dumps(pwl_config, indent=2), encoding="utf-8")
    discrete_manifest_path.write_text(json.dumps({"cases": discrete_entries}, indent=2), encoding="utf-8")

    return SpisCaseArtifacts(
        case_root=case_root,
        output_root=output_root,
        config_pwl_path=config_pwl_path,
        discrete_manifest_path=discrete_manifest_path,
        import_manifest_path=import_manifest_path,
        scan_plan_path=scan_plan_path,
        spis_reference_output_root=reference_output_root,
        discrete_config_paths=discrete_paths,
    )


def build_plasma_wake_documents(case_root: pathlib.Path, output_root: pathlib.Path) -> SpisCaseArtifacts:
    case_root = case_root.resolve()
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    generated_root = output_root / "generated"
    generated_root.mkdir(parents=True, exist_ok=True)

    model_paths = parse_model_paths(case_root)
    global_parameters = parse_global_parameters(model_paths["global_parameters_path"])
    groups = parse_groups(model_paths["groups_path"])
    material, unsupported_fields = _material_from_groups(groups)

    sample_period_s = float(global_parameters.get("instrumentSamplePeriod") or global_parameters.get("plasmaDt") or 2.8e-6)
    duration_s = float(global_parameters.get("duration") or 4.0e-4)
    step_count = max(1, int(round(duration_s / max(sample_period_s, 1.0e-12))))

    electrical_nodes = groups["electrical_nodes"]
    node0 = electrical_nodes[0] if electrical_nodes else {"id": 0, "name": "node_0"}

    reference_output_root = model_paths["output_root"]
    numkernel_output_root = model_paths["run_root"] / "NumKernel/Output"
    source_keys = discover_source_keys(reference_output_root, numkernel_output_root)
    node_indices = [0]
    comparison_targets = [
        {
            "name": "Average_surface_potential_of_node_0",
            "source_family": "DataFieldMonitored",
            "source_pattern": "Average_surface_potential_of_node_0*.nc",
            "unit": "V",
        },
        {
            "name": "Total_current_on_node_0",
            "source_family": "DataFieldMonitored",
            "source_pattern": "Total_current_on_node_0*.nc",
            "unit": "A",
        },
        {
            "name": "Total_current_on_spacecraft",
            "source_family": "DataFieldMonitored",
            "source_pattern": "Total_current_on_spacecraft*.nc",
            "unit": "A",
        },
    ]
    comparison_targets.extend(
        build_source_resolved_comparison_targets(source_keys, node_indices, reference_output_root)
    )
    comparison_target_source_keys = sorted(
        {
            str(item.get("source_key", "")).strip()
            for item in comparison_targets
            if str(item.get("source_key", "")).strip()
        },
        key=_natural_sort_key,
    )
    numkernel_only_source_keys = [
        source_key for source_key in source_keys if source_key not in comparison_target_source_keys
    ]

    import_manifest = {
        "case_name": "Plasma_Wake",
        "case_id": "plasma_wake",
        "case_root": case_root.as_posix(),
        "source_files": {key: value.as_posix() for key, value in model_paths.items()},
        "source_keys": source_keys,
        "comparison_target_source_keys": comparison_target_source_keys,
        "numkernel_only_source_keys": numkernel_only_source_keys,
        "numkernel_output_root": numkernel_output_root.as_posix(),
        "global_parameters_subset": {
            key: global_parameters.get(key)
            for key in (
                "environmentType",
                "electronDistrib",
                "ionDistrib",
                "electronTemperature",
                "ionTemperature",
                "duration",
                "plasmaDt",
                "instrumentSamplePeriod",
            )
        },
        "electrical_nodes": electrical_nodes,
        "unsupported_fields": unsupported_fields,
        "mapping_notes": [
            "SPIS single-node spacecraft is mapped onto a floating surface body node.",
            "Static Plasma_Wake timeline uses a single run config without PWL/discrete sweep artifacts.",
            "Source-resolved compare targets are emitted only when a matching SPIS Java netCDF reference exists.",
        ],
    }

    scan_plan = {"mode": "static_timeline", "timeline_points": [{"time_s": 0.0, "value": 0.0, "label": "static"}], "discrete_points": []}

    config_document = {
        "schema_version": "v1",
        "module": "surface",
        "name": "plasma_wake_spis_import",
        "base_preset": "geo_ecss_kapton_pic_circuit",
        "run": {
            "time_step_s": sample_period_s,
            "steps": step_count,
            "adaptive_time_stepping": False,
        },
        "output_csv": (output_root / "plasma_wake.csv").as_posix(),
        "config": {
            "runtime_route": "unified",
            "surface_pic_runtime": "graph_coupled_shared_surface",
            "surface_instrument_set": "surface_pic_observer_set",
            "current_algorithm_mode": "unified",
            "benchmark_mode": "unified",
            "enable_live_pic_window": True,
            "enable_live_pic_mcc": False,
            "enable_pic_calibration": False,
            "enable_body_patch_circuit": False,
            "body_floating": True,
            "body_initial_potential_v": 0.0,
            "body_capacitance_f": 1.0e-10,
            "surface_area_m2": 1.0,
            "internal_substeps": 1,
            "surface_material": material,
            "surface_nodes": [
                {
                    "name": f"node_{node0['id']}_spacecraft",
                    "area_m2": 1.0,
                    "is_patch": False,
                    "initial_potential_v": 0.0,
                    "capacitance_f": 1.0e-10,
                    "fixed_potential": False,
                    "fixed_value_v": 0.0,
                }
            ],
            "surface_branches": [],
            "spis_import": {
                "case_name": "Plasma_Wake",
                "case_id": "plasma_wake",
                "case_root": case_root.as_posix(),
                "model_xml_path": model_paths["model_xml_path"].as_posix(),
                "study_root": (case_root / "DefaultStudy").as_posix(),
                "run_root": model_paths["run_root"].as_posix(),
                "global_parameters_path": model_paths["global_parameters_path"].as_posix(),
                "groups_path": model_paths["groups_path"].as_posix(),
                "circuit_path": model_paths["circuit_path"].as_posix(),
                "import_manifest_path": (generated_root / "plasma_wake.import_manifest.json").as_posix(),
                "scan_plan_path": (generated_root / "plasma_wake.scan_plan.json").as_posix(),
                "numkernel_output_root": numkernel_output_root.as_posix(),
                "source_keys": source_keys,
                "comparison_target_source_keys": comparison_target_source_keys,
                "numkernel_only_source_keys": numkernel_only_source_keys,
            },
            "comparison_targets": comparison_targets,
            "spis_reference_output": {
                "output_root": reference_output_root.as_posix(),
                "monitored_root": (reference_output_root / "DataFieldMonitored").as_posix(),
                "extracted_root": (reference_output_root / "DataFieldExtracted").as_posix(),
                "comparison_csv_path": (output_root / "plasma_wake.comparison.csv").as_posix(),
                "comparison_summary_path": (output_root / "plasma_wake.comparison.md").as_posix(),
                "comparison_json_path": (output_root / "plasma_wake.comparison.json").as_posix(),
            },
        },
    }
    config_document["run"]["output_csv"] = config_document["output_csv"]

    import_manifest_path = generated_root / "plasma_wake.import_manifest.json"
    scan_plan_path = generated_root / "plasma_wake.scan_plan.json"
    config_path = generated_root / "plasma_wake.surface.json"
    discrete_manifest_path = generated_root / "plasma_wake.discrete_points.json"

    import_manifest_path.write_text(json.dumps(import_manifest, indent=2), encoding="utf-8")
    scan_plan_path.write_text(json.dumps(scan_plan, indent=2), encoding="utf-8")
    config_path.write_text(json.dumps(config_document, indent=2), encoding="utf-8")
    discrete_manifest_path.write_text(json.dumps({"cases": []}, indent=2), encoding="utf-8")

    return SpisCaseArtifacts(
        case_root=case_root,
        output_root=output_root,
        config_pwl_path=config_path,
        discrete_manifest_path=discrete_manifest_path,
        import_manifest_path=import_manifest_path,
        scan_plan_path=scan_plan_path,
        spis_reference_output_root=reference_output_root,
        discrete_config_paths=[],
    )


def build_spis_case_documents(case_root: pathlib.Path, output_root: pathlib.Path) -> SpisCaseArtifacts:
    case_name = case_root.name.lower()
    if "child_langmuir" in case_name:
        return build_child_langmuir_documents(case_root, output_root)
    if "plasma_wake" in case_name:
        return build_plasma_wake_documents(case_root, output_root)
    raise ValueError(f"Unsupported SPIS case for importer: {case_root}")
