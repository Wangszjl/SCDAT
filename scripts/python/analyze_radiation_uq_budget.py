#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple


@dataclass
class ProfileMetrics:
    peak_dose_gy: float
    integrated_deposited_energy_j_per_m2: float


@dataclass
class PerturbationSpec:
    parameter: str
    rel_std: float


PERTURBATION_SPECS: List[PerturbationSpec] = [
    PerturbationSpec("particle_flux_m2_s", 0.05),
    PerturbationSpec("mean_energy_ev", 0.05),
    PerturbationSpec("thickness_m", 0.02),
    PerturbationSpec("attenuation_length_fraction", 0.10),
]

PRESET_NOMINALS: Dict[str, Dict[str, float]] = {
    "geo_electron_belt_dose": {
        "particle_flux_m2_s": 8.0e10,
        "mean_energy_ev": 8.5e4,
        "thickness_m": 4.0e-3,
        "attenuation_length_fraction": 0.20,
    },
    "meo_proton_harness_dose": {
        "particle_flux_m2_s": 1.8e10,
        "mean_energy_ev": 1.6e5,
        "thickness_m": 2.0e-3,
        "attenuation_length_fraction": 0.28,
    },
    "heavy_ion_storm_dose": {
        "particle_flux_m2_s": 2.5e9,
        "mean_energy_ev": 8.0e5,
        "thickness_m": 5.0e-3,
        "attenuation_length_fraction": 0.35,
    },
}


def load_json(path: Path) -> Dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid json at {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"expected object json at {path}")
    return payload


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def run_command(command: List[str], cwd: Path) -> None:
    completed = subprocess.run(
        command,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if completed.stdout.strip():
        print(completed.stdout.strip())
    if completed.stderr.strip():
        print(completed.stderr.strip(), file=sys.stderr)
    if completed.returncode != 0:
        raise RuntimeError(f"command failed ({completed.returncode}): {' '.join(command)}")


def deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    result = json.loads(json.dumps(base))
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def set_nested_parameter(config: Dict[str, Any], key: str, value: float) -> Dict[str, Any]:
    mutated = json.loads(json.dumps(config))
    mutated[key] = value
    return mutated


def build_nominal_config(base_preset: str, config_override: Dict[str, Any]) -> Dict[str, Any]:
    nominal: Dict[str, Any] = {}
    if base_preset in PRESET_NOMINALS:
        nominal = deep_merge(nominal, PRESET_NOMINALS[base_preset])
    nominal = deep_merge(nominal, config_override)
    return nominal


def load_profile(csv_path: Path) -> Tuple[List[float], List[float], List[float]]:
    if not csv_path.is_file():
        raise FileNotFoundError(f"missing csv: {csv_path}")

    lines = csv_path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"empty csv: {csv_path}")

    header = [item.strip() for item in lines[0].split(",")]
    expected = ["depth_m", "deposited_energy_j_per_m3", "dose_gy"]
    if header != expected:
        raise ValueError(f"unexpected header in {csv_path}: {header}")

    depth_m: List[float] = []
    deposited: List[float] = []
    dose: List[float] = []
    for row_index, raw in enumerate(lines[1:], start=2):
        parts = [item.strip() for item in raw.split(",")]
        if len(parts) != 3:
            raise ValueError(f"invalid csv row at {csv_path}:{row_index}")
        values = [float(parts[0]), float(parts[1]), float(parts[2])]
        if not all(math.isfinite(v) for v in values):
            raise ValueError(f"non-finite value at {csv_path}:{row_index}")
        depth_m.append(values[0])
        deposited.append(values[1])
        dose.append(values[2])

    if not depth_m:
        raise ValueError(f"csv has no rows: {csv_path}")

    return depth_m, deposited, dose


def depth_edges(depth_centers_m: List[float]) -> List[float]:
    if len(depth_centers_m) == 1:
        half = max(depth_centers_m[0], 1.0e-9)
        return [max(0.0, depth_centers_m[0] - half), depth_centers_m[0] + half]

    edges = [0.0] * (len(depth_centers_m) + 1)
    for idx in range(1, len(depth_centers_m)):
        edges[idx] = 0.5 * (depth_centers_m[idx - 1] + depth_centers_m[idx])
    first_step = depth_centers_m[1] - depth_centers_m[0]
    last_step = depth_centers_m[-1] - depth_centers_m[-2]
    edges[0] = max(0.0, depth_centers_m[0] - 0.5 * first_step)
    edges[-1] = depth_centers_m[-1] + 0.5 * last_step
    return edges


def metrics_from_csv(csv_path: Path) -> ProfileMetrics:
    depth_m, deposited, dose = load_profile(csv_path)
    edges = depth_edges(depth_m)
    integrated = 0.0
    for idx, dep in enumerate(deposited):
        dz = max(0.0, edges[idx + 1] - edges[idx])
        integrated += dep * dz

    return ProfileMetrics(
        peak_dose_gy=max(dose),
        integrated_deposited_energy_j_per_m2=integrated,
    )


def relative_change(candidate: float, baseline: float) -> float:
    return abs(candidate - baseline) / max(abs(baseline), 1.0e-12)


def run_radiation_case(
    project_root: Path,
    scdat_exe: Path,
    case_id: str,
    tag: str,
    base_preset: str,
    run_object: Dict[str, Any],
    config_object: Dict[str, Any],
    output_root: Path,
) -> Tuple[Path, ProfileMetrics]:
    configs_root = output_root / "generated_configs"
    configs_root.mkdir(parents=True, exist_ok=True)

    csv_path = output_root / f"{case_id}_{tag}.csv"
    config_path = configs_root / f"{case_id}_{tag}.json"

    payload = {
        "base_preset": base_preset,
        "name": f"{case_id}_{tag}",
        "run": deep_merge(run_object, {"output_csv": csv_path.as_posix()}),
        "config": config_object,
    }
    write_json(config_path, payload)

    run_command(
        [str(scdat_exe), "radiation-config", str(config_path), str(csv_path)],
        project_root,
    )

    return csv_path, metrics_from_csv(csv_path)


def format_pct(value: float) -> str:
    return f"{100.0 * value:.3f}%"


def write_markdown(report: Dict[str, Any], output_md: Path) -> None:
    lines: List[str] = []
    lines.append("# Radiation RD-010 UQ Budget Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- matrix_json: {report['matrix_json']}")
    lines.append("")

    for case in report["cases"]:
        lines.append(f"## Case: {case['case_id']}")
        lines.append("")
        base = case["base_metrics"]
        lines.append(f"- base_peak_dose_gy: {base['peak_dose_gy']:.6e}")
        lines.append(
            f"- base_integrated_deposited_energy_j_per_m2: {base['integrated_deposited_energy_j_per_m2']:.6e}"
        )
        combined = case["combined_uncertainty"]
        lines.append(f"- combined_peak_dose_rel_sigma: {format_pct(combined['peak_dose_rel_sigma'])}")
        lines.append(
            f"- combined_integral_rel_sigma: {format_pct(combined['integrated_deposited_energy_rel_sigma'])}"
        )
        lines.append("")
        lines.append(
            "| parameter | rel_std | peak_rel_response | integral_rel_response | peak_sensitivity | integral_sensitivity |"
        )
        lines.append("|---|---:|---:|---:|---:|---:|")
        for row in case["perturbations"]:
            lines.append(
                f"| {row['parameter']} | {row['rel_std']:.3f} | {row['peak_rel_response']:.6e} | {row['integral_rel_response']:.6e} | {row['peak_sensitivity']:.6e} | {row['integral_sensitivity']:.6e} |"
            )
        lines.append("")

    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="RD-010: parameter perturbation based uncertainty budget for radiation presets"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--scdat-exe", required=True)
    parser.add_argument("--matrix-json", required=True)
    parser.add_argument("--output-root", default="build/radiation_uq_budget")
    parser.add_argument("--report-json", default="")
    parser.add_argument("--report-md", default="")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    scdat_exe = Path(args.scdat_exe)
    if not scdat_exe.is_absolute():
        scdat_exe = (project_root / scdat_exe).resolve()

    matrix_json = Path(args.matrix_json)
    if not matrix_json.is_absolute():
        matrix_json = (project_root / matrix_json).resolve()

    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = (project_root / output_root).resolve()

    report_json = Path(args.report_json) if args.report_json else output_root / "radiation_uq_budget.json"
    if not report_json.is_absolute():
        report_json = (project_root / report_json).resolve()
    report_md = Path(args.report_md) if args.report_md else output_root / "radiation_uq_budget.md"
    if not report_md.is_absolute():
        report_md = (project_root / report_md).resolve()

    try:
        if not scdat_exe.is_file():
            raise FileNotFoundError(f"missing scdat executable: {scdat_exe}")
        if not matrix_json.is_file():
            raise FileNotFoundError(f"missing matrix json: {matrix_json}")

        matrix = load_json(matrix_json)
        cases = matrix.get("cases", [])
        if not isinstance(cases, list) or not cases:
            raise ValueError("matrix cases must be non-empty list")

        case_reports: List[Dict[str, Any]] = []
        for case in cases:
            if not isinstance(case, dict):
                raise ValueError("case entry must be object")
            case_id = str(case.get("case_id", "")).strip()
            base_preset = str(case.get("base_preset", "")).strip()
            if not case_id or not base_preset:
                raise ValueError(f"invalid case definition: {case}")

            deterministic_variant = case.get("deterministic", {})
            if not isinstance(deterministic_variant, dict):
                raise ValueError(f"deterministic variant must be object for {case_id}")

            run_object: Dict[str, Any] = {}
            if isinstance(case.get("run"), dict):
                run_object = case["run"]

            base_config_override: Dict[str, Any] = {}
            if isinstance(deterministic_variant.get("config"), dict):
                base_config_override = deterministic_variant["config"]

            nominal_config = build_nominal_config(base_preset, base_config_override)

            _, base_metrics = run_radiation_case(
                project_root,
                scdat_exe,
                case_id,
                "base",
                base_preset,
                run_object,
                nominal_config,
                output_root,
            )

            perturbation_rows: List[Dict[str, Any]] = []
            peak_components: List[float] = []
            integral_components: List[float] = []

            for spec in PERTURBATION_SPECS:
                if spec.parameter not in nominal_config:
                    continue
                base_value = nominal_config[spec.parameter]
                if not isinstance(base_value, (int, float)):
                    continue
                base_value_f = float(base_value)
                if not math.isfinite(base_value_f) or base_value_f <= 0.0:
                    continue

                plus_value = base_value_f * (1.0 + spec.rel_std)
                minus_value = base_value_f * max(1.0 - spec.rel_std, 1.0e-9)

                plus_config = set_nested_parameter(nominal_config, spec.parameter, plus_value)
                minus_config = set_nested_parameter(nominal_config, spec.parameter, minus_value)

                _, plus_metrics = run_radiation_case(
                    project_root,
                    scdat_exe,
                    case_id,
                    f"{spec.parameter}_plus",
                    base_preset,
                    run_object,
                    plus_config,
                    output_root,
                )
                _, minus_metrics = run_radiation_case(
                    project_root,
                    scdat_exe,
                    case_id,
                    f"{spec.parameter}_minus",
                    base_preset,
                    run_object,
                    minus_config,
                    output_root,
                )

                peak_rel_response = 0.5 * (
                    relative_change(plus_metrics.peak_dose_gy, base_metrics.peak_dose_gy)
                    + relative_change(minus_metrics.peak_dose_gy, base_metrics.peak_dose_gy)
                )
                integral_rel_response = 0.5 * (
                    relative_change(
                        plus_metrics.integrated_deposited_energy_j_per_m2,
                        base_metrics.integrated_deposited_energy_j_per_m2,
                    )
                    + relative_change(
                        minus_metrics.integrated_deposited_energy_j_per_m2,
                        base_metrics.integrated_deposited_energy_j_per_m2,
                    )
                )

                peak_components.append(peak_rel_response)
                integral_components.append(integral_rel_response)

                perturbation_rows.append(
                    {
                        "parameter": spec.parameter,
                        "base_value": base_value_f,
                        "rel_std": spec.rel_std,
                        "plus_value": plus_value,
                        "minus_value": minus_value,
                        "peak_rel_response": peak_rel_response,
                        "integral_rel_response": integral_rel_response,
                        "peak_sensitivity": peak_rel_response / spec.rel_std,
                        "integral_sensitivity": integral_rel_response / spec.rel_std,
                        "plus_metrics": {
                            "peak_dose_gy": plus_metrics.peak_dose_gy,
                            "integrated_deposited_energy_j_per_m2": plus_metrics.integrated_deposited_energy_j_per_m2,
                        },
                        "minus_metrics": {
                            "peak_dose_gy": minus_metrics.peak_dose_gy,
                            "integrated_deposited_energy_j_per_m2": minus_metrics.integrated_deposited_energy_j_per_m2,
                        },
                    }
                )

            peak_combined = math.sqrt(sum(value * value for value in peak_components))
            integral_combined = math.sqrt(sum(value * value for value in integral_components))

            case_reports.append(
                {
                    "case_id": case_id,
                    "base_preset": base_preset,
                    "base_metrics": {
                        "peak_dose_gy": base_metrics.peak_dose_gy,
                        "integrated_deposited_energy_j_per_m2": base_metrics.integrated_deposited_energy_j_per_m2,
                    },
                    "perturbations": perturbation_rows,
                    "combined_uncertainty": {
                        "peak_dose_rel_sigma": peak_combined,
                        "integrated_deposited_energy_rel_sigma": integral_combined,
                    },
                }
            )

        status = "PASS"
        report = {
            "status": status,
            "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "project_root": project_root.as_posix(),
            "matrix_json": matrix_json.as_posix(),
            "output_root": output_root.as_posix(),
            "perturbation_specs": [
                {"parameter": spec.parameter, "rel_std": spec.rel_std}
                for spec in PERTURBATION_SPECS
            ],
            "cases": case_reports,
        }

        write_json(report_json, report)
        write_markdown(report, report_md)

        print(f"report_json={report_json.as_posix()}")
        print(f"report_md={report_md.as_posix()}")
        print(f"status={status}")
        return 0

    except Exception as exc:  # pylint: disable=broad-except
        print(f"status=FAIL({exc})")
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
