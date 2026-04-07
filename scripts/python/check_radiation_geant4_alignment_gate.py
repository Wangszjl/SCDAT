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
from typing import Any, Dict, List


@dataclass
class RadiationProfile:
    depth_m: List[float]
    deposited_energy_j_per_m3: List[float]
    dose_gy: List[float]


@dataclass
class ProfileMetrics:
    rows: int
    peak_dose_gy: float
    surface_dose_gy: float
    back_dose_gy: float
    integrated_deposited_energy_j_per_m2: float


@dataclass
class CaseSummary:
    case_id: str
    deterministic_csv: str
    monte_carlo_csv: str
    deterministic: ProfileMetrics
    monte_carlo: ProfileMetrics
    profile_rmse: float
    integral_ratio_mc_to_det: float
    peak_ratio_mc_to_det: float


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
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}"
        )


def load_profile(csv_path: Path) -> RadiationProfile:
    if not csv_path.is_file():
        raise FileNotFoundError(f"missing csv: {csv_path}")

    lines = csv_path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"empty csv: {csv_path}")

    header = [item.strip() for item in lines[0].split(",")]
    expected = ["depth_m", "deposited_energy_j_per_m3", "dose_gy"]
    if header != expected:
        raise ValueError(
            f"unexpected header in {csv_path}: {header}, expected {expected}"
        )

    depth_m: List[float] = []
    deposited: List[float] = []
    dose: List[float] = []
    for row_index, raw in enumerate(lines[1:], start=2):
        parts = [item.strip() for item in raw.split(",")]
        if len(parts) != 3:
            raise ValueError(f"invalid csv row at {csv_path}:{row_index}")
        values = [float(parts[0]), float(parts[1]), float(parts[2])]
        if not all(math.isfinite(v) for v in values):
            raise ValueError(f"non-finite csv row at {csv_path}:{row_index}")
        depth_m.append(values[0])
        deposited.append(values[1])
        dose.append(values[2])

    if not depth_m:
        raise ValueError(f"csv has no data rows: {csv_path}")

    return RadiationProfile(depth_m, deposited, dose)


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


def integrate_deposited_energy(profile: RadiationProfile) -> float:
    edges = depth_edges(profile.depth_m)
    total = 0.0
    for idx, deposited in enumerate(profile.deposited_energy_j_per_m3):
        dz = max(0.0, edges[idx + 1] - edges[idx])
        total += deposited * dz
    return total


def metrics_from_profile(profile: RadiationProfile) -> ProfileMetrics:
    return ProfileMetrics(
        rows=len(profile.depth_m),
        peak_dose_gy=max(profile.dose_gy),
        surface_dose_gy=profile.dose_gy[0],
        back_dose_gy=profile.dose_gy[-1],
        integrated_deposited_energy_j_per_m2=integrate_deposited_energy(profile),
    )


def normalized_profile(values: List[float]) -> List[float]:
    non_negative = [max(0.0, value) for value in values]
    norm = max(sum(non_negative), 1.0e-30)
    return [value / norm for value in non_negative]


def profile_rmse(candidate: List[float], reference: List[float]) -> float:
    if len(candidate) != len(reference):
        raise ValueError("profile length mismatch")
    if not candidate:
        return 0.0
    candidate_norm = normalized_profile(candidate)
    reference_norm = normalized_profile(reference)
    mse = 0.0
    for cand, ref in zip(candidate_norm, reference_norm):
        diff = cand - ref
        mse += diff * diff
    mse /= float(len(candidate))
    return math.sqrt(mse)


def as_dict(metrics: ProfileMetrics) -> Dict[str, float]:
    return {
        "rows": metrics.rows,
        "peak_dose_gy": metrics.peak_dose_gy,
        "surface_dose_gy": metrics.surface_dose_gy,
        "back_dose_gy": metrics.back_dose_gy,
        "integrated_deposited_energy_j_per_m2": metrics.integrated_deposited_energy_j_per_m2,
    }


def deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    result = json.loads(json.dumps(base))
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def case_config_payload(case: Dict[str, Any], variant_override: Dict[str, Any], output_csv: Path) -> Dict[str, Any]:
    payload: Dict[str, Any] = {
        "base_preset": case["base_preset"],
        "name": f"{case['case_id']}_{variant_override['variant_id']}",
        "run": {
            "output_csv": output_csv.as_posix(),
        },
        "config": {},
    }
    if isinstance(case.get("run"), dict):
        payload["run"] = deep_merge(payload["run"], case["run"])

    override_config = variant_override.get("config", {})
    if not isinstance(override_config, dict):
        raise ValueError("variant config must be object")
    payload["config"] = deep_merge(payload["config"], override_config)
    return payload


def relative_delta(candidate: float, baseline: float) -> float:
    return abs(candidate - baseline) / max(abs(baseline), 1.0e-12)


def compare_with_baseline(
    case_summary: CaseSummary,
    baseline_case: Dict[str, Any],
    drift_thresholds: Dict[str, float],
    case_thresholds: Dict[str, float],
) -> List[str]:
    failures: List[str] = []

    det_baseline = baseline_case.get("deterministic", {})
    mc_baseline = baseline_case.get("monte_carlo", {})
    cross_baseline = baseline_case.get("cross_metrics", {})

    det_peak_rel = relative_delta(case_summary.deterministic.peak_dose_gy, float(det_baseline["peak_dose_gy"]))
    if det_peak_rel > drift_thresholds["det_peak_dose_rel"]:
        failures.append(
            f"det_peak_dose_rel={det_peak_rel:.6e}>{drift_thresholds['det_peak_dose_rel']:.6e}"
        )

    det_integral_rel = relative_delta(
        case_summary.deterministic.integrated_deposited_energy_j_per_m2,
        float(det_baseline["integrated_deposited_energy_j_per_m2"]),
    )
    if det_integral_rel > drift_thresholds["det_integral_rel"]:
        failures.append(
            f"det_integral_rel={det_integral_rel:.6e}>{drift_thresholds['det_integral_rel']:.6e}"
        )

    mc_peak_rel = relative_delta(case_summary.monte_carlo.peak_dose_gy, float(mc_baseline["peak_dose_gy"]))
    if mc_peak_rel > drift_thresholds["mc_peak_dose_rel"]:
        failures.append(
            f"mc_peak_dose_rel={mc_peak_rel:.6e}>{drift_thresholds['mc_peak_dose_rel']:.6e}"
        )

    mc_integral_rel = relative_delta(
        case_summary.monte_carlo.integrated_deposited_energy_j_per_m2,
        float(mc_baseline["integrated_deposited_energy_j_per_m2"]),
    )
    if mc_integral_rel > drift_thresholds["mc_integral_rel"]:
        failures.append(
            f"mc_integral_rel={mc_integral_rel:.6e}>{drift_thresholds['mc_integral_rel']:.6e}"
        )

    rmse_rel = relative_delta(case_summary.profile_rmse, float(cross_baseline["profile_rmse"]))
    if rmse_rel > drift_thresholds["profile_rmse_rel"]:
        failures.append(
            f"profile_rmse_rel={rmse_rel:.6e}>{drift_thresholds['profile_rmse_rel']:.6e}"
        )

    if case_summary.profile_rmse > float(case_thresholds["max_profile_rmse"]):
        failures.append(
            f"profile_rmse={case_summary.profile_rmse:.6e}>{float(case_thresholds['max_profile_rmse']):.6e}"
        )

    if case_summary.integral_ratio_mc_to_det < float(case_thresholds["integral_ratio_min"]):
        failures.append(
            f"integral_ratio_mc_to_det={case_summary.integral_ratio_mc_to_det:.6e}<{float(case_thresholds['integral_ratio_min']):.6e}"
        )

    if case_summary.integral_ratio_mc_to_det > float(case_thresholds["integral_ratio_max"]):
        failures.append(
            f"integral_ratio_mc_to_det={case_summary.integral_ratio_mc_to_det:.6e}>{float(case_thresholds['integral_ratio_max']):.6e}"
        )

    return failures


def write_markdown(report: Dict[str, Any], output_md: Path) -> None:
    lines: List[str] = []
    lines.append("# Radiation Geant4 Alignment Gate")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- matrix_json: {report['matrix_json']}")
    lines.append(f"- baseline_json: {report['baseline_json']}")
    lines.append("")
    lines.append("| case_id | profile_rmse | integral_ratio_mc_to_det | peak_ratio_mc_to_det | status |")
    lines.append("|---|---:|---:|---:|---|")
    for case in report["cases"]:
        lines.append(
            f"| {case['case_id']} | {case['profile_rmse']:.6e} | {case['integral_ratio_mc_to_det']:.6e} | {case['peak_ratio_mc_to_det']:.6e} | {case['status']} |"
        )

    if report["failures"]:
        lines.append("")
        lines.append("## Failures")
        lines.append("")
        for failure in report["failures"]:
            lines.append(f"- {failure}")

    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="RD-009: fixed-geometry/material Geant4 alignment benchmark gate"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--scdat-exe", required=True)
    parser.add_argument("--matrix-json", required=True)
    parser.add_argument("--baseline-json", required=True)
    parser.add_argument("--output-root", default="build/radiation_geant4_alignment_gate")
    parser.add_argument("--summary-json", default="")
    parser.add_argument("--summary-md", default="")
    parser.add_argument("--write-baseline", action="store_true")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    scdat_exe = Path(args.scdat_exe)
    if not scdat_exe.is_absolute():
        scdat_exe = (project_root / scdat_exe).resolve()

    matrix_json = Path(args.matrix_json)
    if not matrix_json.is_absolute():
        matrix_json = (project_root / matrix_json).resolve()

    baseline_json = Path(args.baseline_json)
    if not baseline_json.is_absolute():
        baseline_json = (project_root / baseline_json).resolve()

    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = (project_root / output_root).resolve()

    summary_json = Path(args.summary_json) if args.summary_json else output_root / "radiation_geant4_alignment_gate.json"
    if not summary_json.is_absolute():
        summary_json = (project_root / summary_json).resolve()

    summary_md = Path(args.summary_md) if args.summary_md else output_root / "radiation_geant4_alignment_gate.md"
    if not summary_md.is_absolute():
        summary_md = (project_root / summary_md).resolve()

    try:
        if not scdat_exe.is_file():
            raise FileNotFoundError(f"missing scdat executable: {scdat_exe}")
        matrix = load_json(matrix_json)
        cases = matrix.get("cases", [])
        if not isinstance(cases, list) or not cases:
            raise ValueError("matrix cases must be a non-empty array")

        drift_thresholds = matrix.get("drift_thresholds", {})
        required_drift = [
            "det_peak_dose_rel",
            "det_integral_rel",
            "mc_peak_dose_rel",
            "mc_integral_rel",
            "profile_rmse_rel",
        ]
        for key in required_drift:
            if key not in drift_thresholds:
                raise ValueError(f"missing drift_thresholds.{key}")

        output_root.mkdir(parents=True, exist_ok=True)
        generated_config_root = output_root / "generated_configs"
        generated_config_root.mkdir(parents=True, exist_ok=True)

        case_summaries: List[CaseSummary] = []

        for case in cases:
            if not isinstance(case, dict):
                raise ValueError("case must be object")
            case_id = str(case.get("case_id", "")).strip()
            base_preset = str(case.get("base_preset", "")).strip()
            if not case_id or not base_preset:
                raise ValueError(f"invalid case entry: {case}")

            case_thresholds = case.get("thresholds", {})
            for key in ["max_profile_rmse", "integral_ratio_min", "integral_ratio_max"]:
                if key not in case_thresholds:
                    raise ValueError(f"missing thresholds.{key} in case {case_id}")

            det_csv = output_root / f"{case_id}_deterministic.csv"
            mc_csv = output_root / f"{case_id}_monte_carlo.csv"
            det_config = generated_config_root / f"{case_id}_deterministic.json"
            mc_config = generated_config_root / f"{case_id}_monte_carlo.json"

            det_variant = case.get("deterministic", {})
            if not isinstance(det_variant, dict):
                raise ValueError(f"deterministic variant must be object in case {case_id}")
            det_variant = deep_merge(det_variant, {"variant_id": "deterministic"})
            mc_variant = case.get("monte_carlo", {})
            if not isinstance(mc_variant, dict):
                raise ValueError(f"monte_carlo variant must be object in case {case_id}")
            mc_variant = deep_merge(mc_variant, {"variant_id": "monte_carlo"})

            det_payload = case_config_payload(
                deep_merge(case, {"base_preset": base_preset}),
                det_variant,
                det_csv,
            )
            mc_payload = case_config_payload(
                deep_merge(case, {"base_preset": base_preset}),
                mc_variant,
                mc_csv,
            )

            write_json(det_config, det_payload)
            write_json(mc_config, mc_payload)

            run_command([str(scdat_exe), "radiation-config", str(det_config), str(det_csv)], project_root)
            run_command([str(scdat_exe), "radiation-config", str(mc_config), str(mc_csv)], project_root)

            det_profile = load_profile(det_csv)
            mc_profile = load_profile(mc_csv)
            if len(det_profile.depth_m) != len(mc_profile.depth_m):
                raise ValueError(f"profile length mismatch in case {case_id}")

            det_metrics = metrics_from_profile(det_profile)
            mc_metrics = metrics_from_profile(mc_profile)
            rmse = profile_rmse(mc_profile.dose_gy, det_profile.dose_gy)
            integral_ratio = mc_metrics.integrated_deposited_energy_j_per_m2 / max(
                det_metrics.integrated_deposited_energy_j_per_m2, 1.0e-12
            )
            peak_ratio = mc_metrics.peak_dose_gy / max(det_metrics.peak_dose_gy, 1.0e-12)

            case_summaries.append(
                CaseSummary(
                    case_id=case_id,
                    deterministic_csv=det_csv.as_posix(),
                    monte_carlo_csv=mc_csv.as_posix(),
                    deterministic=det_metrics,
                    monte_carlo=mc_metrics,
                    profile_rmse=rmse,
                    integral_ratio_mc_to_det=integral_ratio,
                    peak_ratio_mc_to_det=peak_ratio,
                )
            )

        baseline_payload = {
            "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "matrix_json": matrix_json.as_posix(),
            "cases": [
                {
                    "case_id": case.case_id,
                    "deterministic": as_dict(case.deterministic),
                    "monte_carlo": as_dict(case.monte_carlo),
                    "cross_metrics": {
                        "profile_rmse": case.profile_rmse,
                        "integral_ratio_mc_to_det": case.integral_ratio_mc_to_det,
                        "peak_ratio_mc_to_det": case.peak_ratio_mc_to_det,
                    },
                }
                for case in case_summaries
            ],
        }

        failures: List[str] = []
        cases_report: List[Dict[str, Any]] = []

        if args.write_baseline:
            baseline_json.parent.mkdir(parents=True, exist_ok=True)
            baseline_json.write_text(
                json.dumps(baseline_payload, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
        else:
            if not baseline_json.is_file():
                raise FileNotFoundError(
                    f"missing baseline json: {baseline_json} (run with --write-baseline to create)"
                )
            baseline = load_json(baseline_json)
            baseline_cases = baseline.get("cases", [])
            if not isinstance(baseline_cases, list):
                raise ValueError("baseline cases must be array")
            baseline_by_id = {
                str(case.get("case_id", "")).strip(): case
                for case in baseline_cases
                if isinstance(case, dict)
            }

            matrix_case_by_id = {
                str(case.get("case_id", "")).strip(): case
                for case in cases
                if isinstance(case, dict)
            }

            for case in case_summaries:
                case_failures: List[str] = []
                baseline_case = baseline_by_id.get(case.case_id)
                if baseline_case is None:
                    case_failures.append(f"missing baseline case: {case.case_id}")
                else:
                    thresholds = matrix_case_by_id[case.case_id].get("thresholds", {})
                    case_failures.extend(
                        compare_with_baseline(
                            case,
                            baseline_case,
                            drift_thresholds,
                            thresholds,
                        )
                    )

                status = "PASS" if not case_failures else "FAIL"
                if case_failures:
                    failures.extend([f"{case.case_id}:{item}" for item in case_failures])

                cases_report.append(
                    {
                        "case_id": case.case_id,
                        "status": status,
                        "deterministic_csv": case.deterministic_csv,
                        "monte_carlo_csv": case.monte_carlo_csv,
                        "deterministic": as_dict(case.deterministic),
                        "monte_carlo": as_dict(case.monte_carlo),
                        "profile_rmse": case.profile_rmse,
                        "integral_ratio_mc_to_det": case.integral_ratio_mc_to_det,
                        "peak_ratio_mc_to_det": case.peak_ratio_mc_to_det,
                        "failures": case_failures,
                    }
                )

        if args.write_baseline:
            status = "PASS"
            cases_report = [
                {
                    "case_id": case.case_id,
                    "status": "BASELINE_WRITTEN",
                    "deterministic_csv": case.deterministic_csv,
                    "monte_carlo_csv": case.monte_carlo_csv,
                    "deterministic": as_dict(case.deterministic),
                    "monte_carlo": as_dict(case.monte_carlo),
                    "profile_rmse": case.profile_rmse,
                    "integral_ratio_mc_to_det": case.integral_ratio_mc_to_det,
                    "peak_ratio_mc_to_det": case.peak_ratio_mc_to_det,
                    "failures": [],
                }
                for case in case_summaries
            ]
            failures = []
        else:
            status = "PASS" if not failures else "FAIL"

        report = {
            "status": status,
            "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "project_root": project_root.as_posix(),
            "matrix_json": matrix_json.as_posix(),
            "baseline_json": baseline_json.as_posix(),
            "output_root": output_root.as_posix(),
            "write_baseline": bool(args.write_baseline),
            "cases": cases_report,
            "failures": failures,
        }

        summary_json.parent.mkdir(parents=True, exist_ok=True)
        summary_json.write_text(
            json.dumps(report, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        write_markdown(report, summary_md)

        print(f"summary_json={summary_json.as_posix()}")
        print(f"summary_md={summary_md.as_posix()}")
        print(f"baseline_json={baseline_json.as_posix()}")
        print(f"status={status}")

        if args.write_baseline:
            return 0
        return 0 if status == "PASS" else 2

    except Exception as exc:  # pylint: disable=broad-except
        print(f"status=FAIL({exc})")
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
