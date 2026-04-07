#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import statistics
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class CouplingProfile:
    depth_m: list[float]
    electric_field_v_per_m: list[float]
    volume_charge_density_c_per_m3: list[float]


@dataclass
class GateSummary:
    row_count: int
    min_depth_step_m: float
    min_electric_delta_v_per_m: float
    max_electric_field_v_per_m: float
    final_electric_field_v_per_m: float
    integrated_abs_charge_c_per_m2: float
    ratio_sample_count: int
    ratio_mean: float
    ratio_cv: float
    rho_increase_count: int


def parse_profile(csv_path: Path) -> CouplingProfile:
    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        expected_columns = {
            "depth_m",
            "electric_field_v_per_m",
            "volume_charge_density_c_per_m3",
        }
        field_names = set(reader.fieldnames or [])
        missing = sorted(expected_columns - field_names)
        if missing:
            raise ValueError(
                "CSV missing required columns: " + ", ".join(missing)
            )

        depth_m: list[float] = []
        electric_field_v_per_m: list[float] = []
        volume_charge_density_c_per_m3: list[float] = []

        for line_no, row in enumerate(reader, start=2):
            try:
                depth = float(row["depth_m"])
                electric = float(row["electric_field_v_per_m"])
                rho = float(row["volume_charge_density_c_per_m3"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid numeric data at line {line_no}: {exc}"
                ) from exc

            if not (math.isfinite(depth) and math.isfinite(electric) and math.isfinite(rho)):
                raise ValueError(f"Non-finite numeric data at line {line_no}")

            depth_m.append(depth)
            electric_field_v_per_m.append(electric)
            volume_charge_density_c_per_m3.append(rho)

    if not depth_m:
        raise ValueError(f"CSV has no data rows: {csv_path}")

    return CouplingProfile(
        depth_m=depth_m,
        electric_field_v_per_m=electric_field_v_per_m,
        volume_charge_density_c_per_m3=volume_charge_density_c_per_m3,
    )


def maybe_run_coupled_case(args: argparse.Namespace) -> None:
    if not args.scdat_exe:
        return

    exe_path = Path(args.scdat_exe)
    if not exe_path.is_file():
        raise FileNotFoundError(f"SCDAT executable not found: {exe_path}")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(exe_path),
        "internal-radiation",
        args.internal_preset,
        args.radiation_preset,
        str(args.output_csv),
    ]

    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.stdout.strip():
        print(completed.stdout.strip())
    if completed.stderr.strip():
        print(completed.stderr.strip(), file=sys.stderr)
    if completed.returncode != 0:
        raise RuntimeError(
            "Failed to run coupled command: "
            + " ".join(command)
            + f" (exit code {completed.returncode})"
        )


def evaluate_gate(profile: CouplingProfile, args: argparse.Namespace) -> tuple[GateSummary, list[str]]:
    failures: list[str] = []
    row_count = len(profile.depth_m)

    if row_count < args.min_rows:
        failures.append(
            f"row_count={row_count} is smaller than min_rows={args.min_rows}"
        )

    depth_steps = [
        profile.depth_m[i] - profile.depth_m[i - 1]
        for i in range(1, row_count)
    ]
    min_depth_step = min(depth_steps) if depth_steps else profile.depth_m[0]
    if min_depth_step <= args.min_depth_step_m:
        failures.append(
            f"minimum depth step {min_depth_step:.6e} <= min_depth_step_m={args.min_depth_step_m:.6e}"
        )

    integrated_abs_charge = 0.0
    ratio_samples: list[float] = []
    rho_increase_count = 0
    previous_electric = 0.0
    min_electric_delta = math.inf

    for i in range(row_count):
        electric = profile.electric_field_v_per_m[i]
        rho = profile.volume_charge_density_c_per_m3[i]
        delta_electric = electric - previous_electric
        previous_electric = electric

        if i > 0 and (rho - profile.volume_charge_density_c_per_m3[i - 1]) > args.rho_increase_tolerance:
            rho_increase_count += 1

        if abs(rho) >= args.rho_floor:
            ratio_samples.append(delta_electric / rho)

        if i == 0:
            if row_count > 1:
                layer_width = profile.depth_m[1] - profile.depth_m[0]
            else:
                layer_width = profile.depth_m[0]
        else:
            layer_width = profile.depth_m[i] - profile.depth_m[i - 1]

        integrated_abs_charge += abs(rho) * max(layer_width, 0.0)
        min_electric_delta = min(min_electric_delta, delta_electric)

    ratio_sample_count = len(ratio_samples)
    ratio_mean = statistics.fmean(ratio_samples) if ratio_samples else 0.0
    ratio_cv = (
        statistics.pstdev(ratio_samples) / abs(ratio_mean)
        if ratio_sample_count > 1 and abs(ratio_mean) > 0.0
        else 0.0
    )

    if ratio_sample_count < args.min_ratio_samples:
        failures.append(
            f"ratio sample count {ratio_sample_count} is smaller than min_ratio_samples={args.min_ratio_samples}"
        )

    if abs(ratio_mean) < args.min_ratio_mean:
        failures.append(
            f"ratio mean {ratio_mean:.6e} is smaller than min_ratio_mean={args.min_ratio_mean:.6e}"
        )

    if ratio_cv > args.max_ratio_cv:
        failures.append(
            f"ratio CV {ratio_cv:.6e} exceeds max_ratio_cv={args.max_ratio_cv:.6e}"
        )

    if integrated_abs_charge < args.min_integrated_abs_charge_c_per_m2:
        failures.append(
            "integrated absolute charge "
            f"{integrated_abs_charge:.6e} is smaller than "
            f"min_integrated_abs_charge_c_per_m2={args.min_integrated_abs_charge_c_per_m2:.6e}"
        )

    if profile.electric_field_v_per_m[-1] < args.min_final_electric_field_v_per_m:
        failures.append(
            "final electric field "
            f"{profile.electric_field_v_per_m[-1]:.6e} is smaller than "
            f"min_final_electric_field_v_per_m={args.min_final_electric_field_v_per_m:.6e}"
        )

    if rho_increase_count > args.max_rho_increase_count:
        failures.append(
            f"rho increase count {rho_increase_count} exceeds max_rho_increase_count={args.max_rho_increase_count}"
        )

    if args.require_monotonic_electric_field and min_electric_delta < -args.max_electric_drop_v_per_m:
        failures.append(
            "minimum electric-field delta "
            f"{min_electric_delta:.6e} < -max_electric_drop_v_per_m={args.max_electric_drop_v_per_m:.6e}"
        )

    summary = GateSummary(
        row_count=row_count,
        min_depth_step_m=min_depth_step,
        min_electric_delta_v_per_m=min_electric_delta,
        max_electric_field_v_per_m=max(profile.electric_field_v_per_m),
        final_electric_field_v_per_m=profile.electric_field_v_per_m[-1],
        integrated_abs_charge_c_per_m2=integrated_abs_charge,
        ratio_sample_count=ratio_sample_count,
        ratio_mean=ratio_mean,
        ratio_cv=ratio_cv,
        rho_increase_count=rho_increase_count,
    )
    return summary, failures


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate one-way internal-radiation coupling output with trend and consistency checks."
    )
    parser.add_argument(
        "--csv",
        dest="output_csv",
        type=Path,
        default=Path("results/internal_radiation_geo_electron_belt__geo_electron_belt_dose.csv"),
        help="CSV output to validate.",
    )
    parser.add_argument(
        "--scdat-exe",
        type=Path,
        default=None,
        help="Optional SCDAT executable path. If provided, the script first runs internal-radiation.",
    )
    parser.add_argument(
        "--internal-preset",
        default="geo_electron_belt",
        help="Internal preset used when --scdat-exe is provided.",
    )
    parser.add_argument(
        "--radiation-preset",
        default="geo_electron_belt_dose",
        help="Radiation preset used when --scdat-exe is provided.",
    )

    parser.add_argument("--min-rows", type=int, default=16)
    parser.add_argument("--min-depth-step-m", type=float, default=0.0)

    parser.add_argument("--rho-floor", type=float, default=1.0e-14)
    parser.add_argument("--min-ratio-samples", type=int, default=8)
    parser.add_argument("--min-ratio-mean", type=float, default=1.0e-6)
    parser.add_argument("--max-ratio-cv", type=float, default=1.0e-4)

    parser.add_argument("--rho-increase-tolerance", type=float, default=1.0e-14)
    parser.add_argument("--max-rho-increase-count", type=int, default=0)

    parser.add_argument(
        "--require-monotonic-electric-field",
        action="store_true",
        help="Require electric field profile to be monotonic non-decreasing.",
    )
    parser.add_argument("--max-electric-drop-v-per-m", type=float, default=1.0e-9)

    parser.add_argument("--min-integrated-abs-charge-c-per-m2", type=float, default=1.0e-12)
    parser.add_argument("--min-final-electric-field-v-per-m", type=float, default=1.0e-6)

    return parser


def main() -> int:
    parser = build_argument_parser()
    args = parser.parse_args()

    try:
        maybe_run_coupled_case(args)
        profile = parse_profile(args.output_csv)
        summary, failures = evaluate_gate(profile, args)
    except Exception as exc:  # pylint: disable=broad-except
        print(f"[InternalRadiationCouplingGate] ERROR: {exc}", file=sys.stderr)
        return 1

    print("[InternalRadiationCouplingGate] Summary")
    print(f"- csv: {args.output_csv}")
    print(f"- row_count: {summary.row_count}")
    print(f"- min_depth_step_m: {summary.min_depth_step_m:.6e}")
    print(f"- min_electric_delta_v_per_m: {summary.min_electric_delta_v_per_m:.6e}")
    print(f"- max_electric_field_v_per_m: {summary.max_electric_field_v_per_m:.6e}")
    print(f"- final_electric_field_v_per_m: {summary.final_electric_field_v_per_m:.6e}")
    print(f"- integrated_abs_charge_c_per_m2: {summary.integrated_abs_charge_c_per_m2:.6e}")
    print(f"- ratio_sample_count: {summary.ratio_sample_count}")
    print(f"- ratio_mean: {summary.ratio_mean:.6e}")
    print(f"- ratio_cv: {summary.ratio_cv:.6e}")
    print(f"- rho_increase_count: {summary.rho_increase_count}")

    if failures:
        print("[InternalRadiationCouplingGate] FAIL")
        for item in failures:
            print(f"  - {item}")
        return 1

    print("[InternalRadiationCouplingGate] PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
