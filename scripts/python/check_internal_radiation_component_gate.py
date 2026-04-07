#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class InternalProfile:
    depth_m: list[float]
    electric_field_v_per_m: list[float]
    volume_charge_density_c_per_m3: list[float]


@dataclass
class ProfileMetrics:
    final_electric_field_v_per_m: float
    integrated_abs_charge_c_per_m2: float


def parse_internal_profile(csv_path: Path) -> InternalProfile:
    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "depth_m",
            "electric_field_v_per_m",
            "volume_charge_density_c_per_m3",
        }
        fields = set(reader.fieldnames or [])
        missing = sorted(required - fields)
        if missing:
            raise ValueError(
                f"CSV missing required columns in {csv_path}: {', '.join(missing)}"
            )

        depth_m: list[float] = []
        electric: list[float] = []
        rho: list[float] = []

        for line_no, row in enumerate(reader, start=2):
            try:
                depth = float(row["depth_m"])
                ef = float(row["electric_field_v_per_m"])
                density = float(row["volume_charge_density_c_per_m3"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid numeric data at line {line_no} in {csv_path}: {exc}"
                ) from exc

            if not (math.isfinite(depth) and math.isfinite(ef) and math.isfinite(density)):
                raise ValueError(f"Non-finite numeric data at line {line_no} in {csv_path}")

            depth_m.append(depth)
            electric.append(ef)
            rho.append(density)

    if not depth_m:
        raise ValueError(f"CSV has no rows: {csv_path}")

    return InternalProfile(depth_m=depth_m, electric_field_v_per_m=electric, volume_charge_density_c_per_m3=rho)


def compute_metrics(profile: InternalProfile) -> ProfileMetrics:
    integrated_abs_charge = 0.0
    for i, density in enumerate(profile.volume_charge_density_c_per_m3):
        if i == 0:
            if len(profile.depth_m) > 1:
                layer_width = profile.depth_m[1] - profile.depth_m[0]
            else:
                layer_width = profile.depth_m[0]
        else:
            layer_width = profile.depth_m[i] - profile.depth_m[i - 1]
        integrated_abs_charge += abs(density) * max(layer_width, 0.0)

    return ProfileMetrics(
        final_electric_field_v_per_m=profile.electric_field_v_per_m[-1],
        integrated_abs_charge_c_per_m2=integrated_abs_charge,
    )


def run_command(command: list[str]) -> None:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.stdout.strip():
        print(completed.stdout.strip())
    if completed.stderr.strip():
        print(completed.stderr.strip(), file=sys.stderr)
    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed: " + " ".join(command) + f" (exit code {completed.returncode})"
        )


def compare_baseline_equivalence(
    profile_a: InternalProfile,
    profile_b: InternalProfile,
    args: argparse.Namespace,
) -> tuple[float, float, float]:
    if len(profile_a.depth_m) != len(profile_b.depth_m):
        raise ValueError(
            f"Baseline row mismatch: {len(profile_a.depth_m)} vs {len(profile_b.depth_m)}"
        )

    max_depth_delta = 0.0
    max_electric_delta = 0.0
    max_rho_delta = 0.0

    for a_depth, b_depth, a_e, b_e, a_rho, b_rho in zip(
        profile_a.depth_m,
        profile_b.depth_m,
        profile_a.electric_field_v_per_m,
        profile_b.electric_field_v_per_m,
        profile_a.volume_charge_density_c_per_m3,
        profile_b.volume_charge_density_c_per_m3,
    ):
        max_depth_delta = max(max_depth_delta, abs(a_depth - b_depth))
        max_electric_delta = max(max_electric_delta, abs(a_e - b_e))
        max_rho_delta = max(max_rho_delta, abs(a_rho - b_rho))

    if max_depth_delta > args.max_baseline_depth_delta_m:
        raise ValueError(
            f"Baseline depth mismatch {max_depth_delta:.6e} exceeds {args.max_baseline_depth_delta_m:.6e}"
        )
    if max_electric_delta > args.max_baseline_electric_delta_v_per_m:
        raise ValueError(
            "Baseline electric-field mismatch "
            f"{max_electric_delta:.6e} exceeds {args.max_baseline_electric_delta_v_per_m:.6e}"
        )
    if max_rho_delta > args.max_baseline_charge_delta_c_per_m3:
        raise ValueError(
            "Baseline charge-density mismatch "
            f"{max_rho_delta:.6e} exceeds {args.max_baseline_charge_delta_c_per_m3:.6e}"
        )

    return max_depth_delta, max_electric_delta, max_rho_delta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Internal/Radiation component gate: baseline equivalence with coupling-off and trend enhancement with coupling-on."
    )
    parser.add_argument("--scdat-exe", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("build/internal_radiation_component_gate"))

    parser.add_argument("--internal-preset", default="geo_electron_belt")
    parser.add_argument("--coupling-off-radiation-preset", default="coupling_zero_flux_baseline")
    parser.add_argument("--coupling-on-radiation-preset", default="geo_electron_belt_dose")

    parser.add_argument("--max-baseline-depth-delta-m", type=float, default=1.0e-12)
    parser.add_argument("--max-baseline-electric-delta-v-per-m", type=float, default=1.0e-10)
    parser.add_argument("--max-baseline-charge-delta-c-per-m3", type=float, default=1.0e-14)

    parser.add_argument("--min-final-electric-gain-v-per-m", type=float, default=1.0e-6)
    parser.add_argument("--min-integrated-charge-gain-c-per-m2", type=float, default=1.0e-12)
    parser.add_argument("--min-coupling-on-final-electric-v-per-m", type=float, default=1.0e-6)

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        if not args.scdat_exe.is_file():
            raise FileNotFoundError(f"SCDAT executable not found: {args.scdat_exe}")

        args.output_dir.mkdir(parents=True, exist_ok=True)

        baseline_a_csv = args.output_dir / "internal_baseline_a.csv"
        baseline_b_csv = args.output_dir / "internal_baseline_b.csv"
        coupling_off_csv = args.output_dir / "internal_coupling_off.csv"
        coupling_on_csv = args.output_dir / "internal_coupling_on.csv"

        run_command([
            str(args.scdat_exe),
            "internal",
            args.internal_preset,
            str(baseline_a_csv),
        ])
        run_command([
            str(args.scdat_exe),
            "internal",
            args.internal_preset,
            str(baseline_b_csv),
        ])
        run_command([
            str(args.scdat_exe),
            "internal-radiation",
            args.internal_preset,
            args.coupling_off_radiation_preset,
            str(coupling_off_csv),
        ])
        run_command([
            str(args.scdat_exe),
            "internal-radiation",
            args.internal_preset,
            args.coupling_on_radiation_preset,
            str(coupling_on_csv),
        ])

        baseline_a = parse_internal_profile(baseline_a_csv)
        baseline_b = parse_internal_profile(baseline_b_csv)
        coupling_off = parse_internal_profile(coupling_off_csv)
        coupling_on = parse_internal_profile(coupling_on_csv)

        max_depth_delta, max_electric_delta, max_rho_delta = compare_baseline_equivalence(
            baseline_a, baseline_b, args
        )

        metrics_off = compute_metrics(coupling_off)
        metrics_on = compute_metrics(coupling_on)

        final_electric_gain = (
            metrics_on.final_electric_field_v_per_m - metrics_off.final_electric_field_v_per_m
        )
        integrated_charge_gain = (
            metrics_on.integrated_abs_charge_c_per_m2 - metrics_off.integrated_abs_charge_c_per_m2
        )

        failures: list[str] = []
        if final_electric_gain < args.min_final_electric_gain_v_per_m:
            failures.append(
                "final electric-field gain "
                f"{final_electric_gain:.6e} < min_final_electric_gain_v_per_m={args.min_final_electric_gain_v_per_m:.6e}"
            )
        if integrated_charge_gain < args.min_integrated_charge_gain_c_per_m2:
            failures.append(
                "integrated charge gain "
                f"{integrated_charge_gain:.6e} < min_integrated_charge_gain_c_per_m2={args.min_integrated_charge_gain_c_per_m2:.6e}"
            )
        if metrics_on.final_electric_field_v_per_m < args.min_coupling_on_final_electric_v_per_m:
            failures.append(
                "coupling-on final electric field "
                f"{metrics_on.final_electric_field_v_per_m:.6e} < min_coupling_on_final_electric_v_per_m={args.min_coupling_on_final_electric_v_per_m:.6e}"
            )

        print("[InternalRadiationComponentGate] Summary")
        print(f"- baseline_row_count: {len(baseline_a.depth_m)}")
        print(f"- baseline_max_depth_delta_m: {max_depth_delta:.6e}")
        print(f"- baseline_max_electric_delta_v_per_m: {max_electric_delta:.6e}")
        print(f"- baseline_max_charge_delta_c_per_m3: {max_rho_delta:.6e}")
        print(f"- coupling_off_final_electric_field_v_per_m: {metrics_off.final_electric_field_v_per_m:.6e}")
        print(f"- coupling_on_final_electric_field_v_per_m: {metrics_on.final_electric_field_v_per_m:.6e}")
        print(f"- final_electric_gain_v_per_m: {final_electric_gain:.6e}")
        print(f"- coupling_off_integrated_abs_charge_c_per_m2: {metrics_off.integrated_abs_charge_c_per_m2:.6e}")
        print(f"- coupling_on_integrated_abs_charge_c_per_m2: {metrics_on.integrated_abs_charge_c_per_m2:.6e}")
        print(f"- integrated_charge_gain_c_per_m2: {integrated_charge_gain:.6e}")

        if failures:
            print("[InternalRadiationComponentGate] FAIL")
            for item in failures:
                print(f"  - {item}")
            return 1

        print("[InternalRadiationComponentGate] PASS")
        return 0

    except Exception as exc:  # pylint: disable=broad-except
        print(f"[InternalRadiationComponentGate] ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
