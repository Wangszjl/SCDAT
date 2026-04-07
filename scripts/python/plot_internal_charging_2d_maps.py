#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _load_profile(csv_path: Path):
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    if not rows:
        raise ValueError(f"No rows found in {csv_path}")

    depth_m = np.array([float(row["depth_m"]) for row in rows], dtype=float)
    deposited_energy_j_per_m3 = np.array(
        [float(row["deposited_energy_j_per_m3"]) for row in rows], dtype=float
    )
    electric_field_v_per_m = np.array([float(row["electric_field_v_per_m"]) for row in rows], dtype=float)

    order = np.argsort(depth_m)
    depth_m = depth_m[order]
    deposited_energy_j_per_m3 = deposited_energy_j_per_m3[order]
    electric_field_v_per_m = electric_field_v_per_m[order]
    return depth_m, deposited_energy_j_per_m3, electric_field_v_per_m


def _compute_depth_edges(depth_centers_m: np.ndarray) -> np.ndarray:
    if depth_centers_m.size < 2:
        half_thickness = max(depth_centers_m[0], 1.0e-6)
        return np.array([max(0.0, depth_centers_m[0] - half_thickness), depth_centers_m[0] + half_thickness])

    dz = np.diff(depth_centers_m)
    edges = np.empty(depth_centers_m.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (depth_centers_m[:-1] + depth_centers_m[1:])
    edges[0] = max(0.0, depth_centers_m[0] - 0.5 * dz[0])
    edges[-1] = depth_centers_m[-1] + 0.5 * dz[-1]
    return edges


def _integrate_potential(depth_centers_m: np.ndarray, electric_field_v_per_m: np.ndarray) -> np.ndarray:
    depth_edges_m = _compute_depth_edges(depth_centers_m)
    layer_thickness_m = np.diff(depth_edges_m)
    cumulative_drop_v = np.cumsum(electric_field_v_per_m * layer_thickness_m)
    # Set top-surface potential to 0 V, then integrate inward.
    return -cumulative_drop_v


def _make_2d_map(profile: np.ndarray, lateral_samples: int) -> np.ndarray:
    return np.tile(profile[:, np.newaxis], (1, lateral_samples))


def _plot_map(
    values_2d: np.ndarray,
    lateral_edges_mm: np.ndarray,
    depth_edges_mm: np.ndarray,
    title: str,
    colorbar_label: str,
    output_png: Path,
    cmap: str,
) -> None:
    figure, axis = plt.subplots(figsize=(8.4, 4.8), constrained_layout=True)
    mesh = axis.pcolormesh(
        lateral_edges_mm,
        depth_edges_mm,
        values_2d,
        shading="auto",
        cmap=cmap,
    )
    axis.invert_yaxis()
    axis.set_xlabel("Lateral position x [mm]")
    axis.set_ylabel("Depth z [mm]")
    axis.set_title(title)
    colorbar = figure.colorbar(mesh, ax=axis)
    colorbar.set_label(colorbar_label)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=220)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render 2D deposition and potential maps from internal charging CSV output."
    )
    parser.add_argument(
        "--input",
        default="results/internal_kapton_100kev_beam.csv",
        help="Internal charging CSV path",
    )
    parser.add_argument(
        "--output-dir",
        default="results",
        help="Directory for generated PNG figures",
    )
    parser.add_argument(
        "--lateral-width-mm",
        type=float,
        default=20.0,
        help="Visualized lateral width in mm for the 2D cross-section",
    )
    parser.add_argument(
        "--lateral-samples",
        type=int,
        default=200,
        help="Number of lateral samples for 2D rendering",
    )
    args = parser.parse_args()

    input_csv = Path(args.input)
    output_dir = Path(args.output_dir)

    depth_m, deposited_energy_j_per_m3, electric_field_v_per_m = _load_profile(input_csv)
    potential_v = _integrate_potential(depth_m, electric_field_v_per_m)

    depth_edges_m = _compute_depth_edges(depth_m)
    depth_edges_mm = depth_edges_m * 1.0e3

    half_width_mm = 0.5 * max(args.lateral_width_mm, 1.0e-6)
    lateral_edges_mm = np.linspace(-half_width_mm, half_width_mm, args.lateral_samples + 1)

    deposition_2d = _make_2d_map(deposited_energy_j_per_m3, args.lateral_samples)
    potential_2d = _make_2d_map(potential_v, args.lateral_samples)

    deposition_png = output_dir / "kapton_100keV_deposition_2d.png"
    potential_png = output_dir / "kapton_100keV_potential_2d.png"

    _plot_map(
        deposition_2d,
        lateral_edges_mm,
        depth_edges_mm,
        title="Kapton under 100 keV electron beam: energy deposition map",
        colorbar_label="Deposited energy density [J/m^3]",
        output_png=deposition_png,
        cmap="inferno",
    )
    _plot_map(
        potential_2d,
        lateral_edges_mm,
        depth_edges_mm,
        title="Kapton internal potential map (derived from E-field)",
        colorbar_label="Potential [V]",
        output_png=potential_png,
        cmap="viridis",
    )

    print(f"Saved: {deposition_png}")
    print(f"Saved: {potential_png}")


if __name__ == "__main__":
    main()
