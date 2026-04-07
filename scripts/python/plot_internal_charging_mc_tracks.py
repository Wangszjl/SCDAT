#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_internal_profile(csv_path: Path):
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    if not rows:
        raise ValueError(f"No data rows in {csv_path}")

    depth_m = np.array([float(row["depth_m"]) for row in rows], dtype=float)
    deposited_energy_j_per_m3 = np.array(
        [float(row["deposited_energy_j_per_m3"]) for row in rows], dtype=float
    )

    order = np.argsort(depth_m)
    return depth_m[order], deposited_energy_j_per_m3[order]


def depth_edges_from_centers(depth_centers_m: np.ndarray) -> np.ndarray:
    if depth_centers_m.size == 1:
        half = max(depth_centers_m[0], 1.0e-6)
        return np.array([max(0.0, depth_centers_m[0] - half), depth_centers_m[0] + half], dtype=float)

    dz = np.diff(depth_centers_m)
    edges = np.empty(depth_centers_m.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (depth_centers_m[:-1] + depth_centers_m[1:])
    edges[0] = max(0.0, depth_centers_m[0] - 0.5 * dz[0])
    edges[-1] = depth_centers_m[-1] + 0.5 * dz[-1]
    return edges


def build_stopping_depth_sampler(depth_centers_m: np.ndarray, deposited_energy_j_per_m3: np.ndarray):
    edges = depth_edges_from_centers(depth_centers_m)
    layer_thickness_m = np.diff(edges)

    # Convert volumetric deposition to areal deposition weights.
    weights = np.maximum(deposited_energy_j_per_m3 * layer_thickness_m, 0.0)
    if np.all(weights <= 0.0):
        # Fallback for degenerate input.
        z = np.maximum(depth_centers_m, 1.0e-9)
        weights = np.exp(-z / np.max(z))

    cdf = np.cumsum(weights)
    cdf /= cdf[-1]

    def sample_depth(rng: np.random.Generator, n: int) -> np.ndarray:
        u = rng.random(n)
        idx = np.searchsorted(cdf, u, side="left")
        idx = np.clip(idx, 0, depth_centers_m.size - 1)
        z0 = edges[idx]
        z1 = edges[idx + 1]
        return z0 + (z1 - z0) * rng.random(n)

    return sample_depth, edges


def generate_tracks(
    sample_stop_depth,
    n_particles: int,
    lateral_width_mm: float,
    max_depth_m: float,
    seed: int,
):
    rng = np.random.default_rng(seed)

    x0_mm = rng.uniform(-0.5 * lateral_width_mm, 0.5 * lateral_width_mm, n_particles)
    stop_depth_m = np.clip(sample_stop_depth(rng, n_particles), 0.0, max_depth_m)

    # Forward-focused beam with small angular spread.
    theta_rad = rng.normal(0.0, np.deg2rad(3.0), n_particles)
    theta_rad = np.clip(theta_rad, np.deg2rad(-12.0), np.deg2rad(12.0))

    tracks = []
    for i in range(n_particles):
        z_end_m = float(stop_depth_m[i])
        if z_end_m <= 0.0:
            continue

        points = max(10, int(14 + 2400 * z_end_m))
        z_m = np.linspace(0.0, z_end_m, points)

        # Mean drift from incidence angle.
        x_mm = x0_mm[i] + np.tan(theta_rad[i]) * z_m * 1.0e3

        # Add lightweight multiple-scattering random walk that grows with depth.
        step_sigma_mm = (0.005 + 0.08 * (z_m / max(max_depth_m, 1.0e-12)))
        noise_mm = rng.normal(0.0, step_sigma_mm, size=points)
        noise_mm = np.cumsum(noise_mm) / np.sqrt(points)
        x_mm += noise_mm

        tracks.append((x_mm, z_m * 1.0e3))

    return tracks


def make_deposition_background(
    deposited_energy_j_per_m3: np.ndarray,
    depth_edges_mm: np.ndarray,
    lateral_edges_mm: np.ndarray,
) -> np.ndarray:
    values = np.tile(deposited_energy_j_per_m3[:, np.newaxis], (1, lateral_edges_mm.size - 1))
    return values


def plot_tracks_with_deposition(
    tracks,
    deposited_background,
    depth_edges_mm,
    lateral_edges_mm,
    output_png: Path,
    title: str,
):
    fig, ax = plt.subplots(figsize=(9.2, 5.4), constrained_layout=True)

    mesh = ax.pcolormesh(
        lateral_edges_mm,
        depth_edges_mm,
        deposited_background,
        shading="auto",
        cmap="magma",
        alpha=0.85,
    )

    for x_mm, z_mm in tracks:
        ax.plot(x_mm, z_mm, color="#66FCF1", linewidth=0.55, alpha=0.45)

    ax.invert_yaxis()
    ax.set_xlabel("Lateral position x [mm]")
    ax.set_ylabel("Depth z [mm]")
    ax.set_title(title)

    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Deposited energy density [J/m^3]")

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=240)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Monte Carlo style trajectory diagnostic for internal charging profiles."
    )
    parser.add_argument(
        "--input",
        default="results/internal_kapton_100kev_beam.csv",
        help="Input internal charging csv",
    )
    parser.add_argument(
        "--output",
        default="results/kapton_100keV_mc_tracks.png",
        help="Output track png",
    )
    parser.add_argument("--n-particles", type=int, default=500, help="Number of sampled trajectories")
    parser.add_argument("--lateral-width-mm", type=float, default=20.0, help="Displayed slab width")
    parser.add_argument("--seed", type=int, default=20260403, help="RNG seed")
    args = parser.parse_args()

    input_csv = Path(args.input)
    output_png = Path(args.output)

    depth_m, deposited_energy_j_per_m3 = load_internal_profile(input_csv)
    sample_depth, depth_edges_m = build_stopping_depth_sampler(depth_m, deposited_energy_j_per_m3)

    tracks = generate_tracks(
        sample_depth,
        n_particles=max(1, args.n_particles),
        lateral_width_mm=max(1.0e-6, args.lateral_width_mm),
        max_depth_m=float(depth_edges_m[-1]),
        seed=args.seed,
    )

    depth_edges_mm = depth_edges_m * 1.0e3
    half_width_mm = 0.5 * max(1.0e-6, args.lateral_width_mm)
    lateral_edges_mm = np.linspace(-half_width_mm, half_width_mm, 320)

    deposited_background = make_deposition_background(
        deposited_energy_j_per_m3, depth_edges_mm, lateral_edges_mm
    )
    plot_tracks_with_deposition(
        tracks,
        deposited_background,
        depth_edges_mm,
        lateral_edges_mm,
        output_png,
        title="Kapton 100 keV electron beam: MC-style trajectories over deposition",
    )

    print(f"Saved: {output_png}")


if __name__ == "__main__":
    main()
