#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_profile(path: Path):
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    if not rows:
        raise ValueError(f"No rows in {path}")

    depth_m = np.array([float(row["depth_m"]) for row in rows], dtype=float)
    deposited_energy = np.array([float(row["deposited_energy_j_per_m3"]) for row in rows], dtype=float)
    order = np.argsort(depth_m)
    return depth_m[order], deposited_energy[order]


def load_tracks(path: Path):
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    if not rows:
        return {}

    tracks = {}
    for row in rows:
        track_id = int(row["track_id"])
        tracks.setdefault(track_id, []).append((float(row["x_m"]), float(row["z_m"])))

    for track_id, points in tracks.items():
        points.sort(key=lambda item: item[1])
    return tracks


def depth_edges(depth_centers_m: np.ndarray) -> np.ndarray:
    if depth_centers_m.size == 1:
        half = max(depth_centers_m[0], 1.0e-6)
        return np.array([max(0.0, depth_centers_m[0] - half), depth_centers_m[0] + half], dtype=float)

    dz = np.diff(depth_centers_m)
    edges = np.empty(depth_centers_m.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (depth_centers_m[:-1] + depth_centers_m[1:])
    edges[0] = max(0.0, depth_centers_m[0] - 0.5 * dz[0])
    edges[-1] = depth_centers_m[-1] + 0.5 * dz[-1]
    return edges


def plot_tracks(profile_csv: Path, tracks_csv: Path, output_png: Path, max_tracks: int):
    depth_m, deposited_energy = load_profile(profile_csv)
    tracks = load_tracks(tracks_csv)

    if not tracks:
        raise ValueError(f"No track points in {tracks_csv}")

    depth_edges_mm = depth_edges(depth_m) * 1.0e3

    all_x_mm = []
    for points in tracks.values():
        all_x_mm.extend([p[0] * 1.0e3 for p in points])
    if not all_x_mm:
        all_x_mm = [-10.0, 10.0]
    lateral_max_mm = max(1.0, np.percentile(np.abs(np.array(all_x_mm)), 99.5))
    lateral_edges_mm = np.linspace(-lateral_max_mm, lateral_max_mm, 320)

    deposition_map = np.tile(deposited_energy[:, np.newaxis], (1, lateral_edges_mm.size - 1))

    fig, ax = plt.subplots(figsize=(9.6, 5.6), constrained_layout=True)
    mesh = ax.pcolormesh(
        lateral_edges_mm,
        depth_edges_mm,
        deposition_map,
        shading="auto",
        cmap="magma",
        alpha=0.84,
    )

    selected_ids = sorted(tracks.keys())[: max(1, max_tracks)]
    for track_id in selected_ids:
        points = tracks[track_id]
        x_mm = [point[0] * 1.0e3 for point in points]
        z_mm = [point[1] * 1.0e3 for point in points]
        ax.plot(x_mm, z_mm, color="#66FCF1", linewidth=0.5, alpha=0.5)

    ax.invert_yaxis()
    ax.set_xlabel("Lateral position x [mm]")
    ax.set_ylabel("Depth z [mm]")
    ax.set_title("True Monte Carlo trajectories (exported track CSV)")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Deposited energy density [J/m^3]")

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=240)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot true Monte Carlo trajectories from radiation track CSV")
    parser.add_argument("--profile-csv", required=True, help="Radiation profile csv path")
    parser.add_argument("--tracks-csv", required=True, help="Radiation track csv path")
    parser.add_argument("--output", required=True, help="Output PNG path")
    parser.add_argument("--max-tracks", type=int, default=1500, help="Maximum tracks to draw")
    args = parser.parse_args()

    plot_tracks(Path(args.profile_csv), Path(args.tracks_csv), Path(args.output), args.max_tracks)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
