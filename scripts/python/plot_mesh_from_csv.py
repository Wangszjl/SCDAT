import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


def load_nodes(nodes_csv: Path) -> np.ndarray:
    rows = []
    with nodes_csv.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append((float(row["x"]), float(row["y"]), float(row["z"])))
    if not rows:
        raise ValueError(f"No nodes found in {nodes_csv}")
    return np.array(rows, dtype=float)


def load_triangles(tri_csv: Path):
    if not tri_csv.exists():
        return np.empty((0, 3), dtype=int), np.empty((0,), dtype=int), []
    rows = []
    boundary_ids = []
    boundary_names = []
    with tri_csv.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append((int(row["n1"]), int(row["n2"]), int(row["n3"])))
            boundary_ids.append(int(row.get("boundary_id", "0") or 0))
            boundary_names.append(row.get("boundary_name", "") or "")
    if not rows:
        return np.empty((0, 3), dtype=int), np.empty((0,), dtype=int), []
    return np.array(rows, dtype=int), np.array(boundary_ids, dtype=int), boundary_names


def load_tetrahedra(tet_csv: Path) -> np.ndarray:
    if not tet_csv.exists():
        return np.empty((0, 4), dtype=int)
    rows = []
    with tet_csv.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append((int(row["n1"]), int(row["n2"]), int(row["n3"]), int(row["n4"])))
    if not rows:
        return np.empty((0, 4), dtype=int)
    return np.array(rows, dtype=int)


def boundary_faces_from_tetra(tets: np.ndarray) -> np.ndarray:
    if tets.size == 0:
        return np.empty((0, 3), dtype=int)

    face_counter = Counter()
    oriented = {}

    for tet in tets:
        faces = [
            (tet[0], tet[1], tet[2]),
            (tet[0], tet[1], tet[3]),
            (tet[0], tet[2], tet[3]),
            (tet[1], tet[2], tet[3]),
        ]
        for face in faces:
            key = tuple(sorted(face))
            face_counter[key] += 1
            oriented.setdefault(key, face)

    boundary = [oriented[k] for k, c in face_counter.items() if c == 1]
    if not boundary:
        return np.empty((0, 3), dtype=int)
    return np.array(boundary, dtype=int)


def set_axes_equal(ax, xyz: np.ndarray) -> None:
    mins = xyz.min(axis=0)
    maxs = xyz.max(axis=0)
    centers = (mins + maxs) / 2.0
    radius = (maxs - mins).max() / 2.0

    ax.set_xlim(centers[0] - radius, centers[0] + radius)
    ax.set_ylim(centers[1] - radius, centers[1] + radius)
    ax.set_zlim(centers[2] - radius, centers[2] + radius)


def plot_mesh(
    nodes: np.ndarray,
    tris: np.ndarray,
    tri_boundary_ids: np.ndarray,
    tri_boundary_names,
    out_png: Path,
    max_faces: int = 20000,
) -> None:
    fig = plt.figure(figsize=(10, 8), dpi=150)
    ax = fig.add_subplot(111, projection="3d")

    if tris.size > 0:
        if len(tris) > max_faces:
            idx = np.linspace(0, len(tris) - 1, max_faces, dtype=int)
            tris = tris[idx]
            tri_boundary_ids = tri_boundary_ids[idx] if tri_boundary_ids.size == len(idx) else tri_boundary_ids
            tri_boundary_names = [tri_boundary_names[i] for i in idx] if len(tri_boundary_names) == len(idx) else tri_boundary_names

        unique_ids = sorted(set(tri_boundary_ids.tolist())) if tri_boundary_ids.size > 0 else [0]
        cmap = plt.get_cmap("tab10", max(1, len(unique_ids)))
        color_by_id = {gid: cmap(i) for i, gid in enumerate(unique_ids)}

        legend_items = []
        for gid in unique_ids:
            mask = tri_boundary_ids == gid if tri_boundary_ids.size > 0 else np.ones(len(tris), dtype=bool)
            group_tris = tris[mask]
            if len(group_tris) == 0:
                continue

            polys = nodes[group_tris]
            color = color_by_id[gid]
            mesh = Poly3DCollection(polys, alpha=0.32, edgecolor="k", linewidth=0.15)
            mesh.set_facecolor((color[0], color[1], color[2], 0.30))
            ax.add_collection3d(mesh)

            group_name = ""
            if len(tri_boundary_names) == len(tri_boundary_ids):
                for name, bid in zip(tri_boundary_names, tri_boundary_ids):
                    if bid == gid and name:
                        group_name = name
                        break
            label = f"ID {gid}" if not group_name else f"{group_name} (ID {gid})"
            legend_items.append(Patch(facecolor=(color[0], color[1], color[2], 0.6), edgecolor="k", label=label))

        if legend_items:
            ax.legend(handles=legend_items, loc="upper right", fontsize=8, frameon=True)

    step = max(1, len(nodes) // 6000)
    sampled = nodes[::step]
    ax.scatter(sampled[:, 0], sampled[:, 1], sampled[:, 2], s=0.5, c="tab:red", alpha=0.6)

    set_axes_equal(ax, nodes)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Mesh Visualization (physical groups / boundary IDs)")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot mesh from exported CSV files.")
    parser.add_argument("--input-dir", required=True, help="Directory containing nodes.csv/triangles.csv/tetrahedra.csv")
    parser.add_argument("--output", required=True, help="Output PNG path")
    parser.add_argument("--max-faces", type=int, default=20000, help="Maximum number of faces to draw")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    nodes = load_nodes(input_dir / "nodes.csv")
    triangles, tri_boundary_ids, tri_boundary_names = load_triangles(input_dir / "triangles.csv")

    if triangles.size == 0:
        tets = load_tetrahedra(input_dir / "tetrahedra.csv")
        triangles = boundary_faces_from_tetra(tets)
        tri_boundary_ids = np.zeros((len(triangles),), dtype=int)
        tri_boundary_names = [""] * len(triangles)

    plot_mesh(
        nodes,
        triangles,
        tri_boundary_ids,
        tri_boundary_names,
        Path(args.output),
        max_faces=args.max_faces,
    )


if __name__ == "__main__":
    main()
