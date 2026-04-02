#!/usr/bin/env python3
"""Plot 3D/1D Kapton potential comparison at a target time."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.colors import Normalize


TIME_COLUMN = "time_s"
POTENTIAL_COLUMN = "surface_potential_v"
BODY_COLUMN = "body_potential_v"
FLOATING_COLUMN = "floating_equilibrium_potential_v"


def nearest_row(df: pd.DataFrame, target_time: float) -> pd.Series:
    idx = (df[TIME_COLUMN] - target_time).abs().idxmin()
    return df.loc[idx]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Kapton 2000s comparison figure")
    parser.add_argument(
        "--advanced",
        default="results/surface_geo_ecss_kapton_pic_circuit.csv",
        help="Advanced-framework CSV path",
    )
    parser.add_argument(
        "--reference",
        default="results/surface_geo_ecss_kapton_ref.csv",
        help="1D-reference CSV path",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=2000.0,
        help="Target comparison time in seconds",
    )
    parser.add_argument(
        "--output",
        default="results/kapton_2000s_3d_1d_comparison.png",
        help="Output PNG path",
    )
    parser.add_argument(
        "--csv-output",
        default="results/kapton_2000s_comparison_values.csv",
        help="Output CSV for extracted values",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    advanced_path = Path(args.advanced)
    reference_path = Path(args.reference)
    output_path = Path(args.output)
    csv_output_path = Path(args.csv_output)

    advanced_df = pd.read_csv(advanced_path)
    reference_df = pd.read_csv(reference_path)

    adv_row = nearest_row(advanced_df, args.time)
    ref_row = nearest_row(reference_df, args.time)

    adv_time = float(adv_row[TIME_COLUMN])
    ref_time = float(ref_row[TIME_COLUMN])
    adv_surface = float(adv_row[POTENTIAL_COLUMN])
    ref_surface = float(ref_row[POTENTIAL_COLUMN])
    adv_body = float(adv_row[BODY_COLUMN]) if BODY_COLUMN in adv_row else 0.0
    ref_body = float(ref_row[BODY_COLUMN]) if BODY_COLUMN in ref_row else 0.0
    adv_floating = float(adv_row[FLOATING_COLUMN]) if FLOATING_COLUMN in adv_row else float("nan")
    ref_floating = float(ref_row[FLOATING_COLUMN]) if FLOATING_COLUMN in ref_row else float("nan")

    delta = adv_surface - ref_surface

    fig = plt.figure(figsize=(16, 5.8), dpi=150)

    ax_3d = fig.add_subplot(1, 3, 1, projection="3d")
    x = np.linspace(-0.5, 0.5, 60)
    y = np.linspace(-0.3, 0.3, 36)
    xx, yy = np.meshgrid(x, y)
    zz = np.zeros_like(xx)

    vmin = min(adv_surface, ref_surface, -16000.0)
    vmax = max(0.0, adv_surface, ref_surface)
    norm = Normalize(vmin=vmin, vmax=vmax)
    face_values = np.full_like(xx, adv_surface, dtype=float)

    ax_3d.plot_surface(
        xx,
        yy,
        zz,
        facecolors=cm.coolwarm(norm(face_values)),
        linewidth=0.0,
        antialiased=True,
        shade=False,
    )
    ax_3d.set_title(f"3D Kapton surface potential\nAdvanced @ t={adv_time:.1f}s")
    ax_3d.set_xlabel("x (normalized)")
    ax_3d.set_ylabel("y (normalized)")
    ax_3d.set_zlabel("z")
    ax_3d.set_zlim(-0.2, 0.2)
    ax_3d.view_init(elev=28, azim=-50)
    mappable = cm.ScalarMappable(norm=norm, cmap=cm.coolwarm)
    mappable.set_array([])
    cb = fig.colorbar(mappable, ax=ax_3d, fraction=0.045, pad=0.04)
    cb.set_label("Potential (V)")

    ax_1d = fig.add_subplot(1, 3, 2)
    ax_1d.plot(
        reference_df[TIME_COLUMN],
        reference_df[POTENTIAL_COLUMN],
        color="#1f77b4",
        linewidth=1.7,
        label="1D reference",
    )
    ax_1d.plot(
        advanced_df[TIME_COLUMN],
        advanced_df[POTENTIAL_COLUMN],
        "o-",
        color="#d62728",
        markersize=3.5,
        linewidth=1.2,
        label="Advanced framework",
    )
    ax_1d.scatter([ref_time], [ref_surface], color="#1f77b4", s=36)
    ax_1d.scatter([adv_time], [adv_surface], color="#d62728", s=36)
    ax_1d.axvline(args.time, color="#555555", linestyle="--", linewidth=1.0)
    ax_1d.set_title("1D surface potential history")
    ax_1d.set_xlabel("time (s)")
    ax_1d.set_ylabel("surface potential (V)")
    ax_1d.grid(alpha=0.28)
    ax_1d.legend(loc="best")

    ax_bar = fig.add_subplot(1, 3, 3)
    labels = ["Advanced", "1D ref"]
    values = [adv_surface, ref_surface]
    colors = ["#d62728", "#1f77b4"]
    bars = ax_bar.bar(labels, values, color=colors, alpha=0.9)
    ax_bar.axhline(0.0, color="#444444", linewidth=1.0)
    ax_bar.set_title(f"2000s surface potential comparison\ndelta(advanced-ref) = {delta:.2f} V")
    ax_bar.set_ylabel("surface potential (V)")
    ax_bar.grid(axis="y", alpha=0.25)

    for bar, value in zip(bars, values):
        ax_bar.text(
            bar.get_x() + bar.get_width() / 2.0,
            value,
            f"{value:.1f}",
            ha="center",
            va="bottom" if value >= 0 else "top",
            fontsize=9,
        )

    fig.suptitle("Kapton thin-film potential: 3D distribution + 1D reference comparison", fontsize=13)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)

    summary = pd.DataFrame(
        [
            {
                "target_time_s": args.time,
                "advanced_time_s": adv_time,
                "reference_time_s": ref_time,
                "advanced_surface_potential_v": adv_surface,
                "reference_surface_potential_v": ref_surface,
                "advanced_minus_reference_v": delta,
                "advanced_body_potential_v": adv_body,
                "reference_body_potential_v": ref_body,
                "advanced_floating_equilibrium_potential_v": adv_floating,
                "reference_floating_equilibrium_potential_v": ref_floating,
            }
        ]
    )
    summary.to_csv(csv_output_path, index=False)

    print(f"Saved figure: {output_path}")
    print(f"Saved values: {csv_output_path}")
    print(
        " | ".join(
            [
                f"advanced_surface={adv_surface:.6f}",
                f"reference_surface={ref_surface:.6f}",
                f"delta={delta:.6f}",
                f"advanced_floating={adv_floating:.6f}",
                f"reference_floating={ref_floating:.6f}",
            ]
        )
    )


if __name__ == "__main__":
    main()
