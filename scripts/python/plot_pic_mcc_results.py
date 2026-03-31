import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REQUIRED_COLUMNS = (
    "z_m",
    "potential_v",
    "electric_field_z_v_per_m",
    "electron_density_m3",
    "ion_density_m3",
    "electron_temperature_ev",
)


def load_profile(csv_path: Path) -> dict[str, np.ndarray]:
    rows: list[dict[str, float]] = []
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV header missing in {csv_path}")

        missing = [name for name in REQUIRED_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise ValueError(f"Missing required columns in {csv_path}: {', '.join(missing)}")

        for row in reader:
            rows.append({name: float(row[name]) for name in REQUIRED_COLUMNS})

    if not rows:
        raise ValueError(f"No profile rows found in {csv_path}")

    return {name: np.array([row[name] for row in rows], dtype=float) for name in REQUIRED_COLUMNS}


def positive_floor(values: np.ndarray, floor: float = 1.0e6) -> np.ndarray:
    return np.maximum(values, floor)


def make_title(csv_path: Path, custom_title: str | None) -> str:
    if custom_title:
        return custom_title
    return f"PIC-MCC Profile: {csv_path.stem}"


def plot_profile(profile: dict[str, np.ndarray], csv_path: Path, output_path: Path,
                 title: str | None = None, dpi: int = 160, show: bool = False) -> None:
    z_mm = profile["z_m"] * 1.0e3
    potential = profile["potential_v"]
    electric_field_kv_per_m = profile["electric_field_z_v_per_m"] / 1.0e3
    electron_density = profile["electron_density_m3"]
    ion_density = profile["ion_density_m3"]
    electron_temperature = profile["electron_temperature_ev"]
    net_charge_density = ion_density - electron_density

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), dpi=dpi, constrained_layout=True)
    fig.suptitle(make_title(csv_path, title), fontsize=15, fontweight="bold")

    ax = axes[0, 0]
    ax.plot(z_mm, potential, color="#1f77b4", linewidth=2.0)
    ax.set_title("Potential")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("Potential [V]")
    ax.grid(True, alpha=0.25)

    ax = axes[0, 1]
    ax.plot(z_mm, electric_field_kv_per_m, color="#d62728", linewidth=2.0)
    ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.35)
    ax.set_title("Electric Field")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("Ez [kV/m]")
    ax.grid(True, alpha=0.25)

    ax = axes[1, 0]
    ax.semilogy(z_mm, positive_floor(electron_density), label="Electrons", color="#2ca02c",
                linewidth=2.0)
    ax.semilogy(z_mm, positive_floor(ion_density), label="Ions", color="#ff7f0e",
                linewidth=2.0)
    ax.set_title("Species Density")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("Density [m^-3]")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=True)

    ax = axes[1, 1]
    ax.plot(z_mm, electron_temperature, color="#9467bd", linewidth=2.0,
            label="Electron temperature")
    ax.set_title("Electron Temperature / Net Charge")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("Te [eV]")
    ax.grid(True, alpha=0.25)

    ax2 = ax.twinx()
    ax2.plot(z_mm, net_charge_density, color="#8c564b", linewidth=1.6, linestyle="--",
             label="ni - ne")
    ax2.set_ylabel("Net density [m^-3]")

    lines = ax.get_lines() + ax2.get_lines()
    labels = [line.get_label() for line in lines]
    ax.legend(lines, labels, frameon=True, loc="best")

    summary = (
        f"z range: {z_mm.min():.3f} to {z_mm.max():.3f} mm\n"
        f"phi range: {potential.min():.3e} to {potential.max():.3e} V\n"
        f"|Ez| max: {np.max(np.abs(electric_field_kv_per_m)):.3e} kV/m\n"
        f"ne max: {electron_density.max():.3e} m^-3\n"
        f"ni max: {ion_density.max():.3e} m^-3\n"
        f"Te max: {electron_temperature.max():.3e} eV"
    )
    fig.text(0.015, 0.02, summary, fontsize=9, family="monospace",
             bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8, "edgecolor": "#cccccc"})

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)

    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot PIC-MCC profile CSV results.")
    parser.add_argument("--input", required=True, help="Input profile CSV path")
    parser.add_argument("--output", help="Output image path; defaults to <input>.png")
    parser.add_argument("--title", help="Custom figure title")
    parser.add_argument("--dpi", type=int, default=160, help="Output image DPI")
    parser.add_argument("--show", action="store_true", help="Display the figure interactively")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else input_path.with_suffix(".png")

    profile = load_profile(input_path)
    plot_profile(profile, input_path, output_path, title=args.title, dpi=args.dpi,
                 show=args.show)


if __name__ == "__main__":
    main()
