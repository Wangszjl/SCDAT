#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


def load_csv(path: Path):
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if not fieldnames:
            raise ValueError(f"No columns found in {path}")

        axis_name = fieldnames[0]
        axis_values = []
        series = {name: [] for name in fieldnames[1:]}

        for row in reader:
            axis_values.append(float(row[axis_name]))
            for name in fieldnames[1:]:
                value = row.get(name, "")
                series[name].append(float(value) if value else 0.0)

    return axis_name, axis_values, series


def build_figure(axis_name: str, axis_values, series, title: str):
    names = list(series.keys())
    if not names:
        raise ValueError("No data series available for plotting")

    columns = 2 if len(names) > 1 else 1
    rows = math.ceil(len(names) / columns)
    figure, axes = plt.subplots(rows, columns, figsize=(6.5 * columns, 3.8 * rows), squeeze=False)
    axes_flat = [axis for row in axes for axis in row]

    for axis, name in zip(axes_flat, names):
        axis.plot(axis_values, series[name], linewidth=2.0, color="#0F6CAD")
        axis.set_title(name.replace("_", " "))
        axis.set_xlabel(axis_name)
        axis.grid(True, alpha=0.3)

    for axis in axes_flat[len(names):]:
        axis.axis("off")

    figure.suptitle(title)
    figure.tight_layout()
    return figure


def main():
    parser = argparse.ArgumentParser(description="Plot Toolkit csv outputs.")
    parser.add_argument("--input", required=True, help="Input csv file")
    parser.add_argument("--output", required=True, help="Output png file")
    parser.add_argument("--title", default=None, help="Optional plot title")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    axis_name, axis_values, series = load_csv(input_path)
    title = args.title or input_path.stem
    figure = build_figure(axis_name, axis_values, series, title)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
