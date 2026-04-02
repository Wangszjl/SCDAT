#!/usr/bin/env python3
import argparse
import csv
import math
from collections import defaultdict
from typing import Dict, Iterable, List, Tuple


def rmse(pairs):
    if not pairs:
        return 0.0
    return math.sqrt(sum((a - b) ** 2 for a, b in pairs) / len(pairs))


def mean_signed_delta(pairs):
    if not pairs:
        return 0.0
    return sum(a - b for a, b in pairs) / len(pairs)


def initial_efolding_energy(rows, field):
    if len(rows) < 2:
        return 0.0
    base = rows[0]
    base_current = abs(float(base[f"actual_{field}"]))
    if base_current <= 1.0e-30:
        return 0.0
    base_potential = float(base["actual_potential_v"])
    for row in rows[1:]:
        potential = float(row["actual_potential_v"])
        delta_v = abs(potential - base_potential)
        if delta_v <= 0.0 or potential >= base_potential:
            continue
        current = abs(float(row[f"actual_{field}"]))
        if current <= 1.0e-30 or current >= base_current:
            continue
        ratio = current / base_current
        if ratio <= 0.0 or ratio >= 1.0:
            continue
        return delta_v / abs(math.log(ratio))
    return 0.0


def initial_delta(rows, field):
    if not rows:
        return 0.0
    first = rows[0]
    return float(first[f"actual_{field}"]) - float(first[f"reference_{field}"])


def initial_ratio(rows, field):
    if not rows:
        return 0.0
    first = rows[0]
    denominator = float(first[f"reference_{field}"])
    if abs(denominator) <= 1.0e-30:
        return 0.0
    return float(first[f"actual_{field}"]) / denominator


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def negative_tail_threshold(rows: List[Dict[str, str]]) -> float:
    potentials = [to_float(r.get("reference_potential_v", 0.0)) for r in rows]
    if not potentials:
        return 0.0
    minimum = min(potentials)
    if minimum >= 0.0:
        return 0.0
    return min(-1.0, 0.5 * minimum)


def segment_name(reference_potential_v: float, tail_threshold_v: float) -> str:
    if reference_potential_v <= tail_threshold_v:
        return "negative_tail"
    if reference_potential_v < -1.0:
        return "negative_mid"
    if reference_potential_v <= 1.0:
        return "near_zero"
    return "positive"


def collect_pairs(rows: Iterable[Dict[str, str]], field: str) -> List[Tuple[float, float]]:
    pairs: List[Tuple[float, float]] = []
    for row in rows:
        pairs.append((to_float(row.get(f"actual_{field}", 0.0)),
                      to_float(row.get(f"reference_{field}", 0.0))))
    return pairs


def summarize_segment(rows: List[Dict[str, str]], field: str) -> Dict[str, float]:
    pairs = collect_pairs(rows, field)
    return {
        "count": float(len(pairs)),
        "rmse": rmse(pairs),
        "mean_delta": mean_signed_delta(pairs),
    }


def render_markdown(role: str,
                    tail_threshold_v: float,
                    summary: Dict[str, Dict[str, Dict[str, float]]]) -> str:
    lines: List[str] = []
    lines.append(f"## {role}")
    lines.append("")
    lines.append(f"- negative_tail_threshold_v: `{tail_threshold_v:.6f}`")
    lines.append("")
    lines.append("| segment | field | count | rmse | mean_delta |")
    lines.append("|---|---:|---:|---:|---:|")
    for segment in ["negative_tail", "negative_mid", "near_zero", "positive", "global"]:
        for field in ["je_a_per_m2", "jnet_a_per_m2"]:
            item = summary[segment][field]
            lines.append(
                f"| {segment} | {field} | {int(item['count'])} | "
                f"{item['rmse']:.6e} | {item['mean_delta']:.6e} |"
            )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize current-component mismatch from a legacy benchmark comparison CSV."
    )
    parser.add_argument("csv_path", help="Path to *.benchmark.csv produced by exportResults")
    parser.add_argument(
        "--report-md",
        default="",
        help="Optional path to write a markdown segment-convergence report.",
    )
    args = parser.parse_args()

    grouped = defaultdict(list)
    with open(args.csv_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            grouped[row["role"]].append(row)

    fields = [
        "jnet_a_per_m2",
        "je_a_per_m2",
        "jse_a_per_m2",
        "jb_a_per_m2",
        "ji_a_per_m2",
        "jsi_a_per_m2",
        "jph_a_per_m2",
        "jcond_a_per_m2",
    ]

    markdown_blocks: List[str] = []
    for role, rows in grouped.items():
        tail_threshold_v = negative_tail_threshold(rows)
        segmented_rows: Dict[str, List[Dict[str, str]]] = defaultdict(list)
        for row in rows:
            segmented_rows[
                segment_name(to_float(row.get("reference_potential_v", 0.0)), tail_threshold_v)
            ].append(row)
        segmented_rows["global"] = rows

        print(f"[{role}]")
        print(f"  negative_tail_threshold_v={tail_threshold_v:.6f}")
        for field in fields:
            pairs = collect_pairs(rows, field)
            print(
                f"  {field}: rmse={rmse(pairs):.6e}, mean_delta={mean_signed_delta(pairs):.6e}"
            )

        segment_summary: Dict[str, Dict[str, Dict[str, float]]] = {}
        print("  segmented_convergence:")
        for segment in ["negative_tail", "negative_mid", "near_zero", "positive", "global"]:
            segment_summary[segment] = {}
            segment_rows = segmented_rows.get(segment, [])
            for field in ["je_a_per_m2", "jnet_a_per_m2"]:
                stats = summarize_segment(segment_rows, field)
                segment_summary[segment][field] = stats
                print(
                    f"    {segment}.{field}: count={int(stats['count'])}, "
                    f"rmse={stats['rmse']:.6e}, mean_delta={stats['mean_delta']:.6e}"
                )

        print(
            "  je_initial_efolding_energy_ev="
            f"{initial_efolding_energy(rows, 'je_a_per_m2'):.6f}"
        )
        print(f"  je_initial_delta={initial_delta(rows, 'je_a_per_m2'):.6e}")
        print(f"  je_initial_ratio={initial_ratio(rows, 'je_a_per_m2'):.6f}")

        markdown_blocks.append(render_markdown(role, tail_threshold_v, segment_summary))

    if args.report_md:
        with open(args.report_md, "w", encoding="utf-8") as handle:
            handle.write("# Legacy Benchmark Segment Convergence Report\n\n")
            handle.write("\n\n".join(markdown_blocks))
            handle.write("\n")
        print(f"markdown_report={args.report_md}")


if __name__ == "__main__":
    main()
