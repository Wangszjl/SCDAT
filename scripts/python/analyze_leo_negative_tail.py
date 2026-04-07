#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# Canonical columns used by SC-002 diagnostics.
TIME_COL = "time_s"
VS_COL = "surface_potential_v"
JE_COL = "electron_current_density_a_per_m2"
JNET_COL = "net_current_density_a_per_m2"

# Candidate emission/collection terms used for dominant contribution ranking.
CONTRIB_COLS = [
    "ion_current_density_a_per_m2",
    "photo_emission_density_a_per_m2",
    "secondary_emission_density_a_per_m2",
    "backscatter_emission_density_a_per_m2",
    "ion_secondary_emission_density_a_per_m2",
    "leakage_current_density_a_per_m2",
    "thermionic_emission_density_a_per_m2",
    "field_emission_density_a_per_m2",
    "ram_ion_current_density_a_per_m2",
    "live_pic_net_collection_density_a_per_m2",
]


@dataclass
class Metrics:
    count: int
    mae: float
    rmse: float
    bias: float
    max_abs_error: float

    def as_dict(self) -> Dict[str, float]:
        return {
            "count": self.count,
            "mae": self.mae,
            "rmse": self.rmse,
            "bias": self.bias,
            "max_abs_error": self.max_abs_error,
        }


def read_numeric_csv(path: Path) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            row: Dict[str, float] = {}
            for key, val in raw.items():
                if key is None:
                    continue
                val = (val or "").strip()
                if val == "":
                    continue
                try:
                    row[key] = float(val)
                except ValueError:
                    # Skip non-numeric fields quietly.
                    continue
            if TIME_COL in row:
                rows.append(row)
    if not rows:
        raise ValueError(f"No numeric rows with '{TIME_COL}' found in {path}")
    rows.sort(key=lambda r: r[TIME_COL])
    return rows


def build_time_index(rows: Sequence[Dict[str, float]]) -> Dict[float, Dict[str, float]]:
    return {row[TIME_COL]: row for row in rows}


def aligned_pairs(
    sim_rows: Sequence[Dict[str, float]],
    ref_rows: Sequence[Dict[str, float]],
    tolerance: float,
) -> List[Tuple[Dict[str, float], Dict[str, float]]]:
    ref_idx = build_time_index(ref_rows)
    ref_times = sorted(ref_idx.keys())

    pairs: List[Tuple[Dict[str, float], Dict[str, float]]] = []
    for srow in sim_rows:
        st = srow[TIME_COL]
        if st in ref_idx:
            pairs.append((srow, ref_idx[st]))
            continue

        # Fallback: nearest-time matching under tolerance.
        nearest_t = min(ref_times, key=lambda t: abs(t - st))
        if abs(nearest_t - st) <= tolerance:
            pairs.append((srow, ref_idx[nearest_t]))

    return pairs


def collect_errors(
    pairs: Sequence[Tuple[Dict[str, float], Dict[str, float]]],
    value_col: str,
    mask: Optional[Sequence[bool]] = None,
) -> List[float]:
    errors: List[float] = []
    for idx, (srow, rrow) in enumerate(pairs):
        if mask is not None and not mask[idx]:
            continue
        if value_col in srow and value_col in rrow:
            errors.append(srow[value_col] - rrow[value_col])
    return errors


def calc_metrics(errors: Sequence[float]) -> Optional[Metrics]:
    if not errors:
        return None
    n = len(errors)
    abs_vals = [abs(e) for e in errors]
    sq = [e * e for e in errors]
    mae = sum(abs_vals) / n
    rmse = math.sqrt(sum(sq) / n)
    bias = sum(errors) / n
    max_abs = max(abs_vals)
    return Metrics(count=n, mae=mae, rmse=rmse, bias=bias, max_abs_error=max_abs)


def build_segments(
    pairs: Sequence[Tuple[Dict[str, float], Dict[str, float]]],
    transition_band_abs_v: float,
    transition_window_points: int,
) -> Dict[str, List[bool]]:
    vs_ref: List[float] = []
    for _, rrow in pairs:
        vs_ref.append(rrow.get(VS_COL, float("nan")))

    n = len(vs_ref)
    global_mask = [True] * n

    # Negative-tail: reference potential is negative.
    negative_mask = [(v < 0.0) if math.isfinite(v) else False for v in vs_ref]

    # Transition: around sign crossing if present; otherwise small-|V| band fallback.
    transition_mask = [False] * n
    crossing_idx: List[int] = []
    for i in range(1, n):
        a, b = vs_ref[i - 1], vs_ref[i]
        if math.isfinite(a) and math.isfinite(b) and ((a <= 0.0 < b) or (a >= 0.0 > b)):
            crossing_idx.append(i)

    if crossing_idx:
        for idx in crossing_idx:
            left = max(0, idx - transition_window_points)
            right = min(n, idx + transition_window_points + 1)
            for j in range(left, right):
                transition_mask[j] = True
    else:
        # Fallback for monotonic/no-crossing traces: take near-zero candidates first,
        # but cap to a small window to keep transition diagnostics informative.
        for i, v in enumerate(vs_ref):
            if math.isfinite(v) and abs(v) <= transition_band_abs_v:
                transition_mask[i] = True

        selected = sum(1 for x in transition_mask if x)
        cap = min(n, max(2 * transition_window_points + 1, 5))
        if selected == 0 or selected > cap:
            transition_mask = [False] * n
            scored: List[Tuple[float, int]] = []
            for i, v in enumerate(vs_ref):
                if math.isfinite(v):
                    scored.append((abs(v), i))
            scored.sort(key=lambda x: x[0])
            for _, idx in scored[:cap]:
                transition_mask[idx] = True

    return {
        "global": global_mask,
        "transition": transition_mask,
        "negative_tail": negative_mask,
    }


def top_contributions(
    pairs: Sequence[Tuple[Dict[str, float], Dict[str, float]]],
    mask: Sequence[bool],
    top_n: int,
) -> Dict[str, List[Dict[str, float]]]:
    mags_sim: Dict[str, List[float]] = {c: [] for c in CONTRIB_COLS}
    mags_diff: Dict[str, List[float]] = {c: [] for c in CONTRIB_COLS}

    for idx, (srow, rrow) in enumerate(pairs):
        if not mask[idx]:
            continue
        for col in CONTRIB_COLS:
            if col in srow:
                mags_sim[col].append(abs(srow[col]))
            if col in srow and col in rrow:
                mags_diff[col].append(abs(srow[col] - rrow[col]))

    def aggregate(data: Dict[str, List[float]]) -> List[Dict[str, float]]:
        ranked: List[Tuple[str, float]] = []
        for col, vals in data.items():
            if vals:
                ranked.append((col, sum(vals) / len(vals)))
        ranked.sort(key=lambda x: x[1], reverse=True)
        return [{"column": c, "mean_abs": v} for c, v in ranked[:top_n]]

    return {
        "sim_dominant": aggregate(mags_sim),
        "mismatch_dominant": aggregate(mags_diff),
    }


def format_metric_block(name: str, metric: Optional[Metrics]) -> List[str]:
    if metric is None:
        return [f"- {name}: no aligned data"]
    return [
        f"- {name}: count={metric.count}, mae={metric.mae:.6e}, rmse={metric.rmse:.6e}, bias={metric.bias:.6e}, max_abs_error={metric.max_abs_error:.6e}"
    ]


def write_markdown(report: Dict[str, object], output_md: Path) -> None:
    output_md.parent.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    lines.append("# LEO Negative Tail Segmented Diagnostics")
    lines.append("")
    lines.append(f"- simulation: {report['simulation']}")
    lines.append(f"- reference: {report['reference']}")
    lines.append(f"- aligned_pairs: {report['aligned_pairs']}")
    lines.append(f"- time_tolerance_s: {report['time_tolerance_s']}")
    lines.append("")

    seg_reports = report.get("segments", {})
    for seg_name in ["global", "transition", "negative_tail"]:
        seg = seg_reports.get(seg_name, {}) if isinstance(seg_reports, dict) else {}
        lines.append(f"## Segment: {seg_name}")
        lines.append("")
        lines.append(f"- sample_count: {seg.get('sample_count', 0)}")
        for key in [VS_COL, JE_COL, JNET_COL]:
            metric = seg.get("errors", {}).get(key)
            if metric:
                lines.extend(
                    format_metric_block(
                        key,
                        Metrics(
                            count=int(metric["count"]),
                            mae=float(metric["mae"]),
                            rmse=float(metric["rmse"]),
                            bias=float(metric["bias"]),
                            max_abs_error=float(metric["max_abs_error"]),
                        ),
                    )
                )
            else:
                lines.append(f"- {key}: no aligned data")

        top = seg.get("dominant_contributions", {})
        lines.append("")
        lines.append("### Dominant Contributions (Simulation)")
        sim_dom = top.get("sim_dominant", []) if isinstance(top, dict) else []
        if sim_dom:
            for item in sim_dom:
                lines.append(f"- {item['column']}: mean_abs={item['mean_abs']:.6e}")
        else:
            lines.append("- (none)")

        lines.append("")
        lines.append("### Dominant Contribution Mismatch (Sim-Ref)")
        mismatch = top.get("mismatch_dominant", []) if isinstance(top, dict) else []
        if mismatch:
            for item in mismatch:
                lines.append(f"- {item['column']}: mean_abs={item['mean_abs']:.6e}")
        else:
            lines.append("- (none)")

        lines.append("")

    output_md.write_text("\n".join(lines), encoding="utf-8")


def diagnose(
    simulation: Path,
    reference: Path,
    time_tolerance_s: float,
    transition_band_abs_v: float,
    transition_window_points: int,
    top_n: int,
) -> Dict[str, object]:
    sim_rows = read_numeric_csv(simulation)
    ref_rows = read_numeric_csv(reference)
    pairs = aligned_pairs(sim_rows, ref_rows, tolerance=time_tolerance_s)

    if not pairs:
        raise ValueError("No aligned time pairs between simulation and reference")

    segments = build_segments(
        pairs,
        transition_band_abs_v=transition_band_abs_v,
        transition_window_points=transition_window_points,
    )

    seg_reports: Dict[str, object] = {}
    for seg_name, mask in segments.items():
        sample_count = sum(1 for x in mask if x)
        errors: Dict[str, object] = {}
        for col in [VS_COL, JE_COL, JNET_COL]:
            metric = calc_metrics(collect_errors(pairs, col, mask=mask))
            if metric is not None:
                errors[col] = metric.as_dict()
        seg_reports[seg_name] = {
            "sample_count": sample_count,
            "errors": errors,
            "dominant_contributions": top_contributions(pairs, mask, top_n=top_n),
        }

    return {
        "simulation": simulation.as_posix(),
        "reference": reference.as_posix(),
        "aligned_pairs": len(pairs),
        "time_tolerance_s": time_tolerance_s,
        "transition_band_abs_v": transition_band_abs_v,
        "transition_window_points": transition_window_points,
        "segments": seg_reports,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="SC-002: segmented diagnostics for LEO negative-tail behavior (global/transition/negative_tail)."
    )
    parser.add_argument("--simulation", required=True, help="Simulation CSV path")
    parser.add_argument("--reference", required=True, help="Reference CSV path")
    parser.add_argument(
        "--output-json",
        default="results/analysis/leo_negative_tail_diagnostics.json",
        help="Output JSON report path",
    )
    parser.add_argument(
        "--output-md",
        default="results/analysis/leo_negative_tail_diagnostics.md",
        help="Output Markdown report path",
    )
    parser.add_argument(
        "--time-tolerance-s",
        type=float,
        default=1e-9,
        help="Time alignment tolerance in seconds for nearest match",
    )
    parser.add_argument(
        "--transition-band-abs-v",
        type=float,
        default=2.0,
        help="Fallback transition band when no sign crossing is found",
    )
    parser.add_argument(
        "--transition-window-points",
        type=int,
        default=3,
        help="Window size around each zero-crossing index",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Top N dominant contribution terms to report",
    )
    args = parser.parse_args()

    sim_path = Path(args.simulation).resolve()
    ref_path = Path(args.reference).resolve()
    out_json = Path(args.output_json)
    out_md = Path(args.output_md)

    if not out_json.is_absolute():
        out_json = Path.cwd() / out_json
    if not out_md.is_absolute():
        out_md = Path.cwd() / out_md

    report = diagnose(
        simulation=sim_path,
        reference=ref_path,
        time_tolerance_s=args.time_tolerance_s,
        transition_band_abs_v=args.transition_band_abs_v,
        transition_window_points=args.transition_window_points,
        top_n=args.top_n,
    )

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    write_markdown(report, out_md)

    print(f"output_json={out_json}")
    print(f"output_md={out_md}")
    print(f"aligned_pairs={report['aligned_pairs']}")
    segs = report["segments"]
    for seg_name in ["global", "transition", "negative_tail"]:
        cnt = segs[seg_name]["sample_count"]
        print(f"segment_{seg_name}_count={cnt}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
