#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import pathlib


def parse_float(values: dict[str, str], key: str, fallback: float = 0.0) -> float:
    text = values.get(key, "")
    if not text:
        return fallback
    try:
        return float(text)
    except ValueError:
        return fallback


def evaluate_acceptance(values: dict[str, str]) -> dict[str, object]:
    contract_id = values.get("acceptance_contract_id", "none")
    applicable = contract_id != "none"
    checks: list[tuple[str, str, float, float, bool]] = []
    failures: list[str] = []
    failed_keys: list[str] = []

    def eval_pair(metric_key: str, threshold_key: str) -> None:
        metric = parse_float(values, metric_key, 0.0)
        threshold = parse_float(values, threshold_key, 0.0)
        if not (threshold > 0.0) or not math.isfinite(threshold):
            return
        passed = math.isfinite(metric) and metric <= threshold
        checks.append((metric_key, threshold_key, metric, threshold, passed))
        if not passed:
            failed_keys.append(metric_key)
            failures.append(f"{metric_key}(actual={metric:.6e}, max={threshold:.6e})")

    if applicable:
        eval_pair("patch_rmse_v", "acceptance_patch_rmse_v_max")
        eval_pair("body_rmse_v", "acceptance_body_rmse_v_max")
        eval_pair("body_je_rmse_a_per_m2", "acceptance_body_je_rmse_a_per_m2_max")
        eval_pair("body_jnet_rmse_a_per_m2", "acceptance_body_jnet_rmse_a_per_m2_max")

        tail_samples = parse_float(values, "body_negative_tail_sample_count", 0.0)
        if tail_samples > 0.0:
            eval_pair(
                "body_negative_tail_je_rmse_a_per_m2",
                "acceptance_negative_tail_body_je_rmse_a_per_m2_max",
            )
            eval_pair(
                "body_negative_tail_jnet_rmse_a_per_m2",
                "acceptance_negative_tail_body_jnet_rmse_a_per_m2_max",
            )

    checks_total = len(checks)
    checks_failed = len(failures)
    if not applicable:
        status = "NOT_APPLICABLE"
    else:
        status = "PASS" if checks_failed == 0 else "FAIL"

    return {
        "applicable": applicable,
        "status": status,
        "checks_total": checks_total,
        "checks_failed": checks_failed,
        "failed_keys": "none" if not failed_keys else ",".join(failed_keys),
        "failure_details": "none" if not failures else "; ".join(failures),
        "checks": checks,
    }


def render_markdown(path: pathlib.Path, values: dict[str, str], acceptance: dict[str, object]) -> str:
    lines: list[str] = []
    lines.append("# Legacy Benchmark Acceptance Summary")
    lines.append("")
    lines.append(f"- path: `{path}`")
    lines.append(f"- benchmark_source: `{values.get('benchmark_source', 'unknown')}`")
    lines.append(f"- execution_mode: `{values.get('execution_mode', 'unknown')}`")
    lines.append(f"- acceptance_contract_id: `{values.get('acceptance_contract_id', 'none')}`")
    lines.append(f"- acceptance_status: `{acceptance['status']}`")
    lines.append("")
    lines.append("| metric | actual | max | pass |")
    lines.append("|---|---:|---:|---:|")
    for metric_key, _threshold_key, metric, threshold, passed in acceptance["checks"]:
        lines.append(
            f"| {metric_key} | {metric:.6e} | {threshold:.6e} | {int(bool(passed))} |"
        )
    if not acceptance["checks"]:
        lines.append("| none | 0 | 0 | 1 |")
    lines.append("")
    lines.append(f"- failed_metric_keys: `{acceptance['failed_keys']}`")
    lines.append(f"- failure_details: `{acceptance['failure_details']}`")
    lines.append("")
    return "\n".join(lines)


def parse_report(path: pathlib.Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize a SCDAT legacy benchmark sidecar report.")
    parser.add_argument("report_path")
    parser.add_argument(
        "--enforce-acceptance",
        action="store_true",
        help="Return non-zero when acceptance status is FAIL.",
    )
    parser.add_argument(
        "--report-md",
        default="",
        help="Optional path to write a markdown acceptance summary.",
    )
    args = parser.parse_args()

    path = pathlib.Path(args.report_path)
    values = parse_report(path)
    acceptance = evaluate_acceptance(values)

    print(f"path={path}")
    for key in [
        "runtime_route",
        "benchmark_source",
        "execution_mode",
        "baseline_family",
        "baseline_origin",
        "consistency_status",
        "consistency_authority",
        "patch_rmse_v",
        "body_rmse_v",
        "patch_je_rmse_a_per_m2",
        "body_je_rmse_a_per_m2",
        "reconstructed_body_initial_je_a_per_m2",
        "reference_body_initial_je_a_per_m2",
        "body_initial_je_input_reference_ratio",
        "acceptance_contract_version",
        "acceptance_contract_id",
        "acceptance_focus_segment",
        "acceptance_gate_status",
        "acceptance_gate_checks_failed",
        "acceptance_gate_failed_metric_keys",
    ]:
        if key in values:
            print(f"{key}={values[key]}")

    print(f"computed_acceptance_status={acceptance['status']}")
    print(f"computed_acceptance_checks_total={acceptance['checks_total']}")
    print(f"computed_acceptance_checks_failed={acceptance['checks_failed']}")
    print(f"computed_acceptance_failed_metric_keys={acceptance['failed_keys']}")
    print(f"computed_acceptance_failure_details={acceptance['failure_details']}")

    if args.report_md:
        report_path = pathlib.Path(args.report_md)
        report_path.write_text(render_markdown(path, values, acceptance), encoding="utf-8")
        print(f"markdown_report={report_path}")

    if args.enforce_acceptance and acceptance["status"] == "FAIL":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
