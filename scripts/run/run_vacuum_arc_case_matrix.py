#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple


def deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if (
            isinstance(value, dict)
            and key in result
            and isinstance(result[key], dict)
        ):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def load_json(path: Path) -> Dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid json at {path}: {exc}") from exc


def resolve_scdat_exe(project_root: Path, override: str) -> Path:
    if override:
        exe_path = Path(override)
        return exe_path if exe_path.is_absolute() else (project_root / exe_path)

    candidates = [
        project_root / "build" / "bin" / "SCDAT.exe",
        project_root / "build" / "bin" / "SCDAT",
        project_root / "build_codex" / "bin" / "SCDAT.exe",
        project_root / "build_codex" / "bin" / "SCDAT",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("unable to find SCDAT executable; use --scdat-exe")


def normalize_variants(matrix: Dict[str, Any]) -> List[Dict[str, str]]:
    default = [
        {"id": "arc_pic_aligned", "alignment_mode": "arc_pic_aligned"},
        {"id": "legacy_baseline", "alignment_mode": "legacy_baseline"},
    ]
    raw_variants = matrix.get("variants", default)
    variants: List[Dict[str, str]] = []
    for raw in raw_variants:
        variant_id = str(raw.get("id", "")).strip()
        alignment_mode = str(raw.get("alignment_mode", "")).strip()
        if not variant_id or not alignment_mode:
            raise ValueError(f"invalid variant entry: {raw}")
        variants.append({"id": variant_id, "alignment_mode": alignment_mode})
    if not variants:
        raise ValueError("matrix must contain at least one variant")
    return variants


def materialize_case_config(
    base_config: Dict[str, Any],
    case_entry: Dict[str, Any],
    variant: Dict[str, str],
    case_id: str,
    output_csv: str,
) -> Dict[str, Any]:
    config = copy.deepcopy(base_config)

    if "base_preset" in case_entry:
        config["base_preset"] = case_entry["base_preset"]

    config["name"] = f"{case_id}_{variant['id']}"

    run_object = config.setdefault("run", {})
    run_object["output_csv"] = output_csv

    config_object = config.setdefault("config", {})
    config_object["alignment_mode"] = variant["alignment_mode"]

    if isinstance(case_entry.get("run"), dict):
        config["run"] = deep_merge(run_object, case_entry["run"])
        config["run"]["output_csv"] = output_csv

    if isinstance(case_entry.get("config"), dict):
        config["config"] = deep_merge(config_object, case_entry["config"])
        config["config"]["alignment_mode"] = variant["alignment_mode"]

    return config


def run_one(
    scdat_exe: Path,
    generated_config_path: Path,
    output_csv_path: Path,
) -> Tuple[int, str, str]:
    command = [
        str(scdat_exe),
        "arc-config",
        str(generated_config_path),
        str(output_csv_path),
    ]
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    return completed.returncode, completed.stdout, completed.stderr


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run Vacuum Arc case matrix (case x variant) with normalized output naming."
    )
    parser.add_argument(
        "--matrix-json",
        required=True,
        help="Path to case matrix json.",
    )
    parser.add_argument(
        "--project-root",
        default=".",
        help="Project root path.",
    )
    parser.add_argument(
        "--scdat-exe",
        default="",
        help="Optional explicit path to SCDAT executable.",
    )
    parser.add_argument(
        "--output-root",
        default="results/arc_repro",
        help="Root output directory for generated configs/results.",
    )
    parser.add_argument(
        "--timestamp",
        default="",
        help="Optional fixed timestamp token; default uses UTC yyyymmddTHHMMSSZ.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Generate configs/manifest without executing SCDAT.",
    )
    parser.add_argument(
        "--continue-on-failure",
        action="store_true",
        help="Continue remaining runs when one run fails.",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    matrix_path = Path(args.matrix_json)
    if not matrix_path.is_absolute():
        matrix_path = (project_root / matrix_path).resolve()
    if not matrix_path.is_file():
        print(f"status=FAIL(missing_matrix_json:{matrix_path})")
        return 3

    try:
        matrix = load_json(matrix_path)
    except ValueError as exc:
        print(f"status=FAIL({exc})")
        return 3

    cases = matrix.get("cases", [])
    if not isinstance(cases, list) or not cases:
        print("status=FAIL(matrix_cases_empty)")
        return 2

    try:
        variants = normalize_variants(matrix)
    except ValueError as exc:
        print(f"status=FAIL({exc})")
        return 2

    timestamp = args.timestamp.strip() or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = (project_root / output_root).resolve()
    configs_root = output_root / "generated_configs"

    try:
        scdat_exe = resolve_scdat_exe(project_root, args.scdat_exe)
    except FileNotFoundError as exc:
        print(f"status=FAIL({exc})")
        return 3

    if not args.dry_run and not scdat_exe.is_file():
        print(f"status=FAIL(missing_scdat_exe:{scdat_exe})")
        return 3

    manifest_runs: List[Dict[str, Any]] = []
    has_failures = False

    for raw_case in cases:
        case_id = str(raw_case.get("case_id", "")).strip()
        config_json = str(raw_case.get("config_json", "")).strip()
        if not case_id or not config_json:
            print(f"status=FAIL(invalid_case_entry:{raw_case})")
            return 2

        config_path = Path(config_json)
        if not config_path.is_absolute():
            config_path = (project_root / config_path).resolve()
        if not config_path.is_file():
            print(f"status=FAIL(missing_case_config:{config_path})")
            return 2

        try:
            base_config = load_json(config_path)
        except ValueError as exc:
            print(f"status=FAIL({exc})")
            return 2

        for variant in variants:
            output_csv = (
                output_root
                / case_id
                / f"{case_id}_{variant['id']}_{timestamp}.csv"
            )
            generated_config = (
                configs_root
                / f"{case_id}_{variant['id']}_{timestamp}.json"
            )

            config_payload = materialize_case_config(
                base_config,
                raw_case,
                variant,
                case_id,
                output_csv.as_posix(),
            )

            ensure_parent(generated_config)
            generated_config.write_text(
                json.dumps(config_payload, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )

            run_record: Dict[str, Any] = {
                "case_id": case_id,
                "variant": variant["id"],
                "alignment_mode": variant["alignment_mode"],
                "input_config": str(config_path),
                "generated_config": str(generated_config),
                "output_csv": str(output_csv),
                "status": "DRY_RUN" if args.dry_run else "PENDING",
            }

            if args.dry_run:
                print(
                    f"dry_run case={case_id} variant={variant['id']} output={output_csv}"
                )
                manifest_runs.append(run_record)
                continue

            ensure_parent(output_csv)
            exit_code, stdout_text, stderr_text = run_one(
                scdat_exe,
                generated_config,
                output_csv,
            )
            run_record["exit_code"] = exit_code
            run_record["stdout"] = stdout_text.strip()
            run_record["stderr"] = stderr_text.strip()
            run_record["status"] = "PASS" if exit_code == 0 else "FAIL"

            print(
                f"run case={case_id} variant={variant['id']} exit_code={exit_code} output={output_csv}"
            )

            if exit_code != 0:
                has_failures = True
                if not args.continue_on_failure:
                    manifest_runs.append(run_record)
                    break

            manifest_runs.append(run_record)

        if has_failures and not args.continue_on_failure:
            break

    manifest = {
        "matrix_json": str(matrix_path),
        "project_root": str(project_root),
        "scdat_exe": str(scdat_exe),
        "timestamp": timestamp,
        "output_root": str(output_root),
        "dry_run": args.dry_run,
        "runs": manifest_runs,
        "status": "FAIL" if has_failures else "PASS",
    }

    manifest_path = output_root / f"manifest_{timestamp}.json"
    ensure_parent(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"manifest_json={manifest_path}")

    if has_failures:
        print("status=FAIL")
        return 2

    print("status=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
