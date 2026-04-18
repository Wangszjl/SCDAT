#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys

from spis_surface_import_lib import build_plasma_wake_documents


def run_command(command: list[str], cwd: pathlib.Path) -> None:
    completed = subprocess.run(command, cwd=cwd, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {completed.returncode}: {' '.join(command)}")


def locate_scdat(repo_root: pathlib.Path) -> pathlib.Path:
    candidate = repo_root / "build_codex/bin/SCDAT.exe"
    if not candidate.exists():
        raise FileNotFoundError(f"Unable to find SCDAT executable at {candidate}")
    return candidate


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the full Plasma_Wake SPIS import -> SCDAT -> comparison workflow."
    )
    parser.add_argument(
        "--case-root",
        type=pathlib.Path,
        default=pathlib.Path(
            "ref/SPIS/Plasma_Wake-6.0.5-SNAPSHOT-2019-12-09_16.28.23.695/Plasma_Wake.spis5"
        ),
        help="SPIS Plasma_Wake .spis5 root.",
    )
    parser.add_argument(
        "--output-root",
        type=pathlib.Path,
        default=pathlib.Path("results/spis_import/plasma_wake"),
        help="Workflow artifact output root.",
    )
    parser.add_argument(
        "--repo-root",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parents[3],
        help="Repository root containing build_codex/bin/SCDAT.exe.",
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    case_root = args.case_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    artifacts = build_plasma_wake_documents(case_root, output_root)
    scdat_exe = locate_scdat(repo_root)

    output_csv = output_root / "plasma_wake.csv"
    run_command([str(scdat_exe), "surface-config", str(artifacts.config_pwl_path), str(output_csv)], repo_root)

    augment_script = repo_root / "Toolkit/Surface Charging/scripts/augment_spis_source_resolved_output.py"
    run_command([sys.executable, str(augment_script), str(artifacts.config_pwl_path), str(output_csv)], repo_root)

    compare_script = repo_root / "Toolkit/Surface Charging/scripts/compare_spis_java_output.py"
    run_command([sys.executable, str(compare_script), str(artifacts.config_pwl_path), str(output_csv)], repo_root)

    summary = {
        "config_path": artifacts.config_pwl_path.as_posix(),
        "output_csv": output_csv.as_posix(),
        "comparison_csv_path": (output_root / "plasma_wake.comparison.csv").as_posix(),
        "comparison_summary_path": (output_root / "plasma_wake.comparison.md").as_posix(),
        "comparison_json_path": (output_root / "plasma_wake.comparison.json").as_posix(),
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
