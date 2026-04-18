#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys

from spis_surface_import_lib import build_child_langmuir_documents


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
        description="Run the full Child_Langmuir SPIS import -> SCDAT -> comparison workflow."
    )
    parser.add_argument(
        "--case-root",
        type=pathlib.Path,
        default=pathlib.Path(
            "ref/SPIS/Child_Langmuir-6.0.5-SNAPSHOT-2019-12-09_16.28.23.695/Child_Langmuir.spis5"
        ),
        help="SPIS Child_Langmuir .spis5 root.",
    )
    parser.add_argument(
        "--output-root",
        type=pathlib.Path,
        default=pathlib.Path("results/spis_import/child_langmuir"),
        help="Workflow artifact output root.",
    )
    parser.add_argument(
        "--repo-root",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parents[3],
        help="Repository root containing build_codex/bin/SCDAT.exe.",
    )
    parser.add_argument(
        "--skip-discrete",
        action="store_true",
        help="Only run the timeline PWL case and skip discrete-point cross-validation runs.",
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    case_root = args.case_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    artifacts = build_child_langmuir_documents(case_root, output_root)
    scdat_exe = locate_scdat(repo_root)

    pwl_csv = output_root / "child_langmuir_pwl.csv"
    run_command([str(scdat_exe), "surface-config", str(artifacts.config_pwl_path), str(pwl_csv)], repo_root)

    discrete_manifest = json.loads(artifacts.discrete_manifest_path.read_text(encoding="utf-8"))
    if not args.skip_discrete:
        for entry in discrete_manifest.get("cases", []):
            run_command(
                [str(scdat_exe), "surface-config", entry["config_path"], entry["output_csv"]],
                repo_root,
            )

    augment_script = repo_root / "Toolkit/Surface Charging/scripts/augment_spis_source_resolved_output.py"
    run_command([sys.executable, str(augment_script), str(artifacts.config_pwl_path), str(pwl_csv)], repo_root)

    compare_script = repo_root / "Toolkit/Surface Charging/scripts/compare_spis_java_output.py"
    compare_command = [sys.executable, str(compare_script), str(artifacts.config_pwl_path), str(pwl_csv)]
    if not args.skip_discrete:
        compare_command.extend(["--discrete-manifest", str(artifacts.discrete_manifest_path)])
    run_command(compare_command, repo_root)

    summary = {
        "config_pwl_path": artifacts.config_pwl_path.as_posix(),
        "pwl_csv": pwl_csv.as_posix(),
        "discrete_manifest_path": artifacts.discrete_manifest_path.as_posix(),
        "comparison_csv_path": (output_root / "child_langmuir.comparison.csv").as_posix(),
        "comparison_summary_path": (output_root / "child_langmuir.comparison.md").as_posix(),
        "comparison_json_path": (output_root / "child_langmuir.comparison.json").as_posix(),
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
