#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pathlib

from spis_surface_import_lib import build_spis_case_documents


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Import a SPIS .spis5 case into reusable surface-config artifacts."
    )
    parser.add_argument("case_root", type=pathlib.Path, help="Path to the .spis5 case root.")
    parser.add_argument(
        "--output-root",
        type=pathlib.Path,
        default=pathlib.Path("results/spis_import/spis_case"),
        help="Directory where generated JSON artifacts will be written.",
    )
    args = parser.parse_args()

    artifacts = build_spis_case_documents(args.case_root, args.output_root)
    summary = {
        "case_root": artifacts.case_root.as_posix(),
        "output_root": artifacts.output_root.as_posix(),
        "config_pwl_path": artifacts.config_pwl_path.as_posix(),
        "discrete_manifest_path": artifacts.discrete_manifest_path.as_posix(),
        "import_manifest_path": artifacts.import_manifest_path.as_posix(),
        "scan_plan_path": artifacts.scan_plan_path.as_posix(),
        "spis_reference_output_root": artifacts.spis_reference_output_root.as_posix(),
        "discrete_config_paths": [path.as_posix() for path in artifacts.discrete_config_paths],
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
