#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys
from typing import List


def run_command(command: List[str]) -> None:
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    if completed.stdout:
        print(completed.stdout.strip())
    if completed.stderr:
        print(completed.stderr.strip(), file=sys.stderr)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the SCDAT external field and volume bridge stubs for one surface case."
    )
    parser.add_argument("csv_path", help="Surface case CSV path used to derive sidecar file paths.")
    parser.add_argument(
        "--skip-field",
        action="store_true",
        help="Skip the external field bridge stub.",
    )
    parser.add_argument(
        "--skip-volume",
        action="store_true",
        help="Skip the external volume bridge stub.",
    )
    parser.add_argument(
        "--include-boundary-group-results",
        action="store_true",
        help="Ask the field bridge stub to also emit boundary-group keyed aggregates.",
    )
    args = parser.parse_args()

    csv_path = pathlib.Path(args.csv_path)
    base_dir = pathlib.Path(__file__).resolve().parent
    python_exe = sys.executable

    field_request = csv_path.with_suffix(".field_request.json")
    volume_request = csv_path.with_suffix(".volume_request.json")
    field_stub = base_dir / "run_external_field_bridge_stub.py"
    volume_stub = base_dir / "run_external_volume_bridge_stub.py"

    print(f"csv_path={csv_path}")

    if not args.skip_field:
        if field_request.exists():
            field_command = [python_exe, str(field_stub), str(field_request)]
            if args.include_boundary_group_results:
                field_command.append("--include-boundary-group-results")
            run_command(field_command)
        else:
            print(f"field_request_missing={field_request}")

    if not args.skip_volume:
        if volume_request.exists():
            run_command([python_exe, str(volume_stub), str(volume_request)])
        else:
            print(f"volume_request_missing={volume_request}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
