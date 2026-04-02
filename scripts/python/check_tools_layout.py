#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path


HEADER_EXTENSIONS = {".h", ".hh", ".hpp", ".hxx"}
SOURCE_EXTENSIONS = {".c", ".cc", ".cpp", ".cxx"}
CODE_EXTENSIONS = HEADER_EXTENSIONS | SOURCE_EXTENSIONS


@dataclass
class Violation:
    path: Path
    message: str


def is_code_file(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in CODE_EXTENSIONS


def classify_allowed_location(path: Path, tools_root: Path) -> str | None:
    """Return None when location is valid, otherwise return violation reason."""
    relative_parts = path.relative_to(tools_root).parts
    parent_parts = relative_parts[:-1]
    suffix = path.suffix.lower()

    if suffix in HEADER_EXTENSIONS:
        if "include" in parent_parts or "test" in parent_parts:
            return None
        return "header file must be under include/ or test/"

    if suffix in SOURCE_EXTENSIONS:
        if "src" in parent_parts or "test" in parent_parts:
            return None
        return "source file must be under src/ or test/"

    return None


def find_violations(tools_root: Path) -> tuple[list[Violation], int, int, int]:
    violations: list[Violation] = []
    header_count = 0
    source_count = 0
    test_source_count = 0

    for path in tools_root.rglob("*"):
        if not is_code_file(path):
            continue

        suffix = path.suffix.lower()
        rel_parts = path.relative_to(tools_root).parts

        if suffix in HEADER_EXTENSIONS:
            header_count += 1
        elif suffix in SOURCE_EXTENSIONS:
            source_count += 1
            if "test" in rel_parts[:-1]:
                test_source_count += 1

        reason = classify_allowed_location(path, tools_root)
        if reason is not None:
            violations.append(Violation(path=path, message=reason))

    return violations, header_count, source_count, test_source_count


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate Tools/** layout: headers in include/, sources in src/, tests in test/."
    )
    parser.add_argument(
        "--tools-root",
        default="Tools",
        help="Path to Tools root directory (default: Tools)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with non-zero code when violations are found.",
    )
    args = parser.parse_args()

    tools_root = Path(args.tools_root).resolve()
    if not tools_root.exists() or not tools_root.is_dir():
        print(f"[tools-layout] invalid Tools root: {tools_root}", file=sys.stderr)
        return 2

    violations, header_count, source_count, test_source_count = find_violations(tools_root)

    print(f"[tools-layout] root={tools_root}")
    print(f"[tools-layout] headers={header_count}")
    print(f"[tools-layout] sources={source_count}")
    print(f"[tools-layout] test_sources={test_source_count}")
    print(f"[tools-layout] violations={len(violations)}")

    if violations:
        for violation in violations:
            print(
                f"[tools-layout][violation] {violation.path.relative_to(tools_root)}: {violation.message}",
                file=sys.stderr,
            )
        if args.strict:
            return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
