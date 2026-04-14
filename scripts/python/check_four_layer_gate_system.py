#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple


LAYER_PATTERNS: Dict[str, str] = {
    "unit": (
        r"^(MatrixTest|Point3DTest|Vector3DTest|CouplingTest|DiagnosticsTest|ExtendedFieldSolverTest|"
        r"PoissonSolverTest|MultiBoundaryElectricFieldSolverTest|MaterialDatabaseTest|OutputTest|"
        r"PICCycleBoundaryTest|ParticleDefinitionsTest|ParticleSourceSamplingTest|HighOrderParticlePusherTest)\."
    ),
    "component": (
        r"^(PlasmaAnalysisSmokeTest|RadiationSmokeTest|RadiationPresetsSmokeTest|"
        r"SurfaceChargingSmokeTest|InternalChargingSmokeTest|InternalChargingRadiationCouplingSmokeTest|"
        r"VacuumArcSmokeTest|VacuumArcPresetSmokeTest)\."
    ),
    "system": (
        r"^(ToolsLayoutCheck|LegacyBenchmarkAcceptanceGate|LegacyReplayCsvConsistencyGate|"
        r"SurfaceLeoUnitLintGate|UnifiedConfigSchemaV1Gate|UnifiedOutputContractV1Gate|"
        r"ReproducibleBuildDependencyLockGate|DocumentationClosurePolicyGate|"
        r"VacuumArcCaseMatrixDryRunGate|VacuumArcArcPicCompareSelfGate|"
        r"InternalRadiationCouplingRegressionGate|InternalRadiationComponentTrendGate|"
        r"RadiationGeant4AlignmentGate|RadiationUqBudgetGate|VacuumArcCrossPlatformConsistencyGate)$"
    ),
    "physics": (
        r"^(Kapton2000sRegressionGate|PlasmaMultiscaleRegressionGate|VacuumArcArcPicCompareSelfGate|"
        r"InternalRadiationCouplingRegressionGate|InternalRadiationComponentTrendGate|"
        r"RadiationGeant4AlignmentGate|RadiationUqBudgetGate|VacuumArcCrossPlatformConsistencyGate)$"
    ),
}

TEST_LINE_RE = re.compile(r"^\s*Test\s+#\d+:\s+(?P<name>.+?)\s*$")


def run_command(command: List[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )


def parse_test_names(ctest_n_output: str) -> List[str]:
    names: List[str] = []
    for line in ctest_n_output.splitlines():
        match = TEST_LINE_RE.match(line)
        if match:
            names.append(match.group("name").strip())
    return names


def run_layer(
    layer_name: str,
    regex: str,
    project_root: Path,
    build_dir: Path,
    config: str,
) -> Dict[str, Any]:
    discover_cmd = [
        "ctest",
        "--test-dir",
        str(build_dir),
        "-N",
        "-R",
        regex,
    ]
    discover = run_command(discover_cmd, project_root)
    discovered_names = parse_test_names(discover.stdout)

    layer_report: Dict[str, Any] = {
        "layer": layer_name,
        "regex": regex,
        "discover_exit_code": discover.returncode,
        "selected_tests": discovered_names,
        "selected_count": len(discovered_names),
        "status": "PASS",
        "failures": [],
    }

    if discover.returncode != 0:
        layer_report["status"] = "FAIL"
        layer_report["failures"].append(f"discover_failed_exit_code={discover.returncode}")
        if discover.stderr.strip():
            layer_report["failures"].append(f"discover_stderr={discover.stderr.strip()}")
        return layer_report

    if not discovered_names:
        layer_report["status"] = "FAIL"
        layer_report["failures"].append("no_tests_selected")
        return layer_report

    execute_cmd = [
        "ctest",
        "--test-dir",
        str(build_dir),
        "-C",
        config,
        "-R",
        regex,
        "--output-on-failure",
    ]
    execute = run_command(execute_cmd, project_root)
    layer_report["execute_exit_code"] = execute.returncode

    if execute.returncode != 0:
        layer_report["status"] = "FAIL"
        layer_report["failures"].append(f"execute_failed_exit_code={execute.returncode}")

    return layer_report


def write_markdown(report: Dict[str, Any], output_md: Path) -> None:
    lines: List[str] = []
    lines.append("# Four-Layer Gate Report")
    lines.append("")
    lines.append(f"- status: {report['status']}")
    lines.append(f"- timestamp_utc: {report['timestamp_utc']}")
    lines.append(f"- build_dir: {report['build_dir']}")
    lines.append(f"- config: {report['config']}")
    lines.append("")
    lines.append("| layer | selected_count | status |")
    lines.append("|---|---:|---|")
    for layer in report["layers"]:
        lines.append(f"| {layer['layer']} | {layer['selected_count']} | {layer['status']} |")

    failures = report.get("failures", [])
    if failures:
        lines.append("")
        lines.append("## Failures")
        lines.append("")
        for failure in failures:
            lines.append(f"- {failure}")

    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="ENG-005 four-layer gate aggregator: unit/component/system/physics"
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--build-dir", default="build_codex")
    parser.add_argument("--config", default="Debug")
    parser.add_argument("--report-json", default="")
    parser.add_argument("--report-md", default="")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    build_dir = Path(args.build_dir)
    if not build_dir.is_absolute():
        build_dir = (project_root / build_dir).resolve()

    report_json = Path(args.report_json) if args.report_json else build_dir / "four_layer_gate_report.json"
    if not report_json.is_absolute():
        report_json = (project_root / report_json).resolve()

    report_md = Path(args.report_md) if args.report_md else build_dir / "four_layer_gate_report.md"
    if not report_md.is_absolute():
        report_md = (project_root / report_md).resolve()

    try:
        layers_report: List[Dict[str, Any]] = []
        failures: List[str] = []

        for layer_name in ["unit", "component", "system", "physics"]:
            layer = run_layer(
                layer_name=layer_name,
                regex=LAYER_PATTERNS[layer_name],
                project_root=project_root,
                build_dir=build_dir,
                config=args.config,
            )
            layers_report.append(layer)
            if layer["status"] != "PASS":
                for failure in layer.get("failures", []):
                    failures.append(f"{layer_name}:{failure}")

        status = "PASS" if not failures else "FAIL"
        report: Dict[str, Any] = {
            "status": status,
            "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "project_root": str(project_root),
            "build_dir": str(build_dir),
            "config": args.config,
            "layers": layers_report,
            "failures": failures,
        }

        report_json.parent.mkdir(parents=True, exist_ok=True)
        report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        write_markdown(report, report_md)

        print(f"report_json={report_json}")
        print(f"report_md={report_md}")
        print(f"status={status}")
        return 0 if status == "PASS" else 2

    except Exception as exc:  # pylint: disable=broad-except
        print(f"status=FAIL({exc})")
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
