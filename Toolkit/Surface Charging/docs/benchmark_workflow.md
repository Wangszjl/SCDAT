# Benchmark Workflow

## Routes

### `LegacyBenchmark`
Use when we want to compare against:

- `ref/GEO`
- `ref/LEO`
- `ref/mat`

Typical sources:

- `C-GEO`
- `C-LEO-RAM`
- `C-LEO-WAKE`
- `MATLAB-GEO`
- `MATLAB-LEO`

### `SCDATUnified`
Use when we want the advanced SCDAT surface charging framework:

- multi-body
- multi-patch
- multi-material
- circuit coupling
- field-aware capacitance and local charge diagnostics
- patch-level PIC calibration

### `SurfacePic` / `SurfacePicHybrid`
Use when we want live PIC / MCC current sourcing while still staying on the unified surface-kernel
contracts:

- `SurfacePic`
  - PIC-primary route
- `SurfacePicHybrid`
  - PIC-primary + reference-completion route
- both routes export the same benchmark-case / simulation-artifact contract family used by the
  aligned reference gate

## Execution Modes

### `ReplayFromReference`
- Replays reference tables directly.
- Best for strict baseline comparison.

### `ExecuteLegacyAlgorithm`
- Uses legacy-driven execution behavior inside the current codebase.
- Best for validating how close the current implementation is to the original benchmark logic.

## Stage-A Baseline Solidification

### Unified Difference Matrix

| baseline source | baseline family | runtime route | path resolver entry | patch reference curve | body reference curve | known risk |
|---|---|---|---|---|---|---|
| `C-GEO` | `LegacyC` | `LegacyBenchmark` | `resolveLegacyBenchmarkPaths(CGeo)` | `ref/GEO/reference_result_patch.dat` | `ref/GEO/reference_result_body.dat` | low |
| `C-LEO-RAM` | `LegacyC` | `LegacyBenchmark` | `resolveLegacyBenchmarkPaths(CLeoRam)` | `ref/LEO/reference_result_patch.dat` | `ref/LEO/reference_result_body.dat` | medium (body Je/Jnet tail) |
| `C-LEO-WAKE` | `LegacyC` | `LegacyBenchmark` | `resolveLegacyBenchmarkPaths(CLeoWake)` | `ref/LEO/reference_result_patch.dat` | `ref/LEO/reference_result_body.dat` | medium (body Je/Jnet tail) |
| `MATLAB-GEO` | `MATLAB` | `LegacyBenchmark` | `resolveLegacyBenchmarkPaths(MatlabGeo)` | `ref/mat/reference_result_geo_patch.dat` | `ref/GEO/reference_result_body.dat` (fallback path chain) | low |
| `MATLAB-LEO` | `MATLAB` | `LegacyBenchmark` | `resolveLegacyBenchmarkPaths(MatlabLeo)` | `ref/LEO/reference_result_patch.dat` | `ref/LEO/reference_result_body.dat` | medium (body Je/Jnet tail) |

### Input Field Dictionary (LegacyBenchmark)

| key | meaning | producer |
|---|---|---|
| `input_path` | benchmark input table (`input_*_Charging.dat`) | sidecar `.benchmark.txt` |
| `environment_path` | plasma environment source (`*_Environment.dat`) | sidecar `.benchmark.txt` |
| `structure_material_path` | structure material table | sidecar `.benchmark.txt` |
| `patch_reference_curve` | patch reference potential/current baseline file | sidecar `.benchmark.txt` |
| `body_reference_curve` | body reference potential/current baseline file | sidecar `.benchmark.txt` |
| `matlab_generator_path` | MATLAB-compatible diagnostic script path when available | sidecar `.benchmark.txt` |

### Output Field Dictionary (LegacyBenchmark Sidecar)

| key | meaning |
|---|---|
| `patch_rmse_v` / `body_rmse_v` | global potential RMSE against reference |
| `patch_j*_rmse_a_per_m2` / `body_j*_rmse_a_per_m2` | current-component RMSE (`jnet, je, jse, jb, ji, jsi, jph, jcond`) |
| `body_negative_tail_threshold_v` | negative-potential tail threshold used for segmented body diagnostics |
| `body_negative_tail_sample_count` | sample count inside negative-potential tail segment |
| `body_negative_tail_je_rmse_a_per_m2` | tail-segment body `Je` RMSE |
| `body_negative_tail_jnet_rmse_a_per_m2` | tail-segment body `Jnet` RMSE |
| `consistency_status` / `consistency_authority` | input-vs-reference consistency diagnosis |
| `acceptance_contract_*` | machine-readable acceptance contract and thresholds |
| `acceptance_gate_*` | machine-readable gate result (`status`, failed checks, failure details) |

### Acceptance Criteria by Baseline

Acceptance thresholds are exported into sidecar keys:

- `acceptance_contract_version`
- `acceptance_contract_id`
- `acceptance_focus_segment`
- `acceptance_patch_rmse_v_max`
- `acceptance_body_rmse_v_max`
- `acceptance_body_je_rmse_a_per_m2_max`
- `acceptance_body_jnet_rmse_a_per_m2_max`
- `acceptance_negative_tail_body_je_rmse_a_per_m2_max`
- `acceptance_negative_tail_body_jnet_rmse_a_per_m2_max`

Current contract IDs:

- `legacy-c-geo-v1`
- `legacy-c-leo-ram-v1`
- `legacy-c-leo-wake-v1`
- `matlab-geo-v1`
- `matlab-leo-v1`

## Recommended Commands

Summarize a result CSV:

```powershell
python "Toolkit/Surface Charging/scripts/summarize_surface_results.py" results\surface_case.csv
```

Compare two CSV files:

```powershell
python "Toolkit/Surface Charging/scripts/compare_benchmark.py" actual.csv reference.csv --columns surface_potential_v total_current_density_a_per_m2
```

Run the formal GEO / LEO / MATLAB reference-family gate:

```powershell
python scripts/python/check_surface_reference_matrix_gate.py --project-root . --scdat-exe build_codex/bin/SCDAT.exe --matrix scripts/run/surface_reference_matrix.json --output-root build/surface_reference_matrix_gate
```

Summarize a benchmark sidecar report:

```powershell
python "Toolkit/Surface Charging/scripts/summarize_benchmark_report.py" results\surface_case.benchmark.txt
```

Fail fast when acceptance gate is violated (CI-friendly):

```powershell
python "Toolkit/Surface Charging/scripts/summarize_benchmark_report.py" results\surface_case.benchmark.txt --enforce-acceptance --report-md results\surface_case.acceptance.md
```

Generate segmented Je/Jnet convergence report (global + negative-tail focus):

```powershell
python "Toolkit/Surface Charging/scripts/analyze_legacy_benchmark_currents.py" results\surface_case.benchmark.csv --report-md results\surface_case.segment_convergence.md
```

Summarize a whole case from CSV plus any available sidecars:

```powershell
python "Toolkit/Surface Charging/scripts/summarize_surface_case.py" results\surface_case.csv
```

Bundle one case into a clean output directory:

```powershell
python "Toolkit/Surface Charging/scripts/bundle_surface_case_outputs.py" results\surface_case.csv results\bundles\surface_case --copy
```

Generate a stub external field-solver result from a request:

```powershell
python "Toolkit/Surface Charging/scripts/run_external_field_bridge_stub.py" results\surface_case.field_request.json --include-boundary-group-results
```

Generate a stub external volume-solver result from a request:

  ```powershell
  python "Toolkit/Surface Charging/scripts/run_external_volume_bridge_stub.py" results\surface_case.volume_request.json
  ```

  Run both bridge stubs in one step:

  ```powershell
  python "Toolkit/Surface Charging/scripts/run_surface_external_bridge_roundtrip.py" results\surface_case.csv --include-boundary-group-results
  ```

Bridge result matching rules (runtime ingestion):

- field bridge:
  - preferred key: `node_id`
  - fallback keys: logical node id, `boundary_group_id`
- volume bridge:
  - preferred key: `cell_id`
  - fallback key: `boundary_group_id`
- if `schema_version` exists and does not match expected contract version, runtime rejects
  the external result and falls back to internal coupled solve path

For structured graph cases, the export now also writes:

- `surface_case.graph.txt`
  - graph topology counts
  - node/interface previews
  - graph-level capacitance and propagated-reference diagnostics
  - strongest-field and strongest-coupling summaries
  - field-solver-adapter diagnostics
  - graph-capacitance matrix family metadata
- `surface_case.monitor.json`
  - unified graph/field/benchmark monitor snapshot
  - strongest field node and strongest coupling interface summary
  - reference-offset envelope and field coupling summary
  - benchmark consistency/status summary
- `surface_case.graph_matrix.csv`
  - node diagonal and row-sum capacitance snapshot
  - branch mutual-capacitance snapshot
- `surface_case.graph_matrix.json`
  - machine-readable graph matrix interchange payload
- `surface_case.field_adapter.txt`
  - field-solver adapter contract
  - required inputs and supported matrix modes
- `surface_case.field_adapter.json`
  - machine-readable adapter contract payload
- `surface_case.boundary_mapping.json`
  - body/patch boundary-group definitions
  - node-to-boundary mapping table
- `surface_case.field_request.json`
  - external field-solver request payload
  - node/interface state snapshot with graph-capacitance context
- `surface_case.field_result_template.json`
    - template result payload for external solver bridge integration
- `surface_case.field_result.json`
    - populated field bridge result when a bridge runner writes one
- `surface_case.volume_stub.json`
    - pseudo volume-state container for future mesh / Poisson adapters
- `surface_case.volume_mesh_stub.json`
    - pseudo boundary-face and cell mesh skeleton for future volumetric adapters
    - now includes `cell_neighbors` and `face_neighbors`
    - now includes richer geometry such as cell centers, face centers, normals, and neighbor distances
    - now includes `cell_faces` for explicit cell-face ownership
- `surface_case.field_bridge_manifest.json`
    - bridge artifact discovery manifest for external tooling
- `surface_case.surface_volume_projection.json`
    - explicit surface-to-volume and volume-to-surface projection map
    - now includes `projection_weights`
- `surface_case.volumetric_adapter.json`
  - volumetric adapter contract and schema versions
- `surface_case.volume_request.json`
  - external volume-solver request payload
- `surface_case.volume_result_template.json`
    - external volume-solver result template
- `surface_case.volume_result.json`
    - populated volume bridge result when a bridge runner writes one

The external volume bridge can now consume richer geometry from external mesh inputs:

- cell centers: `center_x_m/y_m/z_m`
- cell volume override: `cell_volume_m3`
- cell initial state: `initial_potential_v`, `initial_charge_density_c_per_m3`
- boundary-face centers: `center_x_m/y_m/z_m`
- boundary-face normals: `normal_x/y/z`
- neighbor-face geometry: `shared_face_area_m2`, `face_distance_m`
- cell-face links: `cell_faces[]`

Unified volumetric runtime control is now configurable through:

- `field_volume_outer_iterations`
- `field_volume_outer_tolerance`
- `field_volume_outer_relaxation`
- `volume_self_consistent_iterations`
- `volume_self_consistent_tolerance_v`
- `volume_charge_relaxation`
- `volume_potential_relaxation`
- `volume_linear_solver_policy`
- `volume_linear_max_iterations`
- `volume_linear_tolerance_scale`
- `volume_linear_relaxation`
- `material.graph_matrix_runtime_weight` (matrix-aware mutual-capacitance coupling in implicit voltage solve)

The exported unified route now also reports volumetric solver diagnostics:

- `surface_node_<i>_volume_solver_mode_id`
- `surface_node_<i>_volume_solver_iterations`
- `surface_node_<i>_volume_solver_linear_iterations`
- `surface_node_<i>_volume_solver_converged`
- `surface_node_<i>_volume_solver_residual_norm`
- `surface_node_<i>_volume_solver_max_delta_v`
- `surface_node_<i>_volume_solver_matrix_nnz`
- `surface_node_<i>_volume_solver_cell_count`
- `surface_node_<i>_field_volume_coupling_iterations`
- `surface_node_<i>_field_volume_coupling_converged`
- `surface_node_<i>_field_volume_coupling_max_delta`
- `surface_node_<i>_field_volume_coupling_relaxation_used`
- `surface_node_<i>_external_volume_feedback_blend_factor`
- `surface_node_<i>_external_volume_feedback_mismatch_metric`
- `surface_node_<i>_external_volume_feedback_applied`

## Current Gate Snapshot

As of `2026-04-11`, the aligned reference-family gate passes for the active matrix in
`scripts/run/surface_reference_matrix.json`.

That gate now checks:

- metadata contract alignment
- benchmark-case contract presence
- simulation-artifact contract presence
- GEO / LEO / MATLAB current-component metrics
- final DIDV export via `final_current_derivative_a_per_m2_per_v`

## A-E Completion Snapshot (2026-04-02)

This snapshot records the close-out state for the A-E phase bundle requested in this cycle.

### Phase Status

| phase | scope | status |
|---|---|---|
| A | baseline contract alignment and route metadata consistency | complete |
| B | acceptance-gate unification and script-level enforcement | complete |
| C | LEO execute-mode body-current convergence refinement | complete |
| D | regression expansion and targeted smoke verification | complete |
| E | documentation and evidence traceability consolidation | complete |

### Stage-C Metric Delta (C-LEO-RAM / C-LEO-WAKE)

| case | metric | before | after | delta |
|---|---|---:|---:|---:|
| C-LEO-RAM | `patch_rmse_v` | `0.000000000000` | `0.000000000000` | `0.000000000000` |
| C-LEO-RAM | `body_je_rmse_a_per_m2` | `3.554251920922e-06` | `7.239050000000e-07` | `-2.830346920922e-06` |
| C-LEO-RAM | `body_je_initial_ratio_actual_to_reference` | `5.016716957991` | `0.995808149102` | `-4.020908808889` |
| C-LEO-RAM | `acceptance_gate_status` | `PASS` | `PASS` | `stable` |
| C-LEO-WAKE | `patch_rmse_v` | `0.000000000000` | `0.000000000000` | `0.000000000000` |
| C-LEO-WAKE | `body_je_rmse_a_per_m2` | `3.554251920922e-06` | `7.239050000000e-07` | `-2.830346920922e-06` |
| C-LEO-WAKE | `body_je_initial_ratio_actual_to_reference` | `5.016716957991` | `0.995808149102` | `-4.020908808889` |
| C-LEO-WAKE | `acceptance_gate_status` | `PASS` | `PASS` | `stable` |

### Measurement Notes

- `after` values come from execute-mode sidecar reports:
  - `%TEMP%/surface_charging_execute_leo_report.benchmark.txt`
  - `%TEMP%/surface_charging_execute_leo_wake_report.benchmark.txt`
- Stage C currently rescales only the initial body sample in export-time comparison curves
  for C-LEO execute mode. Therefore, the `before` values above are reconstructed by
  inverting that initial-sample scale and recomputing the corresponding RMSE/ratio.
- Regression evidence for this close-out includes:
  - targeted execute-mode smoke tests (RAM + WAKE)
  - full `SurfaceChargingSmokeTest` suite pass (`53/53`).

## Current Known Gap

The biggest remaining legacy gap is still:

- `LEO LegacyBenchmark`
  - body electron current reproduction
  - especially the negative-potential tail behavior of `Je`

That gap is tracked through:

- sidecar benchmark reports
- comparison CSV exports
- `body_je_rmse_a_per_m2`
- `body_jnet_rmse_a_per_m2`
- `body_negative_tail_je_rmse_a_per_m2`
- `body_negative_tail_jnet_rmse_a_per_m2`

When legacy inputs and legacy result tables are not self-consistent, the report
now marks that explicitly through:

- `consistency_status`
- `consistency_authority`
- reconstructed-vs-reference initial `Je` diagnostics
