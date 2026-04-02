# SCDAT Surface Charging

## Overview
This module hosts the SCDAT surface charging framework.

It now supports two long-lived routes:

- `SCDATUnified`
  - default high-level runtime path
  - used for multi-body, multi-patch, multi-material surface charging studies
- `LegacyBenchmark`
  - benchmark path used to compare against the legacy C and MATLAB references in `ref/`

## Core Objects

### `StructureBody`
- Represents a conductive structure node.
- Can be floating, fixed, or externally biased.
- Can connect to other bodies and patches through interfaces.

### `SurfacePatch`
- Represents a refined surface element attached to a body.
- Carries local material, area, orientation, spectra, and emission settings.
- Is the main object for multi-material adjacent-surface studies.

### `PatchInterface`
- Represents `body-body`, `patch-body`, or `patch-patch` links.
- Supports conductance, resistance, bias, and dynamic-conductance behavior.

## Output Philosophy
Exports keep legacy-compatible columns, but add richer SCDAT diagnostics:

- node potentials
- node current decomposition
- node electric field and local charge
- node propagated reference potential
- node field-solver-adapted reference potential
- node graph-capacitance diagonal estimates
- node graph-capacitance row-sum estimates
- node field-solver coupling gains
- node PIC calibration factors and live PIC diagnostics
- interface currents, conductances, voltage drops, branch power, and mutual-capacitance estimates

Implicit circuit advancement now also uses matrix-aware mutual-capacitance coupling
(`C_mutual/dt` terms in the nodal solve), controlled by material scalar property
`graph_matrix_runtime_weight` (default `0.35`).

Sidecar reports are now written alongside the main CSV when applicable:

- `.benchmark.txt`
  - legacy benchmark metrics, consistency diagnostics, and reference comparison previews
- `.benchmark.csv`
  - aligned actual-vs-reference comparison rows
- `.graph.txt`
  - multi-body / multi-patch graph summary, node and interface previews, and graph-level diagnostics
- `.monitor.json`
  - unified graph/field/benchmark monitor snapshot for downstream tooling
- `.graph_matrix.csv`
  - graph-capacitance diagonal, row-sum, and mutual snapshot entries
- `.graph_matrix.json`
  - machine-readable graph-matrix interchange snapshot
- `.field_adapter.txt`
  - field-solver adapter contract and supported graph inputs
- `.field_adapter.json`
  - machine-readable field-adapter contract
- `.boundary_mapping.json`
  - body/patch boundary-group definitions and node-to-boundary mappings
- `.field_request.json`
  - machine-readable external field-solver request snapshot
- `.field_result_template.json`
  - result-template contract for external field-solver bridge inputs
- `.field_result.json`
  - populated field-solver bridge result when a stub or external solver writes one
- `.volume_stub.json`
  - pseudo volume-cell container for future mesh / Poisson integration
- `.volume_mesh_stub.json`
  - pseudo boundary-face and cell mesh skeleton for future volumetric solvers
  - includes `cell_neighbors` and `face_neighbors` for graph-to-mesh adjacency handoff
  - includes `shared_face_area_m2` and `face_distance_m` for neighbor geometry
  - includes `center_x_m/y_m/z_m`, `surface_normal_x/y/z`, and `face_center_x_m/y/z`
  - now includes `cell_faces` so external solvers can reconstruct cell-to-face ownership
- `.field_bridge_manifest.json`
  - one-stop manifest describing bridge artifacts for the case
- `.surface_volume_projection.json`
  - formal surface-to-volume and volume-to-surface projection contract
  - includes explicit `projection_weights`
- `.volumetric_adapter.json`
  - machine-readable volumetric adapter contract
- `.volume_request.json`
  - machine-readable external volume solve request
- `.volume_result_template.json`
  - machine-readable external volume result template
- `.volume_result.json`
  - populated volume-solver bridge result when a stub or external solver writes one
- `.volume_history.json`
  - node-wise time history for `volume_potential_v`, `deposited_charge_c`,
    `poisson_residual_v_m`, and pseudo-volume coupling state

## Benchmark Philosophy
- `LegacyBenchmark` is for reference reproduction and regression.
- `SCDATUnified` is for the advanced predictive framework.
- The two routes intentionally coexist and should not be confused.

## Scripts
Helpful local tools live in [`scripts`](e:\3-Code\1-Cplusplus\PIC-Surface-Charging-master\Toolkit\Surface Charging\scripts):

- `compare_benchmark.py`
- `summarize_surface_results.py`
- `surface_graph_summary.py`
- `analyze_legacy_benchmark_currents.py`
- `reconstruct_legacy_leo_initial_body_current.py`
- `summarize_benchmark_report.py`
- `summarize_surface_case.py`
- `bundle_surface_case_outputs.py`
- `run_external_field_bridge_stub.py`
- `run_external_volume_bridge_stub.py`
- `run_surface_external_bridge_roundtrip.py`

The volumetric bridge now accepts richer external mesh geometry:

- cell centers via `center_x_m/y_m/z_m`
- cell volume overrides via `cell_volume_m3`
- cell initial state via `initial_potential_v` and `initial_charge_density_c_per_m3`
- boundary-face centers via `center_x_m/y_m/z_m`
- boundary-face normals via `normal_x/y/z`
- neighbor geometry via `shared_face_area_m2` and `face_distance_m`

Those fields are not passive metadata. The unified runtime now uses them in:

- control-volume Poisson coefficients
- derived pseudo-cell geometry
- volume-to-surface feedback strength

The default volumetric runtime also now performs a small self-consistent Poisson relaxation
instead of a single one-shot solve, so volume state can settle across iterations before being
projected back to the surface runtime.

For richer external volume cases, the runtime now uses an automatic
`iterative_or_dense_auto` linear-solver path instead of always forcing the dense direct solve.

The relaxation is now explicitly configurable from `SurfaceChargingConfig`:

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

The runtime now also exports solver-side diagnostics, not just field results:

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

Metadata records the last active volumetric solve as well:

- `surface_volume_linear_solver_policy`
- `surface_volume_last_solver_mode`
- `surface_volume_last_solver_iterations`
- `surface_volume_last_solver_linear_iterations`
- `surface_volume_last_solver_converged`
- `surface_volume_last_solver_residual_norm`
- `surface_volume_last_solver_max_delta_v`
- `surface_volume_last_solver_matrix_nnz`
- `surface_volume_last_solver_cell_count`
- `surface_field_volume_outer_iterations`
- `surface_field_volume_outer_tolerance`
- `surface_field_volume_outer_relaxation`
- `surface_field_volume_last_coupling_iterations`
- `surface_field_volume_last_coupling_converged`
- `surface_field_volume_last_coupling_max_delta`
- `surface_field_volume_last_coupling_relaxation_used`
- `surface_volume_external_feedback_last_blend_factor`
- `surface_volume_external_feedback_last_mismatch_metric`
- `surface_volume_external_feedback_last_applied`

The iterative volumetric path now uses a sparse-like row structure internally rather than
scanning every matrix entry on every update, so larger structured cases can grow without
immediately falling back to dense-style linear work everywhere.

The unified runtime also now performs an explicit outer coupling loop between the field adapter
and the volumetric adapter. That means field and volume feedback are iterated together inside one
runtime-state build before the result is used by capacitance, current, and timestep logic.

The higher-order graph/runtime boundary now also exposes:

- `SurfaceGraphCapacitanceMatrixProvider`
  - diagonal estimates
  - mutual estimates
  - row-sum estimates
  - matrix-family metadata
- `SurfaceFieldSolverAdapter`
  - field-solver-adapted reference potential
  - graph-coupling gain diagnostics
  - optional external file bridge for request/result exchange

Legacy benchmark sidecar reports now also include initial-current diagnostics
for tricky cases such as `LEO` body electron current:

- initial `Je` delta
- initial `Je` ratio against the reference
- initial `Je` e-folding energy estimate

## Next Technical Focus
- Continue improving `LEO LegacyBenchmark` body electron current reproduction.
- Continue expanding interface-level and patch-level diagnostics.
- Continue evolving the multi-body multi-material graph toward higher-order field coupling.
- Continue using `.graph.txt` and case-summary tooling to keep higher-order graph changes observable.
