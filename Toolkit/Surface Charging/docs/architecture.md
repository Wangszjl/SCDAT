# SCDAT Surface Charging Architecture

## Reading Guide

This document focuses on runtime architecture and kernel layering.

For the native project-facing workflow, read these first:

- [native_surface_architecture_overview.md](./native_surface_architecture_overview.md)
- [native_surface_mainline_workflow.md](./native_surface_mainline_workflow.md)

For the optional SPIS adapter and comparison path, read:

- [spis_import_compute_workflow.md](./spis_import_compute_workflow.md)

That distinction matters:

- native modeling / configuration / simulation / output is the primary project workflow
- `SPIS import` is an auxiliary adapter and benchmarking path, not the architectural center

## Alignment State
The runtime is now organized around a neutral-named unified surface kernel:

- `surface distribution`
- `material interaction`
- `surface DIDV assembly`
- `surface-to-circuit mapping`
- `implicit circuit kernel`

`SCDATUnified`, `SurfacePic`, and `SurfacePicHybrid` now all terminate in that same kernel path.
`LegacyBenchmark` remains separate as the reference reproduction / replay boundary.

## Runtime Layers

### `SurfaceRuntimePlan`
- compiles scenario-time inputs before runtime object construction
- currently freezes:
  - normalized `SurfaceChargingConfig`
  - route / strategy / legacy execution selection
  - solver policy flags
  - bridge-enable decisions
  - spectrum-to-plasma derived moments
- must not grow step-time mutable state, history containers, or exporter-owned side effects

### `SurfaceTransitionEngine`
- evaluates step-time transition decisions before entering the main advance pipeline
- currently resolves:
  - legacy replay short-circuit
  - reference-circuit mainline selection
  - shared surface runtime enablement
  - shared global coupled solve / live-PIC refresh / particle-transport coupling gates
  - solver-policy-derived iteration behavior for the coupled runtime
- should remain a decision layer, not become a second runtime state store

### `SurfaceScenarioOrchestrator`
- assembles the runtime route
- creates current, voltage, capacitance, reference, field, charge, and circuit models

### `SurfaceCurrentModel`
- computes local current decomposition
- owns benchmark or unified current logic
- now also exposes patch-level calibration and live PIC diagnostics
- exports `SurfaceKernelSnapshot` state so runtime routes can hand off distribution/current/DIDV
  as one aligned kernel input object

### `SurfaceCircuitModel`
- owns the body/patch/interface graph
- supports node and branch metadata
- supports implicit multi-node advancement
- consumes explicit kernel-input linearization rather than ad hoc local derivative scattering

### `SurfaceCapacitanceModel`
- owns capacitance-per-area logic
- can be combined with `BubbleCapacitanceEstimator`

### Providers
- `PotentialReferenceModel`
- `ElectricFieldProvider`
- `VolumeChargeProvider`
- `BubbleCapacitanceEstimator`
- `SurfaceReferenceStateProvider`

These providers make the current codebase field-aware without requiring a full volumetric Poisson solver.

## Current Graph Capabilities
- multiple bodies
- multiple patches
- body-body interfaces
- patch-body interfaces
- patch-patch interfaces
- node-level material ownership
- node-level current decomposition export
- node-level PIC calibration export
- interface current and conductance export
- graph-propagated reference-potential diagnostics
- interface-aware graph capacitance diagnostics

## Planned Higher-Order Direction
The current framework is still a reduced-order surface graph.

The next higher-order extension points are:

- interface-aware capacitance matrices
- graph-level reference propagation
- stronger neighbor-aware electric field estimates
- external field-solver bridge contracts
- future adapters for higher-order field solvers

The current codebase now includes lightweight default implementations for:

- `SurfaceReferenceGraphPropagator`
- `SurfaceGraphCapacitanceMatrixProvider`
- `SurfaceFieldSolverAdapter`

These are diagnostic and extension-oriented by design: they expose graph-level
state without forcing the runtime to become a full volumetric field solver.

The latest bridge boundary now also includes:

- `BodyBoundaryGroup`
- `PatchBoundaryGroup`
- `SurfaceBoundaryMapping`
- external request/result JSON sidecars
- adapter-fed `field_solver_capacitance_scale`

Those graph-level diagnostics are now exported through:

- node-level CSV series such as
  - `surface_node_<i>_propagated_reference_potential_v`
  - `surface_node_<i>_field_solver_reference_potential_v`
  - `surface_node_<i>_graph_capacitance_diagonal_f`
  - `surface_node_<i>_graph_capacitance_row_sum_f`
  - `surface_node_<i>_field_solver_coupling_gain`
- interface-level CSV series such as
  - `surface_interface_<i>_mutual_capacitance_f`
  - `surface_interface_<i>_voltage_drop_v`
  - `surface_interface_<i>_power_w`
- `.graph.txt` sidecar summaries for quick topology inspection
- `.monitor.json` unified monitor snapshot for graph/field/benchmark governance
- `.boundary_mapping.json` for surface-to-boundary ownership
- `.field_request.json` for external solver handoff
- `.field_result_template.json` for bridge result schema
- `.field_result.json` for populated field bridge outputs
- `.graph_matrix.json` for machine-readable matrix interchange
- `.field_adapter.json` for machine-readable adapter contract
- `.volume_stub.json` for future volume-state handoff
- `.volume_mesh_stub.json` for pseudo boundary-face and cell mesh handoff
  - includes graph-derived `cell_neighbors` and `face_neighbors`
  - includes `shared_face_area_m2` and `face_distance_m`
  - includes cell centers, surface normals, and face centers for richer geometry-aware exchange
- includes `cell_faces` for explicit cell-to-face ownership
- supports external cell volume and initial cell-state overrides
- `.field_bridge_manifest.json` for bridge artifact discovery
- `.surface_volume_projection.json` for explicit projection contracts
  - includes explicit `projection_weights`
- `.volumetric_adapter.json` for volumetric adapter capabilities
- `.volume_request.json` for external volume-solve requests
- `.volume_result_template.json` for external volume-solve result schema
- `.volume_result.json` for populated volume bridge outputs
- `.volume_history.json` for time-resolved volume-state history

External bridge result ingestion is now contract-aware:

- result payloads honor `schema_version` when present
- unsupported schema versions are rejected and runtime falls back to internal adapters
- field results can map by `node_id`, logical node id, or `boundary_group_id`
- volume results can map by `cell_id` or `boundary_group_id`

The default volumetric runtime now does more than export bridge artifacts:

- deposits surface charge into pseudo control volumes
- solves a dense discrete Poisson system on the control-volume graph
- projects volume potential and field back to node runtime state
- carries volume potential and reconstructed charge density forward across steps
- can prefer externally supplied mesh/projection geometry over internally generated pseudo mesh
- uses externally supplied cell centers, face centers, face normals, and neighbor distances
  directly inside the control-volume Poisson assembly when present
- performs a short self-consistent relaxation loop so reconstructed volume charge and
  volume potential feed the next internal Poisson iteration before surface writeback
- automatically switches between dense direct and iterative linear solves depending on
  volumetric case richness and cell count
- can be forced to `dense_only`, `iterative_only`, or left at `iterative_or_dense_auto`
- can tune iterative linear solve behavior through:
  - `volume_linear_max_iterations`
  - `volume_linear_tolerance_scale`
  - `volume_linear_relaxation`
- exposes explicit control parameters for that loop:
  - `field_volume_outer_iterations`
  - `field_volume_outer_tolerance`
  - `field_volume_outer_relaxation`
  - `volume_self_consistent_iterations`
  - `volume_self_consistent_tolerance_v`
  - `volume_charge_relaxation`
  - `volume_potential_relaxation`
- exports per-node solver diagnostics so convergence behavior is observable:
  - `volume_solver_mode_id`
  - `volume_solver_iterations`
  - `volume_solver_linear_iterations`
  - `volume_solver_converged`
  - `volume_solver_residual_norm`
  - `volume_solver_max_delta_v`
  - `volume_solver_matrix_nnz`
  - `volume_solver_cell_count`
  - `field_volume_coupling_iterations`
  - `field_volume_coupling_converged`
  - `field_volume_coupling_max_delta`
  - `field_volume_coupling_relaxation_used`
  - `external_volume_feedback_blend_factor`
  - `external_volume_feedback_mismatch_metric`
  - `external_volume_feedback_applied`

The iterative path is now sparse-like at the row level: off-diagonal couplings are compacted
before iteration, so iterative volumetric solves no longer pay a full dense scan cost for each
row update.

Above that, the runtime now also performs an outer field-volume coupling loop. Field-adapter
updates and volumetric-adapter updates are iterated together until the coupled state settles or
the configured outer iteration budget is exhausted.

The implicit circuit voltage runtime is now matrix-aware for graph mutual capacitance:

- branch-level mutual capacitance is converted to `C_mutual / dt` coupling terms
- those terms are injected into the implicit nodal matrix (diagonal + off-diagonal + rhs)
- this makes graph capacitance participate directly in voltage advancement instead of only
  scaling local lumped capacitance
- runtime weight is controlled by material scalar property:
  - `graph_matrix_runtime_weight` (default `0.35`)

The graph sidecar is intended to make higher-order graph behavior easy to audit.
It now reports:

- strongest-field node
- strongest-coupling interface
- maximum neighbor potential delta
- maximum neighbor field contrast
- maximum node reference offset
- maximum field-solver reference offset
- matrix family and solver-adapter hint

The monitor sidecar (`.monitor.json`) provides a stable machine-readable monitor
surface for tooling and dashboards:

- `graph`: topology/coupling and strongest-node/interface summaries
- `field`: reference-offset envelope and strongest-field diagnostics
- `benchmark`: route/mode/source plus consistency status and RMSE snapshot

## Benchmark Boundary
- `LegacyBenchmark` remains the reference reproduction route
- `SCDATUnified` remains the advanced framework route
- `SurfacePic` and `SurfacePicHybrid` now share the same benchmark/export contracts as the unified
  surface-kernel route
- benchmark replay and execute modes are intentionally separate from unified prediction

## Regression Contracts
The regression boundary is now explicit rather than implicit:

- `.benchmark_case.json`
  - machine-readable benchmark-case contract
  - carries reference datasets and validation metrics
- `.simulation_artifact.json`
  - machine-readable runtime artifact contract
  - carries field / particle / surface metric families
- `surface_reference_matrix.json`
  - declares required metadata, benchmark metrics, and artifact metrics for aligned cases
- `check_surface_reference_matrix_gate.py`
  - executes the formal GEO / LEO / MATLAB reference-family gate

The aligned benchmark metrics now cover:

- patch/body potential RMSE
- current-component RMSE families (`jnet`, `je`, `jse`, `jb`, `ji`, `jsi`, `jph`, `jcond`)
- segmented transient / midrise / plateau metrics where available
- negative-tail body-current metrics where available
- final DIDV export via `final_current_derivative_a_per_m2_per_v`
