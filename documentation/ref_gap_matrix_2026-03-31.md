# `ref` Content Gap Matrix After Surface-Charging Multiscale Update

## Current baseline

- Active runtime architecture is now `Tools/* + Toolkit/* + Main/*`; `ref` is no longer part of the build.
- The current executable chain is complete for:
  - geometry / mesh core
  - particle / collision / boundary
  - field solver / PIC core
  - material / diagnostics / output / coupling base
  - plasma analysis / surface charging / internal charging / vacuum arc toolkit scenarios
- `Main/main.cpp` is the scenario dispatcher; `Main/main_PICcore.cpp` remains the PIC-MCC regression and benchmark entry.
- Current verified baseline is `78/78` passing tests in `build_codex`.

## What is already covered, even if the names differ from `ref`

| `ref` area | current landing place | status |
|---|---|---|
| geometry | `Tools/Geometry` | covered |
| collision / MCC | `Tools/Interactions/Collisions` | covered |
| Poisson / multiboundary field solve | `Tools/FieldSolver` | covered |
| PIC cycle / particle push / runtime boundary handling | `Tools/PICcore` + `Tools/Boundary` | covered |
| materials, loaders, surface interaction | `Tools/Material` | mostly covered |
| diagnostics | `Tools/Diagnostics` | covered |
| csv / vtk / lightweight h5 export | `Tools/Output` | covered |
| toolkit-level plasma / surface / internal / arc cases | `Toolkit/*` | covered |

## Missing content that is still real, not just a rename

### P1: the main functional gaps relative to `ref`

| `ref` module | missing or partial content in current repo | note |
|---|---|---|
| `mesh` | `AdaptiveMesh`, `MeshLoader`, `SimpleMesh`, `SpatialIndex`, `MeshTypes` compatibility layer | current mesh core runs, but these named extensions are still absent |
| `field` | `EnhancedElectricFieldSolver` | current Poisson and multiboundary solvers are usable, but there is no dedicated compatibility facade |
| `coupling` | `MultiPhysicsCoupling` | current `Tools/Coupling` has base manager pieces, not the exact high-level facade |
| `FluidAlgorithms` | `MultiScaleSpatialDecomposer` | current plasma-analysis toolkit runs without this spatial decomposition layer |
| `material` | `UnifiedMaterialDatabase`, `SPISXMLParser`, `SemiconductorMaterial`, `MaterialProperties` compatibility layer | current import path is csv/text oriented and focused on active cases |
| `solver` | `HighOrderSolver`, `PICPusher`, `SolverTypes` | some adjacent functionality exists, but not these first-class compatibility modules |
| `core` | `AlgorithmSwitchingManager`, `ModuleRegistry`, `MultiScaleCouplingValidator`, `PerformanceMonitor`, `SystemIntegrationOptimizer` | orchestration and integration meta-layer still thin |
| `particle` | `AdvancedBoundaryConditions`, `CollisionModel`, `ModernParticle`, `ParallelParticlePusher` | current particle stack is sufficient for the active CPU chain, not parity-complete |

### P2: useful but currently deferred

| `ref` module | missing content | note |
|---|---|---|
| `config` | `ConfigurationManager`, `ConfigurationUI`, `ParameterValidator` | no centralized configuration subsystem yet |
| `Geant4Algorithms` | `Geant4Integrator`, `Geant4ParticleTransport`, `MonteCarloCollisionModel` | intentionally deferred to avoid external dependency pressure |
| `SPISAlgorithms` | `SpacecraftAnalysisTools` | main dispatcher exists, but not the full high-level SPIS-style orchestration layer |
| `benchmarks` | `GPUPerformanceBenchmark` | current benchmark focus is CPU PIC-MCC regression |
| `collision` | direct parity wrappers such as `CoulombCollision`, `MonteCarloCollision` | active collision path exists under the reorganized interaction module |

### P3: intentionally postponed subsystems

| `ref` module | current status |
|---|---|
| `gpu` | not ported |
| `performance` | not ported |
| `python` | plotting scripts only, no embedded binding layer |
| `ui` / `ui_demo` | not ported |
| `visualization` | replaced for now by csv/vtk export and Python plotting |
| `io` | no `HighPerformanceIO` parity layer |

## What is complete enough to call the project structurally whole

- The active PIC-MCC runtime is no longer blocked by placeholder solver or boundary code.
- The toolkit layer is no longer empty shell code: all four domains are buildable and runnable.
- Surface charging now has distinct operating regimes:
  - `LEO` flowing-plasma branch
  - `GEO` short-window PIC calibration + long-horizon kinetic charging branch
  - `thruster plume` directed-flow branch
- Result export, summary, plotting, and regression testing are all live and connected.

## Recommended next ports if the goal is fuller `ref` parity

1. `Tools/Mesh`: add `AdaptiveMesh`, `MeshLoader`, and `SpatialIndex`.
2. `Tools/FieldSolver`: add `EnhancedElectricFieldSolver` as a facade over current solvers.
3. `Tools/Coupling`: add `MultiPhysicsCoupling` on top of the current manager base.
4. `Toolkit/Plasma Analysis`: add `MultiScaleSpatialDecomposer`.
5. `Tools/Material`: add `SPISXMLParser` and a `UnifiedMaterialDatabase` compatibility layer.

## Placement rule that still holds

- Reusable numerics, solvers, material handling, diagnostics, and exporters belong in `Tools/*`.
- Plasma simulation, surface charging, internal charging, and discharge / arc workflows belong in their corresponding `Toolkit/*` directories.
- Further `ref` ports should continue extending the active modules rather than recreating `ref` as a parallel source tree.
