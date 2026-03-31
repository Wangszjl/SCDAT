# `ref` Gap Matrix After Phase-1 Refactor

## Current baseline

- Active architecture now centers on `Tools/* + Toolkit/* + Main/*`.
- `Tools/Material`, `Tools/Diagnostics`, `Tools/Output`, and `Tools/Coupling` are already live.
- `Toolkit/Plasma Analysis`, `Toolkit/Surface Charging`, `Toolkit/Internal Charging`, and `Toolkit/Vacuum Arc` are all runnable libraries with stable `initialize / advance / getStatus / reset / exportResults` entry points.
- The current test baseline is `69/69` green in `build_codex`.
- `Main/main.cpp` is now a scenario dispatcher, while `Main/main_PICcore.cpp` remains the PIC-MCC benchmark entry.

## Remaining `ref` gaps by module

| `ref` module | current status | current target in this repo | remaining gap | priority |
|---|---|---|---|---|
| `material` | mostly implemented | `Tools/Material` | exact `SPISXMLParser` API is still missing; importer path is csv/text oriented | P1 |
| `diagnostics` | implemented | `Tools/Diagnostics` | no blocking gap for phase 1 | covered |
| `output` | implemented | `Tools/Output` | native binary HDF5 backend is still absent; current `.h5` is a lightweight archive format | P1 |
| `coupling` | partially implemented | `Tools/Coupling` + `Toolkit/Vacuum Arc/SurfaceArcCoupling.*` | exact `MultiPhysicsCoupling` compatibility layer is still missing | P1 |
| `FluidAlgorithms` | mostly implemented | `Toolkit/Plasma Analysis` | `MultiScaleSpatialDecomposer` and `MultiTimeScaleSynchronizer` are not present yet | P1 |
| `algorithms` | mostly implemented | `Toolkit/Surface Charging`, `Toolkit/Internal Charging`, `Toolkit/Vacuum Arc`, `Tools/FieldSolver` | `AdvancedNumericalMethods` and `TransportSolvers` compatibility wrappers are still missing | P2 |
| `field` | partially implemented | `Tools/FieldSolver` | exact `EnhancedElectricFieldSolver` API is still missing | P1 |
| `solver` | partially implemented | `Tools/Solver`, `Tools/FieldSolver`, `Tools/PICcore` | `HighOrderSolver`, `ParallelSolver`, `SubcyclingController`, `TimeStepController` are still missing as first-class modules | P1 |
| `solvers` | implemented | `Tools/FieldSolver` | no blocking gap for the current phase | covered |
| `mesh` | partially implemented | `Tools/Mesh` | `AdaptiveMesh`, `MeshLoader`, and `SpatialIndex` are still missing as named modules | P1 |
| `particle` | partially implemented | `Tools/Particle`, `Tools/Boundary`, `Tools/PICcore` | `ParticleSoA`, `ParallelParticlePusher`, and architecture-specific AVX/CUDA pushers are still absent | P2 |
| `collision` | covered | `Tools/Interactions/Collisions` | no phase-1 blocking gap | covered |
| `benchmarks` | partially implemented | `Main/main_PICcore.cpp` | `GPUPerformanceBenchmark` is still absent; `TurnerBenchmark` is represented by the current PIC-MCC benchmark driver rather than a direct class port | P2 |
| `ArcPICAlgorithms` | mostly implemented | `Toolkit/Vacuum Arc` | current implementation is CPU-only and does not mirror every `ArcPICTypes` helper from `ref` | P1 |
| `Geant4Algorithms` | not implemented | deferred | external dependency path intentionally postponed | P2 |
| `SPISAlgorithms` | not implemented | deferred | `SpacecraftAnalysisTools` style high-level orchestration still missing | P2 |
| `core` | lightly covered | `Tools/Basic` and module-local interfaces | `AlgorithmSwitchingManager`, `ModuleRegistry`, `InterfaceStandardization`, `PerformanceMonitor`, `SystemIntegrationOptimizer`, `MultiScaleCouplingValidator` are still missing | P2 |
| `config` | not implemented | deferred | `ConfigurationManager`, `ConfigurationUI`, `ParameterValidator` still missing | P2 |
| `io` | not implemented | deferred | `HighPerformanceIO` not present | P2 |
| `gpu` | not implemented | deferred | CUDA / multi-GPU path intentionally postponed | P3 |
| `performance` | not implemented | deferred | thread-pool and optimization layer from `ref` not migrated | P3 |
| `python` | not implemented | deferred | no embedded Python or binding layer | P3 |
| `ui` | not implemented | deferred | intentionally postponed | P3 |
| `visualization` | not implemented | deferred | replaced for now by csv/vtk export plus Python plotting scripts | P3 |

## Most important next ports

### P1: directly useful to the current architecture

1. `Tools/Mesh`
   Add `AdaptiveMesh`, `MeshLoader`, and `SpatialIndex` in the current mesh module instead of rebuilding a parallel mesh stack.
2. `Tools/FieldSolver`
   Add an `EnhancedElectricFieldSolver` compatibility layer that wraps the current Poisson and multiboundary solvers.
3. `Tools/Coupling`
   Add a `MultiPhysicsCoupling` facade on top of `CouplingBase` and `MultiPhysicsManager`.
4. `Toolkit/Plasma Analysis`
   Add `MultiScaleSpatialDecomposer` and `MultiTimeScaleSynchronizer` so the fluid path can genuinely coordinate different resolution and cadence regions.
5. `Tools/Output`
   Replace the lightweight `.h5` archive with a real HDF5 backend only if external dependencies are acceptable.

### P2: useful, but not blocking the current runtime chain

1. `Tools/Solver`
   Port `HighOrderSolver`, `ParallelSolver`, `TimeStepController`, and `SubcyclingController`.
2. `Toolkit/Internal Charging`
   Add `Geant4` transport adapters as optional plugins, not as core build dependencies.
3. `Tools/Material`
   Add `SPISXMLParser` if you want direct import parity with older SPIS-style material datasets.
4. `Main`
   Add a higher-level `SpacecraftAnalysisTools` dispatcher only after the toolkit cases stabilize.

### P3: postpone until physics and workflow settle

1. `gpu`
2. `performance`
3. `python`
4. `ui`
5. `visualization`

## Placement rule for all further ports

- Reusable numerics, data structures, import/export, and diagnostics go to `Tools/*`.
- Plasma, surface charging, internal charging, and vacuum arc scenario logic stay under their respective `Toolkit/*` folders.
- New ports must extend the current active modules instead of reviving `ref` as a parallel build tree.

## Related docs

- `documentation/ref_audit_2026-03-30.md`
- `documentation/toolkit_case_presets_2026-03-30.md`
