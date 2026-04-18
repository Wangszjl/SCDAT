# `Top/Top + Top/Simulation` shell-to-runtime 对照

## 1. 目的

本文档用于落实 `Round 49 / Wave AD`，把 `ref/SPIS/num-master/src/main/java/spis/Top/Top/**` 与 `ref/SPIS/num-master/src/main/java/spis/Top/Simulation/**` 中的壳层 / 入口层 class，与当前项目中的 `SurfaceScenarioCatalog`、`SurfaceScenarioLoader`、`SurfaceSimulationRunner`、`SurfaceRuntimePlan`、`DensePlasmaSurfaceCharging` 之间的吸收关系显式写清。

本文档回答的问题是：

- `Top/Top` 里的入口壳层 class 在当前项目里由谁承接；
- `Top/Simulation` 里的 simulation 壳层 class 在当前项目里由谁承接；
- 哪些类属于 `聚合等价实现`，哪些仍只达到 `文档化近似实现`。

> 口径说明：`top_top` 与 `top_simulation` family 已封板并通过 sealoff。  
> 本文档只处理 class-level 壳层吸收关系，不改变当前 family `completed` 结论。

## 2. 当前入口骨架

当前项目的 `Top` 壳层入口主要收口到以下路径：

| 当前入口层 | 角色 | 主要文件 |
| --- | --- | --- |
| preset catalog | 列举可运行 scenario、解析别名、生成默认 preset | `Toolkit/Surface Charging/include/SurfaceScenarioCatalog.h`；`Toolkit/Surface Charging/src/SurfaceScenarioCatalog.cpp` |
| config loader | 把 `surface-config` JSON 解析为 `SurfaceChargingScenarioPreset` | `Main/SurfaceScenarioLoader.cpp` |
| CLI shell | 暴露 `surface` / `surface-config` 命令，并把 preset / JSON 输入统一导入运行主线 | `Main/main.cpp` |
| simulation runner | 用统一 runner 承接运行生命周期、步进、导出 | `Toolkit/Surface Charging/include/SurfaceSimulationRunner.h`；`Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp` |
| runtime plan | 将 preset/config 编译为运行时 route、solver policy、time-step 计划 | `Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp` |
| charging engine | 真正初始化、推进、导出 metadata 与 sidecar | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |

统一执行链是：

`CLI / JSON -> SurfaceScenarioCatalog or SurfaceScenarioLoader -> SurfaceSimulationRunner -> compileSurfaceRuntimePlan -> DensePlasmaSurfaceCharging`

## 3. 逐 class 对照

### 3.1 `Top/Top`

| SPIS class | 当前承载落点 | 当前承载类型 | 主要证据 | 说明 |
| --- | --- | --- | --- | --- |
| `Scenario` | `SurfaceScenarioCatalog` 的 preset 目录、`SurfaceChargingScenarioPreset`、`surface` CLI 命令 | 聚合等价实现 | `Toolkit/Surface Charging/include/SurfaceScenarioCatalog.h`；`Toolkit/Surface Charging/src/SurfaceScenarioCatalog.cpp`；`Main/main.cpp` | SPIS 中的 scenario 壳层，在当前项目中由 preset catalog + preset struct 统一承接。 |
| `PotentialSweep` | `SurfaceChargingScenarioPreset` 中的步数 / 时间步长 / duration 计划，以及 `SurfaceSimulationRunner` 的步进主线 | 聚合等价实现 | `Toolkit/Surface Charging/include/SurfaceChargingCases.h`；`Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp` | 当前没有独立 `PotentialSweep` 类型，其职责已并入 preset 执行计划与 adaptive stepping。 |
| `UIInvokable` | `Main/main.cpp` 中的 `surface` / `surface-config` CLI 壳层入口 | 文档化近似实现 | `Main/main.cpp`；`Main/SurfaceScenarioLoader.cpp` | SPIS 的 UI 可调用入口在当前项目中被 CLI / JSON 命令接口吸收。 |
| `SpisTopMenu` | `main()` 中的命令分发、`presets` 列表、模块命令入口 | 文档化近似实现 | `Main/main.cpp` | 顶层菜单在当前项目中由命令行分发承担，不存在独立菜单对象。 |
| `NumTopFromUI` | `surface-config` 路径中 `SurfaceScenarioLoader::loadFromJson(...)` + `runSurface(...)` | 聚合等价实现 | `Main/SurfaceScenarioLoader.cpp`；`Main/main.cpp`；`Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp` | “来自 UI/配置的数值入口”被 JSON loader + runner 组合承接。 |

### 3.2 `Top/Simulation`

| SPIS class | 当前承载落点 | 当前承载类型 | 主要证据 | 说明 |
| --- | --- | --- | --- | --- |
| `Simulation` | `SurfaceSimulationRunner::run(...)` + `DensePlasmaSurfaceCharging::initialize/advance/exportResults` | 已等价实现 | `Toolkit/Surface Charging/include/SurfaceSimulationRunner.h`；`Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp` | 当前已经存在明确的 simulation 主执行对象。 |
| `PlasmaScSimulation` | `SurfaceSimulationRunner` 驱动的 surface charging 主线 | 已等价实现 | `Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | plasma + spacecraft 联合模拟在当前项目中由统一 surface charging runner 承载。 |
| `SimulationFromUIParams` | `SurfaceScenarioLoader::loadFromJson(...)` + `compileSurfaceRuntimePlan(...)` | 已等价实现 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp` | UI 参数构造 simulation 的职责已被 JSON loader + runtime plan 编译链吸收。 |
| `SimulationListener` | runner 的错误回传、transition observer、metadata / smoke 可观测路径 | 文档化近似实现 | `Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 当前没有单独 listener interface，监听/观察职责被 runner 结果对象与 observer/metadata 路径吸收。 |
| `BundleAPI_SimulationInit` | `compileSurfaceRuntimePlan(...)` + `DensePlasmaSurfaceCharging::initialize(...)` 初始化链 | 聚合等价实现 | `Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 初始化 bundle API 的职责被 runtime plan 编译与 charging engine 初始化共同承担。 |

## 4. 证据交叉

### 4.1 代码入口

- `Main/main.cpp`
  - `surface` 命令走 preset catalog
  - `surface-config` 命令走 JSON loader
  - 最终统一进入 `runSurface(...)`
- `Main/SurfaceScenarioLoader.cpp`
  - JSON -> `SurfaceChargingScenarioPreset`
  - 从 `base_preset`、`run`、`config` 等字段构造 scenario
- `Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`
  - 调用 `compileSurfaceRuntimePlan(...)`
  - 驱动 `DensePlasmaSurfaceCharging::initialize/advance/exportResults`
- `Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`
  - 把 preset/config 编译为 route、solver policy、time-step 计划
- `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`
  - 导出 `surface_top_top_*`、`surface_top_simulation_*` metadata

### 4.2 测试与 smoke

- `Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`
  - `SurfaceScenarioCatalogTest.ListsPresetsAndResolvesAliases`
  - `SurfaceSimulationRunnerTest.RunsDefaultPresetWithZeroSteps`
  - `SurfaceRuntimePlanCompilerTest.CompilesPresetMetadataAndSolverPolicy`
- `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`
  - `SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsTopTransitionFamilyMetadata`
  - `SurfaceChargingSmokeTest.PresetCatalogRouteExportsTopTopFamilyMetadata`

其中：

- `PresetCatalogRouteExportsTopTopFamilyMetadata` 直接证明 `Top/Top` family signature 在 preset catalog 路径上可观测；
- `ScdatUnifiedRouteExportsTopTransitionFamilyMetadata` 同时给出 `surface_top_simulation_*` metadata，可作为 `Top/Simulation` family 的 smoke 证据；
- object-layer tests 证明 catalog、runner、runtime plan 三层入口实际可执行。

## 5. 汇总结论

### 5.1 `Top/Top`

- `Scenario`、`PotentialSweep`、`NumTopFromUI` 当前已达到 `聚合等价实现`。
- `UIInvokable`、`SpisTopMenu` 当前仍属于 `文档化近似实现`。
- 其根本原因不是功能缺失，而是当前项目选择用 CLI / JSON / preset catalog 替代 SPIS UI 壳层对象。

### 5.2 `Top/Simulation`

- `Simulation`、`PlasmaScSimulation`、`SimulationFromUIParams` 当前已有明确执行落点，可视为 `已等价实现`。
- `BundleAPI_SimulationInit` 当前属于 `聚合等价实现`，由 runtime plan 编译 + charging engine 初始化共同承担。
- `SimulationListener` 当前仍属于 `文档化近似实现`，因为监听语义被 runner 结果对象与 observer/metadata 路径吸收，未保留独立 listener interface。

### 5.3 对主计划的影响

- `Round 49` 完成后，`Top/Top + Top/Simulation` 已从“class audit 摘要中的壳层弱映射”升级为“有专门文档、能回链到代码与测试的 shell-to-runtime 解释”。
- 当前 class-level 主线将转入：
  - `Round 50 / Wave AE`：`Util/Func`
  - `Round 51 / Wave AF`：`Util/Instrument`

## 6. 证据索引

- 代码
  - `Main/main.cpp`
  - `Main/SurfaceScenarioLoader.cpp`
  - `Toolkit/Surface Charging/include/SurfaceScenarioCatalog.h`
  - `Toolkit/Surface Charging/src/SurfaceScenarioCatalog.cpp`
  - `Toolkit/Surface Charging/include/SurfaceSimulationRunner.h`
  - `Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`
  - `Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`
  - `Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`
  - `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`
  - `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`
- 测试
  - `Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`
  - `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`
- 审计基线
  - `Toolkit/Surface Charging/docs/surface_spis_class_audit.json`
  - `Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md`
  - `Toolkit/Surface Charging/docs/surface_num_alignment_refactor_plan.md`
