# `Top/Default` class-to-runtime 映射

## 1. 目的

本文档用于落实 `Round 48 / Wave AC`，把 `ref/SPIS/num-master/src/main/java/spis/Top/Default/**` 中的 class 与当前项目中的 runtime / scenario / config / metadata 承载落点显式对应起来。

本文档回答的问题不是“`top_default` family 是否已经封板”，而是：

- `Top/Default` 下每个 SPIS class 现在由哪个对象层或配置层承载；
- 这些承载是“单类等价”还是“聚合承载”；
- 当前有哪些类仍然只达到“文档化近似实现”。

> 口径说明：`top_default` family 已在 `Round 38` family 扩边后进入 `completed`，并已通过 `surface_reference_matrix gate` 与 sealoff。  
> 本文档只处理 class-level 显式映射，不推翻现有 family sealoff。

## 2. 当前承载骨架

`Top/Default` 在当前项目中主要由以下 4 层共同承载：

| 当前承载层 | 角色 | 主要文件 |
| --- | --- | --- |
| preset / catalog 层 | 管理默认 preset、别名、说明文本、主线入口 | `Toolkit/Surface Charging/include/SurfaceScenarioCatalog.h`；`Toolkit/Surface Charging/src/SurfaceScenarioCatalog.cpp`；`Toolkit/Surface Charging/include/SurfaceChargingCases.h` |
| config / loader 层 | 将 JSON / CLI 输入装配为 `SurfaceChargingScenarioPreset` 与 `SurfaceChargingConfig` | `Main/SurfaceScenarioLoader.cpp` |
| runtime plan 层 | 把 preset/config 编译为运行时可执行计划、求解策略、步进参数 | `Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp` |
| runtime / metadata 层 | 初始化 surface kernel、导出 `surface_top_default_*` metadata、保存 source/entrypoint/solver/observer 结果 | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |

辅助证据主要来自：

- `Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`
  - `SurfaceScenarioCatalogTest.ListsPresetsAndResolvesAliases`
  - `SurfaceRuntimePlanCompilerTest.CompilesPresetMetadataAndSolverPolicy`
- `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`
  - `SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsTopTransitionFamilyMetadata`
  - `SurfaceChargingSmokeTest.PresetCatalogRouteExportsTopTopFamilyMetadata`

## 3. 映射规则

本文档对每个 SPIS class 使用以下结论词：

- `文档化近似实现`
  - 当前已能指出明确承载层和证据；
  - 但未在 C++ 中保留与 SPIS Java 对象同构的独立类型。
- `聚合承载`
  - 该 class 的职责被多个 runtime/config 对象共同承担。

对 `Top/Default` 来说，当前 25 个 class 均属于“文档化近似实现”或“聚合承载”，没有要求补成 25 个独立 C++ 类型。

## 4. 逐 class 映射

| SPIS class | 当前承载落点 | 当前承载类型 | 主要证据 | 说明 |
| --- | --- | --- | --- | --- |
| `ActiveZone` | `SurfaceChargingConfig.boundary_mappings`、structured-topology / body-patch 分组、runtime 中的边界激活路径 | 文档化近似实现 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`（boundary mapping smoke） | 当前没有单独 `ActiveZone` 类型，活跃区域由边界映射和 patch/body 选择承载。 |
| `Common` | `SurfaceScenarioCatalog::makeDefaultPreset()` 与 `SurfaceChargingScenarioPreset` 默认字段集合 | 聚合承载 | `Toolkit/Surface Charging/src/SurfaceScenarioCatalog.cpp`；`Toolkit/Surface Charging/include/SurfaceChargingCases.h` | 通用默认值被折叠到 preset catalog 和 `SurfaceChargingCases`。 |
| `Credits` | `SurfaceChargingScenarioPreset.description`、`SurfaceRuntimePlan.source_description`、导出 metadata 中的 source description | 文档化近似实现 | `Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | SPIS 中的说明性对象被当前 preset 描述文本与 sidecar metadata 吸收。 |
| `Distribution` | `SurfaceChargingConfig.distribution_model`、谱参数、`SurfaceDistributionFunction` 路由 | 聚合承载 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Tools/Particle/include/SurfaceDistributionFunction.h` | 不保留 `Distribution` 壳层对象，直接落到 distribution 配置和函数族。 |
| `Duplicable` | `SurfaceChargingScenarioPreset` 的值语义复制、catalog alias 解析与 preset 派生 | 文档化近似实现 | `Toolkit/Surface Charging/src/SurfaceScenarioCatalog.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp` | “可复制对象”语义由 preset struct 与 alias 解析承担。 |
| `Field` | `SurfaceChargingConfig` 中 field/solver 相关参数、native-volume / bridge route、runtime metadata | 聚合承载 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 当前没有单独 `Field` 壳层对象，field 语义被 solver/bridge/runtime route 直接消费。 |
| `GlobalParameter` | `SurfaceRuntimePlan` 中的 route、strategy、scheduled steps、时间步长与 source 元信息 | 文档化近似实现 | `Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`SurfaceRuntimePlanCompilerTest.CompilesPresetMetadataAndSolverPolicy` | 全局参数被收口到 runtime plan。 |
| `InstrumentParameter` | `SurfaceChargingConfig.surface_instrument_set_kind`、observer/instrument metadata、loader 中的 `surface_instrument_set` 解析 | 文档化近似实现 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 仪器参数未作为单独对象存在，而是作为 instrument set 与 observer route 的配置字段。 |
| `Interactor` | `config.material` 中的 interactor family 选择与 DensePlasmaSurfaceCharging 中的 interactor route | 聚合承载 | `Main/SurfaceScenarioLoader.cpp`；`Tools/Material/include/SurfaceInteraction.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `Top/Default/Interactor` 在当前项目中表现为对 `SurfInteract` family 的顶层配置引用。 |
| `Listener` | 运行时观察与回调需求由 `SurfaceSimulationRunner`、transition observer、smoke/object-layer 测试吸收 | 文档化近似实现 | `Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 没有单独 listener 接口；监听语义已并入 runner/observer 路径。 |
| `LocalParameter` | `MaterialProperty` 局部标量、patch/body 局部属性、material scalar property bags | 聚合承载 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 局部参数通过 material/body_material 和 patch/boundary 配置携带。 |
| `MaterialProperty` | `Material::MaterialProperty`、`config.material`、`config.body_material` | 聚合承载 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 这是 `Top/Default` 中与当前 C++ 命名最接近的一层，但仍不是 SPIS `Top/Default` 原类的直接镜像。 |
| `Mesh` | structured topology、boundary mapping、grid/runtime 装配与 volume bridge | 聚合承载 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | mesh 组织层在当前项目中由 loader + runtime glue 承载。 |
| `NamedObject` | `SurfaceChargingScenarioPreset.name`、`SurfaceRuntimePlan.source_name`、`surface_top_top_source_name` metadata | 文档化近似实现 | `Toolkit/Surface Charging/include/SurfaceChargingCases.h`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 名称语义已稳定存在，但作为 preset/runtime metadata 字段存在。 |
| `ObjectWithBoundaries` | `boundary_mappings`、body/patch/interface 分组、surface/grid metadata | 聚合承载 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | “带边界对象”被当前 body/patch/interface 组织与 boundary mapping 吸收。 |
| `ObjectWithTemperature` | `photoelectron_temperature_ev`、`body_photoelectron_temperature_ev`、plasma 温度与 thermal sidecar | 聚合承载 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 温度属性由 config 字段和 material scalar property 承载。 |
| `Parameter` | `SurfaceChargingConfig` 顶层配置字段与 JSON/CLI loader 键解析 | 文档化近似实现 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h` | 参数对象层在当前项目中被具体配置字段取代。 |
| `PlasmaIO` | JSON loader、CSV/metadata sidecar 导出、reference matrix case I/O | 聚合承载 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 输入输出职责拆分到 loader 与 exporter/metadata。 |
| `ScalableObject` | material scalar properties、solver scaling、barrier/current scaling、runtime metadata 中的 family/scalar 导出 | 聚合承载 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Tools/FieldSolver/include/SurfaceBarrierModels.h` | “可缩放对象”没有统一壳层，而是分散到 material/solver/barrier 标量属性。 |
| `SimulationObject` | `SurfaceRuntimePlan`、`DensePlasmaSurfaceCharging::initialize/advance/exportResults`、`SurfaceSimulationRunner::run` | 聚合承载 | `Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp` | 模拟对象被统一收口为 runtime plan + charging engine。 |
| `Solver` | `solver_config`、`SolverPolicyFlags`、`SurfaceSolverFacade` 路由 | 聚合承载 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Tools/Solver/include/SurfaceSolverFacade.h` | 当前的 `Solver` 壳层直接映射到 solver policy / facade route。 |
| `SpisDefaultPartTypes` | 默认粒子/发射物种选择、surface emission / distribution 默认路由 | 文档化近似实现 | `Toolkit/Surface Charging/include/SurfaceChargingCases.h`；`Tools/Particle/include/SurfaceDistributionFunction.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 默认粒子类型不再作为单独 Top/Default 对象，而是作为 preset+distribution 约定存在。 |
| `SpisDefaultSampling` | `sampling_policy`、deterministic/default sampling 路由、sampler metadata | 文档化近似实现 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 默认采样策略被当前 config 字段和 sampler metadata 承载。 |
| `ThermalElement` | 温度相关 config、photoelectron/body thermal 参数、thermal current 相关 runtime 输出 | 聚合承载 | `Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 当前 thermal 语义没有独立 Top 层对象，而是嵌入材料/发射/温度字段。 |
| `ThreadManager` | `SurfaceSimulationRunner` 的执行协调、adaptive stepping 循环与 step budget 管理 | 文档化近似实现 | `Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`；`SurfaceSimulationRunnerTest.RunsDefaultPresetWithZeroSteps` | 当前没有线程管理器类；执行协调被 runner 顺序循环吸收。 |

## 5. 汇总结论

### 5.1 结论

- `Top/Default` 当前 **没有功能缺口型 `尚无实现` 项**。
- 当前 25 个 SPIS class 都已能回链到明确的当前承载层，因此 `Round 48` 之后它们不再只是“审计表里的一句聚合描述”，而是具备文档化的逐类解释。
- 这些 class 仍然大多属于：
  - `文档化近似实现`
  - 或 `聚合承载`
- 这是因为当前项目把 `Top/Default` 作为“配置/对象壳层”而不是“运行算法壳层”处理，天然更适合压缩进 preset/config/runtime plan，而不是保留 1:1 Java 类树。

### 5.2 对主计划的影响

- `top_default` family 的 sealoff 结论保持不变。
- `Round 48` 的完成，意味着 `Top/Default` 已从“class-level 近似实现但未文档化”升级为“class-level 近似实现且已文档化、可审计、可回链证据”。
- 下一步重点转向：
  - `Top/Top + Top/Simulation` shell-to-runtime 对照
  - `Util/Func` class-to-helper / formula 映射
  - `Util/Instrument` catalog / observer-route 映射

## 6. 证据索引

- 代码
  - `Toolkit/Surface Charging/include/SurfaceChargingCases.h`
  - `Toolkit/Surface Charging/include/SurfaceScenarioCatalog.h`
  - `Toolkit/Surface Charging/src/SurfaceScenarioCatalog.cpp`
  - `Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`
  - `Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`
  - `Toolkit/Surface Charging/include/SurfaceSimulationRunner.h`
  - `Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`
  - `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`
  - `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`
  - `Main/SurfaceScenarioLoader.cpp`
- 测试
  - `Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`
  - `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`
- 审计基线
  - `Toolkit/Surface Charging/docs/surface_spis_class_audit.json`
  - `Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md`
  - `Toolkit/Surface Charging/docs/surface_num_alignment_refactor_plan.md`
