# Surface Charging vs SPIS 数值内核差距分析与 Surface 主线改进计划

## 1. 目标与使用规则

### 1.1 文档定位

- 本文档是 `Toolkit/Surface Charging` 对齐 `ref/SPIS/num-master` 数值内核的唯一长期维护主计划。
- 本文档既是差距分析文档，也是执行看板；每一轮推进结束后，必须直接在本文档中更新状态。
- 本文档当前只服务 `Surface Charging` 主线，不是整个仓库的总计划，也不负责管理其他 toolkit 的独立路线。
- 本文档对比范围限定为 Surface 主线直接依赖的模块：
  - `Toolkit/Surface Charging`
  - `Tools/Material`
  - `Tools/Particle`
  - `Tools/FieldSolver`
  - `Tools/Coupling`
  - `Tools/Solver`
  - `Tools/PICcore`
  - 与上述模块直接关联的 `scripts/`、`documentation/contracts/`
- `Toolkit/Internal Charging`、`Toolkit/Radiation`、`Toolkit/Vacuum Arc` 等模块只在需要说明复用边界时被引用，不进入本文档正文主表。

### 1.2 权威来源

- 架构与职责对标源：
  - `ref/SPIS/num-master/src/main/java/spis/Surf/**`
  - `ref/SPIS/num-master/src/main/java/spis/Circ/**`
  - `ref/SPIS/num-master/src/main/java/spis/Solver/**`
  - `ref/SPIS/num-master/src/main/java/spis/Vol/**`
  - `ref/SPIS/num-master/src/main/java/spis/Util/**`
  - `ref/SPIS/num-master/src/main/java/spis/Top/**`
- 回归与验收源：
  - `ref/GEO/**`
  - `ref/LEO/**`
  - `ref/mat/**`
- 当前项目事实来源：
  - `Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`
  - `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`
  - `Toolkit/Surface Charging/include/SurfaceScenarioCatalog.h`
  - `Toolkit/Surface Charging/include/SurfaceSimulationRunner.h`
  - `Tools/Material/include/SurfaceMaterialModel.h`
  - `Tools/Particle/include/SurfaceDistributionFunction.h`
  - `Tools/FieldSolver/include/SurfaceBarrierModels.h`
  - `Tools/Coupling/include/SurfaceCircuitCoupling.h`
  - `Tools/Coupling/include/SurfaceCurrentLinearizationKernel.h`
  - `Tools/Solver/include/SurfaceSolverFacade.h`
  - `Tools/PICcore/include/PICCycle.h`
  - `scripts/run/surface_reference_matrix.json`
  - `documentation/contracts/**`

### 1.3 使用规则

- 每一轮执行结束后，必须更新：
  - `状态`
  - `最后更新`
  - `证据`
  - `下一步`
- 任何任务不得只写大方向，必须能在一轮执行后明确改一次状态。
- 任何被冻结、受阻或转入历史归档的任务都必须写明原因。
- 本文档优先描述职责边界、模块落点、验收方式与兼容要求，不做逐文件 1:1 迁移清单。

### 1.4 外部兼容边界

- 必须保持以下外部接口和契约可持续兼容：
  - `SurfaceChargingConfig`
  - `SurfaceChargingScenarioPreset`
  - `SurfaceRuntimeRoute`
  - `SurfacePicStrategy`
  - `LegacyBenchmarkExecutionMode`
  - 现有主 CSV 时序输出
  - 现有 benchmark / graph / monitor / field / volume / mapping / bridge sidecar 家族
- 允许扩展元数据和 sidecar 字段。
- 不允许删除已有消费者依赖的主字段或 route/preset 名称。

## 2. Surface 主线六域差距分析

### 2.1 总览

| 域 | SPIS 主要落点 | 当前项目主要落点 | 现状判断 | 当前关键问题 |
| --- | --- | --- | --- | --- |
| `Surf` | `spis/Surf/*` | `Tools/Material`、`Tools/Particle`、`Tools/FieldSolver`、`Tools/Coupling`、`Toolkit/Surface Charging` | Surface 主线核心物理能力已基本收口到 `Tools/*` | 模型家族深度仍明显少于 SPIS Surf（侵蚀材料、丰富分布族、多 scaler 家族） |
| `Circ` | `spis/Circ/*` | `Tools/Coupling`、`Toolkit/Surface Charging` | Surface 节点网络、DIDV、动态激励入口已具备 | `R11-D` 已补齐 `PWL/EXP` 与 weighted/multi-surface DIDV；剩余差距集中在更高阶组合类（如 reducable/多层包装） |
| `Solver` | `spis/Solver/*` | `Tools/Solver`、`Tools/FieldSolver`、`Tools/Particle`、`Tools/PICcore` | facade 与求解策略路由已建立 | `initialize/advance` 仍承担较重 runtime 编排 |
| `Vol` | `spis/Vol/*` | `Tools/Boundary`、`Tools/Mesh`、`Tools/FieldSolver`、`Tools/Particle`、`Toolkit/Surface Charging` bridge/runtime | surface bridge 契约与 pseudo-volume 路径已系统化 | 当前仍不是 native volume kernel，对标口径需要纠偏 |
| `Util` | `spis/Util/*` | `Tools/Basic`、`Tools/Geometry`、`Tools/Diagnostics`、`Tools/Output`、`Tools/Particle`、`Tools/Solver` | Surface 主线依赖的 helper 已完成首轮归位 | 仍需防止新 helper 回流到 toolkit 巨型文件 |
| `Top` | `spis/Top/*` | `SurfaceScenarioCatalog`、`SurfaceScenarioLoader`、`SurfaceSimulationRunner`、`SurfaceRuntimePlan`、`SurfaceTransitionEngine`、`DensePlasmaSurfaceCharging`、`Main/*` | 顶层对象层第二阶段（plan + transition）已落地 | `R11-E` 已接入 `LocalTime/Spinning/SunFlux` 最小事件链；剩余差距集中在 `ConductivityEvolution` 等扩展事件族 |

### 2.2 Surf 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | `Material`、`SurfInteract`、`SurfField`、`SurfDistrib` 与 `Util/DistribFunc` 共同构成表面材料响应、入射/发射分布、barrier/recollection、DIDV 前置物理内核 |
| 当前项目落点 | `Tools/Material`、`Tools/Particle`、`Tools/FieldSolver`、`Tools/Coupling` 已沉淀核心能力；`Toolkit/Surface Charging` 负责场景装配与运行时 glue |
| 现状判断 | Surface 主线数值物理核已完成首轮公共化，单测与回归门禁强于 SPIS 旧 Java 结构 |
| 核心缺口 | toolkit 侧仍承担较重的 config 适配与 runtime glue，同时 `Tools/*` 中 Surf 模型家族深度仍偏浅（当前以 `BasicSurfaceMaterialModel` + 有限 distrib/scaler 为主） |
| 改进方向 | 继续压缩 toolkit 中的物理装配与状态编译逻辑，并补齐 SPIS Surf 关键模型族的最小可用子集（erosion/distrib/scaler） |

### 2.3 Circ 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | `Circ/Circ`、`Circ/DIDV`、`CircField` 定义电路元件、电流分布、DIDV 组合与电路求解前置数据 |
| 当前项目落点 | `Tools/Coupling` 已有通用节点/支路/动态激励/DIDV 组合接口；`Toolkit/Surface Charging` 消费这些接口 |
| 现状判断 | 对 surface 主线而言，`Circ` 已不是当前主阻塞项 |
| 核心缺口 | `R11-D` 后核心波形与 DIDV 语义最小集已补齐；当前缺口是 SPIS 中更高阶 DIDV 组合包装类是否需要工程化引入 |
| 改进方向 | 维持当前统一组合接口，按真实场景需求再决定是否引入 `ReducableDIDV` 风格包装层，避免过度同构迁移 |

### 2.4 Solver 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | `Circuit`、`ElectroMag`、`Matter` 负责电路、泊松/电磁、粒子推进等求解任务 |
| 当前项目落点 | `Tools/Solver`、`Tools/FieldSolver`、`Tools/Particle`、`Tools/PICcore` |
| 现状判断 | solver facade 已建立，线性求解、策略解析和路由入口已有公共承载点 |
| 核心缺口 | `DensePlasmaSurfaceCharging::initialize(...)` 与 `advance(...)` 仍直接承担 solver 相关阶段编排 |
| 改进方向 | 将 solver 域在 surface 主线上的后续工作聚焦到 runtime 编排解耦，而不是继续扩展 solver 能力面 |

### 2.5 Vol 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | `BC`、`Geom`、`VolMesh`、`VolDistrib`、`VolField`、`VolInteract` 共同组成原生体域数值内核 |
| 当前项目落点 | `Tools/Boundary`、`Tools/Mesh`、`Tools/FieldSolver`、`Tools/Particle` 的一部分能力 + `Toolkit/Surface Charging` 中的 field/volume bridge、pseudo-volume runtime |
| 现状判断 | surface bridge 契约、外部体域反馈和 pseudo-volume 运行路径已清晰很多 |
| 核心缺口 | 当前“Volume 已完成”的结论只适用于 surface bridge 边界冻结，不适用于 native volume kernel 对标完成 |
| 改进方向 | 在正文中冻结“surface bridge 已成型、native volume 未完成”的口径，避免后续计划误判优先级 |

### 2.6 Util 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | 分布函数、数学函数、积分、异常、参数/IO 工具等基础设施 |
| 当前项目落点 | `Tools/Basic`、`Tools/Geometry`、`Tools/Diagnostics`、`Tools/Output`，以及 `Tools/Particle`、`Tools/Solver` 中的 surface 公共 helper |
| 现状判断 | Surface 主线直接依赖的 helper 已完成首轮公共化 |
| 核心缺口 | runtime/top 层新增辅助逻辑若没有门禁，仍可能回流到 toolkit 巨型文件 |
| 改进方向 | 将 Util 域后续工作转化为“禁止回流”规则与 gate，而不是继续扩写大而全的 util 计划 |

### 2.7 Top 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | `Simulation`、`Scenario`、`Transition`、`Plasma`、`SC` 等统一管理场景、运行时生命周期与对象关系 |
| 当前项目落点 | `SurfaceScenarioCatalog`、`SurfaceScenarioLoader`、`SurfaceSimulationRunner`、`SurfaceRuntimePlan`、`SurfaceTransitionEngine`、`DensePlasmaSurfaceCharging`、`Main/main.cpp` |
| 现状判断 | `SurfaceScenarioCatalog`、`SurfaceScenarioLoader`、`SurfaceSimulationRunner`、`SurfaceRuntimePlan`、`SurfaceTransitionEngine` 已落地，对象层第二阶段完成 |
| 核心缺口 | `R11-E` 已覆盖 `LocalTime/Spinning/SunFlux` 最小事件链；尚缺 `ConductivityEvolution`、`SourceFluxUpdater`、`SimulationParamUpdater` 等扩展事件族 |
| 改进方向 | 继续在 `SurfaceTransitionEngine` 扩展事件描述层，但保持“最小可用 + gate 驱动”策略，优先补齐会影响主线数值结果的事件 |

### 2.8 本轮重新审计结论

| 主题 | 文档原结论 | 重新审计后结论 | 正文处理 |
| --- | --- | --- | --- |
| `Solver` | 已完成 solver 编排统一 | facade 与求解路由保持稳定，当前不再是主阻塞域 | 保留为已取得成果，转入按需维护 |
| `Vol` | 已完成 volume 原生化轮 | 当前只完成 surface bridge/native 边界冻结与契约下沉，不能等同于 native volume kernel 完成 | 在正文中纠偏，并转为背景结论 |
| `Top` | 已完成 Top 对象层规划与首批实现 | `RuntimePlan + TransitionEngine` 已完成；当前缺口转为“事件级 transition 家族” | 作为 `Round 11` 主线之一 |
| `Surf` | 已完成 Surface 通用能力下沉 | 公共边界已建立，但模型家族深度与 SPIS 仍有明显差距 | 作为 `Round 11` 主线之一 |
| `Circ` | 已完成 Circuit 通用化轮 | 通用电路/激励边界已建立，但 waveform 与 DIDV 语义族仍不完整 | 作为 `Round 11` 主线之一 |
| `Util` | 已完成 | 对 surface 主线而言已完成首轮收口；后续主要是防回流维护 | 保留为已完成背景工作，转入历史归档口径 |

### 2.9 二次详细对比（SPIS class-family 级）

| 能力维度 | SPIS 代表族（证据） | 当前实现（证据） | 对齐判断 | 后续动作 |
| --- | --- | --- | --- | --- |
| Surf Material 模型族 | `spis/Surf/Material/BasicMaterialModel.java`、`spis/Surf/Material/ErosionMaterialModel.java` | `Tools/Material/include/SurfaceMaterialModel.h`（已含 `BasicSurfaceMaterialModel` 与 `ErosionSurfaceMaterialModel`，并提供 `resolveSurfaceMaterialModelVariant` 路由） | `已对齐（R11-A 已完成最小集）` | 维持默认 Basic 兼容，按需扩展 erosion 参数与产额曲线 |
| Surf Distrib 模型族 | `spis/Surf/SurfDistrib/AxisymTabulatedVelocitySurfDistrib.java`、`LocalModifiedPearsonIVSurfDistrib.java`、`NonLocalizedSurfDistrib.java`、`MultipleSurfDistrib.java` | `Tools/Particle/include/SurfaceDistributionFunction.h`（已含 `TabulatedVelocity`、`NonLocalizedHybrid`，并保留 `MaxwellianProjected/WakeAnisotropic/MultiPopulationHybrid`） | `部分对齐（R11-B 已完成最小集）` | 保持现有模型默认兼容；将 `MultipleSurfDistrib`/`LocalModifiedPearsonIV` 纳入 Round 12 候选 |
| Surf Field Scaler 模型族 | `spis/Surf/SurfField/OMLCurrentScaler.java`、`LTE_OML_CurrentScaler.java`、`FowlerNordheimCurrentScaler.java`、`VariableBarrierCurrentScaler.java` | `Tools/FieldSolver/include/SurfaceBarrierModels.h`（已含 `SecondaryRecollection/VariableBarrier/OML/LTE/FN`） | `部分对齐（R11-C 已完成核心集）` | 按场景需求评估 `AutomaticBarrier`/`MultipleCurrentScaler` 类扩展是否进入主线 |
| Circ 波形族 | `spis/Circ/Circ/SIN.java`、`PULSE.java`、`PWL.java`、`EXP.java` | `Tools/Coupling/include/SurfaceCircuitCoupling.h`（已含 `Constant/Sinusoidal/Pulse/PiecewiseLinear/Exponential`） | `已对齐（R11-D 已完成）` | 维持采样回归与 gate 监测，暂不扩展额外 waveform 语法糖 |
| Circ DIDV 语义组合族 | `spis/Circ/DIDV/WeightedSurfDIDV.java`、`MultipleSurfDIDV.java`、`SingleScalerMultipleSurfDIDV.java` | `Tools/Coupling/include/SurfaceCircuitCoupling.h` + `composeCircuitDidv(...)`（已含 `WeightedTerms`、`MultiSurface`） | `已对齐（R11-D 已完成）` | 保持组合接口统一；仅在出现真实需求时新增 reducable 包装层 |
| Top Transition 事件族 | `spis/Top/Transition/LocalTimeTransition.java`、`SpinningSpacecraft.java`、`SunFluxUpdater.java`、`ConductivityEvolution.java` | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`（已接入 `LocalTime/Spinning/SunFlux`） | `部分对齐（R11-E 已完成最小集）` | 下一步补齐 `ConductivityEvolution` 与参数更新类事件，形成扩展事件最小闭环 |
| Vol 原生体域内核 | `spis/Vol/VolMesh/ThreeDUnstructVolMesh.java`、`ThreeDUnstructVolMeshWithViewFactor.java`、`spis/Vol/VolField/EMField.java` | `Toolkit/Surface Charging` pseudo-volume + external bridge | `工程替代` | 维持 bridge 口径，不在当前主线轮次内推进 native volume kernel |

### 2.10 再次对比后的后续改进计划与清单（2026-04-14）

| 优先级 | 任务ID | 任务 | 目标落点 | 验证方式 | 完成判据 |
| --- | --- | --- | --- | --- | --- |
| `P0` | `R11-A1` | 引入 `ErosionSurfaceMaterialModel` 最小实现 | `Tools/Material/include/SurfaceMaterialModel.h`、`Tools/Material/src/SurfaceMaterialModel.cpp` | `Material_test` 新增 erosion 族用例 + surface smoke | `SurfaceMaterialModel` 支持 family 选择，erosion 电流/产额路径可被调用且默认不改旧行为 |
| `P0` | `R11-A2` | 打通 material family 配置路由与导出元数据 | `Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | object-layer + smoke + 导出 sidecar 检查 | metadata 明确导出 family，`surface`/`surface-config` 在未配置时保持 `BasicSurfaceMaterialModel` |
| `P1` | `R11-A3` | 建立 R11-A gate 证据并封板 Round 11 | `scripts/run/surface_reference_matrix_round11_subset.json`、`build/round11_surface_reference_matrix_gate/*`、`build/surface_reference_matrix_gate/*` | Round11 子集 gate + 全量 gate | 子集与全量 gate 均 `PASS`，Round11 状态更新为 `已完成` |
| `P2` | `R12-A` | Top 事件族扩展到 `ConductivityEvolution` 最小闭环 | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`、`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp` | object-layer test + smoke | 事件触发能改变导电参数路径，并具备门控与默认兼容 |
| `P2` | `R12-B` | Surf distrib 深度扩展候选评估（`MultipleSurf`/`PearsonIV`） | `Tools/Particle/include/SurfaceDistributionFunction.h`、`Tools/Particle/src/SurfaceDistributionFunction.cpp` | Particle 单测 + 对比 smoke | 至少一项候选实现达到“可路由+可回归验证” |
| `P2` | `R12-C` | Surf scaler 扩展候选评估（`AutomaticBarrier`/`MultipleCurrentScaler`） | `Tools/FieldSolver/include/SurfaceBarrierModels.h`、`Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | FieldSolver 单测 + 目标场景 smoke | 明确进入/不进入主线的工程决策，并有证据记录 |

- 执行结果（`2026-04-14`）：`R11-A1`、`R11-A2`、`R11-A3` 已完成。
- 证据：`Tools/Material/include/SurfaceMaterialModel.h`、`Tools/Material/src/SurfaceMaterialModel.cpp`、`Tools/Material/test/Material_test.cpp`、`Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`build/round11_surface_reference_matrix_gate/surface_reference_matrix_gate.json`、`build/surface_reference_matrix_gate.json`。

### 2.11 对标完成预期声明（执行前置）

- 本文档后续全量对标任务执行完毕后，可判定当前项目与 SPIS 在“表面充电分析能力”上达到**基本相近**。
- “基本相近”采用可审计判定口径，必须同时满足下表条件：

| 判定项 | 必须满足条件 |
| --- | --- |
| Top 事件族 | `Round 12` 任务全部标记 `已完成`，并具备对象层 + smoke + sidecar 证据 |
| Surf 分布族 | `Round 13` 任务全部标记 `已完成`，新增分布族可路由、可测试、可回归 |
| Surf scaler 族 | `Round 14` 任务全部标记 `已完成`，自动/组合策略具备可观测元数据 |
| Circ DIDV 高阶组合 | `Round 15` 任务全部标记 `已完成`，并通过 Coupling 单测与场景 smoke |
| 回归封板 | `Round 18` 的子集 gate 与全量 gate 均 `PASS`，并完成 contract 校验 |

- 若以上条件满足，则本项目在 Surface 主线的能力、模式与回归证据上与 SPIS 保持同级可用；若 `Round 17` 同时完成，则可进一步覆盖 native volume parity 轨道。

## 3. 改进总原则

- 新的通用数值能力优先进入 `Tools/*`，不得默认继续堆叠到 `Toolkit/Surface Charging`。
- `Toolkit/Surface Charging` 只保留：
  - 场景装配
  - 兼容适配
  - 产物导出
  - route/preset/runtime 组织
- 凡是未来可能被 surface 主线中多个运行路径复用的能力，都不得继续私有化在巨型 runtime 文件中。
- 回归门禁与物理等价必须持续维持，但“回归通过”不等于“边界已经合理”。
- 文档中的改进任务优先围绕“职责收口”和“复用边界”展开，而不是无节制加新功能。
- 已完成但不再属于当前主线的工作，进入历史归档，不继续占用正文轮次总表。

## 4. 分轮执行计划

### 4.1 状态词典

| 状态 | 含义 |
| --- | --- |
| `未开始` | 尚未启动 |
| `已盘点` | 已完成现状摸底，但未进入设计或实现 |
| `设计中` | 正在确定接口、边界、拆分方案 |
| `实现中` | 已开始代码级推进 |
| `已完成` | 代码、测试、文档、验收均已收口 |
| `历史归档` | 已形成稳定结论或已有阶段性落地，但不再属于当前正文主线 |
| `已冻结` | 确认保留现状，暂不推进，但结论已记录 |
| `受阻` | 存在外部依赖或关键分歧，当前无法继续 |

### 4.2 轮次总表

| 轮次 | 主题 | 当前状态 | 最后更新 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- |
| `Round 0` | 计划基线轮 | `已完成` | `2026-04-13` | 本文档已建立六域差距表、状态规则、分轮模板 | 进入 `Round 1` 盘点 Surface 私有残留逻辑 |
| `Round 1` | Surface 边界收口轮 | `已完成` | `2026-04-13` | 本文档 4.4 节残留矩阵；`Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 进入 `Round 2`，按矩阵逐类抽离剩余通用能力 |
| `Round 2` | Surface 通用能力下沉轮 | `已完成` | `2026-04-13` | `Tools/Material/include/SurfaceMaterialModel.h`；`Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/FieldSolver/include/SurfaceBarrierModels.h`；`Tools/Coupling/include/SurfaceCurrentLinearizationKernel.h`；`Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 进入 `Round 3`，开始 Circuit 通用化设计与实现 |
| `Round 3` | Circuit 通用化轮 | `已完成` | `2026-04-13` | `Tools/Coupling/include/SurfaceCircuitCoupling.h`；`Tools/Coupling/src/SurfaceCircuitCoupling.cpp`；`Tools/Coupling/test/Coupling_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 进入 `Round 9`，推进 `RuntimePlan + TransitionEngine` 拆分 |
| `Round 9` | RuntimePlan + TransitionEngine 拆分轮 | `已完成` | `2026-04-14` | `Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`build/surface_reference_matrix_gate.json` | 进入 `Round 10` 等价验收收口 |
| `Round 10` | Surface 功能等价验收轮 | `已完成` | `2026-04-14` | `Toolkit/Surface Charging/docs/surface_round10_equivalence_matrix.md`；`build/surface_reference_matrix_gate.json`；`build/surface_reference_matrix_gate.md`；`scripts/run/surface_reference_matrix.json` | 进入 `Round 11` 深差距补齐 |
| `Round 11` | SPIS 深差距补齐轮（Surf/Circ/Top） | `已完成` | `2026-04-14` | 本文档 2.9 节对照矩阵；`Tools/Material/include/SurfaceMaterialModel.h`；`Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Coupling/include/SurfaceCircuitCoupling.h`；`Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`build/round11_surface_reference_matrix_gate/surface_reference_matrix_gate.json`；`build/surface_reference_matrix_gate.json` | 进入 `Round 12-18` 全量对标执行清单（4.11 节） |
| `Round 12` | Top Transition 全量对标轮 | `实现中` | `2026-04-14` | 见 4.11 节 `R12-A*` 清单；`R12-A1`、`R12-A2`、`R12-A3`、`R12-A4` 已完成 | 继续推进 `R12-A5` 到 `R12-A8` |
| `Round 13` | Surf Distrib 全量对标轮 | `未开始` | `2026-04-14` | 见 4.11 节 `R13-B*` 清单 | 完成分布族扩展并进入 `Round 14` |
| `Round 14` | Surf Scaler 全量对标轮 | `未开始` | `2026-04-14` | 见 4.11 节 `R14-C*` 清单 | 完成 scaler 族扩展并进入 `Round 15` |
| `Round 15` | Circ DIDV 高阶组合对标轮 | `未开始` | `2026-04-14` | 见 4.11 节 `R15-D*` 清单 | 完成高阶 DIDV 组合并进入 `Round 16` |
| `Round 16` | Material 参数资产全量对标轮 | `未开始` | `2026-04-14` | 见 4.11 节 `R16-E*` 清单 | 完成参数资产闭环并进入 `Round 17` |
| `Round 17` | Vol Native Parity 轨道轮 | `未开始` | `2026-04-14` | 见 4.11 节 `R17-F*` 清单 | 完成 native parity 最小闭环并进入 `Round 18` |
| `Round 18` | 全量回归封板与文档闭环轮 | `未开始` | `2026-04-14` | 见 4.11 节 `R18-G*` 清单 | 子集+全量 gate `PASS` 后封板 |
| `Round 4-8` | Solver/Volume/Util/Top 历史工作 | `历史归档` | `2026-04-14` | 见 4.10 节历史归档 | 仅保留历史证据，不再作为正文主线轮次 |

### 4.3 Round 0：计划基线轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 完成六域差距盘点，冻结对比边界和状态规则 |
| 输入 | `ref/SPIS/num-master`、`Toolkit/Surface Charging`、`Tools/*`、`scripts/run/surface_reference_matrix.json`、`documentation/contracts/**` |
| 实施项 | 建立六域总表；标出已有模块、缺口模块、私有化热点模块；建立统一状态词典和轮次更新规范 |
| 输出 | 本文档主框架；六域差距表；状态规则表；轮次模板 |
| 依赖 | 无 |
| 验收 | 后续执行者不需要再决定“比什么、怎么记状态” |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-13` |
| 证据 | 本文档第 2、4、5、6、7 节 |

### 4.4 Round 1：Surface 边界收口轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 把 Surface 已有成果从“能力存在”推进到“边界清晰” |
| 输入 | `Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、当前 `Tools/Material`、`Tools/Particle`、`Tools/FieldSolver`、`Tools/Coupling` 接口 |
| 实施项 | 盘点仍未下沉到 `Tools/*` 的通用逻辑；按 `Material`、`Particle`、`FieldSolver`、`Coupling` 四类归档；建立 Surface 私有逻辑残留矩阵 |
| 输出 | Surface 私有逻辑残留矩阵；下一轮抽离任务清单 |
| 依赖 | `Round 0` |
| 验收 | 每一项残留逻辑都能指向唯一目标落点 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-13` |
| 证据 | `Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |

#### Round 1 任务卡

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R1-A` | 盘点 `Material` 类残留通用逻辑 | `已完成` | `2026-04-13` | `Round 0` | 残留材料逻辑清单 | 代码检视 | `SurfaceChargingKernelFramework.cpp` | 建立 `Tools/Material` 迁移项 |
| `R1-B` | 盘点 `Particle` 类残留通用逻辑 | `已完成` | `2026-04-13` | `Round 0` | 残留分布/采样逻辑清单 | 代码检视 | `SurfaceChargingKernelFramework.cpp` | 建立 `Tools/Particle` 迁移项 |
| `R1-C` | 盘点 `FieldSolver` 类残留通用逻辑 | `已完成` | `2026-04-13` | `Round 0` | 残留 barrier/recollection/field helper 清单 | 代码检视 | `SurfaceChargingKernelFramework.cpp`；`DensePlasmaSurfaceCharging.cpp` | 建立 `Tools/FieldSolver` 迁移项 |
| `R1-D` | 盘点 `Coupling` 类残留通用逻辑 | `已完成` | `2026-04-13` | `Round 0` | 残留 DIDV/电路装配 helper 清单 | 代码检视 | `SurfaceChargingKernelFramework.cpp`；`DensePlasmaSurfaceCharging.cpp` | 建立 `Tools/Coupling` 迁移项 |

### 4.5 Round 2：Surface 通用能力下沉轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 把可复用的 surface 数值能力继续抽离到 `Tools/*` |
| 输入 | `Round 1` 的 Surface 私有逻辑残留矩阵 |
| 实施项 | 抽离剩余分布函数、材料响应、barrier/recollection、线性化 helper；清理 toolkit 私有重复实现与 wrapper；为每个抽离点补对应单测 |
| 输出 | `Tools/*` 新增或完善的公共接口；toolkit 侧减少私有重复实现 |
| 依赖 | `Round 1` |
| 验收 | Surface 主线运行时不再依赖 toolkit 私有公式实现完成核心数值步骤 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-13` |
| 证据 | `Tools/Material/include/SurfaceMaterialModel.h`；`Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/FieldSolver/include/SurfaceBarrierModels.h`；`Tools/Coupling/include/SurfaceCurrentLinearizationKernel.h`；`cmake --build build --config Debug`；`ctest --test-dir build -C Debug -R "Material\|SurfaceDistributionFunction\|SurfaceFieldVolumeBridge\|SurfaceBarrierModels" --output-on-failure` |

#### Round 2 任务卡

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R2-A` | 抽离剩余 `Material` 通用 helper | `已完成` | `2026-04-13` | `R1-A` | `Tools/Material` 公共接口补齐 | 单测 + toolkit 替换点检查 | `Tools/Material/include/SurfaceMaterialModel.h` | 保留适配层仅承载场景语义 |
| `R2-B` | 抽离剩余 `Particle` 通用 helper | `已完成` | `2026-04-13` | `R1-B` | `Tools/Particle` 公共接口补齐 | 单测 + toolkit 替换点检查 | `Tools/Particle/include/SurfaceDistributionFunction.h` | toolkit 仅做 config->request 适配 |
| `R2-C` | 抽离剩余 `FieldSolver` 通用 helper | `已完成` | `2026-04-13` | `R1-C` | `Tools/FieldSolver` 公共接口补齐 | 单测 + toolkit 替换点检查 | `Tools/FieldSolver/include/SurfaceBarrierModels.h` | 继续减少场状态细节重复持有 |
| `R2-D` | 抽离剩余 `Coupling` 通用 helper | `已完成` | `2026-04-13` | `R1-D` | `Tools/Coupling` 公共接口补齐 | 单测 + toolkit 替换点检查 | `Tools/Coupling/include/SurfaceCurrentLinearizationKernel.h` | 进入 `Round 3` |
| `R2-E` | 清理 toolkit 残留重复实现和 wrapper | `已完成` | `2026-04-13` | `R2-A`、`R2-B`、`R2-C`、`R2-D` | toolkit 侧重复实现减少 | diff + smoke test | `SurfaceChargingKernelFramework.cpp` | `Round 2` 收口完成 |

### 4.6 Round 3：Circuit 通用化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 从 surface-specific circuit 走向 SPIS `Circ` 对标层 |
| 输入 | `Tools/Coupling` 当前接口、`Toolkit/Surface Charging` 电路装配热点、`ref/SPIS/num-master/src/main/java/spis/Circ/**` |
| 实施项 | 建立通用节点/支路/元件抽象；建立通用 DIDV 组合接口；增加波形源/外加激励模型入口；让 `Surface Charging` 消费公共接口而非内嵌代数逻辑 |
| 输出 | `Tools/Coupling` 的 surface 主线公共化边界 |
| 依赖 | `Round 2` |
| 验收 | surface 主线中的电路代数不再依赖 toolkit 手工组装 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-13` |
| 证据 | `Tools/Coupling/include/SurfaceCircuitCoupling.h`；`Tools/Coupling/src/SurfaceCircuitCoupling.cpp`；`Tools/Coupling/test/Coupling_test.cpp`；`ctest --test-dir build -C Debug -R "CouplingTest\\." --output-on-failure` |

#### Round 3 任务卡

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R3-A` | 定义通用节点/支路/元件抽象 | `已完成` | `2026-04-13` | `Round 2` | `CircuitAssembly` 与 `assembleSurfaceCircuit` 适配入口 | 编译 + Coupling 单测 | `Tools/Coupling/include/SurfaceCircuitCoupling.h` | 已完成 |
| `R3-B` | 定义通用 DIDV 组合接口 | `已完成` | `2026-04-13` | `Round 2` | `composeCircuitDidv` 与线性化路径接入 | 编译 + Coupling 单测 + Surface smoke | `SurfaceCurrentLinearizationKernel.cpp` | 已完成 |
| `R3-C` | 增加波形源/外加激励模型入口 | `已完成` | `2026-04-13` | `R3-A` | `CircuitExcitationDescriptor` 与动态覆盖入口 | 编译 + Coupling 单测 | `SurfaceCircuitCoupling.cpp` | 已完成 |
| `R3-D` | 将 `Surface Charging` 电路代数迁移到公共接口 | `已完成` | `2026-04-13` | `R3-A`、`R3-B` | `composeCircuitKernelInput` 与主线路径接入 | 编译 + Coupling 单测 + Surface smoke | `DensePlasmaSurfaceCharging.cpp` | 进入 `Round 9` |

### 4.7 Round 9：RuntimePlan + TransitionEngine 拆分轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `DensePlasmaSurfaceCharging` 中的顶层编排拆成 `SurfaceRuntimePlan + SurfaceTransitionEngine + runtime core` 三层，使对象层真正落地到初始化与步进主路径 |
| 输入 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`、`Toolkit/Surface Charging/include/SurfaceSimulationRunner.h`、`Main/SurfaceScenarioLoader.cpp`、当前 `SurfaceScenarioCatalog/Loader/Runner` 首批对象层 |
| 实施项 | 冻结 `scenario -> runtime plan -> transition engine -> runtime core` 四层边界；从 `initialize(...)` 中抽离编译态数据与策略矩阵；从 `advance(...)` 中抽离 route/strategy/legacy/feature gate 转场；补齐对象层单测与 smoke/gate 验证 |
| 输出 | `SurfaceRuntimePlan`、`SurfaceRuntimeContext`、`SurfaceTransitionEngine` 的落地方案与详细任务清单；正文状态从“仅有规划”推进到“可执行实现任务” |
| 依赖 | `Round 3`；历史归档中的 solver/top 首批成果 |
| 验收 | `DensePlasmaSurfaceCharging` 不再直接承担场景编译、阶段转场、feature gate 规则汇总三类职责；现有 CLI、preset、route、CSV、sidecar 契约保持兼容 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-14` |
| 证据 | `Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/docs/architecture.md`；`build/bin/SurfaceCharging_object_layer_test.exe`；`build/bin/SurfaceCharging_smoke_test.exe`；`build/surface_reference_matrix_gate.json` |

#### Round 9 任务卡

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R9-A` | 冻结 `RuntimePlan` 边界与字段清单 | `已完成` | `2026-04-14` | `Round 8` 历史成果 | `SurfaceRuntimePlan` 已落地为 scenario 编译态对象，当前冻结字段包括：normalized config、route / strategy / legacy execution、solver policy flags、bridge enable、spectrum->plasma 派生态、preset metadata | 代码检视 + 文档复核 | `Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/docs/architecture.md` | 进入 `R9-B` / `R9-C` |
| `R9-B` | 将 `initialize(...)` 拆为 scenario 编译与 runtime context 组装两阶段 | `已完成` | `2026-04-14` | `R9-A` | `compileSurfaceRuntimePlan(...)` 已接入 `SurfaceSimulationRunner` 与 `DensePlasmaSurfaceCharging::initialize(...)`；`initialize(plan)` 仅消费编译态配置并组装 runtime core/context | 编译 + object-layer 单测 + targeted smoke | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`SurfaceChargingSmokeTest.ToolkitInitializesAdvancesAndExports` | 已完成 |
| `R9-C` | 抽离 `SurfaceTransitionEngine` 与 route/strategy/legacy 转场矩阵 | `已完成` | `2026-04-14` | `R9-A` | `SurfaceAdvanceTransition` / `evaluateSurfaceAdvanceTransition(...)` 已落地；`advance(...)` 入口不再直接持有 legacy replay、reference-circuit、shared runtime、solver policy 主决策 | object-layer 单测 + targeted smoke | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`SurfaceTransitionEngineTest.LegacyReplayRequiresLegacyRouteAndReplayMode`；`SurfaceChargingSmokeTest.LegacyBenchmarkExecuteModeAdvancesWithoutReplay` | 已完成 |
| `R9-D` | 将 `advance(...)` 拆为 prepare / solve / commit 阶段管线 | `已完成` | `2026-04-14` | `R9-B`、`R9-C` | `advance(...)` 已形成 prepare/solve/commit 三阶段边界，求解路径与导出契约保持兼容 | 编译 + smoke + 回归 gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`SurfaceChargingSmokeTest.LegacyBenchmarkExecuteModeAdvancesWithoutReplay`；`SurfaceChargingSmokeTest.ScdatUnifiedRouteIgnoresLegacyFallbackModeOutsideLegacyRoute` | 已完成 |
| `R9-E` | 收敛 runtime state / telemetry / sidecar 触发边界 | `已完成` | `2026-04-14` | `R9-D` | runtime status/history 写回统一收敛到 commit 阶段；legacy replay 触发边界收敛到 transition 决策（legacy route + replay mode） | object-layer + smoke + artifact contract 检查 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`SurfaceChargingSmokeTest.SurfaceReferenceAndRuntimeContractsExportMetadata` | 已完成 |
| `R9-F` | 建立 `Round 9` 验收门禁与详细回归清单 | `已完成` | `2026-04-14` | `R9-B`、`R9-C`、`R9-D` | transition replay 边界门禁、对象层与关键 smoke 门禁、reference matrix gate 脚本化记录已固定 | object-layer + CLI smoke + reference matrix gate | `Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`SurfaceTransitionEngineTest.LegacyReplayRequiresLegacyRouteAndReplayMode`；`SurfaceChargingSmokeTest.ScdatUnifiedRouteIgnoresLegacyFallbackModeOutsideLegacyRoute`；`build/surface_reference_matrix_gate.json`；`build/surface_reference_matrix_gate.md` | 进入 `Round 10` 收口 |

#### Round 9 详细任务拆解

| 子系统 | 当前热点 | Round 9 必须完成的收口动作 | 不在本轮做的事 |
| --- | --- | --- | --- |
| `Scenario -> Runtime` 编译边界 | `Main/SurfaceScenarioLoader.cpp`、`DensePlasmaSurfaceCharging::initialize(...)` | 把 JSON/preset/overlay 的结果编译成 `SurfaceRuntimePlan`，并冻结字段归属 | 不重写外部 JSON schema；不改 preset 名称与 alias |
| `Transition / Stage` 编排边界 | `DensePlasmaSurfaceCharging::advance(...)` | 把 route/strategy/legacy/feature gate 判断与阶段切换抽成 `SurfaceTransitionEngine` | 不在本轮把所有物理求解实现再下沉到新模块 |
| `Runtime Context` ownership | `DensePlasmaSurfaceCharging` 成员集合、orchestrator 组件实例 | 明确 context 持有组件实例、共享状态与求解器句柄，runtime core 只消费 context | 不在本轮做全量状态仓库化 |
| `Runner / Export` 边界 | `SurfaceSimulationRunner`、主 CSV 与 sidecar 导出路径 | runner 只负责生命周期驱动与导出触发，不再回收 runtime 内部判断职责 | 不在本轮迁移全部 sidecar 导出实现 |
| `Regression Gate` | object-layer test、surface smoke、reference matrix | 为新边界建立可重复的门禁组合，保证拆分不破坏现有行为 | 不新增与 surface 主线无关的跨 toolkit gate |

#### Round 9 当前已冻结的 `RuntimePlan` 字段归属

| 类别 | 当前归属到 `SurfaceRuntimePlan` 的字段/信息 | 禁止回流到 giant runtime 的内容 |
| --- | --- | --- |
| route / mode | `runtime_route`、`surface_pic_strategy`、`legacy_input_adapter_kind`、`surface_pic_runtime_kind`、`surface_instrument_set_kind`、`legacy_benchmark_execution_mode`、`benchmark_source` | 在 `initialize(...)` 中重新判定 route / strategy 主分支 |
| solver policy | `solver_policy_flags`；由策略编译得到的 `internal_substeps` 下限与 shared global coupled policy material scalars | 在 `advance(...)` 中重新解析 coupling mode / convergence policy token |
| bridge contract | `enable_external_field_solver_bridge`、`enable_external_volume_solver_bridge` 与其编译后的 normalized config 结果 | 在 runtime core 中重新改写 bridge enable / disable 语义 |
| plasma-derived state | `has_electron_spectrum` / `has_ion_spectrum` 派生的 plasma density / temperature / ion mass 同步结果 | 在 step-time 路径重新做 spectrum -> plasma moments 写回 |
| source metadata | preset `name`、`description` | 把 preset / loader 元数据再次塞回 runtime mutable state |

#### Round 9 当前禁止回流清单

- 不再向 `DensePlasmaSurfaceCharging::initialize(...)` 回填新的 scenario-token 解析逻辑。
- 不再把 solver policy token 归一化、legacy execute config 展开、spectrum-to-plasma 派生写回直接堆回 giant runtime。
- 不把 history 容器、telemetry、export side effects 并入 `SurfaceRuntimePlan`。
- 不在 `SurfaceSimulationRunner` 中绕过 plan 直接操纵 runtime route / legacy mode。
- 不把 step 前的 route / legacy / shared runtime / solver policy 总判断重新散回 `advance(...)` 巨型入口。

### 4.8 Round 10：Surface 功能等价验收轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 在 `Round 9` 完成后，系统性确认当前项目表面充电主线与 SPIS 在功能覆盖、模式行为、关键输出和参考场景结果上是否已达到“基本一致” |
| 输入 | `ref/GEO/**`、`ref/LEO/**`、`ref/mat/**`、`ref/SPIS/num-master/src/main/java/spis/Surf/**`、`ref/SPIS/num-master/src/main/java/spis/Circ/**`、`ref/SPIS/num-master/src/main/java/spis/Top/**`、`scripts/run/surface_reference_matrix.json`、现有 smoke/object-layer 测试 |
| 实施项 | 定义“基本一致”的验收口径；把 surface 主线功能拆成能力矩阵、运行模式矩阵、输出矩阵和参考场景矩阵；逐项标注已一致、近似一致、待补齐、仅工程替代四种状态；把未一致项转成后续实现任务 |
| 输出 | 一份可审计的 surface 功能等价矩阵；一组回归门禁；一批后续补齐任务，不再用泛化表述代替对标结论 |
| 依赖 | `Round 9` |
| 验收 | 能明确回答“当前项目与 SPIS 在表面充电功能上哪些已经基本一致，哪些还没有”，且结论有测试、artifact 或参考场景证据支撑 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-14` |
| 证据 | `Toolkit/Surface Charging/docs/surface_round10_equivalence_matrix.md`；`build/surface_reference_matrix_gate.json`；`build/surface_reference_matrix_gate.md`；`scripts/run/surface_reference_matrix.json`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`documentation/contracts/**` |

#### Round 10 任务卡

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R10-A` | 定义“基本一致”验收口径 | `已完成` | `2026-04-14` | `Round 9` | 已固定四类结论：`已一致`、`近似一致`、`待补齐`、`工程替代`，并明确证据门槛 | 文档审查 | `Toolkit/Surface Charging/docs/surface_round10_equivalence_matrix.md`；`scripts/run/surface_reference_matrix.json` | 进入 `R10-B` |
| `R10-B` | 建立 surface 功能能力矩阵 | `已完成` | `2026-04-14` | `R10-A` | 已形成材料/粒子/场势垒/电路/legacy/shared runtime/bridge 七类能力矩阵并完成证据闭环 | 代码检视 + 对标审查 + smoke | `Toolkit/Surface Charging/docs/surface_round10_equivalence_matrix.md`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build/surface_reference_matrix_gate.json` | 已完成 |
| `R10-C` | 建立运行模式与输出矩阵 | `已完成` | `2026-04-14` | `R10-A` | 已覆盖 `surface`、`surface-config`、preset/replay、legacy execute/replay、shared-coupled、bridge sidecar 的输入/行为/输出结论并通过门禁验证 | CLI smoke + artifact 检查 | `Toolkit/Surface Charging/docs/surface_round10_equivalence_matrix.md`；`Main/main.cpp`；`SurfaceSimulationRunner`；`build/surface_reference_matrix_gate.json`；`documentation/contracts/**` | 已完成 |
| `R10-D` | 建立参考场景等价矩阵 | `已完成` | `2026-04-14` | `R10-B`、`R10-C` | GEO/LEO/material 与 surface-pic 参考场景矩阵已刷新为全量 `PASS`，并形成可审计 gate 报告 | reference matrix gate + smoke + 人工审查 | `Toolkit/Surface Charging/docs/surface_round10_equivalence_matrix.md`；`build/surface_reference_matrix_gate.json`；`build/surface_reference_matrix_gate.md`；`scripts/run/surface_reference_matrix.json` | 已完成 |
| `R10-E` | 将未一致项转为实现 backlog | `已完成` | `2026-04-14` | `R10-D` | 本轮未一致项已收敛，历史 backlog 条目关闭并归档为“已收口”记录 | 文档审查 + gate 复核 | `Toolkit/Surface Charging/docs/surface_round10_equivalence_matrix.md`；`build/surface_reference_matrix_gate.json` | 转入 `Round 11+` 按需维护 |

#### Round 10 详细验收维度

| 维度 | 必须回答的问题 | 合格证据 |
| --- | --- | --- |
| 能力覆盖 | 当前 surface 主线是否具备 SPIS 对应的核心表面充电能力 | 代码路径 + smoke/test + contract |
| 模式兼容 | preset、route、strategy、legacy replay/execute 是否语义兼容 | CLI smoke + object-layer test + artifact |
| 输出兼容 | 主 CSV、graph/monitor/field/volume sidecar、benchmark artifact 是否可持续兼容 | `documentation/contracts/**` + 导出结果检查 |
| 场景等价 | GEO/LEO/material 参考场景是否达到通过或可解释的近似一致 | `surface_reference_matrix` + smoke/gate + 偏差说明 |
| 剩余缺口 | 仍未一致的点是否已经转成明确任务，而不是口头结论 | 本文档后续 backlog |

### 4.9 Round 11：SPIS 深差距补齐轮（Surf/Circ/Top）

| 字段 | 内容 |
| --- | --- |
| 目标 | 在保持 `Round 10` 兼容结论不回退的前提下，补齐 SPIS 对照中最关键的 Surf/Circ/Top 深度能力缺口 |
| 输入 | 本文档 2.9 节 class-family 对照；`ref/SPIS/num-master/src/main/java/spis/Surf/**`；`ref/SPIS/num-master/src/main/java/spis/Circ/**`；`ref/SPIS/num-master/src/main/java/spis/Top/Transition/**`；当前 `Tools/*` 与 `SurfaceTransitionEngine` |
| 实施项 | 引入可插拔 material/distrib/scaler 家族扩展点；补齐 `PWL/EXP` 波形与 DIDV 语义组合；增加 Top 事件级 transition 描述与执行链路；将新增能力纳入 reference matrix gate |
| 输出 | `Round 11` 差距闭环任务卡与最小可用实现路径；更新后的等价矩阵证据 |
| 依赖 | `Round 10` |
| 验收 | 新增能力通过单测 + smoke + reference matrix 子集门禁，且现有 route/preset/CSV/sidecar 契约不回归 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-14` |
| 证据 | `Tools/Material/include/SurfaceMaterialModel.h`；`Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/FieldSolver/include/SurfaceBarrierModels.h`；`Tools/Coupling/include/SurfaceCircuitCoupling.h`；`Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`scripts/run/surface_reference_matrix.json` |

#### Round 11 任务卡

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R11-A` | 扩展 Surf Material 家族（含 erosion 路径） | `已完成` | `2026-04-14` | `Round 10` | 已增加 `ErosionSurfaceMaterialModel` 与 family 选择路由，默认 Basic 兼容保持不变 | Material 单测 + smoke + matrix gate | `Tools/Material/include/SurfaceMaterialModel.h`；`Tools/Material/src/SurfaceMaterialModel.cpp`；`Tools/Material/test/Material_test.cpp`；`Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build/round11_surface_reference_matrix_gate/surface_reference_matrix_gate.json`；`build/surface_reference_matrix_gate.json` | `Round 11` 已封板 |
| `R11-B` | 扩展 Surf Distrib 家族（tabulated-velocity / non-localized） | `已完成` | `2026-04-14` | `R11-A` | 已新增 `TabulatedVelocity` 与 `NonLocalizedHybrid` 分布族，并接入 build/synthesis 路由（默认兼容） | Particle 单测 + smoke | `Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`SurfaceDistributionFunctionTest.BuildAndSynthesisSupportTabulatedVelocityAndNonLocalizedModels`；`SurfaceChargingSmokeTest.DistributionModelSwitchRebalancesReferencePatchCurrents` | 进入 `R11-C` |
| `R11-C` | 扩展 Surf Field Scaler 家族（OML/LTE/FN） | `已完成` | `2026-04-14` | `R11-A` | 已补齐 OML/LTE/FN 三类 scaler 与材料参数路由（默认兼容） | FieldSolver 单测 + smoke | `Tools/FieldSolver/include/SurfaceBarrierModels.h`；`Tools/FieldSolver/src/SurfaceBarrierModels.cpp`；`Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`；`SurfaceChargingSmokeTest.ReferenceModelBodyPatchBarrierContrastStaysWithinExpectedRanges` | 进入 `R11-D` |
| `R11-D` | 补齐 Circ 波形与 DIDV 语义族 | `已完成` | `2026-04-14` | `Round 10` | 波形 `PWL/EXP` 与 DIDV weighted/multi-surface 语义组合均已接入并通过回归 | Coupling 单测 + surface smoke | `Tools/Coupling/include/SurfaceCircuitCoupling.h`；`Tools/Coupling/src/SurfaceCircuitCoupling.cpp`；`Tools/Coupling/test/Coupling_test.cpp`；`SurfaceChargingSmokeTest.ReferenceDidvMatchesFiniteDifferenceAcrossGeoLeoRoutes` | 进入 `R11-E` |
| `R11-E` | 增加 Top 事件级 Transition 机制 | `已完成` | `2026-04-14` | `R11-D` | `SurfaceTransitionEngine` 已接入 LocalTime/Spinning/SunFlux 事件链，并对 pic recalibration 增加 SunFlux 门控 | object-layer + smoke | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`SurfaceTransitionEngineTest.AppliesLocalTimeSpinningAndSunFluxEventChain`；`SurfaceChargingSmokeTest.LegacyBenchmarkExecuteModeAdvancesWithoutReplay` | 进入 `R11-F` |
| `R11-F` | 更新 reference matrix 与等价矩阵证据 | `已完成` | `2026-04-14` | `R11-B`、`R11-C`、`R11-D`、`R11-E` | `Round 11` 子集矩阵 gate 与全量 `surface_reference_matrix` gate 均已 `PASS`，并完成等价矩阵文档回填 | matrix gate + 文档审查 | `scripts/run/surface_reference_matrix_round11_subset.json`；`build/round11_surface_reference_matrix_gate/surface_reference_matrix_gate.json`；`build/round11_surface_reference_matrix_gate/surface_reference_matrix_gate.md`；`scripts/run/surface_reference_matrix.json`；`build/surface_reference_matrix_gate/surface_reference_matrix_gate.json`；`build/surface_reference_matrix_gate/surface_reference_matrix_gate.md`；`Toolkit/Surface Charging/docs/surface_round10_equivalence_matrix.md` | 待 `R11-A` 收口后归档 Round 11 |

### 4.10 历史归档：已完成但不再属于当前正文主线的轮次

| 轮次 | 归档主题 | 归档原因 | 核心结论 | 主要证据 |
| --- | --- | --- | --- | --- |
| `Round 4` | Solver 编排统一轮 | solver facade 与边界重画已形成稳定结论，但当前主线问题已转为 runtime 编排拆分 | `Tools/Solver` 已承担求解策略解析、路由与线性求解公共入口；surface 主线剩余问题在 `RuntimePlan + TransitionEngine` | `Tools/Solver/include/SurfaceSolverFacade.h`；`Tools/Solver/src/SurfaceSolverFacade.cpp`；`Tools/Solver/test/SurfaceSolverFacade_test.cpp`；`DensePlasmaSurfaceCharging.cpp` |
| `Round 5` | Volume 边界冻结轮 | 现阶段只需要保留 bridge/native 纠偏结论，不再继续在正文展开 volume 子域任务 | 当前已完成 surface bridge 契约和边界冻结，但不能等同于 native volume kernel 完成 | `Tools/Mesh/include/SurfaceVolumeBridgeTypes.h`；`Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/Particle/include/SurfaceVolumeDistributionContracts.h` |
| `Round 6` | Util 公共基础能力轮 | Surface 主线所需 helper 首轮归位已完成，后续主要是防回流治理 | Numeric aggregation、string token、array math、distribution helper 已完成公共化 | `Tools/Basic/include/NumericAggregation.h`；`Tools/Basic/include/StringTokenUtils.h`；`Tools/Basic/include/ArrayMath3.h`；`Tools/Particle/include/SurfaceDistributionFunction.h` |
| `Round 7` | Top 运行时对象层规划轮 | 规划已被 `Round 8` 首批实现吸收，正文无需继续保留规划轮次 | 已形成 `ScenarioCatalog / ScenarioLoader / RuntimePlan / RuntimeContext / TransitionEngine / SimulationRunner` 对象层蓝图 | `Toolkit/Surface Charging/docs/architecture.md`；`Toolkit/Surface Charging/docs/surface_num_alignment_refactor_plan.md` |
| `Round 8` | Top 对象层首批实现轮 | 入口降耦已完成，当前主线进入第二阶段 runtime 拆分 | `SurfaceScenarioCatalog`、`SurfaceScenarioLoader`、`SurfaceSimulationRunner` 已落地，CLI 与 preset 兼容性已验证 | `Toolkit/Surface Charging/include/SurfaceScenarioCatalog.h`；`Toolkit/Surface Charging/include/SurfaceSimulationRunner.h`；`Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp` |

#### 历史归档对当前主线的约束

- 保持 `SurfaceChargingConfig` 外部字段语义和默认行为不变。
- 保持 preset canonical name 与现有 alias 映射不变。
- 保持 `SurfaceRuntimeRoute`、`SurfacePicStrategy`、`LegacyBenchmarkExecutionMode` 外部行为兼容。
- 保持主 CSV 与 benchmark / graph / monitor / field / volume sidecar 合同兼容。

### 4.11 Round 12-18：SPIS 全量对标执行清单（可逐项标记完成）

#### 4.11.1 使用说明

- 本节任务默认初始状态为 `未开始`；执行后按 5.2 节规则更新为 `设计中` / `实现中` / `已完成` / `受阻`。
- 任何任务标记 `已完成` 时，必须填写 `完成时间` 与 `证据` 两列。
- 本节任务全部完成并通过 `Round 18` 封板后，按 2.11 节口径判定为“与 SPIS 表面充电分析能力基本相近”。

#### 4.11.2 任务清单

| 阶段 | 任务ID | 任务 | 目标落点 | 当前状态 | 完成时间 | 证据 | 完成判据 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `Round 12` | `R12-A1` | 统一 Top 事件描述结构与门控顺序 | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h` | `已完成` | `2026-04-14` | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp` | 新增事件描述结构可承载扩展事件族 |
| `Round 12` | `R12-A2` | 实现 `ConductivityEvolution` 事件 | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp` | `已完成` | `2026-04-14` | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 导电参数动态更新生效且默认兼容 |
| `Round 12` | `R12-A3` | 实现 `SourceFluxUpdater` 事件 | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp` | `已完成` | `2026-04-14` | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 入射源参数可按事件动态更新 |
| `Round 12` | `R12-A4` | 实现 `SimulationParamUpdater` 事件 | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp` | `已完成` | `2026-04-14` | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 求解参数可按事件动态更新 |
| `Round 12` | `R12-A5` | 实现 `SheathOrPresheathPoissonBCUpdater` 事件 | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp` | `未开始` | `N/A` | `-` | 鞘层边界策略可切换并可回归 |
| `Round 12` | `R12-A6` | 实现 `RCCabsSCUpdater` / `VcrossBfieldUpdater` / `BasicEclipseExit` / `TransientArtificialSources` / `LangmuirProbeTransition` | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp` | `未开始` | `N/A` | `-` | 事件可独立启停并进入执行链 |
| `Round 12` | `R12-A7` | 补齐对象层事件测试 | `Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp` | `未开始` | `N/A` | `-` | 新增事件用例通过 |
| `Round 12` | `R12-A8` | 补齐事件 smoke 与 sidecar 可观测字段 | `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `未开始` | `N/A` | `-` | smoke 通过且导出字段可审计 |
| `Round 13` | `R13-B1` | 实现 `MultipleSurfDistrib` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `未开始` | `N/A` | `-` | 分布族可路由且默认兼容 |
| `Round 13` | `R13-B2` | 实现 `LocalModifiedPearsonIV` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `未开始` | `N/A` | `-` | 可配置、可回归、可解释 |
| `Round 13` | `R13-B3` | 实现 `LocalTabulatedSurfDistrib` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `未开始` | `N/A` | `-` | 局地制表分布可用 |
| `Round 13` | `R13-B4` | 实现 `TwoAxesTabulatedVelocitySurfDistrib` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `未开始` | `N/A` | `-` | 双轴速度制表分布可用 |
| `Round 13` | `R13-B5` | 实现 `UniformVelocitySurfDistrib` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `未开始` | `N/A` | `-` | 均匀速度分布可用 |
| `Round 13` | `R13-B6` | 实现 `FluidSurfDistrib` / `FowlerNordheimSurfDistrib` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `未开始` | `N/A` | `-` | 特化分布族可路由与回归 |
| `Round 13` | `R13-B7` | 补齐分布族单测与 smoke | `Tools/Particle/test/SurfaceDistributionFunction_test.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | `未开始` | `N/A` | `-` | 每个新增分布族具备单测证据 |
| `Round 14` | `R14-C1` | 实现 `AutomaticBarrierCurrentScaler` | `Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | `未开始` | `N/A` | `-` | 自动选择规则可复现 |
| `Round 14` | `R14-C2` | 实现 `MultipleCurrentScaler` | `Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | `未开始` | `N/A` | `-` | 组合策略可配置 |
| `Round 14` | `R14-C3` | 实现 `GlobalTempCurrentScaler` / `SmoothedGlobalTempCurrentScaler` | `Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | `未开始` | `N/A` | `-` | 温度族 scaler 可回归 |
| `Round 14` | `R14-C4` | 实现显式 `CurrentVariation` / `LocalVariation` 策略类 | `Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | `未开始` | `N/A` | `-` | 显式策略与现有线性化兼容 |
| `Round 14` | `R14-C5` | 导出自动/组合 scaler 决策元数据 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `未开始` | `N/A` | `-` | 决策链可观测可追溯 |
| `Round 14` | `R14-C6` | 补齐 scaler 单测与 smoke | `Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | `未开始` | `N/A` | `-` | 新增 scaler 用例通过 |
| `Round 15` | `R15-D1` | 实现 `ReducableDIDV` 基础框架 | `Tools/Coupling/include/SurfaceCircuitCoupling.h`、`Tools/Coupling/src/SurfaceCircuitCoupling.cpp` | `未开始` | `N/A` | `-` | 可规约组合与现有聚合共存 |
| `Round 15` | `R15-D2` | 实现 `RedDIDVFromSurfDIDV` / `RedDIDVfromRegDIDV` | `Tools/Coupling/src/SurfaceCircuitCoupling.cpp` | `未开始` | `N/A` | `-` | 规约路径可运行 |
| `Round 15` | `R15-D3` | 实现 `DIDVfromSurfDIDV` / `SurfDIDVFromMatrices` 显式路由 | `Tools/Coupling/src/SurfaceCircuitCoupling.cpp` | `未开始` | `N/A` | `-` | DIDV 构造链完整 |
| `Round 15` | `R15-D4` | 补齐 DIDV 高阶组合测试 | `Tools/Coupling/test/Coupling_test.cpp` | `未开始` | `N/A` | `-` | 高阶 DIDV 单测 + smoke 通过 |
| `Round 16` | `R16-E1` | 补齐 `ErosionParamSet` 对应参数资产 | `Tools/Material/src/SurfaceMaterialModel.cpp` | `未开始` | `N/A` | `-` | GEO/LEO 主线材料覆盖完整 |
| `Round 16` | `R16-E2` | 补齐背散射二维表（能量 x 角度）路径 | `Tools/Material/src/SurfaceMaterialModel.cpp` | `未开始` | `N/A` | `-` | 二维表驱动生效且兼容 |
| `Round 16` | `R16-E3` | 复核 FN 参数链全路径一致性 | `Tools/FieldSolver/include/SurfaceBarrierModels.h`、`Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | `未开始` | `N/A` | `-` | 参数从配置到求解一致 |
| `Round 16` | `R16-E4` | 补齐材料参数回归测试 | `Tools/Material/test/Material_test.cpp` | `未开始` | `N/A` | `-` | 材料参数测试与场景 smoke 通过 |
| `Round 17` | `R17-F1` | 建立 native volume parity 子计划与接口边界 | `Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h` | `未开始` | `N/A` | `-` | native parity 设计评审通过 |
| `Round 17` | `R17-F2` | 实现最小 `VolMesh` native 路径 | `Tools/Mesh/*` | `未开始` | `N/A` | `-` | 原生体网格最小能力可运行 |
| `Round 17` | `R17-F3` | 实现最小 `VolField` native 路径 | `Tools/FieldSolver/*` | `未开始` | `N/A` | `-` | 原生体场步进可运行 |
| `Round 17` | `R17-F4` | 实现最小 `VolDistrib` / `VolInteract` native 路径 | `Tools/Particle/*` | `未开始` | `N/A` | `-` | 体分布与体交互可运行 |
| `Round 17` | `R17-F5` | 完成 bridge 与 native 并存兼容 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `未开始` | `N/A` | `-` | 双路径并存且不回归 |
| `Round 18` | `R18-G1` | 执行关键单测集合封板 | `build/*` | `未开始` | `N/A` | `-` | 关键测试集全绿 |
| `Round 18` | `R18-G2` | 执行子集矩阵 gate 封板 | `scripts/run/surface_reference_matrix_round11_subset.json` | `未开始` | `N/A` | `-` | 子集 gate `PASS` |
| `Round 18` | `R18-G3` | 执行全量矩阵 gate 封板 | `scripts/run/surface_reference_matrix.json` | `未开始` | `N/A` | `-` | 全量 gate `PASS` |
| `Round 18` | `R18-G4` | 执行 contract / schema 校验 | `documentation/contracts/config_schema_v1.json` | `未开始` | `N/A` | `-` | schema 与输出契约一致 |
| `Round 18` | `R18-G5` | 文档闭环与状态回填封板 | `Toolkit/Surface Charging/docs/surface_num_alignment_refactor_plan.md`、`Toolkit/Surface Charging/docs/surface_round10_equivalence_matrix.md` | `未开始` | `N/A` | `-` | 任务状态、证据、完成时间全部回填 |

## 5. 状态跟踪与更新规范

### 5.1 必填字段

- 每次更新任务状态时，至少补齐：
  - `状态`
  - `最后更新`
  - `完成时间`（未完成可填 `N/A`）
  - `证据`
  - `下一步`

### 5.2 更新规则

- 若任务刚开始摸底，状态从 `未开始` 改为 `已盘点`。
- 若正在明确接口、模块边界、拆分方案，状态改为 `设计中`。
- 若已经开始代码改动或测试改动，状态改为 `实现中`。
- 若代码、测试、文档、验收均收口，状态改为 `已完成`。
- 若状态改为 `已完成`，必须同时填写 `完成时间` 与最小证据集合（代码 + 测试 + gate/contract）。
- 若该轮次已有稳定结论但不再属于当前正文主线，状态改为 `历史归档`。
- 若确认本轮不做但结论已经稳定，状态改为 `已冻结`。
- 若受到依赖阻塞或存在关键分歧，状态改为 `受阻`，并写明原因。

### 5.3 证据记录格式

- 证据优先写成以下形式之一：
  - 代码路径，例如 `Tools/Coupling/include/SurfaceCircuitCoupling.h`
  - 测试路径，例如 `Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`
  - gate 或脚本，例如 `scripts/run/surface_reference_matrix.json`
  - 文档路径，例如 `documentation/contracts/**`
- 若任务已完成，建议补：
  - 影响文件
  - 测试名称
  - 对应 artifact 或 gate

### 5.4 更新示例

| 轮次 | 任务ID | 任务 | 状态 | 最后更新 | 完成时间 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `Round 9` | `R9-A` | 冻结 `RuntimePlan` 边界与字段清单 | `已完成` | `2026-04-14` | `2026-04-14` | `Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/docs/architecture.md` | 进入 `R9-B` / `R9-C` |

## 6. 验收标准

### 6.1 架构验收

- 六个 SPIS 域都必须包含：
  - SPIS 模块职责
  - 当前项目落点
  - 现状判断
  - 核心缺口
  - 改进方向
- 不允许把 surface bridge 契约冻结误写成 native volume kernel 对标完成。
- 不允许把 facade 落地误写成 runtime 编排拆分完成。

### 6.2 计划验收

- 每一轮都必须具备：
  - 目标
  - 输入
  - 实施项
  - 输出
  - 依赖
  - 验收
  - 当前状态
  - 最后更新
  - 证据
- 当前正文主线必须始终清晰可见；历史归档不得与活跃轮次混写。

### 6.3 执行验收

- 任意一轮结束后，本文档都能准确显示：
  - 哪些任务已完成
  - 哪些任务进行中
  - 哪些任务被冻结
  - 哪些任务受阻
  - 哪些任务已经转入历史归档

### 6.4 兼容验收

- 必须维持以下边界可持续兼容：
  - `SurfaceChargingConfig`
  - preset 名称与别名
  - `SurfaceRuntimeRoute`
  - `SurfacePicStrategy`
  - 主 CSV 输出
  - 现有 sidecar 合同

### 6.5 回归验收

- Surface 主线继续以以下路径为主 gate：
  - `ref/GEO/**`
  - `ref/LEO/**`
  - `ref/mat/**`
  - `scripts/run/surface_reference_matrix.json`
- `Round 9` 必须额外补齐：
  - `SurfaceScenarioCatalogTest`
  - `SurfaceSimulationRunnerTest`
  - `SurfaceChargingSmokeTest` 中与 route/legacy/shared-coupled 直接相关的用例
  - 需要时增加新的 object-layer / transition-engine 单测入口
- `Round 10` 必须额外补齐：
  - surface 功能能力矩阵
  - 运行模式与输出矩阵
  - GEO / LEO / material 参考场景等价矩阵
  - 未一致项的明确 backlog
- `Round 11` 必须额外补齐：
  - Surf material/distrib/scaler 新增族的单测覆盖
  - Circ `PWL/EXP` 波形与 weighted/multi-surface DIDV 组合测试
  - Top 事件级 transition（`LocalTime/Spinning/SunFlux` 最小集）对象层与 smoke 验证
  - reference matrix 子集门禁与等价矩阵状态刷新

## 7. 附录：当前已完成的 Surface 基础对齐成果摘要

### 7.1 已沉淀到 `Tools/*` 的核心表面能力

- `Tools/Material`
  - `SurfaceMaterialModel`
  - `BasicSurfaceMaterialModel`
  - 表面材料响应与 legacy yield/helper 的公共化接口
- `Tools/Particle`
  - `SurfaceDistributionFunction`
  - `IsotropicMaxwellianFluxDistribution`
  - `BiMaxwellianFluxDistribution`
  - `TabulatedSurfaceDistributionFunction`
  - `SurfaceFluxSampler`
  - `SurfaceEmissionSampler`
- `Tools/FieldSolver`
  - `SurfaceBarrierModels`
  - `BarrierCurrentScaler`
  - `SecondaryRecollectionScaler`
  - `VariableBarrierScaler`
- `Tools/Coupling`
  - `SurfaceCircuitCoupling`
  - `SurfaceCurrentLinearizationKernel`
- `Tools/Solver`
  - `SurfaceSolverFacade`

### 7.2 已建立的 Surface 主线运行时与回归基础

- `SCDATUnified`、`SurfacePic`、`SurfacePicHybrid` 已共享统一 surface kernel 路径。
- `LegacyBenchmark` 已被限定为回归与参考 reproduction 边界。
- 现有 Surface 域已具备：
  - benchmark artifact
  - graph / monitor / field / volume sidecar
  - machine-readable contract
  - reference matrix gate
- `SurfaceScenarioCatalog`、`SurfaceScenarioLoader`、`SurfaceSimulationRunner` 已完成首批对象层入口降耦。

### 7.3 当前仍需持续关注的巨型热点

- `Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`
- `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`

这两个文件当前仍是后续边界收口与职责拆分的首要热点。
