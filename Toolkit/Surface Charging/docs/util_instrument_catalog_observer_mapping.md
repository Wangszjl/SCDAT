# `Util/Instrument` catalog / observer-route 映射

## 1. 目的

本文用于落实 `Round 51 / Wave AF`，把 `ref/SPIS/num-master/src/main/java/spis/Util/Instrument/**` 下的 `51` 个 SPIS class，逐个映射到当前项目中的 instrument catalog、observer route、monitor sidecar、probe/current API 与 diagnostics 承载层。

本文回答的不是 “`util_instrument` family 是否已经封板”，而是：

- 哪些 instrument class 在当前项目里已经有明确的运行态行为或导出物证；
- 哪些 class 仅保留为 instrument metadata / catalog / route 语义；
- 哪些 class 目前还没有独立可观测物证，只能视作聚合承载下的近似解释。

> 口径说明：`util_instrument` family 已在 family 级 gate / reference matrix / sealoff 中封板。  
> 本文只处理 class-level 证据细化，不改变现有 `completed` 结论。

## 2. 当前承载骨架

| 当前承载层 | 角色 | 主要文件 |
| --- | --- | --- |
| instrument set 路由层 | 解析 `surface_instrument_set`、把观测器集合并入 preset / config / runtime plan | `Main/main.cpp`、`Main/SurfaceScenarioLoader.cpp`、`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp` |
| runtime metadata 层 | 导出 `surface_instrument_set`、`surface_instrument_contract_id`、observer artifact 路径与 transition observer 最新状态 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |
| observer / monitor artifact 层 | 写出 `.metadata.json`、`.monitor.json`、`.shared_runtime_observer.json` | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |
| probe / current runtime 层 | 提供 live probe、电流响应、Langmuir probe transition 与电流导出 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/src/PicMccSurfaceCurrentSampler.cpp` |
| diagnostics 基础层 | 提供通用 `Monitor` / `MonitorManager` 承载监测器树 | `Tools/Diagnostics/include/Monitor.h`、`Tools/Diagnostics/include/MonitorManager.h` |
| 测试证据层 | 对 instrument set、observer sidecar、Langmuir probe 与 shared runtime observer 做 smoke / object-layer 断言 | `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp` |

## 3. 判定规则

- `已有 runtime 行为`
  - 当前 class 能回链到显式运行态字段、artifact、probe API 或 smoke 断言；
  - 不要求存在与 SPIS 同名的独立 C++ class，但要有直接物证。
- `仅 metadata route`
  - 当前保留了 catalog / factory / instrument set / route / sidecar 层语义；
  - 但没有单独的、可点名复核的运行态实体。
- `尚未单独可观测`
  - family 语义仍在，但当前只能以更大粒度的 observer / monitor / distribution / detector 聚合层解释；
  - 还没有单独 artifact 或 smoke 专门证明这一类对象。

## 4. 逐 class 映射

| SPIS class | 当前主要承载落点 | 结论 | 说明 |
| --- | --- | --- | --- |
| `AbstractParticleDetector` | `computeSurfaceCurrents()`、shared runtime observer、surface node histories | `尚未单独可观测` | 当前没有抽象粒子探测器 class 树，探测语义被电流/粒子响应聚合吸收。 |
| `AvSCPotentialMonitor` | `DensePlasmaSurfaceCharging` 电位历史、metadata sidecar、monitor json | `已有 runtime 行为` | 航天器平均电位在运行时历史和导出 sidecar 中可观测。 |
| `BundleAPI_DefaultInstruments` | `SurfaceInstrumentSetKind`、preset catalog、`surface_instrument_set` route | `仅 metadata route` | 默认仪器集合由 instrument set 枚举和 preset 吸收。 |
| `BundleAPI_InstrumentsCatalogue` | `parseSurfaceInstrumentSetKind(...)`、scenario loader、preset catalog | `仅 metadata route` | 当前由 CLI/JSON 入口和 catalog 路由承担。 |
| `CreateDefaultMonitoringInstruments` | `SurfaceChargingCases.cpp` 默认 preset、observer set 默认值 | `仅 metadata route` | 默认监测器集合创建被 preset 构造过程折叠承载。 |
| `DensityPS` | observer artifact 与分布/粒子导出聚合层 | `尚未单独可观测` | 没有独立 density instrument artifact。 |
| `DistFuncOfOneVariable` | 分布函数输出与 sampling helper 聚合层 | `尚未单独可观测` | 当前没有专门的一元分布仪器对象。 |
| `ESNAbstractMonitor` | `Monitor` / `MonitorManager`、surface observer metadata | `仅 metadata route` | 抽象监测器层已由通用 diagnostics 承担。 |
| `ESNBulkConductivityMonitor` | transition conductivity metadata、monitor json | `已有 runtime 行为` | 电导演化相关观测已进入 sidecar 和 transition observer 状态。 |
| `ESNCurrentMonitor` | surface current histories、`computeSurfaceCurrents()`、CSV/metadata 导出 | `已有 runtime 行为` | 电流监测语义有直接运行态输出。 |
| `ESNIndividualCurrentMonitor` | surface node current histories 与 node-level series 导出 | `已有 runtime 行为` | patch/node 级电流已可在运行导出中点名追踪。 |
| `ESNPotentialDifferenceMonitor` | patch/body/reference potential histories、shared runtime observer | `已有 runtime 行为` | 电位差可通过 current probe 与 observer artifact 间接复核。 |
| `ESNPotentialMonitor` | potential histories、metadata sidecar、monitor json | `已有 runtime 行为` | 电位观测在主导出与 sidecar 中有直接物证。 |
| `ESNPotentialSCFrameMonitor` | SC frame potential metadata 与 observer artifact | `已有 runtime 行为` | 航天器参考框架下的电位观测已聚合到 sidecar。 |
| `ESNPotentialVariationMonitor` | potential histories、transition observer、monitor json | `已有 runtime 行为` | 电位变化可由历史序列和 observer 复核。 |
| `EnergyDistFuncPS` | distribution/sampler 导出聚合层 | `尚未单独可观测` | 当前没有单独 energy distribution instrument artifact。 |
| `Instrument` | `surface_instrument_set`、`surface_instrument_contract_id` | `仅 metadata route` | 基础 instrument 概念由 instrument set 和 contract id 承担。 |
| `InstrumentSerial` | metadata/contract 序列化路径 | `仅 metadata route` | 序列化语义被 sidecar/JSON 导出统一吸收。 |
| `InstrumentsCatalogue` | `SurfaceScenarioCatalog`、`parseSurfaceInstrumentSetKind(...)` | `仅 metadata route` | 仪器目录概念由 preset catalog 与 parser 承担。 |
| `InteractorInstrumentContainer` | material/interactor route + instrument set metadata | `仅 metadata route` | interactor 与 instrument 的组合语义仍停留在路由层。 |
| `KineticEnergyPS` | 粒子/分布统计导出聚合层 | `尚未单独可观测` | 当前没有专门 kinetic-energy post-sensor artifact。 |
| `LangmuirProbe` | `transition_langmuir_probe_*`、`live_pic_probe_delta_v`、probe current API、smoke 断言 | `已有 runtime 行为` | Langmuir probe 是当前最明确的 instrument runtime 之一。 |
| `LinePS` | 空间采样形状由 observer/export 聚合层吸收 | `尚未单独可观测` | 线型 post-sensor 当前没有独立 artifact。 |
| `MonitorInstrument` | `Monitor` / `MonitorManager`、`.monitor.json` artifact | `已有 runtime 行为` | 通用监视型 instrument 已有显式 diagnostics 落点。 |
| `MultiPlasmaSensor` | instrument set + plasma diagnostics route | `仅 metadata route` | 多传感器组合当前主要体现在 route/config 层。 |
| `MultipleESNBulkConductivityMonitor` | conductivity observer metadata、monitor json | `已有 runtime 行为` | 多通道电导监测被聚合到 observer/monitor artifact。 |
| `MultipleESNPotentialSensors` | node/patch potential histories、shared observer json | `已有 runtime 行为` | 多电位传感器语义可由多节点电位序列复核。 |
| `MultipleESNPotentialSensorsSCFrame` | SC frame observer metadata 与 shared runtime observer | `已有 runtime 行为` | SC frame 多电位观测已被 artifact 聚合承载。 |
| `MultipleESNTotalCurrentMonitor` | total current histories、shared observer json、CSV 导出 | `已有 runtime 行为` | 多通道总电流监测具备直接导出物证。 |
| `NumericsMonitor` | `.monitor.json`、数值状态/graph consistency metadata | `已有 runtime 行为` | numerics/consistency/graph coupling 指标有专门 monitor artifact。 |
| `ParticleDetector` | probe/current API 与粒子响应聚合层 | `尚未单独可观测` | 没有单独的 ParticleDetector C++ 类型或 artifact。 |
| `ParticleDetectorInterface` | detector 概念由 probe API 与采样器吸收 | `尚未单独可观测` | 接口树未在当前 C++ 项目保留。 |
| `ParticleSPS` | 粒子统计后处理导出聚合层 | `尚未单独可观测` | 当前没有专门 particle post-sensor artifact。 |
| `PlasmaSensor` | diagnostic set / instrument set route | `仅 metadata route` | plasma sensor 概念主要由 route 与 preset 维护。 |
| `PointPS` | 空间采样形状聚合层 | `尚未单独可观测` | 点采样后处理没有单独 artifact。 |
| `PotentialLPS` | potential histories + probe API 聚合层 | `尚未单独可观测` | 当前没有单独 `PotentialLPS` 对象或专门 smoke。 |
| `PotentialPS` | potential histories + metadata sidecar 聚合层 | `尚未单独可观测` | 电位 post-sensor 仍通过更大粒度 artifact 承载。 |
| `ScalarPlasmaSensor` | plasma diagnostics route、instrument set metadata | `仅 metadata route` | 标量等离子体传感器未作为独立 runtime class 保留。 |
| `SimulationControl` | `SurfaceSimulationRunner`、transition observer、runtime plan | `已有 runtime 行为` | 运行控制语义在 runner/transition observer 中有直接行为。 |
| `SphericalPS` | 空间采样形状聚合层 | `尚未单独可观测` | 球形 post-sensor 当前未形成专门导出物证。 |
| `SurfaceFluxDistFunctionPS` | surface flux sampling 与 distribution helper 聚合层 | `尚未单独可观测` | surface flux distribution 语义仍停留在更底层分布/采样模块。 |
| `ThreadInstrument` | 顺序 runner + observer checkpoints | `仅 metadata route` | 当前没有线程级 instrument 对象，只有 runner/observer 协调。 |
| `TopInstrumentFactory` | CLI/loader/parser/preset catalog | `仅 metadata route` | factory 语义由 parser 和 preset 构造承担。 |
| `TotalCurrentOnSC` | total current histories、CSV 导出、`computeSurfaceCurrents()` | `已有 runtime 行为` | 航天器总电流在当前实现中有直接运行态输出。 |
| `TotalEnergyDistFuncPS` | energy distribution 聚合导出层 | `尚未单独可观测` | 当前没有专门总能量分布 instrument artifact。 |
| `TotalEnergyPS` | 能量统计聚合导出层 | `尚未单独可观测` | 总能量 post-sensor 仍无独立可观测对象。 |
| `TotalSuperParticlePS` | 粒子统计聚合导出层 | `尚未单独可观测` | 超粒子数统计没有单独 artifact。 |
| `VelocityDistFunctionPS` | velocity distribution/sampler 聚合层 | `尚未单独可观测` | 速度分布仍由 sampler/distribution 家族聚合承载。 |
| `VirtualInstrument` | metadata-only instrument set、sidecar contract | `仅 metadata route` | 虚拟仪器概念由 metadata-only 路由承担。 |
| `VirtualParticleDetector` | probe/current API 聚合层 | `尚未单独可观测` | 虚拟粒子探测器没有专门 artifact。 |
| `VolDistribMomentPS` | volume distribution / moment 聚合导出层 | `尚未单独可观测` | 体分布矩统计仍停留在更大粒度导出层。 |

## 5. 当前最重要的可观测证据

### 5.1 instrument route

- `Main/main.cpp`
  - `parseSurfaceInstrumentSetKind(...)`
  - 支持 `metadata_only` 与 `surface_pic_observer_set`
- `Main/SurfaceScenarioLoader.cpp`
  - 解析 `surface_instrument_set`
- `Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`
  - 将 instrument set 编译进 runtime plan

### 5.2 metadata / observer artifact

- `DensePlasmaSurfaceCharging.cpp`
  - 导出 `surface_instrument_set`
  - 导出 `surface_instrument_contract_id = surface-instrument-observer-v1`
  - 导出 `surface_pic_runtime_boundary_observer_contract_id = surface-pic-boundary-observer-v1`
  - 写出 `.monitor.json`
  - 写出 `.shared_runtime_observer.json`

### 5.3 runtime behavior

- `DensePlasmaSurfaceCharging.cpp`
  - `computeSurfaceCurrents()` 提供 probe/current 直接查询
  - `transition_langmuir_probe_derivative_scale`
  - `transition_observer_*` 一组运行态字段
- `SurfaceTransitionEngine.cpp`
  - `LangmuirProbeTransition`
  - `transition_observer_enabled`
  - observer checkpoint 与 active transition 计算

### 5.4 smoke / object-layer tests

- `SurfaceCharging_smoke_test.cpp`
  - `SurfaceReferenceAndRuntimeContractsExportMetadata`
  - `surface_instrument_set` 与 `surface_instrument_contract_id` 断言
  - `surface_transition_observer_*` 断言
  - `Langmuir probe` 相关断言
  - `.shared_runtime_observer.json` 存在性与内容断言
- `SurfaceCharging_object_layer_test.cpp`
  - `transition.observer.active`
  - `transition.observer.checkpoint_dt_s`
  - `transition.extended_events.langmuir_probe_transition_active`

## 6. 汇总结论

- `Util/Instrument` 当前没有 family 级功能缺口，核心差异仍是 class-level 对象层没有保留 SPIS 的完整 class 树。
- 最稳定的一组 instrument runtime 语义已经比较清楚：
  - `LangmuirProbe`
  - potential/current/observer monitors
  - `NumericsMonitor`
  - `TotalCurrentOnSC`
- 最大的剩余弱项不是“没有 instrument 能力”，而是大量 `PS` / `ParticleDetector` / distribution 型 instrument 还缺独立 artifact 或单独 smoke。

## 7. 后续建议

如果后面还要继续增强 `Util/Instrument`，优先顺序建议是：

1. 为 `PotentialPS / PotentialLPS / PointPS / LinePS / SphericalPS` 补一个空间采样 artifact 口径。
2. 为 `ParticleDetector / VirtualParticleDetector / ParticleSPS` 补一个粒子探测输出索引。
3. 为 `EnergyDistFuncPS / VelocityDistFunctionPS / SurfaceFluxDistFunctionPS / VolDistribMomentPS` 补一个 distribution-observer 对照页。

## 8. 证据索引

- route / config
  - `Main/main.cpp`
  - `Main/SurfaceScenarioLoader.cpp`
  - `Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`
  - `Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`
- runtime / artifact
  - `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`
  - `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`
  - `Toolkit/Surface Charging/include/PicMccSurfaceCurrentSampler.h`
  - `Toolkit/Surface Charging/src/PicMccSurfaceCurrentSampler.cpp`
  - `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`
  - `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`
- diagnostics
  - `Tools/Diagnostics/include/Monitor.h`
  - `Tools/Diagnostics/include/MonitorManager.h`
- tests / audit
  - `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`
  - `Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`
  - `Toolkit/Surface Charging/docs/surface_spis_class_audit.json`
  - `Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md`
  - `Toolkit/Surface Charging/docs/surface_num_alignment_refactor_plan.md`
