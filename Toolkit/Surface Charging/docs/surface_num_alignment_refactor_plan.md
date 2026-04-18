# Surface Charging 与 SPIS 功能 1:1 对齐总计划

## 1. 目标与使用规则

### 1.1 文档定位

- 本文档是 `Toolkit/Surface Charging` 对齐 `ref/SPIS/num-master` 的唯一长期维护主计划。
- 本文档既是差异盘点文档，也是执行看板；后续每一轮推进结束后，都必须直接在本文档中回填状态、证据和下一步。
- 本文档当前只服务 `Surface Charging` 主线及其直接依赖模块，不负责整个仓库的总规划。
- 本文档当前覆盖的代码范围限定为：
  - `Toolkit/Surface Charging`
  - `Main`
  - `Tools/Basic`
  - `Tools/Geometry`
  - `Tools/Diagnostics`
  - `Tools/Output`
  - `Tools/Boundary`
  - `Tools/Material`
  - `Tools/Mesh`
  - `Tools/Particle`
  - `Tools/FieldSolver`
  - `Tools/Coupling`
  - `Tools/Solver`
  - `Tools/PICcore`
  - 与上述模块直接关联的 `scripts/` 与 `documentation/contracts/`
- `Toolkit/Internal Charging`、`Toolkit/Radiation`、`Toolkit/Vacuum Arc` 等模块只在需要说明复用边界时被引用，不进入本文档正文主表。

### 1.2 权威来源

- 功能与结构对标唯一基线：
  - `ref/SPIS/num-master/src/main/java/spis/Surf/**`
  - `ref/SPIS/num-master/src/main/java/spis/Circ/**`
  - `ref/SPIS/num-master/src/main/java/spis/Solver/**`
  - `ref/SPIS/num-master/src/main/java/spis/Vol/**`
  - `ref/SPIS/num-master/src/main/java/spis/Util/**`
  - `ref/SPIS/num-master/src/main/java/spis/Top/**`
- 回归与验收参考：
  - `ref/GEO/**`
  - `ref/LEO/**`
  - `ref/mat/**`
- 当前项目事实来源：
  - `Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`
  - `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`
  - `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`
  - `Tools/Material/include/SurfaceInteraction.h`
  - `Tools/Material/include/SurfaceMaterialModel.h`
  - `Tools/Particle/include/SurfaceDistributionFunction.h`
  - `Tools/FieldSolver/include/SurfaceBarrierModels.h`
  - `Tools/FieldSolver/include/SurfaceFieldVolumeBridge.h`
  - `Tools/Coupling/include/SurfaceCircuitCoupling.h`
  - `Tools/Solver/include/SurfaceSolverFacade.h`
  - `scripts/run/surface_reference_matrix.json`
  - `documentation/contracts/**`

### 1.3 1:1 对齐判定规则

- 本文档中的“1:1 对齐”专指与 `ref/SPIS/num-master` 在 `Surface Charging` 相关功能上的一一映射，不再使用“最小集可用”代替“功能对齐”。
- 只有同时满足以下条件，某个 SPIS family 才能标记为 `已等价实现`：
  - 在当前项目中有明确代码落点
  - 有可达的配置或运行路由
  - 运行时语义与 SPIS 对应 family 一致或已明确建立适配关系
  - 有单测、smoke、gate 或参考矩阵证据
- 若当前项目仅有近似实现、最小闭环或统一抽象承载，必须标记为 `近似实现`，不得写成 `已对齐`。
- 若当前项目没有实现，必须标记为 `尚无实现` 并进入后续 backlog。
- 若某项经确认不属于本文档边界，必须标记为 `不纳入当前目标范围`，并给出理由与证据。

### 1.4 使用规则

- 每一轮执行结束后，必须更新：
  - `状态`
  - `最后更新`
  - `证据`
  - `下一步`
- 任何任务不得只写大方向，必须能在一轮执行后明确改一次状态。
- 任何被冻结、受阻或转入历史归档的任务都必须写明原因。
- 本文档允许保留历史轮次记录，但历史事实与当前 1:1 对齐结论必须分开书写，避免把“已完成某轮任务”等价写成“已与 SPIS 对齐”。

### 1.5 外部兼容边界

- 默认保持以下外部接口与契约持续兼容：
  - `SurfaceChargingConfig`
  - `SurfaceChargingScenarioPreset`
  - `SurfaceRuntimeRoute`
  - `SurfacePicStrategy`
  - `LegacyBenchmarkExecutionMode`
  - 现有主 CSV 时序输出
  - 现有 benchmark / graph / monitor / field / volume / mapping / bridge sidecar 家族
- 允许扩展内部 family、对象层、元数据与 sidecar 字段。
- 不允许删除已有消费者依赖的主字段或 route/preset 名称，除非文档先明确批准兼容性破坏。

## 2. Surface 主线六域 1:1 差距分析

### 2.1 总览

| 域 | SPIS 主要落点 | 当前项目主要落点 | 当前状态 | 当前关键差异 |
| --- | --- | --- | --- | --- |
| `Surf` | `spis/Surf/*` | `Tools/Material`、`Tools/Particle`、`Tools/FieldSolver`、`Tools/Coupling`、`Toolkit/Surface Charging` | `已完成` | `Surf/Material`、`SurfInteract`、`SurfDistrib`、`SurfMesh`、`SurfField` 已全部具备 direct route / config / test / matrix 证据，转入长期回归维护 |
| `Circ` | `spis/Circ/*` | `Tools/Coupling`、`Toolkit/Surface Charging` | `已完成` | `Circ/Circ`、`CircField`、`Circ/DIDV` 已补齐 direct route / smoke / matrix 证据，转入长期回归维护 |
| `Solver` | `spis/Solver/*` | `Tools/Solver`、`Tools/FieldSolver`、`Tools/Particle`、`Tools/PICcore` | `已完成` | `Circuit`、`Util`、`ElectroMag`、`Matter` 已完成 direct-SPIS semantic evidence / route / test / matrix 证据闭环，转入长期回归维护 |
| `Vol` | `spis/Vol/*` | `Tools/Boundary`、`Tools/Mesh`、`Tools/FieldSolver`、`Tools/Particle`、`Toolkit/Surface Charging` bridge/runtime | `已完成` | `BC` 抽象长尾、`Field`、`Distrib`、`Interact` 四组 family 已具备显式实现与 direct evidence；`Geom`、`VolMesh` 已升级为独立 family 级 contract / gate |
| `Util` | `spis/Util/*` | `Tools/Basic`、`Tools/Geometry`、`Tools/Diagnostics`、`Tools/Output`、`Tools/Particle`、`Tools/Solver`、`Tools/Mesh` | `已完成` | `Util/**` 已以 `util_full_stack` 形式纳入 target families，并新增 `util_distrib_func`、`util_exception`、`util_instrument`、`util_func`、`util_list`、`util_octree`、`util_octree_mesh`、`util_part`、`util_phys`、`util_sampler`、`util_monitor`、`util_io`、`util_table`、`util_matrix`、`util_vect` 十五个子包级 direct code / test / config-route / reference-matrix 证据 |
| `Top` | `spis/Top/*` | `SurfaceScenarioCatalog`、`SurfaceScenarioLoader`、`SurfaceSimulationRunner`、`SurfaceRuntimePlan`、`SurfaceTransitionEngine`、`DensePlasmaSurfaceCharging`、`Main/*` | `已完成` | 除 `Transition` 外，`Top/Top`、`Simulation`、`Plasma`、`SC`、`Grid`、`Default` 已纳入逐族盘点与独立 gate |

### 2.2 Surf 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | `Material`、`SurfInteract`、`SurfField`、`SurfDistrib` 与 `Util/DistribFunc` 共同构成表面材料响应、入射/发射分布、barrier/recollection、DIDV 前置物理内核 |
| 当前项目落点 | `Tools/Material`、`Tools/Particle`、`Tools/FieldSolver`、`Tools/Coupling` 已沉淀核心能力；`Toolkit/Surface Charging` 负责场景装配与运行时 glue |
| 当前状态 | `已完成` |
| 关键差异 | 当前无阻塞性缺口；`Surf/Material`、`SurfInteract`、`SurfDistrib`、`SurfMesh`、`SurfField` 已全部具备 direct route / config / test / matrix 证据 |
| 对齐方向 | 转入长期回归维护，保留 family-to-code / test / config-route 与 reference matrix gate |

### 2.3 Circ 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | `Circ/Circ`、`Circ/DIDV`、`CircField` 定义电路元件、电流分布、DIDV 组合与电路求解前置数据 |
| 当前项目落点 | `Tools/Coupling` 已有通用节点/支路/动态激励/DIDV 组合接口；`Toolkit/Surface Charging` 消费这些接口 |
| 当前状态 | `已完成` |
| 关键差异 | 当前无阻塞性缺口；`Circ/Circ`、`CircField`、`Circ/DIDV` 已完成 direct route、运行态断言与 matrix 证据闭环 |
| 对齐方向 | 转入长期回归维护，保留 circuit / didv 组合与 multi-patch 场景回归 gate |

### 2.4 Solver 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | `Circuit`、`ElectroMag`、`Matter` 负责电路、泊松/电磁、粒子推进等求解任务 |
| 当前项目落点 | `Tools/Solver`、`Tools/FieldSolver`、`Tools/Particle`、`Tools/PICcore` |
| 当前状态 | `已完成` |
| 关键差异 | solver facade 仍承担 runtime 编排，但已不再构成 1:1 缺口；`Circuit`、`Util`、`ElectroMag`、`Matter` family 均已具备 direct-SPIS semantic evidence |
| 对齐方向 | 转入长期回归维护，保留 solver family metadata 与 reference matrix gate |

### 2.5 Vol 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | `BC`、`Geom`、`VolMesh`、`VolDistrib`、`VolField`、`VolInteract` 共同组成原生体域数值内核 |
| 当前项目落点 | `Tools/Boundary`、`Tools/Mesh`、`Tools/FieldSolver`、`Tools/Particle` 的一部分能力 + `Toolkit/Surface Charging` 中的 field/volume bridge、pseudo-volume runtime |
| 当前状态 | `已完成` |
| 关键差异 | 当前无阻塞性缺口；`BC` 抽象长尾、`Field`、`Distrib`、`Interact` 四组 family 已完成显式 route、单测、smoke 与 matrix case 闭环；`Geom`、`VolMesh` 已补独立 metadata / contract / gate |
| 对齐方向 | 转入长期回归维护；bridge 继续保留为兼容路径，但不再作为独立 remediation 目标 |

### 2.6 Util 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | 分布函数、数学函数、积分、异常、参数/IO 工具等基础设施 |
| 当前项目落点 | `Tools/Basic`、`Tools/Geometry`、`Tools/Diagnostics`、`Tools/Output`、`Tools/Mesh`，以及 `Tools/Particle`、`Tools/Solver` 中的 surface 公共 helper |
| 当前状态 | `已完成` |
| 关键差异 | 当前无阻塞性缺口；`Util/**` 已按 `util_full_stack` 纳入 target families，并新增 `surface_util_distrib_func_*`、`surface_util_exception_*`、`surface_util_instrument_*`、`surface_util_func_*`、`surface_util_list_*`、`surface_util_octree_*`、`surface_util_octree_mesh_*`、`surface_util_part_*`、`surface_util_phys_*`、`surface_util_sampler_*`、`surface_util_monitor_*`、`surface_util_io_*`、`surface_util_table_*`、`surface_util_matrix_*`、`surface_util_vect_*` metadata，形成四批子包级 gate |
| 对齐方向 | 转入长期回归维护，保留 family-to-code / test / config-route 与 reference matrix 证据链 |

### 2.7 Top 域

| 项目 | 内容 |
| --- | --- |
| SPIS 模块职责 | `Simulation`、`Scenario`、`Transition`、`Plasma`、`SC` 等统一管理场景、运行时生命周期与对象关系 |
| 当前项目落点 | `SurfaceScenarioCatalog`、`SurfaceScenarioLoader`、`SurfaceSimulationRunner`、`SurfaceRuntimePlan`、`SurfaceTransitionEngine`、`DensePlasmaSurfaceCharging`、`Main/main.cpp` |
| 当前状态 | `已完成` |
| 关键差异 | 当前无阻塞性缺口；除 `Transition` 外，`Top/Top`、`Simulation`、`Plasma`、`SC`、`Grid`、`Default` 六组 family 已具备 route evidence 与独立 metadata / gate |
| 对齐方向 | 转入长期回归维护，保留 `top_*` family metadata 与 route smoke / reference matrix gate |

### 2.8 当前审计结论（按 1:1 对齐口径）

| 主题 | 当前事实 | 1:1 对齐判断 | 正文处理 |
| --- | --- | --- | --- |
| `Surf` | `Surf/Material`、`Surf/SurfInteract`、`Surf/SurfDistrib`、`Surf/SurfMesh`、`Surf/SurfField` 已全部具备 direct route / config / test / matrix 证据 | `已完成` | 转入长期回归维护 |
| `Circ` | `Circ/Circ`、`CircField`、`Circ/DIDV` 的 direct route、运行态断言与 matrix 证据已补齐 | `已完成` | 转入长期回归维护 |
| `Solver` | `Circuit`、`Util`、`ElectroMag`、`Matter` family 的 direct-SPIS semantic evidence 已补齐 | `已完成` | 转入长期回归维护 |
| `Vol` | `BC`、`Field`、`Distrib`、`Interact` 四组 family 已全部闭环，`vol_full_stack` 已升级为 `completed`；`vol_bc_abstract`、`vol_field_abstract`、`vol_distrib_long_tail`、`vol_interact_long_tail`、`vol_geom`、`vol_mesh` 已具备独立 contract / gate | `已完成` | 转入长期回归维护 |
| `Util` | `Util/**` 子包级拆分已完成；`util_distrib_func`、`util_exception`、`util_instrument`、`util_func`、`util_list`、`util_octree`、`util_octree_mesh`、`util_part`、`util_phys`、`util_sampler`、`util_monitor`、`util_io`、`util_table`、`util_matrix`、`util_vect` 均已升级为独立 target family，且全部为 `completed` | `已完成` | 转入长期回归维护 |
| `Top` | `top_transition` 之外，`top_top`、`top_simulation`、`top_plasma`、`top_sc`、`top_grid`、`top_default` 已具备独立 route evidence 与 gate | `已完成` | 转入长期回归维护 |

### 2.9 SPIS -> SCDAT 全量映射总表（一级总表）

| 域 | SPIS 源路径 | 当前 SCDAT 落点 | 当前状态 | 主要缺口 | 后续任务组 |
| --- | --- | --- | --- | --- | --- |
| `Surf/Material` | `ref/SPIS/num-master/src/main/java/spis/Surf/Material/**` | `Tools/Material/include/SurfaceMaterialModel.h`、`Tools/Material/src/SurfaceMaterialModel.cpp` | `已完成` | 当前无阻塞性缺口；仅保留参数资产与 reference matrix 回归维护 | `R32-F1（已完成）` |
| `Surf/SurfInteract` | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfInteract/**` | `Tools/Material/include/SurfaceInteraction.h`、`Tools/Material/src/SurfaceInteraction.cpp` | `已完成` | 已补齐显式 family 入口、active signature 与 direct evidence；转入回归维护 | `R27-A（已完成）` |
| `Surf/SurfDistrib` | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfDistrib/**` | `Tools/Particle/include/SurfaceDistributionFunction.h`、`Tools/Particle/src/SurfaceDistributionFunction.cpp` | `已完成（Surface 范围）` | 全家族已接入主线路由、sidecar 与 matrix gate；后续仅保留回归维护 | `R22-K（已完成）` |
| `Surf/SurfMesh` | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfMesh/**` | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`、`Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Main/*` | `已完成` | `SurfMesh`、`ThreeDUnstructSurfMesh` 已形成独立 target family、structured-topology route、metadata 与 matrix 证据链；转入回归维护 | `R39-M（已完成）` |
| `Surf/SurfField` | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfField/**` | `Tools/FieldSolver/include/SurfaceBarrierModels.h`、`Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | `已完成` | 已补齐 scaler family 的 direct route / config / smoke / matrix evidence；转入回归维护 | `R31-E1（已完成）` |
| `Util/Func` | `ref/SPIS/num-master/src/main/java/spis/Util/Func/**` | `Tools/Basic/include/MathUtils.h`、`Tools/Basic/include/NumericAggregation.h`、`Tools/Basic/include/StringTokenUtils.h`、`Tools/Basic/include/ModelRegistry.h`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_func` target family，并具备 `surface_util_func_*` metadata、Basic 单测、smoke 与 reference matrix 证据 | `R44-Y（已完成）` |
| `Util/DistribFunc` | `ref/SPIS/num-master/src/main/java/spis/Util/DistribFunc/**` | `Tools/Particle/include/SurfaceDistributionFunction.h`、`Tools/Particle/src/SurfaceDistributionFunction.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_distrib_func` target family，并具备 `surface_util_distrib_func_*` metadata、Particle 单测、smoke 与 reference matrix 证据 | `R45-Z（已完成）` |
| `Util/Exception` | `ref/SPIS/num-master/src/main/java/spis/Util/Exception/**` | `Tools/Basic/include/ErrorCode.h`、`Tools/Basic/include/ErrorHandling.h`、`Tools/Basic/include/VoidResult.h`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_exception` target family，并具备 `surface_util_exception_*` metadata、Basic 单测、smoke 与 reference matrix 证据 | `R45-Z（已完成）` |
| `Util/Instrument` | `ref/SPIS/num-master/src/main/java/spis/Util/Instrument/**` | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Main/SurfaceScenarioLoader.cpp`、`Main/main.cpp` | `已完成` | 已形成独立 `util_instrument` target family，并具备 `surface_util_instrument_*` metadata、observer/instrument smoke 与 reference matrix 证据 | `R46-AA（已完成）` |
| `Util/List` | `ref/SPIS/num-master/src/main/java/spis/Util/List/**` | `Tools/Particle/include/ParticleSource.h`、`Tools/Particle/src/ParticleSource.cpp`、`Tools/Particle/include/SurfaceDistributionFunction.h`、`Tools/Particle/src/SurfaceDistributionFunction.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_list` target family，并具备 `surface_util_list_*` metadata、Particle 单测、smoke 与 reference matrix 证据 | `R47-AB（已完成）` |
| `Util/OcTree` | `ref/SPIS/num-master/src/main/java/spis/Util/OcTree/**` | `Tools/Mesh/include/MeshPartitioning.h`、`Tools/Mesh/src/MeshPartitioning.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_octree` target family，并具备 `surface_util_octree_*` metadata、Mesh 单测、smoke 与 reference matrix 证据 | `R47-AB（已完成）` |
| `Util/OcTreeMesh` | `ref/SPIS/num-master/src/main/java/spis/Util/OcTreeMesh/**` | `Tools/Mesh/include/MeshAlgorithms.h`、`Tools/Mesh/include/MeshPartitioning.h`、`Tools/Mesh/src/MeshAlgorithms.cpp`、`Tools/Mesh/src/MeshPartitioning.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_octree_mesh` target family，并具备 `surface_util_octree_mesh_*` metadata、Mesh 单测、smoke 与 reference matrix 证据 | `R47-AB（已完成）` |
| `Util/Sampler` | `ref/SPIS/num-master/src/main/java/spis/Util/Sampler/**` | `Tools/Particle/include/SurfaceDistributionFunction.h`、`Tools/Particle/src/SurfaceDistributionFunction.cpp`、`Toolkit/Surface Charging/include/PicMccSurfaceCurrentSampler.h`、`Toolkit/Surface Charging/src/PicMccSurfaceCurrentSampler.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_sampler` target family，并具备 `surface_util_sampler_*` metadata、单测、smoke 与 reference matrix 证据 | `R44-Y（已完成）` |
| `Util/Monitor` | `ref/SPIS/num-master/src/main/java/spis/Util/Monitor/**` | `Tools/Diagnostics/include/Monitor.h`、`Tools/Diagnostics/include/FieldMonitor.h`、`Tools/Diagnostics/include/MonitorManager.h`、`Tools/Diagnostics/src/MonitorManager.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_monitor` target family，并具备 `surface_util_monitor_*` metadata、Diagnostics 单测、smoke 与 reference matrix 证据 | `R44-Y（已完成）` |
| `Util/io` | `ref/SPIS/num-master/src/main/java/spis/Util/io/**` | `Tools/Output/include/ResultExporter.h`、`Tools/Output/include/CSVExporter.h`、`Tools/Output/src/ResultExporter.cpp`、`Main/SurfaceScenarioLoader.cpp`、`Main/main.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_io` target family，并具备 `surface_util_io_*` metadata、Output 单测、smoke 与 reference matrix 证据 | `R44-Y（已完成）` |
| `Util/Part` | `ref/SPIS/num-master/src/main/java/spis/Util/Part/**` | `Tools/Particle/include/ParticleDefinitions.h`、`Tools/Particle/include/ParticleSource.h`、`Tools/Particle/src/ParticleDefinitions.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_part` target family，并具备 `surface_util_part_*` metadata、ParticleDefinitions 单测、smoke 与 reference matrix 证据 | `R45-Z（已完成）` |
| `Util/Phys` | `ref/SPIS/num-master/src/main/java/spis/Util/Phys/**` | `Tools/Basic/include/Constants.h`、`Tools/Basic/include/MathUtils.h`、`Tools/Particle/include/ParticleDefinitions.h`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_phys` target family，并具备 `surface_util_phys_*` metadata、Basic/Particle 单测、smoke 与 reference matrix 证据 | `R45-Z（已完成）` |
| `Util/Table` | `ref/SPIS/num-master/src/main/java/spis/Util/Table/**` | `Tools/Particle/include/ParticleSource.h`、`Tools/Particle/include/SurfaceDistributionFunction.h`、`Tools/Particle/src/SurfaceDistributionFunction.cpp`、`Tools/Output/include/DataAnalyzer.h`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_table` target family，并具备 `surface_util_table_*` metadata、Particle/Output 单测、smoke 与 reference matrix 证据 | `R45-Z（已完成）` |
| `Util/Matrix` | `ref/SPIS/num-master/src/main/java/spis/Util/Matrix/**` | `Tools/Solver/include/LinearSystem.h`、`Tools/Solver/include/SparseMatrix.h`、`Tools/Solver/include/BandMatrix.h`、`Tools/Solver/src/LinearSystem.cpp`、`Tools/Solver/src/SurfaceSolverFacade.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_matrix` target family，并具备 `surface_util_matrix_*` metadata、Solver 单测、smoke 与 reference matrix 证据 | `R44-Y（已完成）` |
| `Util/Vect` | `ref/SPIS/num-master/src/main/java/spis/Util/Vect/**` | `Tools/Geometry/include/Vector3D.h`、`Tools/Geometry/include/Point3D.h`、`Tools/Geometry/src/Vector3D.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | 已形成独立 `util_vect` target family，并具备 `surface_util_vect_*` metadata、Geometry 单测、smoke 与 reference matrix 证据 | `R45-Z（已完成）` |
| `Circ/**` | `ref/SPIS/num-master/src/main/java/spis/Circ/**` | `Tools/Coupling/include/SurfaceCircuitCoupling.h`、`Tools/Coupling/src/SurfaceCircuitCoupling.cpp` | `已完成` | `Circ/Circ`、`CircField`、`Circ/DIDV` 已完成 direct route 与 matrix 证据闭环；转入回归维护 | `R31-E3`~`R31-E5（已完成）` |
| `Solver/Util` | `ref/SPIS/num-master/src/main/java/spis/Solver/Util/**` | `Tools/Solver/include/LinearSystem.h`、`Tools/Solver/src/LinearSystem.cpp`、`Tools/FieldSolver/include/PoissonSolver.h`、`Tools/FieldSolver/src/PoissonSolver.cpp`、`Tools/PICcore/src/PICCycle.cpp` | `已完成` | `SolverUtil`、`LeastSquare` 已形成独立 target family、runtime metadata 与 matrix 证据链；转入回归维护 | `R37-K（已完成）` |
| `Solver/**` | `ref/SPIS/num-master/src/main/java/spis/Solver/**` | `Tools/Solver/**`、`Tools/FieldSolver/**`、`Tools/Particle/**`、`Tools/PICcore/**` | `已完成` | `Circuit`、`ElectroMag`、`Matter` 已完成 direct-SPIS semantic evidence 闭环；转入回归维护 | `R32-F2`~`R32-F4（已完成）` |
| `Vol/**` | `ref/SPIS/num-master/src/main/java/spis/Vol/**` | `Tools/FieldSolver/include/SurfaceFieldVolumeBridge.h`、`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp` 及相关 native volume 路径 | `已完成` | `BC`、`Field`、`Distrib`、`Interact` 已完成显式实现、smoke 与 matrix 证据闭环；转入回归维护 | `R28-B`~`R30-D`、`R33-G（已完成）` |
| `Top/Top` | `ref/SPIS/num-master/src/main/java/spis/Top/Top/**` | `Toolkit/Surface Charging/include/SurfaceScenarioCatalog.h`、`Toolkit/Surface Charging/src/SurfaceScenarioCatalog.cpp`、`Toolkit/Surface Charging/include/SurfaceSimulationRunner.h`、`Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`、`Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`、`Main/*` | `已完成` | `Scenario`、`PotentialSweep`、`UIInvokable`、`SpisTopMenu`、`NumTopFromUI` 已形成独立 target family 与 route / matrix 证据链；转入回归维护 | `R38-L（已完成）` |
| `Top/**` | `ref/SPIS/num-master/src/main/java/spis/Top/**` | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`、`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Main/*` | `已完成` | transition/lifecycle 对象层级、route evidence 与 lifecycle metadata 已闭环；转入回归维护 | `R31-E2（已完成）` |
| `Vol/BC(Abstract)` | `ref/SPIS/num-master/src/main/java/spis/Vol/BC/BC.java`、`DirichletPoissonBC.java`、`MatterBC.java`、`TestableMatterBC.java`、`VoltageGenerator.java` | `Tools/FieldSolver/include/BoundaryCondition.h`、`Tools/FieldSolver/include/PoissonSolver.h`、`Tools/Solver/include/LinearSystem.h`、`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | `BC`、`DirichletPoissonBC`、`MatterBC`、`TestableMatterBC`、`VoltageGenerator` 已形成独立 target family 与 metadata / gate；转入回归维护 | `R40-N（已完成）` |
| `Vol/VolField(Abstract)` | `ref/SPIS/num-master/src/main/java/spis/Vol/VolField/VolField.java`、`ScalVolField.java`、`VectVolField.java`、`EField.java`、`DirVectVolField.java`、`DirEField.java`、`PotEField.java`、`PotVectVolField.java`、`BField.java` | `Tools/FieldSolver/include/VolField.h`、`Tools/FieldSolver/src/VolField.cpp`、`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | `VolField` 抽象层已形成独立 target family 与 `surface_native_volume_field_*` metadata / gate；转入回归维护 | `R41-O（已完成）` |
| `Vol/VolDistrib(LongTail)` | `ref/SPIS/num-master/src/main/java/spis/Vol/VolDistrib/**` | `Tools/FieldSolver/include/VolDistrib.h`、`Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`、`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`、`Main/SurfaceScenarioLoader.cpp`、`Main/main.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | `VolDistrib`、`NonPICVolDistrib`、`AnalyticVolDistrib`、`MultipleVolDistrib`、`LocalMaxwellVolDistrib`、`LocalMaxwellBoltzmannVolDistrib`、`GlobalMaxwellBoltzmannVolDistrib`、`PICBoltzmannVolDistrib`、`SteadyMaxwellBoltzmannVolDistrib`、`UnlimitedGlobalMaxwellBoltzmannVolDistrib`、`SurfaceLimitedGlobalMaxwellBoltzmannVolDistrib`、`TrunckatedGlobalMaxwellBoltzmannVolDistrib`、`ImplicitableVolDistrib`、`Updatable` 已形成独立 target family 与 `surface_native_volume_distribution_*` metadata / gate；转入回归维护 | `R42-P（已完成）` |
| `Vol/VolInteract(LongTail)` | `ref/SPIS/num-master/src/main/java/spis/Vol/VolInteract/VolInteractor.java`、`VolInteractModel.java`、`VolInteractionAlongTrajectory.java`、`ElasticCollisions.java`、`ConstantIonizationInteractor.java` | `Tools/FieldSolver/include/VolInteract.h`、`Tools/FieldSolver/src/VolInteract.cpp`、`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | `VolInteractor`、`VolInteractModel`、`VolInteractionAlongTrajectory`、`ElasticCollisions`、`ConstantIonizationInteractor` 已形成独立 target family 与 `surface_native_volume_interaction_long_tail_*` metadata / gate；转入回归维护 | `R43-Q（已完成）` |

### 2.10 SPIS class-family 盘点基线（首批明确点名）

| 能力维度 | SPIS family | 当前现状 | 对齐判断 | 后续动作 |
| --- | --- | --- | --- | --- |
| Surf Material | `BasicMaterialModel`、`ErosionMaterialModel` | 已补齐 `BasicSurfaceMaterialModel`、`ErosionSurfaceMaterialModel` 的参数资产、默认映射与 direct evidence | `已完成` | 转入长期回归维护 |
| Surf Interact | `MaterialInteractor`、`MaterialDFInteractor`、`GenericDFInteractor`、`MaxwellianInteractor`、`MaxwellianInteractorWithRecollection`、`MultipleInteractor`、`MultipleMaxwellianInteractor`、`YieldInteractor`、`ReflectionInteractor`、`ErosionInteractor`、`ImprovedPhotoEmInteractor`、`TabulatedSEYModel`、`RecollectedSEYModel`、`DefaultPEEModel`、`DefaultSEEEYModel`、`DefaultSEEPModel`、`DefaultErosionModel`、`BasicInducedConductInteractor`、`Device` | 已完成显式 family 入口、active signature、smoke 与 matrix evidence 闭环 | `已完成` | 转入长期回归维护 |
| Surf Distrib | `AxisymTabulatedVelocitySurfDistrib`、`LocalModifiedPearsonIVSurfDistrib`、`NonLocalizedSurfDistrib`、`MultipleSurfDistrib`、`LocalTabulatedSurfDistrib`、`TwoAxesTabulatedVelocitySurfDistrib`、`FowlerNordheimSurfDistrib`、`GlobalMaxwellBoltzmannSurfDistrib`、`GlobalMaxwellBoltzmannSurfDistrib2`、`GlobalMaxwellSurfDistrib`、`LocalMaxwellSurfDistrib`、`LocalMaxwellSurfDistrib2`、`RecollMaxwellSurfDistrib`、`PICSurfDistrib`、`NonPICSurfDistrib`、`GenericSurfDistrib`、`GlobalSurfDistrib`、`LocalGenericSurfDistrib`、`TestableSurfDistrib`、`TestableForASurfDistrib`、`MaxwellianThruster` | 当前已覆盖上述全量 family，并完成 `distribution_model` 解析、`surface_distribution_family` sidecar 导出与 matrix gate 覆盖 | `已完成` | 转入 `Round 23`，仅保留回归 gate 维护 |
| Surf Field | `OMLCurrentScaler`、`LTE_OML_CurrentScaler`、`FowlerNordheimCurrentScaler`、`VariableBarrierCurrentScaler`、`AutomaticBarrierCurrentScaler`、`MultipleCurrentScaler`、`CurrentScalerFromCurrentVariation`、`CurrentScalerFromLocalCurrentVariation` | 当前已补齐全部 scaler family 的 direct route / config / smoke / matrix evidence | `已完成` | 转入长期回归维护 |
| Top Transition | `LocalTimeTransition`、`SpinningSpacecraft`、`SunFluxUpdater`、`ConductivityEvolution`、`SourceFluxUpdater`、`SimulationParamUpdater`、`SheathOrPresheathPoissonBCUpdater`、`RCCabsSCUpdater`、`VcrossBfieldUpdater`、`BasicEclipseExit`、`TransientArtificialSources`、`LangmuirProbeTransition`、`TransitionObserver`、`Finalization`、`SunFluxIntensityUpdater` | 当前主线事件、extended events、lifecycle hooks 与 route evidence 已闭环 | `已完成` | 转入长期回归维护 |
| Vol BC / Field / Distrib / Interact | `VoltageDependentMBC`、`MixedDirichletFourierPoissonBC`、`FourierPoissonBC`、`SurfDistribMatterBC`、`OneSurfDistribTestableMatterBC`、`CapacitiveVoltageGenerator`、`UniformBField`、`SolenoidBField`、`DipolarBField`、`MultipleEField`、`PICVolDistrib`、`PICVolDistribNoAcc`、`PICVolDistribUpdatable`、`SmartPICVolDistrib`、`CompositeVolDistrib`、`BackTrackingVolDistrib`、`BacktrackingPICCompositeVolDistrib`、`BacktrackingBoltzmannCompositeVolDistrib`、`MCCInteractor`、`CEXInteractor`、`PhotoIonization`、`TrajectoryInteractionFromField`、`SpinningSpacecraftTrajectory` | 当前已补齐四组 family 的显式实现、active signature、smoke 与 matrix evidence | `已完成` | 转入长期回归维护 |

#### 2.10.1 `Round 24` Top/Transition 对象层映射（2026-04-16）

| SPIS 对象层 | SPIS family / 角色 | 当前 SCDAT 落点 | 当前证据 |
| --- | --- | --- | --- |
| interface layer | `TransitionInterface`、`Transition`、`SimulationParamUpdater` | `SurfaceTransitionObjectLayerState.interface_layer_family_signature` / `active_interface_layer_family_signature` | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`SurfaceTransitionEngineTest.ExposesObserverFinalizationAndSunFluxIntensityObjectHooks` |
| transition family | `LocalTimeTransition`、`SpinningSpacecraft`、`ConductivityEvolution`、`BasicEclipseExit`、`TransientArtificialSources`、`LangmuirProbeTransition` | `SurfaceAdvanceTransition` 事件求值 + `object_layer.active_family_signature` | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp` |
| updater family | `SunFluxUpdater`、`SourceFluxUpdater`、`SimulationParamUpdater`、`SheathOrPresheathPoissonBCUpdater`、`RCCabsSCUpdater`、`VcrossBfieldUpdater` | `SurfaceAdvanceTransition.extended_events` + `object_layer.active_family_signature` | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`SurfaceTransitionEngineTest.EvaluatesExtendedEventBundleWhenSunFluxIsDisabled` |
| lifecycle family | `TransitionObserver`、`Finalization`、`SunFluxIntensityUpdater` | `SurfaceAdvanceTransition.observer` / `finalization` / `sun_flux_intensity` + `object_layer.lifecycle_family_signature` | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp` |
| sidecar / route evidence | `Top` family 可观测字段与 route 证据链 | metadata `surface_top_transition_*` 系列字段 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsTopTransitionFamilyMetadata` |

#### 2.10.2 `Round 25` Circ/Solver family 映射（2026-04-16）

| 域 | SPIS family / 层次 | 当前 SCDAT 落点 | 当前证据 |
| --- | --- | --- | --- |
| `Circ/Circ` | `Circuit`、`ElecComponent`、`SurfaceComponent`、`RCCabsCirc`、`RLCCirc`、`SIN`、`PULSE`、`PWL`、`EXP` | `Tools/Coupling/include/SurfaceCircuitCoupling.h` 中的 `CircuitAssembly` / `SurfaceCircuitCoupling` / `CircuitExcitationDescriptor`，以及 `resolveCircuitFamilyView(...)` 输出的 family 签名 | `Tools/Coupling/src/SurfaceCircuitCoupling.cpp`；`Tools/Coupling/test/Coupling_test.cpp`；`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsCircFamilyMetadata` |
| `Circ/CircField` | `CircField`、`CurrentDistribution`、`DirCircField` | `SurfaceCircuitLinearization` / `SurfaceCircuitKernelInput` 与 `resolveCircuitFamilyView(...)` 的 `circfield` family 视图 | `Tools/Coupling/include/SurfaceCircuitCoupling.h`；`Tools/Coupling/test/Coupling_test.cpp`；`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsCircFamilyMetadata` |
| `Circ/DIDV` | `DIDV`、`SurfDIDV*`、`WeightedSurfDIDV`、`MultipleSurfDIDV`、`ReducableDIDV`、`DIDVOnCirc` | `composeCircuitDidv(...)`、`composeSurfDidvFromMatrices(...)`、`reduceCircuitDidvFrom*` 与 `resolveCircuitFamilyView(...)` 的 `didv` family 视图 | `Tools/Coupling/src/SurfaceCircuitCoupling.cpp`；`CouplingTest.SurfDidvFromMatricesBuildsExplicitNodeAndOffDiagonalRoute`；`CouplingTest.ReducedSurfDidvAggregatesMappedSurfaceMatrices`；`SurfaceChargingSmokeTest.ReferenceDidvMatchesFiniteDifferenceAcrossGeoLeoRoutes` |
| `Solver/Circuit` | `CircSolve` | `Tools/Solver/include/SurfaceSolverFacade.h` 中的 `SurfaceSolverFamilyView` + surface runtime circuit solve 路径 | `Tools/Solver/src/SurfaceSolverFacade.cpp`；`SurfaceSolverFacadeTest.ResolveSurfaceSolverFamilyViewBuildsSupportedAndActiveSignatures`；`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsSolverFamilyMetadata` |
| `Solver/ElectroMag` | `AbstractEMSolver`、`PoissonSolver`、`PotPoissonSolver`、`ConjGrad3DUnstructPoissonSolver`、`PoissonMagnetostaticSolver` | `resolveVolumeLinearSolverRouting(...)`、`solveDenseLinearSystemWithResidual(...)`、`solveIterativeLinearSystem(...)` 与 `SurfaceSolverFamilyView.electromag_*` | `Tools/Solver/src/SurfaceSolverFacade.cpp`；`Tools/Solver/test/SurfaceSolverFacade_test.cpp`；`SurfaceChargingSmokeTest.SolverConfigImplicitCouplingDrivesSharedGlobalCoupledPolicy`；`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsSolverFamilyMetadata` |
| `Solver/Matter` | `ParticlePusher`、`PICPusher`、`CrossTetra*` 家族 | live-PIC / shared-surface route 与 `SurfaceSolverFamilyView.matter_*` | `Tools/Solver/include/SurfaceSolverFacade.h`；`Tools/Solver/src/SurfaceSolverFacade.cpp`；`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsSolverFamilyMetadata`；`SurfaceChargingSmokeTest.PicCircuitPresetsKeepLivePicMccEnabled` |

#### 2.10.3 `Round 26` family coverage artifact 基线（2026-04-16）

| artifact | 作用 | 当前证据 |
| --- | --- | --- |
| `documentation/contracts/surface_family_coverage_matrix_v1.json` | 冻结 `family coverage matrix` artifact 的 section / field / tracked family 契约 | `documentation/contracts/surface_family_coverage_matrix_v1.json` |
| `scripts/run/surface_spis_family_coverage_matrix.json` | 统一承载 `family-to-code`、`family-to-test`、`family-to-config-route` 三张 machine-readable matrix | `scripts/run/surface_spis_family_coverage_matrix.json` |
| `scripts/python/check_surface_family_coverage_gate.py` | 验证每个 tracked family 是否同时具备代码、测试、config/route 三类证据，并生成 gate report | `scripts/python/check_surface_family_coverage_gate.py`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.md` |

#### 2.10.4 `Post-R26` 状态升级准则（2026-04-16）

| 当前口径 / 状态 | 是否视为 1:1 已达成 | 必须补齐的附加条件 | 后续归属 |
| --- | --- | --- | --- |
| `near_equivalent` / `近似实现` | `否` | 必须补齐显式代码入口、显式 config/route 入口、单测、smoke 或 matrix case，并在 `surface_spis_family_coverage_matrix.json` 中回填为 `completed` | 进入 `Round 27`~`Round 33` remediation 主线 |
| `已完成（family 可观测）` | `否` | 必须从“仅 parity/signature 可观测”升级到“显式实现 + direct evidence”；仅 metadata / sidecar / parity signature 不足以升为 `completed` | 进入 `Round 28`~`Round 30` remediation 主线 |
| `近似可达` | `否` | 必须补齐独立 family route、主线路由证据与至少一个 direct smoke / matrix case | 进入 `Round 28`~`Round 30` remediation 主线 |
| `未显式实现` | `否` | 必须先补显式 family 入口，再补 supported/active family signature、单测和 matrix 证据 | 进入 `Round 28`~`Round 30` remediation 主线 |

- `Round 22`~`Round 26` 中出现的 `near_equivalent`、`近似可达`、`未显式实现`、`已完成（family 可观测）` 一律不得视为历史已收口项。
- 这些条目在 `Post-R26` 阶段自动纳入 `Round 27`~`Round 33` 的 remediation 主线，直到其上层 family 在 `scripts/run/surface_spis_family_coverage_matrix.json` 中升级为 `completed`。

### 2.11 1:1 对齐达成标准

- 本文档中的“达成 1:1 对齐”仅在 `Surface Charging` 相关范围内成立，不含 `ref/SPIS/ui-master/**`、`ref/SPIS/instruments-master/**`、`ref/SPIS/num-plugins-master/**` 与 `ref/SPIS/num-master/**/osgi/**` 等非当前目标范围。
- 只有同时满足下表全部条件，才允许宣告“与 SPIS 功能 1:1 对齐完成”：

| 判定项 | 必须满足条件 |
| --- | --- |
| family 覆盖 | `Surf/Circ/Solver/Vol/Util/Top` 目标 family 全部归类为 `已等价实现`、`近似实现`、`尚无实现`、`不纳入当前目标范围` 之一，且不存在未归类项 |
| family 实现状态 | 除明确排除项外，所有目标 family 均达到 `已等价实现` |
| 配置/路由 | 每个目标 family 都有明确配置、路由或运行入口证据 |
| 测试/回归 | 每个目标 family 都至少挂接到单测、smoke、gate 或参考矩阵中的一种 |
| 文档一致性 | 本文档不再把 `近似实现` 或 `最小闭环` 误写成 `已对齐` |
| 外部兼容 | 现有 `SurfaceChargingConfig`、preset、route、sidecar 契约不发生未记录破坏 |

### 2.12 当前结论（2026-04-17）

- `Round 33` 封板之后，`Round 34`~`Round 47` 已继续把 `Util/**`、`Util/DistribFunc+Exception+Instrument+Func+List+OcTree+OcTreeMesh+Part+Phys+Sampler+Monitor+io+Table+Matrix+Vect`、`Top/Top+Simulation+Plasma+SC+Grid+Default`、`Vol/BC(Abstract)+VolField(Abstract)+VolDistrib(LongTail)+VolInteract(LongTail)+Geom+VolMesh`、`Solver/Util`、`Surf/SurfMesh` 纳入扩展 target families，并重新跑通 `family coverage gate`、全量 `surface_reference_matrix gate` 与 sealoff。
- `build_codex/surface_round47_sealoff/surface_round26_sealoff.json` 当前给出的最终 verdict 为 `achieved`，因此按本文档 2.11 节的 1:1 判定标准，当前项目**已达到**与 `ref/SPIS/num-master` 的功能 1:1 对齐。
- 为避免把 family 封板误读成 class-level 逐对象镜像已全部完成，当前已新增：
  - `Toolkit/Surface Charging/docs/surface_spis_class_audit.json`
  - `Toolkit/Surface Charging/docs/surface_spis_class_audit.csv`
  - `Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md`
- 上述逐 class 审计产物当前覆盖 `ref/SPIS/num-master/src/main/java/spis/**` 下 `602` 个 class；其用途是补充回答“family 已封板之后，哪些 SPIS class 在当前项目中由单类承载、哪些由聚合对象层承载、哪些仍只达到近似实现”。该审计**不推翻**当前 family sealoff，而是作为下一阶段对象层收口与结构优化基线。
- 当前目标 family 中的状态分布已经更新为：
  - `completed`：`surf_material`、`surf_interact`、`surf_distrib`、`surf_mesh`、`surf_field`、`util_full_stack`、`util_distrib_func`、`util_exception`、`util_instrument`、`util_func`、`util_list`、`util_octree`、`util_octree_mesh`、`util_part`、`util_phys`、`util_sampler`、`util_monitor`、`util_io`、`util_table`、`util_matrix`、`util_vect`、`top_transition`、`top_top`、`top_simulation`、`top_plasma`、`top_sc`、`top_grid`、`top_default`、`vol_full_stack`、`vol_bc_abstract`、`vol_field_abstract`、`vol_distrib_long_tail`、`vol_interact_long_tail`、`vol_geom`、`vol_mesh`、`circ_circuit`、`circ_circfield`、`circ_didv`、`solver_circuit`、`solver_util`、`solver_electromag`、`solver_matter`
  - `unresolved`：无
- 当前结论应表述为：“**Round 31 / Round 32 已完成，Round 33 已完成最终重封板；Round 34 / Round 35 / Round 36 / Round 37 / Round 38 / Round 39 / Round 40 / Round 41 / Round 42 / Round 43 / Round 44 / Round 45 / Round 46 / Round 47 已完成 scope expansion 收口；Surface/Circ/Solver/Vol/Util/Top 目标 family 已全部升级为 `completed`，项目在本文档约束范围内已达成与 SPIS 的 1:1 对齐。**”
- 当前 class-level 审计结论应表述为：“**family 层 1:1 对齐已达成，但 class-level 对象层仍存在聚合承载差异；下一阶段工作的目标不再是新增 in-scope family，而是补 Top/Default、Top/Top+Simulation、Util/Func、Util/Instrument 的 class-to-runtime / class-to-helper 显式映射与可审计证据。**”

### 2.13 当前阶段性目标

- `Wave A` 目标：完成文档重构与全量差异盘点封板。
- `Wave B` 目标：完成 `SurfInteract` family 1:1 对齐。
- `Wave C` 目标：完成 `SurfDistrib` family 1:1 对齐。
- `Wave D` 目标：完成 `Vol` full stack 1:1 对齐。
- `Wave E` 目标：完成 `Top` lifecycle / transition 结构对齐。
- `Wave F` 目标：完成 `Circ/Solver` 剩余语义复核与补差。
- `Wave G` 目标：建立 1:1 对齐 gate、family coverage matrix 与最终封板证据链。
  当前状态：`已完成最终封板`；`family coverage gate`、`surface_reference_matrix gate` 与 sealoff 全部 `PASS`，最终 verdict 为 `achieved`。
- `Wave O` 目标：把 `spis/Util/**` 从“Surface 刚需支撑域”升级为扩展 target family，并补 `util_full_stack` contract / gate。
  当前状态：`已完成`；`util_full_stack` 已进入 `surface_spis_family_coverage_matrix.json`，并通过 smoke / reference matrix / coverage gate。
- `Wave P` 目标：把 `Top/**` 从以 `Transition` 为主的 gate 扩展到 `Simulation`、`Plasma`、`SC`、`Grid`、`Default` 逐族盘点。
  当前状态：`已完成`；`top_simulation`、`top_plasma`、`top_sc`、`top_grid`、`top_default` 已具备独立 metadata / contract / gate。
- `Wave Q` 目标：把 `Vol/Geom`、`Vol/VolMesh` 从“能力闭环”升级为独立 family 级 contract / gate。
  当前状态：`已完成`；`vol_geom`、`vol_mesh` 已进入 target families，并通过 native-volume metadata / reference matrix / coverage gate。
- `Wave R` 目标：把 `Solver/Util` 从 `Solver/**` 的隐含支撑层升级为独立 target family，并补 `SolverUtil` / `LeastSquare` 的 runtime metadata / tests / gate。
  当前状态：`已完成`；`solver_util` 已进入 target families，并通过 smoke / reference matrix / coverage gate / sealoff。
- `Wave S` 目标：把 `Top/Top` 从“Top 的剩余入口层”升级为独立 target family，并补 `Scenario` / `PotentialSweep` / `UIInvokable` / `SpisTopMenu` / `NumTopFromUI` 的 runtime metadata / tests / gate。
  当前状态：`已完成`；`top_top` 已进入 target families，并通过 smoke / reference matrix / coverage gate / sealoff。
- `Wave T` 目标：把 `Surf/SurfMesh` 从 structured topology / boundary mapping 的隐含能力升级为独立 target family，并补 `SurfMesh` / `ThreeDUnstructSurfMesh` 的 runtime metadata / tests / gate。
  当前状态：`已完成`；`surf_mesh` 已进入 target families，并通过 smoke / reference matrix / coverage gate / sealoff。
- `Wave U` 目标：把 `Vol/BC` 抽象长尾从 `vol_full_stack` 的隐含层升级为独立 target family，并补 `BC` / `DirichletPoissonBC` / `MatterBC` / `TestableMatterBC` / `VoltageGenerator` 的 runtime metadata / tests / gate。
  当前状态：`已完成`；`vol_bc_abstract` 已进入 target families，并通过 reference matrix / coverage gate / sealoff。
- `Wave V` 目标：把 `Vol/VolField` 抽象长尾从 `vol_full_stack` 的隐含层升级为独立 target family，并补 `VolField` / `ScalVolField` / `VectVolField` / `EField` / `DirVectVolField` / `DirEField` / `PotEField` / `PotVectVolField` / `BField` 的 runtime metadata / tests / gate。
  当前状态：`已完成`；`vol_field_abstract` 已进入 target families，并通过 reference matrix / coverage gate / sealoff。
- `Wave X` 目标：把 `Vol/VolInteract` 长尾抽象/兼容层升级为独立 target family，并补 `VolInteractor` / `VolInteractModel` / `VolInteractionAlongTrajectory` / `ElasticCollisions` / `ConstantIonizationInteractor` 的 runtime metadata / tests / gate。
  当前状态：`已完成`；`vol_interact_long_tail` 已进入 target families，并通过 reference matrix / coverage gate / sealoff。
- `Wave Y` 目标：把 `Util/**` 中最直接影响 surface 主线的子包升级为独立 target family，并补 `util_func` / `util_sampler` / `util_monitor` / `util_io` / `util_matrix` 的 runtime metadata / tests / gate。
  当前状态：`已完成`；五个 util 子包级 target family 已进入 `surface_spis_family_coverage_matrix.json`，并通过 smoke / reference matrix / coverage gate / sealoff。
- `Wave Z` 目标：继续把 `Util/**` 中具备稳定代码与测试落点的子包升级为独立 target family，并补 `util_distrib_func` / `util_exception` / `util_part` / `util_phys` / `util_table` / `util_vect` 的 runtime metadata / tests / gate。
  当前状态：`已完成`；六个 util 子包级 target family 已进入 `surface_spis_family_coverage_matrix.json`，并通过 smoke / reference matrix / coverage gate / sealoff。
- `Wave AA` 目标：把 `Util/Instrument` 从 `util_full_stack` 聚合层升级为独立 target family，并补 `util_instrument` 的 runtime metadata / tests / gate。
  当前状态：`已完成`；`util_instrument` 已进入 `surface_spis_family_coverage_matrix.json`，并通过 observer/instrument smoke / reference matrix / coverage gate / sealoff。
- `Wave AB` 目标：把 `Util/List`、`Util/OcTree`、`Util/OcTreeMesh` 从 `util_full_stack` 聚合层升级为独立 target family，并补 `util_list`、`util_octree`、`util_octree_mesh` 的 runtime metadata / tests / gate。
  当前状态：`已完成`；三个 util 子包级 target family 已进入 `surface_spis_family_coverage_matrix.json`，并通过 Particle/Mesh 单测、smoke、reference matrix、coverage gate 与 sealoff。
- `Wave AC` 目标：在 family sealoff 不变的前提下，补 `Top/Default` 的 class-to-runtime 对象层职责索引，并将 `top_default` 中 class-level `近似实现` 项从“仅审计可见”升级为“文档可追踪”。
  当前状态：`已完成`；`Top/Default` class-to-runtime 映射文档已落地，class-level 差异已进入文档化、可回链状态。
- `Wave AD` 目标：补 `Top/Top` 与 `Top/Simulation` 的 shell-to-runtime 对照，说明 `UIInvokable`、`SpisTopMenu`、`SimulationListener` 等壳层 class 在当前项目中如何由 runtime plan / runner / scenario 入口吸收。
  当前状态：`已完成`；`Top/Top + Top/Simulation` shell-to-runtime 对照文档已落地，class-level 壳层吸收关系已文档化。
- `Wave AE` 目标：为 `Util/Func` 建立 class-to-helper / class-to-formula 映射索引，区分“已由基础函数库聚合承载”“仍缺显式命名桥接”“不建议做逐类镜像”的三类结论。
  当前状态：`已完成`；`Util/Func` 映射文档已落地，`98` 个 class 已逐项归入三类之一，并抽出了桥接候选函数族。
- `Wave AF` 目标：为 `Util/Instrument` 建立 instrument catalog / observer 路由 / smoke 断言映射，区分“已有 runtime 行为”“仅 metadata route”“尚未单独可观测”的三类仪器对象。
  当前状态：`已完成`；`Util/Instrument` 映射文档已落地，`51` 个 class 已按三类完成逐项归档，class-level 主线收口。
- `Post-R26 remediation` 目标：将 `Round 22`~`Round 26` 中所有 `near_equivalent`、`近似可达`、`未显式实现`、`已完成（family 可观测）` 条目逐项拉平为 `completed`，并在 `Round 33` 重新生成 `final_alignment_verdict = achieved`。
  当前进展：`Round 27`~`Round 33` 已收口，`Round 34`~`Round 47` 已完成扩展 scope 收口；全部 target family 已升级为 `completed`；项目转入长期回归维护。

### 2.14 当前仍未实现 / 未纳入的后续缺口（Round 47 执行后）

当前 in-scope backlog 已清零。下表仅保留最后一批已收口项与明确排除项，用于回答“若后续继续扩边界，还剩哪些不在当前 verdict 内的内容”。这些条目不是当前 sealoff target families 的一部分，因此**不影响** 2.12 节现有 `achieved` 结论。

| 主题 | SPIS 源路径 | 当前状态 | 当前差异 | 后续处理 |
| --- | --- | --- | --- | --- |
| `Vol/VolField` 抽象层长尾 | `ref/SPIS/num-master/src/main/java/spis/Vol/VolField/**` 中 `VolField`、`ScalVolField`、`VectVolField`、`EField`、`DirVectVolField`、`DirEField`、`PotEField`、`PotVectVolField`、`BField` 等 | `已完成` | 已提升为独立 `vol_field_abstract` target family，并具备 `surface_native_volume_field_*` metadata、smoke、reference matrix 与 gate 证据 | 转入长期回归维护 |
| `Vol/VolDistrib` 长尾 | `ref/SPIS/num-master/src/main/java/spis/Vol/VolDistrib/**` 中 `AnalyticVolDistrib`、`NonPICVolDistrib`、`MultipleVolDistrib`、`LocalMaxwellVolDistrib`、`LocalMaxwellBoltzmannVolDistrib`、`GlobalMaxwellBoltzmannVolDistrib`、`PICBoltzmannVolDistrib`、`SteadyMaxwellBoltzmannVolDistrib`、`UnlimitedGlobalMaxwellBoltzmannVolDistrib` 等 | `已完成` | 已在 `Round 42` 升级为独立 target family，并补齐 `surface_native_volume_distribution_*` metadata、config-route、smoke 与 reference-matrix 证据 | 转入长期回归维护 |
| `Vol/VolInteract` 长尾 | `ref/SPIS/num-master/src/main/java/spis/Vol/VolInteract/**` 中 `VolInteractor`、`VolInteractModel`、`VolInteractionAlongTrajectory`、`ElasticCollisions`、`ConstantIonizationInteractor` 等 | `已完成` | 已在 `Round 43` 升级为独立 target family，并补齐 `surface_native_volume_interaction_long_tail_*` metadata、smoke、reference matrix 与 gate 证据 | 转入长期回归维护 |
| `Util/**` 子包级拆分 | `ref/SPIS/num-master/src/main/java/spis/Util/**` | `已完成` | `Round 44`~`Round 47` 已独立纳入 `util_distrib_func`、`util_exception`、`util_instrument`、`util_func`、`util_list`、`util_octree`、`util_octree_mesh`、`util_part`、`util_phys`、`util_sampler`、`util_monitor`、`util_io`、`util_table`、`util_matrix`、`util_vect`；全部剩余 util 子包已脱离 `util_full_stack` 聚合层 | 转入长期回归维护 |
| 非当前边界项 | `ref/SPIS/ui-master/**`、`ref/SPIS/instruments-master/**`、`ref/SPIS/num-plugins-master/**`、`ref/SPIS/num-master/**/osgi/**` | `不纳入当前目标范围` | 当前文档边界仍锁定 `num-master` 数值内核六个域及其直接依赖，不含 UI / Instruments / Plugins / OSGi 工程；这些内容只做边界说明，不进入当前主线 verdict | 保持排除，除非后续变更边界 |

#### 2.14.1 下一阶段 backlog 总结

- 已明确但尚未纳入实现主线的 family / family-group：无。
- `Solver/Util` 已在 `Round 37` 中完成收口，说明“先文档化缺口，再新增 target family、metadata、tests、gate”的扩展模式可行。
- `Top/Top` 已在 `Round 38` 中完成收口。
- `Surf/SurfMesh` 已在 `Round 39` 中完成收口。
- `Vol/BC` 抽象长尾已在 `Round 40` 中完成收口。
- `Vol/VolField` 抽象长尾已在 `Round 41` 中完成收口。
- `Vol/VolDistrib` 长尾已在 `Round 42` 中完成收口。
- `Vol/VolInteract` 长尾已在 `Round 43` 中完成收口。
- `Util/Func`、`Util/Sampler`、`Util/Monitor`、`Util/io`、`Util/Matrix` 已在 `Round 44` 中完成收口。
- `Util/DistribFunc`、`Util/Exception`、`Util/Part`、`Util/Phys`、`Util/Table`、`Util/Vect` 已在 `Round 45` 中完成收口。
- `Util/Instrument` 已在 `Round 46` 中完成收口。
- `Util/List`、`Util/OcTree`、`Util/OcTreeMesh` 已在 `Round 47` 中完成收口。
- 当前 in-scope backlog 已清零；若后续继续扩边界，仅剩 UI / Instruments / Plugins / OSGi 工程等非当前目标范围项。

#### 2.14.2 Class-level 审计后续主线（Round 47 执行后）

`surface_spis_class_audit_summary.md` 当前给出的逐 class 统计为：

| 域 | class 数 | 已等价实现 | 聚合等价实现 | 近似实现 | 尚无实现 | 不纳入范围 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `Surf` | `100` | `55` | `45` | `0` | `0` | `0` |
| `Vol` | `69` | `54` | `15` | `0` | `0` | `0` |
| `Top` | `65` | `17` | `20` | `28` | `0` | `0` |
| `Circ` | `29` | `21` | `8` | `0` | `0` | `0` |
| `Solver` | `23` | `7` | `16` | `0` | `0` | `0` |
| `Util` | `315` | `0` | `166` | `149` | `0` | `0` |
| `osgi` | `1` | `0` | `0` | `0` | `0` | `1` |

该统计的含义是：

- `Surf`、`Vol`、`Circ`、`Solver` 当前不再存在 class-level `尚无实现` 项；后续重点是防回退，而不是继续扩 family。
- `Top` 当前最大的 class-level 弱区集中在 `Top/Default`，其次是 `Top/Top` 与 `Top/Simulation` 的壳层对象。
- `Util` 当前最大的 class-level 弱区集中在 `Util/Func` 与 `Util/Instrument`；这两块更适合补“显式映射索引 / catalog / helper-bridge”，而不是强行做 C++ 单类镜像。
- 这些 class-level `近似实现` / `聚合等价实现` 结论属于“family 已封板后的对象层差异盘点”，**不影响** 2.12 节现有 `achieved` verdict。

下一阶段 class-level backlog 统一收口为下表：

| 优先级 | 主题 | 当前 class-level 差异 | 目标产出 | 处理策略 |
| --- | --- | --- | --- | --- |
| `P0` | `Top/Default` | `25` 个 class 仍为 `近似实现`，对象层职责被 runtime structs / catalog / glue 折叠承载 | `Top/Default` class-to-runtime 映射文档 + 任务索引 + 引用证据 | 先补对象层职责索引，不急于新增实现 |
| `P1` | `Top/Top + Top/Simulation` | `UIInvokable`、`SpisTopMenu`、`SimulationListener` 等壳层 class 在当前项目中被入口层吸收 | shell-to-runtime 映射文档 + runner/scenario 引用证据 | 用文档与测试入口解释吸收路径 |
| `P2` | `Util/Func` | `98` 个 class 为 `近似实现`，大量 Java 函数对象被合并到 `MathUtils` / `Constants` / helper 中 | class-to-helper / formula 映射索引 | 按“已聚合承载 / 需显式命名桥接 / 不建议逐类镜像”三层归类 |
| `P3` | `Util/Instrument` | `51` 个 class 为 `近似实现`，当前更多体现为 instrument metadata / observer route | instrument catalog / observer-route 映射文档 + 测试矩阵索引 | 优先补可观测性与 catalog，不先做大规模实现 |

`Round 51` 完成后，上表中的 `P0`、`P1`、`P2`、`P3` 已全部从“待落地”升级为“已文档化收口”；当前 class-level backlog 已完成本轮计划范围内的收口，转入长期回归维护。

#### 2.14.3 后续代码增强 backlog（非必需，但高收益）

在当前 `sealoff = achieved` 与 class-level 文档化收口保持不变的前提下，后续代码层面最值得继续推进的增强项如下：

| 优先级 | 主题 | 目标 | 首选文件 | 收益 | 风险 |
| --- | --- | --- | --- | --- | --- |
| `E1` | `Util/Instrument` 可观测性增强 | 把当前仅停留在 observer route / metadata 的一部分 instrument class 提升为“有直接 runtime 物证” | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 最高；最直接提升 class-level 说服力 | 低 |
| `E2` | `Util/Func` 显式命名 helper | 为 `DebyeSurfSheath`、`SimpleSurfSheath`、`OMLFactor`、`Maxwellian1D`、`Kappa1_*`、`PowerLaw` 等补显式 helper 名称 | `Tools/Material/src/SurfaceMaterialModel.cpp`、`Tools/Particle/src/ParticleSource.cpp`、`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp` | 高；提升物理公式可审计性 | 中低 |
| `E3` | runtime 热点拆分 | 把 instrument/export、transition observer、shared runtime observer 从巨型文件中继续剥离 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp` | 高；长期维护收益大 | 中 |

当前执行建议顺序：

1. 先完成 `E1`，优先补 metadata / sidecar / smoke，不急于新增对象类。
2. 再做 `E2`，只补少量高价值 helper，避免追求 `Util/Func` 全量镜像。
3. 最后做 `E3`，作为中期结构优化主线。

当前进展补记：

- `E1` 已落第二批增强：instrument observability metadata 已抽离为独立 helper，并新增 `monitor / particle_detector / metadata_catalog` artifact 与 observability 口径，smoke 断言已同步补齐。
- `E2` 已启动第一批增强：`SurfaceMaterialModel` 已补 `Debye length / simple sheath / OML-like / Maxwellian weight` 命名 helper，并新增对应单测，现阶段仍属于增强显式性而非改动主求解语义。
- `E3` 已启动第一步拆分：`DensePlasmaSurfaceCharging::exportResults` 中的 instrument observability metadata 组装已迁移到 `SurfaceInstrumentObservability` 辅助单元，先从低风险导出路径拆分开始。
- `E3` 已继续拆第二步：`monitor / shared_runtime_observer` artifact 的导出编排已迁移到 `SurfaceObserverArtifactExport`，主文件保留数据收集与 writer 回调，进一步降低 `exportResults` 的职责密度。
- `E3` 已继续拆第三步：`shared_runtime_consistency / particle_transport_domain / global_particle_domain / global_particle_repository / global_sheath_field_solve` 这组 shared-runtime artifact 的导出编排已并入统一 helper，`exportResults` 进一步收缩为数据准备层。
- `E3` 已继续拆第四步：`graph_report / graph_matrix / field_adapter / boundary_mapping / field_request / field_result / volume_stub / projection` 这组 surface-bridge artifact 的导出编排已迁移到统一 helper，`exportResults` 只保留导出数据和 writer 回调。
- `E3` 已继续拆第五步：`volumetric_adapter / volume_request / volume_result_template / volume_result / volume_history` 这组 volumetric artifact 的导出编排已迁移到统一 helper，`exportResults` 的主流程导出编排基本收口完成。
- `E3` 已继续拆第六步：`writeExternalFieldSolveResultTemplateJson / writeExternalFieldSolveResultJson / writeExternalFieldSolveRequestJson / writeSurfaceVolumeProjectionJson` 的 writer 定义本体已迁移到 `SurfaceBridgeArtifactWriters` 单元，开始从“只抽导出编排”进入“削减巨型实现文件中的 writer 实现本体”阶段。
- `E3` 已继续拆第七步：`writeExternalVolumeMeshStubJson / writeVolumeMeshSkeletonJson / writeExternalFieldBridgeManifestJson` 的 writer 定义本体也已迁移到 `SurfaceBridgeArtifactWriters`，`DensePlasmaSurfaceCharging.cpp` 进一步收缩，`surface-bridge / volume-mesh-stub` 这组导出逻辑已经从 orchestration 拆分推进到 writer-implementation 拆分。
- `E3` 已继续拆第八步：`writeVolumeHistoryJson / writeVolumetricSolverAdapterContractJson / writeExternalVolumeSolveRequestJson / writeExternalVolumeSolveResultTemplateJson / writeExternalVolumeSolveResultJson` 这一组 volumetric writer 定义本体已迁移到 `SurfaceBridgeArtifactWriters`，`DensePlasmaSurfaceCharging.cpp` 中 volume bridge 的请求/结果/history 序列化已基本脱离主文件。
- `E3` 已继续拆第九步：`writeSurfaceGraphMatrixSnapshotCsv / writeSurfaceGraphMatrixSnapshotJson / writeSurfaceFieldSolverAdapterContractReport / writeSurfaceFieldSolverAdapterContractJson / writeSurfaceBoundaryMappingJson` 已迁移到 `SurfaceBridgeArtifactWriters`，`DensePlasmaSurfaceCharging.cpp` 又收掉一批纯合同/矩阵/边界映射导出实现。
- `E3` 已继续拆第十步：`writeSurfaceGraphReport` 已迁移到 `SurfaceBridgeArtifactWriters`，并把 graph report 依赖的 patch override / node owner / node material / branch topology 描述 helper 一并局部下沉到 writer 单元；`DensePlasmaSurfaceCharging.cpp` 再减少一段高体量文本报告导出实现。
- `E3` 已继续拆第十一步：`writeSharedSurfaceRuntimeObserverJson` 与 `writeSharedSurfaceRuntimeConsistencyJson` 已迁移到 `SurfaceBridgeArtifactWriters`；其中 `observer` 的旧实现已从主文件移除，`consistency` 的旧实现已从主文件编译路径摘除，shared-runtime artifact 开始进入独立 writer 单元收口阶段。
- `E3` 已继续拆第十二步：`writeSharedSurfaceParticleTransportDomainJson / writeGlobalSurfaceParticleDomainJson / writeGlobalParticleRepositoryJson / writeGlobalSurfaceSheathFieldSolveJson` 已全部迁移到 `SurfaceBridgeArtifactWriters`；`DensePlasmaSurfaceCharging.cpp` 中对应旧实现已整体从编译路径摘除，shared/global particle + sheath field 这组 writer 收口完成，并已通过 `SurfaceChargingSmokeTest.SurfaceReferenceAndRuntimeContractsExportMetadata` 验证。
- `E3` 已继续拆第十三步：`SurfaceGraphMonitorSnapshot / SurfaceFieldMonitorSnapshot / SurfaceBenchmarkMonitorSnapshot` 与 `writeSurfaceMonitorJson` 已提升到 `SurfaceBridgeArtifactWriters` 接口层；`DensePlasmaSurfaceCharging.cpp` 只保留 monitor snapshot 的构建逻辑，monitor artifact 序列化正式脱离主文件，并已通过 `SurfaceChargingSmokeTest.SurfaceReferenceAndRuntimeContractsExportMetadata` 验证。
- `E3` 已继续拆第十四步：`exportResults` 中内联构造 `SurfaceBenchmarkMonitorSnapshot` 的逻辑已提炼为独立 helper，`DensePlasmaSurfaceCharging.cpp` 的 monitor/export orchestration 再缩短一段；该步未改变导出合同，仅减少主流程内联装配代码，并已通过 `SurfaceChargingSmokeTest.SurfaceReferenceAndRuntimeContractsExportMetadata` 验证。
- `E3` 已继续拆第十五步：legacy benchmark sidecar 的 `report/csv` 写出路径与错误收口已提炼为 `writeLegacyBenchmarkArtifacts(...)` helper，`exportResults` 中 benchmark reporting 内联流程继续缩短；该步仅收口 orchestration，不改变 benchmark artifact 合同，并已通过 `SurfaceChargingSmokeTest.SurfaceReferenceAndRuntimeContractsExportMetadata` 验证。
- `E3` 已继续拆第十六步：`exportResults` 中 `data_set` 导出前的 scalar series 规范化/校验已统一收口为 `normalizeExportScalarSeries(...) + validateExportScalarSeries(...)`，同时 legacy benchmark display curve 的初始样本补齐流程已成组提炼为 `.cpp` 内 free helper，`exportResults` 主流程一次性缩掉了导出前校验与 benchmark curve 预处理两段内联逻辑；该批修改不改变 CSV/sidecar 合同，并已通过 `SurfaceChargingSmokeTest.SurfaceReferenceAndRuntimeContractsExportMetadata` 验证。
- `E3` 已继续拆第十七步：`benchmark_case` artifact 的 `case_id / relative RMSE / tolerance profile / validation metrics / reference datasets / json 写出收口` 已整体提炼为 `.cpp` 内 helper（`BenchmarkCaseArtifactInputs` + `writeSurfaceBenchmarkCaseArtifact(...)` 及其派生函数），`exportResults` 末尾不再内联堆叠整段 benchmark-case 组装逻辑；该步不改变 benchmark case 合同，并已通过 `SurfaceChargingSmokeTest.SurfaceReferenceAndRuntimeContractsExportMetadata` 验证。

## 3. 改进总原则

- 新的通用数值能力优先进入 `Tools/*`，不得默认继续堆叠到 `Toolkit/Surface Charging`。
- `Toolkit/Surface Charging` 只保留：
  - 场景装配
  - 兼容适配
  - 产物导出
  - route/preset/runtime 组织
- 任何“最小插件层”“最小闭环”“bridge 可用”“gate 通过”都不得直接当作“与 SPIS family 已 1:1 对齐”的证据。
- family 对齐优先于口径乐观化；若实现未完成，应明确记录为缺口，而不是通过改写定义掩盖缺口。
- family sealoff 与 class-level 显式对象镜像是两个不同层级：前者用于回答“当前主线功能是否已达成 1:1”，后者用于回答“逐对象职责是否已可审计”。后续不得把 class-level `近似实现` 误写成 family 未封板，也不得把 family 已封板误写成所有 class 都已文件级 1:1 镜像。
- 历史 round 记录保留事实，但不能替代新的 1:1 对齐审计。
- 每个新增 family 最终都应沉淀到：
  - 代码落点
  - 配置/路由落点
  - 测试落点
  - 文档状态落点

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
| `Round 12` | Top Transition 全量对标轮 | `已完成` | `2026-04-14` | 见 4.11 节 `R12-A*` 清单；`R12-A1` 到 `R12-A8` 已完成；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 进入 `Round 13`，推进 Surf distrib 全量对标 |
| `Round 13` | Surf Distrib 全量对标轮 | `已完成` | `2026-04-14` | 见 4.11 节 `R13-B*` 清单；`R13-B1` 到 `R13-B7` 已完成；`Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 进入 `Round 14`，推进 Surf scaler 全量对标 |
| `Round 14` | Surf Scaler 全量对标轮 | `已完成` | `2026-04-14` | 见 4.11 节 `R14-C*` 清单；`R14-C1` 到 `R14-C6` 已完成；`Tools/FieldSolver/include/SurfaceBarrierModels.h`；`Tools/FieldSolver/src/SurfaceBarrierModels.cpp`；`Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 进入 `Round 15`，推进 DIDV 高阶组合 |
| `Round 15` | Circ DIDV 高阶组合对标轮 | `已完成` | `2026-04-14` | 见 4.11 节 `R15-D*` 清单；`R15-D1` 到 `R15-D4` 已完成；`Tools/Coupling/include/SurfaceCircuitCoupling.h`；`Tools/Coupling/src/SurfaceCircuitCoupling.cpp`；`Tools/Coupling/test/Coupling_test.cpp` | 进入 `Round 16`，推进材料参数资产与 FN 参数链复核 |
| `Round 16` | Material 参数资产全量对标轮 | `已完成` | `2026-04-14` | 见 4.11 节 `R16-E*` 清单；`R16-E1` 到 `R16-E4` 已完成；`Tools/Material/include/SurfaceMaterialModel.h`；`Tools/Material/src/SurfaceMaterialModel.cpp`；`Tools/Material/test/Material_test.cpp`；`Tools/FieldSolver/src/SurfaceBarrierModels.cpp`；`Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp` | 进入 `Round 17`，推进 native volume parity 子计划 |
| `Round 17` | Vol Native Parity 轨道轮 | `已完成` | `2026-04-14` | 见 4.11 节 `R17-F*` 清单；`R17-F1` 到 `R17-F5` 已完成；`Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/FieldSolver/include/SurfaceFieldVolumeBridge.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Tools/FieldSolver/test/SurfaceFieldVolumeBridge_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 进入 `Round 18`，推进后续 native volume 扩展或收口 |
| `Round 18` | 全量回归封板与文档闭环轮 | `已完成` | `2026-04-14` | 见 4.11 节 `R18-G*` 清单；`R18-G1` 到 `R18-G5` 已完成；`build/surface_reference_matrix_gate_subset/surface_reference_matrix_gate.json`；`build/surface_reference_matrix_gate_full/surface_reference_matrix_gate.json`；`build/config_schema_v1_gate/config_schema_v1_gate.json`；`build/kernel_contract_catalog_gate.json`；`build/benchmark_case_matrix_gate.json`；`Toolkit/Surface Charging/docs/surface_round10_equivalence_matrix.md` | Surface 主线完成封板，可按 2.11 节口径收口 |
| `Round 19` | 标准B强约束执行轮（SurfInteract/Axisym/Top 生命周期/Vol phase-2/gate 扩展） | `已完成` | `2026-04-14` | 本文档 2.12-2.13 节；`Tools/Material/include/SurfaceInteraction.h`；`Tools/Particle/include/SurfaceDistributionFunction.h`；`Main/SurfaceScenarioLoader.cpp`；`Main/main.cpp`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build/surface_reference_matrix_gate_full/surface_reference_matrix_gate.json`；`build/config_schema_v1_gate/config_schema_v1_gate.json`；`build/kernel_contract_catalog_gate.json`；`build/benchmark_case_matrix_gate.json` | 标准B硬门禁达成，转入持续回归监测 |
| `Round 27` | SurfInteract 口径校正与显式实现轮 | `已完成` | `2026-04-16` | `Tools/Material/include/SurfaceInteraction.h`；`Tools/Material/src/SurfaceInteraction.cpp`；`Tools/Material/test/Material_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_reference_matrix_gate_round27/surface_reference_matrix_gate.json`；`build_codex/surface_round26_sealoff/surface_round26_sealoff.json` | 进入 `Round 28`，推进 `Vol BC/Field` 显式实现 |
| `Round 28` | Vol BC/Field 显式实现轮 | `已完成` | `2026-04-16` | `Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Tools/FieldSolver/test/SurfaceFieldVolumeBridge_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `Round 29`（已完成），并汇入 `vol_full_stack` 最终封板 |
| `Round 29` | Vol Distrib 显式实现轮 | `已完成` | `2026-04-16` | `Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Tools/FieldSolver/test/SurfaceFieldVolumeBridge_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `Round 30`（已完成），推进 `Vol Interact` 显式实现 |
| `Round 30` | Vol Interact 显式实现轮 | `已完成` | `2026-04-16` | 本文档 2.10、2.10.4、4.14 节；`R23-L1` 明细表；`Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/src/VolInteract.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_family_coverage_gate_round30/surface_family_coverage_gate.json` | 进入 `Round 31`，推进 `near-equivalent` 家族提升；`vol_full_stack` 暂维持 `not_fully_implemented` |
| `Round 31` | `SurfField + Top + Circ` near-equivalent 提升轮 | `已完成` | `2026-04-16` | `Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Tools/Coupling/test/Coupling_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json`；`scripts/run/surface_spis_family_coverage_matrix.json` | 进入 `Round 32`（已完成），并转入 `Round 33` 重封板复核 |
| `Round 32` | `SurfMaterial + Solver` near-equivalent 提升轮 | `已完成` | `2026-04-16` | `Tools/Material/test/Material_test.cpp`；`Tools/Solver/test/SurfaceSolverFacade_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json`；`scripts/run/surface_spis_family_coverage_matrix.json` | 进入 `Round 33`（已完成），执行最终重封板复核 |
| `Round 33` | 最终重封板轮 | `已完成` | `2026-04-16` | `build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json`；`build_codex/surface_round33_sealoff/surface_round26_sealoff.json`；本文档 2.12、4.14 节 | `final_alignment_verdict = achieved`；项目转入长期回归维护 |
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
| `Round 12` | `R12-A5` | 实现 `SheathOrPresheathPoissonBCUpdater` 事件 | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp` | `已完成` | `2026-04-14` | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 鞘层边界策略可切换并可回归 |
| `Round 12` | `R12-A6` | 实现 `RCCabsSCUpdater` / `VcrossBfieldUpdater` / `BasicEclipseExit` / `TransientArtificialSources` / `LangmuirProbeTransition` | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp` | `已完成` | `2026-04-14` | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 事件可独立启停并进入执行链 |
| `Round 12` | `R12-A7` | 补齐对象层事件测试 | `Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp` | `已完成` | `2026-04-14` | `Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`SurfaceTransitionEngineTest.EvaluatesExtendedEventBundleWhenSunFluxIsDisabled`；`SurfaceTransitionEngineTest.ExtendedEventBundleDefaultsToDisabled`；`SurfaceTransitionEngineTest.EnablesSelectedExtendedEventsIndependently` | 新增事件用例通过 |
| `Round 12` | `R12-A8` | 补齐事件 smoke 与 sidecar 可观测字段 | `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | `2026-04-14` | `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`SurfaceChargingSmokeTest.SheathOrPresheathPoissonBCUpdaterEventSwitchesBoundaryScaleWhenEnabled`；`SurfaceChargingSmokeTest.ExtendedTransitionEventsA6CanBeEnabledIndependently`；`SurfaceChargingSmokeTest.SurfaceReferenceAndRuntimeContractsExportMetadata` | smoke 通过且导出字段可审计 |
| `Round 13` | `R13-B1` | 实现 `MultipleSurfDistrib` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `已完成` | `2026-04-14` | `Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`SurfaceDistributionFunctionTest.MultipleSurfDistributionCombinesComponentMoments`；`SurfaceDistributionFunctionTest.BuildAndSynthesisSupportMultipleSurfAndLocalModifiedPearsonModels`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`SurfaceChargingSmokeTest.DistributionModelSwitchRebalancesReferencePatchCurrents` | 分布族可路由且默认兼容 |
| `Round 13` | `R13-B2` | 实现 `LocalModifiedPearsonIV` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `已完成` | `2026-04-14` | `Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`SurfaceDistributionFunctionTest.LocalModifiedPearsonDistributionProducesFiniteShapeAwareMoments`；`SurfaceDistributionFunctionTest.BuildAndSynthesisSupportMultipleSurfAndLocalModifiedPearsonModels`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`SurfaceChargingSmokeTest.DistributionModelSwitchRebalancesReferencePatchCurrents` | 可配置、可回归、可解释 |
| `Round 13` | `R13-B3` | 实现 `LocalTabulatedSurfDistrib` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `已完成` | `2026-04-14` | `Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`SurfaceDistributionFunctionTest.LocalTabulatedDistributionAppliesLocalShiftAndScale`；`SurfaceDistributionFunctionTest.BuildAndSynthesisSupportLocalTabulatedAndTwoAxesVelocityModels`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`SurfaceChargingSmokeTest.DistributionModelSwitchRebalancesReferencePatchCurrents` | 局地制表分布可用 |
| `Round 13` | `R13-B4` | 实现 `TwoAxesTabulatedVelocitySurfDistrib` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `已完成` | `2026-04-14` | `Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`SurfaceDistributionFunctionTest.TwoAxesTabulatedVelocityDistributionCombinesNormalAndTangentialEnergy`；`SurfaceDistributionFunctionTest.BuildAndSynthesisSupportLocalTabulatedAndTwoAxesVelocityModels`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`SurfaceChargingSmokeTest.DistributionModelSwitchRebalancesReferencePatchCurrents` | 双轴速度制表分布可用 |
| `Round 13` | `R13-B5` | 实现 `UniformVelocitySurfDistrib` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `已完成` | `2026-04-14` | `Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`SurfaceDistributionFunctionTest.UniformVelocityDistributionProducesSingleModeMoments`；`SurfaceDistributionFunctionTest.BuildAndSynthesisSupportUniformFluidAndFowlerNordheimModels`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`SurfaceChargingSmokeTest.DistributionModelSwitchRebalancesReferencePatchCurrents` | 均匀速度分布可用 |
| `Round 13` | `R13-B6` | 实现 `FluidSurfDistrib` / `FowlerNordheimSurfDistrib` | `Tools/Particle/src/SurfaceDistributionFunction.cpp` | `已完成` | `2026-04-14` | `Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`SurfaceDistributionFunctionTest.FluidDistributionTracksDirectedFlowAndCompressibility`；`SurfaceDistributionFunctionTest.FowlerNordheimDistributionProducesFiniteEmissionMoments`；`SurfaceDistributionFunctionTest.BuildAndSynthesisSupportUniformFluidAndFowlerNordheimModels`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`SurfaceChargingSmokeTest.DistributionModelSwitchRebalancesReferencePatchCurrents` | 特化分布族可路由与回归 |
| `Round 13` | `R13-B7` | 补齐分布族单测与 smoke | `Tools/Particle/test/SurfaceDistributionFunction_test.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | `已完成` | `2026-04-14` | `Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`SurfaceDistributionFunctionTest.BuildAndSynthesisSupportMultipleSurfAndLocalModifiedPearsonModels`；`SurfaceDistributionFunctionTest.BuildAndSynthesisSupportLocalTabulatedAndTwoAxesVelocityModels`；`SurfaceDistributionFunctionTest.BuildAndSynthesisSupportUniformFluidAndFowlerNordheimModels`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`SurfaceChargingSmokeTest.DistributionModelSwitchRebalancesReferencePatchCurrents`；`build/bin/SurfaceDistributionFunction_test.exe`；`build/bin/SurfaceCharging_smoke_test.exe` | 每个新增分布族具备单测证据 |
| `Round 14` | `R14-C1` | 实现 `AutomaticBarrierCurrentScaler` | `Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | `已完成` | `2026-04-14` | `Tools/FieldSolver/include/SurfaceBarrierModels.h`、`Tools/FieldSolver/src/SurfaceBarrierModels.cpp`、`Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`、`SurfaceBarrierModelsTest.AutomaticBarrierCurrentScalerSelectsExpectedInnerFamilies` | 自动选择规则可复现 |
| `Round 14` | `R14-C2` | 实现 `MultipleCurrentScaler` | `Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | `已完成` | `2026-04-14` | `Tools/FieldSolver/include/SurfaceBarrierModels.h`、`Tools/FieldSolver/src/SurfaceBarrierModels.cpp`、`Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`、`SurfaceBarrierModelsTest.MultipleCurrentScalerBlendsConfiguredInnerScalers` | 组合策略可配置 |
| `Round 14` | `R14-C3` | 实现 `GlobalTempCurrentScaler` / `SmoothedGlobalTempCurrentScaler` | `Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | `已完成` | `2026-04-14` | `Tools/FieldSolver/include/SurfaceBarrierModels.h`、`Tools/FieldSolver/src/SurfaceBarrierModels.cpp`、`Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`、`SurfaceBarrierModelsTest.GlobalTempCurrentScalerUsesConfiguredGlobalTemperature`、`SurfaceBarrierModelsTest.SmoothedGlobalTempCurrentScalerBlendsLocalAndGlobalTemperature` | 温度族 scaler 可回归 |
| `Round 14` | `R14-C4` | 实现显式 `CurrentVariation` / `LocalVariation` 策略类 | `Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | `已完成` | `2026-04-14` | `Tools/FieldSolver/include/SurfaceBarrierModels.h`、`Tools/FieldSolver/src/SurfaceBarrierModels.cpp`、`Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`、`SurfaceBarrierModelsTest.CurrentVariationScalerUsesConfiguredGlobalReferenceState`、`SurfaceBarrierModelsTest.LocalVariationScalerTracksLocalReferenceState` | 显式策略与现有线性化兼容 |
| `Round 14` | `R14-C5` | 导出自动/组合 scaler 决策元数据 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | `2026-04-14` | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsUnifiedKernelMetadata`、`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsAutomaticAndMultipleBarrierMetadata` | 决策链可观测可追溯 |
| `Round 14` | `R14-C6` | 补齐 scaler 单测与 smoke | `Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | `已完成` | `2026-04-14` | `Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`build/bin/SurfaceBarrierModels_test.exe`、`build/bin/SurfaceCharging_smoke_test.exe` | 新增 scaler 用例通过 |
| `Round 15` | `R15-D1` | 实现 `ReducableDIDV` 基础框架 | `Tools/Coupling/include/SurfaceCircuitCoupling.h`、`Tools/Coupling/src/SurfaceCircuitCoupling.cpp` | `已完成` | `2026-04-14` | `Tools/Coupling/include/SurfaceCircuitCoupling.h`、`Tools/Coupling/src/SurfaceCircuitCoupling.cpp`、`Tools/Coupling/test/Coupling_test.cpp`、`CouplingTest.ReducableDidvFrameworkKeepsComposedKernelInputUsable` | 可规约组合与现有聚合共存 |
| `Round 15` | `R15-D2` | 实现 `RedDIDVFromSurfDIDV` / `RedDIDVfromRegDIDV` | `Tools/Coupling/src/SurfaceCircuitCoupling.cpp` | `已完成` | `2026-04-14` | `Tools/Coupling/include/SurfaceCircuitCoupling.h`、`Tools/Coupling/src/SurfaceCircuitCoupling.cpp`、`Tools/Coupling/test/Coupling_test.cpp`、`CouplingTest.ReducedRegularDidvAggregatesMappedNodesAndInternalCouplings`、`CouplingTest.ReducedSurfDidvAggregatesMappedSurfaceMatrices` | 规约路径可运行 |
| `Round 15` | `R15-D3` | 实现 `DIDVfromSurfDIDV` / `SurfDIDVFromMatrices` 显式路由 | `Tools/Coupling/src/SurfaceCircuitCoupling.cpp` | `已完成` | `2026-04-14` | `Tools/Coupling/include/SurfaceCircuitCoupling.h`、`Tools/Coupling/src/SurfaceCircuitCoupling.cpp`、`Tools/Coupling/test/Coupling_test.cpp`、`CouplingTest.SurfDidvFromMatricesBuildsExplicitNodeAndOffDiagonalRoute` | DIDV 构造链完整 |
| `Round 15` | `R15-D4` | 补齐 DIDV 高阶组合测试 | `Tools/Coupling/test/Coupling_test.cpp` | `已完成` | `2026-04-14` | `Tools/Coupling/test/Coupling_test.cpp`、`build/bin/Coupling_test.exe` | 高阶 DIDV 单测 + smoke 通过 |
| `Round 16` | `R16-E1` | 补齐 `ErosionParamSet` 对应参数资产 | `Tools/Material/src/SurfaceMaterialModel.cpp` | `已完成` | `2026-04-14 15:24` | `resolveErosionParamSet(...)` 已补齐 kapton / ptfe(teflon) / aluminum 默认参数与 scalar override 链 | GEO/LEO 主线材料覆盖完整 |
| `Round 16` | `R16-E2` | 补齐背散射二维表（能量 x 角度）路径 | `Tools/Material/src/SurfaceMaterialModel.cpp` | `已完成` | `2026-04-14 15:24` | `backscatter_table_energy_count` / `backscatter_table_angle_count` / `backscatter_table_value_i_j` 双轴插值路径已接通，并保留 legacy fallback | 二维表驱动生效且兼容 |
| `Round 16` | `R16-E3` | 复核 FN 参数链全路径一致性 | `Tools/FieldSolver/include/SurfaceBarrierModels.h`、`Tools/FieldSolver/src/SurfaceBarrierModels.cpp` | `已完成` | `2026-04-14 15:24` | `resolveSurfaceFieldEmissionParameters(...)` 与 `evaluateFowlerNordheimBarrier(...)` 已统一到同一组 FN 键名/默认值/threshold/sheath-floor 语义 | 参数从配置到求解一致 |
| `Round 16` | `R16-E4` | 补齐材料参数回归测试 | `Tools/Material/test/Material_test.cpp` | `已完成` | `2026-04-14 15:24` | `Material_test.exe` 34/34 通过；`SurfaceBarrierModels_test.exe` 20/20 通过；新增 erosion/defaults、2D backscatter、FN alias、FN barrier regression 用例 | 材料参数测试与场景 smoke 通过 |
| `Round 17` | `R17-F1` | 建立 native volume parity 子计划与接口边界 | `Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h` | `已完成` | `2026-04-14 15:32` | 已新增 `NativeVolumeParityRoute`、mesh/field/distribution/interaction/snapshot contracts 与 route name 约定 | native parity 设计评审通过 |
| `Round 17` | `R17-F2` | 实现最小 `VolMesh` native 路径 | `Tools/Mesh/*` | `已完成` | `2026-04-14 15:32` | `summarizeNativeVolumeMesh(...)` 基于 `Mesh::VolMesh` 输出 node/element/bbox/volume 摘要，并在 native parity snapshot 中闭环 | 原生体网格最小能力可运行 |
| `Round 17` | `R17-F3` | 实现最小 `VolField` native 路径 | `Tools/FieldSolver/*` | `已完成` | `2026-04-14 15:32` | `solveNativeVolumePotential(...)` 已复用 `PoissonSolver` 建立最小 Dirichlet 体场求解路径 | 原生体场步进可运行 |
| `Round 17` | `R17-F4` | 实现最小 `VolDistrib` / `VolInteract` native 路径 | `Tools/Particle/*` | `已完成` | `2026-04-14 15:32` | `evaluateNativeVolumeDistributions(...)` 与 `executeNativeVolumeInteractions(...)` 已接通 `VolDistrib` / `VolInteract` 最小链路 | 体分布与体交互可运行 |
| `Round 17` | `R17-F5` | 完成 bridge 与 native 并存兼容 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已完成` | `2026-04-14 15:32` | metadata 已导出 `surface_external_volume_solver_bridge_enabled`、`surface_native_volume_parity_mode`、`surface_native_volume_parity_bridge_compatible`，并通过 smoke | 双路径并存且不回归 |
| `Round 18` | `R18-G1` | 执行关键单测集合封板 | `build/*` | `已完成` | `2026-04-14 15:40` | `Material_test.exe` 34/34、`SurfaceDistributionFunction_test.exe` 37/37、`SurfaceBarrierModels_test.exe` 20/20、`Coupling_test.exe` 39/39、`SurfaceFieldVolumeBridge_test.exe` 5/5、关键 `SurfaceCharging_smoke_test` 3/3 全部通过 | 关键测试集全绿 |
| `Round 18` | `R18-G2` | 执行子集矩阵 gate 封板 | `scripts/run/surface_reference_matrix_round11_subset.json` | `已完成` | `2026-04-14 15:40` | `build/surface_reference_matrix_gate_subset/surface_reference_matrix_gate.json` 状态 `PASS`；2/2 case `PASS` | 子集 gate `PASS` |
| `Round 18` | `R18-G3` | 执行全量矩阵 gate 封板 | `scripts/run/surface_reference_matrix.json` | `已完成` | `2026-04-14 15:40` | `build/surface_reference_matrix_gate_full/surface_reference_matrix_gate.json` 状态 `PASS`；7/7 case `PASS` | 全量 gate `PASS` |
| `Round 18` | `R18-G4` | 执行 contract / schema 校验 | `documentation/contracts/config_schema_v1.json` | `已完成` | `2026-04-14 15:40` | `build/config_schema_v1_gate/config_schema_v1_gate.json`、`build/kernel_contract_catalog_gate.json`、`build/benchmark_case_matrix_gate.json` 均为 `PASS` | schema 与输出契约一致 |
| `Round 18` | `R18-G5` | 文档闭环与状态回填封板 | `Toolkit/Surface Charging/docs/surface_num_alignment_refactor_plan.md`、`Toolkit/Surface Charging/docs/surface_round10_equivalence_matrix.md` | `已完成` | `2026-04-14 15:40` | 本文档 `Round 18` 与 `R18-G*` 状态回填完成；等价矩阵补充 Round 18 封板说明 | 任务状态、证据、完成时间全部回填 |

### 4.12 Round 19：标准B强约束执行轮（新增，锁定版）

| 字段 | 内容 |
| --- | --- |
| 目标 | 在不回退 `Round 18` 封板结论的前提下，完成标准B（结构同构级）所需的全部剩余实现（SurfInteract、Axisym 显式族、Top 生命周期对象层、Vol phase-2、长尾 gate 覆盖） |
| 输入 | 本文档 2.12 节复核结论；`ref/SPIS/num-master/src/main/java/spis/Surf/SurfInteract/**`；`ref/SPIS/num-master/src/main/java/spis/Top/Transition/**`；`ref/SPIS/num-master/src/main/java/spis/Vol/**`；`scripts/run/surface_reference_matrix.json` |
| 实施项 | 对 `R19-H1`~`R19-H5` 逐项强制实现并补齐测试与 gate；禁止“冻结差异”结题路径 |
| 输出 | 标准B达成证据链（代码/测试/gate/文档）与 2.9 对齐状态回填 |
| 依赖 | `Round 18` |
| 验收 | 2.13 节全部硬门禁满足；且主线 `surface_reference_matrix` 与 contract gate 不回归 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-14` |
| 证据 | 本文档 2.12-2.13 节；`build/surface_reference_matrix_gate_full/surface_reference_matrix_gate.json`；`build/config_schema_v1_gate/config_schema_v1_gate.json`；`build/kernel_contract_catalog_gate.json`；`build/benchmark_case_matrix_gate.json` |

#### Round 19 任务卡

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R19-H1` | 实现最小 `SurfInteract` 插件层（强制） | `已完成` | `2026-04-14` | `Round 18` | 独立 `SurfInteract` 最小对象层接口 + 主线路由接入 + 默认兼容 | 代码检视 + Material/FieldSolver 单测 + smoke | `Tools/Material/include/SurfaceInteraction.h`；`Tools/Material/src/SurfaceInteraction.cpp`；`Tools/Material/test/Material_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`SurfaceInteractionTest.PluggableInteractorOverridesEvaluationAndCanResetToDefault`；`SurfaceChargingSmokeTest.SurfaceReferenceAndRuntimeContractsExportMetadata` | 进入 `R19-H2`，实现 axisym 显式族 |
| `R19-H2` | 实现 `AxisymTabulatedVelocity` 显式族并接入主线（强制） | `已完成` | `2026-04-14` | `Round 18` | axisym 显式族、参数约束、主线分布路由与元数据导出 | Particle 单测 + 分布切换 smoke + matrix gate | `Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Main/SurfaceScenarioLoader.cpp`；`Main/main.cpp`；`build/surface_reference_matrix_gate_full/surface_reference_matrix_gate.json` | 转入回归监测 |
| `R19-H3` | 实现 Top 生命周期对象层最小钩子（强制） | `已完成` | `2026-04-14` | `Round 12` | `TransitionObserver`/`Finalization`/`SunFluxIntensity` 最小对象层接口 + 运行时可观测字段 | object-layer test + smoke + sidecar 可观测字段复核 | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 转入回归监测 |
| `R19-H4` | 实现 Vol parity phase-2 最小闭环（强制） | `已完成` | `2026-04-14` | `Round 17` | 多 BC / 多分布 / 多交互最小闭环与桥接兼容报告 | FieldSolver/Particle 单测 + volume 相关 smoke + matrix gate | `Tools/FieldSolver/include/SurfaceFieldVolumeBridge.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Tools/FieldSolver/test/SurfaceFieldVolumeBridge_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build/surface_reference_matrix_gate_full/surface_reference_matrix_gate.json` | 转入回归监测 |
| `R19-H5` | 扩展 `surface_reference_matrix` 长尾场景并补 gate（强制） | `已完成` | `2026-04-14` | `Round 18` | 扩展后的矩阵 case（`12` 个）与 gate 报告 | matrix gate + contract gate | `scripts/run/surface_reference_matrix.json`；`build/surface_reference_matrix_gate_full/surface_reference_matrix_gate.json`；`build/config_schema_v1_gate/config_schema_v1_gate.json`；`build/kernel_contract_catalog_gate.json`；`build/benchmark_case_matrix_gate.json` | Round 19 封板完成 |

#### Round 19 执行顺序（一次性细化，禁止反复改口径）

1. `R19-H1`：先补 `SurfInteract` 最小对象层，确保主线兼容不回归。
2. `R19-H2`：接入 axisym 显式族与参数约束，完成分布族语义闭环。
3. `R19-H3`：补齐 Top 生命周期对象层钩子与可观测字段。
4. `R19-H4`：推进 Vol parity phase-2 最小闭环（多 BC / 多分布 / 多交互）。
5. `R19-H5`：最后扩展长尾矩阵并执行全量 gate 封板。

- 执行锁定规则：除非出现阻断式编译/契约冲突，不得新增“冻结”分支；若发生阻断，必须先修复阻断后继续执行，不得降级标准B口径。

### 4.13 Round 20-26：SPIS 功能 1:1 对齐执行轮

#### Round 20：Wave A 文档重构与全量差异盘点封板

| 字段 | 内容 |
| --- | --- |
| 目标 | 将本文档从“主线最小集对齐”口径重构为“与 `ref/SPIS/num-master` 做功能 1:1 对齐”的总控文档，并完成首轮全量 family 盘点 |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/**`；本文档当前版本；当前仓库对应模块实现 |
| 实施项 | 执行 `R20-I1`~`R20-I5`，重写 1~3 节、补齐 2.9 映射表、建立 family 分类与 backlog |
| 输出 | 1:1 对齐口径文档、全量 family 映射总表、后续任务组与任务 ID |
| 依赖 | `Round 19` 历史事实保留 |
| 验收 | 不存在把“最小闭环”写成“已对齐”的条目；不存在未归类的 SPIS family |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-15` |
| 证据 | 本文档 1~3 节与 4.13 节；`ref/SPIS/num-master/src/main/java/spis/**` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R20-I1` | 重写导言、权威来源与 1:1 判定规则 | `已完成` | `2026-04-15` | `Round 19` | 1:1 对齐文档口径 | 文档检视 | 本文档 1.1~1.5 | 进入 `R20-I2` |
| `R20-I2` | 重写六域总览与当前审计结论 | `已完成` | `2026-04-15` | `R20-I1` | 2.1~2.8 改为 1:1 差距分析 | 文档检视 | 本文档 2.1~2.8 | 进入 `R20-I3` |
| `R20-I3` | 建立 SPIS -> SCDAT 一级映射总表 | `已完成` | `2026-04-15` | `R20-I2` | 2.9 一级映射表 | 文档检视 | 本文档 2.9 | 进入 `R20-I4` |
| `R20-I4` | 建立首批 class-family 盘点基线 | `已完成` | `2026-04-15` | `R20-I3` | 2.10 首批 family 清单 | 文档检视 | 本文档 2.10 | 进入 `R20-I5` |
| `R20-I5` | 为未对齐 family 分配后续任务组 | `已完成` | `2026-04-15` | `R20-I4` | `R21-J*`~`R26-O*` 任务骨架 | 文档检视 | 本文档 4.13 节 | 进入 `Round 21` |

#### Round 21：Wave B `SurfInteract` 全家族补齐

| 字段 | 内容 |
| --- | --- |
| 目标 | 让 `Tools/Material` 与 `Toolkit/Surface Charging` 覆盖 `spis/Surf/SurfInteract/**` 的 family 级功能语义，而不是只保留最小插件层 |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfInteract/**`；`Tools/Material/include/SurfaceInteraction.h`；`Tools/Material/src/SurfaceInteraction.cpp` |
| 实施项 | 执行 `R21-J1`~`R21-J6`，建立多 family 选择器、组合器、默认适配器与配置/元数据导出 |
| 输出 | `SurfInteract` family 路由层、配置入口、测试与 smoke 证据 |
| 依赖 | `Round 20` |
| 验收 | SPIS `SurfInteract` family 不再存在未映射项；每个 family 至少具备代码落点与测试落点 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-15` |
| 证据 | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfInteract/**`；`Tools/Material/include/SurfaceInteraction.h`；`Tools/Material/src/SurfaceInteraction.cpp`；`Tools/Material/test/Material_test.cpp`；`Main/main.cpp`；`Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R21-J1` | 建立 `SurfInteract` family-to-family 映射表 | `已完成` | `2026-04-15` | `Round 20` | 映射表与实现优先级 | 文档检视 | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfInteract/**`；本文档 2.10；`Tools/Material/include/SurfaceInteraction.h` | 进入 `R21-J2` |
| `R21-J2` | 将 `SurfaceInteraction` 扩展为多 family 选择器与组合器 | `已完成` | `2026-04-15` | `R21-J1` | family 路由层骨架 | `Material_test` | `Tools/Material/include/SurfaceInteraction.h`；`Tools/Material/src/SurfaceInteraction.cpp`；`SurfaceInteractionTest.ResolvesConfiguredInteractorFamiliesFromMaterialFlagsAndCodes` | 进入 `R21-J3` |
| `R21-J3` | 补齐 `MaterialInteractor`、`GenericDFInteractor`、`MaterialDFInteractor`、`MaxwellianInteractor` 家族 | `已完成` | `2026-04-15` | `R21-J2` | 首批 interactor family + `surface_interactor_family` 配置入口 | 单测 + smoke | `Tools/Material/src/SurfaceInteraction.cpp`；`Tools/Material/test/Material_test.cpp`；`Main/main.cpp`；`Main/SurfaceScenarioLoader.cpp`；`SurfaceInteractionTest.ResolvesConfiguredInteractorFamiliesFromMaterialFlagsAndCodes`；`SurfaceChargingSmokeTest.SurfaceConfigCliAcceptsSurfaceInteractorFamilyOverrideKeys`；`SurfaceChargingSmokeTest.SurfaceConfigCliRejectsUnsupportedSurfaceInteractorFamilyOverride` | 进入 `R21-J4` |
| `R21-J4` | 补齐 `ReflectionInteractor`、`YieldInteractor`、`ErosionInteractor`、`ImprovedPhotoEmInteractor` 家族 | `已完成` | `2026-04-15` | `R21-J3` | 二批 interactor family | 单测 + smoke | `Tools/Material/src/SurfaceInteraction.cpp`；`Tools/Material/test/Material_test.cpp`；`SurfaceInteractionTest.ExtendedInteractorFamiliesResolveCanonicalMetadata`；`SurfaceChargingSmokeTest.SurfaceReferenceAndRuntimeContractsExportMetadata` | 进入 `R21-J5` |
| `R21-J5` | 补齐 `MultipleInteractor`、`MultipleMaxwellianInteractor`、`MaxwellianInteractorWithRecollection`、`RecollectedSEYModel`、`TabulatedSEYModel` | `已完成` | `2026-04-15` | `R21-J4` | 组合与回收语义 family | 单测 + smoke | `Tools/Material/src/SurfaceInteraction.cpp`；`SurfaceInteractionTest.BuiltinTabulatedSeyInteractorUsesConfiguredLookupTable`；`SurfaceInteractionTest.BuiltinMultipleInteractorCombinesYieldPhotoAndConductionSemantics`；`SurfaceInteractionTest.BuiltinRecollectedSeyInteractorReducesEmissionWithStrongerBarrier`；`SurfaceInteractionTest.BuiltinMultipleMaxwellianInteractorRespondsToEmissionTemperature`；`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsTabulatedSeyInteractorFamilyMetadata`；`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsMultipleInteractorFamilyMetadata` | 进入 `R21-J6` |
| `R21-J6` | 补齐默认模型族与 `Device` 语义导出 | `已完成` | `2026-04-15` | `R21-J5` | `DefaultPEEModel`、`DefaultSEEEYModel`、`DefaultSEEPModel`、`DefaultErosionModel`、`BasicInducedConductInteractor`、`Device` 映射 | 单测 + sidecar 检查 | `Tools/Material/include/SurfaceInteraction.h`；`Tools/Material/src/SurfaceInteraction.cpp`；`Main/main.cpp`；`Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`SurfaceChargingSmokeTest.SurfaceConfigCliAcceptsSurfaceInteractorFamilyOverrideKeys`；`SurfaceChargingSmokeTest.SurfaceConfigCliAcceptsCanonicalDefaultInteractorFamilyAlias` | 进入 `Round 22` / `R22-K1` |

#### Round 22：Wave C `SurfDistrib` 全家族补齐

| 字段 | 内容 |
| --- | --- |
| 目标 | 完成 `spis/Surf/SurfDistrib/**` 全家族 1:1 对齐，而不是当前主线子集 |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfDistrib/**`；`Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp` |
| 实施项 | 执行 `R22-K1`~`R22-K6`，扩展 `SurfaceDistributionModel`、build/synthesis/environment/role 路由与 family 测试 |
| 输出 | `SurfDistrib` 全家族路由、测试与 matrix gate 证据 |
| 依赖 | `Round 20`，建议串联 `Round 21` 中的 interactor 路径 |
| 验收 | SPIS `SurfDistrib` family 无未映射项；每个新增 family 至少具备单测与主线路由 smoke |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-15` |
| 证据 | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfDistrib/**`；`Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`Main/main.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate_round22/surface_reference_matrix_gate.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R22-K1` | 建立 `SurfDistrib` family-to-family 映射表 | `已完成` | `2026-04-15` | `Round 20` | 映射表与实现顺序 | 文档检视 | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfDistrib/**`；本文档 2.10（Surf Distrib 行） | 进入 `R22-K2` |
| `R22-K2` | 补齐 `GlobalMaxwellBoltzmann*`、`GlobalMaxwellSurfDistrib`、`LocalMaxwell*` family | `已完成` | `2026-04-15` | `R22-K1` | Maxwell family 实现 | `SurfaceDistributionFunction_test` | `Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp` | 进入 `R22-K3` |
| `R22-K3` | 补齐 `RecollMaxwellSurfDistrib` 与回收语义 | `已完成` | `2026-04-15` | `R22-K2` | 回收型分布 family | 单测 + smoke | `Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`scripts/run/surface_reference_matrix.json` | 进入 `R22-K4` |
| `R22-K4` | 补齐 `PICSurfDistrib`、`NonPICSurfDistrib`、`GenericSurfDistrib`、`GlobalSurfDistrib`、`LocalGenericSurfDistrib` | `已完成` | `2026-04-15` | `R22-K3` | PIC/Generic family | 单测 + smoke | `Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`scripts/run/surface_reference_matrix.json` | 进入 `R22-K5` |
| `R22-K5` | 补齐 `TestableSurfDistrib`、`TestableForASurfDistrib`、`MaxwellianThruster` | `已完成` | `2026-04-15` | `R22-K4` | 测试型与特化 family | 单测 + smoke | `Tools/Particle/src/SurfaceDistributionFunction.cpp`；`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`；`scripts/run/surface_reference_matrix.json` | 进入 `R22-K6` |
| `R22-K6` | 将全部 `SurfDistrib` family 接入 matrix gate 与 sidecar 元数据 | `已完成` | `2026-04-15` | `R22-K5` | family coverage gate、route 证据 | matrix gate + smoke | `Main/main.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`（`SurfaceChargingSmokeTest.SurfaceConfigCliAcceptsCanonicalSurfDistribFamilyAlias`）；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate_round22/surface_reference_matrix_gate.json`（`PASS`，`26/26`） | 进入 `Round 23` |

#### Round 23：Wave D `Vol` 全家族补齐

| 字段 | 内容 |
| --- | --- |
| 目标 | 让当前项目从 native parity minimal/phase-2 扩展到 `spis/Vol/**` 功能全覆盖 |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Vol/**`；`Tools/FieldSolver/**`；`Tools/Particle/**`；相关 mesh/boundary 模块 |
| 实施项 | 执行 `R23-L1`~`R23-L7`，补齐 `BC`、`Geom`、`VolMesh`、`VolField`、`VolDistrib`、`VolInteract` 家族 |
| 输出 | `Vol` full stack 路由、实现、测试与 parity/matrix 证据 |
| 依赖 | `Round 20` |
| 验收 | `Vol` family 不再依赖“minimal parity”冒充已对齐；排除项若有必须明确记录 |
| 当前状态 | `已完成（family 可观测）` |
| 最后更新 | `2026-04-15` |
| 证据 | `ref/SPIS/num-master/src/main/java/spis/Vol/**`；`Tools/FieldSolver/include/VolField.h`；`Tools/FieldSolver/include/VolDistrib.h`；`Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/FieldSolver/include/SurfaceFieldVolumeBridge.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate.json`；本文档 `R23-L1` 映射明细 |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R23-L1` | 建立 `Vol` 全族映射表 | `已完成` | `2026-04-15` | `Round 20` | `Vol` family 盘点总表 | 文档检视 | `spis/Vol/**`；`Tools/FieldSolver/include/VolField.h`；`Tools/FieldSolver/include/VolDistrib.h`；`Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；本文档 `R23-L1` 映射明细 | 进入 `R23-L2` |
| `R23-L2` | 补齐 `BC` 家族：`VoltageDependentMBC`、`MixedDirichletFourierPoissonBC`、`FourierPoissonBC`、`SurfDistribMatterBC`、`OneSurfDistribTestableMatterBC`、`CapacitiveVoltageGenerator` | `已完成（family 可观测）` | `2026-04-15` | `R23-L1` | `Vol/BC` family route + signature | `SurfaceFieldVolumeBridge_test` + smoke | `Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 进入 `R23-L3` |
| `R23-L3` | 补齐 `VolField` 家族：`UniformBField`、`SolenoidBField`、`DipolarBField`、`MultipleEField`、`EMField` 组合语义 | `已完成（family 可观测）` | `2026-04-15` | `R23-L2` | `Vol/Field` family route + signature | 单测 + smoke | `Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 进入 `R23-L4` |
| `R23-L4` | 补齐 `VolDistrib` 家族：`PICVolDistrib`、`PICVolDistribNoAcc`、`PICVolDistribUpdatable`、`SmartPICVolDistrib` | `已完成（family 可观测）` | `2026-04-15` | `R23-L3` | PIC volume family route + signature | 单测 + smoke | `Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 进入 `R23-L5` |
| `R23-L5` | 补齐 `CompositeVolDistrib`、`BackTrackingVolDistrib`、`BacktrackingPICCompositeVolDistrib`、`BacktrackingBoltzmannCompositeVolDistrib` | `已完成（family 可观测）` | `2026-04-15` | `R23-L4` | composite/backtracking family route + signature | 单测 + smoke | `Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`scripts/run/surface_reference_matrix.json` | 进入 `R23-L6` |
| `R23-L6` | 补齐 `VolInteract` 家族：`MCCInteractor`、`CEXInteractor`、`PhotoIonization`、`TrajectoryInteractionFromField`、`SpinningSpacecraftTrajectory` | `已完成（显式实现）` | `2026-04-16` | `R23-L5` | interaction family route + signature | 单测 + smoke + matrix gate | `Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/src/VolInteract.cpp`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`scripts/run/surface_reference_matrix.json` | 进入 `R23-L7` |
| `R23-L7` | 将 `Vol` family 纳入 volume 专项 smoke 与 parity/reference gate | `已完成` | `2026-04-15` | `R23-L6` | volume gate 证据链 | gate + smoke + snapshot | `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate.json` | 进入 `Round 24` |

##### `R23-L1` `Vol` 全族映射明细（2026-04-15）

| SPIS `Vol` family | 当前 SCDAT 落点 | 当前可达状态 | 主要缺口 | 对应后续任务 |
| --- | --- | --- | --- | --- |
| `VoltageDependentMBC` | `Tools/FieldSolver/include/BoundaryCondition.h`，`Tools/FieldSolver/include/PoissonSolver.h`，`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp` (`solveNativeVolumePotential`) | `近似可达` | 当前仅有静态 Dirichlet 组装，尚缺显式“电压依赖”边界语义入口 | `R23-L2` |
| `MixedDirichletFourierPoissonBC` | 同上（`PoissonSolver` + `BoundaryCondition`） | `近似可达` | 尚缺 mixed/fourier 家族显式路由与可观测元数据 | `R23-L2` |
| `FourierPoissonBC` | 同上 | `近似可达` | 尚缺 Fourier 家族命名路由与 smoke 级证据 | `R23-L2` |
| `SurfDistribMatterBC` | `Tools/FieldSolver/include/VolDistrib.h` + `SurfaceFieldVolumeBridge.cpp` (`evaluateNativeVolumeDistributions`) | `近似可达` | 边界条件与分布耦合关系尚未显式封装为独立 family | `R23-L2` |
| `OneSurfDistribTestableMatterBC` | 同上 | `近似可达` | 测试型单分布边界语义尚未独立路由 | `R23-L2` |
| `CapacitiveVoltageGenerator` | `SurfaceFieldVolumeBridge.cpp` (`solveNativeVolumePotential`) | `近似可达` | 尚缺显式电容型激励器 family 入口与参数面 | `R23-L2` |
| `UniformBField` | `Tools/FieldSolver/include/VolField.h` (`MagneticField::setUniformField`) | `可达` | 尚未纳入 native volume parity 汇总与 gate 元数据 | `R23-L3` |
| `SolenoidBField` | `VolField` 组合路径（可用轴向场近似） | `近似可达` | 缺显式 solenoid family 封装与测试命名 | `R23-L3` |
| `DipolarBField` | `VolField` (`MagneticField::addDipoleField`) | `可达` | 尚未纳入 parity 报告可观测字段 | `R23-L3` |
| `MultipleEField` | `VolField` (`ElectricField` 叠加) | `可达` | 缺显式多电场组合 family 路由 | `R23-L3` |
| `EMField` 组合语义 | `VolField` (`ElectricField` + `MagneticField`) | `近似可达` | 缺统一组合导出与 smoke 证据链 | `R23-L3` |
| `PICVolDistrib` | `VolDistrib` (`ParticleDensityDistrib`/`ChargeDensityDistrib`) | `近似可达` | 缺命名 family 路由与 PIC 变体区分 | `R23-L4` |
| `PICVolDistribNoAcc` | 同上 | `未显式实现` | 尚无 no-acc 专用语义入口 | `R23-L4` |
| `PICVolDistribUpdatable` | 同上（`updateFromParticles`) | `近似可达` | 缺 updatable family 元数据与回归断言 | `R23-L4` |
| `SmartPICVolDistrib` | 同上 | `未显式实现` | 尚无 smart 策略分支与可观测指标 | `R23-L4` |
| `CompositeVolDistrib` | `VolDistrib` 聚合计算路径（当前在 bridge 内部拼装） | `近似可达` | 缺独立 composite family 封装 | `R23-L5` |
| `BackTrackingVolDistrib` | 暂无独立实现 | `未显式实现` | 缺 backtracking 轨迹追溯分布实现 | `R23-L5` |
| `BacktrackingPICCompositeVolDistrib` | 暂无独立实现 | `未显式实现` | 缺 PIC+composite+backtracking 组合实现 | `R23-L5` |
| `BacktrackingBoltzmannCompositeVolDistrib` | 暂无独立实现 | `未显式实现` | 缺 Boltzmann+composite+backtracking 组合实现 | `R23-L5` |
| `MCCInteractor` | `Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/src/VolInteract.cpp`（`MCCInteractor` 显式类） | `已完成（显式实现）` | 已具备显式 family 类、factory 入口与 active signature | `R30-D1` |
| `CEXInteractor` | `Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/src/VolInteract.cpp`（`CEXInteractor` 显式类） | `已完成（显式实现）` | 已具备显式 family 类、factory 入口与 active signature | `R30-D1` |
| `PhotoIonization` | `Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/src/VolInteract.cpp`（`PhotoIonizationInteractor` 显式类） | `已完成（显式实现）` | 已具备显式 family 类、config route 与 matrix 证据 | `R30-D1` |
| `TrajectoryInteractionFromField` | `Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/src/VolInteract.cpp`（`TrajectoryInteractionFromFieldInteractor` 显式类） | `已完成（显式实现）` | 已具备显式 family 类、active signature 与 unit/smoke 证据 | `R30-D1` |
| `SpinningSpacecraftTrajectory` | `Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/src/VolInteract.cpp`（`SpinningSpacecraftTrajectoryInteractor` 显式类） | `已完成（显式实现）` | 已具备显式 family 类、config route 与 active signature | `R30-D1` |

##### `R23-L2`~`R23-L7` 实施与验证闭环（2026-04-15）

- 已在 `SurfaceFieldVolumeContracts` 与 `SurfaceFieldVolumeBridge` 建立 `BC/Field/Distribution/Interaction` 四类 family 的枚举、解析与 signature 导出。
- `DensePlasmaSurfaceCharging` native-volume sidecar 现已导出四类 family 的 `*_family_count` 与 `*_family_signature` 元数据。
- `SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsNativeVolumeParityMetadata` 已补充四类 family 元数据断言。
- `surface_reference_matrix.json` 的 native-volume `hybrid`/`external` case 已挂接上述新增元数据断言。
- `Round 30` 已在 `VolInteract` 中补齐 `MCCInteractor`、`CEXInteractor`、`PhotoIonizationInteractor`、`TrajectoryInteractionFromFieldInteractor`、`SpinningSpacecraftTrajectoryInteractor` 五个显式交互 family，并通过 `buildNativeVolumeInteractionManager(...)` 接入 `SurfaceFieldVolumeBridge`。
- `DensePlasmaSurfaceCharging` 已新增 `surface_native_volume_active_interaction_route`、`surface_native_volume_active_interaction_family_*` 元数据，以及 `surface_native_volume_interaction_families` / `native_volume_interaction_families` 配置入口。
- `surface_reference_matrix.json` 已新增 `geo_ecss_kapton_native_volume_interaction_override` case，`build_codex/surface_family_coverage_gate_round30/surface_family_coverage_gate.json` 通过。
- 验证闭环：
  - 构建：`cmake --build build_codex --config Debug --target SCDAT SurfaceFieldVolumeBridge_test SurfaceCharging_smoke_test` 通过。
  - 目标测试：`SurfaceFieldVolumeBridgeTest.NativeVolumeParitySnapshotBuildsMinimalClosedLoop`、`SurfaceFieldVolumeBridgeTest.HybridBlendParityAddsPhase2BoundaryAndInteractionCoverage`、`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsNativeVolumeParityMetadata` 全部 `PASS`。
  - 全量 gate：`ctest --test-dir build_codex -C Debug -R '^SurfaceReferenceMatrixGate$' --output-on-failure` 通过；`build_codex/surface_reference_matrix_gate.json` 显示 `status=PASS`、`26/26` case 通过、`0` 失败。

#### Round 24：Wave E `Top` 生命周期与 transition 结构对齐

| 字段 | 内容 |
| --- | --- |
| 目标 | 在保留现有外部接口前提下，使 `Top` 层的 transition/lifecycle 组织方式与 SPIS 更接近并完成明确映射 |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Top/**`；`Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |
| 实施项 | 执行 `R24-M1`~`R24-M5`，对齐事件族、observer/finalization/interface 层级与可观测字段 |
| 输出 | `Top` family 映射、对象层级对齐实现与测试证据 |
| 依赖 | `Round 20` |
| 验收 | `Top` 不再只有“事件存在”，还要有 family 级对象组织与语义证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `ref/SPIS/num-master/src/main/java/spis/Top/Transition/**`；`Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R24-M1` | 建立 `Top/Transition` family 映射表 | `已完成` | `2026-04-16` | `Round 20` | family 映射与层级差异表 | 文档检视 | `ref/SPIS/num-master/src/main/java/spis/Top/Transition/**`；本文档 2.10.1 | 进入 `R24-M2` |
| `R24-M2` | 对齐 `TransitionObserver`、`Finalization`、`SunFluxIntensityUpdater` 的对象层级 | `已完成` | `2026-04-16` | `R24-M1` | lifecycle family | `SurfaceCharging_object_layer_test` | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`SurfaceTransitionEngineTest.ExposesObserverFinalizationAndSunFluxIntensityObjectHooks` | 进入 `R24-M3` |
| `R24-M3` | 复核并对齐 `LocalTimeTransition`、`SpinningSpacecraft`、`SunFluxUpdater`、`ConductivityEvolution` 等事件语义 | `已完成` | `2026-04-16` | `R24-M2` | 主线 transition family | object-layer test + smoke | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`SurfaceTransitionEngineTest.EvaluatesExtendedEventBundleWhenSunFluxIsDisabled` | 进入 `R24-M4` |
| `R24-M4` | 对齐 transition interface / lifecycle 组织关系 | `已完成` | `2026-04-16` | `R24-M3` | 更接近 SPIS 的对象层组织 | object-layer test | `Toolkit/Surface Charging/include/SurfaceTransitionEngine.h`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`SurfaceTransitionEngineTest.ExposesObserverFinalizationAndSunFluxIntensityObjectHooks`；`SurfaceTransitionEngineTest.ExtendedEventBundleDefaultsToDisabled` | 进入 `R24-M5` |
| `R24-M5` | 将 `Top` family 纳入 sidecar 与 route 证据链 | `已完成` | `2026-04-16` | `R24-M4` | metadata 与 smoke 证据 | smoke + sidecar 检查 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsTopTransitionFamilyMetadata` | 进入 `Round 25` |

##### `R24-M1`~`R24-M5` 实施与验证闭环（2026-04-16）
- 已在 `SurfaceTransitionEngine` 中补充 `SurfaceTransitionObjectLayerState`，把 `TransitionInterface / Transition / SimulationParamUpdater` 的对象层签名与 active family 视图固定下来。
- `DensePlasmaSurfaceCharging` metadata sidecar 已新增 `surface_top_transition_route`、`surface_top_transition_*family_count`、`surface_top_transition_*family_signature` 系列字段，形成 `Top` family 的 route 证据链。
- 目标测试：`SurfaceTransitionEngineTest.ExposesObserverFinalizationAndSunFluxIntensityObjectHooks`、`SurfaceTransitionEngineTest.EvaluatesExtendedEventBundleWhenSunFluxIsDisabled`、`SurfaceTransitionEngineTest.ExtendedEventBundleDefaultsToDisabled`、`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsTopTransitionFamilyMetadata` 全部 `PASS`。
- 验证方式：
  - 语法校验：使用 `C:\mingw64\bin\c++.exe -fsyntax-only` 对 `SurfaceTransitionEngine.cpp`、`DensePlasmaSurfaceCharging.cpp`、`SurfaceCharging_object_layer_test.cpp`、`SurfaceCharging_smoke_test.cpp` 逐文件通过。
  - 目标执行：直接运行 `build_codex/bin/SurfaceCharging_object_layer_test.exe` 与 `build_codex/bin/SurfaceCharging_smoke_test.exe` 对上述 gtest 过滤子集验证通过。

#### Round 25：Wave F `Circ/Solver` 剩余语义复核与补差

| 字段 | 内容 |
| --- | --- |
| 目标 | 完成 `spis/Circ/**` 与 `spis/Solver/**` 的全量 family 复核，补齐近似实现与缺失语义 |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Circ/**`；`ref/SPIS/num-master/src/main/java/spis/Solver/**`；`Tools/Coupling/**`；`Tools/Solver/**`；相关 field/particle 模块 |
| 实施项 | 执行 `R25-N1`~`R25-N4`，补全 family 映射、组合语义与测试 |
| 输出 | `Circ/Solver` family 级对齐证据 |
| 依赖 | `Round 20`，建议在 `Round 21`~`Round 24` 后执行 |
| 验收 | `Circ` 与 `Solver` 不再依赖“主线够用”作为对齐依据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `ref/SPIS/num-master/src/main/java/spis/Circ/**`；`ref/SPIS/num-master/src/main/java/spis/Solver/**`；`Tools/Coupling/include/SurfaceCircuitCoupling.h`；`Tools/Coupling/src/SurfaceCircuitCoupling.cpp`；`Tools/Solver/include/SurfaceSolverFacade.h`；`Tools/Solver/src/SurfaceSolverFacade.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R25-N1` | 建立 `Circ` family 映射表 | `已完成` | `2026-04-16` | `Round 20` | `Circ/Circ`、`CircField`、`Circ/DIDV` 映射表 | 文档检视 | `ref/SPIS/num-master/src/main/java/spis/Circ/**`；本文档 2.10.2 | 进入 `R25-N2` |
| `R25-N2` | 补齐 `Circ` 剩余组合语义与测试 | `已完成` | `2026-04-16` | `R25-N1` | `Circ` family 代码与测试 | `Coupling_test` + smoke | `Tools/Coupling/include/SurfaceCircuitCoupling.h`；`Tools/Coupling/src/SurfaceCircuitCoupling.cpp`；`Tools/Coupling/test/Coupling_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 进入 `R25-N3` |
| `R25-N3` | 建立 `Solver` family 映射表 | `已完成` | `2026-04-16` | `R25-N2` | `Circuit`、`ElectroMag`、`Matter` 映射表 | 文档检视 | `ref/SPIS/num-master/src/main/java/spis/Solver/**`；本文档 2.10.2 | 进入 `R25-N4` |
| `R25-N4` | 补齐 `Solver` 剩余语义、路由与测试 | `已完成` | `2026-04-16` | `R25-N3` | solver family 代码与测试 | solver 单测 + circuit/surface 联动 smoke | `Tools/Solver/include/SurfaceSolverFacade.h`；`Tools/Solver/src/SurfaceSolverFacade.cpp`；`Tools/Solver/test/SurfaceSolverFacade_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 进入 `Round 26` |

##### `R25-N1`~`R25-N4` 实施与验证闭环（2026-04-16）
- 已在 `Tools/Coupling` 中补充 `CircuitFamilyRouteInput` / `CircuitFamilyView`，把 `Circ/Circ`、`CircField`、`Circ/DIDV` 的 supported / active family 签名固定下来。
- 已在 `Tools/Solver` 中补充 `SurfaceSolverFamilyRouteInput` / `SurfaceSolverFamilyView`，把 `Solver/Circuit`、`ElectroMag`、`Matter` 的 supported / active family 签名固定下来。
- `DensePlasmaSurfaceCharging` metadata sidecar 已新增 `surface_circ_*` 与 `surface_solver_*family_*` 系列字段，形成 `Circ/Solver` family 的 route 证据链。
- 目标测试：
  - `CouplingTest.ResolveCircuitFamilyViewBuildsSupportedAndActiveSignatures`
  - `SurfaceSolverFacadeTest.ResolveSurfaceSolverFamilyViewBuildsSupportedAndActiveSignatures`
  - `SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsCircFamilyMetadata`
  - `SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsSolverFamilyMetadata`
  - `SurfaceChargingSmokeTest.ReferenceDidvMatchesFiniteDifferenceAcrossGeoLeoRoutes`
  - `SurfaceChargingSmokeTest.SolverConfigImplicitCouplingDrivesSharedGlobalCoupledPolicy`
  均已 `PASS`。
- 验证方式：
  - 语法校验：使用 `C:\mingw64\bin\c++.exe -fsyntax-only` 对 `SurfaceCircuitCoupling.cpp`、`SurfaceSolverFacade.cpp`、`DensePlasmaSurfaceCharging.cpp`、相关测试文件逐文件通过。
  - 目标执行：直接运行 `build_codex/bin/Coupling_test.exe`、`build_codex/bin/SurfaceSolverFacade_test.exe`、`build_codex/bin/SurfaceCharging_smoke_test.exe` 的对应 gtest 子集验证通过。

#### Round 26：Wave G 1:1 对齐 gate 与封板

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 1:1 对齐变为可审计状态，建立 family coverage matrix、family-to-test matrix、family-to-config-route matrix 与最终封板证据 |
| 输入 | 前述所有轮次产出；`scripts/run/surface_reference_matrix.json`；contract/schema gates |
| 实施项 | 执行 `R26-O1`~`R26-O5`，建立 coverage gate、长尾参考矩阵与最终封板文档 |
| 输出 | 1:1 对齐证据链、family coverage gate、最终封板结论 |
| 依赖 | `Round 20`~`Round 25` |
| 验收 | 每个目标 family 都同时具备代码落点、测试落点、配置/路由证据与文档状态 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `documentation/contracts/surface_family_coverage_matrix_v1.json`；`documentation/contracts/kernel_contract_catalog_v1.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/python/check_surface_family_coverage_gate.py`；`scripts/python/check_surface_round26_sealoff.py`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate_round26_full/surface_reference_matrix_gate.json`；`build_codex/surface_round26_sealoff/surface_round26_sealoff.json`；本文档 |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R26-O1` | 建立 SPIS family coverage matrix | `已完成` | `2026-04-16` | `Round 25` | family-to-code 覆盖表 | 文档检视 + gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R26-O2` |
| `R26-O2` | 建立 family-to-test matrix | `已完成` | `2026-04-16` | `R26-O1` | family-to-test 覆盖表 | 文档检视 + gate | `scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/python/check_surface_family_coverage_gate.py`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R26-O3` |
| `R26-O3` | 建立 family-to-config-route matrix | `已完成` | `2026-04-16` | `R26-O2` | family-to-route 覆盖表 | 文档检视 + gate | `scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/python/check_surface_family_coverage_gate.py`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R26-O4` |
| `R26-O4` | 扩展 long-tail reference matrix 与 direct-SPIS semantic comparison cases | `已完成` | `2026-04-16` | `R26-O3` | 长尾 matrix case 与对比 artifact | matrix gate + contract gate | `scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_round26_o4_subset.json`；`build_codex/surface_reference_matrix_gate_round26_o4/surface_reference_matrix_gate.json` | 进入 `R26-O5` |
| `R26-O5` | 1:1 对齐封板与最终达成声明 | `已完成` | `2026-04-16` | `R26-O4` | 最终封板文档与证据链 | 全量 gate | `scripts/python/check_surface_round26_sealoff.py`；`build_codex/surface_reference_matrix_gate_round26_full/surface_reference_matrix_gate.json`；`build_codex/surface_round26_sealoff/surface_round26_sealoff.json`；本文档 2.12 | 项目转入长期回归维护，并保留 `SurfInteract` / `Vol` 主缺口 backlog |

##### `R26-O1`~`R26-O5` 实施与验证闭环（2026-04-16）
- 已新增 `surface-spis-family-coverage-matrix-v1` 契约文件，冻结 `target_families`、`family_to_code`、`family_to_test`、`family_to_config_route` 四段结构。
- 已新增 `scripts/run/surface_spis_family_coverage_matrix.json`，覆盖 `Surf / Top / Vol / Circ / Solver` 共 `12` 个 tracked family，并把 `Round 21`~`Round 25` 已沉淀的代码、测试、route/metadata 证据统一收口到 machine-readable artifact。
- 已新增 `scripts/python/check_surface_family_coverage_gate.py`，并生成 `build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` 与 `.md`；当前 gate 结果为 `PASS`，`12/12` family 的 code / test / config_route 三类证据均完整。
- 已扩展 `scripts/run/surface_reference_matrix.json` 中的 `leo_pic_circuit_ram_facing_hybrid_multi_patch` 长尾 case，把 `surface_circ_*` 与 `surface_solver_*family_signature` 系列 metadata 纳入 direct-SPIS semantic comparison 预期。
- 已使用 `build_codex/surface_reference_matrix_round26_o4_subset.json` 运行 targeted matrix gate；`build_codex/surface_reference_matrix_gate_round26_o4/surface_reference_matrix_gate.json` 当前为 `PASS`，新增 `Circ/Solver` metadata 预期已被 reference matrix gate 验证。
- 已重新执行更新后的全量 `surface_reference_matrix` gate；`build_codex/surface_reference_matrix_gate_round26_full/surface_reference_matrix_gate.json` 当前为 `PASS`，`25/25` parity-applicable case 全部通过。
- 已新增 `scripts/python/check_surface_round26_sealoff.py`，并生成 `build_codex/surface_round26_sealoff/surface_round26_sealoff.json` 与 `.md`；当前 sealoff gate 自身状态为 `PASS`，但最终 alignment verdict 为 `not_yet_achieved`。
- `Round 26` 的最终封板声明因此固定为：审计化封板完成、证据链完整、全量 gate 通过，但按 2.11 节标准当前项目仍未达到功能 1:1，对齐主缺口保留在 `SurfInteract` 与 `Vol` full stack。 

### 4.14 Round 27-33：Post-R26 explicit-implementation remediation

- 本节是 `Post-R26` 的唯一正文主线，用来处理 `Round 22`~`Round 26` 中所有 `near_equivalent`、`近似可达`、`未显式实现`、`已完成（family 可观测）` 但仍未满足 2.11 节 1:1 判定标准的项。
- 本节任务完成前，`surface_spis_family_coverage_matrix.json` 中除 `surf_distrib`、`surf_interact` 与已在后续轮次完成回填的 family 外，非 `completed` family 一律不得视为历史已收口项。

#### Round 27：Wave H `SurfInteract` 口径校正与显式实现轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 解决 `Round 21` 已完成但 sealoff 仍为 `surf_interact = not_fully_implemented` 的矛盾，补齐全部 `SurfInteract` 家族的显式实现与 direct evidence |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfInteract/**`；`Tools/Material/include/SurfaceInteraction.h`；`Tools/Material/src/SurfaceInteraction.cpp`；`scripts/run/surface_spis_family_coverage_matrix.json` |
| 实施项 | 执行 `R27-A1`~`R27-A4`，重新盘点已实现/缺失家族，补显式 family resolver、active signature、主线路由与测试 |
| 输出 | `SurfInteract` 三态盘点表、显式 family 路由、更新后的 coverage 状态 |
| 依赖 | `Round 26` |
| 验收 | `surf_interact` 不再只依赖“主线路由已接入”口径；只有显式实现 + direct evidence 完整时才允许改为 `completed` |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfInteract/**`；`Tools/Material/include/SurfaceInteraction.h`；`Tools/Material/src/SurfaceInteraction.cpp`；`Tools/Material/test/Material_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_reference_matrix_gate_round27/surface_reference_matrix_gate.json`；`build_codex/surface_round26_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R27-A1` | 重新盘点 `SurfInteract` 已实现家族与缺失家族 | `已完成` | `2026-04-16` | `Round 26` | “已显式实现 / 仍缺显式实现 / 仅元数据可达” 三态结论：当前 SPIS `SurfInteract` 家族均已有显式落点，缺口集中在 family view / direct route evidence 而非缺少枚举项 | 文档检视 + 代码检视 | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfInteract/**`；`Tools/Material/include/SurfaceInteraction.h`；`Tools/Material/src/SurfaceInteraction.cpp` | 进入 `R27-A2` |
| `R27-A2` | 逐项核对 `MaterialInteractor`、`MaterialDFInteractor`、`GenericDFInteractor`、`MaxwellianInteractor`、`MaxwellianInteractorWithRecollection`、`MultipleInteractor`、`MultipleMaxwellianInteractor`、`YieldInteractor`、`ReflectionInteractor`、`ErosionInteractor`、`ImprovedPhotoEmInteractor`、`TabulatedSEYModel`、`RecollectedSEYModel`、`DefaultPEEModel`、`DefaultSEEEYModel`、`DefaultSEEPModel`、`DefaultErosionModel`、`BasicInducedConductInteractor`、`Device` 的显式落点 | `已完成` | `2026-04-16` | `R27-A1` | 全家族显式落点总表 | 代码检视 | `Tools/Material/include/SurfaceInteraction.h`；`Tools/Material/src/SurfaceInteraction.cpp`；`Tools/Material/test/Material_test.cpp` | 进入 `R27-A3` |
| `R27-A3` | 为 `SurfInteract` 家族补显式 family resolver、active family signature、主线路由与 direct test evidence | `已完成` | `2026-04-16` | `R27-A2` | `surface_interactor_supported_family_signature` / `surface_interactor_active_family_signature` 与 direct metadata / matrix evidence | 单测 + smoke + matrix case | `Tools/Material/include/SurfaceInteraction.h`；`Tools/Material/src/SurfaceInteraction.cpp`；`Tools/Material/test/Material_test.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json` | 进入 `R27-A4` |
| `R27-A4` | 在 `surf_interact.alignment_status = completed` 后封板 `Round 27` | `已完成` | `2026-04-16` | `R27-A3` | 更新后的 coverage 状态与 sealoff 输入 | family coverage gate + sealoff 复核 | `scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_round26_sealoff/surface_round26_sealoff.json` | 进入 `Round 28` |

#### Round 28：Wave I `Vol BC/Field` 显式实现轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `R23-L2`、`R23-L3` 中标为 `近似可达` 或仅 `family 可观测` 的 `Vol BC / Field` 家族升级为显式原生实现 |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Vol/**`；`Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/FieldSolver/include/SurfaceFieldVolumeBridge.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；本文档 `R23-L1` 明细表 |
| 实施项 | 执行 `R28-B1`~`R28-B4`，补 config/route 入口、supported/active family signature、工具层单测、surface smoke 与 matrix case |
| 输出 | `Vol BC / Field` 显式 family 路由与 direct evidence |
| 依赖 | `Round 27` 可并行摸底，但封板依赖 `Round 27` 结果回填 |
| 验收 | 以下 family 不再只是 `近似可达` 或 `family 可观测`：`VoltageDependentMBC`、`MixedDirichletFourierPoissonBC`、`FourierPoissonBC`、`SurfDistribMatterBC`、`OneSurfDistribTestableMatterBC`、`CapacitiveVoltageGenerator`、`UniformBField`、`SolenoidBField`、`DipolarBField`、`MultipleEField`、`EMField` |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | 本文档 `R23-L1` 明细表；`Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Tools/FieldSolver/test/SurfaceFieldVolumeBridge_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R28-B1` | 为 `VoltageDependentMBC`、`MixedDirichletFourierPoissonBC`、`FourierPoissonBC`、`SurfDistribMatterBC`、`OneSurfDistribTestableMatterBC`、`CapacitiveVoltageGenerator`、`UniformBField`、`SolenoidBField`、`DipolarBField`、`MultipleEField`、`EMField` 建立独立 config/route 入口 | `已完成` | `2026-04-16` | `Round 26` | 显式 config/route 表 | 代码检视 + config 检视 + CLI smoke | `Main/main.cpp`；`Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`documentation/contracts/config_schema_v1.json`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`（`SurfaceChargingSmokeTest.SurfaceConfigCliAcceptsNativeVolumeFamilyOverrideKeys`） | 进入 `R28-B2` |
| `R28-B2` | 在 `SurfaceFieldVolumeContracts / SurfaceFieldVolumeBridge` 中补 supported/active family signature，不再只输出总 signature | `已完成` | `2026-04-16` | `R28-B1` | active family signature | 单测 + sidecar 检查 | `Tools/FieldSolver/include/SurfaceFieldVolumeBridge.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Tools/FieldSolver/test/SurfaceFieldVolumeBridge_test.cpp`（`SurfaceFieldVolumeBridgeTest.ExplicitBoundaryAndFieldFamiliesNormalizeDuplicates`）；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`（`SurfaceChargingSmokeTest.ScdatUnifiedRouteExportsNativeVolumeParityMetadata`） | 进入 `R28-B3` |
| `R28-B3` | 为每个 `Vol BC / Field` family 增加至少一个工具层单测和一个 surface smoke / matrix case | `已完成` | `2026-04-16` | `R28-B2` | direct evidence | 单测 + smoke + matrix gate | `Tools/FieldSolver/test/SurfaceFieldVolumeBridge_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `R28-B4` |
| `R28-B4` | 仅在 `Vol BC / Field` 家族全部具备显式 route + test + smoke + matrix 证据时，才允许进入 `vol_full_stack` 提升候选 | `已完成` | `2026-04-16` | `R28-B3` | `vol_full_stack` 子域提升候选记录 | family coverage gate + 文档检视 | `scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；本文档 2.10.4 | 进入 `Round 29` |

#### Round 29：Wave J `Vol Distrib` 显式实现轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `R23-L4`、`R23-L5` 中标为 `近似可达` 或 `未显式实现` 的 `Vol Distrib` 家族升级为显式原生实现 |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Vol/**`；`Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；本文档 `R23-L1` 明细表 |
| 实施项 | 执行 `R29-C1`~`R29-C3`，为分布家族补 route / metadata / 单测 / smoke / matrix case |
| 输出 | `Vol Distrib` 显式家族闭环 |
| 依赖 | `Round 28` |
| 验收 | 以下 family 不再只是 `近似可达 / 未显式实现`：`PICVolDistrib`、`PICVolDistribNoAcc`、`PICVolDistribUpdatable`、`SmartPICVolDistrib`、`CompositeVolDistrib`、`BackTrackingVolDistrib`、`BacktrackingPICCompositeVolDistrib`、`BacktrackingBoltzmannCompositeVolDistrib` |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | 本文档 `R23-L1` 明细表；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Tools/FieldSolver/test/SurfaceFieldVolumeBridge_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R29-C1` | 为 `PICVolDistrib`、`PICVolDistribNoAcc`、`PICVolDistribUpdatable`、`SmartPICVolDistrib`、`CompositeVolDistrib`、`BackTrackingVolDistrib`、`BacktrackingPICCompositeVolDistrib`、`BacktrackingBoltzmannCompositeVolDistrib` 建立独立 route 或显式 strategy 分支 | `已完成` | `2026-04-16` | `Round 28` | 显式 route / strategy 表 | 代码检视 | `Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | 进入 `R29-C2` |
| `R29-C2` | 为 no-acc / updatable / smart / backtracking 家族增加可观测 metadata，禁止继续只靠 bridge 内部拼装 | `已完成` | `2026-04-16` | `R29-C1` | 运行态 metadata 证据 | smoke + sidecar 检查 | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json` | 进入 `R29-C3` |
| `R29-C3` | 为每个 `Vol Distrib` family 增加最小单测 + smoke + reference matrix case | `已完成` | `2026-04-16` | `R29-C2` | direct evidence | 单测 + smoke + matrix gate | `Tools/FieldSolver/test/SurfaceFieldVolumeBridge_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `Round 30` |

#### Round 30：Wave K `Vol Interact` 显式实现轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `R23-L6` 中标为 `近似可达` 或 `未显式实现` 的 `Vol Interact` 家族升级为显式原生实现 |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Vol/**`；`Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；本文档 `R23-L1` 明细表 |
| 实施项 | 执行 `R30-D1`~`R30-D3`，补显式 family 入口、active signature、route 证据，并在四组齐全后重新评估 `vol_full_stack` |
| 输出 | `Vol Interact` 显式家族闭环与 `vol_full_stack` 重评估 |
| 依赖 | `Round 29` |
| 验收 | 以下 family 不再只是 `近似可达 / 未显式实现`：`MCCInteractor`、`CEXInteractor`、`PhotoIonization`、`TrajectoryInteractionFromField`、`SpinningSpacecraftTrajectory`；若 `BC / Field / Distrib / Interact` 四组都闭环，才允许把 `vol_full_stack.alignment_status` 改为 `completed` |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | 本文档 `R23-L1` 明细表；`Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/src/VolInteract.cpp`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_family_coverage_gate_round30/surface_family_coverage_gate.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R30-D1` | 将 `MCCInteractor`、`CEXInteractor`、`PhotoIonization`、`TrajectoryInteractionFromField`、`SpinningSpacecraftTrajectory` 从“框架可承载”升级为显式 family 入口 + active signature + route 证据 | `已完成` | `2026-04-16` | `Round 29` | 显式 family 入口 | 代码检视 + 单测 | `Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/src/VolInteract.cpp`；`Tools/FieldSolver/include/SurfaceFieldVolumeBridge.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Tools/FieldSolver/test/SurfaceFieldVolumeBridge_test.cpp` | 进入 `R30-D2` |
| `R30-D2` | 让 `surface_native_volume_active_interaction_family_signature` 反映真实运行态，而不是仅 parity 可观测 | `已完成` | `2026-04-16` | `R30-D1` | 运行态 active signature | smoke + matrix gate | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Main/SurfaceScenarioLoader.cpp`；`documentation/contracts/config_schema_v1.json`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json` | 进入 `R30-D3` |
| `R30-D3` | 对 `vol_full_stack` 执行整体重评估；若 `BC / Field / Distrib / Interact` 四组全部闭环，才允许把 `vol_full_stack.alignment_status` 改为 `completed` | `已完成` | `2026-04-16` | `R30-D2` | 更新后的 `vol_full_stack` 状态 | family coverage gate + sealoff 预检查 | `scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate_round30/surface_family_coverage_gate.json` | 进入 `Round 31`；本轮仅完成重评估，`vol_full_stack.alignment_status` 维持 `not_fully_implemented` |

#### Round 31：Wave L `SurfField + Top + Circ` near-equivalent 提升轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `Round 22`~`Round 25` 中仍为 `near_equivalent` 的 `surf_field`、`top_transition`、`circ_circuit`、`circ_circfield`、`circ_didv` 全部提升到 `completed` |
| 输入 | `scripts/run/surface_spis_family_coverage_matrix.json`；`Tools/FieldSolver/**`；`Tools/Coupling/**`；`Toolkit/Surface Charging/**` |
| 实施项 | 执行 `R31-E1`~`R31-E5`，按 family 组逐项补 direct route / config / smoke / matrix evidence |
| 输出 | 更新后的 `surf_field`、`top_transition`、`circ_*` 状态 |
| 依赖 | `Round 30` |
| 验收 | 不再允许把“只存在 family view / object-layer signature / metadata signature”当作 `completed` 依据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Tools/Coupling/test/Coupling_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json`；`scripts/run/surface_spis_family_coverage_matrix.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R31-E1` | `surf_field`：为 `OMLCurrentScaler`、`LTE_OML_CurrentScaler`、`FowlerNordheimCurrentScaler`、`VariableBarrierCurrentScaler`、`AutomaticBarrierCurrentScaler`、`MultipleCurrentScaler`、`CurrentScalerFromCurrentVariation`、`CurrentScalerFromLocalCurrentVariation` 补 direct route / config / smoke evidence | `已完成` | `2026-04-16` | `Round 30` | `surf_field` direct evidence | 单测 + smoke + matrix case | `Tools/FieldSolver/test/SurfaceBarrierModels_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `R31-E2` |
| `R31-E2` | `top_transition`：让 `LocalTimeTransition`、`SpinningSpacecraft`、`SunFluxUpdater`、`ConductivityEvolution`、`SourceFluxUpdater`、`SimulationParamUpdater`、`SheathOrPresheathPoissonBCUpdater`、`RCCabsSCUpdater`、`VcrossBfieldUpdater`、`BasicEclipseExit`、`TransientArtificialSources`、`LangmuirProbeTransition`、`TransitionObserver`、`Finalization`、`SunFluxIntensityUpdater` 全部具备 route evidence，而非仅 object-layer signature | `已完成` | `2026-04-16` | `R31-E1` | `top_transition` direct evidence | object-layer test + smoke + matrix case | `Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `R31-E3` |
| `R31-E3` | `circ_circuit`：除 current metadata 外，补 `Circuit`、`ElecComponent`、`SurfaceComponent`、`RCCabsCirc`、`RLCCirc`、`SIN`、`PULSE`、`PWL`、`EXP` 在 surface 主线中的 direct route / smoke 证据 | `已完成` | `2026-04-16` | `R31-E2` | `circ_circuit` direct evidence | `Coupling_test` + smoke + matrix case | `Tools/Coupling/test/Coupling_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `R31-E4` |
| `R31-E4` | `circ_circfield`：补 `CircField`、`CurrentDistribution`、`DirCircField` 的 active 路由和运行态断言 | `已完成` | `2026-04-16` | `R31-E3` | `circ_circfield` direct evidence | 单测 + smoke | `Tools/Coupling/test/Coupling_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` | 进入 `R31-E5` |
| `R31-E5` | `circ_didv`：对 `SurfDIDVFromMatrices`、`WeightedSurfDIDV`、`MultipleSurfDIDV`、`ReducableDIDV`、`RedDIDVFromSurfDIDV`、`RedDIDVFromRegDIDV`、`DIDVOnCirc` 等家族补 direct route evidence，避免只依赖 helper / metadata | `已完成` | `2026-04-16` | `R31-E4` | `circ_didv` direct evidence | 单测 + smoke + matrix case | `Tools/Coupling/test/Coupling_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `Round 32` |

#### Round 32：Wave M `SurfMaterial + Solver` near-equivalent 提升轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 sealoff 中仍为 `near_equivalent` 的 `surf_material`、`solver_circuit`、`solver_electromag`、`solver_matter` 提升到 `completed` |
| 输入 | `scripts/run/surface_spis_family_coverage_matrix.json`；`Tools/Material/**`；`Tools/Solver/**`；`Toolkit/Surface Charging/**` |
| 实施项 | 执行 `R32-F1`~`R32-F4`，按 family 补齐 route 条件、direct test evidence 与 matrix case |
| 输出 | 更新后的 `surf_material` 与 `solver_*` 状态 |
| 依赖 | `Round 31` |
| 验收 | `surf_material` 与 `solver_*` 不再只依赖 family view / metadata 推断，必须具备 direct-SPIS semantic comparison 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `Tools/Material/test/Material_test.cpp`；`Tools/Solver/test/SurfaceSolverFacade_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json`；`scripts/run/surface_spis_family_coverage_matrix.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R32-F1` | `surf_material`：为 `BasicMaterialModel`、`ErosionMaterialModel` 补参数资产、默认映射、family-to-config 与 family-to-test 的最后闭环 | `已完成` | `2026-04-16` | `Round 31` | `surf_material` completed evidence | `Material_test` + smoke + matrix case | `Tools/Material/test/Material_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `R32-F2` |
| `R32-F2` | `solver_circuit`：将 `CircSolve` 的运行态证据从“存在 family view”提升到“direct-SPIS semantic comparison 可单独验证” | `已完成` | `2026-04-16` | `R32-F1` | `solver_circuit` direct evidence | solver 单测 + smoke + matrix case | `Tools/Solver/test/SurfaceSolverFacade_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `R32-F3` |
| `R32-F3` | `solver_electromag`：对 `AbstractEMSolver`、`PoissonSolver`、`PotPoissonSolver`、`ConjGrad3DUnstructPoissonSolver`、`PoissonMagnetostaticSolver` 分别补 route 条件与 smoke / matrix evidence | `已完成` | `2026-04-16` | `R32-F2` | `solver_electromag` direct evidence | solver 单测 + smoke + matrix case | `Tools/Solver/test/SurfaceSolverFacade_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `R32-F4` |
| `R32-F4` | `solver_matter`：对 `ParticlePusher`、`PICPusher`、`CrossTetra*` 家族补 route 触发条件与 direct evidence，不再只由 `SurfaceSolverFamilyView` 推断 | `已完成` | `2026-04-16` | `R32-F3` | `solver_matter` direct evidence | solver 单测 + smoke + matrix case | `Tools/Solver/test/SurfaceSolverFacade_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `Round 33` |

#### Round 33：Wave N 最终重封板轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 重新执行 `R26-O1`~`R26-O5` 的全部 gate，并仅在全部 family 为 `completed` 时生成真正的 `1:1 achieved` verdict |
| 输入 | 更新后的 `scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/run/surface_reference_matrix.json`；`scripts/python/check_surface_family_coverage_gate.py`；`scripts/python/check_surface_round26_sealoff.py` |
| 实施项 | 执行 `R33-G1`~`R33-G4`，重跑 family coverage gate、全量 reference matrix gate、sealoff gate，并回填文档结论 |
| 输出 | 最终 1:1 verdict 或明确的 `未达成 / 受阻` 声明 |
| 依赖 | `Round 32` |
| 验收 | 只有当所有 target family 的 `alignment_status = completed` 且 `final_alignment_verdict = achieved` 时，才允许宣告“1:1 已完成” |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json`；`build_codex/surface_round33_sealoff/surface_round26_sealoff.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/run/surface_reference_matrix.json`；本文档 2.12 |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R33-G1` | 重跑 `surface_family_coverage_gate`，要求所有 `target_families[].alignment_status = completed` | `已完成` | `2026-04-16` | `Round 32` | 更新后的 coverage gate report | gate | `build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R33-G2` |
| `R33-G2` | 重跑全量 `surface_reference_matrix gate` | `已完成` | `2026-04-16` | `R33-G1` | 更新后的 full matrix gate report | gate | `build_codex/surface_reference_matrix_gate_round33_full/surface_reference_matrix_gate.json` | 进入 `R33-G3` |
| `R33-G3` | 重跑 sealoff gate，要求 `final_alignment_verdict = achieved` | `已完成` | `2026-04-16` | `R33-G2` | 最终 sealoff report | gate | `build_codex/surface_round33_sealoff/surface_round26_sealoff.json` | 进入 `R33-G4`；`final_alignment_verdict = achieved` |
| `R33-G4` | 回填 2.12 / 2.13 / 4.2 / 4.14 文档结论；若仍有 family 不是 `completed`，则必须记为 `未达成` 或 `受阻` | `已完成` | `2026-04-16` | `R33-G3` | 最终文档结论 | 文档检视 + gate 复核 | 本文档；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_round33_sealoff/surface_round26_sealoff.json` | 项目进入长期回归维护 |

#### Round 34：Wave O `Util/**` 扩展 target family 轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `spis/Util/**` 从“Surface 主线刚需支撑域”升级为扩展 target family，并补 `util_full_stack` contract / smoke / matrix / coverage gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Util/**`；`Tools/Basic/**`；`Tools/Geometry/**`；`Tools/Diagnostics/**`；`Tools/Output/**`；`Tools/Particle/**`；`Tools/Solver/**` |
| 实施项 | 执行 `R34-H1`~`R34-H3`，扩展 `surface_family_coverage_matrix_v1`、`surface_spis_family_coverage_matrix.json`、runtime metadata、smoke 与 reference matrix |
| 输出 | `util_full_stack` target family、`surface_util_*` metadata、扩展后的 gate 报告 |
| 依赖 | `Round 33` |
| 验收 | `Util/**` 不再只作为“部分纳入”描述；必须进入 target families，并通过 code / test / config-route / reference-matrix 证据链 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `documentation/contracts/surface_family_coverage_matrix_v1.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_reference_matrix_gate_round34_full/surface_reference_matrix_gate.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R34-H1` | 扩展 `surface_family_coverage_matrix_v1` 与 `surface_spis_family_coverage_matrix.json`，新增 `util_full_stack` target family | `已完成` | `2026-04-16` | `Round 33` | 扩展后的 contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R34-H2` |
| `R34-H2` | 在 runtime metadata 中新增 `surface_util_*` 家族签名，并补 smoke / reference matrix 期望值 | `已完成` | `2026-04-16` | `R34-H1` | `surface_util_route` / `surface_util_*family_*` 字段 | smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate_round34_full/surface_reference_matrix_gate.json` | 进入 `R34-H3` |
| `R34-H3` | 重新跑 `SurfaceFamilyCoverageGate` / `SurfaceAlignmentSealoffGate`，要求扩展后 `util_full_stack.alignment_status = completed` | `已完成` | `2026-04-16` | `R34-H2` | 扩展后的 gate / sealoff 报告 | gate | `build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_round26_sealoff/surface_round26_sealoff.json` | 进入 `Round 35` |

#### Round 35：Wave P `Top/**` 深化 family 盘点轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `Top/**` 从当前以 `Transition` 为主的 gate 扩展到 `Simulation`、`Plasma`、`SC`、`Grid`、`Default` 逐族盘点 |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Top/**`；`Toolkit/Surface Charging/**`；`Main/SurfaceScenarioLoader.cpp`；`Tools/Coupling/**` |
| 实施项 | 执行 `R35-I1`~`R35-I3`，扩展 target family、runtime metadata、object-layer / smoke / reference-matrix 证据 |
| 输出 | `top_simulation`、`top_plasma`、`top_sc`、`top_grid`、`top_default` 目标 family 与 gate 证据 |
| 依赖 | `Round 34` |
| 验收 | `Top/**` 不再只以 `top_transition` 为唯一 target family；五组新 family 必须具备独立 metadata / config-route / test / matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_reference_matrix_gate_round34_full/surface_reference_matrix_gate.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R35-I1` | 在 coverage matrix 中新增 `top_simulation`、`top_plasma`、`top_sc`、`top_grid`、`top_default` 五个 target family | `已完成` | `2026-04-16` | `Round 34` | 扩展后的 Top family contract | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R35-I2` |
| `R35-I2` | 在 `DensePlasmaSurfaceCharging` 中导出 `surface_top_{simulation,plasma,sc,grid,default}_*` metadata，并补 smoke / matrix 断言 | `已完成` | `2026-04-16` | `R35-I1` | Top 五组 family metadata | object-layer + smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate_round34_full/surface_reference_matrix_gate.json` | 进入 `R35-I3` |
| `R35-I3` | 重新跑 sealoff，要求扩展后 Top family 全部 `completed` | `已完成` | `2026-04-16` | `R35-I2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round26_sealoff/surface_round26_sealoff.json` | 进入 `Round 36` |

#### Round 36：Wave Q `Vol/Geom + Vol/VolMesh` 独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `Vol/Geom`、`Vol/VolMesh` 从“能力闭环”升级为独立 family 级 contract / gate，而不是继续隐含在 `vol_full_stack` 口径里 |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Vol/Geom/**`；`ref/SPIS/num-master/src/main/java/spis/Vol/VolMesh/**`；`Tools/Mesh/**`；`Tools/FieldSolver/**`；`Toolkit/Surface Charging/**` |
| 实施项 | 执行 `R36-J1`~`R36-J3`，扩展 target family、runtime metadata、native-volume smoke 与 reference-matrix 证据 |
| 输出 | `vol_geom`、`vol_mesh` family contract / route metadata / gate 结果 |
| 依赖 | `Round 35` |
| 验收 | `Vol/Geom`、`Vol/VolMesh` 必须拥有独立 target family、独立 metadata key、独立 config-route 映射与 reference-matrix case，不再只靠 `vol_full_stack` 代管 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `Tools/Mesh/include/MeshParsing.h`；`Tools/Mesh/src/MeshParsing.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_reference_matrix_gate_round34_full/surface_reference_matrix_gate.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R36-J1` | 在 coverage matrix 中新增 `vol_geom`、`vol_mesh` 两个 target family，并补 family-to-code / test / config-route 映射 | `已完成` | `2026-04-16` | `Round 35` | 扩展后的 Vol family contract | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R36-J2` |
| `R36-J2` | 在 runtime metadata 中新增 `surface_native_volume_parity_geom_*`、`surface_native_volume_parity_mesh_*` 字段，并补 native-volume smoke / matrix 断言 | `已完成` | `2026-04-16` | `R36-J1` | Geom / VolMesh family metadata | smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate_round34_full/surface_reference_matrix_gate.json` | 进入 `R36-J3` |
| `R36-J3` | 重新跑 sealoff，要求扩展后 `vol_geom`、`vol_mesh` 均为 `completed` | `已完成` | `2026-04-16` | `R36-J2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round26_sealoff/surface_round26_sealoff.json` | 项目回到长期回归维护 |

#### Round 37：Wave R `Solver/Util` 独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `spis/Solver/Util/**` 从 `Solver/**` 的隐含支撑层升级为独立 target family，并补 `SolverUtil` / `LeastSquare` 的 runtime metadata / tests / gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Solver/Util/**`；`Tools/Solver/**`；`Tools/FieldSolver/**`；`Tools/PICcore/**`；`Toolkit/Surface Charging/**` |
| 实施项 | 执行 `R37-K1`~`R37-K3`，扩展 target family、runtime metadata、solver smoke / reference-matrix 证据与 sealoff |
| 输出 | `solver_util` target family、`surface_solver_util_*` metadata、扩展后的 gate / sealoff 报告 |
| 依赖 | `Round 36` |
| 验收 | `Solver/Util` 不再隐含在 `solver_electromag` 或工具层实现中；必须具备独立 target family、独立 metadata key、独立 test / config-route / reference-matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/run/surface_reference_matrix.json`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Tools/FieldSolver/test/PoissonSolver_test.cpp`；`Tools/PICcore/test/PICCycleBoundary_test.cpp`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_reference_matrix_gate_round37_full/surface_reference_matrix_gate.json`；`build_codex/surface_round26_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R37-K1` | 在 coverage matrix / contract 中新增 `solver_util` target family，并补 family-to-code / test / config-route 映射 | `已完成` | `2026-04-16` | `Round 36` | `solver_util` contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R37-K2` |
| `R37-K2` | 在 runtime metadata 中新增 `surface_solver_util_route`、`surface_solver_util_*family_*` 字段，并补 smoke / reference matrix 断言 | `已完成` | `2026-04-16` | `R37-K1` | solver util family metadata | smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate_round37_full/surface_reference_matrix_gate.json` | 进入 `R37-K3` |
| `R37-K3` | 重新跑 `SurfaceFamilyCoverageGate` / `SurfaceAlignmentSealoffGate`，要求扩展后 `solver_util.alignment_status = completed` | `已完成` | `2026-04-16` | `R37-K2` | 扩展后的 gate / sealoff 报告 | gate | `build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_round26_sealoff/surface_round26_sealoff.json` | 下一轮优先评估 `Top/Top` 或 `Surf/SurfMesh` |

#### Round 38：Wave S `Top/Top` 独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `spis/Top/Top/**` 从 Top 剩余入口层升级为独立 target family，并补 `Scenario`、`PotentialSweep`、`UIInvokable`、`SpisTopMenu`、`NumTopFromUI` 的 runtime metadata / tests / gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Top/Top/**`；`Toolkit/Surface Charging/**`；`Main/main.cpp`；`Main/SurfaceScenarioLoader.cpp` |
| 实施项 | 执行 `R38-L1`~`R38-L3`，扩展 target family、runtime-plan entrypoint metadata、smoke / reference-matrix 证据与 sealoff |
| 输出 | `top_top` target family、`surface_top_top_*` metadata、扩展后的 gate / sealoff 报告 |
| 依赖 | `Round 37` |
| 验收 | `Top/Top` 不再留在 backlog；必须具备独立 target family、独立 metadata key、独立 test / config-route / reference-matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `Toolkit/Surface Charging/include/SurfaceChargingCases.h`；`Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Main/SurfaceScenarioLoader.cpp`；`Main/main.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_reference_matrix_gate_round38_full/surface_reference_matrix_gate.json`；`build_codex/surface_round38_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R38-L1` | 在 coverage matrix / contract 中新增 `top_top` target family，并补 family-to-code / test / config-route 映射 | `已完成` | `2026-04-16` | `Round 37` | `top_top` contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R38-L2` |
| `R38-L2` | 在 runtime metadata 中新增 `surface_top_top_route`、`surface_top_top_entrypoint_mode`、`surface_top_top_*family_*` 字段，并补 smoke / reference matrix 断言 | `已完成` | `2026-04-16` | `R38-L1` | top-top family metadata | smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate_round38_full/surface_reference_matrix_gate.json` | 进入 `R38-L3` |
| `R38-L3` | 重新跑 sealoff，要求扩展后 `top_top.alignment_status = completed` | `已完成` | `2026-04-16` | `R38-L2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round38_sealoff/surface_round26_sealoff.json` | 下一轮优先进入 `Surf/SurfMesh` |

#### Round 39：Wave T `Surf/SurfMesh` 独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `spis/Surf/SurfMesh/**` 从隐含的 structured topology / boundary mapping 能力升级为独立 target family，并补 `SurfMesh`、`ThreeDUnstructSurfMesh` 的 runtime metadata / tests / gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Surf/SurfMesh/**`；`Toolkit/Surface Charging/**`；`Main/SurfaceScenarioLoader.cpp`；`Main/main.cpp`；`Tools/Mesh/**` |
| 实施项 | 执行 `R39-M1`~`R39-M3`，扩展 target family、surf-mesh metadata、structured-topology / boundary-mapping smoke 与 reference-matrix 证据 |
| 输出 | `surf_mesh` target family、`surface_surf_mesh_*` metadata、扩展后的 gate / sealoff 报告 |
| 依赖 | `Round 38` |
| 验收 | `Surf/SurfMesh` 不再停留在 backlog；必须具备独立 target family、独立 metadata key、独立 test / config-route / reference-matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-16` |
| 证据 | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Main/SurfaceScenarioLoader.cpp`；`Main/main.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_reference_matrix_gate_round39_full/surface_reference_matrix_gate.json`；`build_codex/surface_round39_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R39-M1` | 在 coverage matrix / contract 中新增 `surf_mesh` target family，并补 family-to-code / test / config-route 映射 | `已完成` | `2026-04-16` | `Round 38` | `surf_mesh` contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R39-M2` |
| `R39-M2` | 在 runtime metadata 中新增 `surface_surf_mesh_route`、`surface_surf_mesh_*family_*` 字段，并补 structured-topology / boundary-mapping smoke / reference matrix 断言 | `已完成` | `2026-04-16` | `R39-M1` | surf-mesh family metadata | smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate_round39_full/surface_reference_matrix_gate.json` | 进入 `R39-M3` |
| `R39-M3` | 重新跑 sealoff，要求扩展后 `surf_mesh.alignment_status = completed` | `已完成` | `2026-04-16` | `R39-M2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round39_sealoff/surface_round26_sealoff.json` | 下一轮继续 `Vol/BC` 抽象长尾 |

#### Round 40：Wave U `Vol/BC(Abstract)` 独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `spis/Vol/BC` 抽象长尾从 `vol_full_stack` 的隐含层升级为独立 target family，并补 `BC`、`DirichletPoissonBC`、`MatterBC`、`TestableMatterBC`、`VoltageGenerator` 的 runtime metadata / tests / gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Vol/BC/BC.java`、`DirichletPoissonBC.java`、`MatterBC.java`、`TestableMatterBC.java`、`VoltageGenerator.java`；`Tools/FieldSolver/**`；`Tools/Solver/**`；`Toolkit/Surface Charging/**` |
| 实施项 | 执行 `R40-N1`~`R40-N3`，扩展 target family、BC-abstract metadata、native-volume smoke / reference-matrix 证据与 sealoff |
| 输出 | `vol_bc_abstract` target family、`surface_native_volume_bc_abstract_*` metadata、扩展后的 gate / sealoff 报告 |
| 依赖 | `Round 39` |
| 验收 | `Vol/BC` 抽象长尾不再停留在 backlog；必须具备独立 target family、独立 metadata key、独立 test / config-route / reference-matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Tools/FieldSolver/include/BoundaryCondition.h`；`Tools/FieldSolver/include/PoissonSolver.h`；`Tools/Solver/include/LinearSystem.h`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_reference_matrix_gate_round40_full/surface_reference_matrix_gate.json`；`build_codex/surface_round40_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R40-N1` | 在 coverage matrix / contract 中新增 `vol_bc_abstract` target family，并补 family-to-code / test / config-route 映射 | `已完成` | `2026-04-17` | `Round 39` | `vol_bc_abstract` contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R40-N2` |
| `R40-N2` | 在 runtime metadata 中新增 `surface_native_volume_bc_abstract_route`、`surface_native_volume_bc_abstract_*family_*` 字段，并补 native-volume smoke / reference matrix 断言 | `已完成` | `2026-04-17` | `R40-N1` | BC abstract family metadata | smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate_round40_full/surface_reference_matrix_gate.json` | 进入 `R40-N3` |
| `R40-N3` | 重新跑 sealoff，要求扩展后 `vol_bc_abstract.alignment_status = completed` | `已完成` | `2026-04-17` | `R40-N2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round40_sealoff/surface_round26_sealoff.json` | 进入 `Round 41` |

#### Round 41：Wave V `Vol/VolField(Abstract)` 独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `spis/Vol/VolField` 抽象长尾从 `vol_full_stack` 的隐含层升级为独立 target family，并补 `VolField`、`ScalVolField`、`VectVolField`、`EField`、`DirVectVolField`、`DirEField`、`PotEField`、`PotVectVolField`、`BField` 的 runtime metadata / tests / gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Vol/VolField/**`；`Tools/FieldSolver/include/VolField.h`；`Tools/FieldSolver/src/VolField.cpp`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |
| 实施项 | 执行 `R41-O1`~`R41-O3`，扩展 target family、VolField-abstract metadata、native-volume smoke / reference-matrix 证据与 sealoff |
| 输出 | `vol_field_abstract` target family、`surface_native_volume_field_*` metadata、扩展后的 gate / sealoff 报告 |
| 依赖 | `Round 40` |
| 验收 | `Vol/VolField` 不再停留在 backlog；必须具备独立 target family、独立 metadata key、独立 test / config-route / reference-matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Tools/FieldSolver/include/VolField.h`；`Tools/FieldSolver/src/VolField.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_spis_family_coverage_matrix.json`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json`；`build_codex/surface_reference_matrix_gate_round41_full/surface_reference_matrix_gate.json`；`build_codex/surface_round41_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R41-O1` | 在 coverage matrix / contract 中新增 `vol_field_abstract` target family，并补 family-to-code / test / config-route 映射 | `已完成` | `2026-04-17` | `Round 40` | `vol_field_abstract` contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`；`scripts/run/surface_spis_family_coverage_matrix.json`；`build_codex/surface_family_coverage_gate/surface_family_coverage_gate.json` | 进入 `R41-O2` |
| `R41-O2` | 在 runtime metadata 中新增 `surface_native_volume_field_route`、`surface_native_volume_field_*family_*` 字段，并补 native-volume smoke / reference matrix 断言 | `已完成` | `2026-04-17` | `R41-O1` | VolField abstract family metadata | smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`build_codex/surface_reference_matrix_gate_round41_full/surface_reference_matrix_gate.json` | 进入 `R41-O3` |
| `R41-O3` | 重新跑 sealoff，要求扩展后 `vol_field_abstract.alignment_status = completed` | `已完成` | `2026-04-17` | `R41-O2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round41_sealoff/surface_round26_sealoff.json` | 进入 `Round 42` |

#### Round 42：Wave W `Vol/VolDistrib(LongTail)` 独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `spis/Vol/VolDistrib/**` 中的长尾 family 从 `vol_full_stack` 的隐含层升级为独立 target family，并补齐 `VolDistrib`、`NonPICVolDistrib`、`AnalyticVolDistrib`、`MultipleVolDistrib`、`LocalMaxwellVolDistrib`、`LocalMaxwellBoltzmannVolDistrib`、`GlobalMaxwellBoltzmannVolDistrib`、`PICBoltzmannVolDistrib`、`SteadyMaxwellBoltzmannVolDistrib`、`UnlimitedGlobalMaxwellBoltzmannVolDistrib`、`SurfaceLimitedGlobalMaxwellBoltzmannVolDistrib`、`TrunckatedGlobalMaxwellBoltzmannVolDistrib`、`ImplicitableVolDistrib`、`Updatable` 的 runtime metadata / tests / gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Vol/VolDistrib/**`、`Tools/FieldSolver/include/VolDistrib.h`、`Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`、`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`、`Main/SurfaceScenarioLoader.cpp`、`Main/main.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |
| 实施项 | 执行 `R42-P1`~`R42-P3`，扩展 target family、VolDistrib-long-tail metadata、config override / smoke / reference-matrix 证据与 sealoff |
| 输出 | `vol_distrib_long_tail` target family、`surface_native_volume_distribution_*` metadata、扩展后的 gate / sealoff 报告 |
| 依赖 | `Round 41` |
| 验收 | `Vol/VolDistrib` 长尾不再停留在候选 backlog；必须具备独立 target family、独立 metadata key、独立 test / config-route / reference-matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Tools/FieldSolver/include/VolDistrib.h`、`Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h`、`Tools/FieldSolver/include/SurfaceFieldVolumeBridge.h`、`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`、`Main/SurfaceScenarioLoader.cpp`、`Main/main.cpp`、`Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`Tools/FieldSolver/test/SurfaceFieldVolumeBridge_test.cpp`、`scripts/run/surface_spis_family_coverage_matrix.json`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_family_coverage_gate_round42/surface_family_coverage_gate.json`、`build_codex/surface_reference_matrix_gate_round42_full/surface_reference_matrix_gate.json`、`build_codex/surface_round42_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R42-P1` | 在 coverage matrix / contract 中新增 `vol_distrib_long_tail` target family，并补齐 family-to-code / test / config-route 映射 | `已完成` | `2026-04-17` | `Round 41` | `vol_distrib_long_tail` contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`、`scripts/run/surface_spis_family_coverage_matrix.json`、`build_codex/surface_family_coverage_gate_round42/surface_family_coverage_gate.json` | 进入 `R42-P2` |
| `R42-P2` | 在 runtime metadata 中新增 `surface_native_volume_distribution_route`、`surface_native_volume_distribution_*family_*` 字段，并补 config override / smoke / reference matrix 断言 | `已完成` | `2026-04-17` | `R42-P1` | VolDistrib long-tail family metadata | smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_reference_matrix_gate_round42_full/surface_reference_matrix_gate.json` | 进入 `R42-P3` |
| `R42-P3` | 重新跑 sealoff，要求扩展后 `vol_distrib_long_tail.alignment_status = completed` | `已完成` | `2026-04-17` | `R42-P2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round42_sealoff/surface_round26_sealoff.json` | 进入 `Round 43` |

#### Round 43：Wave X `Vol/VolInteract(LongTail)` 独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `spis/Vol/VolInteract/**` 中的长尾抽象/兼容 family 从既有 `Vol Interact` 显式实现层中升级为独立 target family，并补齐 `VolInteractor`、`VolInteractModel`、`VolInteractionAlongTrajectory`、`ElasticCollisions`、`ConstantIonizationInteractor` 的 runtime metadata / tests / gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Vol/VolInteract/VolInteractor.java`、`VolInteractModel.java`、`VolInteractionAlongTrajectory.java`、`ElasticCollisions.java`、`ConstantIonizationInteractor.java`；`Tools/FieldSolver/include/VolInteract.h`；`Tools/FieldSolver/src/VolInteract.cpp`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |
| 实施项 | 执行 `R43-Q1`~`R43-Q3`，扩展 target family、VolInteract-long-tail metadata、FieldSolver 单测 / smoke / reference-matrix 证据与 sealoff |
| 输出 | `vol_interact_long_tail` target family、`surface_native_volume_interaction_long_tail_*` metadata、扩展后的 gate / sealoff 报告 |
| 依赖 | `Round 42` |
| 验收 | `Vol/VolInteract` 长尾不再停留在候选 backlog；必须具备独立 target family、独立 metadata key、独立 test / config-route / reference-matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Tools/FieldSolver/include/VolInteract.h`、`Tools/FieldSolver/src/VolInteract.cpp`、`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`、`Tools/FieldSolver/test/VolInteract_test.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`documentation/contracts/surface_family_coverage_matrix_v1.json`、`scripts/run/surface_spis_family_coverage_matrix.json`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_family_coverage_gate_round43/surface_family_coverage_gate.json`、`build_codex/surface_reference_matrix_gate_round43_full/surface_reference_matrix_gate.json`、`build_codex/surface_round43_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R43-Q1` | 在 coverage matrix / contract 中新增 `vol_interact_long_tail` target family，并补齐 family-to-code / test / config-route 映射 | `已完成` | `2026-04-17` | `Round 42` | `vol_interact_long_tail` contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`、`scripts/run/surface_spis_family_coverage_matrix.json`、`build_codex/surface_family_coverage_gate_round43/surface_family_coverage_gate.json` | 进入 `R43-Q2` |
| `R43-Q2` | 在 runtime metadata 中新增 `surface_native_volume_interaction_long_tail_route`、`surface_native_volume_interaction_long_tail_*family_*` 字段，并补 FieldSolver 单测 / smoke / reference matrix 断言 | `已完成` | `2026-04-17` | `R43-Q1` | VolInteract long-tail family metadata | 单测 + smoke + matrix gate | `Tools/FieldSolver/test/VolInteract_test.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_reference_matrix_gate_round43_full/surface_reference_matrix_gate.json` | 进入 `R43-Q3` |
| `R43-Q3` | 重新跑 sealoff，要求扩展后 `vol_interact_long_tail.alignment_status = completed` | `已完成` | `2026-04-17` | `R43-Q2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round43_sealoff/surface_round26_sealoff.json` | 下一轮优先评估 `Util/**` 子包级拆分 |

#### Round 44：Wave Y `Util/**` 第一批子包独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `spis/Util/**` 中最直接影响 surface 主线的子包从 `util_full_stack` 聚合层升级为独立 target family，并补齐 `Util/Func`、`Util/Sampler`、`Util/Monitor`、`Util/io`、`Util/Matrix` 的 runtime metadata / tests / gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Util/Func/**`、`Util/Sampler/**`、`Util/Monitor/**`、`Util/io/**`、`Util/Matrix/**`；`Tools/Basic/**`；`Tools/Diagnostics/**`；`Tools/Output/**`；`Tools/Particle/**`；`Tools/Solver/**`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |
| 实施项 | 执行 `R44-Y1`~`R44-Y3`，扩展 target family、util 子包 metadata、Basic/Particle/Diagnostics/Output/Solver 单测与 reference-matrix / sealoff 证据 |
| 输出 | `util_func`、`util_sampler`、`util_monitor`、`util_io`、`util_matrix` target families，以及 `surface_util_*` 子包级 metadata 与扩展后的 gate / sealoff 报告 |
| 依赖 | `Round 43` |
| 验收 | 上述 5 个 util 子包不再只停留在 `util_full_stack` 聚合层；必须具备独立 target family、独立 metadata key、独立 test / config-route / reference-matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Tools/Basic/test/Basic_test.cpp`、`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`、`Tools/Diagnostics/test/Diagnostics_test.cpp`、`Tools/Output/test/Output_test.cpp`、`Tools/Solver/test/SurfaceSolverFacade_test.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`documentation/contracts/surface_family_coverage_matrix_v1.json`、`scripts/run/surface_spis_family_coverage_matrix.json`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_family_coverage_gate_round44/surface_family_coverage_gate.json`、`build_codex/surface_reference_matrix_gate_round44_full/surface_reference_matrix_gate.json`、`build_codex/surface_round44_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R44-Y1` | 在 coverage matrix / contract 中新增 `util_func`、`util_sampler`、`util_monitor`、`util_io`、`util_matrix` target family，并补齐 family-to-code / test / config-route 映射 | `已完成` | `2026-04-17` | `Round 43` | util 子包 contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`、`scripts/run/surface_spis_family_coverage_matrix.json`、`build_codex/surface_family_coverage_gate_round44/surface_family_coverage_gate.json` | 进入 `R44-Y2` |
| `R44-Y2` | 在 runtime metadata 中新增 `surface_util_func_*`、`surface_util_sampler_*`、`surface_util_monitor_*`、`surface_util_io_*`、`surface_util_matrix_*` 字段，并补 util smoke / reference matrix 断言 | `已完成` | `2026-04-17` | `R44-Y1` | util 子包 family metadata | 单测 + smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_reference_matrix_gate_round44_full/surface_reference_matrix_gate.json` | 进入 `R44-Y3` |
| `R44-Y3` | 重新跑 sealoff，要求扩展后 `util_func`、`util_sampler`、`util_monitor`、`util_io`、`util_matrix` 的 `alignment_status = completed` | `已完成` | `2026-04-17` | `R44-Y2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round44_sealoff/surface_round26_sealoff.json` | 下一轮优先推进剩余 `Util/**` 子包 |

#### Round 45：Wave Z `Util/**` 第二批子包独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 继续将 `spis/Util/**` 中具备稳定代码与测试落点的子包从 `util_full_stack` 聚合层升级为独立 target family，并补齐 `Util/DistribFunc`、`Util/Exception`、`Util/Part`、`Util/Phys`、`Util/Table`、`Util/Vect` 的 runtime metadata / tests / gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Util/DistribFunc/**`、`Util/Exception/**`、`Util/Part/**`、`Util/Phys/**`、`Util/Table/**`、`Util/Vect/**`；`Tools/Basic/**`；`Tools/Geometry/**`；`Tools/Particle/**`；`Tools/Output/**`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |
| 实施项 | 执行 `R45-Z1`~`R45-Z3`，扩展 target family、util 子包 metadata、Basic/Geometry/Particle/Output 单测与 reference-matrix / sealoff 证据 |
| 输出 | `util_distrib_func`、`util_exception`、`util_part`、`util_phys`、`util_table`、`util_vect` target families，以及扩展后的 gate / sealoff 报告 |
| 依赖 | `Round 44` |
| 验收 | 上述 6 个 util 子包不再停留在 `util_full_stack` 聚合层；必须具备独立 target family、独立 metadata key、独立 test / config-route / reference-matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Tools/Basic/test/Basic_test.cpp`、`Tools/Particle/test/ParticleDefinitions_test.cpp`、`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`、`Tools/Geometry/test/Vector3D_test.cpp`、`Tools/Output/test/Output_test.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`documentation/contracts/surface_family_coverage_matrix_v1.json`、`scripts/run/surface_spis_family_coverage_matrix.json`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_family_coverage_gate_round45/surface_family_coverage_gate.json`、`build_codex/surface_reference_matrix_gate_round45_full/surface_reference_matrix_gate.json`、`build_codex/surface_round45_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R45-Z1` | 在 coverage matrix / contract 中新增 `util_distrib_func`、`util_exception`、`util_part`、`util_phys`、`util_table`、`util_vect` target family，并补齐 family-to-code / test / config-route 映射 | `已完成` | `2026-04-17` | `Round 44` | util 子包 contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`、`scripts/run/surface_spis_family_coverage_matrix.json`、`build_codex/surface_family_coverage_gate_round45/surface_family_coverage_gate.json` | 进入 `R45-Z2` |
| `R45-Z2` | 在 runtime metadata 中新增 `surface_util_distrib_func_*`、`surface_util_exception_*`、`surface_util_part_*`、`surface_util_phys_*`、`surface_util_table_*`、`surface_util_vect_*` 字段，并补 util smoke / reference matrix 断言 | `已完成` | `2026-04-17` | `R45-Z1` | util 子包 family metadata | 单测 + smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_reference_matrix_gate_round45_full/surface_reference_matrix_gate.json` | 进入 `R45-Z3` |
| `R45-Z3` | 重新跑 sealoff，要求扩展后 `util_distrib_func`、`util_exception`、`util_part`、`util_phys`、`util_table`、`util_vect` 的 `alignment_status = completed` | `已完成` | `2026-04-17` | `R45-Z2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round45_sealoff/surface_round26_sealoff.json` | 进入 `Round 46` |

#### Round 46：Wave AA `Util/Instrument` 独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `spis/Util/Instrument/**` 从 `util_full_stack` 聚合层升级为独立 target family，并补齐 `util_instrument` 的 runtime metadata / tests / gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Util/Instrument/**`；`Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Main/SurfaceScenarioLoader.cpp`；`Main/main.cpp` |
| 实施项 | 执行 `R46-AA1`~`R46-AA3`，扩展 target family、instrument metadata、observer/instrument smoke 与 reference-matrix / sealoff 证据 |
| 输出 | `util_instrument` target family，以及 `surface_util_instrument_*` metadata 与扩展后的 gate / sealoff 报告 |
| 依赖 | `Round 45` |
| 验收 | `Util/Instrument` 不再停留在 `util_full_stack` 聚合层；必须具备独立 target family、独立 metadata key、独立 test / config-route / reference-matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`Main/SurfaceScenarioLoader.cpp`、`Main/main.cpp`、`documentation/contracts/surface_family_coverage_matrix_v1.json`、`scripts/run/surface_spis_family_coverage_matrix.json`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_family_coverage_gate_round46/surface_family_coverage_gate.json`、`build_codex/surface_reference_matrix_gate_round46_full/surface_reference_matrix_gate.json`、`build_codex/surface_round46_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R46-AA1` | 在 coverage matrix / contract 中新增 `util_instrument` target family，并补齐 family-to-code / test / config-route 映射 | `已完成` | `2026-04-17` | `Round 45` | util instrument contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`、`scripts/run/surface_spis_family_coverage_matrix.json`、`build_codex/surface_family_coverage_gate_round46/surface_family_coverage_gate.json` | 进入 `R46-AA2` |
| `R46-AA2` | 在 runtime metadata 中新增 `surface_util_instrument_*` 字段，并补 observer/instrument smoke / reference matrix 断言 | `已完成` | `2026-04-17` | `R46-AA1` | util instrument family metadata | 单测 + smoke + matrix gate | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_reference_matrix_gate_round46_full/surface_reference_matrix_gate.json` | 进入 `R46-AA3` |
| `R46-AA3` | 重新跑 sealoff，要求扩展后 `util_instrument.alignment_status = completed` | `已完成` | `2026-04-17` | `R46-AA2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round46_sealoff/surface_round26_sealoff.json` | 进入 `Round 47` |

#### Round 47：Wave AB `Util/List + Util/OcTree + Util/OcTreeMesh` 独立 family 化轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 将 `spis/Util/List/**`、`spis/Util/OcTree/**`、`spis/Util/OcTreeMesh/**` 从 `util_full_stack` 聚合层升级为独立 target family，并补齐 `util_list`、`util_octree`、`util_octree_mesh` 的 runtime metadata / tests / gate |
| 输入 | `ref/SPIS/num-master/src/main/java/spis/Util/List/**`、`Util/OcTree/**`、`Util/OcTreeMesh/**`；`Tools/Particle/**`；`Tools/Mesh/**`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json` |
| 实施项 | 执行 `R47-AB1`~`R47-AB3`，扩展 target family、util 子包 metadata、Particle/Mesh 单测与 reference-matrix / sealoff 证据 |
| 输出 | `util_list`、`util_octree`、`util_octree_mesh` target families，以及 `surface_util_list_*`、`surface_util_octree_*`、`surface_util_octree_mesh_*` metadata 与扩展后的 gate / sealoff 报告 |
| 依赖 | `Round 46` |
| 验收 | 三个 util 子包不再停留在 `util_full_stack` 聚合层；必须具备独立 target family、独立 metadata key、独立 test / config-route / reference-matrix 证据 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Tools/Particle/test/ParticleSource_test.cpp`、`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`、`Tools/Mesh/test/MeshPartitioning_test.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`documentation/contracts/surface_family_coverage_matrix_v1.json`、`scripts/run/surface_spis_family_coverage_matrix.json`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_family_coverage_gate_round47/surface_family_coverage_gate.json`、`build_codex/surface_reference_matrix_gate_round47_full/surface_reference_matrix_gate.json`、`build_codex/surface_round47_sealoff/surface_round26_sealoff.json` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R47-AB1` | 在 coverage matrix / contract 中新增 `util_list`、`util_octree`、`util_octree_mesh` target family，并补齐 family-to-code / test / config-route 映射 | `已完成` | `2026-04-17` | `Round 46` | util 子包 contract / matrix | contract 检视 + coverage gate | `documentation/contracts/surface_family_coverage_matrix_v1.json`、`scripts/run/surface_spis_family_coverage_matrix.json`、`build_codex/surface_family_coverage_gate_round47/surface_family_coverage_gate.json` | 进入 `R47-AB2` |
| `R47-AB2` | 在 runtime metadata 中新增 `surface_util_list_*`、`surface_util_octree_*`、`surface_util_octree_mesh_*` 字段，并补 Particle/Mesh 单测、smoke 与 reference matrix 断言 | `已完成` | `2026-04-17` | `R47-AB1` | util 子包 family metadata | 单测 + smoke + matrix gate | `Tools/Particle/test/ParticleSource_test.cpp`、`Tools/Particle/test/SurfaceDistributionFunction_test.cpp`、`Tools/Mesh/test/MeshPartitioning_test.cpp`、`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`、`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`、`scripts/run/surface_reference_matrix.json`、`build_codex/surface_reference_matrix_gate_round47_full/surface_reference_matrix_gate.json` | 进入 `R47-AB3` |
| `R47-AB3` | 重新跑 sealoff，要求扩展后 `util_list`、`util_octree`、`util_octree_mesh` 的 `alignment_status = completed` | `已完成` | `2026-04-17` | `R47-AB2` | 扩展后的 sealoff 结果 | gate | `build_codex/surface_round47_sealoff/surface_round26_sealoff.json` | `Util/**` 子包级拆分全部收口，转入长期回归维护 |

#### Round 48：Wave AC `Top/Default` class-to-runtime 对象层索引轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 以 `surface_spis_class_audit.json` 为基线，补齐 `Top/Default/**` 的 class-to-runtime 对象层职责索引，解释 `NamedObject`、`Parameter`、`SimulationObject`、`GlobalParameter`、`LocalParameter`、`InstrumentParameter` 等 class 在当前项目中的承载落点 |
| 输入 | `Toolkit/Surface Charging/docs/surface_spis_class_audit.json`；`ref/SPIS/num-master/src/main/java/spis/Top/Default/**`；`Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp`；`Toolkit/Surface Charging/include/SurfaceScenarioCatalog.h`；`Toolkit/Surface Charging/src/SurfaceScenarioCatalog.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` |
| 实施项 | 执行 `R48-AC1`~`R48-AC3`，新增 `Top/Default` class-to-runtime 索引文档、回填主文档 class-level 结论，并补一组 object-layer / smoke 引用证据 |
| 输出 | `Top/Default` class-to-runtime 映射文档、更新后的逐 class 审计条目与主计划文档引用 |
| 依赖 | `Round 47`；`surface_spis_class_audit_summary.md` |
| 验收 | `Top/Default` 当前 `25` 个 `近似实现` class 都必须在文档中具备“SPIS class -> 当前对象层/配置层/运行时承载落点 -> 证据”的显式索引；不得再只停留在聚合口径 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Toolkit/Surface Charging/docs/top_default_class_runtime_mapping.md`；`Toolkit/Surface Charging/docs/surface_spis_class_audit.json`；`Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R48-AC1` | 为 `Top/Default` 建立逐 class 对象层索引，按 “SPIS class / 当前 runtime struct / scenario object / config carrier / 证据” 五列输出 | `已完成` | `2026-04-17` | `Round 47` | `Top/Default` 映射文档 | 文档检视 + 代码交叉核对 | `Toolkit/Surface Charging/docs/top_default_class_runtime_mapping.md`；`Toolkit/Surface Charging/include/SurfaceRuntimePlan.h`；`Toolkit/Surface Charging/src/SurfaceRuntimePlan.cpp` | 进入 `R48-AC2` |
| `R48-AC2` | 在主文档 2.14 节回填 `Top/Default` class-level 差异与收口结论，明确其不影响 family sealoff | `已完成` | `2026-04-17` | `R48-AC1` | 文档化 class-level 结论 | 文档检视 | 本文档；`Toolkit/Surface Charging/docs/top_default_class_runtime_mapping.md`；`Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md` | 进入 `R48-AC3` |
| `R48-AC3` | 补 object-layer / smoke 引用证据，使 `Top/Default` 文档映射可回链到测试或运行时入口 | `已完成` | `2026-04-17` | `R48-AC2` | `Top/Default` 证据链 | 测试检视 + smoke 引用复核 | `Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Toolkit/Surface Charging/docs/top_default_class_runtime_mapping.md` | 进入 `Round 49` |

#### Round 49：Wave AD `Top/Top + Top/Simulation` shell-to-runtime 对照轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 为 `Top/Top/**` 与 `Top/Simulation/**` 建立 shell-to-runtime 对照，说明 `UIInvokable`、`SpisTopMenu`、`SimulationListener`、`SimulationFromUIParams` 等壳层/入口层 class 如何被当前 `SurfaceScenarioCatalog`、`SurfaceSimulationRunner`、`SurfaceRuntimePlan` 吸收 |
| 输入 | `Toolkit/Surface Charging/docs/surface_spis_class_audit.json`；`ref/SPIS/num-master/src/main/java/spis/Top/Top/**`；`ref/SPIS/num-master/src/main/java/spis/Top/Simulation/**`；`Toolkit/Surface Charging/include/SurfaceSimulationRunner.h`；`Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`；`Main/SurfaceScenarioLoader.cpp`；`Main/main.cpp` |
| 实施项 | 执行 `R49-AD1`~`R49-AD3`，新增 shell-to-runtime 对照文档、补 scenario/runner 入口证据，并将 `top_top` / `top_simulation` 的 class-level 审计结论回填主文档 |
| 输出 | `Top/Top + Top/Simulation` shell-to-runtime 映射文档、主文档更新、scenario/runner 引用证据 |
| 依赖 | `Round 48` |
| 验收 | `Top/Top` 中 `UIInvokable` / `SpisTopMenu` 与 `Top/Simulation` 中 `SimulationListener` 等 class 必须具备明确的当前承载解释；不能只写“已被吸收”而不给出实际入口文件 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Toolkit/Surface Charging/docs/top_shell_simulation_runtime_mapping.md`；`Toolkit/Surface Charging/docs/surface_spis_class_audit.json`；`Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Main/main.cpp`；`Main/SurfaceScenarioLoader.cpp` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R49-AD1` | 为 `Top/Top` 中 `Scenario`、`PotentialSweep`、`UIInvokable`、`SpisTopMenu`、`NumTopFromUI` 建立 shell-to-runtime 对照 | `已完成` | `2026-04-17` | `Round 48` | `Top/Top` 映射文档 | 文档检视 + 代码检视 | `Toolkit/Surface Charging/docs/top_shell_simulation_runtime_mapping.md`；`Toolkit/Surface Charging/include/SurfaceScenarioCatalog.h`；`Toolkit/Surface Charging/src/SurfaceScenarioCatalog.cpp`；`Main/SurfaceScenarioLoader.cpp` | 进入 `R49-AD2` |
| `R49-AD2` | 为 `Top/Simulation` 中 `Simulation`、`PlasmaScSimulation`、`SimulationFromUIParams`、`SimulationListener` 建立 runner-to-runtime 对照 | `已完成` | `2026-04-17` | `R49-AD1` | `Top/Simulation` 映射文档 | 文档检视 + 代码检视 | `Toolkit/Surface Charging/docs/top_shell_simulation_runtime_mapping.md`；`Toolkit/Surface Charging/include/SurfaceSimulationRunner.h`；`Toolkit/Surface Charging/src/SurfaceSimulationRunner.cpp`；`Toolkit/Surface Charging/include/SurfaceRuntimePlan.h` | 进入 `R49-AD3` |
| `R49-AD3` | 在主文档与逐 class 审计结果中回填 `Top/Top + Top/Simulation` class-level 结论与引用证据 | `已完成` | `2026-04-17` | `R49-AD2` | 文档化 shell-to-runtime 结论 | 文档检视 | 本文档；`Toolkit/Surface Charging/docs/top_shell_simulation_runtime_mapping.md`；`Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md` | 进入 `Round 50` |

#### Round 50：Wave AE `Util/Func` class-to-helper / formula 映射轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 为 `Util/Func/**` 建立 class-to-helper / formula 映射索引，将当前 `98` 个 `近似实现` class 分为“已由基础函数库聚合承载”“需显式命名桥接”“不建议逐类镜像”三类 |
| 输入 | `Toolkit/Surface Charging/docs/surface_spis_class_audit.json`；`ref/SPIS/num-master/src/main/java/spis/Util/Func/**`；`Tools/Basic/include/MathUtils.h`；`Tools/Basic/include/Constants.h`；`Tools/Basic/include/NumericAggregation.h`；`Tools/Basic/include/StringTokenUtils.h`；`Tools/Basic/include/ModelRegistry.h` |
| 实施项 | 执行 `R50-AE1`~`R50-AE3`，新增 `Util/Func` class-to-helper 映射文档、按三类归档结果回填审计结论，并列出需要后续显式命名桥接的函数族 |
| 输出 | `Util/Func` class-to-helper / formula 映射索引、桥接候选清单、更新后的主文档引用 |
| 依赖 | `Round 49` |
| 验收 | `Util/Func` 下所有 class 都必须被归入三类之一；不得再以“整体聚合到 MathUtils”一语带过 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Toolkit/Surface Charging/docs/util_func_class_helper_mapping.md`；`Tools/Basic/include/MathUtils.h`；`Tools/Basic/include/Constants.h`；`Tools/Basic/include/NumericAggregation.h`；`Tools/Basic/test/Basic_test.cpp`；`Tools/Particle/src/ParticleSource.cpp`；`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`；`Toolkit/Surface Charging/docs/surface_spis_class_audit.json`；`Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R50-AE1` | 为 `Util/Func/**` 生成逐 class 映射索引，并标记 “已聚合承载 / 需显式命名桥接 / 不建议逐类镜像” | `已完成` | `2026-04-17` | `Round 49` | `Util/Func` 映射索引 | 文档检视 + 代码检视 | `Toolkit/Surface Charging/docs/util_func_class_helper_mapping.md`；`Tools/Basic/include/MathUtils.h`；`Tools/Basic/include/Constants.h`；`Toolkit/Surface Charging/docs/surface_spis_class_audit.json` | 进入 `R50-AE2` |
| `R50-AE2` | 输出函数族桥接候选清单，优先标出物理相关、边界模型相关、分布函数相关 helper | `已完成` | `2026-04-17` | `R50-AE1` | 桥接候选清单 | 文档检视 | `Toolkit/Surface Charging/docs/util_func_class_helper_mapping.md`；`Tools/Basic/include/NumericAggregation.h`；`Tools/Basic/include/ModelRegistry.h`；`Tools/Particle/include/SurfaceDistributionFunction.h`；`Tools/Particle/src/ParticleSource.cpp` | 进入 `R50-AE3` |
| `R50-AE3` | 在主文档 2.14 节回填 `Util/Func` class-level 结论，明确哪些差异属于结构优化而非功能缺口 | `已完成` | `2026-04-17` | `R50-AE2` | 主文档更新 | 文档检视 | 本文档；`Toolkit/Surface Charging/docs/util_func_class_helper_mapping.md`；`Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md` | 进入 `Round 51` |

#### Round 51：Wave AF `Util/Instrument` catalog / observer-route 映射轮

| 字段 | 内容 |
| --- | --- |
| 目标 | 为 `Util/Instrument/**` 建立 instrument catalog / observer-route / smoke 断言映射，将当前 `51` 个 `近似实现` class 分为“已有 runtime 行为”“仅 metadata route”“尚未单独可观测”的三类 |
| 输入 | `Toolkit/Surface Charging/docs/surface_spis_class_audit.json`；`ref/SPIS/num-master/src/main/java/spis/Util/Instrument/**`；`Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Main/SurfaceScenarioLoader.cpp`；`Main/main.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp` |
| 实施项 | 执行 `R51-AF1`~`R51-AF3`，新增 instrument catalog 映射文档、补 observer/instrument 测试矩阵索引，并在主文档中回填 class-level 结论 |
| 输出 | `Util/Instrument` catalog / route 映射文档、测试矩阵索引、更新后的主文档引用 |
| 依赖 | `Round 50` |
| 验收 | `Util/Instrument` 当前所有 class 都必须被归入三类之一，并给出 `Main` / `DensePlasmaSurfaceCharging` / smoke 证据中的至少一种；不得再只停留在 family metadata 口径 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-17` |
| 证据 | `Toolkit/Surface Charging/docs/util_instrument_catalog_observer_mapping.md`；`Main/main.cpp`；`Main/SurfaceScenarioLoader.cpp`；`Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/src/SurfaceTransitionEngine.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`Toolkit/Surface Charging/docs/surface_spis_class_audit.json`；`Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R51-AF1` | 为 `Util/Instrument/**` 建立逐 class catalog / observer-route 映射，并按 “已有 runtime 行为 / 仅 metadata route / 尚未单独可观测” 归类 | `已完成` | `2026-04-17` | `Round 50` | `Util/Instrument` 映射文档 | 文档检视 + 代码检视 | `Toolkit/Surface Charging/docs/util_instrument_catalog_observer_mapping.md`；`Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Main/SurfaceScenarioLoader.cpp` | 进入 `R51-AF2` |
| `R51-AF2` | 补 observer/instrument smoke 与引用测试矩阵索引，使 class-level 结论可回链到现有测试 | `已完成` | `2026-04-17` | `R51-AF1` | 测试矩阵索引 | 测试检视 + smoke 引用复核 | `Toolkit/Surface Charging/docs/util_instrument_catalog_observer_mapping.md`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`Toolkit/Surface Charging/test/SurfaceCharging_object_layer_test.cpp`；`scripts/run/surface_reference_matrix.json` | 进入 `R51-AF3` |
| `R51-AF3` | 在主文档 2.14 节回填 `Util/Instrument` class-level 结论，明确其为“后续结构增强主线”而非 “family remediation 主线” | `已完成` | `2026-04-17` | `R51-AF2` | 主文档更新 | 文档检视 | 本文档；`Toolkit/Surface Charging/docs/util_instrument_catalog_observer_mapping.md`；`Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md` | 转入 class-level 长期回归维护 |

#### Round 52：严格 per-source 粒子账本推广轮
| 字段 | 内容 |
| --- | --- |
| 目标 | 将当前 source-resolved 输出从“单 source 兼容 + 导出层补充”升级为“多 source 同时存在时的严格 per-source 粒子账本”，并保证输出结果与 SPIS 多 source 语义一致 |
| 输入 | `Toolkit/Surface Charging/include/PicMccSurfaceCurrentSampler.h`；`Toolkit/Surface Charging/src/PicMccSurfaceCurrentSampler.cpp`；`Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/scripts/spis_surface_import_lib.py`；`Toolkit/Surface Charging/scripts/compare_spis_java_output.py`；`Toolkit/Surface Charging/scripts/augment_spis_source_resolved_output.py`；`Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json` |
| 实施项 | 执行 `R52-AG1`~`R52-AG4`，统一 sampler/source slot 记账、补 strict per-source runtime ledger、推广 SPIS import/compare 为多 source 自动生成，并为 source/node/source/spacecraft/all-source 聚合关系补 gate |
| 输出 | 严格 per-source 账本文档、统一 runtime/source 输出、无 `source1` 硬编码 compare/import 逻辑、多 source smoke 与 matrix 证据 |
| 依赖 | `Round 51`；`results/spis_import/child_langmuir/*`；`Toolkit/Surface Charging/docs/strict_per_source_particle_ledger_plan.md` |
| 验收 | 单 source 与多 source 共享同一条主路径；运行时不再依赖 `source_keys.front()` 或 `source_keys.size()==1`；`surface_source_*` 与 `surface_node_*_source_*` 由原子账本推导；compare/import 不再写死 `source1`；gate 能验证 per-source 聚合守恒与与 SPIS 多 source 输出一致 |
| 当前状态 | `已完成` |
| 最后更新 | `2026-04-18` |
| 证据 | `Toolkit/Surface Charging/docs/strict_per_source_particle_ledger_plan.md`；`Toolkit/Surface Charging/src/PicMccSurfaceCurrentSampler.cpp`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/scripts/spis_surface_import_lib.py`；`Toolkit/Surface Charging/scripts/compare_spis_java_output.py`；`Toolkit/Surface Charging/scripts/augment_spis_source_resolved_output.py`；`Toolkit/Surface Charging/scripts/run_child_langmuir_spis_workflow.py`；`results/spis_import/child_langmuir_regen_workflow/child_langmuir.comparison.md` |

| 任务ID | 任务 | 状态 | 最后更新 | 依赖 | 产出 | 验证方式 | 证据 | 下一步 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R52-AG1` | 重构 `PicMccSurfaceCurrentSampler` 的 source slot 与 `source_resolved_samples` 生成逻辑，消除 `config.source_keys.front()` 单 source 假设 | `已完成` | `2026-04-18` | `Round 51` | 多 source sampler 原子样本 | 单测 + 代码检视 | `Toolkit/Surface Charging/include/PicMccSurfaceCurrentSampler.h`；`Toolkit/Surface Charging/src/PicMccSurfaceCurrentSampler.cpp`；`Toolkit/Surface Charging/docs/strict_per_source_particle_ledger_plan.md` | 进入 `R52-AG2` |
| `R52-AG2` | 在 `DensePlasmaSurfaceCharging` 中建立严格 per-source 原子账本，统一 `node/source`、`source/spacecraft` 与 `all-source` 聚合路径，并补 `interactor_current` 严格归因或显式未归因语义 | `已完成` | `2026-04-18` | `R52-AG1` | strict runtime ledger + metadata mode | 单测 + smoke + 输出契约检视 | `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`；`Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`；`Toolkit/Surface Charging/docs/strict_per_source_particle_ledger_plan.md` | 进入 `R52-AG3` |
| `R52-AG3` | 推广 SPIS import / compare 脚本为多 source 通用生成逻辑，移除 `source1` 专用 comparison target 与 mapping table | `已完成` | `2026-04-18` | `R52-AG2` | 多 source compare/import 通用脚本 | 脚本回归 + `child_langmuir` 对比复核 | `Toolkit/Surface Charging/scripts/spis_surface_import_lib.py`；`Toolkit/Surface Charging/scripts/compare_spis_java_output.py`；`results/spis_import/child_langmuir/child_langmuir.comparison.csv` | 进入 `R52-AG4` |
| `R52-AG4` | 为严格 per-source 账本补 smoke / reference matrix / output contract gate，校验 source 聚合守恒与 metadata 完整性 | `已完成` | `2026-04-18` | `R52-AG3` | strict per-source bookkeeping gate | smoke + matrix gate + contract gate | `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`；`scripts/run/surface_reference_matrix.json`；`scripts/python/check_surface_reference_matrix_gate.py`；`scripts/python/check_output_contract_v1_gate.py`；`Toolkit/Surface Charging/docs/strict_per_source_particle_ledger_plan.md` | 转入实现轮 |

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
  - 当前状态
  - 关键差异
  - 对齐方向
- 不允许把“最小闭环”“最小插件层”“主线可用”误写成“与 SPIS family 已 1:1 对齐”。
- 不允许把 surface bridge 契约冻结误写成 native volume kernel 全家族对标完成。
- 不允许把 facade 落地误写成 solver family 语义已全部对齐。

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
- `Round 20` 之后的正文主线必须以 `Wave A`~`Wave G` 为核心，不得再使用“标准A/标准B 已完成”替代 family 级盘点。

### 6.3 执行验收

- 任意一轮结束后，本文档都能准确显示：
  - 哪些任务已完成
  - 哪些任务进行中
  - 哪些任务被冻结
  - 哪些任务受阻
  - 哪些任务已经转入历史归档
- 任意一个 SPIS family 在被标记为 `已等价实现` 前，本文档都必须同时给出：
  - 代码落点
  - 配置或运行路由落点
  - 测试或 gate 证据
  - 与 SPIS 对应源路径
- 任意一个 SPIS class 若在逐 class 审计中被标记为 `已等价实现` 或 `聚合等价实现`，至少必须能回链到以下证据之一：
  - 代码落点
  - config/route 落点
  - 测试或 smoke
  - gate / reference matrix / sealoff
- 若某个 SPIS class 被标记为 `近似实现`，本文档必须明确说明它属于：
  - family 已封板后的聚合承载差异
  - 还是需要重新触发 family 级复审的功能缺口

### 6.4 兼容验收

- 必须维持以下边界可持续兼容：
  - `SurfaceChargingConfig`
  - preset 名称与别名
  - `SurfaceRuntimeRoute`
  - `SurfacePicStrategy`
  - 主 CSV 输出
  - 现有 sidecar 合同
- 若为了 1:1 对齐新增内部 family、对象层或元数据，必须保持上述外部边界默认兼容，除非文档先记录破坏式变更批准。

### 6.5 回归验收

- Surface 主线继续以以下路径为主 gate：
  - `ref/GEO/**`
  - `ref/LEO/**`
  - `ref/mat/**`
  - `scripts/run/surface_reference_matrix.json`
- `Round 20` 必须额外补齐：
  - SPIS family coverage matrix
  - family-to-code 映射
  - family-to-test 映射
  - family-to-config-route 映射
- `Round 21`~`Round 25` 必须额外补齐：
  - 对应 family 的单测覆盖
  - 主线路由 smoke
  - 必要时的 sidecar/metadata 可观测字段
  - 若影响参考矩阵，则同步补齐 matrix gate
- `Round 26` 必须额外补齐：
  - long-tail reference matrix 场景
  - direct-SPIS semantic comparison cases
  - family coverage gate
  - 最终 1:1 对齐封板证据

## 7. 附录：当前已完成的 Surface 基础对齐成果摘要

### 7.1 已沉淀到 `Tools/*` 的核心表面能力

- `Tools/Material`
  - `SurfaceMaterialModel`
  - `BasicSurfaceMaterialModel`
  - `ErosionSurfaceMaterialModel`
  - `ErosionParamSet` 与材料参数资产路由
  - 表面材料响应与 legacy yield/helper 的公共化接口
- `Tools/Particle`
  - `SurfaceDistributionFunction`
  - `TabulatedVelocity`、`NonLocalizedHybrid`、`MultipleSurf`、`LocalModifiedPearsonIV`
  - `LocalTabulated`、`TwoAxesTabulatedVelocity`、`UniformVelocity`、`Fluid`、`FowlerNordheim`
  - `SurfaceFluxSampler`
  - `SurfaceEmissionSampler`
- `Tools/FieldSolver`
  - `SurfaceBarrierModels`
  - `BarrierCurrentScaler`
  - `SecondaryRecollection/VariableBarrier/OML/LTE/FN/Automatic/Multiple`
  - `GlobalTemp/SmoothedGlobalTemp/CurrentVariation/LocalVariation`
  - `SurfaceFieldVolumeBridge`（native volume parity 最小闭环）
- `Tools/Coupling`
  - `SurfaceCircuitCoupling`
  - `SurfaceCurrentLinearizationKernel`
  - DIDV 高阶组合与规约链（`Reducable` + matrix/reduction route）
- `Tools/Solver`
  - `SurfaceSolverFacade`
- 当前仍未独立沉淀的对象层：`SurfInteract` 插件家族（见 2.12 与 `Round 19`）。

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
