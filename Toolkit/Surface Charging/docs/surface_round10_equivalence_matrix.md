# Surface Round 10 功能等价矩阵（执行版）

## 1. R10-A 基本一致口径

| 结论 | 定义 | 最低证据要求 |
| --- | --- | --- |
| `已一致` | 能力/模式在当前实现中可稳定复现，且与 SPIS 对标语义一致 | 至少 1 条代码路径 + 1 条 smoke/object-layer 用例 + 1 个 artifact/contract 证据 |
| `近似一致` | 主行为一致，但在边界条件、收敛策略或观测粒度上仍有可解释差异 | 至少 1 条 reference matrix 或 smoke 证据 + 偏差说明 |
| `待补齐` | 对标能力尚未形成可重复 gate，或关键指标未通过 | 失败 gate 或缺失项清单 |
| `工程替代` | 保持外部契约兼容，但内部实现为工程替代而非 SPIS 同构实现 | 兼容 contract + 差异声明 |

### 1.1 口径补充（与 Round 11 衔接）

- 本文中的 `已一致` 默认指 **Round 10 场景级/契约级口径**：当前主线在既定 preset、route、artifact contract 与 reference matrix gate 上可稳定复现。
- `已一致` 不等同于 SPIS class-family 全量同构；Surf/Circ/Top 的模型家族深度差距在 `surface_num_alignment_refactor_plan.md` 的 2.9 与 Round 11 任务卡中继续跟踪。
- 后续若 Round 11 新增能力进入 gate，本矩阵应同步把“场景级一致”与“家族级补齐”证据分开记录。

## 2. R10-B Surface 能力矩阵

| 能力域 | SPIS 对标源 | 当前实现落点 | 当前结论 | 证据 |
| --- | --- | --- | --- | --- |
| 材料响应（yield/导电/电容） | `ref/SPIS/num-master/src/main/java/spis/Surf/**` | `Tools/Material/include/SurfaceMaterialModel.h` | `已一致` | `SurfaceChargingSmokeTest.SurfaceMaterialInteractionScalesReachCapacitanceAndCurrentMainline`；`SurfaceChargingSmokeTest.SurfaceMaterialLibraryImportOverridesPrimaryMaterial` |
| 粒子分布与采样 | `ref/SPIS/num-master/src/main/java/spis/Surf/**` | `Tools/Particle/include/SurfaceDistributionFunction.h` | `已一致` | `SurfaceChargingSmokeTest.DistributionModelSwitchRebalancesReferencePatchCurrents`；`SurfaceChargingSmokeTest.SurfacePicRouteExportsUnifiedKernelSnapshotMetadata` |
| 场/势垒与回收缩放 | `ref/SPIS/num-master/src/main/java/spis/Surf/**` | `Tools/FieldSolver/include/SurfaceBarrierModels.h` | `已一致` | `SurfaceChargingSmokeTest.ReferenceModelBodyPatchBarrierContrastStaysWithinExpectedRanges` |
| 电路耦合与 DIDV | `ref/SPIS/num-master/src/main/java/spis/Circ/**` | `Tools/Coupling/include/SurfaceCircuitCoupling.h` | `已一致` | `SurfaceChargingSmokeTest.ReferenceDidvMatchesFiniteDifferenceAcrossGeoLeoRoutes`；`SurfaceChargingSmokeTest.GeoLeoMatBenchmarkArtifactsExposeDidvMetric` |
| legacy benchmark execute/replay | `ref/SPIS/num-master/src/main/java/spis/Top/**` | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已一致` | `SurfaceChargingSmokeTest.LegacyBenchmarkExecuteModeAdvancesWithoutReplay`；`SurfaceTransitionEngineTest.LegacyReplayRequiresLegacyRouteAndReplayMode` |
| shared runtime 全局耦合一致性 | `ref/SPIS/num-master/src/main/java/spis/Top/**` | `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp` | `已一致` | `SurfaceChargingSmokeTest.SharedSurfacePicRuntimeGlobalCoupledSolveConvergesWithinSubstep`；`build/surface_reference_matrix_gate.json`（cases 全量 `PASS`） |
| bridge/sidecar 契约 | `ref/SPIS/num-master/src/main/java/spis/Top/**` | `documentation/contracts/**`、`SurfaceSimulationRunner` | `已一致` | `SurfaceChargingSmokeTest.ExternalFieldSolverBridgeIngestsNodeResult`；`SurfaceChargingSmokeTest.ExternalVolumeSolverBridgeIngestsCellResult`；`SurfaceChargingSmokeTest.SurfaceReferenceAndRuntimeContractsExportMetadata` |

## 3. R10-C 运行模式与输出矩阵

| 运行模式 | 输入入口 | 关键行为 | 关键输出 | 当前结论 | 证据 |
| --- | --- | --- | --- | --- | --- |
| `surface`（统一主线） | preset + runner | 初始化/步进/导出一体化 | 主 CSV + sidecar | `已一致` | `SurfaceChargingSmokeTest.ToolkitInitializesAdvancesAndExports` |
| `surface-config`（显式参数覆盖） | CLI 覆盖字段 | body/photo/material 覆盖生效 | 输出元数据与主序列同步 | `已一致` | `SurfaceChargingSmokeTest.SurfaceConfigCliAppliesBodyAndPhotoOverrides` |
| preset + alias | 场景名/别名解析 | 主线 preset 与 replay preset 分离 | metadata 中 source 名称稳定 | `已一致` | `SurfaceScenarioCatalogTest.ListsPresetsAndResolvesAliases`；`SurfaceChargingSmokeTest.MainlinePresetLookupExcludesLegacyDeckReplayPresets` |
| legacy replay | Legacy route + replay mode | 复现参考曲线，不进入 execute path | benchmark/report sidecar | `已一致` | `SurfaceTransitionEngineTest.LegacyReplayRequiresLegacyRouteAndReplayMode` |
| legacy execute | Legacy route + execute mode | 用当前 runtime 执行，不依赖 replay 数据 | benchmark/report + metrics | `已一致` | `SurfaceChargingSmokeTest.LegacyBenchmarkExecuteModeAdvancesWithoutReplay` |
| shared-coupled surface-pic | shared runtime + global solve | patch 间耦合求解与收敛统计 | shared runtime observer/consistency | `已一致` | `SurfaceChargingSmokeTest.SharedSurfacePicRuntimeGlobalCoupledSolveConvergesWithinSubstep`；`build/surface_reference_matrix_gate.json` |
| bridge sidecar | external field/volume 输入 | 读取 bridge 输入并映射至边界/体域 | field/volume sidecar 与 contract | `已一致` | `SurfaceChargingSmokeTest.ExternalFieldSolverBridgeIngestsBoundaryGroupResult`；`SurfaceChargingSmokeTest.ExternalVolumeSolverBridgeIngestsBoundaryGroupResult` |

## 4. R10-D 参考场景等价矩阵（基于当前 gate 资产）

数据来源：`build/surface_reference_matrix_gate.json`、`scripts/run/surface_reference_matrix.json`

| 参考场景/族 | 代表 preset | 当前状态 | 说明 |
| --- | --- | --- | --- |
| GEO C 参考族 | `geo_ecss_kapton_ref` | `已一致` | command/metadata/benchmark/observer/consistency/global domain/repository/sheath 全项 `PASS` |
| GEO SPIS-PIC direct 族 | `geo_ecss_kapton_surface_pic_direct` | `已一致` | 参考矩阵 case 全项 `PASS`，一致性/全局域/仓库/鞘层子项均通过 |
| LEO RAM C 参考族 | `leo_ref_ram_facing` | `已一致` | 关键 benchmark 与 observer 子项 `PASS` |
| LEO Wake C 参考族 | `leo_ref_wake_facing` | `已一致` | 关键 benchmark 与 observer 子项 `PASS` |
| LEO PIC-hybrid 族 | `leo_pic_circuit_ram_facing_hybrid` | `已一致` | shared-runtime 一致性、global domain/repository/sheath 子项均 `PASS` |
| LEO PIC-hybrid multi-patch 族 | `leo_pic_circuit_ram_facing_hybrid_multi_patch` | `已一致` | 多 patch shared-runtime 门禁子项 `PASS`，总状态 `PASS` |
| MATLAB 兼容参考族 | `geo_ecss_kapton_ref_matlab_compatible` | `已一致` | MATLAB 兼容参考路径与 benchmark contract 子项 `PASS` |

## 5. R10-E 未一致项收口记录

| 历史条目 | 本轮状态 | 收口证据 | 备注 |
| --- | --- | --- | --- |
| `R11-SHARED-001` | `已关闭` | `build/surface_reference_matrix_gate.json`（`geo_ecss_kapton_surface_pic_direct` 全项 `PASS`） | 由 Round 10 gate 刷新直接收口 |
| `R11-SHARED-002` | `已关闭` | `build/surface_reference_matrix_gate.json`（`leo_pic_circuit_ram_facing_hybrid` 全项 `PASS`） | global domain/repository/sheath 全部通过 |
| `R11-SHARED-003` | `已关闭` | `build/surface_reference_matrix_gate.json`（`leo_pic_circuit_ram_facing_hybrid_multi_patch` 全项 `PASS`） | 多 patch 一致性门禁通过 |
| `R11-GATE-001` | `已关闭` | `build/surface_reference_matrix_gate.md`（report 总状态与 case 状态一致） | 汇总语义与 case 子状态一致 |

## 6. Round 11 衔接清单（差距转任务）

| 差距主题 | Round 10 当前口径 | Round 11 任务锚点 |
| --- | --- | --- |
| Surf material 家族深度（含 erosion） | 场景级 `已一致`（`R11-A` 已完成并通过子集+全量 matrix gate） | `R11-A` |
| Surf distrib 家族深度（tabulated velocity / non-localized） | 场景级 `已一致`（`R11-B` 已完成） | `R11-B` |
| Surf scaler 家族深度（OML/LTE/FN） | 场景级 `已一致`（`R11-C` 已完成） | `R11-C` |
| Circ 波形族（PWL/EXP）与 DIDV 语义族 | 场景级 `已一致`（`R11-D` 已完成） | `R11-D` |
| Top 事件级 transition（LocalTime/Spinning/SunFlux） | 场景级 `已一致`（`R11-E` 已完成） | `R11-E` |
| Top 扩展事件族（ConductivityEvolution/参数更新） | 场景级 `已一致`（最小事件链已完成） | `Round 12` 候选 |
| 新增能力 gate 证据刷新 | `R11` 子集 gate 与全量 `surface_reference_matrix` gate 均已 `PASS`（见 `build/round11_surface_reference_matrix_gate/*` 与 `build/surface_reference_matrix_gate/*`） | `R11-F` |

- `Round 11` 封板结论：`R11-A/B/C/D/E/F` 全部完成，且 `build/round11_surface_reference_matrix_gate/surface_reference_matrix_gate.json` 与 `build/surface_reference_matrix_gate.json` 均为 `PASS`。

## 7. 本文档更新策略

- 每次刷新 `surface_reference_matrix` 后，仅允许在本文件更新：
  - 场景状态
  - 偏差说明
  - 收口条目（若存在）
- 若 gate 全量 `PASS`，保持 backlog 为空并记录“已收口”；若出现回归失败，再新增下一轮 backlog 条目。
