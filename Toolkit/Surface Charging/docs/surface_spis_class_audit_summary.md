# Surface/SPIS 逐 Class 审计摘要

## 1. 基线

- 生成时间：`2026-04-17T07:03:33+00:00`
- 审计范围：`Surface Charging 主线边界内六大域逐 class 审计`
- SPIS class 总数：`602`
- sealoff 基线：`tracked_family_count = 42`，`final_alignment_verdict = achieved`
- 参考主文档：`Toolkit/Surface Charging/docs/surface_num_alignment_refactor_plan.md`

## 2. 产物

- 机器可读主表：`Toolkit/Surface Charging/docs/surface_spis_class_audit.json`
- 可筛选 CSV：`Toolkit/Surface Charging/docs/surface_spis_class_audit.csv`
- family matrix：`scripts/run/surface_spis_family_coverage_matrix.json`
- coverage gate：`build_codex/surface_family_coverage_gate_round47/surface_family_coverage_gate.json`
- sealoff：`build_codex/surface_round47_sealoff/surface_round26_sealoff.json`

## 3. 域级统计

| 域 | class 数 | 已等价实现 | 聚合等价实现 | 近似实现 | 尚无实现 | 不纳入范围 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Surf | 100 | 55 | 45 | 0 | 0 | 0 |
| Vol | 69 | 54 | 15 | 0 | 0 | 0 |
| Top | 65 | 17 | 20 | 28 | 0 | 0 |
| Circ | 29 | 21 | 8 | 0 | 0 | 0 |
| Solver | 23 | 7 | 16 | 0 | 0 | 0 |
| Util | 315 | 0 | 166 | 149 | 0 | 0 |
| osgi | 1 | 0 | 0 | 0 | 0 | 1 |

## 4. 解读规则

- `已等价实现`：当前项目里有较明确的单类或强命名落点，且证据链可回溯。
- `聚合等价实现`：语义已闭环，但由多个 C++ 模块共同承载，不是 Java 单类镜像。
- `近似实现`：family 已有 route/metadata/gate 证据，但 class 级对象层或独立语义承载偏弱。
- `尚无实现`：当前边界内未找到足够证据。
- `不纳入当前目标范围`：当前审计不纳入 Surface Charging 主线。

> 说明：如果某个 class 被标为 `聚合等价实现` 或 `近似实现`，并不自动推翻当前 family sealoff；这表示“family 已封板，但 class 级对象层仍有聚合承载差异”。

## 5. 后续任务摘要

- `top_default`: `25` 个 class 仍是 `近似实现/尚无实现`。Top/Default 在 family gate 上已封板，但 class 级对象层镜像最弱，建议单独补一份对象层职责索引。
- `top_simulation`: `1` 个 class 仍是 `近似实现/尚无实现`。优先补类级职责索引、显式对象层映射和更直接的单测/route 证据。
- `top_top`: `2` 个 class 仍是 `近似实现/尚无实现`。Top/Top 的 UI/top-shell 类在当前项目中被 runtime plan 和 scenario runner 吸收，适合补一份 shell-to-runtime 对照。
- `util_func`: `98` 个 class 仍是 `近似实现/尚无实现`。Util/Func 大量 Java 函数对象已折叠为基础数学/物理 helper，若需要 class 级可信度，应补命名桥接与函数索引。
- `util_instrument`: `51` 个 class 仍是 `近似实现/尚无实现`。Util/Instrument 目前更接近 family route 完成，尚未形成逐 instrument 的 C++ 类型与专门单测矩阵。

## 6. 高风险热点

- `Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`
- `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`
- `Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`

这些文件仍是 class-level 语义被集中承载的主要热点，后续如果继续追求更强的逐 class 可审计性，应优先从这里拆对象层索引。
