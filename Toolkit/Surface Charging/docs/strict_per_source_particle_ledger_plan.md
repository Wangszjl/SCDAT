# 多 source 严格 per-source 粒子账本改造计划

## 1. 目标

- 将当前 `Surface Charging` 中“source 只在导出层可见”的弱 source-resolved 路径，升级为运行时核心就严格按 `source_key` 记账的粒子账本。
- 保证单 source 与多 source 场景都遵循同一套账本语义，不再依赖 `source_keys.front()`、`source_keys.size() == 1` 等单 source 特判。
- 使输出结果能够与 SPIS 多 source 输出保持一一对应，并能够通过 machine-readable gate 验证：
  - `per-node + per-source` 原子账本成立；
  - `per-source -> spacecraft` 聚合成立；
  - `all sources -> total` 聚合成立；
  - SCDAT 与 SPIS 的 source-resolved 监测量命名、数量、语义一致。

## 2. 当前问题

### 2.1 运行时 source 维度尚未成为核心一等对象

- `PicMccCurrentSample` 已暴露 `source_resolved_samples`，但 `PicMccSurfaceCurrentSampler.cpp` 当前只写入 `config.source_keys.front()`，本质仍是单 source 采样结果外包一层 source 标签。
- `DensePlasmaSurfaceCharging.cpp` 已经能消费 `config_.spis_import.source_keys` 并导出部分 `surface_source_*` 序列，但其上游采样并没有真正提供严格的多 source 原子账本。

### 2.2 interactor 电流没有严格 per-source 归因

- `history_surface_source_spacecraft_interactor_currents_` 当前只有在 `source_keys.size() == 1` 时，才用总二次/背散射等发射电流回填。
- 这意味着一旦出现多个 source，`interactor_current` 既不守恒，也不能与 SPIS 的 per-source interactor 输出严格对齐。

### 2.3 SPIS 导入与 compare 仍然写死 `source1`

- `spis_surface_import_lib.py` 当前只生成 `source1` 对应的 comparison targets。
- `compare_spis_java_output.py` 当前映射表也只支持 `Individual_current_on_spacecraft_Collected_source1`、`Number_of_source1` 等固定字段。
- 这导致现在的 `child_langmuir` 结果更像“单一案例适配”，还不是可推广的通用 source-resolved 契约。

### 2.4 缺少严格账本不变量 gate

- 当前输出虽有 `surface_source_*` 序列，但尚未形成“原子账本 -> 聚合账本 -> 对外输出”的严格不变量校验链。
- 缺少对以下关系的统一 gate：
  - `sum(node_source_collected_current_a) == source_spacecraft_collected_current_a`
  - `sum(source_spacecraft_collected_current_a) == spacecraft_total_collected_current_a`
  - `sum(source_superparticle_count) == global/live repository total`
  - `sum(source interactor) == all-source interactor aggregate`

## 3. 目标语义

### 3.1 原子账本维度

后续所有 source-resolved 输出，必须从以下原子账本维度推导：

- `time_step`
- `node_index`
- `source_key`
- `collected_current_a`
- `interactor_current_a`
- `alive_superparticle_count`

如实现阶段能稳定获得更多可审计量，可继续扩展：

- `injected_superparticle_count`
- `absorbed_superparticle_count`
- `emitted_superparticle_count`
- `escaped_superparticle_count`
- `net_charge_c`

### 3.2 归因规则

- `collected_current_a`
  - 按入射并被表面收集的粒子自身 `source_key` 记账。
- `alive_superparticle_count`
  - 按时间步结束后仍存活于 source-resolved 粒子仓库中的粒子 `source_key` 记账。
- `interactor_current_a`
  - 按“触发表面相互作用的入射 source_key”归因，而不是按发射粒子种类重新分类。
  - 不能再把总 secondary/backscatter/ion-secondary 电流在多 source 情形下模糊并账。
- 若某一类相互作用暂时无法严格归因：
  - 必须显式输出 `unattributed_*` 或 `unsupported_*` 语义；
  - 不允许静默并入某个 source；
  - 不允许继续维持“单 source 时正确，多 source 时空缺”的隐式行为。

### 3.3 聚合规则

- 原子账本是唯一事实来源。
- 一切聚合量必须由原子账本求和得到：
  - `source x spacecraft`
  - `all-source x node`
  - `all-source x spacecraft`
- 禁止单独走旁路公式生成看似同名的 source-resolved 对外字段。

## 4. 设计原则

### 4.1 内部账本优先于导出兼容

- 第一优先级是建立严格的运行时 source 账本。
- 导出字段、compare 脚本、SPIS 适配应全部建立在该账本之上。
- 不接受“输出层拼接出看似正确的 source 字段，但底层运行时并无严格归因”的方案。

### 4.2 单 source 与多 source 同构

- 单 source 场景必须视为多 source 的特例。
- 任何使用 `front()`、`size()==1`、`if only one source then backfill` 的逻辑都应从主路径移除。

### 4.3 保持既有外部接口兼容

- 保留既有 `surface_source_*` 命名风格与 `surface_spis_import_source_key_*` metadata。
- 允许新增更严格的 metadata 与 gate 字段，但不破坏现有 CSV / sidecar 主字段。

## 5. 实施分期

### Phase A：账本模型收口

目标：

- 将 source 维度前移到 sampler / live PIC 运行时核心。
- 明确原子账本的数据结构和归因规则。

实施项：

1. 为 `PicMccSourceResolvedSample` 增补严格账本语义注释与必要字段。
2. 梳理 `PicMccSurfaceCurrentSampler` 的 source 输入来源，建立稳定的 `source_key -> source slot` 映射。
3. 消灭 `config.source_keys.front()` 的单 source 写法，改为真正逐 source 采样/汇总。
4. 为 `interactor_current` 设计严格归因路径，或显式落地 `unattributed` 兼容字段。

阶段验收：

- `PicMccCurrentSample.source_resolved_samples` 在多 source 输入下包含全部 source。
- 不再依赖单 source 特判。
- 单测能够证明 source slot 顺序稳定、缺失 source 返回零样本而不是丢字段。

### Phase B：DensePlasmaSurfaceCharging 运行时账本统一

目标：

- 用统一的原子账本驱动 node/source history、spacecraft/source history 与总量 history。

实施项：

1. 将 `history_surface_node_source_collected_currents_`、`history_surface_source_*` 收敛到同一账本来源。
2. 新增 `history_surface_node_source_interactor_currents_`，避免 spacecraft 层只有 source 聚合、节点层却没有原子事实。
3. 为 superparticle count 明确“瞬时存活数”还是“累计注入数”，并统一命名。
4. 在 runtime metadata 中新增账本模式标识，例如：
   - `surface_source_resolved_bookkeeping_mode = strict_per_source_particle_ledger_v1`

阶段验收：

- source/node、source/spacecraft、total 三层输出都由统一账本生成。
- 可以自动检查聚合一致性。

### Phase C：SPIS 导入与 comparison 通用化

目标：

- 从“固定支持 source1”推广为“自动发现全部 source_keys 并生成 comparison targets”。

实施项：

1. 更新 `spis_surface_import_lib.py`：
   - 自动扫描 SPIS case 中可用 source 列表；
   - 生成每个 source 的 `Collected / Interactor / Number` comparison target；
   - 生成每个 source x node 的 node-level collected/interactor target。
2. 更新 `compare_spis_java_output.py`：
   - 按 `target.source_key`、`target.source_family`、`scdat_series_hint` 动态解析；
   - 移除硬编码的 `source1` mapping table。
3. 保持 `child_langmuir` 单 source case 输出不回退。

阶段验收：

- compare 脚本可直接消费任意数量 source。
- `child_langmuir` 继续 PASS 当前 source1 对齐。

### Phase D：Gate 与 reference matrix 收口

目标：

- 将严格 per-source 粒子账本纳入现有 matrix / contract / smoke 验收链。

实施项：

1. 在 smoke test 中增加多 source 场景：
   - 至少两个 source；
   - 至少两个节点；
   - 校验原子账本与聚合账本一致。
2. 在 reference matrix 中新增 per-source bookkeeping case。
3. 在 output contract / reference matrix gate 中新增以下断言：
   - source key metadata 完整；
   - required `surface_source_*`、`surface_node_*_source_*` 字段完整；
   - per-source sum 与 total sum 一致；
   - strict bookkeeping mode 与期望值一致。

阶段验收：

- 新 gate 能独立判定账本是否严格成立。
- reference matrix 可以明确给出 `strict_per_source_particle_ledger_v1` 证据。

## 6. 任务拆分建议

### Round 52A：Sampler/source slot 核心改造

- 责任范围：
  - `Toolkit/Surface Charging/include/PicMccSurfaceCurrentSampler.h`
  - `Toolkit/Surface Charging/src/PicMccSurfaceCurrentSampler.cpp`
- 交付：
  - 多 source `source_resolved_samples`
  - 单测补齐

### Round 52B：Dense runtime ledger 统一

- 责任范围：
  - `Toolkit/Surface Charging/include/DensePlasmaSurfaceCharging.h`
  - `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`
- 交付：
  - 严格原子账本
  - source/node/source/spacecraft 输出统一
  - metadata 模式升级

### Round 52C：SPIS import / compare 通用化

- 责任范围：
  - `Toolkit/Surface Charging/scripts/spis_surface_import_lib.py`
  - `Toolkit/Surface Charging/scripts/compare_spis_java_output.py`
  - `Toolkit/Surface Charging/scripts/augment_spis_source_resolved_output.py`
- 交付：
  - 多 source comparison target 自动生成
  - 无 `source1` 专用映射

### Round 52D：Smoke / matrix / gate 收口

- 责任范围：
  - `Toolkit/Surface Charging/test/SurfaceCharging_smoke_test.cpp`
  - `scripts/run/surface_reference_matrix.json`
  - `scripts/python/check_surface_reference_matrix_gate.py`
  - `scripts/python/check_output_contract_v1_gate.py`
- 交付：
  - 严格账本 gate
  - 多 source smoke / matrix evidence

## 7. 风险与注意事项

- 最大风险不是代码改不出来，而是“语义看似完整、实际归因不严格”。
- `interactor_current` 的 source 归因必须优先澄清，否则多 source 对齐会继续停留在表面字段一致。
- `superparticle_count` 必须区分“瞬时存活数”和“累计产生数”，不要让 SPIS compare 与内部仓库统计使用不同口径。
- 若 SPIS 某些 source-resolved 场当前只存在于 ASCII numkernel sidecar，而不在标准 nc 输出中，需要在文档中明确这是“导入适配边界”，不能混写成运行时天然对齐。

## 8. 完成判据

- 运行时主路径中不存在单 source 特判账本逻辑。
- `source_resolved_samples` 真正反映所有 source，而不是仅包装 `front()`。
- `surface_source_*` 与 `surface_node_*_source_*` 全部由原子账本求和得到。
- compare / import 脚本不再出现 `source1` 硬编码。
- 至少一个多 source smoke case 和一个 matrix case 通过 gate。
- 文档、测试、gate 三者都明确声明并验证：
  - `strict_per_source_particle_ledger_v1`

