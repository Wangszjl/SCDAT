# SPIS Case Onboarding Guide

## 1. 目的

本文档说明：如果你要把你自己的 SPIS 算例接入当前项目，并复用现有

- importer
- `surface-config` 计算主线
- source-resolved augment
- SPIS compare

应该按什么顺序修改、检查和验证。

本文档不是只针对 `Child_Langmuir` 或 `Plasma_Wake`。

这两个案例只是已经跑通的参考样例。你的新算例应尽量复用同一条总流程。

## 2. 总原则

新算例接入时，优先修改 importer 层，而不是修改计算主线。

也就是说，优先改这些地方：

- [build_spis_case_documents()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:783)
- 在 [spis_surface_import_lib.py](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py) 中新增一个 `build_<your_case>_documents()` 函数

尽量不要为了某个新算例去特化这些通用层：

- [PicMccSurfaceCurrentSampler](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/include/PicMccSurfaceCurrentSampler.h:81)
- [DensePlasmaSurfaceCharging](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/include/DensePlasmaSurfaceCharging.h:1083)
- [augment_spis_source_resolved_output.py](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/augment_spis_source_resolved_output.py)
- [compare_spis_java_output.py](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/compare_spis_java_output.py)

只有当新算例暴露出真正的通用能力缺口时，才应扩展这些公共层。

## 3. 接入步骤

### 3.1 准备输入

先确认你的 SPIS 算例目录至少包含：

- `model.xml`
- `DefaultStudy/Preprocessing/Groups/groups.xml`
- `DefaultStudy/Simulations/Run1/GlobalParameters/globalParameters.xml`
- `DefaultStudy/Simulations/Run1/OutputFolder/DataFieldMonitored`
- `DefaultStudy/Simulations/Run1/OutputFolder/DataFieldExtracted`
- `DefaultStudy/Simulations/Run1/NumKernel/Output`

如果这些路径不齐，现有 importer 很难直接复用。

### 3.2 先回答 6 个问题

在写 builder 之前，先把这 6 件事搞清楚：

1. 这个 case 是静态运行，还是 sweep / PWL / 多工况扫描
2. 有几个 electrical node
3. `surface_nodes` 在当前项目里应该如何映射
4. 是否需要 `surface_branches`
5. 哪些 source 是真实 source
6. 哪些 source 有 Java `nc` 参考，哪些 source 只有 NumKernel 文本

如果这 6 件事没有想清楚，后面的 `comparison_targets` 很容易生成错。

### 3.3 新增 builder

在 [spis_surface_import_lib.py](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py) 中新增：

```python
def build_your_case_documents(case_root: pathlib.Path, output_root: pathlib.Path) -> SpisCaseArtifacts:
    ...
```

建议按下面结构写：

1. `parse_model_paths()`
2. `parse_global_parameters()`
3. `parse_groups()`
4. 如有必要，解析 sweep / bias / circuit
5. `discover_source_keys()`
6. 构建 `comparison_targets`
7. 构建 `import_manifest`
8. 构建 `scan_plan`
9. 构建 `surface-config` JSON
10. 返回 `SpisCaseArtifacts`

### 3.4 更新 case 分派入口

在 [build_spis_case_documents()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:783) 里注册你的 case：

```python
if "your_case_name" in case_name:
    return build_your_case_documents(case_root, output_root)
```

如果你的 case 命名模式不稳定，建议改成更稳妥的分派条件，不要只依赖目录名字符串包含关系。

### 3.5 决定 workflow 形式

如果你的 case 是：

- 单个静态运行

那么可以照 `Plasma_Wake` 的模式，新增一个最小 workflow：

- [run_plasma_wake_spis_workflow.py](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/run_plasma_wake_spis_workflow.py)

如果你的 case 是：

- PWL / discrete sweep / 多阶段扫描

那么可以照 `Child_Langmuir` 的模式：

- [run_child_langmuir_spis_workflow.py](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/run_child_langmuir_spis_workflow.py)

## 4. `comparison_targets` 怎么设计

这是最容易出错的一层。

### 4.1 不要按 `source_keys` 盲目全量生成 compare target

必须区分：

- `source_keys`
- `comparison_target_source_keys`
- `numkernel_only_source_keys`

原因是：

- 有些 source 在 SPIS Java `nc` 输出里存在
- 有些 source 只在 NumKernel 文本里存在

如果你不区分，compare 结果里会出现很多假 `unsupported`。

### 4.2 推荐最小 target 集

通常先从这几类开始：

- `Average_surface_potential_of_node_*`
- `Total_current_on_node_*`
- `Total_current_on_spacecraft`
- `Individual_current_on_spacecraft_Collected_<source>`
- `Collected_current_<source>_node_<n>`
- `Interactor_current_<source>_node_<n>`：仅当 Java 输出真的存在
- `Number_of_<source>`

### 4.3 target 的原则

只有当 SPIS 参考输出中真的存在匹配文件时，才把它加入 `comparison_targets`。

当前这条逻辑已经在这里实现：

- [build_source_resolved_comparison_targets()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:255)

## 5. source-resolved 接入规则

### 5.1 真实 source 发现

优先使用：

- Java `OutputFolder` 中的 source-resolved 输出
- NumKernel `Number_of_*` 文本

对应函数：

- [discover_source_keys()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:212)

### 5.2 聚合伪 source 要过滤

例如：

- `all_interactions`
- `all_populations`

这类不是严格 per-source 账本里的真实 source，不应该作为真实 source key 写进 compare / importer 主集合。

### 5.3 NumKernel-only source 也是合法 source

即使某个 source 没有 Java `nc` 参考输出，只要它真实存在于 NumKernel 账本中，它也应该：

- 出现在 `source_keys`
- 进入 augment 后的 source-resolved metadata
- 在 CSV 中拥有自己的 source-resolved 列

只是它不一定出现在 Java compare targets 里。

## 6. 验证顺序

推荐按这个顺序做，不要一上来就跑全流程。

### 6.1 先跑 importer

命令示意：

```powershell
python "Toolkit/Surface Charging/scripts/import_spis_surface_case.py" "<your_case_root>" --output-root "<your_output_root>"
```

先检查：

- `generated/*.surface.json`
- `generated/*.import_manifest.json`
- `generated/*.scan_plan.json`

重点看：

- `source_keys`
- `comparison_target_source_keys`
- `numkernel_only_source_keys`
- `comparison_targets`

### 6.2 再跑 `surface-config`

确认 `SCDAT.exe` 能把你的配置跑完，并成功导出：

- CSV
- metadata sidecar

### 6.3 再跑 augment

确认：

- source-resolved 列是否被补齐
- `surface_source_resolved_*` metadata 是否完整

### 6.4 最后跑 compare

看这几个指标：

- `aligned`
- `adapted`
- `unsupported`
- `missing_required_series_count`

理想状态下：

- `unsupported = 0`
- `missing_required_series_count = 0`

## 7. 成功接入的判据

你的新算例接入成功，至少应满足以下 7 条：

1. importer 能稳定生成配置
2. `surface-config` 能正常计算完成
3. CSV / metadata / sidecar 正常导出
4. source-resolved augment 正常执行
5. compare 正常生成 `comparison.csv/.md/.json`
6. `comparison_targets` 不含明显伪 target
7. source-resolved checks 中没有缺列

## 8. 常见失败点

### 8.1 把 case 特性写死在公共脚本里

例如只按 `Child_Langmuir` 的 PWL 假设写 importer，会导致静态 case 无法复用。

### 8.2 把所有 source 都当成 Java 可比对 source

这会制造大量假 `unsupported`。

### 8.3 只补 compare，不补 augment

这样 source-resolved CSV 列不完整，对比再准也只是部分对齐。

### 8.4 为个别 case 去改运行时主线

除非确认是通用能力缺口，否则不要轻易改：

- `PicMccSurfaceCurrentSampler`
- `DensePlasmaSurfaceCharging`

更常见的正确做法，是在 importer 中把 case 映射清楚。

## 9. 推荐落地方式

如果你后面要接自己的算例，我建议按这个节奏：

1. 先复制一个最接近的样例 builder
2. 把 importer 跑通
3. 看 manifest / targets 是否合理
4. 再跑 workflow
5. 最后根据 compare 结果修 importer 和 target 映射

不要一开始就试图把所有 source、所有 monitor、所有 extracted 字段一次性接满。

## 10. 一句话原则

新 SPIS 算例接入时，优先做“case 到通用链路的映射”，而不是为新 case 重新发明一条计算链路。
