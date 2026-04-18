# SPIS Import Compute Workflow

> 说明：这份文档描述的是 `SPIS -> 本项目` 的辅助适配与对标链路，不是 surface 模块的唯一主流程。
> 原生建模入口、属性配置结构、仿真启动主线，请优先参见 [native_surface_mainline_workflow.md](./native_surface_mainline_workflow.md)。

## 1. 目的

本文档整理当前项目中“从 `ref/SPIS` 算例导入，到 `surface-config` 数值计算，再到与 SPIS 结果对比”的实际计算链路。

这条链路的目标不是只服务 `Child_Langmuir` 和 `Plasma_Wake` 两个案例。

- `Child_Langmuir`
- `Plasma_Wake`

这两个案例只是当前已经跑通的验证样例，用来证明这条链路在现有工程中是可执行的。

对你自己的新算例，目标也是复用同一条主链路：

1. 从 SPIS case 解析结构与参数
2. 生成本项目可运行的 `surface-config` JSON
3. 调用 `SCDAT.exe surface-config` 完成计算
4. 导出 CSV / metadata / sidecar
5. 用 NumKernel 输出补齐严格 per-source 粒子账本
6. 与 SPIS Java `nc` 输出做字段级对比

## 2. 总体流程

```mermaid
flowchart TD
    A[SPIS case 根目录\nmodel.xml / groups.xml / globalParameters.xml / NumKernel / OutputFolder]
    B[Importer 生成 surface-config JSON\nbuild_spis_case_documents()\nbuild_child_langmuir_documents()\nbuild_plasma_wake_documents()]
    C[解析 SPIS 参数\nparse_model_paths()\nparse_global_parameters()\nparse_groups()\nparse_pwl()\ndiscover_source_keys()]
    D[生成 comparison_targets / import_manifest / scan_plan\nbuild_source_resolved_comparison_targets()]
    E[Workflow 调用 SCDAT.exe\nrun_child_langmuir_spis_workflow.main()\nrun_plasma_wake_spis_workflow.main()]
    F[SCDAT surface-config 主计算]
    G[运行时采样\nPicMccSurfaceCurrentSampler.sample()\nsampleWithDerivative()]
    H[per-source 账本聚合\nDensePlasmaSurfaceCharging\nsource_resolved_samples -> history]
    I[导出 CSV + metadata + sidecar\nDensePlasmaSurfaceCharging.exportResults()]
    J[NumKernel 文本增强 source-resolved 列\naugment()\nbuild_augmented_series()]
    K[Java nc vs SCDAT CSV 对比\ncompare()\nresolve_scdat_series()]
    L[comparison.csv / comparison.md / comparison.json]

    A --> B
    B --> C
    C --> D
    D --> E
    E --> F
    F --> G
    G --> H
    H --> I
    I --> J
    J --> K
    K --> L
```

## 3. 分阶段链路

### 3.1 SPIS Case 导入层

职责：

- 从 `.spis5` 目录中提取模型路径、组定义、全局参数、NumKernel 输出根目录、Java 输出根目录
- 识别 source keys
- 识别哪些 source 有 Java `nc` 参考输出，哪些 source 只有 NumKernel 文本输出
- 生成本项目使用的 `surface-config` JSON、`import_manifest.json`、`scan_plan.json`

核心入口：

- [build_spis_case_documents()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:783)
- [build_child_langmuir_documents()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:328)
- [build_plasma_wake_documents()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:601)
- [main()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/import_spis_surface_case.py:11)

关键解析函数：

- [parse_model_paths()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:57)
- [parse_global_parameters()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:46)
- [parse_groups()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:82)
- [parse_pwl()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:136)
- [_material_from_groups()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:156)
- [discover_source_keys()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:212)
- [build_source_resolved_comparison_targets()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:255)

这一层的输出不是最终结果，而是“可运行输入”。

## 3.2 Workflow 执行层

职责：

- 调用 importer
- 定位 `SCDAT.exe`
- 执行 `surface-config`
- 对结果做 augment
- 对结果做 compare

当前两个样例入口：

- [run_child_langmuir_spis_workflow.py](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/run_child_langmuir_spis_workflow.py)
- [run_plasma_wake_spis_workflow.py](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/run_plasma_wake_spis_workflow.py)

关键函数：

- `Child_Langmuir`：[main()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/run_child_langmuir_spis_workflow.py:26)
- `Plasma_Wake`：[main()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/run_plasma_wake_spis_workflow.py:26)
- 公共调用器：[run_command()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/run_child_langmuir_spis_workflow.py:13)
- 可执行定位：[locate_scdat()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/run_child_langmuir_spis_workflow.py:19)

这一层只是 orchestration，不负责物理计算本身。

## 3.3 SCDAT 运行时计算层

职责：

- 按 importer 生成的 `surface-config` 进行数值推进
- 对 live PIC / source-resolved 电流进行采样
- 在 `DensePlasmaSurfaceCharging` 内聚合为历史序列
- 导出 CSV 和 metadata

关键入口对象：

- [PicMccSurfaceCurrentSampler](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/include/PicMccSurfaceCurrentSampler.h:81)
- [sample()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/src/PicMccSurfaceCurrentSampler.cpp:635)
- [sampleWithDerivative()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/src/PicMccSurfaceCurrentSampler.cpp:641)
- source-resolved 样本字段 [source_resolved_samples](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/include/PicMccSurfaceCurrentSampler.h:78)

账本聚合与导出：

- source/node interactor 历史 [history_surface_node_source_interactor_currents_](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/include/DensePlasmaSurfaceCharging.h:1357)
- bookkeeping 完整性标志 [source_resolved_bookkeeping_slots_complete_](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/include/DensePlasmaSurfaceCharging.h:1361)
- 每步 source-resolved 消费主逻辑 [DensePlasmaSurfaceCharging.cpp:5904](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/src/DensePlasmaSurfaceCharging.cpp:5904)
- 导出主入口 [exportResults()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/src/DensePlasmaSurfaceCharging.cpp:11348)
- node/source 系列导出 [DensePlasmaSurfaceCharging.cpp:11660](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/src/DensePlasmaSurfaceCharging.cpp:11660)
- `spis_import` / source-resolved metadata 导出 [DensePlasmaSurfaceCharging.cpp:11810](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/src/DensePlasmaSurfaceCharging.cpp:11810)

这一层是项目内部真正的数值主线。

## 3.4 Source-Resolved 增强层

职责：

- 读取 SCDAT 导出的 CSV
- 读取 SPIS NumKernel 文本输出
- 将 source-resolved `collected_current`、`interactor_current`、`superparticle_count` 回填到 CSV
- 同步写入 source-resolved metadata

关键函数：

- [augment()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/augment_spis_source_resolved_output.py:215)
- [build_augmented_series()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/augment_spis_source_resolved_output.py:163)
- [parse_current_matrix()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/augment_spis_source_resolved_output.py:76)
- [load_count_series_for_source()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/augment_spis_source_resolved_output.py:137)

这一层的意义是：把“运行时导出的 source-resolved 主账本”和“SPIS NumKernel 的 source-resolved 参考账本”合并到统一结果文件中。

## 3.5 对比层

职责：

- 读取 SPIS Java `DataFieldMonitored` / `DataFieldExtracted` 的 `nc` 文件
- 读取当前项目导出的 CSV
- 按 target 配置做字段级映射
- 输出 `comparison.csv`、`comparison.md`、`comparison.json`

关键函数：

- [compare()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/compare_spis_java_output.py:328)
- [resolve_target_file()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/compare_spis_java_output.py:103)
- [resolve_scdat_series()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/compare_spis_java_output.py:189)
- [build_source_resolved_checks()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/compare_spis_java_output.py:151)

这一层输出的是最终诊断结果，而不是计算结果本身。

## 4. 当前已跑通样例在链路中的位置

### 4.1 Child_Langmuir

特点：

- 有 PWL 偏压扫描
- 有 discrete point cross-check
- source-resolved Java `nc` 输出主要集中在 `source1`
- 多个其他 source 只在 NumKernel 文本中可见

对应入口：

- [build_child_langmuir_documents()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:328)
- [run_child_langmuir_spis_workflow.py](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/run_child_langmuir_spis_workflow.py)

当前已跑通结果：

- [results/spis_import/child_langmuir_regen_workflow](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/results/spis_import/child_langmuir_regen_workflow)

### 4.2 Plasma_Wake

特点：

- 单节点静态 timeline
- 无 PWL / discrete sweep
- source 主要是 `elec1` 与 `ions1`
- Java `nc` 中已经有两类 source 的 collected / number 输出

对应入口：

- [build_plasma_wake_documents()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:601)
- [run_plasma_wake_spis_workflow.py](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/run_plasma_wake_spis_workflow.py)

当前已跑通结果：

- [results/spis_import/plasma_wake_workflow](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/results/spis_import/plasma_wake_workflow)

## 5. 对你自己的算例意味着什么

最重要的一点：

`Child_Langmuir` 和 `Plasma_Wake` 不是“专属流程”，而是“已验证流程”。

你的新算例接入时，原则上仍然走同一条总链路：

1. 新算例 `.spis5` 进入 importer
2. importer 为该算例生成本项目可运行配置
3. workflow 调 `SCDAT.exe surface-config`
4. `DensePlasmaSurfaceCharging` 导出结果
5. augment 用 NumKernel 补 source-resolved 列
6. compare 用 Java `nc` 做结果对比

真正需要为新算例适配的，通常只有 importer 侧的“case 映射层”，而不是后面的计算主线。

也就是说：

- `PicMccSurfaceCurrentSampler`
- `DensePlasmaSurfaceCharging`
- `augment_spis_source_resolved_output.py`
- `compare_spis_java_output.py`

这些层应该尽量保持通用。

新算例最常见的改动位置是：

- [build_spis_case_documents()](/E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/scripts/spis_surface_import_lib.py:783)
- 新增一个类似 `build_<your_case>_documents()` 的 builder

## 6. 新算例接入的推荐检查单

如果后面要接你的新算例，建议按这个检查单推进。

### 6.1 SPIS 输入是否齐全

- `model.xml`
- `groups.xml`
- `globalParameters.xml`
- `OutputFolder/DataFieldMonitored`
- `OutputFolder/DataFieldExtracted`
- `NumKernel/Output`

### 6.2 importer 是否能回答这几个问题

- 这个 case 是静态还是 sweep
- 有几个 electrical node
- 哪些 source 是真实 source
- 哪些 source 有 Java `nc` 参考
- 哪些 source 只有 NumKernel 文本
- surface node / branch / bias 应如何映射到本项目配置

### 6.3 compare target 是否来自真实参考输出

不能把 importer 写成“按 source_keys 全量盲目生成 compare target”。

必须区分：

- 真实 source
- Java-reference-backed source
- NumKernel-only source

否则新算例一旦 source 更多，就会产生大量假 `unsupported`。

### 6.4 成功接入的判据

一个新算例被认为“接入成功”，至少应满足：

1. importer 能生成合法配置
2. `surface-config` 能跑完
3. CSV / metadata / sidecar 正常导出
4. augment 能补齐 source-resolved 列
5. compare 能生成 `comparison.csv/.md/.json`
6. `unsupported = 0` 或者所有 unsupported 都有明确原因

## 7. 一句话总结

当前项目的 SPIS 对接流程，本质上是：

“SPIS case 解析层” + “SCDAT 通用计算主线” + “source-resolved 增强层” + “SPIS 结果对比层”。

`Child_Langmuir` 和 `Plasma_Wake` 只是这条通用流程已经跑通的两个样例，不应被理解成流程只能服务这两个案例。
