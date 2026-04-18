# 原生 Surface 构建与验证注意事项

## 目的

这份文档记录本次为验证原生 surface 主链路时遇到的构建问题、定位过程和可复用处理方式。

它不是架构文档，而是偏工程操作的补充说明。

## 1. 本次实际验证对象

本次验证的目标是：

- [SurfaceCharging_object_layer_test.cpp](</E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/test/SurfaceCharging_object_layer_test.cpp:25>)

验证结果：

- `10/10` 测试通过

## 2. 首个阻塞点

最初 `SurfaceCharging_object_layer_test` 卡在链接阶段，缺失的主要是 surface artifact/export 相关符号，例如：

- `writeSurfaceMonitorJson(...)`
- `writeSharedSurfaceRuntimeObserverJson(...)`
- `exportSurfaceObserverArtifacts(...)`
- `exportSurfaceBridgeArtifactSuite(...)`

对应实现文件是：

- [SurfaceObserverArtifactExport.cpp](</E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/src/SurfaceObserverArtifactExport.cpp:1>)
- [SurfaceBridgeArtifactWriters.cpp](</E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/src/SurfaceBridgeArtifactWriters.cpp:1>)

## 3. 关键定位结论

定位后发现：

1. 这两个 `.cpp` 文件本身是存在的
2. 重新 `cmake -S . -B build` 后，`build.ninja` 已经包含这两个源文件
3. 但旧的 `libToolkit_Surface_Charging.a` 仍然没有把对应对象归档进去

因此，第一个问题不是“代码缺实现”，而是：

- **静态库归档陈旧**

## 4. 本次有效处理方式

### 4.1 先重新配置

```powershell
cmake -S . -B build
```

### 4.2 强制重建 surface 模块静态库

```powershell
cmake --build build --target Toolkit_Surface_Charging --config Debug --clean-first
```

这一步会让 `SurfaceObserverArtifactExport.cpp` 和
`SurfaceBridgeArtifactWriters.cpp` 真正重新进入 `Toolkit_Surface_Charging`。

## 5. 第二个阻塞点

在继续构建过程中，又遇到更底层的问题：

- `Tools_Output.a`
- `Tools_Particle.a`

这两个静态库里出现了坏掉的归档成员，典型现象是：

- `nm` 报 `file format not recognized`
- linker / `ar` / `ranlib` 报 `file truncated`

受影响的对象包括：

- `HDF5Exporter.cpp.obj`
- `ParticlePusher.cpp.obj`

这说明问题已经下沉到：

- **底层静态库归档稳定性**

而不再是 surface 模块本身。

## 6. 本次有效恢复步骤

### 6.1 确保输出目录存在

```powershell
New-Item -ItemType Directory -Force -Path 'build/lib','build/bin' | Out-Null
```

### 6.2 删除损坏的静态库与关键对象

```powershell
Remove-Item -Force 'build/lib/libTools_Output.a','build/lib/libTools_Particle.a' -ErrorAction SilentlyContinue
Remove-Item -Force 'build/Tools/CMakeFiles/Tools_Output.dir/Output/src/HDF5Exporter.cpp.obj' -ErrorAction SilentlyContinue
Remove-Item -Force 'build/Tools/CMakeFiles/Tools_Particle.dir/Particle/src/ParticlePusher.cpp.obj' -ErrorAction SilentlyContinue
```

### 6.3 单线程重建受影响基础库

```powershell
cmake --build build --target Tools_Output Tools_Particle --config Debug -j 1
```

### 6.4 再构建 object-layer test

```powershell
cmake --build build --target SurfaceCharging_object_layer_test --config Debug -j 1
```

### 6.5 运行测试

```powershell
build/bin/SurfaceCharging_object_layer_test.exe
```

## 7. 本次最终状态

最终 `SurfaceCharging_object_layer_test` 已成功：

- 构建
- 链接
- 运行

并得到：

- `10/10` passing

## 8. 对项目状态的解释

这次验证后的结论应当区分两层：

### A. 对 surface 原生主链路的结论

- 原生主链路已经存在
- object-layer 自动化验证已经通过

### B. 对构建系统稳定性的结论

- 当前 workspace 的静态库归档在个别构建轮次下不够稳定
- 问题主要体现在：
  - 旧归档未及时刷新
  - 个别 `.a` 成员损坏

所以，这次发现的是：

- **主链路实现没问题**
- **构建稳定性还有治理空间**

## 9. 后续建议

如果后面要把 object-layer test 升级成稳定 gate，建议继续做两件事：

1. 把“重新配置 + 强制重建”纳入标准验证步骤
2. 继续排查 MinGW `ar/ranlib` 与当前 Ninja 增量构建下的静态库损坏问题

## 10. 推荐配合阅读

- [native_surface_architecture_overview.md](./native_surface_architecture_overview.md)
- [native_surface_mainline_workflow.md](./native_surface_mainline_workflow.md)

## 11. 2026-04 CMake Hardening Follow-up

The original recovery above proved the native object-layer path was healthy, but it still relied on
manual rebuild discipline. We have now added build-system hardening so the recovery is no longer
just an ad hoc operator procedure.

Applied changes:

- top-level [CMakeLists.txt](</E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/CMakeLists.txt:31>)
  now creates shared `bin/` and `lib/` output directories at configure time
- the same file also defines `scdat_output_dirs` and `scdat_harden_build_target(...)`
- Ninja builds now route link/archive work through a single-slot `scdat_link_pool`
- [Tools/CMakeLists.txt](</E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Tools/CMakeLists.txt:63>)
  hardens generated `Tools_*` static libraries
- [Toolkit/CMakeLists.txt](</E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/CMakeLists.txt:63>)
  hardens generated `Toolkit_*` static libraries
- [Main/CMakeLists.txt](</E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Main/CMakeLists.txt:42>)
  hardens `SCDAT` and `main_*` executables
- [Toolkit/Surface Charging/test/CMakeLists.txt](</E:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/Toolkit/Surface%20Charging/test/CMakeLists.txt:17>)
  hardens the surface smoke/object-layer test targets

Post-hardening verification:

```powershell
cmake -S . -B build
cmake --build build --target SurfaceCharging_object_layer_test --config Debug --clean-first
build/bin/SurfaceCharging_object_layer_test.exe
```

Observed result:

- clean-first rebuild completed successfully without the previous manual directory bootstrap
- `SurfaceCharging_object_layer_test` linked successfully
- runtime result remained `10/10` passing

This does not prove every possible target is immune to future MinGW archive issues, but it does
move the surface native verification path from "manual recovery required" to "CMake-hardened and
repeatably rebuildable in the current workspace".
