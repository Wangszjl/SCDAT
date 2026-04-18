# `Util/Func` class-to-helper / formula 映射

## 1. 目的

本文用于落实 `Round 50 / Wave AE`，把 `ref/SPIS/num-master/src/main/java/spis/Util/Func/**` 下的 `98` 个 SPIS class，逐个映射到当前项目中的基础函数库、物理公式承载层、运行时桥接层或结构性替代层。

本文回答的不是 “`util_func` family 是否已经封板”，而是：

- 哪些 `Util/Func` class 已经由 `Tools/Basic` 与现有 helper 聚合承载；
- 哪些 class 仍值得后续补“显式命名桥接”；
- 哪些 class 本质上是 Java 风格函数接口 / 组合器 / 可调用包装器，不建议在当前 C++ 项目中逐类镜像。

> 口径说明：`util_func` family 已在 family 级 gate / reference matrix / sealoff 中封板。  
> 本文只处理 class-level 映射解释，不推翻现有 `completed` 结论。

## 2. 当前承载骨架

| 当前承载层 | 角色 | 主要文件 |
| --- | --- | --- |
| 数学基础层 | 三角函数、指数对数、插值、积分、梯度、特殊函数等通用数学 helper | `Tools/Basic/include/MathUtils.h` |
| 常量与单位层 | 数学常量、物理常量、单位换算、容差、宏函数 | `Tools/Basic/include/Constants.h` |
| 数值聚合层 | 求和、均值、梯形积分等简单数值统计 | `Tools/Basic/include/NumericAggregation.h` |
| 运行时注册/别名层 | 文本别名、模型注册与轻量包装 | `Tools/Basic/include/ModelRegistry.h` |
| 分布与粒子层 | Maxwellian、Kappa、幂律、采样与表面分布相关显式实现 | `Tools/Particle/include/SurfaceDistributionFunction.h`、`Tools/Particle/src/ParticleSource.cpp` |
| 场/边界/几何层 | Dipole、Solenoid、距离、网格/Octree 相关几何或场 helper | `Tools/FieldSolver/src/VolField.cpp`、`Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`、`Tools/Mesh/src/MeshPartitioning.cpp`、`Tools/Geometry/include/Point3D.h` |
| 材料/鞘层 | OML、sheath、Maxwellian yield、边界物理公式 | `Tools/Material/include/SurfaceMaterialModel.h`、`Tools/Material/src/SurfaceMaterialModel.cpp`、`Tools/Boundary/src/BoundaryConditions.cpp` |

## 3. 判定规则

- `已聚合承载`
  - 当前已有明确 helper、公式族或测试入口承载语义；
  - 不要求保留与 SPIS 同名的独立 C++ class。
- `需显式命名桥接`
  - 当前已有相近运行时语义或公式散落在 domain helper 中；
  - 但 class 名义与证据链还不够显式，后续值得补桥接文档或命名包装。
- `不建议逐类镜像`
  - 主要是函数接口、组合器、谓词、UI-callable 包装器、单位检查接口等结构壳层；
  - 当前 C++ 更适合用模板、`std::function`、配置字段、注册表和测试替代。

## 4. 逐 class 映射

| SPIS class | 当前主要承载落点 | 结论 | 说明 |
| --- | --- | --- | --- |
| `Abs` | `MathUtils.h` + `std::abs` | `已聚合承载` | 绝对值函数已由基础数学层直接承载。 |
| `And` | 条件组合由配置分支、普通布尔表达式承担 | `不建议逐类镜像` | 属于布尔组合壳层，不值得单独保留 class。 |
| `ArcCos` | `MathUtils.h` + `std::acos` | `已聚合承载` | 反余弦能力已在通用数学层具备。 |
| `Atan` | `MathUtils.h` + `std::atan` | `已聚合承载` | 反正切能力已在通用数学层具备。 |
| `BoundedFunctionOf3Scal` | `std::function` + clamp/range 约束 | `不建议逐类镜像` | 是函数签名/边界约束壳层，C++ 里不必建同名类。 |
| `BVect` | `VolField.cpp`、`SurfaceFieldVolumeBridge.cpp` | `需显式命名桥接` | 磁场向量语义存在，但尚未形成与 SPIS 同名的显式 helper。 |
| `CombinationOfTwoScalFunctionOfScal` | `std::function` 组合 + runtime policy | `不建议逐类镜像` | 组合器语义可由 lambda 和策略组合完成。 |
| `CombineInt` | `NumericAggregation.h` + 局部组合逻辑 | `不建议逐类镜像` | 更像整型组合辅助器，而非独立数值模型。 |
| `ComputeIntegral` | `MathUtils::integrate1D`、`gaussLegendreIntegrate`、`adaptiveIntegrate` | `已聚合承载` | 通用积分器已显式存在。 |
| `ComputeIntegralMoment` | `MathUtils` 积分 + `NumericAggregation` | `已聚合承载` | moment 类积分由现有积分/聚合能力组合承担。 |
| `ComputeIntegralVectorMoment` | `MathUtils` 积分 + 几何向量运算 | `已聚合承载` | 向量矩积分语义已能组合表达。 |
| `ConstantFunction` | `Constants.h` + lambda 常值函数 | `已聚合承载` | 常值函数无需单独对象层。 |
| `ConvertibleUnitScalFunctionOfScal` | `Constants.h` 单位换算宏与常量 | `已聚合承载` | 单位转换与标量函数已被常量/宏体系承载。 |
| `Correlation` | `NumericAggregation.h` + 本地统计逻辑 | `已聚合承载` | 相关性统计属聚合数学能力。 |
| `Cos` | `MathUtils.h` + `std::cos` | `已聚合承载` | 已由基础数学层直接承载。 |
| `CosMultiD` | `MathUtils` + 多维向量/坐标运算 | `已聚合承载` | 多维余弦场可由现有数学与几何组合表示。 |
| `DebyeSurfSheath` | `Constants.h` Debye 常量、材料/鞘层 runtime | `需显式命名桥接` | Debye 鞘层语义散落在 sheath/runtime 字段中，建议后续显式命名。 |
| `DipoleVect` | `VolField.cpp` 的 dipole field 计算 | `需显式命名桥接` | 偶极场已实现，但还不是清晰的 `DipoleVect` 入口。 |
| `Dirac` | 采样/离散极限逻辑 + piecewise helper | `已聚合承载` | Dirac 型近似以离散化和采样方式承载。 |
| `DistToPoint` | `Point3D::distanceTo`、`Vector3D` 距离函数 | `需显式命名桥接` | 距离到点的几何语义已存在，但 class 名义尚未显式。 |
| `DriftingMaxwellian3VFunction` | `ParticleSource.cpp`、`SurfaceDistributionFunction.h` | `需显式命名桥接` | 漂移 Maxwellian 在粒子/分布层有运行语义，值得补名称桥接。 |
| `Erf` | `MathUtils::errorFunction` | `已聚合承载` | 误差函数已存在。 |
| `ErfFromMath3` | `MathUtils::errorFunction` | `已聚合承载` | SPIS 对外部数学库的绑定已由本地 helper 吸收。 |
| `Exp` | `MathUtils.h` + `std::exp` | `已聚合承载` | 指数函数已存在。 |
| `ExpPower` | `MathUtils` + 幂律/指数组合公式 | `已聚合承载` | 指数幂组合可由现有公式直接表达。 |
| `FractionHandler` | 数值归一化、比例裁剪、`safeDivide` | `不建议逐类镜像` | 更像局部数值处理模式，不值得独立类。 |
| `FunctionOfpointDistance` | `Point3D` / `Vector3D` 距离计算 + 自定义函数 | `需显式命名桥接` | 距离函数组合已有能力，但缺少显式文档名义。 |
| `GammaFromMath3` | `MathUtils::gammaFunction` | `已聚合承载` | Gamma 函数已存在。 |
| `GenericOcTreeSplittingHeuristic` | `MeshPartitioning.cpp`、`util_octree_mesh` 证据链 | `需显式命名桥接` | Octree 切分启发式存在相近实现，但 class-level 仍需桥接索引。 |
| `GFunction` | 材料/场/分布层经验公式 | `需显式命名桥接` | 已有类似经验公式承载，但 `GFunction` 名义不够显式。 |
| `GradCosMultiD` | `MathUtils::gradient` + 多维余弦函数 | `已聚合承载` | 梯度与多维余弦已可组合得到。 |
| `GradPowerLaw` | `ParticleSource.cpp` 幂律拟合与梯度相关逻辑 | `需显式命名桥接` | 幂律导数语义存在，但应补显式 helper 名称。 |
| `Heaviside` | 阈值分段逻辑 + clamp/比较 | `已聚合承载` | 阶跃函数通常以分支或分段函数实现。 |
| `HistogramScalFunctionOfScal` | `NumericAggregation` + 采样统计逻辑 | `已聚合承载` | 直方图型标量函数可由统计层承担。 |
| `Identity` | lambda 恒等映射 | `已聚合承载` | 恒等函数不必单独建类。 |
| `IndexedFunctionof4Scal` | `std::function` + 索引访问 | `不建议逐类镜像` | 属于多参数函数接口壳层。 |
| `InvertableScalFunctionOfScal` | 反函数语义由局部求解器承担 | `不建议逐类镜像` | 可逆接口是抽象壳层，而非缺功能。 |
| `InvertibleScalFunctionOfScalFromFunctionOf4` | 多参函数包装成单参函数 | `不建议逐类镜像` | 典型 Java 包装器，C++ 无需直译。 |
| `IsDifferent` | `Constants.h` 容差比较、普通谓词 | `不建议逐类镜像` | 谓词接口不值得独立 class。 |
| `IsEqual` | `Constants.h` 的 `IS_EQUAL` | `不建议逐类镜像` | 已由容差比较宏显式承担。 |
| `Kappa1_1D` | `ParticleSource.cpp` 的 kappa 拟合与采样 | `需显式命名桥接` | Kappa 1D 语义已存在于粒子源公式中。 |
| `Kappa1_3D` | `ParticleSource.cpp` 的 kappa 拟合与热速度修正 | `需显式命名桥接` | 3D Kappa 语义也已有实现痕迹。 |
| `LimitedExp` | `MathUtils` + clamp 限幅指数 | `已聚合承载` | 限幅指数可由现有数学 helper 实现。 |
| `Log` | `MathUtils.h` + `std::log` | `已聚合承载` | 自然对数已存在。 |
| `log10` | `MathUtils.h` + `std::log10` | `已聚合承载` | 十进对数已存在。 |
| `LowerGammaIncFunction` | `MathUtils` 特殊函数扩展点 | `已聚合承载` | 归入特殊函数聚合承载。 |
| `MaskedIsDifferent` | 谓词 + mask 配置 | `不建议逐类镜像` | 属于条件包装壳层。 |
| `Max` | `std::max`、`MathUtils::clamp` | `已聚合承载` | 最大值操作已存在。 |
| `Maxwellian1D` | `BoundaryConditions.cpp`、`ParticleSource.cpp`、`SurfaceDistributionFunction.h` | `需显式命名桥接` | Maxwellian 1D 采样/分布已有明确运行时承载。 |
| `Min` | `std::min`、`MathUtils::clamp` | `已聚合承载` | 最小值操作已存在。 |
| `ModifInt` | 积分/整数修正辅助逻辑 | `不建议逐类镜像` | 更像内部修正包装器而非主要数学模型。 |
| `OMLFactor` | `SurfaceMaterialModel.cpp`、sheath/current 相关公式 | `需显式命名桥接` | OML 因子语义值得后续显式命名。 |
| `Or` | 条件组合由普通布尔逻辑承担 | `不建议逐类镜像` | 逻辑组合壳层无需单独类。 |
| `PanelEdge` | 几何边界/面片相关 helper | `需显式命名桥接` | 面板边界语义已分散在几何/mesh 层。 |
| `Polynomial` | 多项式公式由通用数学层直接表达 | `已聚合承载` | 无需单独 Java 风格 class。 |
| `PolynomialRatio` | 多项式比值由通用公式直接表达 | `已聚合承载` | 已可通过普通公式与安全除法表达。 |
| `Power` | `std::pow` | `已聚合承载` | 幂函数已存在。 |
| `PowerLaw` | `ParticleSource.cpp` 幂律拟合与源分布 | `需显式命名桥接` | 幂律在粒子源层已有显式 runtime 语义。 |
| `ProductOfTwoScalFunctionOfScal` | lambda 组合、普通乘积表达式 | `不建议逐类镜像` | 纯组合壳层。 |
| `ProductOfUICallableFunction` | CLI/config route + callable 组合 | `不建议逐类镜像` | 属于 SPIS UI 风格包装器。 |
| `ReciprocalFunction` | `safeDivide` + `1/x` 表达 | `已聚合承载` | 倒数函数已可直接表达。 |
| `ReverseUnitCheckable` | 单位检查接口的反向包装 | `不建议逐类镜像` | 当前项目不保留这类接口树。 |
| `ScalFuncOfScalFromTestOfInt` | 谓词转函数包装 | `不建议逐类镜像` | 属于适配壳层。 |
| `ScalFunctionOf2Scal` | `std::function<double(double,double)>` | `不建议逐类镜像` | 纯函数签名接口。 |
| `ScalFunctionOf3Scal` | `std::function` / lambda | `不建议逐类镜像` | 纯函数签名接口。 |
| `ScalFunctionOf4Scal` | `std::function` / lambda | `不建议逐类镜像` | 纯函数签名接口。 |
| `ScalFunctionOf5Scal` | `std::function` / lambda | `不建议逐类镜像` | 纯函数签名接口。 |
| `ScalFunctionOfNothing` | 零参 lambda | `不建议逐类镜像` | C++ 不必用 class 包装。 |
| `ScalFunctionOfOctree` | Octree 查询函数由 mesh/octree helper 承担 | `不建议逐类镜像` | 函数接口壳层，不是功能缺口。 |
| `ScalFunctionOfScal` | `std::function<double(double)>` | `不建议逐类镜像` | 纯函数签名接口。 |
| `ScalFunctionOfScalFromFunctionOf4` | 多参到单参包装 | `不建议逐类镜像` | 结构适配壳层。 |
| `ScalFunctionOfScalFromFunctionOf5` | 多参到单参包装 | `不建议逐类镜像` | 结构适配壳层。 |
| `ScalFunctionOfVect` | 向量到标量的 lambda / 几何 helper | `不建议逐类镜像` | 当前更适合函数对象而非 class 树。 |
| `SeparableFunctionOf2Scal` | 分离变量公式可由普通乘积表达 | `不建议逐类镜像` | 不值得单独类。 |
| `SeparableFunctionOf3Scal` | 分离变量公式可由普通乘积表达 | `不建议逐类镜像` | 不值得单独类。 |
| `Set` | 条件集 / 谓词集由普通容器与分支承担 | `不建议逐类镜像` | 结构壳层。 |
| `Sgn` | `Constants.h` 的 `SIGN` 宏 | `已聚合承载` | 符号函数已显式存在。 |
| `ShiftSet` | 区间/集合平移逻辑 | `不建议逐类镜像` | 更像结构包装器。 |
| `SimpleSurfSheath` | 材料/边界/transition 中的 sheath 长度与电容逻辑 | `需显式命名桥接` | 表面鞘层已有大量 runtime 字段，后续可补独立 helper 名称。 |
| `Sin` | `MathUtils.h` + `std::sin` | `已聚合承载` | 正弦函数已存在。 |
| `SmartTabulatedScalFunctionOf2Scal` | `MathUtils::interpolate2D` | `已聚合承载` | 智能二维表函数已由插值能力承担。 |
| `SolenoidVect` | `SurfaceFieldVolumeBridge.cpp` 的 `SolenoidBField` | `需显式命名桥接` | 螺线管场已有运行时 route。 |
| `Sqrt` | `std::sqrt` | `已聚合承载` | 平方根已存在。 |
| `StepWiseScalFunctionOfScal` | 分段函数、阈值路由、piecewise helper | `已聚合承载` | 阶梯函数可由分段逻辑承担。 |
| `TabulatedScalFunctionOf2Scal` | `MathUtils::interpolate2D` | `已聚合承载` | 二维表函数已显式存在。 |
| `TabulatedScalFunctionOf3Scal` | `MathUtils::interpolate3D` | `已聚合承载` | 三维表函数已显式存在。 |
| `TabulatedScalFunctionOfScal` | `MathUtils::interpolate1D`、`splineInterpolate` | `已聚合承载` | 一维表函数已显式存在。 |
| `Test` | 普通布尔测试接口 | `不建议逐类镜像` | Java 谓词接口壳层。 |
| `TestOfInt` | 整型谓词接口 | `不建议逐类镜像` | Java 谓词接口壳层。 |
| `ThinWireEnd` | 几何/边界近场语义 | `需显式命名桥接` | 细导线端点语义可能存在于几何/场近似中，值得单独桥接索引。 |
| `ThinWireLog` | 细导线对数项近似 | `需显式命名桥接` | 细导线相关对数近似值得后续显式命名。 |
| `TriLinearTabulatedScalFunctionOf3Scal` | `MathUtils::interpolate3D` | `已聚合承载` | 三线性插值已显式存在。 |
| `TruncatedLimitedExp` | 限幅+截断指数逻辑 | `已聚合承载` | 属于指数函数的聚合变体。 |
| `UICallableFunction` | CLI/config/loader 路由 | `不建议逐类镜像` | 当前项目以 CLI/JSON route 替代 SPIS UI-callable。 |
| `UnitCheckable` | 单位检查接口 | `不建议逐类镜像` | 当前以常量、命名和测试而非接口树承载。 |
| `UnitUtil` | `Constants.h` 单位换算宏、常量与稳定换算测试 | `已聚合承载` | 单位工具已直接存在。 |
| `VectFunctionOfVect` | 几何/场向量变换 helper | `需显式命名桥接` | 向量到向量函数广泛存在，但尚无统一显式命名桥接。 |
| `YZPermutation` | 坐标重排由向量/网格处理代码承担 | `需显式命名桥接` | 坐标置换语义有实现空间，但缺专门文档落点。 |

## 5. 桥接候选清单

当前最值得后续补“显式命名桥接”的函数族有四组：

1. 鞘层/边界物理
   - `DebyeSurfSheath`
   - `SimpleSurfSheath`
   - `OMLFactor`
   - 主要落点：`SurfaceMaterialModel`、`DensePlasmaSurfaceCharging`、`SurfaceCurrentLinearizationKernel`

2. 分布函数/粒子源
   - `DriftingMaxwellian3VFunction`
   - `Maxwellian1D`
   - `Kappa1_1D`
   - `Kappa1_3D`
   - `PowerLaw`
   - `GradPowerLaw`
   - 主要落点：`ParticleSource.cpp`、`SurfaceDistributionFunction.h/.cpp`

3. 场与几何向量
   - `BVect`
   - `DipoleVect`
   - `SolenoidVect`
   - `VectFunctionOfVect`
   - `YZPermutation`
   - 主要落点：`VolField.cpp`、`SurfaceFieldVolumeBridge.cpp`、`Point3D.h`、`Vector3D.h`

4. 几何近场 / 网格启发式
   - `DistToPoint`
   - `FunctionOfpointDistance`
   - `PanelEdge`
   - `ThinWireEnd`
   - `ThinWireLog`
   - `GenericOcTreeSplittingHeuristic`
   - 主要落点：`Point3D.h`、`MeshPartitioning.cpp`

## 6. 汇总结论

- `Util/Func` 当前不存在 family 级功能缺口，问题主要是 class-level 命名显式度不足。
- 大量 SPIS `Util/Func` class 实际属于 Java 函数接口和包装器；这类对象在当前 C++ 项目中不值得 1:1 镜像。
- 后续如果继续推进，优先级不是“补 98 个同名 class”，而是只对少数物理含义强、后续仍会被审计追问的函数族补显式桥接。

## 7. 证据索引

- 基础层
  - `Tools/Basic/include/MathUtils.h`
  - `Tools/Basic/include/Constants.h`
  - `Tools/Basic/include/NumericAggregation.h`
  - `Tools/Basic/include/ModelRegistry.h`
- domain helper
  - `Tools/Particle/include/SurfaceDistributionFunction.h`
  - `Tools/Particle/src/SurfaceDistributionFunction.cpp`
  - `Tools/Particle/src/ParticleSource.cpp`
  - `Tools/FieldSolver/src/VolField.cpp`
  - `Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`
  - `Tools/Mesh/src/MeshPartitioning.cpp`
  - `Tools/Material/include/SurfaceMaterialModel.h`
  - `Tools/Material/src/SurfaceMaterialModel.cpp`
  - `Tools/Boundary/src/BoundaryConditions.cpp`
- 测试/审计
  - `Tools/Basic/test/Basic_test.cpp`
  - `Toolkit/Surface Charging/docs/surface_spis_class_audit.json`
  - `Toolkit/Surface Charging/docs/surface_spis_class_audit_summary.md`
  - `Toolkit/Surface Charging/docs/surface_num_alignment_refactor_plan.md`
