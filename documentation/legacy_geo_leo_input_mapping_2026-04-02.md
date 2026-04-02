# Legacy GEO/LEO 输入字段映射（执行基线）

本文件把 legacy `.dat` 输入到当前 `SurfaceChargingConfig` 的映射关系显式化，用于回放验证与后续重构对齐。

## 1. 映射代码基准

- 入口函数：`translateLegacyBenchmarkConfig(...)`
- 文件：`Toolkit/Surface Charging/src/LegacyBenchmarkSupport.cpp`
- 关键行为：
  - 强制启用 legacy 回放路由 `LegacyBenchmark`
  - 关闭 live PIC/mcc 和 body-patch 电路
  - 使用 `raw_values` 按索引读取 GEO/LEO 输入参数

## 2. GEO 输入索引（`ref/GEO/input_GEO_Charging.dat`）

`raw_values` 顺序映射（对应 `is_geo == true` 分支）：

| 索引 | legacy 含义 | 目标配置字段 / 逻辑 |
|---:|---|---|
| 0 | `jsun` | `patch_incidence_angle_deg = (jsun == 0 ? 90 : theta_deg)` |
| 1 | `jenv` | 环境 ID，匹配 `LegacyEnvironmentRecord.case_id` |
| 2 | `theta_deg` | 入射角（与 `jsun` 共同决定 patch 入射角） |
| 3 | `jmb` | 结构材料 ID，优先在结构材料表匹配 |
| 4 | `jsee` | SEE 模型映射：`1->Whipple,2->Sims,3->Katz` |
| 5 | `jp` | patch 类型：`1` 介质 patch，`2` 金属 patch |
| 6 | `jmp` | patch 材料 ID（在 `jp` 对应材料表中查找） |
| 7 | `thickness_m` | `dielectric_thickness_m`（最小 1e-9 m） |
| 8 | `epsilon_r` | 当 `jp==2` 时覆盖层介电常数 |
| 9 | `sigma` | 当 `jp==2` 时覆盖层电导率 |

## 3. LEO 输入索引（`ref/LEO/input_LEO_Charging.dat`）

`raw_values` 顺序映射（对应 `is_geo == false` 分支）：

| 索引 | legacy 含义 | 目标配置字段 / 逻辑 |
|---:|---|---|
| 0 | `jsun` | `patch_incidence_angle_deg` 规则同 GEO |
| 1 | `jram` | RAM/WAKE 逻辑开关 |
| 2 | `altitude_km` | 用于轨道速度估算 `legacyLeoOrbitalSpeedMPerS` |
| 3 | `flow_angle_deg` | `patch_flow_angle_deg`，并参与 `flow_alignment_cosine` |
| 4 | `jenv` | 环境 ID |
| 5 | `theta_deg` | 入射角 |
| 6 | `jmb` | 结构材料 ID |
| 7 | `jsee` | SEE 模型 |
| 8 | `jp` | patch 类型 |
| 9 | `jmp` | patch 材料 ID |
| 10 | `thickness_m` | `dielectric_thickness_m` |
| 11 | `epsilon_r` | 当 `jp==2` 时覆盖层介电常数 |
| 12 | `sigma` | 当 `jp==2` 时覆盖层电导率 |

`jram` 规则：

- `jram == 1`：无定向流（`bulk_flow_velocity_m_per_s = 0`, `flow_alignment_cosine = 0`）
- `jram == 2`：RAM 方向，`flow_alignment_cosine = max(0, cos(flow_angle))`
- `jram == 3`：WAKE 方向，`flow_alignment_cosine = -abs(cos(flow_angle))`

## 4. 环境文件热参数映射（GEO/LEO 通用）

环境解析来源：`parseLegacyEnvironmentFile(...)` 与 `applyLegacyEnvironmentToConfig(...)`。

`thermal_values` 长度至少 15 时，索引映射如下：

| 索引 | 含义 | 目标 |
|---:|---|---|
| 0,2,4 | `Ne1/2/3 (1/cm^3)` | 电子密度分量，转换到 `m^-3`（乘 1e6） |
| 1,3,5 | `Te1/2/3 (eV)` | 电子温度分量 |
| 6,9,12 | `Ni1/2/3 (1/cm^3)` | 离子密度分量，转换到 `m^-3` |
| 7,10,13 | `Ti1/2/3 (eV)` | 离子温度分量 |
| 8,11,14 | `Mi1/2/3 (amu)` | 离子质量分量 |

聚合策略：

- 总密度：各分量直接求和
- 温度/离子质量：按分量密度加权平均
- 电子/离子能谱：从环境文件 `Electron energy` / `Ion energy` 与后续 `Flux(...)` 行装载到 TABULATED 频谱

## 5. 材料选择与电容

- 结构材料选择优先级：`structure` -> `metal_patch` -> `dielectric_patch` -> fallback
- patch 材料来源：
  - `jp==1`：`dielectric_patch_materials[jmp]`
  - `jp==2`：`metal_patch_materials[jmp]`
- patch 电容：
  - `jp==1`：`epsilon0 * epsilon_r / thickness`
  - `jp==2`：沿用输入或 fallback 的 `capacitance_per_area_f_per_m2`

## 6. 维护约定

- 本文件与 `validate_legacy_replay_csv.py` 共同构成 P0 回放基线。
- 若修改 `raw_values` 解释或字段命名，需同步更新：
  1. 本映射文档
  2. 对应回放校验脚本阈值或列映射
  3. CTest 用例参数
