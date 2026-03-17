# Spacecraft Charging and Discharging Analysis Toolkit (SCDAT)

## 空间飞行器-等离子体相互作用套件（SCDAT）

本项目为计算空间飞行器与等离子体相互作用的工具包，提供了多种模型和方法来分析飞行器在空间环境中的充电和放电现象。该工具包旨在帮助研究人员和工程师更好地理解和预测飞行器在不同空间环境下的行为，从而提高飞行器的设计和运行安全性。

---

## 环境要求

| 工具 | 最低版本 | 备注 |
|------|---------|------|
| CMake | 3.20 | 建议使用 4.x |
| C++ 编译器 | GCC 13 / Clang 16 / MSVC 19.38 | 须支持 C++20 |
| Ninja | 任意 | 推荐，也可使用 MinGW Make |

> Windows 下推荐使用 [MinGW-W64 GCC 13](https://github.com/brechtsanders/winlibs_mingw) + Ninja，已通过验证（GCC 13.2.0 + CMake 4.0.2 + Ninja，Windows 11）。

---

## 编译步骤

### 1. 克隆仓库

```bash
git clone https://github.com/Wangszjl/SCDAT.git
cd SCDAT
```

### 2. 配置（以 Ninja + GCC 为例）

```bash
cmake -S . -B build -G "Ninja" -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Debug
```

如需 Release 构建，将 `Debug` 替换为 `Release`。

### 3. 编译

```bash
cmake --build build
```

编译成功后，可执行文件位于：
- `build/bin/SCDAT.exe` — 正式主程序
- `build/bin/SCDAT_PICcore.exe` — PICcore 独立验证入口（来自 `Main/main_PICcore.cpp`）

### 4. 运行

```bash
./build/bin/SCDAT.exe
```

---

## 项目结构

```
SCDAT/
├── CMakeLists.txt          # 根配置（C++20，链接三个子目录）
├── Main/
│   ├── CMakeLists.txt      # 主可执行目标，支持 main_*.cpp 自动注册
│   ├── main.cpp            # 正式主程序入口 → SCDAT
│   └── main_PICcore.cpp    # PICcore 测试入口 → SCDAT_PICcore
├── Tools/
│   ├── CMakeLists.txt      # 自动 GLOB，每个子目录注册为独立静态库
│   ├── Basic/
│   ├── Interactions/
│   ├── Matrix/
│   ├── Mesh/
│   ├── Parallel/
│   ├── PICcore/
│   ├── Scripts/
│   ├── UI/
│   ├── Vector/
│   └── Visualization/
└── Toolkit/
    ├── CMakeLists.txt      # 自动 GLOB，支持目录名含空格
    ├── Internal Charging/
    ├── Plasma Analysis/
    ├── Surface Charging/
    └── Vaccum Arc/
```

---

## 新增源文件的工作流

本项目的 CMakeLists 使用 **`file(GLOB_RECURSE)`** + **`CMAKE_CONFIGURE_DEPENDS`** 自动追踪文件变化：

- 在 `Tools/<模块>/` 或 `Toolkit/<模块>/` 下新建 `.cpp` / `.h` 文件
- 在 `Main/` 下新建 `main_<名称>.cpp`（含独立 `main` 函数）即可自动生成新的可执行目标 `SCDAT_<名称>`
- 直接运行 `cmake --build build` — CMake 会检测到目录内容变化并自动重新配置，**无需手动修改任何 CMakeLists.txt**

