# ==============================================================================
# Spacecraft Charging and Discharging Analysis Toolkit (SCDAT)
# 空间飞行器充放电分析工具包
#
# 根 CMakeLists.txt — 项目顶层配置
# ==============================================================================

cmake_minimum_required(VERSION 3.20)

project(SCDAT
    VERSION     0.1.0
    DESCRIPTION "Spacecraft Charging and Discharging Analysis Toolkit"
    LANGUAGES   CXX
)

include(CTest)

# ------------------------------------------------------------------------------
# 全局 C++ 标准
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------------------------------------------------------
# 输出目录
# ------------------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ------------------------------------------------------------------------------
# 全局头文件搜索路径（让所有子模块互相可见）
# ------------------------------------------------------------------------------
include_directories(
    ${CMAKE_SOURCE_DIR}/Tools
    ${CMAKE_SOURCE_DIR}/Toolkit
    ${CMAKE_SOURCE_DIR}/Main
)

# ------------------------------------------------------------------------------
# 子目录
# ------------------------------------------------------------------------------
add_subdirectory(Tools)
add_subdirectory(Toolkit)
add_subdirectory(Main)

if(BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)

    add_test(
        NAME ToolsLayoutCheck
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_tools_layout.py
            --tools-root ${CMAKE_SOURCE_DIR}/Tools
            --strict
    )

    add_test(
        NAME LegacyBenchmarkAcceptanceGate
        COMMAND
            ${Python3_EXECUTABLE}
            "${CMAKE_SOURCE_DIR}/Toolkit/Surface Charging/scripts/benchmark_regression_gate.py"
            --glob "${CMAKE_SOURCE_DIR}/results/*.benchmark.txt"
            --allow-empty
    )

    add_test(
        NAME LegacyReplayCsvConsistencyGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/validate_legacy_replay_csv.py
            --csv
                ${CMAKE_SOURCE_DIR}/results/tmp_geo_legacy.csv
                ${CMAKE_SOURCE_DIR}/results/tmp_leo_ram_legacy.csv
                ${CMAKE_SOURCE_DIR}/results/tmp_leo_wake_legacy.csv
            --strict-schema
            --max-abs-delta 1e-12
    )

    add_test(
        NAME Kapton2000sRegressionGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_kapton_2000s_regression.py
            --advanced-csv ${CMAKE_SOURCE_DIR}/results/surface_geo_ecss_kapton_pic_circuit.csv
            --reference-csv ${CMAKE_SOURCE_DIR}/results/surface_geo_ecss_kapton_ref.csv
            --target-time-s 2000
            --max-time-offset-s 1e-6
            --max-abs-surface-delta-v 3000
            --max-equilibrium-gap-v 1e-6
            --require-advanced-more-negative
    )

    add_subdirectory(Tools/Geometry/test)
    add_subdirectory(Tools/Coupling/test)
    add_subdirectory(Tools/Diagnostics/test)
    add_subdirectory(Tools/FieldSolver/test)
    add_subdirectory(Tools/Material/test)
    add_subdirectory(Tools/Output/test)
    add_subdirectory(Tools/PICcore/test)
    add_subdirectory(Tools/Particle/test)
    add_subdirectory("Toolkit/Plasma Analysis/test")
    add_subdirectory("Toolkit/Surface Charging/test")
    add_subdirectory("Toolkit/Internal Charging/test")
    add_subdirectory("Toolkit/Vacuum Arc/test")
endif()
