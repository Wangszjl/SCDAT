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
    add_subdirectory(Tools/Geometry/test)
    add_subdirectory(Tools/Particle/test)
endif()
