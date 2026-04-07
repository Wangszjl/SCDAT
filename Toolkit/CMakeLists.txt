# ==============================================================================
# SCDAT — Toolkit 应用模块目录
#
# 策略：每个子目录自动扫描 .cpp/.h 并注册为独立静态库。
# 目录名含空格，CMake target 名使用下划线版本。
# 新增源文件后只需重新运行 CMake（无需手动修改本文件）。
# ==============================================================================

# 格式："target名称|目录名称"（以 | 分隔，支持目录名含空格）
set(TOOLKIT_MODULES
    "Internal_Charging|Internal Charging"
    "Plasma_Analysis|Plasma Analysis"
    "Radiation|Radiation"
    "Surface_Charging|Surface Charging"
    "Vacuum_Arc|Vacuum Arc"
)

foreach(entry IN LISTS TOOLKIT_MODULES)
    # ── 拆分 target 名 与 目录名 ──────────────────────────────────────────────
    string(REPLACE "|" ";" _parts "${entry}")
    list(GET _parts 0 _target)
    list(GET _parts 1 _dirname)

    set(_dir "${CMAKE_CURRENT_SOURCE_DIR}/${_dirname}")

    # ── 监听目录：新增 / 删除文件时触发重新配置 ──────────────────────────────
    set_property(
        DIRECTORY "${CMAKE_SOURCE_DIR}"
        APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_dir}"
    )

    # ── 自动收集源文件与头文件（include/src 分层） ────────────────────────────
    file(GLOB_RECURSE _sources
        "${_dir}/src/*.cpp"
        "${_dir}/src/*.cc"
        "${_dir}/src/*.cxx"
        "${_dir}/src/*.c"
    )
    file(GLOB_RECURSE _headers
        "${_dir}/include/*.h"
        "${_dir}/include/*.hpp"
        "${_dir}/include/*.hxx"
        "${_dir}/include/*.hh"
    )

    # 兼容旧布局：若仍存在根目录头源文件，继续纳入构建
    file(GLOB _legacy_sources
        "${_dir}/*.cpp"
        "${_dir}/*.cc"
        "${_dir}/*.cxx"
        "${_dir}/*.c"
    )
    file(GLOB _legacy_headers
        "${_dir}/*.h"
        "${_dir}/*.hpp"
        "${_dir}/*.hxx"
        "${_dir}/*.hh"
    )
    list(APPEND _sources ${_legacy_sources})
    list(APPEND _headers ${_legacy_headers})

    # ── 注册库目标 ────────────────────────────────────────────────────────────
    if(_sources)
        # 目录非空：创建真实静态库
        add_library(Toolkit_${_target} STATIC ${_sources} ${_headers})
        target_include_directories(Toolkit_${_target}
            PUBLIC
                $<BUILD_INTERFACE:${_dir}/include>
                $<BUILD_INTERFACE:${_dir}>
                $<INSTALL_INTERFACE:include/SCDAT/Toolkit/${_target}>
        )
        message(STATUS "[SCDAT] Toolkit::${_target} — STATIC (${_dir})")
    else()
        # 目录暂无源文件：创建 INTERFACE 占位库
        add_library(Toolkit_${_target} INTERFACE)
        target_include_directories(Toolkit_${_target}
            INTERFACE
                $<BUILD_INTERFACE:${_dir}/include>
                $<BUILD_INTERFACE:${_dir}>
                $<INSTALL_INTERFACE:include/SCDAT/Toolkit/${_target}>
        )
        message(STATUS "[SCDAT] Toolkit::${_target} — INTERFACE (empty, ${_dir})")
    endif()

    # 链接到对应的 Tools 模块（示例通用依赖，可按需定制）
    if(_sources)
        target_link_libraries(Toolkit_${_target}
            PUBLIC
                SCDAT::Tools::Basic
                SCDAT::Tools::Boundary
                SCDAT::Tools::Coupling
                SCDAT::Tools::Diagnostics
                SCDAT::Tools::FieldSolver
                SCDAT::Tools::Geometry
                SCDAT::Tools::Interactions
                SCDAT::Tools::Material
                SCDAT::Tools::Matrix
                SCDAT::Tools::Mesh
                SCDAT::Tools::Output
                SCDAT::Tools::Particle
                SCDAT::Tools::PICcore
                SCDAT::Tools::Solver
                SCDAT::Tools::Vector
        )
    endif()

    # ── 命名空间别名：target_link_libraries(... SCDAT::Toolkit::Surface_Charging)
    add_library(SCDAT::Toolkit::${_target} ALIAS Toolkit_${_target})

    if(_target STREQUAL "Vacuum_Arc" AND NOT TARGET SCDAT::Toolkit::Vaccum_Arc)
        add_library(SCDAT::Toolkit::Vaccum_Arc ALIAS Toolkit_${_target})
    endif()

    unset(_parts)
    unset(_target)
    unset(_dirname)
    unset(_dir)
    unset(_sources)
    unset(_headers)
    unset(_legacy_sources)
    unset(_legacy_headers)
endforeach()
