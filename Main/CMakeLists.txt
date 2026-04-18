# ==============================================================================
# SCDAT — Main 可执行入口
#
# 命名约定：
#   main.cpp          → 可执行目标 SCDAT（正式主程序）
#   main_<name>.cpp   → 可执行目标 SCDAT_<name>（独立入口 / 演示 / 测试）
#
# 新增 main_xxx.cpp 后只需重新运行 CMake，无需修改本文件。
# ==============================================================================

# ── 监听目录变化 ──────────────────────────────────────────────────────────────
set_property(
    DIRECTORY "${CMAKE_SOURCE_DIR}"
    APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}"
)

# ── 公共链接依赖（所有入口共享）──────────────────────────────────────────────
set(_SCDAT_LINK_LIBS
    SCDAT::Tools::Basic
    SCDAT::Tools::Coupling
    SCDAT::Tools::Diagnostics
    SCDAT::Tools::Interactions
    SCDAT::Tools::Material
    SCDAT::Tools::Matrix
    SCDAT::Tools::Mesh
    SCDAT::Tools::Output
    SCDAT::Tools::Parallel
    SCDAT::Tools::PICcore
    SCDAT::Tools::Scripts
    SCDAT::Tools::UI
    SCDAT::Tools::Vector
    SCDAT::Tools::Visualization
    SCDAT::Toolkit::Internal_Charging
    SCDAT::Toolkit::Plasma_Analysis
    SCDAT::Toolkit::Radiation
    SCDAT::Toolkit::Surface_Charging
    SCDAT::Toolkit::Vacuum_Arc
)

# ── 正式主程序：main.cpp → 可执行 SCDAT ─────────────────────────────────────
add_executable(SCDAT
    "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SurfaceScenarioLoader.cpp"
)
scdat_harden_build_target(SCDAT)

target_link_libraries(SCDAT PRIVATE ${_SCDAT_LINK_LIBS})
target_include_directories(SCDAT PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# ── 附加入口：main_<name>.cpp → 可执行 SCDAT_<name> ─────────────────────────
file(GLOB _extra_mains "${CMAKE_CURRENT_SOURCE_DIR}/main_*.cpp")

foreach(_src IN LISTS _extra_mains)
    get_filename_component(_stem "${_src}" NAME_WE)   # e.g. main_PICcore
    string(REGEX REPLACE "^main_" "" _suffix "${_stem}")  # e.g. PICcore
    set(_tgt "SCDAT_${_suffix}")

    add_executable(${_tgt} "${_src}")
    scdat_harden_build_target(${_tgt})
    target_link_libraries(${_tgt} PRIVATE ${_SCDAT_LINK_LIBS})
    if(_suffix STREQUAL "mesh_export")
        # MinGW static link order is sensitive; append Geometry after Mesh consumers.
        target_link_libraries(${_tgt} PRIVATE SCDAT::Tools::Geometry)
    endif()
    target_include_directories(${_tgt} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

    message(STATUS "[SCDAT] Extra entry: ${_tgt} <- ${_src}")
    unset(_tgt)
    unset(_suffix)
    unset(_stem)
endforeach()

message(STATUS "[SCDAT] Main executable targets registered")
