# ==============================================================================
# SCDAT — Main 可执行入口
#
# 自动扫描本目录下所有 .cpp 源文件，生成主可执行文件 SCDAT。
# 新增源文件后只需重新运行 CMake。
# ==============================================================================

# ── 监听目录变化 ──────────────────────────────────────────────────────────────
set_property(
    DIRECTORY "${CMAKE_SOURCE_DIR}"
    APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}"
)

# ── 自动收集主程序源文件 ──────────────────────────────────────────────────────
file(GLOB_RECURSE MAIN_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

# ── 可执行目标 ────────────────────────────────────────────────────────────────
add_executable(SCDAT ${MAIN_SOURCES})

# ── 链接所有 Tools 模块 ───────────────────────────────────────────────────────
target_link_libraries(SCDAT
    PRIVATE
        # Tools 基础模块
        SCDAT::Tools::Basic
        SCDAT::Tools::Interactions
        SCDAT::Tools::Matrix
        SCDAT::Tools::Mesh
        SCDAT::Tools::Parallel
        SCDAT::Tools::PICcore
        SCDAT::Tools::Scripts
        SCDAT::Tools::UI
        SCDAT::Tools::Vector
        SCDAT::Tools::Visualization
        # Toolkit 应用模块
        SCDAT::Toolkit::Internal_Charging
        SCDAT::Toolkit::Plasma_Analysis
        SCDAT::Toolkit::Surface_Charging
        SCDAT::Toolkit::Vaccum_Arc
)

target_include_directories(SCDAT
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

message(STATUS "[SCDAT] Main executable target registered")
