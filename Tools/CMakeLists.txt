# ==============================================================================
# SCDAT — Tools 子库目录
#
# 策略：每个子目录自动扫描 .cpp/.h 并注册为独立静态库。
# 新增源文件后只需重新运行 CMake（无需手动修改本文件），
# 原理：set_property CMAKE_CONFIGURE_DEPENDS 监听目录变化，
#       文件增减时 CMake 自动重新配置。
# ==============================================================================

# 所有 Tools 模块目录（名称 = CMake target 后缀 = 子目录名）
set(TOOLS_MODULES
    Basic
    Boundary
    Coupling
    Diagnostics
    FieldSolver
    Geometry
    Interactions
    Material
    Mesh
    Output
    Particle
    Parallel
    PICcore
    Scripts
    Solver
    UI
    Vector
    Visualization
)

foreach(module IN LISTS TOOLS_MODULES)
    set(_dir "${CMAKE_CURRENT_SOURCE_DIR}/${module}")

    # ── 监听目录：新增 / 删除文件时触发重新配置 ──────────────────────────────
    set_property(
        DIRECTORY "${CMAKE_SOURCE_DIR}"
        APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_dir}"
    )

    # ── 自动收集源文件与头文件 ────────────────────────────────────────────────
    file(GLOB_RECURSE _sources
        "${_dir}/*.cpp"
        "${_dir}/*.cc"
        "${_dir}/*.cxx"
    )
    file(GLOB_RECURSE _headers
        "${_dir}/*.h"
        "${_dir}/*.hpp"
        "${_dir}/*.hxx"
    )

    # 测试代码由独立测试目标管理，避免混入模块库目标。
    list(FILTER _sources EXCLUDE REGEX ".*/test/.*")
    list(FILTER _headers EXCLUDE REGEX ".*/test/.*")

    # Mesh模块系统性重构：仅保留三类聚合编译入口文件
    if(module STREQUAL "Mesh")
        list(FILTER _sources INCLUDE REGEX ".*[\\/](MeshAlgorithms|MeshPartitioning|MeshParsing)\\.cpp$")
    endif()

    # ── 注册库目标 ────────────────────────────────────────────────────────────
    if(_sources)
        # 目录非空：创建真实静态库
        add_library(Tools_${module} STATIC ${_sources} ${_headers})
        target_include_directories(Tools_${module}
            PUBLIC
                $<BUILD_INTERFACE:${_dir}/include>
                $<BUILD_INTERFACE:${_dir}>
                $<INSTALL_INTERFACE:include/SCDAT/Tools/${module}>
        )
        message(STATUS "[SCDAT] Tools::${module} — STATIC (${_dir})")
    else()
        # 目录暂无源文件：创建 INTERFACE 占位库，头文件仍可用
        add_library(Tools_${module} INTERFACE)
        target_include_directories(Tools_${module}
            INTERFACE
                $<BUILD_INTERFACE:${_dir}/include>
                $<BUILD_INTERFACE:${_dir}>
                $<INSTALL_INTERFACE:include/SCDAT/Tools/${module}>
        )
        message(STATUS "[SCDAT] Tools::${module} — INTERFACE (empty, ${_dir})")
    endif()

    # ── 命名空间别名：target_link_libraries(... SCDAT::Tools::Basic) ─────────
    add_library(SCDAT::Tools::${module} ALIAS Tools_${module})

    unset(_dir)
    unset(_sources)
    unset(_headers)
endforeach()

# 兼容旧的链接名：历史上Main模块链接到SCDAT::Tools::Matrix。
if(TARGET Tools_Geometry AND NOT TARGET SCDAT::Tools::Matrix)
    add_library(SCDAT::Tools::Matrix ALIAS Tools_Geometry)
endif()

if(TARGET Tools_Particle AND TARGET Tools_Geometry)
    target_link_libraries(Tools_Particle PUBLIC Tools_Geometry)
endif()

if(TARGET Tools_Solver)
    target_link_libraries(Tools_Solver
        PUBLIC
            Tools_Basic
            Tools_Geometry
            Tools_Mesh
    )

    if(ENABLE_SUPERLU)
        find_package(SuperLU QUIET)

        set(_superlu_linked FALSE)

        if(TARGET SuperLU::SuperLU)
            target_link_libraries(Tools_Solver PUBLIC SuperLU::SuperLU)
            set(_superlu_linked TRUE)
        elseif(TARGET SuperLU::superlu)
            target_link_libraries(Tools_Solver PUBLIC SuperLU::superlu)
            set(_superlu_linked TRUE)
        endif()

        if(NOT _superlu_linked)
            set(_superlu_libs "")
            set(_superlu_includes "")

            if(DEFINED SuperLU_LIBRARIES)
                set(_superlu_libs ${SuperLU_LIBRARIES})
            elseif(DEFINED SUPERLU_LIBRARIES)
                set(_superlu_libs ${SUPERLU_LIBRARIES})
            endif()

            if(DEFINED SuperLU_INCLUDE_DIRS)
                set(_superlu_includes ${SuperLU_INCLUDE_DIRS})
            elseif(DEFINED SUPERLU_INCLUDE_DIRS)
                set(_superlu_includes ${SUPERLU_INCLUDE_DIRS})
            endif()

            if(NOT _superlu_libs)
                find_path(_superlu_header_dir slu_ddefs.h)
                find_library(_superlu_library NAMES superlu)
                if(_superlu_header_dir AND _superlu_library)
                    set(_superlu_includes ${_superlu_header_dir})
                    set(_superlu_libs ${_superlu_library})
                endif()
            endif()

            if(NOT _superlu_libs)
                message(FATAL_ERROR "ENABLE_SUPERLU=ON but SuperLU was not found. Provide CMAKE_PREFIX_PATH/SuperLU_DIR or install development headers/libs.")
            endif()

            if(_superlu_includes)
                target_include_directories(Tools_Solver PUBLIC ${_superlu_includes})
            endif()
            target_link_libraries(Tools_Solver PUBLIC ${_superlu_libs})
        endif()

        target_compile_definitions(Tools_Solver PUBLIC HAVE_SUPERLU)

        message(STATUS "[SCDAT] Tools::Solver - SuperLU ENABLED")
    else()
        message(STATUS "[SCDAT] Tools::Solver - SuperLU DISABLED")
    endif()
endif()

if(TARGET Tools_Material)
    target_link_libraries(Tools_Material
        PUBLIC
            Tools_Basic
            Tools_Geometry
            Tools_Mesh
    )
endif()

if(TARGET Tools_Diagnostics)
    target_link_libraries(Tools_Diagnostics
        PUBLIC
            Tools_Basic
            Tools_Geometry
            Tools_Particle
    )
endif()

if(TARGET Tools_Output)
    target_link_libraries(Tools_Output
        PUBLIC
            Tools_Basic
    )
endif()

if(TARGET Tools_Coupling)
    target_link_libraries(Tools_Coupling
        PUBLIC
            Tools_Basic
    )
endif()

if(TARGET Tools_FieldSolver)
    target_link_libraries(Tools_FieldSolver
        PUBLIC
            Tools_Basic
            Tools_Geometry
            Tools_Material
            Tools_Mesh
            Tools_Solver
    )
endif()

if(TARGET Tools_Particle)
    target_link_libraries(Tools_Particle
        PUBLIC
            Tools_Basic
            Tools_Geometry
            Tools_Mesh
    )
endif()

if(TARGET Tools_Boundary)
    target_link_libraries(Tools_Boundary
        PUBLIC
            Tools_Geometry
            Tools_Particle
    )
endif()

if(TARGET Tools_Parallel)
    target_link_libraries(Tools_Parallel
        PUBLIC
            Tools_Solver
    )
endif()

if(TARGET Tools_PICcore)
    target_link_libraries(Tools_PICcore
        PUBLIC
            Tools_Basic
            Tools_Boundary
            Tools_FieldSolver
            Tools_Geometry
            Tools_Mesh
            Tools_Particle
            Tools_Solver
    )
endif()
