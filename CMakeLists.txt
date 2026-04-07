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

option(ENABLE_SUPERLU "Enable SuperLU-backed sparse direct solver path" OFF)

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

    add_test(
        NAME SurfaceLeoUnitLintGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/lint_surface_leo_units.py
            --config ${CMAKE_SOURCE_DIR}/scripts/run/surface_config_full_example.json
            --strict
            --report-json ${CMAKE_BINARY_DIR}/surface_leo_unit_lint_report.json
    )

    add_test(
        NAME UnifiedConfigSchemaV1Gate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_config_schema_v1_gate.py
            --project-root ${CMAKE_SOURCE_DIR}
            --schema ${CMAKE_SOURCE_DIR}/documentation/contracts/config_schema_v1.json
            --scdat-exe $<TARGET_FILE:SCDAT>
            --output-root ${CMAKE_BINARY_DIR}/config_schema_v1_gate
            --report-json ${CMAKE_BINARY_DIR}/config_schema_v1_gate.json
            --report-md ${CMAKE_BINARY_DIR}/config_schema_v1_gate.md
    )

    add_test(
        NAME UnifiedOutputContractV1Gate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_output_contract_v1_gate.py
            --project-root ${CMAKE_SOURCE_DIR}
            --scdat-exe $<TARGET_FILE:SCDAT>
            --output-root ${CMAKE_BINARY_DIR}/output_contract_v1_gate
            --report-json ${CMAKE_BINARY_DIR}/output_contract_v1_gate.json
            --report-md ${CMAKE_BINARY_DIR}/output_contract_v1_gate.md
    )

    add_test(
        NAME FourLayerGateSystemGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_four_layer_gate_system.py
            --project-root ${CMAKE_SOURCE_DIR}
            --build-dir ${CMAKE_BINARY_DIR}
            --config Debug
            --report-json ${CMAKE_BINARY_DIR}/four_layer_gate_report.json
            --report-md ${CMAKE_BINARY_DIR}/four_layer_gate_report.md
    )

    add_test(
        NAME DriftMonitoringDashboardGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_drift_monitoring_dashboard_gate.py
            --project-root ${CMAKE_SOURCE_DIR}
            --build-dir ${CMAKE_BINARY_DIR}
            --config Debug
            --auto-run-missing
            --schema-report ${CMAKE_BINARY_DIR}/config_schema_v1_gate.json
            --output-contract-report ${CMAKE_BINARY_DIR}/output_contract_v1_gate.json
            --unit-lint-report ${CMAKE_BINARY_DIR}/surface_leo_unit_lint_report.json
            --four-layer-report ${CMAKE_BINARY_DIR}/four_layer_gate_report.json
            --history-json ${CMAKE_BINARY_DIR}/drift_monitor_dashboard_history.json
            --report-json ${CMAKE_BINARY_DIR}/drift_monitor_dashboard.json
            --report-md ${CMAKE_BINARY_DIR}/drift_monitor_dashboard.md
    )
    set_tests_properties(
        DriftMonitoringDashboardGate
        PROPERTIES DEPENDS "SurfaceLeoUnitLintGate;UnifiedConfigSchemaV1Gate;UnifiedOutputContractV1Gate;FourLayerGateSystemGate"
    )

    add_test(
        NAME ReproducibleBuildDependencyLockGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_reproducible_build_dependency_lock_gate.py
            --project-root ${CMAKE_SOURCE_DIR}
            --lock-json ${CMAKE_SOURCE_DIR}/scripts/reproducible/toolchain_lock_v1.json
            --report-json ${CMAKE_BINARY_DIR}/reproducible_build_lock_gate.json
            --report-md ${CMAKE_BINARY_DIR}/reproducible_build_lock_gate.md
    )

    add_test(
        NAME PrecisionPerformanceParetoGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_precision_performance_pareto_gate.py
            --project-root ${CMAKE_SOURCE_DIR}
            --template ${CMAKE_SOURCE_DIR}/documentation/templates/precision_performance_pareto_template_v1.md
            --four-layer-report ${CMAKE_BINARY_DIR}/four_layer_gate_report.json
            --drift-report ${CMAKE_BINARY_DIR}/drift_monitor_dashboard.json
            --output-root ${CMAKE_BINARY_DIR}/pareto_reports
            --report-json ${CMAKE_BINARY_DIR}/pareto_report_workflow_gate.json
            --report-md ${CMAKE_BINARY_DIR}/pareto_report_workflow_gate.md
    )
    set_tests_properties(
        PrecisionPerformanceParetoGate
        PROPERTIES DEPENDS "FourLayerGateSystemGate;DriftMonitoringDashboardGate"
    )

    add_test(
        NAME DocumentationClosurePolicyGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_documentation_closure_gate.py
            --project-root ${CMAKE_SOURCE_DIR}
            --template-dir ${CMAKE_SOURCE_DIR}/documentation/templates
            --runlog-glob "documentation/*_stage_runlog_*.md"
            --min-runlog-count 10
            --report-json ${CMAKE_BINARY_DIR}/documentation_closure_gate.json
            --report-md ${CMAKE_BINARY_DIR}/documentation_closure_gate.md
    )

    add_test(
        NAME ReleaseCandidateSignatureGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_release_candidate_signature_gate.py
            --project-root ${CMAKE_SOURCE_DIR}
            --build-dir ${CMAKE_BINARY_DIR}
            --candidate-label m6-rc-2026-04-07
            --documentation-closure-report ${CMAKE_BINARY_DIR}/documentation_closure_gate.json
            --pareto-report ${CMAKE_BINARY_DIR}/pareto_report_workflow_gate.json
            --drift-report ${CMAKE_BINARY_DIR}/drift_monitor_dashboard.json
            --reproducible-lock-report ${CMAKE_BINARY_DIR}/reproducible_build_lock_gate.json
            --release-matrix ${CMAKE_SOURCE_DIR}/documentation/release_acceptance_matrix_2026-04-07.md
            --execution-plan ${CMAKE_SOURCE_DIR}/documentation/full_alignment_execution_plan_2026-04-03.md
            --legacy-warning-budget 0
            --signature-json ${CMAKE_BINARY_DIR}/release_candidate_signature.json
            --warning-inventory-json ${CMAKE_BINARY_DIR}/documentation_legacy_warning_inventory.json
            --warning-inventory-md ${CMAKE_BINARY_DIR}/documentation_legacy_warning_inventory.md
            --report-json ${CMAKE_BINARY_DIR}/release_candidate_signature_gate.json
            --report-md ${CMAKE_BINARY_DIR}/release_candidate_signature_gate.md
    )
    set_tests_properties(
        ReleaseCandidateSignatureGate
        PROPERTIES DEPENDS "DocumentationClosurePolicyGate;PrecisionPerformanceParetoGate;ReproducibleBuildDependencyLockGate;DriftMonitoringDashboardGate"
    )

    add_test(
        NAME VacuumArcCaseMatrixDryRunGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/run/run_vacuum_arc_case_matrix.py
            --matrix-json ${CMAKE_SOURCE_DIR}/scripts/run/arc_case_matrix_example.json
            --project-root ${CMAKE_SOURCE_DIR}
            --dry-run
            --output-root ${CMAKE_BINARY_DIR}/arc_repro_gate
    )

    add_test(
        NAME VacuumArcArcPicCompareSelfGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_vacuum_arc_matrix_gate.py
            --project-root ${CMAKE_SOURCE_DIR}
            --scdat-exe $<TARGET_FILE:SCDAT>
            --matrix-json ${CMAKE_SOURCE_DIR}/scripts/run/arc_case_matrix_example.json
            --output-root ${CMAKE_BINARY_DIR}/arc_compare_matrix_gate
            --columns discharge_current_a current_density_a_per_m2 channel_conductivity_s_per_m cathode_temperature_k anode_temperature_k surface_potential_v surface_charge_density_c_per_m2 stability_substeps_used stability_rollbacks
            --breakdown-column discharge_current_a
            --breakdown-threshold 1e-9
            --max-abs-delta 1e-2
            --max-rel-rmse 1e-3
            --max-rel-peak-delta 1e-3
            --max-rel-final-delta 1e-3
            --max-rel-integral-delta 1e-3
            --max-rel-steady-mean-delta 1e-3
            --max-breakdown-time-delta-s 1e-8
            --summary-json ${CMAKE_BINARY_DIR}/arc_compare_matrix_gate.json
            --summary-md ${CMAKE_BINARY_DIR}/arc_compare_matrix_gate.md
    )

    add_test(
        NAME InternalRadiationCouplingRegressionGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_internal_radiation_coupling_gate.py
            --scdat-exe $<TARGET_FILE:SCDAT>
            --internal-preset geo_electron_belt
            --radiation-preset geo_electron_belt_dose
            --csv ${CMAKE_BINARY_DIR}/internal_radiation_gate.csv
            --min-rows 20
            --max-ratio-cv 1e-4
            --max-rho-increase-count 0
            --require-monotonic-electric-field
    )

    add_test(
        NAME InternalRadiationComponentTrendGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_internal_radiation_component_gate.py
            --scdat-exe $<TARGET_FILE:SCDAT>
            --output-dir ${CMAKE_BINARY_DIR}/internal_radiation_component_gate
            --internal-preset geo_electron_belt
            --coupling-off-radiation-preset coupling_zero_flux_baseline
            --coupling-on-radiation-preset geo_electron_belt_dose
            --max-baseline-electric-delta-v-per-m 1e-10
            --max-baseline-charge-delta-c-per-m3 1e-14
            --min-final-electric-gain-v-per-m 1e-6
            --min-integrated-charge-gain-c-per-m2 1e-12
    )

    add_test(
        NAME RadiationGeant4AlignmentGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_radiation_geant4_alignment_gate.py
            --project-root ${CMAKE_SOURCE_DIR}
            --scdat-exe $<TARGET_FILE:SCDAT>
            --matrix-json ${CMAKE_SOURCE_DIR}/scripts/run/radiation_geant4_alignment_matrix.json
            --baseline-json ${CMAKE_SOURCE_DIR}/results/baselines/radiation/radiation_geant4_alignment_baseline_2026-04-04.json
            --output-root ${CMAKE_BINARY_DIR}/radiation_geant4_alignment_gate
            --summary-json ${CMAKE_BINARY_DIR}/radiation_geant4_alignment_gate.json
            --summary-md ${CMAKE_BINARY_DIR}/radiation_geant4_alignment_gate.md
    )

    add_test(
        NAME RadiationUqBudgetGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/analyze_radiation_uq_budget.py
            --project-root ${CMAKE_SOURCE_DIR}
            --scdat-exe $<TARGET_FILE:SCDAT>
            --matrix-json ${CMAKE_SOURCE_DIR}/scripts/run/radiation_geant4_alignment_matrix.json
            --output-root ${CMAKE_BINARY_DIR}/radiation_uq_budget
            --report-json ${CMAKE_BINARY_DIR}/radiation_uq_budget.json
            --report-md ${CMAKE_BINARY_DIR}/radiation_uq_budget.md
    )

    add_test(
        NAME VacuumArcCrossPlatformConsistencyGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_vacuum_arc_cross_platform_consistency_gate.py
            --synthesize-platform-summaries
            --project-root ${CMAKE_SOURCE_DIR}
            --scdat-exe $<TARGET_FILE:SCDAT>
            --matrix-json ${CMAKE_SOURCE_DIR}/scripts/run/arc_case_matrix_example.json
            --synth-output-root ${CMAKE_BINARY_DIR}/arc_cross_platform_consistency_gate
            --platform-a-label windows
            --platform-b-label linux
            --max-rel-rmse-delta 1e-3
            --max-rel-peak-delta-delta 1e-3
            --max-rel-final-delta-delta 1e-3
            --max-rel-integral-delta-delta 1e-3
            --max-rel-steady-mean-delta-delta 1e-3
            --max-breakdown-time-delta-delta-s 1e-8
            --report-json ${CMAKE_BINARY_DIR}/vacuum_arc_cross_platform_consistency.json
            --report-md ${CMAKE_BINARY_DIR}/vacuum_arc_cross_platform_consistency.md
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
    add_subdirectory("Toolkit/Radiation/test")
    add_subdirectory("Toolkit/Surface Charging/test")
    add_subdirectory("Toolkit/Internal Charging/test")
    add_subdirectory("Toolkit/Vacuum Arc/test")
endif()
