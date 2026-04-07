# Plasma Analysis smoke tests

include(FetchContent)
if(NOT TARGET GTest::gtest_main)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    )
    if(MSVC)
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    endif()
    FetchContent_MakeAvailable(googletest)
endif()

include(GoogleTest)

add_executable(PlasmaAnalysis_smoke_test PlasmaAnalysis_smoke_test.cpp)
target_link_libraries(PlasmaAnalysis_smoke_test
    PRIVATE
        SCDAT::Toolkit::Plasma_Analysis
        SCDAT::Tools::Output
        GTest::gtest_main
)
gtest_discover_tests(PlasmaAnalysis_smoke_test)

if(BUILD_TESTING)
    if(NOT Python3_EXECUTABLE)
        find_package(Python3 COMPONENTS Interpreter REQUIRED)
    endif()

    add_test(
        NAME PlasmaMultiscaleRegressionGate
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/python/check_plasma_multiscale_regression_gate.py
            --project-root ${CMAKE_SOURCE_DIR}
            --plasma-test-exe $<TARGET_FILE:PlasmaAnalysis_smoke_test>
            --matrix-json ${CMAKE_SOURCE_DIR}/scripts/run/plasma_multiscale_regression_matrix.json
            --baseline-json ${CMAKE_SOURCE_DIR}/results/baselines/plasma/plasma_multiscale_regression_baseline_2026-04-04.json
            --summary-json ${CMAKE_BINARY_DIR}/plasma_multiscale_regression_gate.json
            --summary-md ${CMAKE_BINARY_DIR}/plasma_multiscale_regression_gate.md
    )
endif()
