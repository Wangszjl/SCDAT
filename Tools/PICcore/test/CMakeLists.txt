# PICcore unit tests

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

set(PICCORE_TEST_TARGETS
    PICCycleBoundary_test
)

foreach(test_target IN LISTS PICCORE_TEST_TARGETS)
    add_executable(${test_target} ${test_target}.cpp)
    target_link_libraries(${test_target}
        PRIVATE
            SCDAT::Tools::Boundary
            SCDAT::Tools::FieldSolver
            SCDAT::Tools::Geometry
            SCDAT::Tools::Mesh
            SCDAT::Tools::PICcore
            SCDAT::Tools::Particle
            GTest::gtest_main
    )
    gtest_discover_tests(${test_target})
endforeach()
