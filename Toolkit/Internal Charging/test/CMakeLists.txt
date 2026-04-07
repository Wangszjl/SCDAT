# Internal Charging smoke tests

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

add_executable(InternalCharging_smoke_test InternalCharging_smoke_test.cpp)
target_link_libraries(InternalCharging_smoke_test
    PRIVATE
        SCDAT::Toolkit::Internal_Charging
        SCDAT::Tools::Output
        GTest::gtest_main
)
gtest_discover_tests(InternalCharging_smoke_test)

add_executable(InternalCharging_radiation_coupling_smoke_test
    InternalCharging_radiation_coupling_smoke_test.cpp
)
target_link_libraries(InternalCharging_radiation_coupling_smoke_test
    PRIVATE
        SCDAT::Toolkit::Internal_Charging
        SCDAT::Toolkit::Radiation
    SCDAT::Tools::Coupling
        SCDAT::Tools::Output
        GTest::gtest_main
)
gtest_discover_tests(InternalCharging_radiation_coupling_smoke_test)
