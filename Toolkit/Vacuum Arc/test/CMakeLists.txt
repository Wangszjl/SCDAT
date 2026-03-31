# Vacuum Arc smoke tests

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

add_executable(VacuumArc_smoke_test VacuumArc_smoke_test.cpp)
target_link_libraries(VacuumArc_smoke_test
    PRIVATE
        SCDAT::Toolkit::Vacuum_Arc
        SCDAT::Tools::Output
        GTest::gtest_main
)
gtest_discover_tests(VacuumArc_smoke_test)
