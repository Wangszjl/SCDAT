# Diagnostics unit tests

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

add_executable(Diagnostics_test Diagnostics_test.cpp)
target_link_libraries(Diagnostics_test
    PRIVATE
        SCDAT::Tools::Diagnostics
        SCDAT::Tools::Geometry
        SCDAT::Tools::Particle
        GTest::gtest_main
)
gtest_discover_tests(Diagnostics_test)
