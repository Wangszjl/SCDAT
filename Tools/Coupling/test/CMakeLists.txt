# Coupling unit tests

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

add_executable(Coupling_test Coupling_test.cpp)
target_link_libraries(Coupling_test
    PRIVATE
        SCDAT::Tools::Coupling
        GTest::gtest_main
)
gtest_discover_tests(Coupling_test)

add_executable(BenchmarkContracts_test BenchmarkContracts_test.cpp)
target_link_libraries(BenchmarkContracts_test
    PRIVATE
        SCDAT::Tools::Coupling
        GTest::gtest_main
)
gtest_discover_tests(BenchmarkContracts_test)
