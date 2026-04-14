# Material unit tests

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

add_executable(Material_test Material_test.cpp)
target_link_libraries(Material_test
    PRIVATE
        SCDAT::Tools::Geometry
        SCDAT::Tools::Material
        SCDAT::Tools::FieldSolver
        GTest::gtest_main
)
gtest_discover_tests(Material_test)
