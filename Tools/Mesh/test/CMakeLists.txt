# Mesh unit tests

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

add_executable(MeshPartitioning_test MeshPartitioning_test.cpp)
target_link_libraries(
    MeshPartitioning_test
    PRIVATE
        SCDAT::Tools::Mesh
        SCDAT::Tools::Geometry
        GTest::gtest_main
)
gtest_discover_tests(MeshPartitioning_test)
