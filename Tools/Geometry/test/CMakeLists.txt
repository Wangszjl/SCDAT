# Geometry单元测试

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
if(MSVC)
    # Windows下使用与项目一致的CRT，避免链接冲突。
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

set(GEOMETRY_TEST_TARGETS
    Matrix_test
    Point3D_test
    Vector3D_test
)

foreach(test_target IN LISTS GEOMETRY_TEST_TARGETS)
    add_executable(${test_target} ${test_target}.cpp)
    target_link_libraries(${test_target} PRIVATE SCDAT::Tools::Geometry GTest::gtest_main)
    gtest_discover_tests(${test_target})
endforeach()
