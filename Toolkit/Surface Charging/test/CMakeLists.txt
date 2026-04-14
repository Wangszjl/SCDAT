# Surface Charging smoke tests

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

add_executable(SurfaceCharging_smoke_test SurfaceCharging_smoke_test.cpp)
target_link_libraries(SurfaceCharging_smoke_test
    PRIVATE
        SCDAT::Toolkit::Surface_Charging
        SCDAT::Tools::Output
        GTest::gtest_main
)
gtest_discover_tests(SurfaceCharging_smoke_test)

add_executable(SurfaceCharging_object_layer_test SurfaceCharging_object_layer_test.cpp)
target_link_libraries(SurfaceCharging_object_layer_test
    PRIVATE
        SCDAT::Toolkit::Surface_Charging
        GTest::gtest_main
)
gtest_discover_tests(SurfaceCharging_object_layer_test)
