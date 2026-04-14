# Particle unit tests

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

add_executable(ParticleDefinitions_test ParticleDefinitions_test.cpp)
target_link_libraries(ParticleDefinitions_test PRIVATE SCDAT::Tools::Particle SCDAT::Tools::Geometry GTest::gtest_main)
gtest_discover_tests(ParticleDefinitions_test)

add_executable(ParticleSource_test ParticleSource_test.cpp)
target_link_libraries(ParticleSource_test PRIVATE SCDAT::Tools::Particle SCDAT::Tools::Geometry GTest::gtest_main)
gtest_discover_tests(ParticleSource_test)

add_executable(HighOrderParticlePusher_test HighOrderParticlePusher_test.cpp)
target_link_libraries(HighOrderParticlePusher_test PRIVATE SCDAT::Tools::Particle SCDAT::Tools::Geometry GTest::gtest_main)
gtest_discover_tests(HighOrderParticlePusher_test)

add_executable(SurfaceDistributionFunction_test SurfaceDistributionFunction_test.cpp)
target_link_libraries(SurfaceDistributionFunction_test PRIVATE SCDAT::Tools::Particle GTest::gtest_main)
gtest_discover_tests(SurfaceDistributionFunction_test)

