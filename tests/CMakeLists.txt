# Unit tests for header-only utilities and (later) the refactored ns_core.

find_package(GTest 1.14 CONFIG QUIET)

if(NOT GTest_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/releases/download/v1.15.2/googletest-1.15.2.tar.gz
        URL_HASH SHA256=7b42b4d6ed48810c5362c265a17faebe90dc2373c885e5216439d37927f02926
    )
    set(gtest_force_shared_crt ON CACHE INTERNAL "")
    FetchContent_MakeAvailable(googletest)
endif()

add_executable(ns_tests
    test_math_utils.cpp
    test_config.cpp
    test_legacy_import.cpp
    test_field.cpp
    test_mesh.cpp
    test_io.cpp
    test_physics.cpp
    test_numerics.cpp
)
target_include_directories(ns_tests PRIVATE "${CMAKE_SOURCE_DIR}/inc")
target_link_libraries(ns_tests PRIVATE GTest::gtest_main ns_config ns_core)
target_compile_definitions(ns_tests PRIVATE NS_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(ns_tests PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(GoogleTest)
gtest_discover_tests(ns_tests DISCOVERY_TIMEOUT 60)
