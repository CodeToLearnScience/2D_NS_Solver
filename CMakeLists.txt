cmake_minimum_required(VERSION 3.23)
project(TwoDNSolver VERSION 0.2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(NS_BUILD_TESTS "Build unit tests (Google Test)" ON)
option(NS_ENABLE_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

function(ns_default_warnings target)
    target_compile_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra>
        $<$<BOOL:${NS_ENABLE_WARNINGS_AS_ERRORS}>:-Werror>)
endfunction()

# ---------------------------------------------------------------------------
# Legacy solver, frozen as the migration baseline (retired at end of Phase 7).
# Built unchanged from src/ + inc/ so regression goldens stay reproducible.
# ---------------------------------------------------------------------------
file(GLOB NS_LEGACY_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_executable(ns_legacy ${NS_LEGACY_SOURCES})
target_include_directories(ns_legacy PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/inc")
target_compile_options(ns_legacy PRIVATE $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall>)
set_target_properties(ns_legacy PROPERTIES OUTPUT_NAME 2D_NS_Solver)

# ---------------------------------------------------------------------------
# Config layer (Phase 2): typed TOML config + legacy importer.
# ---------------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(
    tomlplusplus
    URL https://github.com/marzer/tomlplusplus/archive/refs/tags/v3.4.0.tar.gz
    URL_HASH SHA256=8517f65938a4faae9ccf8ebb36631a38c1cadfb5efa85d9a72e15b9e97d25155
    EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(tomlplusplus)

add_library(ns_config STATIC
    src/config/config.cpp
    src/config/legacy_import.cpp
)
target_include_directories(ns_config PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(ns_config PUBLIC tomlplusplus::tomlplusplus)
target_compile_features(ns_config PUBLIC cxx_std_23)
ns_default_warnings(ns_config)

add_executable(cfg2toml tools/cfg2toml/main.cpp)
target_link_libraries(cfg2toml PRIVATE ns_config)
ns_default_warnings(cfg2toml)

# ---------------------------------------------------------------------------
# Core infrastructure (Phase 3): fields, mesh/metrics, IO.
# ---------------------------------------------------------------------------
add_library(ns_core STATIC
    src/mesh/mesh.cpp
    src/io/vtk.cpp
    src/io/restart.cpp
    src/physics/transport.cpp
    src/physics/dependent_variables.cpp
    src/physics/gradients.cpp
    src/physics/viscous_flux.cpp
)
target_include_directories(ns_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_compile_features(ns_core PUBLIC cxx_std_23)
ns_default_warnings(ns_core)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(NS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
