cmake_minimum_required(VERSION 3.20)
project(vela-tcad VERSION 0.1.0 LANGUAGES CXX)

# ---------------------------------------------------------------------------
# C++20
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# Dependencies (installed system packages)
# ---------------------------------------------------------------------------

# Eigen 3 (header-only linear algebra)
find_package(Eigen3 REQUIRED NO_MODULE)

# nlohmann/json (header-only JSON)
find_package(nlohmann_json REQUIRED)

# Catch2 v3 (unit test framework)
find_package(Catch2 3 REQUIRED)

# ---------------------------------------------------------------------------
# Library: vela_core
# ---------------------------------------------------------------------------
add_library(vela_core STATIC
    src/core/ScalingSystem.cpp
    src/material/Material.cpp
    src/material/MaterialDatabase.cpp
    src/mesh/DeviceMesh.cpp
    src/mesh/BoxGeometryBuilder.cpp
    src/io/MeshReader.cpp
    src/io/VTKWriter.cpp
    src/io/CSVWriter.cpp
    src/physics/DopingModel.cpp
    src/physics/CarrierStatistics.cpp
    src/physics/SRHRecombination.cpp
    src/physics/MobilityModel.cpp
    src/physics/RecombinationModel.cpp
    src/physics/BandgapNarrowing.cpp
    src/equation/PoissonAssembler.cpp
    src/equation/DDAssembler.cpp
    src/discretization/Bernoulli.cpp
    src/discretization/ScharfetterGummel.cpp
    src/solver/LinearSolver.cpp
    src/solver/GummelSolver.cpp
    src/post/ContactCurrent.cpp
    src/solver/NewtonSolver.cpp
    src/equation/CoupledDDAssembler.cpp
    src/numerics/ResidualNorm.cpp
    src/numerics/LineSearch.cpp
    src/simulation/PoissonSimulation.cpp
    src/simulation/DCSweep.cpp
)

target_include_directories(vela_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(vela_core PUBLIC
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
enable_testing()
include(CTest)

add_executable(test_scaling tests/test_scaling.cpp)
target_link_libraries(test_scaling PRIVATE vela_core Catch2::Catch2WithMain)

add_executable(test_mesh tests/test_mesh.cpp)
target_link_libraries(test_mesh PRIVATE vela_core Catch2::Catch2WithMain)

add_executable(test_poisson tests/test_poisson.cpp)
target_link_libraries(test_poisson PRIVATE vela_core Catch2::Catch2WithMain)

add_executable(test_box_geometry tests/test_box_geometry.cpp)
target_link_libraries(test_box_geometry PRIVATE vela_core Catch2::Catch2WithMain)

add_executable(test_bernoulli tests/test_bernoulli.cpp)
target_link_libraries(test_bernoulli PRIVATE vela_core Catch2::Catch2WithMain)

add_executable(test_sg_flux tests/test_sg_flux.cpp)
target_link_libraries(test_sg_flux PRIVATE vela_core Catch2::Catch2WithMain)

add_executable(test_dd_gummel tests/test_dd_gummel.cpp)
target_link_libraries(test_dd_gummel PRIVATE vela_core Catch2::Catch2WithMain)

add_executable(test_newton_solver tests/test_newton_solver.cpp)
target_link_libraries(test_newton_solver PRIVATE vela_core Catch2::Catch2WithMain)

add_executable(test_mobility tests/test_mobility.cpp)
target_link_libraries(test_mobility PRIVATE vela_core Catch2::Catch2WithMain)

add_executable(test_recombination tests/test_recombination.cpp)
target_link_libraries(test_recombination PRIVATE vela_core Catch2::Catch2WithMain)

add_executable(test_dc_sweep tests/test_dc_sweep.cpp)
target_link_libraries(test_dc_sweep PRIVATE vela_core Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(test_scaling)
catch_discover_tests(test_mesh)
catch_discover_tests(test_poisson)
catch_discover_tests(test_box_geometry)
catch_discover_tests(test_bernoulli)
catch_discover_tests(test_sg_flux)
catch_discover_tests(test_dd_gummel)
catch_discover_tests(test_newton_solver)
catch_discover_tests(test_mobility)
catch_discover_tests(test_recombination)
catch_discover_tests(test_dc_sweep)
