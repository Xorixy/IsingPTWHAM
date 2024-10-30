
cmake_minimum_required(VERSION 3.24)

### Setup project name and language
project(Ising VERSION 1.2
        DESCRIPTION "A CMake project template"
        HOMEPAGE_URL "https://github.com/DavidAce/CMakeTemplate"
        LANGUAGES C CXX)


### Set up options
option(ISING_ENABLE_OPENMP   "Enable OpenMP flags such as -fopenmp"                 OFF)
option(ISING_ENABLE_MPI      "Enables use of MPI (work in progress)"                OFF)
option(ISING_ENABLE_COVERAGE "Adds the --coverage flag (used with github actions)"  OFF)
option(ISING_BUILD_BENCH     "Builds benchmarks in ./bench"                         OFF)
option(ISING_BUILD_EXAMPLES  "Builds examples in ./examples"                        OFF)
option(ISING_BUILD_TESTS     "Builds unit tests in ./tests"                         OFF)
option(ISING_CMAKE_DEBUG     "Print info during CMake configuration"                OFF)



if(ISING_DEBUG_CMAKE)
    ### Print operating system details
    include(cmake/PrintHostInfo.cmake)
endif()

### Add all source files
add_executable(Ising
        source/main.cpp
        source/ising.cpp
        source/rnd.cpp
        source/io.cpp
        source/stats.cpp
        source/pt.cpp
        # ...
        )

target_include_directories(Ising PRIVATE source)


### Find dependencies
find_package(fmt        10.1.0  REQUIRED)
find_package(h5pp       1.11.1  REQUIRED)
find_package(Eigen3     3.4.0   REQUIRED)
find_package(spdlog     1.13.0  REQUIRED)
find_package(pcg-cpp            REQUIRED)
find_package(doctest    2.4.1   REQUIRED)
find_package(nanobench  4.3.11  REQUIRED)
find_package(CLI11      2.3.2   REQUIRED)
find_package(nlohmann_json      3.11.3  REQUIRED)



### Create a helper target that links dependencies (this target can be used to build tests and benchmarks)
add_library(ising-libs INTERFACE)

# Link the libraries to the helper target
target_link_libraries(ising-libs INTERFACE fmt::fmt)
target_link_libraries(ising-libs INTERFACE h5pp::h5pp)
target_link_libraries(ising-libs INTERFACE Eigen3::Eigen )
target_link_libraries(ising-libs INTERFACE spdlog::spdlog )
target_link_libraries(ising-libs INTERFACE pcg-cpp::pcg-cpp )
target_link_libraries(ising-libs INTERFACE doctest::doctest)
target_link_libraries(ising-libs INTERFACE nanobench::nanobench)
target_link_libraries(ising-libs INTERFACE CLI11::CLI11)
target_link_libraries(ising-libs INTERFACE nlohmann_json::nlohmann_json)

if(ISING_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS CXX REQUIRED)
    target_link_libraries(ising-libs INTERFACE OpenMP::OpenMP_CXX)
endif()

if(ISING_ENABLE_MPI)
    find_package(MPI COMPONENTS CXX REQUIRED)
    target_link_libraries(ising-libs INTERFACE MPI::MPI_CXX)
endif()


### Link targets to the main executable
target_link_libraries(Ising PUBLIC ising-libs)


################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git identifiers     ###
################################################################
include(cmake/gitversion.cmake)

add_subdirectory(bench)    # Builds benchmarks if ISING_BUILD_BENCH==ON
add_subdirectory(examples) # Builds examples if ISING_BUILD_EXAMPLES==ON


### CTest
if(ISING_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
    if(ISING_ENABLE_COVERAGE)
        # Coverage is used to check what fraction of the source code is tested
        # We use this option in the github actions when compiling on Ubuntu with gcc
        target_compile_options(ising-libs INTERFACE --coverage)
        target_link_options(ising-libs INTERFACE --coverage)
    endif()
endif()


if(ISING_CMAKE_DEBUG)
    # Print summary of CMake configuration
    include(cmake/PrintTargetInfo.cmake)
    print_and_write_project_summary(Ising)
endif()
