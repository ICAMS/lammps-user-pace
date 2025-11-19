cmake_minimum_required(VERSION 3.15)
project(lammps_pyace_python)

# Find Python
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Find or fetch pybind11
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found, fetching from GitHub")
    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.11.1
    )
    FetchContent_MakeAvailable(pybind11)
endif()

# Find or build parent ML-PACE library
find_package(pace CONFIG QUIET)
if(NOT pace_FOUND)
    set(PACE_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/ml-pace-build)
    if(NOT TARGET pace::pace)
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PACE_BUILD_DIR})
    endif()
endif()


pybind11_add_module(_basis
    ace_c_basis_binding.cpp
    helpers/ace_bbasis_spec_helper.cpp
    helpers/ace_c_basis_helper.cpp
    helpers/ace_c_basisfunction_helper.cpp
    helpers/ace_radial_helper.cpp
)

target_link_libraries(_basis PRIVATE pace::pace)

target_include_directories(_basis PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/helpers
    ${CMAKE_CURRENT_SOURCE_DIR}/../ML-PACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../wigner-cpp/include
)

# Honor the directory provided by setup.py
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lammps_pyace")
endif()

if(APPLE)
    set_target_properties(_basis PROPERTIES SUFFIX ".so")
endif()
