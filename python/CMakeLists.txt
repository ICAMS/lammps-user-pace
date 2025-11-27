cmake_minimum_required(VERSION 3.15)

# ============================================================
#   Enforce Debug builds and full DWARF on macOS and Linux
# ============================================================

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_MACOSX_RPATH ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -fno-omit-frame-pointer")
set(CMAKE_C_FLAGS_DEBUG   "-g3 -O0 -fno-omit-frame-pointer")

# Disable stripping globally
set(CMAKE_STRIP ":")

add_compile_options("$<$<CONFIG:Debug>:-g3>")
add_compile_options("$<$<CONFIG:Debug>:-fno-omit-frame-pointer>")

# ============================================================
# Continue with project() and targets
# ============================================================

project(lammps_pyace_python LANGUAGES CXX)

# Find Python
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# -------------------------------- yaml-cpp --------------------------------

find_package(yaml-cpp QUIET)

if(NOT yaml-cpp_FOUND)
    # enforce building libyaml-cpp as static library and turn off optional features
    set(YAML_BUILD_SHARED_LIBS OFF)
    set(YAML_CPP_BUILD_CONTRIB OFF)
    set(YAML_CPP_BUILD_TOOLS OFF)
    add_subdirectory(yaml-cpp build-yaml-cpp)
    add_library(yaml-cpp::yaml-cpp ALIAS yaml-cpp-pace)

    install(DIRECTORY yaml-cpp/include/yaml-cpp
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
      FILES_MATCHING
      PATTERN *.h
    )
    set(BUILD_YAML_CPP TRUE)
else()
    find_package(yaml-cpp REQUIRED)
    set(BUILD_YAML_CPP FALSE)
endif()

# -------------------------------- cnpy --------------------------------

find_package(cnpy QUIET)

if(NOT cnpy_FOUND)
    add_library(cnpy-static STATIC ../cnpy/cnpy.cpp)
    target_include_directories(cnpy-static PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/pace>
    )
    set_target_properties(cnpy-static PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME cnpy)
    add_library(cnpy::cnpy ALIAS cnpy-static)

    install(DIRECTORY cnpy
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pace
      FILES_MATCHING
      PATTERN *.h
    )
    set(BUILD_CNPY TRUE)
else()
    find_package(cnpy REQUIRED)
    set(BUILD_CNPY FALSE)
endif()

# -------------------------------- wigner-cpp --------------------------------
# this is header-only library

set(WIGNER_PATH ../wigner-cpp)
set(WIGNER_INCLUDE_PATH ${WIGNER_PATH}/include)

install(DIRECTORY wigner
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pace
  FILES_MATCHING
  PATTERN *.hpp
)

# -------------------------------- pybind --------------------------------

include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)

# -------------------------------- basis --------------------------------

message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

set(ML_PACE "../ML-PACE")
set(ACE ${ML_PACE}/ace)
set(ACE_EVALUATOR ${ML_PACE}/ace-evaluator)
set(HELPERS "helpers")

include_directories(${ML_PACE})
include_directories(${HELPERS})

set(SOURCES_BASIS "${ACE_EVALUATOR}/ace_c_basis.cpp"
        "${ACE_EVALUATOR}/ace_radial.cpp"
        "${ACE_EVALUATOR}/ace_spherical_cart.cpp"
        "${ACE_EVALUATOR}/ace_abstract_basis.cpp"
        "${ACE_EVALUATOR}/ace_flatten_basis.cpp"
        "${ACE_EVALUATOR}/ships_radial.cpp"
        "${ACE}/ace_b_basis.cpp"
        "${ACE}/ace_b_basisfunction.cpp"
        "${ACE}/ace_clebsch_gordan.cpp"
        "${ACE}/ace_couplings.cpp"
        "${ACE}/ace_yaml_input.cpp"
        "${HELPERS}/ace_c_basis_binding.cpp"
        "${HELPERS}/ace_bbasis_spec_helper.cpp"
        "${HELPERS}/ace_radial_helper.cpp"
        "${HELPERS}/ace_c_basisfunction_helper.cpp"
        "${HELPERS}/ace_c_basis_helper.cpp"
  )
  
pybind11_add_module(basis ${SOURCES_BASIS} )
target_include_directories(basis PRIVATE ${WIGNER_INCLUDE_PATH})
target_link_libraries(basis PRIVATE yaml-cpp::yaml-cpp cnpy::cnpy)
target_compile_options(basis PRIVATE -Wno-deprecated-declarations)
set(deps_targets "")

if(BUILD_YAML_CPP)
    list(APPEND deps_targets yaml-cpp-pace)
else()
    target_compile_definitions(basis PUBLIC YAML_PACE=YAML)
    set_target_properties(basis PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
    )
endif()

if(BUILD_CNPY)
    list(APPEND deps_targets cnpy-static)
else()
    file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/ML-PACE/cnpy)
    file(WRITE ${PROJECT_SOURCE_DIR}/ML-PACE/cnpy/cnpy.h "#include <cnpy.h>")
endif()





# Honor the directory provided by setup.py
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lammps_pyace")
endif()
