cmake_minimum_required(VERSION 3.10)

project(libpace CXX)

# set policy to silence warnings about ignoring <PackageName>_ROOT but use it
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()
# set policy to silence warnings about missing executable permissions in
# pythonx.y-config when cross-compiling. review occasionally if it may be set to NEW
if(POLICY CMP0109)
    cmake_policy(SET CMP0109 OLD)
endif()


##yaml
# enforce building libyaml-cpp as static library and turn off optional features
set(YAML_BUILD_SHARED_LIBS OFF)
set(YAML_CPP_BUILD_CONTRIB OFF)
set(YAML_CPP_BUILD_TOOLS OFF)
add_subdirectory(yaml-cpp build-yaml-cpp)
set(YAML_CPP_INCLUDE_DIR yaml-cpp/include)

## cnpy
set(CNPY_PATH cnpy)
set(CNPY_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CNPY_SRC ${CNPY_PATH}/cnpy.cpp)
add_library(cnpy-static STATIC ${CNPY_SRC})
set_target_properties(cnpy-static PROPERTIES LINKER_LANGUAGE CXX)

## winger-cpp
# this is header-only library
set(WIGNER_PATH wigner-cpp)
set(WIGNER_INCLUDE_PATH ${WIGNER_PATH}/include)

# ML-PACE includes
file(GLOB PACE_INCLUDE_DIR ML-PACE)

# ML-PACE sources
file(GLOB PACE_EVALUATOR_SOURCES ML-PACE/ace-evaluator/*.cpp)
file(GLOB PACE_SOURCES ML-PACE/ace/*.cpp)
list(FILTER PACE_EVALUATOR_SOURCES EXCLUDE REGEX pair_pace*.cpp)
list(FILTER PACE_SOURCES EXCLUDE REGEX pair_pace*.cpp)

add_library(pace STATIC ${PACE_EVALUATOR_SOURCES} ${PACE_SOURCES})
#set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})
target_include_directories(pace PUBLIC ${PACE_INCLUDE_DIR} ${YAML_CPP_INCLUDE_DIR})
target_include_directories(pace PRIVATE ${CNPY_INCLUDE_PATH} ${WIGNER_INCLUDE_PATH})
target_compile_definitions(pace PUBLIC EXTRA_C_PROJECTIONS) # important for B-projections and extrapolation grade calculations
target_link_libraries(pace PRIVATE yaml-cpp-pace cnpy-static)
