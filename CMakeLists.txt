cmake_minimum_required(VERSION 3.10)

project(libpace VERSION 20240911 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

##yaml
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

## cnpy
find_package(cnpy QUIET)

if(NOT cnpy_FOUND)
    add_library(cnpy-static STATIC cnpy/cnpy.cpp)
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

## winger-cpp
# this is header-only library
set(WIGNER_PATH wigner-cpp)
set(WIGNER_INCLUDE_PATH ${WIGNER_PATH}/include)

install(DIRECTORY wigner
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pace
  FILES_MATCHING
  PATTERN *.hpp
)

# ML-PACE sources
file(GLOB PACE_EVALUATOR_SOURCES ML-PACE/ace-evaluator/*.cpp)
file(GLOB PACE_SOURCES ML-PACE/ace/*.cpp)

add_library(pace STATIC ${PACE_EVALUATOR_SOURCES} ${PACE_SOURCES})
target_include_directories(pace PRIVATE ${WIGNER_INCLUDE_PATH})
target_compile_definitions(pace PUBLIC EXTRA_C_PROJECTIONS) # important for B-projections and extrapolation grade calculations
target_compile_definitions(pace PUBLIC COMPUTE_B_GRAD) # important for gradients of B-projections and ctilde functions
target_link_libraries(pace PUBLIC yaml-cpp::yaml-cpp cnpy::cnpy)

set(deps_targets "")

if(BUILD_YAML_CPP)
    list(APPEND deps_targets yaml-cpp-pace)
else()
    target_compile_definitions(pace PUBLIC YAML_PACE=YAML)
    set_target_properties(pace PROPERTIES
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



include(GNUInstallDirs)

# ML-PACE includes
target_include_directories(pace
PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/ML-PACE>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/pace>
)

install(DIRECTORY ML-PACE/ace ML-PACE/ace-evaluator
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pace
  FILES_MATCHING
  PATTERN *.h
)

install(TARGETS pace ${deps_targets} EXPORT pace_Targets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(EXPORT pace_Targets FILE pace-targets.cmake
        NAMESPACE pace:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pace
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/pace-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/pace-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pace
)

write_basic_package_version_file("pace-config-version.cmake"
                                 VERSION ${PROJECT_VERSION} COMPATIBILITY ExactVersion)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/pace-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/pace-config-version.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Findcnpy.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pace)

add_library(pace::pace ALIAS pace)
