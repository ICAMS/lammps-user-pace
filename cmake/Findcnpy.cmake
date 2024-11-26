find_path(CNPY_INCLUDE_DIR cnpy.h PATH_SUFFIXES include)
find_library(CNPY_LIBRARY NAMES cnpy PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cnpy REQUIRED_VARS CNPY_LIBRARY CNPY_INCLUDE_DIR)

if(cnpy_FOUND)
  set(CNPY_LIBRARIES ${CNPY_LIBRARY})
  set(CNPY_INCLUDE_DIRS ${CNPY_INCLUDE_DIR})

  mark_as_advanced(CNPY_LIBRARIES CNPY_INCLUDE_DIRS)

  if(NOT TARGET cnpy::cnpy)
    add_library(cnpy::cnpy UNKNOWN IMPORTED)
    set_target_properties(cnpy::cnpy PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${CNPY_INCLUDE_DIR}")

    set_property(TARGET cnpy::cnpy APPEND PROPERTY
        IMPORTED_LOCATION "${CNPY_LIBRARY}")
  endif()
endif()
