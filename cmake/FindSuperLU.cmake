
# Umfpack lib usually requires linking to a blas library.

find_package(BLAS)

if (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)
  set(SUPERLU_FIND_QUIETLY TRUE)
endif (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)

find_path(SUPERLU_INCLUDES
  NAMES
  supermatrix.h
  PATHS
  $ENV{SUPERLUDIR}
  ${INCLUDE_INSTALL_DIR}
  ${SUPERLU_INSTALL_DIR}
  PATH_SUFFIXES
  superlu
  include/superlu
  SRC
)

if (SuperLU_FIND_VERSION_EXACT)
  set(_SUPERLU_KNOWN_VERSIONS ${SuperLU_FIND_VERSION})
else()
  set(_SUPERLU_KNOWN_VERSIONS ${SuperLU_FIND_VERSION} "5.2.1" "5.2" "5.1.1" "5.1" "5.0" "4.3" "4.2" "4.2" "4.0" "3.1" "3.0")
endif()

set(_SUPERLU_LIBRARY_NAMES)
foreach(_VERSION ${_SUPERLU_KNOWN_VERSIONS})
  list(APPEND _SUPERLU_LIBRARY_NAMES "superlu_${_VERSION}" "superlu-${_VERSION}")
endforeach()
list(APPEND _SUPERLU_LIBRARY_NAMES "superlu")
unset(_VERSION)
unset(_SUPERLU_KNOWN_VERSIONS)

find_library(SUPERLU_LIBRARIES
  NAMES ${_SUPERLU_LIBRARY_NAMES}
  PATHS $ENV{SUPERLUDIR} ${LIB_INSTALL_DIR}
  PATH_SUFFIXES lib)

if(SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state()

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${SUPERLU_INCLUDES})

# check whether struct mem_usage_t is globally defined
check_cxx_source_compiles("
typedef int int_t;
#include <supermatrix.h>
#include <slu_util.h>
int main() {
  mem_usage_t mem;
  return 0;
}"
SUPERLU_HAS_GLOBAL_MEM_USAGE_T)


check_cxx_source_compiles("
typedef int int_t;
#include <supermatrix.h>
#include <superlu_enum_consts.h>
int main() {
  return SLU_SINGLE;
}"
SUPERLU_HAS_CLEAN_ENUMS)

check_cxx_source_compiles("
typedef int int_t;
#include <supermatrix.h>
#include <slu_util.h>
int main(void)
{
  GlobalLU_t glu;
  return 0;
}"
SUPERLU_HAS_GLOBALLU_T)

if(SUPERLU_HAS_GLOBALLU_T)
  # at least 5.0
  set(SUPERLU_VERSION_VAR "5.0")
elseif(SUPERLU_HAS_CLEAN_ENUMS)
  # at least 4.3
  set(SUPERLU_VERSION_VAR "4.3")
elseif(SUPERLU_HAS_GLOBAL_MEM_USAGE_T)
  # at least 4.0
  set(SUPERLU_VERSION_VAR "4.0")
else()
  set(SUPERLU_VERSION_VAR "3.0")
endif()

cmake_pop_check_state()

if(SuperLU_FIND_VERSION)
  if(${SUPERLU_VERSION_VAR} VERSION_LESS ${SuperLU_FIND_VERSION})
    set(SUPERLU_VERSION_OK FALSE)
  else()
    set(SUPERLU_VERSION_OK TRUE)
  endif()
else()
  set(SUPERLU_VERSION_OK TRUE)
endif()

endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SuperLU
                                  REQUIRED_VARS SUPERLU_INCLUDES SUPERLU_LIBRARIES SUPERLU_VERSION_OK
                                  VERSION_VAR SUPERLU_VERSION_VAR)

mark_as_advanced(SUPERLU_INCLUDES SUPERLU_LIBRARIES)

if(SUPERLU_FOUND)
  if (NOT TARGET superlu)
    if (NOT CMAKE_CFG_INTDIR STREQUAL .)
        string(TOUPPER ${CMAKE_CFG_INTDIR} _CURRENT_BUILD_TYPE)
    elseif(CMAKE_BUILD_TYPE)
        string(TOUPPER ${CMAKE_BUILD_TYPE} _CURRENT_BUILD_TYPE)
    else()
        set(_CURRENT_BUILD_TYPE NOCONFIG)
    endif()
    set(_FIRST_LIB ${SUPERLU_LIBRARIES})
    # I'm not aware of a framework for SuperLU but just in case.
    # Treat apple frameworks separate
    # See http://stackoverflow.com/questions/12547624/cant-link-macos-frameworks-with-cmake
    if(APPLE AND ${_FIRST_LIB} MATCHES ".framework$")
        string(REGEX REPLACE ".*/([A-Za-z0-9.]+).framework$" "\\1" FW_NAME ${_FIRST_LIB})
        #message(STATUS "Matched '${FW_NAME}' to ${LIB}")
        set(_FIRST_LIB "${_FIRST_LIB}/${FW_NAME}")
    endif()

    add_library(superlu UNKNOWN IMPORTED)
    if (TARGET blas)
      append_link_library(superlu blas)
    endif ()
    set_target_properties(superlu PROPERTIES
      IMPORTED_LOCATION_${_CURRENT_BUILD_TYPE} ${_FIRST_LIB}
      IMPORTED_CONFIGURATIONS ${_CURRENT_BUILD_TYPE}
      INTERFACE_INCLUDE_DIRECTORIES ${SUPERLU_INCLUDES}
    )

    unset(_CURRENT_BUILD_TYPE)
    unset(_FIRST_LIB)
  endif()
endif()
