

cmake_minimum_required(VERSION 3.2)
project(Aether VERSION 0.1.0 LANGUAGES Fortran C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Required for creating documentation from restructured text files
find_package(Sphinx)

# Create additional build types
include(buildtypes)

# Build shared or static library
if(NOT BUILD_SHARED)
  set(BUILD_SHARED ON)
endif()
set(AETHER_BUILD_SHARED ${BUILD_SHARED} CACHE BOOL "Build shared library (true) or static library (false)")
set(BUILD_SHARED_LIBS ${AETHER_BUILD_SHARED} CACHE INTERNAL "Internalise BUILD_SHARED_LIBS, manipulate via AETHER_BUILD_SHARED" FORCE)

# Set the build type Release, Debug plus more if you define them
if(NOT BUILD_TYPE)
  set(BUILD_TYPE Release)
endif()
set(AETHER_BUILD_TYPE ${BUILD_TYPE} CACHE STRING "Set the build type; [Release], Debug, Pedantic (GNU Fortran only)")
set(CMAKE_BUILD_TYPE ${AETHER_BUILD_TYPE} CACHE INTERNAL "Internalise CMAKE_BUILD_TYPE, manipulate via AETHER_BUILD_TYPE" FORCE)

# Set whether to build bindgins or not, may require SWIG
option(AETHER_BUILD_BINDINGS "Build bindings for ${PROJECT}, may require SWIG" YES)

# Set base install location
if(NOT INSTALL_PREFIX)
  set(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()
set(AETHER_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE STRING  "Set the base install prefix")
set(CMAKE_INSTALL_PREFIX ${AETHER_INSTALL_PREFIX} CACHE INTERNAL "Internlaise CMAKE_INSTALL_PREFIX, manipulate via AETHER_INSTALL_PREFIX" FORCE)
mark_as_advanced(AETHER_INSTALL_PREFIX)

# Turn on the ability to create folders to organize projects (.vcproj)
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# The list of all source files that are required to build the library
set(LIB_SRCS
  src/arrays.f90
  src/diagnostics.f90
  src/exports.f90
  src/filenames.f90
  src/gasmix.f90
  src/geometry.f90
  src/indices.f90
  src/math_utilities.f90
  src/other_consts.f90
  src/pressure_resistance_flow.f90
  src/solve.f90
  src/ventilation.f90
)

# Create an object library so we can use the objects in two different libraries
add_library(aether_obj OBJECT ${LIB_SRCS})
set_target_properties(aether_obj PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Modules"
  POSITION_INDEPENDENT_CODE TRUE)

if (MSVC)
  # Visual Studio requires a dummy file as it will not build a library from only objects.
  set(dummy_file "${CMAKE_CURRENT_BINARY_DIR}/dummy.f90")
  file(WRITE ${dummy_file} "module dummy
implicit none
contains
end module dummy\n")
endif ()

# Create the aether fortran library
add_library(aether $<TARGET_OBJECTS:aether_obj> ${dummy_file})

set_target_properties(aether PROPERTIES
  Fortran_VISIBILITY_PRESET hidden
  C_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN 1)

if(AETHER_BUILD_BINDINGS)
  add_subdirectory(bindings)
endif()

if(SPHINX_FOUND)
  if(NOT DEFINED SPHINX_THEME)
    set(SPHINX_THEME default)
  endif()

  if(NOT DEFINED SPHINX_THEME_DIR)
    set(SPHINX_THEME_DIR)
  endif()

  # Sphinx cache with pickled ReST documents
  set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

  # HTML output directory
  set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

  # Sphinx source directory
  set(SPHINX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Documentation")

  add_custom_target(docs
    ${SPHINX_EXECUTABLE}
      -q -b html
      -d "${SPHINX_CACHE_DIR}"
      "${SPHINX_SOURCE_DIR}"
      "${SPHINX_HTML_DIR}"
    COMMENT "Building HTML documentation with Sphinx")

endif()

