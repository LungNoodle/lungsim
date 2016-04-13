

cmake_minimum_required(VERSION 3.2)
project(Aether Fortran CXX)

# Build shared or static library
if(NOT BUILD_SHARED)
  set(BUILD_SHARED ON)
endif()
set(VENTILATION_BUILD_SHARED ${BUILD_SHARED} CACHE BOOL "Build shared library (true) or static library (false)")
set(BUILD_SHARED_LIBS ${VENTILATION_BUILD_SHARED} CACHE INTERNAL "Internalise BUILD_SHARED_LIBS, manipulate via VENTILATION_BUILD_SHARED" FORCE)

# Set the build type Release, Debug plus more if you define them
if(NOT BUILD_TYPE)
  set(BUILD_TYPE Release)
endif()
set(VENTILATION_BUILD_TYPE ${BUILD_TYPE} CACHE STRING "Set the build type; [Release], Debug, Pedantic (GNU Fortran only)")
set(CMAKE_BUILD_TYPE ${VENTILATION_BUILD_TYPE} CACHE INTERNAL "Internalise CMAKE_BUILD_TYPE, manipulate via VENTILATION_BUILD_TYPE" FORCE)

# Set base install location
if(NOT INSTALL_PREFIX)
  set(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()
set(VENTILATION_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE STRING  "Set the base install prefix")
set(CMAKE_INSTALL_PREFIX ${VENTILATION_INSTALL_PREFIX} CACHE INTERNAL "Internlaise CMAKE_INSTALL_PREFIX, manipulate via VENTILATION_INSTALL_PREFIX" FORCE)

set(LIB_SRCS
  src/arrays.f90
  src/diagnostics.f90
  src/exports.f90
  src/filenames.f90
  src/gasmixmod.f90
  src/geometry.f90
  src/indices.f90
  src/other_consts.f90
  src/solve.f90
  src/ventilation.f90
)

add_library(aether ${LIB_SRCS})#Aether to be replaced with whatever the final name of the code is

set_target_properties(aether PROPERTIES 
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Modules"
  Fortran_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN 1)
