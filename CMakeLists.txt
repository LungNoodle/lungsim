

cmake_minimum_required(VERSION 3.2)
project(Aether VERSION 0.1.0 LANGUAGES Fortran C CXX)

# Policy settings

if (POLICY CMP0078)
  cmake_policy(SET CMP0078 OLD)
endif ()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Required for creating documentation from restructured text files
find_package(Sphinx)

# Helpful functions
include(miscfunctions)

# Create additional build types
include(buildtypes)

# Mark some variables as advanced hiding them from basic users
tidy_gui()

# Build shared or static library
if(NOT BUILD_SHARED)
  set(BUILD_SHARED ON)
endif()
set(AETHER_BUILD_SHARED ${BUILD_SHARED} CACHE BOOL "Build shared library (true) or static library (false)")
set(BUILD_SHARED_LIBS ${AETHER_BUILD_SHARED} CACHE INTERNAL "Internalise BUILD_SHARED_LIBS, manipulate via AETHER_BUILD_SHARED" FORCE)

# Set the build type Release, Debug plus more if you define them
if (NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT BUILD_TYPE)
    set(BUILD_TYPE Release)
  endif()
  set(AETHER_BUILD_TYPE ${BUILD_TYPE} CACHE STRING "Set the build type; [Release], Debug, Pedantic (GNU Fortran only)")
  set(CMAKE_BUILD_TYPE ${AETHER_BUILD_TYPE} CACHE INTERNAL "Internalise CMAKE_BUILD_TYPE, manipulate via AETHER_BUILD_TYPE" FORCE)
endif()

# Set whether to build bindings or not, may require SWIG
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

# All the source files for the library and bindings are located under this directory.
add_subdirectory(src)

# Must enable testing at the lowest level in CMake.
enable_testing()

# All the tests for the library are located under this directory.
add_subdirectory(tests)

# If Sphinx was found we can add some targets for building the documentation.
set(ADDITIONAL_GENERATED_FILES)
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
  set(SPHINX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/documentation")

  # Add a target 'docs' that will build the documentation inside the build directory on the local computer.
  add_custom_target(docs
    ${SPHINX_EXECUTABLE}
      -q -b html
      -d "${SPHINX_CACHE_DIR}"
      "${SPHINX_SOURCE_DIR}"
      "${SPHINX_HTML_DIR}"
    WORKING_DIRECTORY "${SPHINX_SOURCE_DIR}"
    COMMENT "Building HTML documentation with Sphinx")

  # Remove Sphinx generated files when the 'clean' target is used.
  list(APPEND ADDITIONAL_GENERATED_FILES html _doctrees)

endif()

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${ADDITIONAL_GENERATED_FILES}")
