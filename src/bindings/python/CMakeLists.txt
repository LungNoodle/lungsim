
include(UseSWIG)

if (match_python_deployment)
  execute_process(
    COMMAND ${Python_EXECUTABLE} -c
            "import sysconfig; print(sysconfig.get_config_var('MACOSX_DEPLOYMENT_TARGET') or '')"
    OUTPUT_VARIABLE PYTHON_MACOS_TARGET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()

if(PYTHON_MACOS_TARGET)
  set(CMAKE_OSX_DEPLOYMENT_TARGET ${PYTHON_MACOS_TARGET} CACHE STRING "" FORCE)
endif()


# SWIG options (target-local, not global)
set(SWIG_COMMON_FLAGS "-relativeimport")

# OS classifier
if(WIN32)
  set(SETUP_PY_OS_CLASSIFIER "Operating System :: Microsoft :: Windows")
elseif(APPLE)
  set(SETUP_PY_OS_CLASSIFIER "Operating System :: MacOS :: MacOS X")
else()
  set(SETUP_PY_OS_CLASSIFIER "Operating System :: Unix")
endif()

# Interface files
set(INTERFACE_SRCS
  ../interface/arrays.i
  ../interface/diagnostics.i
  ../interface/exports.i
  ../interface/filenames.i
  ../interface/field_utilities.i
  ../interface/gas_exchange.i
  ../interface/gasmix.i
  ../interface/geometry.i
  ../interface/growtree.i
  ../interface/indices.i
  ../interface/imports.i
  ../interface/lymphatics.i
  ../interface/parameter_types.i
  ../interface/pressure_resistance_flow.i
  ../interface/species_transport.i
  ../interface/surface_fitting.i
  ../interface/ventilation.i
  ../interface/wave_transmission.i
)

if(IS_MULTI_CONFIG)
  set(_CFG_DIR "/$<CONFIG>")
endif()

set(PYTHON_PACKAGE_DIR ${CMAKE_CURRENT_BINARY_DIR}${_CFG_DIR}/package)

# Output directories
set(SWIG_OUT_DIR ${PYTHON_PACKAGE_DIR}/aether)
set(SWIG_GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/swig)

file(MAKE_DIRECTORY ${SWIG_OUT_DIR})
file(MAKE_DIRECTORY ${SWIG_GEN_DIR})

set(SWIG_TARGETS)
get_target_property(AETHER_C_INCLUDES aether_c INTERFACE_INCLUDE_DIRECTORIES)

foreach(interface ${INTERFACE_SRCS})
  get_filename_component(module ${interface} NAME_WE)

  swig_add_library(${module}
    TYPE MODULE
    LANGUAGE python
    OUTPUT_DIR ${SWIG_OUT_DIR}
    OUTFILE_DIR ${SWIG_GEN_DIR}
    SOURCES ${interface}
  )

  # Apply SWIG flags per target
  set_target_properties(${module} PROPERTIES
    SWIG_COMPILE_OPTIONS "${SWIG_COMMON_FLAGS}"
    SWIG_INCLUDE_DIRECTORIES ${AETHER_C_INCLUDES}
  )

  target_link_libraries(${module}
    PRIVATE
      Python::Module
      Python::NumPy
      aether_c
  )

  set_target_properties(${module} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${SWIG_OUT_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${SWIG_OUT_DIR}
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
  )

  list(APPEND SWIG_TARGETS ${module})
  set(_last_swig_module ${module})

endforeach()

add_custom_command(TARGET ${_last_swig_module} POST_BUILD
  COMMENT "Copying aether_c to Python package dir"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${SWIG_OUT_DIR}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_FILE:aether_c>
    ${SWIG_OUT_DIR}
)

if(WIN32)
  set(SETUP_PY_OPERATING_SYSTEM_CLASSIFIER "Operating System :: Microsoft :: Windows")
elseif(DARWIN)
  set(SETUP_PY_OPERATING_SYSTEM_CLASSIFIER "Operating System :: MacOS :: MacOS X")
else()
  set(SETUP_PY_OPERATING_SYSTEM_CLASSIFIER "Operating System :: Unix")
endif()

set(CLASSIFIERS
"  \"Programming Language :: Python :: 3\",
  \"Development Status :: 4 - Beta\",
  \"Intended Audience :: Developers\",
  \"Intended Audience :: Education\",
  \"Intended Audience :: Science/Research\",
  \"Programming Language :: Python\",
  \"Programming Language :: Python :: ${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}\",
  \"${SETUP_PY_OPERATING_SYSTEM_CLASSIFIER}\",
  \"Topic :: Scientific/Engineering :: Medical Science Apps.\",
  \"Topic :: Software Development :: Libraries :: Python Modules\",
  \"Operating System :: Unix\""
)

set(PACKAGE_PYPROJECT_TOML_FILE ${PYTHON_PACKAGE_DIR}/pyproject.toml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pyproject.in.toml ${PACKAGE_PYPROJECT_TOML_FILE} @ONLY)
set(PACKAGE_INIT_FILE ${SWIG_OUT_DIR}/__init__.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${PACKAGE_INIT_FILE} @ONLY)
set(PACKAGE_SETUP_PY_FILE ${PYTHON_PACKAGE_DIR}/setup.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${PACKAGE_SETUP_PY_FILE} COPYONLY)
set(PACKAGE_README_FILE ${PYTHON_PACKAGE_DIR}/README.rst)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../../README.rst ${PACKAGE_README_FILE} COPYONLY)
set(PACKAGE_LICENCE_FILE ${PYTHON_PACKAGE_DIR}/LICENCE.txt)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../../LICENCE.txt ${PACKAGE_LICENCE_FILE} COPYONLY)

# Aggregate target
add_custom_target(pybindings ALL DEPENDS ${SWIG_TARGETS})

# Install compiled modules
install(TARGETS ${SWIG_TARGETS}
  LIBRARY DESTINATION lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/aether
  RUNTIME DESTINATION lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/aether
)

# Install Python support files
install(
  DIRECTORY ${PYTHON_PACKAGE_DIR}/
  DESTINATION lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/aether
  FILES_MATCHING PATTERN "*.py"
)
install( FILES
  ${PACKAGE_LICENCE_FILE}
  ${PACKAGE_README_FILE}
  ${PACKAGE_PYPROJECT_TOML_FILE}
  ${PACKAGE_INIT_FILE}
  DESTINATION lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/aether
)

# Copy core library into package
install(FILES $<TARGET_FILE:aether_c>
  DESTINATION lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/aether
)
