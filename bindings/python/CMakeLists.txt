
project(AetherPythonBindings LANGUAGES C)

# SWIG has already been found at this point.
include(${SWIG_USE_FILE})

string(SUBSTRING ${PYTHONLIBS_VERSION_STRING} 0 1 PYTHONLIBS_MAJOR_VERSION)
if(${PYTHONLIBS_MAJOR_VERSION} STREQUAL "3")
  set(PYTHONLIBS_SWIG_FLAGS -py3)
endif()

set(INTERFACE_SRCS
  ../interface/arrays.i
  ../interface/diagnostics.i
  ../interface/exports.i
  ../interface/filenames.i
  ../interface/gasmix.i
  ../interface/geometry.i
  ../interface/indices.i
  ../interface/ventilation.i
)

include_directories(${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../c/src)

set(PYTHON_PACKAGE_DIR ${CMAKE_CURRENT_BINARY_DIR}/aether)
set(CMAKE_SWIG_OUTDIR "${PYTHON_PACKAGE_DIR}")

foreach(SWIG_INTERFACE ${INTERFACE_SRCS})
  get_filename_component(MODULE_NAME ${SWIG_INTERFACE} NAME_WE)
  swig_add_module(${MODULE_NAME} python ${SWIG_INTERFACE})
  swig_link_libraries(${MODULE_NAME} ${PYTHON_LIBRARIES} aether_c)
  set(MODULE_TARGET ${SWIG_MODULE_${MODULE_NAME}_REAL_NAME})
  set_target_properties(${MODULE_TARGET} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PYTHON_PACKAGE_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${PYTHON_PACKAGE_DIR}
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1)
endforeach()


