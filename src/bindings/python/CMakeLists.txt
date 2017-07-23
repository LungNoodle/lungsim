
project(AetherPythonBindings LANGUAGES C)

find_package(NumPy REQUIRED)

set(ADDITIONAL_GENERATED_FILES)
# SWIG has already been found at this point.
include(${SWIG_USE_FILE})

string(SUBSTRING ${PYTHONLIBS_VERSION_STRING} 0 1 PYTHONLIBS_MAJOR_VERSION)
if(${PYTHONLIBS_MAJOR_VERSION} STREQUAL "3")
  set(PYTHONLIBS_SWIG_FLAGS -py3 -relativeimport)
endif()

set(INTERFACE_SRCS
  ../interface/arrays.i
  ../interface/diagnostics.i
  ../interface/exports.i
  ../interface/filenames.i
  ../interface/gasmix.i
  ../interface/geometry.i
  ../interface/indices.i
  ../interface/pressure_resistance_flow.i
  ../interface/tree_wave_propagation.i
  ../interface/ventilation.i
)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../c/src)

if(MSVC)
  set(PYTHON_PACKAGE_DIR ${CMAKE_CURRENT_BINARY_DIR}/$(Configuration)/aether)
  set(PYTHON_PACKAGE_INIT_FILE ${CMAKE_CURRENT_BINARY_DIR}/$(Configuration)/aether/__init__.py)
  set(PYTHON_SETUP_PY_FILE_DIR ${CMAKE_CURRENT_BINARY_DIR}/$(Configuration))
else()
  set(PYTHON_PACKAGE_DIR ${CMAKE_CURRENT_BINARY_DIR}/aether)
  set(PYTHON_PACKAGE_INIT_FILE ${CMAKE_CURRENT_BINARY_DIR}/aether/__init__.py)
  set(PYTHON_SETUP_PY_FILE_DIR ${CMAKE_CURRENT_BINARY_DIR})
endif()
list(APPEND ADDITIONAL_GENERATED_FILES ${PYTHON_PACKAGE_INIT_FILE})
set(CMAKE_SWIG_OUTDIR "${PYTHON_PACKAGE_DIR}")

set(CMAKE_SWIG_FLAGS ${PYTHONLIBS_SWIG_FLAGS})

set(MODULE_TARGETS)
foreach(SWIG_INTERFACE ${INTERFACE_SRCS})
  get_filename_component(MODULE_NAME ${SWIG_INTERFACE} NAME_WE)
  swig_add_module(${MODULE_NAME} python ${SWIG_INTERFACE})
  # The Python debug libraries on windows also require the Python release libraries
  # (this appears to be the case for the libraries built with Visual Studio 14 2015),
  # so we add the release libraries here for all platforms.  This probably means that
  # when we try and link to the debug version of the Python libraires on OSX and GNU/Linux 
  # we will not succeed.  But as this is rarely done it hopefully won't become an issue.
  swig_link_libraries(${MODULE_NAME} ${PYTHON_LIBRARY_RELEASE} aether_c)
  set(MODULE_TARGET ${SWIG_MODULE_${MODULE_NAME}_REAL_NAME})
  set_target_properties(${MODULE_TARGET} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PYTHON_PACKAGE_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${PYTHON_PACKAGE_DIR}
    FOLDER bindings/python
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1)
  list(APPEND MODULE_TARGETS ${MODULE_TARGET})

  if (MSVC)
    if (${PYTHON_LIBRARY_DEBUG})
      # Again from the above comment we need the debug Python library and release Python library when
      # linking with link.exe on Windows.
      target_link_libraries(${MODULE_TARGET} $<$<CONFIG:Debug>:PYTHON_LIBRARY_DEBUG>)
    else()
      # Disable use of Debug Python libraries when not available and using MSVC.
      # See http://stackoverflow.com/questions/11311877/creating-a-dll-from-a-wrapped-cpp-file-with-swig
      target_compile_definitions(${MODULE_TARGET} PRIVATE $<$<CONFIG:Debug>:SWIG_PYTHON_INTERPRETER_NO_DEBUG>)
    endif()
    set_target_properties(${MODULE_TARGET} PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PYTHON_PACKAGE_DIR}
      LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PYTHON_PACKAGE_DIR})
  endif()
endforeach()

# Copy over setup.py
add_custom_command(
  OUTPUT ${PYTHON_SETUP_PY_FILE_DIR}/setup.py
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${PYTHON_SETUP_PY_FILE_DIR}/setup.py
  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)  

set_property(SOURCE ../interface/arrays.i APPEND PROPERTY OBJECT_DEPENDS ${PYTHON_SETUP_PY_FILE_DIR}/setup.py)
# Alternative option is to configure the setup.py file to the appropriate location.
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py 

add_custom_target(pybindings ALL)
add_dependencies(pybindings ${MODULE_TARGETS})
# Create the Python __init__.py file.
add_custom_command( TARGET pybindings POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E touch ${PYTHON_PACKAGE_INIT_FILE})

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${ADDITIONAL_GENERATED_FILES}")
