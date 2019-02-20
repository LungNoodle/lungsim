
project(AetherPythonBindings LANGUAGES C)

find_package(NumPy REQUIRED)

set(ADDITIONAL_GENERATED_FILES)
# SWIG has already been found at this point.
include(${SWIG_USE_FILE})

string(REPLACE "." ";" PYTHONLIBS_VERSION_LIST ${PYTHONLIBS_VERSION_STRING})
list(GET PYTHONLIBS_VERSION_LIST 0 PYTHONLIBS_MAJOR_VERSION)
list(GET PYTHONLIBS_VERSION_LIST 1 PYTHONLIBS_MINOR_VERSION)
if(${PYTHONLIBS_MAJOR_VERSION} STREQUAL "3")
  set(PYTHONLIBS_SWIG_FLAGS -py3 -relativeimport)
endif()

if(WIN32)
  set(SETUP_PY_OPERATING_SYSTEM_CLASSIFIER "Operating System :: Microsoft :: Windows")
elseif(DARWIN)
  set(SETUP_PY_OPERATING_SYSTEM_CLASSIFIER "Operating System :: MacOS :: MacOS X")
else()
  set(SETUP_PY_OPERATING_SYSTEM_CLASSIFIER "Operating System :: Unix")
endif()

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
  ../interface/pressure_resistance_flow.i
  ../interface/species_transport.i
  ../interface/ventilation.i
)
include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../c)

if(MSVC)
  set(PYTHON_PACKAGE_DIR ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/aether)
  set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/aether)
  set(PYTHON_PACKAGE_INIT_FILE ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/aether/__init__.py)
  set(PYTHON_SETUP_PY_FILE_DIR ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR})
  set(PYTHON_SETUP_PY_FILE_CONFIGURATION_DIR ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>)
  set(MSVC_INTDIR /$<CONFIG>)
else()
  set(PYTHON_PACKAGE_DIR ${CMAKE_CURRENT_BINARY_DIR}/aether)
  set(CMAKE_SWIG_OUTDIR "${PYTHON_PACKAGE_DIR}")
  set(PYTHON_PACKAGE_INIT_FILE ${CMAKE_CURRENT_BINARY_DIR}/aether/__init__.py)
  set(PYTHON_SETUP_PY_FILE_DIR ${CMAKE_CURRENT_BINARY_DIR})
  set(PYTHON_SETUP_PY_FILE_CONFIGURATION_DIR ${PYTHON_SETUP_PY_FILE_DIR})
endif()
set(PYTHON_MODULE_INSTALL_PREFIX lib/python${PYTHONLIBS_MAJOR_VERSION}.${PYTHONLIBS_MINOR_VERSION}${MSVC_INTDIR}/aether)

list(APPEND ADDITIONAL_GENERATED_FILES ${PYTHON_PACKAGE_INIT_FILE})

set(CMAKE_SWIG_FLAGS ${PYTHONLIBS_SWIG_FLAGS})

set(MODULE_TARGETS)
set(PYTHON_MODULE_FILES)
foreach(SWIG_INTERFACE ${INTERFACE_SRCS})
  get_filename_component(MODULE_NAME ${SWIG_INTERFACE} NAME_WE)
  if (${CMAKE_VERSION} VERSION_LESS 3.8)
    swig_add_module(${MODULE_NAME} python ${SWIG_INTERFACE})
  else ()
    swig_add_library(${MODULE_NAME}
      LANGUAGE python
      SOURCES ${SWIG_INTERFACE})
  endif ()
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
  get_target_property(_SOURCES_LIST ${MODULE_TARGET} SOURCES)
  foreach(_SOURCE in ${_SOURCES_LIST})
    string(REGEX MATCH "^.+\.py$" _MATCH ${_SOURCE})
    if(_MATCH)
      string(REPLACE "$(Configuration)" "$<CONFIG>" _MATCH ${_MATCH})
      list(APPEND PYTHON_MODULE_FILES ${_MATCH})
    endif()
  endforeach()
endforeach()

# Alternative option is to configure the setup.py file to the appropriate location.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${PYTHON_SETUP_PY_FILE_DIR}/setup.py @ONLY)

if (CMAKE_CONFIGURATION_TYPES)
  file(GENERATE OUTPUT ${PYTHON_SETUP_PY_FILE_CONFIGURATION_DIR}/setup.py
    INPUT ${PYTHON_SETUP_PY_FILE_DIR}/setup.py)
endif ()

add_custom_target(pybindings ALL)
add_dependencies(pybindings ${MODULE_TARGETS})
# Create the Python __init__.py file.
add_custom_command(TARGET pybindings POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E touch ${PYTHON_PACKAGE_INIT_FILE})
if (MSVC)
  add_custom_command(TARGET ${MODULE_TARGETS} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:aether_c> $<$<CONFIG:Debug>:$<TARGET_PDB_FILE:aether_c>> ${PYTHON_PACKAGE_DIR})
endif ()

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${ADDITIONAL_GENERATED_FILES}")

install(TARGETS ${MODULE_TARGETS}
  RUNTIME DESTINATION ${PYTHON_MODULE_INSTALL_PREFIX}
  LIBRARY DESTINATION ${PYTHON_MODULE_INSTALL_PREFIX})

install(FILES ${PYTHON_PACKAGE_INIT_FILE} ${PYTHON_MODULE_FILES} $<TARGET_FILE:aether_c> 
  DESTINATION ${PYTHON_MODULE_INSTALL_PREFIX})
if (MSVC)
  install(FILES $<$<CONFIG:Debug>:$<TARGET_PDB_FILE:aether_c>> DESTINATION ${PYTHON_MODULE_INSTALL_PREFIX})
endif ()

install(FILES ${PYTHON_SETUP_PY_FILE_CONFIGURATION_DIR}/setup.py ${Aether_SOURCE_DIR}/LICENCE.txt
  ${Aether_SOURCE_DIR}/README.rst
  DESTINATION ${PYTHON_MODULE_INSTALL_PREFIX}/..)
