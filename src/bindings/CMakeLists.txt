
add_subdirectory(c)

find_package(SWIG QUIET)
mark_as_advanced(SWIG_EXECUTABLE)

if(SWIG_FOUND)
  option(AETHER_BUILD_PYTHON_BINDINGS "Build Python bindings for ${PROJECT_NAME}" YES)

  if(AETHER_BUILD_PYTHON_BINDINGS)
    add_project_config_parameter(PREFERRED_PYTHON_VERSION "" STRING "Set the preferred Python version for the bindings.")
    add_project_config_parameter(Python_ROOT_DIR "" PATH "Define the root directory of a Python installation.")
    set(Python_ROOT_DIR ${AETHER_Python_ROOT_DIR})

    find_package(Python ${AETHER_PREFERRED_PYTHON_VERSION} COMPONENTS Interpreter Development NumPy)
    if(Python_FOUND)
      add_subdirectory(python)
    else()
      message(WARNING "Python bindings requested but a Python interpreter, library, and numpy were not found.")
    endif()
  endif()
endif()
