
add_subdirectory(c)

find_package(SWIG QUIET)
mark_as_advanced(SWIG_EXECUTABLE)

if(SWIG_FOUND)
  option(AETHER_BUILD_PYTHON_BINDINGS "Build Python bindings for ${PROJECT_NAME}" YES)

  if(AETHER_BUILD_PYTHON_BINDINGS)
    find_package(Python COMPONENTS Interpreter Development NumPy)
    if(Python_FOUND)
      add_subdirectory(python)
    else()
      message(WARNING "Python bindings requested but a Python interpreter, library, and numpy were not found.")
    endif()
  endif()
endif()
