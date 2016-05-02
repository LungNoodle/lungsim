
add_subdirectory(c)

find_package(SWIG QUIET)
mark_as_advanced(SWIG_EXECUTABLE)

if(SWIG_FOUND)
  option(AETHER_BUILD_PYTHON_BINDINGS "Build Python bindings for ${PROJECT_NAME}" YES)


  if(AETHER_BUILD_PYTHON_BINDINGS)
    find_package(PythonInterp)
    find_package(PythonLibs)
    if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
      add_subdirectory(python)
    else()
      message(WARNING "Python bindings requested but a Python interpreter and library were not found.")
    endif()
  endif()
endif()
