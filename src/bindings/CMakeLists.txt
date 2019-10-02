
add_subdirectory(c)

find_package(SWIG QUIET)
mark_as_advanced(SWIG_EXECUTABLE)

if(SWIG_FOUND)
  option(AETHER_BUILD_PYTHON_BINDINGS "Build Python bindings for ${PROJECT_NAME}" YES)

  if(AETHER_BUILD_PYTHON_BINDINGS)
    if (${CMAKE_VERSION} VERSION_LESS 3.13)
      find_package(PythonInterp)
      find_package(PythonLibs)
      if (PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND)
        find_package(NumPy REQUIRED)
        set(Python_FOUND TRUE)
        set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
        set(Python_VERSION_MAJOR ${PYTHON_VERSION_MAJOR})
        set(Python_VERSION_MAJOR ${PYTHON_VERSION_MINOR})
        if (NOT TARGET Python::Python)
          add_library(Python::Python SHARED IMPORTED)
        endif()
        set_property(TARGET Python::Python
                      PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${PYTHON_INCLUDE_DIRS}")
        if (PYTHON_LIBRARY_RELEASE AND PYTHON_LIBRARY_DEBUG)
          set_property (TARGET Python::Python PROPERTY IMPORTED_CONFIGURATIONS RELEASE DEBUG)
          set_target_properties(Python::Python
                                 PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
                                            IMPORTED_LOCATION_RELEASE "${PYTHON_LIBRARY_RELEASE}")
          set_target_properties(Python::Python
                                 PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
                                            IMPORTED_LOCATION_DEBUG "${PYTHON_LIBRARY_DEBUG}")
        else()
          set_target_properties(Python::Python
                                 PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "C"
                                            IMPORTED_LOCATION "${PYTHON_LIBRARY}")
        endif()
      endif()
      if (NOT TARGET Python::NumPy)
        add_library(Python::NumPy SHARED IMPORTED)
      endif()
      set_property(TARGET Python::NumPy
                    PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${NUMPY_INCLUDE_DIRS}")

    else()
      find_package(Python COMPONENTS Interpreter Development NumPy)
    endif()
    if(Python_FOUND)
      add_subdirectory(python)
    else()
      message(WARNING "Python bindings requested but a Python interpreter, library, and numpy were not found.")
    endif()
  endif()
endif()
