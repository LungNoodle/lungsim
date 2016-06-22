# Find the Python NumPy package
# PYTHON_NUMPY_INCLUDE_DIR
# NUMPY_FOUND
# will be set by this script

cmake_minimum_required(VERSION 2.6)

if(NOT PYTHON_EXECUTABLE)
  if(NumPy_FIND_QUIETLY)
    find_package(PythonInterp QUIET)
  else()
    find_package(PythonInterp)
    set(_numpy_out 1)
  endif()
endif()

if (PYTHON_EXECUTABLE)
  # write a python script that finds the numpy path
  file(WRITE ${PROJECT_BINARY_DIR}/FindNumpyPath.py
      "try: import numpy; print(numpy.get_include())\nexcept:pass\n")

  # execute the find script
  exec_program("${PYTHON_EXECUTABLE}" ${PROJECT_BINARY_DIR}
    ARGS "FindNumpyPath.py"
    OUTPUT_VARIABLE NUMPY_PATH)
    unset (NUMPY_LIBRARIES)

  execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import sys, site; map(lambda x: sys.stdout.write('%s\\n' % x), site.getsitepackages())"
    RESULT_VARIABLE _PYTHON_SUCCESS
    OUTPUT_VARIABLE _PYTHON_SITE_PACKAGES
    ERROR_VARIABLE _PYTHON_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  set(_NUMPY_LIB_PATHS)
  if (_PYTHON_SUCCESS STREQUAL "0")
    string(REPLACE "\n" ";" _PYTHON_SITE_PACKAGES ${_PYTHON_SITE_PACKAGES})
    foreach(d ${_PYTHON_SITE_PACKAGES})
      list(APPEND _NUMPY_LIB_PATHS "${d}/numpy/core")
    endforeach()
    find_library (NUMPY_NPYMATH_LIBRARY npymath
      HINTS ${_NUMPY_LIB_PATHS}
      PATH_SUFFIXES lib)
  endif()
elseif(_numpy_out)
  message(STATUS "Python executable not found.")
endif(PYTHON_EXECUTABLE)

find_path(PYTHON_NUMPY_INCLUDE_DIR numpy/arrayobject.h
  HINTS "${NUMPY_PATH}" "${PYTHON_INCLUDE_PATH}")

if(PYTHON_NUMPY_INCLUDE_DIR)
  set(PYTHON_NUMPY_FOUND 1 CACHE INTERNAL "Python numpy found")
endif(PYTHON_NUMPY_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NumPy DEFAULT_MSG NUMPY_NPYMATH_LIBRARY PYTHON_NUMPY_INCLUDE_DIR)

if(NUMPY_FOUND)
  set(NUMPY_INCLUDE_DIRS ${PYTHON_NUMPY_INCLUDE_DIR})
  set(NUMPY_LIBRARIES ${NUMPY_NPYMATH_LIBRARY})
endif()
