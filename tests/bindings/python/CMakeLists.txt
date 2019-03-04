
# Create a Python virtual environment to use for testing.  This is only done once.
include(PythonVirtualEnv)

set(TEST_SRCS
  geometry_test.py
)

foreach(_TEST ${TEST_SRCS})
  # We'll actually take a copy of the test and run it from the build directory and not the source directory.
  configure_file(${_TEST} ${_TEST} COPYONLY)

  # Set up the test.
  get_filename_component(_TEST_NAME ${_TEST} NAME_WE)
  set(_TEST_NAME python_${_TEST_NAME})
  add_test(NAME ${_TEST_NAME} COMMAND ${VIRTUALENV_PYTHON_EXECUTABLE} ${_TEST})
  if (WIN32)
    set_tests_properties(${_TEST_NAME} PROPERTIES
      ENVIRONMENT "PATH=$<TARGET_FILE_DIR:aether>\;%PATH%;PYTHONPATH=${PROJECT_BINARY_DIR}/src/bindings/python;TEST_RESOURCES_DIR=${CMAKE_CURRENT_SOURCE_DIR}/resources")
  else ()
    set_tests_properties(${_TEST_NAME} PROPERTIES
      ENVIRONMENT "PYTHONPATH=${PROJECT_BINARY_DIR}/src/bindings/python;TEST_RESOURCES_DIR=${CMAKE_CURRENT_SOURCE_DIR}/resources")
  endif ()

endforeach()
