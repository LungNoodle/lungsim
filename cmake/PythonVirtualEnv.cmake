

set(VIRTUALENV_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/venv_for_tests"
    CACHE INTERNAL "Path to virtualenv" )
if (WIN32)
  set(PYTHON_PLATFORM_BIN_DIR Scripts)
else ()
  set(PYTHON_PLATFORM_BIN_DIR bin)
endif ()

find_program(VIRTUALENV_PYTHON_EXECUTABLE python PATHS "${VIRTUALENV_DIRECTORY}/${PYTHON_PLATFORM_BIN_DIR}" NO_DEFAULT_PATH)
mark_as_advanced(VIRTUALENV_PYTHON_EXECTUABLE)

if (NOT EXISTS "${VIRTUALENV_PYTHON_EXECTUABLE}")
  function(_create_virtualenv_from_exec call)
    execute_process(COMMAND
      ${call} --python=${PYTHON_EXECUTABLE} ${VIRTUALENV_DIRECTORY}
      RESULT_VARIABLE RESULT
      ERROR_VARIABLE ERROR
      OUTPUT_VARIABLE OUTPUT
    )
    if(NOT "${RESULT}" STREQUAL "0")
      message(STATUS "${RESULT}")
      message(STATUS "${OUTPUT}")
      message(STATUS "${ERROR}")
      message(FATAL_ERROR "Could not create virtual environment.")
    endif()
  endfunction()

  get_filename_component(_PYTHON_BIN "${PYTHON_EXECUTABLE}" PATH)
  find_program(VIRTUALENV_EXECUTABLE virtualenv)
  mark_as_advanced(VIRTUALENV_EXECUTABLE)

  # Could also check for 'venv' and 'virtualenv' Python packages as we could use these instead of
  # the virtualenv executable.
  if (VIRTUALENV_EXECUTABLE)
    message(STATUS "Creating virtual environment ...")
    _create_virtualenv_from_exec("${VIRTUALENV_EXECUTABLE}")
  else ()
    message(FATAL_ERROR "Could not find virtualenv.")
  endif ()

  find_program(VIRTUALENV_PYTHON_EXECUTABLE python PATHS "${VIRTUALENV_DIRECTORY}/${PYTHON_PLATFORM_BIN_DIR}" NO_DEFAULT_PATH)

  if (VIRTUALENV_PYTHON_EXECUTABLE)
    message(STATUS "Creating virtual environment ... success")
  else ()
    message(FATAL_ERROR "Virtual environment Python executable does not exist.")
  endif ()
endif ()
