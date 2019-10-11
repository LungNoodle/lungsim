
function(TIDY_GUI)
    # Hide some variables in the GUI
    mark_as_advanced(CMAKE_CODEBLOCKS_EXECUTABLE)
    mark_as_advanced(CMAKE_CODEBLOCKS_EXECUTABLE)
    mark_as_advanced(QT_QMAKE_EXECUTABLE)
    if (APPLE)
        mark_as_advanced(CMAKE_OSX_ARCHITECTURES)
        mark_as_advanced(CMAKE_OSX_DEPLOYMENT_TARGET)
        mark_as_advanced(CMAKE_OSX_SYSROOT)
    endif ()
endfunction()

# Appends a library to the list of interface_link_libraries
function(append_link_library TARGET LIB)
    get_target_property(CURRENT_ILL
        ${TARGET} INTERFACE_LINK_LIBRARIES)
    if (NOT CURRENT_ILL)
        set(CURRENT_ILL )
    endif()
    # Treat framework references different
    if(APPLE AND ${LIB} MATCHES ".framework$")
        string(REGEX REPLACE ".*/([A-Za-z0-9.]+).framework$" "\\1" FW_NAME ${LIB})
        #message(STATUS "Matched '${FW_NAME}' to ${LIB}")
        set(LIB "-framework ${FW_NAME}")
    endif()
    set_target_properties(${TARGET} PROPERTIES
        INTERFACE_LINK_LIBRARIES "${CURRENT_ILL};${LIB}")
endfunction()

function(add_project_config_parameter _VARIABLE_NAME _DEFAULT_VALUE _VARIABLE_TYPE _VARIABLE_DOCS)
  string(TOUPPER ${PROJECT_NAME} _UPPER_PREFIX)
  if (DEFINED ${_VARIABLE_NAME})
    set(FORCE_VARIABLE FORCE)
  elseif (NOT ${_VARIABLE_NAME})
    set(${_VARIABLE_NAME} ${_DEFAULT_VALUE})
  endif()
  set(${_UPPER_PREFIX}_${_VARIABLE_NAME} ${${_VARIABLE_NAME}} CACHE ${_VARIABLE_TYPE} "${_VARIABLE_DOCS}" ${FORCE_VARIABLE})
  unset(${_VARIABLE_NAME} CACHE)
endfunction()

