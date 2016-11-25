
function(TIDY_GUI)
    # Hide some variables in the GUI
    mark_as_advanced(CMAKE_CODEBLOCKS_EXECUTABLE)
    mark_as_advanced(CMAKE_CODEBLOCKS_EXECUTABLE)
    if (APPLE)
        mark_as_advanced(CMAKE_OSX_ARCHITECTURES)
        mark_as_advanced(CMAKE_OSX_DEPLOYMENT_TARGET)
        mark_as_advanced(CMAKE_OSX_SYSROOT)
    endif ()
endfunction()
