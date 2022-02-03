
include(FindBLAS.kitware)

if (NOT TARGET blas)

  message(STATUS "Creating target 'blas' ...")
  set(LIBS ${BLAS_LIBRARIES})

  set(INCS )
  foreach(DIRSUFF _INCLUDE_DIRS _INCLUDES _INCLUDE_PATH _INCLUDE_DIR)
    if (DEFINED BLAS${DIRSUFF})
     list(APPEND INCS ${BLAS${DIRSUFF}})
    endif()
  endforeach()

  if (NOT CMAKE_CFG_INTDIR STREQUAL .)
    string(TOUPPER ${CMAKE_CFG_INTDIR} CURRENT_BUILD_TYPE)
  elseif(CMAKE_BUILD_TYPE)
    string(TOUPPER ${CMAKE_BUILD_TYPE} CURRENT_BUILD_TYPE)
  else()
    set(CURRENT_BUILD_TYPE NOCONFIG)
  endif()

  # Creating target 'blas'
  list(GET LIBS 0 _FIRST_LIB)
  add_library(blas UNKNOWN IMPORTED)

  # Treat apple frameworks separate
  # See http://stackoverflow.com/questions/12547624/cant-link-macos-frameworks-with-cmake
  if(APPLE AND ${_FIRST_LIB} MATCHES ".framework$")
    string(REGEX REPLACE ".*/([A-Za-z0-9.]+).framework$" "\\1" FW_NAME ${_FIRST_LIB})
    #message(STATUS "Matched '${FW_NAME}' to ${LIB}")
    set(_FIRST_LIB "${_FIRST_LIB}/${FW_NAME}")
  endif()
  set_target_properties(blas PROPERTIES
    IMPORTED_LOCATION_${CURRENT_BUILD_TYPE} ${_FIRST_LIB}
    IMPORTED_CONFIGURATIONS ${CURRENT_BUILD_TYPE}
    INTERFACE_INCLUDE_DIRECTORIES "${INCS}"
    INTERFAVE_COMPILE_DEFINITIONS "${BLAS_DEFINITIONS}"
  )

  list(REMOVE_AT LIBS 0)
  # Add non-matched libraries as link libraries so nothing gets forgotten
  foreach(LIB ${LIBS})
    append_link_library(blas ${LIB})
  endforeach()
  message(STATUS "Creating target 'blas' ... done.")

endif()
