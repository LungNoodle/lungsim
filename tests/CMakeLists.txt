
add_subdirectory(pFUnit-3.2.9)

# Include the config file so we can use the test creation function directly.
include("${CMAKE_CURRENT_BINARY_DIR}/pFUnit-3.2.9/pFUnitConfig.cmake")

add_subdirectory(lib)
