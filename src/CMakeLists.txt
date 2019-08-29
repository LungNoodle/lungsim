
# Our version for FindSuperLU creates an imported target for us to use.
find_package(SuperLU)

add_subdirectory(lib)

if(AETHER_BUILD_BINDINGS)
  add_subdirectory(bindings)
endif()

