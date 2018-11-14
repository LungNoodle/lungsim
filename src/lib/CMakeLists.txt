
# The list of all source files that are required to build the library
set(LIB_SRCS
  arrays.f90
  capillaryflow.f90
  diagnostics.f90
  exports.f90
  filenames.f90
  field_utilities.f90
  gas_exchange.f90
  gasmix.f90
  geometry.f90
  growtree.f90
  indices.f90
  imports.f90
  math_utilities.f90
  mesh_utilities.f90
  other_consts.f90
  pressure_resistance_flow.f90
  solve.f90
  species_transport.f90
  ventilation.f90

)

set(AETHER_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Modules")

# Create an object library so we can use the objects in two different libraries
add_library(aether_obj OBJECT ${LIB_SRCS})
set_target_properties(aether_obj PROPERTIES
  Fortran_MODULE_DIRECTORY "${AETHER_MODULE_DIRECTORY}"
  POSITION_INDEPENDENT_CODE TRUE)

if (MSVC)
  target_compile_options(aether_obj PRIVATE /heap-arrays0)
  # Visual Studio requires a dummy file as it will not build a library from only objects.
  set(dummy_file "${CMAKE_CURRENT_BINARY_DIR}/dummy.f90")
  file(WRITE ${dummy_file} "module dummy
implicit none
contains
end module dummy\n")
endif ()

# Create the aether fortran library
add_library(aether $<TARGET_OBJECTS:aether_obj> ${dummy_file})

set_target_properties(aether PROPERTIES
  Fortran_VISIBILITY_PRESET hidden
  C_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN 1)

target_include_directories(aether PUBLIC "${AETHER_MODULE_DIRECTORY}")

install(TARGETS aether
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
