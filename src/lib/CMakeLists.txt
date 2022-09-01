
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
  precision.f90
  pressure_resistance_flow.f90
  solve.f90
  species_transport.f90
  surface_fitting.f90
  ventilation.f90
  wave_transmission.f90
)
if (SUPERLU_FOUND)
  set_source_files_properties(capillaryflow.f90 PROPERTIES COMPILE_DEFINITIONS HAVE_SUPERLU)
  list(APPEND LIB_SRCS c_fortran_dgssv.c)
  if (SUPERLU_HAS_GLOBALLU_T)
    set_source_files_properties(c_fortran_dgssv.c PROPERTIES COMPILE_DEFINITIONS SUPERLU_HAS_GLOBALLU_T)
  endif()
endif ()

if (NOT CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  set_source_files_properties(capillaryflow.f90 PROPERTIES COMPILE_FLAGS -cpp)
endif ()

set(AETHER_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Modules")

# Create an object library so we can use the objects in two different libraries
add_library(aether STATIC ${LIB_SRCS})
set_target_properties(aether PROPERTIES
  Fortran_MODULE_DIRECTORY "${AETHER_MODULE_DIRECTORY}"
  POSITION_INDEPENDENT_CODE TRUE)
if (MSVC)
  target_compile_options(aether PRIVATE /heap-arrays0)
endif ()

if (SUPERLU_FOUND)
  target_link_libraries(aether PUBLIC superlu)
endif ()
