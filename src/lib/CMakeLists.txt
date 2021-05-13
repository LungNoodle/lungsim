
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
endif ()
set_source_files_properties(capillaryflow.f90 PROPERTIES COMPILE_FLAGS -cpp)

set(AETHER_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Modules")

# Create an object library so we can use the objects in two different libraries
add_library(aether_obj OBJECT ${LIB_SRCS})
set_target_properties(aether_obj PROPERTIES
  Fortran_MODULE_DIRECTORY "${AETHER_MODULE_DIRECTORY}"
  POSITION_INDEPENDENT_CODE TRUE)

if (MSVC)
  target_compile_options(aether_obj PRIVATE /heap-arrays0)
endif ()

if (SUPERLU_FOUND)
  set(SUPERLU_FORTRAN_SRCS
    c_fortran_dgssv.c
    #c_fortran_cgssv.c
  )
  if (SUPERLU_HAS_GLOBALLU_T)
    set_source_files_properties(c_fortran_dgssv.c PROPERTIES COMPILE_DEFINITIONS SUPERLU_HAS_GLOBALLU_T)
  endif()
  # Create an object library for the SuperLU C Fortran cgssv file.
  add_library(aether_superlu_obj OBJECT ${SUPERLU_FORTRAN_SRCS})
  # Our object library only needs the include directory information carried by the 'superlu' target
  if (CMAKE_VERSION VERSION_LESS "3.12")
    target_include_directories(aether_superlu_obj PRIVATE $<TARGET_PROPERTY:superlu,INTERFACE_INCLUDE_DIRECTORIES>)
  else ()
    target_link_libraries(aether_superlu_obj PRIVATE superlu)
  endif ()
  set_property(TARGET aether_superlu_obj PROPERTY POSITION_INDEPENDENT_CODE ON)
endif ()

# Create the aether library from the object libraries.
add_library(aether $<TARGET_OBJECTS:aether_obj> ${dummy_file})
if (SUPERLU_FOUND)
  target_sources(aether PRIVATE $<TARGET_OBJECTS:aether_superlu_obj>)
  # The aether library actually needs to link against the SuperLU libraries.
  target_link_libraries(aether PUBLIC superlu)
endif ()

set_target_properties(aether PROPERTIES
  Fortran_VISIBILITY_PRESET hidden
  C_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN 1)

target_include_directories(aether PUBLIC "${AETHER_MODULE_DIRECTORY}")

install(TARGETS aether
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
