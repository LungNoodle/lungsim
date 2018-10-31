
project(AetherCBindings LANGUAGES C Fortran)

set(C_FORTRAN_LIB_SRCS
    arrays.f90
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
    pressure_resistance_flow.f90
    species_transport.f90
    utils.f90
    ventilation.f90
)
set(C_C_LIB_SRCS
    arrays.c
    diagnostics.c
    exports.c
    filenames.c
    field_utilities.c
    gas_exchange.c
    gasmix.c
    geometry.c
    growtree.c
    indices.c
    imports.c
    pressure_resistance_flow.c
    species_transport.c
    utils.c
    ventilation.c
)

set(C_LIB_HDRS
  symbol_export.h
  arrays.h
  diagnostics.h
  exports.h
  filenames.h
  field_utilities.h
  gas_exchange.h
  gasmix.h
  geometry.h
  growtree.h
  indices.h
  imports.h
  pressure_resistance_flow.h
  species_transport.h
  utils.h
  ventilation.h
)

add_library(aether_c_fortran STATIC ${C_FORTRAN_LIB_SRCS} $<TARGET_OBJECTS:aether_obj>)
set_target_properties(aether_c_fortran PROPERTIES
    POSITION_INDEPENDENT_CODE True
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Modules"
    FOLDER bindings/c)
target_include_directories(aether_c_fortran PRIVATE $<TARGET_PROPERTY:aether_obj,Fortran_MODULE_DIRECTORY>)
add_library(aether_c SHARED ${C_C_LIB_SRCS} ${C_LIB_HDRS}) #  Aether to be replaced with whatever the final name of the code is
set_target_properties(aether_c PROPERTIES
    FOLDER bindings/c
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1)
# Have to set this manually, as it is only set automatically if the file extension has an
# uppercase 'F' and not lowercase as we have here for GNU compiler.
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    set_source_files_properties(${C_FORTRAN_LIB_SRCS} PROPERTIES COMPILE_FLAGS -cpp)
endif ()
target_link_libraries(aether_c PUBLIC aether_c_fortran)

add_custom_target(cbindings ALL)
add_dependencies(cbindings aether_c)

install(TARGETS aether_c aether_c_fortran
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

