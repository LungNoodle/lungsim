
project(AetherCBindings LANGUAGES C Fortran)

set(C_FORTRAN_LIB_SRCS
    src/arrays.f90
    src/diagnostics.f90
    src/exports.f90
    src/filenames.f90
    src/gasmix.f90
    src/geometry.f90
    src/indices.f90
    src/pressure_resistance_flow.f90
    src/utils.f90
    src/ventilation.f90
)
set(C_C_LIB_SRCS
    src/arrays.c
    src/diagnostics.c
    src/exports.c
    src/filenames.c
    src/gasmix.c
    src/geometry.c
    src/indices.c
    src/pressure_resistance_flow.c
    src/utils.c
    src/ventilation.c
)

set(C_LIB_HDRS
  src/symbol_export.h
  src/arrays.h
  src/diagnostics.h
  src/exports.h
  src/filenames.h
  src/gasmix.h
  src/geometry.h
  src/indices.h
  src/pressure_resistance_flow.h
  src/utils.h
  src/ventilation.h
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

