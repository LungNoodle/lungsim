
project(AetherCBindings LANGUAGES C Fortran)

set(C_LIB_SRCS
    src/arrays.f90
    src/arrays.c
    src/diagnostics.c
    src/diagnostics.f90
    src/exports.c
    src/exports.f90
    src/filenames.c
    src/filenames.f90
    src/gasmix.c
    src/gasmix.f90
    src/geometry.c
    src/geometry.f90
    src/indices.c
    src/indices.f90
    src/utils.f90
    src/utils.c
#  src/solve.c
    src/ventilation.c
    src/ventilation.f90
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
  src/utils.h
  src/ventilation.h
)

add_library(aether_c SHARED ${C_LIB_SRCS} ${C_LIB_HDRS} $<TARGET_OBJECTS:aether_obj>) #  Aether to be replaced with whatever the final name of the code is
set_target_properties(aether_c PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Modules"
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1)
target_include_directories(aether_c PRIVATE $<TARGET_PROPERTY:aether_obj,Fortran_MODULE_DIRECTORY>)
#target_sources(aether_c PUBLIC aether)
#target_link_libraries(aether_c aether)
#target_link_libraries(aether_c PUBLIC aether)
