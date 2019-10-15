
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # Add new build type PEDANTIC
  SET(CMAKE_Fortran_FLAGS_PEDANTIC
    "-g3 -Wall -Wextra -Wconversion -pedantic -fbounds-check -ffpe-trap=zero,overflow,underflow -Wimplicit-procedure"
    CACHE STRING "Flags used by the Fortran compiler during pedantic builds."
    FORCE )
  SET(CMAKE_EXE_LINKER_FLAGS_PEDANTIC
    ""
    CACHE STRING "Flags used for linking binaries during pedantic builds."
    FORCE )
  SET(CMAKE_SHARED_LINKER_FLAGS_PEDANTIC
    ""
    CACHE STRING "Flags used by the shared libraries linker during pedantic builds."
    FORCE )
  MARK_AS_ADVANCED(
    CMAKE_Fortran_FLAGS_PEDANTIC
    CMAKE_EXE_LINKER_FLAGS_PEDANTIC
    CMAKE_SHARED_LINKER_FLAGS_PEDANTIC)
endif()

