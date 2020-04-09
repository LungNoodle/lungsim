# ------------------------------------------------------------------------ #
# DOCUMENTATION:
# ------------------------------------------------------------------------ #
#
# Command line options:
#     MPI=YES                                 ! defaults to NO
#     MAX_RANK=<max array rank for generated code>
#                                             ! defaults to 5 or
#                                             !   ${PFUNIT_MAX_ARRAY_RANK}
#
# Usage:
#   cmake -DMPI=YES <path to source>
#
# ------------------------------------------------------------------------ #
cmake_minimum_required(VERSION 3.3)
cmake_policy(SET CMP0048 NEW)
project (pFUnit
    VERSION 3.1.1
    LANGUAGES Fortran)

if (DEFINED MAX_RANK)
    set(PFUNIT_MAX_RANK ${MAX_RANK} CACHE STRING "Maximum array rank for generated code." FORCE)
endif ()
set(PFUNIT_MAX_RANK 5 CACHE STRING "Maximum array rank for generated code.")
if (NOT PFUNIT_MAX_RANK)
   # Remove PFUNIT_MAX_RANK in version 4.
   if (DEFINED ENV{PFUNIT_MAX_RANK})
      set(PFUNIT_MAX_RANK "$ENV{PFUNIT_MAX_RANK}")
   # Promote following to primary in version 4.
   elseif (DEFINED ENV{PFUNIT_MAX_ARRAY_RANK})
      set(PFUNIT_MAX_RANK "$ENV{PFUNIT_MAX_ARRAY_RANK}")
   else ()
      set(PFUNIT_MAX_RANK 5)
   endif ()
endif ()

if (DEFINED MPI)
    set(PFUNIT_MPI ${MPI} CACHE BOOL "Use MPI for parallel runs" FORCE)
endif ()
option(PFUNIT_MPI "Use MPI for parallel runs" NO)
option(PFUNIT_FIND_MPI "Attempt to automatically find MPI information" YES)
option(PFUNIT_MPI_USE_MPIEXEC "Use CMake's FindMPI to find mpiexec" NO)
option(PFUNIT_OPENMP "Use OPENMP for parallel runs" NO)

# The following sets robust to "ON".
option(PFUNIT_ROBUST "Build robust test runner" YES)

# If the user specifies this on the command line (or equivalently, in the
# CMake cache), we take them at their word and do not attempt any kind of
# autodetection.
set(PFUNIT_OPENMP_FLAGS "OPENMP_FLAGS-NOTFOUND"
  CACHE STRING "OpenMP flags for the Fortran compiler.")


if (PFUNIT_MPI)
    add_definitions(-DUSE_MPI)
    if (PFUNIT_FIND_MPI)
        set(MPI_USE_MPIEXEC ${PFUNIT_MPI_USE_MPIEXEC})
        find_package(MPI REQUIRED)
        add_definitions(${MPI_Fortran_COMPILE_FLAGS})
        include_directories(${MPI_Fortran_INCLUDE_PATH})
    endif()
    message( STATUS "MPI enabled")
endif()

# 2014-0807 MLR MAX_RANK needed for code generation, not here.
#if (MAX_RANK)
#   add_definitions(-DMAX_RANK=${MAX_RANK})
#endif()

if (PFUNIT_ROBUST)
    add_definitions(-DBUILD_ROBUST)
endif()

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

    if (UNIX)
       execute_process(COMMAND ifort --version
                        OUTPUT_VARIABLE INTEL_COMPILER_VERSION)
                               string (REGEX MATCH "([0-9]+)"
                        NUM_VERSION ${INTEL_COMPILER_VERSION})
      if ("${NUM_VERSION}" STREQUAL "13")
            add_definitions(-DINTEL_13)
      endif ()
      if ("${NUM_VERSION}" STREQUAL "16")
            add_definitions(-DINTEL_16)
      endif ()
    endif ()

    #Use ifort
    #There seems to be no easy way to determine the compiler version on Windows
    #since it isn't necessary defined in the PATH. As a workaround, assume version
    #different to 13 if exectuing ifort doesn't give any results
    if (NOT "${INTEL_COMPILER_VERSION}" STREQUAL "")
      string (REGEX MATCH "([0-9]+)"
              NUM_VERSION ${INTEL_COMPILER_VERSION})
      if ("${NUM_VERSION}" STREQUAL "13")
          add_definitions(-DINTEL_13)
      endif ()
      if ("${NUM_VERSION}" STREQUAL "16")
          add_definitions(-DINTEL_16)
      endif ()
    endif ()

    add_definitions(-DIntel)
    if (WIN32)
        set(CMAKE_Fortran_FLAGS
            "-check:uninit -check:bounds -traceback -assume:realloc_lhs ${CMAKE_Fortran_FLAGS}"
            )
    else (WIN32)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -check uninit -check bounds -check pointers -check stack -traceback -assume realloc_lhs ${CMAKE_Fortran_FLAGS}"
            )
    endif (WIN32)

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)

    # Use Gfortran
    add_definitions(-DGNU)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fbounds-check")

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL NAG)

    # use nagfor
    add_definitions(-DNAG)
# Strong checking (-C=all breaks 6.0, but later releases will include it.
#    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -f2003 -w=uda -gline -fpp -C=all -fpp")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -f2003 -w=uda -gline -fpp -C=present -fpp")

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL XL)

    # Use XLF.
    add_definitions(-DIBM)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -WF,-qfpp -C")

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)

    # Use PGI.
    add_definitions(-DPGI)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -g -traceback -Mallocatable=03 -Mbounds -Mchkfpstk -Mchkstk")

else()

    message( FATAL_ERROR "Unrecognized compiler. Please use ifort, gfortran, gfortran-mp-4.8, PGI, or nagfor" )

endif()

# Fortran OpenMP support is not yet integrated into the CMake distribution.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/include/cmake/Modules/")
if (PFUNIT_OPENMP)
#    find_package(OpenMP_Fortran REQUIRED)
#    if(OPENMP_FORTRAN_FOUND)
#      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
#      set(CMAKE_Fortran_LINKER_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
#      message( STATUS "OpenMP enabled")
#    endif()
   if (PFUNIT_OPENMP_FLAGS)
     set(OpenMP_Fortran_FLAGS "${PFUNIT_OPENMP_FLAGS}")
   else()
     find_package(OpenMP_Fortran REQUIRED)
   endif()
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
   set(CMAKE_Fortran_LINKER_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
   message( STATUS "OpenMP enabled")
endif()

set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
set(CMAKE_SKIP_RPATH ON)

add_subdirectory (source)
# add_subdirectory (tests)
add_subdirectory (include)
add_subdirectory (bin)

# Packaging
include (cmake/packaging.cmake)

# Create the pFUnitConfig.cmake and pFUnitConfigVersion files
configure_file(cmake/pFUnitConfig.cmake.in
  "${PROJECT_BINARY_DIR}/pFUnitConfig.cmake" @ONLY)
configure_file(cmake/pFUnitConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/pFUnitConfigVersion.cmake" @ONLY)

# Install the pFUnitConfig.cmake and pFUnitConfigVersion.cmake
#install(FILES
#  "${PROJECT_BINARY_DIR}/pFUnitConfig.cmake"
#  "${PROJECT_BINARY_DIR}/pFUnitConfigVersion.cmake"
#  DESTINATION . )
