
find_package(Python COMPONENTS Interpreter REQUIRED) #Sets ${Python_EXECUTABLE}

set(srcs Assert.F90)

set(GENERATED_ASSERT_FILES_DIR "${CMAKE_CURRENT_BINARY_DIR}/tmp_assert_files_dir")
if (NOT EXISTS "${GENERATED_ASSERT_FILES_DIR}")
    make_directory("${GENERATED_ASSERT_FILES_DIR}")
endif ()
list(APPEND _ADDTIONAL_SRCS_TO_CLEAN ${GENERATED_ASSERT_FILES_DIR})

# Get list of source files
execute_process(
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GenerateAssertsOnArrays.py --maxRank ${PFUNIT_MAX_RANK}
    WORKING_DIRECTORY "${GENERATED_ASSERT_FILES_DIR}"
    OUTPUT_VARIABLE generate_asserts_on_array_output
)
string(REPLACE "\n" ";" GENERATED_ASSERT_FILES ${generate_asserts_on_array_output})

# Generate AssertArray files:  AssertArrays.fh, generated.inc, and AssertXYZ?.F90
add_custom_command(OUTPUT ${GENERATED_ASSERT_FILES} AssertArrays.fh generated.inc
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GenerateAssertsOnArrays.py --maxRank ${PFUNIT_MAX_RANK} --quiet True
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Generating assert files."
)
add_custom_target(generate_assert_files
    DEPENDS ${GENERATED_ASSERT_FILES}
)

list(APPEND srcs "${GENERATED_ASSERT_FILES}")
#list(APPEND _ADDTIONAL_SRCS_TO_CLEAN ${GENERATED_ASSERTS})

list(APPEND srcs AssertBasic.F90)

list(APPEND srcs AbstractTestResult.F90)
list(APPEND srcs BaseTestRunner.F90)
list(APPEND srcs DebugListener.F90)

list(APPEND srcs DynamicTestCase.F90)
list(APPEND srcs Exception.F90)
list(APPEND srcs Expectation.F90)
list(APPEND srcs MockCall.F90)
list(APPEND srcs MockRepository.F90)

list(APPEND srcs ParallelContext.F90)
list(APPEND srcs ParallelException.F90)
list(APPEND srcs AbstractTestParameter.F90)
list(APPEND srcs ParameterizedTestCase.F90)
list(APPEND srcs Params.F90)

list(APPEND srcs XmlPrinter.F90)
list(APPEND srcs ResultPrinter.F90)

if ( PFUNIT_ROBUST )
  list(APPEND srcs UnixPipeInterfaces.F90)
  list(APPEND srcs UnixProcess.F90)
  list(APPEND srcs RobustRunner.F90)
  list(APPEND srcs RemoteProxyTestCase.F90)
endif()

list(APPEND srcs SerialContext.F90)
list(APPEND srcs SourceLocation.F90)
list(APPEND srcs StringConversionUtilities.F90)
list(APPEND srcs SubsetRunner.F90)
list(APPEND srcs SurrogateTestCase.F90)
list(APPEND srcs Test.F90)
list(APPEND srcs TestCase.F90)
list(APPEND srcs TestFailure.F90)
list(APPEND srcs TestListener.F90)
list(APPEND srcs TestMethod.F90)
list(APPEND srcs TestResult.F90)
list(APPEND srcs TestRunner.F90)
list(APPEND srcs TestSuite.F90)
list(APPEND srcs ThrowFundamentalTypes.F90)
list(APPEND srcs pFUnit.F90)
list(APPEND srcs pFUnitPackage.F90)

if (MPI)
  list(APPEND srcs MpiContext.F90)
  list(APPEND srcs MpiStubs.F90)
  list(APPEND srcs MpiTestParameter.F90)
  list(APPEND srcs MpiTestCase.F90)
  list(APPEND srcs MpiTestMethod.F90)
endif()

list(APPEND pysrcs CodeUtilities.py)
list(APPEND pysrcs GenerateAssertsOnArrays.py)
list(APPEND pysrcs Utilities.py)

set(SOURCES ${srcs})

add_library(pfunit STATIC ${SOURCES})
add_dependencies(pfunit generate_assert_files)

set(MODULES_DIR "${CMAKE_CURRENT_BINARY_DIR}/Modules")
set_target_properties(pfunit PROPERTIES
    Fortran_MODULE_DIRECTORY "${MODULES_DIR}")
target_include_directories(pfunit PUBLIC
    "${CMAKE_CURRENT_BINARY_DIR}"
    "${MODULES_DIR}")

#install(TARGETS pfunit DESTINATION lib)
#install(DIRECTORY "${MODULES_DIR}" DESTINATION .)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${_ADDTIONAL_SRCS_TO_CLEAN}")
