
# Add all the files that make a single test, we could have multiple files testing
# the same module.  Don't add test files into the same test that test different modules.
# These are all .pf files.
set(DIAGNOSTICS_TEST_SRCS
    test_diagnostics.pf)

# Make use of the pFUnit helper function to create a test.
# Arguments    : - test_package_name: Name of the test package
#                - test_sources     : List of pf-files to be compiled
#                - extra_sources    : List of extra Fortran source code used for testing (if none, input empty string "")
#                - extra_sources_c  : List of extra C/C++ source code used for testing (if none, input empty string "")
add_pfunit_test(diagnostics_test ${DIAGNOSTICS_TEST_SRCS} "" "")
# Link the test to the aether library target.
target_link_libraries(diagnostics_test aether)
target_include_directories(diagnostics_test PRIVATE $<TARGET_PROPERTY:aether,Fortran_MODULE_DIRECTORY>)
