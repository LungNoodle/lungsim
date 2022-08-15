
=======
Testing
=======

Testing is an integral part of developing software and validating the code.
It is also builds trust in users of the software that the software will work as expected.
When adding new code to the library itself a test must also be added.  Ideally the tests test every single line of code so that we have 100% code coverage.
When we have good coverage of tests over the library we inherently get regression testing.
We don't want to have new code changes breaking code that is already considered to be working.

For the testing of the Fortran code the pFUnit testing framework has been chosen.
The pFUnit testing framework uses Python to manage some of the test generation, without Python we cannot build the tests.

How to add a test
=================

All tests live under the *tests* tree and mirror what is in the source tree.
In the following example we are going to add a new testing module for the *diagnostics* module in the *lib* directory from the *src* tree.

Write test
----------

To start we are first going to make sure we have the correct structure that matches the *src* tree.
Starting from the root directory of the lungsim repository we need to make sure that the directory::

   tests/lib

exists and if not create it, from the command line on UNIX based oses this can be done with the *mkdir* command::

   mkdir tests/lib

Once the directory structure is correct we then create the testing module.
Because we want to test the diagnostics module from the library we will create a test file named *test_diagnostics.pf* in the *tests/lib* directory.
The *pf* extension indicates that this file is a hybrid Python fortran file, this file is a preprocessor input file which is Fortran free format file with preprocessor directives added.
To create the test a Python script will generate a valid Fortran file from directives written into this file.
With your favourite text editor create a file named *test_diagnostics.pf*.
We could choose *vi* for this task as shown below but any text editor will work::

   vi tests/lib/test_diagnostics.pf

Into this file we will write our first test for the module.
This test will check that the diagnositcs flag has been set when using the *set_diagnostics_on* subroutine::

   @test
   subroutine testSetDiagnostics()
      use pfunit_mod
      use diagnostics, only: get_diagnostics_on, set_diagnostics_on
      implicit none

      logical :: state

      call get_diagnostics_on(state)
      @assertFalse(state)
      call set_diagnostics_on(.true.)
      call get_diagnostics_on(state)
      @assertTrue(state)

   end subroutine testSetDiagnostics

With our test written we now need to add this into the CMake build generation system.


Add test to CMake
-----------------

The first task to do when adding a test to the CMake files is to check that a CMake file exists.
When adding a test to a new directory, as we are doing here, there won't be a CMake file for us to use.
To fix this we first need to tell CMake that a new subdirectory is available.
We do this by adding a *sub_directory* command into an existing *CMakeLists.txt* file in a parent directory of the directory we have just added a test to.
In our example we would edit the file (any text editor will do, don't feel you need to use *vi*)::

   vi tests/CMakeLists.txt

and add the line at the bottom of the file::

   add_subdirectory(lib)

Then we need to create a new *CMakeLists.txt* (the capitalisation of this file is important) file in the *tests/lib* directory (any text editor will do, don't feel you need to use *vi*)::

   vi tests/lib/CMakeLists.txt

and add the following to create an executable test that will work with CTest (we will also be able to execute this test directly)::

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
   target_link_libraries (diagnostics_test aether)
   target_include_directories(diagnostics_test PRIVATE $<TARGET_PROPERTY:aether,Fortran_MODULE_DIRECTORY>)

With our test added to the test framework we can now build and run our test.

Build and run test
------------------

The test we have just completed will be built when we build the configuration from the build directory by default.
That is if we execute the *BUILD_ALL* build target for IDEs like Visual Studio or on *Makefile* generation builds we would simple issue the command *make* in the build directory.
We can also build our test directly be building the target *diagnostics_test*, for *Makefile* generation builds we would issue the command::

   make diagnostics_test

To run the test we can execute the ctest command from the command line in the build directory with the following arguments::

   ctest -R diagnostics_test

we will also execute all tests if we execute the command::

   ctest

A handy flag to add to both of these commands is the *--verbose* flag.
This gives us the details output from each test and not just the summary statement.

