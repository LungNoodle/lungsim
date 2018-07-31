
=======
Testing
=======

Testing is an integral part of developing software and validating the code.  It is also builds trust in users of the software that the software will work as expected.  When adding new code to the library itself a test must also be added.  Ideally the tests test every single line of code so that we have 100% code coverage.  When we have good coverage of tests over the library we inherently get regression testing.  We don't want to have new code changes breaking code that is already considered to be working.

For the testing of the Fortran code the pFUnit testing framework has been chosen.  The pFUnit testing framework uses Python to manage some of the test generation, therefore without Python we cannot build the tests.  For testing the code from Python we use the *unittest* testing module that is part of the Python distribution.

How to add a Fortran test
=========================

All tests live under the *tests* tree and mirror what is in the source tree.  In the following example we are going to add a new testing module for the *diagnostics* module in the *lib* directory from the *src* tree.

Write test
----------

To start we are first going to make sure we have the correct structure that matches the *src* tree.  Starting from the root directory of the lungsim repository we need to make sure that the directory::

   tests/lib

exists and if not create it, from the command line on UNIX based oses this can be done with the *mkdir* command (of course this should already be done, so there are big problems if you have to run this command!)::

   mkdir -p tests/lib

Once the directory structure is correct we then create the testing module.  Because we want to test the *diagnostics* module from the library we will create a test file named *test_diagnostics.pf* in the *tests/lib* directory.  The *pf* extension indicates that this file is a hybrid Python fortran file, this file is a preprocessor input file which is Fortran free format file with preprocessor directives added.  To create the test a Python script will generate a valid Fortran file from directives written into this file.  With your favourite text editor create a file named *test_diagnostics.pf*.  We could choose *vi* for this task as shown below but any text editor will work::

   vi tests/lib/test_diagnostics.pf

Into this file we will write our first test for the module.  This test will check that the diagnositcs flag has been set when using the *set_diagnostics_on* subroutine::

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

The first task to do when adding a test to the CMake files is to check that a CMake file exists.  When adding a test to a new directory, as we are doing here, there won't be a CMake file for us to use.  To fix this we first need to tell CMake that a new subdirectory is available.  We do this by adding a *sub_directory* command into an existing *CMakeLists.txt* file in a parent directory of the directory we have just added a test to.  In our example we would edit the file (any text editor will do, don't feel you need to use *vi*)::

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

With our test added to the test framework we can now build and run our test.

Build and run test
------------------

The test we have just completed will be built when we build the configuration from the build directory by default.  That is if we execute the *BUILD_ALL* build target for IDEs like Visual Studio or on *Makefile* generation builds we would simple issue the command *make* in the build directory.  We can also build our test directly be building the target *diagnostics_test*, for *Makefile* generation builds we would issue the command::

   make diagnostics_test

To run the test we can execute the ctest command from the command line in the build directory with the following arguments::

   ctest -R diagnostics_test

we will also execute all tests if we execute the command::

   ctest

A handy flag to add to both of these commands is the *--verbose* flag.  This gives us the details output from each test and not just the summary statement.


How to add a Python test
========================

In the following example we are going to add a new testing module for the *geometry* module for the *Python* bindings in the *src* tree.

Write test
----------

To start we are first going to make sure we have the correct structure that matches the *src* tree.  Starting from the root directory of the lungsim repository we need to make sure that the directory::

   tests/bindings/python

exists and if not create it, from the command line on UNIX based oses this can be done with the *mkdir* command (of course this should already be done, so there are big problems if you have to run this command!)::

   mkdir -p tests/bindings/python

Once the directory structure is correct we then create the testing module.  Because we want to test the *geometry* module from the library we will create a test file named *geometry_test.py* in the *tests/bindings/python* directory.  We could choose *vi* for this task as shown below but any text editor will work::

   vi tests/bindings/python/geometry_test.py

This file is going to be a standard Python file that makes use of the *unittest* unit testing framework.  In this file we are going to write our first test for the module.  This test will check that the *define_node_geometry_2d* method correctly sets the value of the nodes read from the *square.ipnode* file::

   import os
   import unittest

   from aether.diagnostics import set_diagnostics_on
   from aether.geometry import define_node_geometry_2d
   from aether.arrays import check_node_xyz_2d

   # Look to see if the 'TEST_RESOURCES_DIR' is set otherwise fallback to a path
   # relative to this file.
   if 'TEST_RESOURCES_DIR' in os.environ:
       resources_dir = os.environ['TEST_RESOURCES_DIR']
   else:
       here = os.path.abspath(os.path.dirname(__file__))
       resources_dir = os.path.join(here, 'resources')


   class GeometryTestCase(unittest.TestCase):

       def test_read_square(self):
           set_diagnostics_on(False)
           define_node_geometry_2d(os.path.join(resources_dir, 'square.ipnode'))
           value = check_node_xyz_2d(1, 1, 10)
           self.assertEqual(10, value)


   if __name__ == '__main__':
       unittest.main()


The first thing to note is that when we are handling external resources like files we should be explicit in where they are coming from.  This is reason for the following statement::

   if 'TEST_RESOURCES_DIR' in os.environ:
       resources_dir = os.environ['TEST_RESOURCES_DIR']
   else:
       here = os.path.abspath(os.path.dirname(__file__))
       resources_dir = os.path.join(here, 'resources')

When we are using ctest to run the tests the *TEST_RESOURCES_DIR* environment variable defines the location of the resources directory.  If this environment variable is not found then the fallback is to set the resources directory relative to the file itself.  This allows us to run the test file in different environments.  The end result is that we can explicitly state (in a relative sense) the location of the *square.ipnode* file resource for this test.

The second point to note is that our test is defined within a **class** that derives from *unittest.TestCase*.  Any class deriving from *unittest.TestCase* found in this file will be run as a test.  This means we are free to have multiple classes derived from *unittest.testcase* if we so choose.  The benefit of this is that we can group our tests by some heuristic.

The third point is on the test itself.  In this test we are testing to make sure that the correct value for the node is set when reading in a node geometry file.  We use the *unittest* framework to assert that the value of the node is indeed 10.

The final point is about the last two lines of the file.  It is these two lines that get executed by ctest if they are missing the test will actually pass this is because no tests will have been run therefor it is important that these two lines are present at the bottom of every file that defines classes derived from *unittest.TestCase*.

With our test written we now need to add this into the CMake build generation system.

Add test to CMake
-----------------

To make CMake aware of the new test we need to add an new entry in to the *TEST_SRCS* CMake variable in the file *tests/bindings/python/CMakeLists.txt*.  For this example we need to add *geometry_test.py*.  With this being our first Python test our *TEST_SRCS* variable will look like the following::

   set(TEST_SRCS
     geometry_test.py
   )

Over time we should have a list of test files defined here.

Run test
--------

The Python tests do not need building, as such, they do however require a little preprocessing for ctest to run them.  We can make sure the preprocessing is done in one of two ways.

1. Execute the test command::

   make test

2. Make CMake perform the preprocessing, and then run the tests with ctest::

   cmake .
   ctest

All of these commands must of course be executed from within the build directory.

We can of course run just some of the tests using the *-R* flag to ctest.  To run just the Python tests we could execute the following command::

   ctest -R python_

Don't forget about the verbose flag::

   ctest -R python_ -V

This command will show us more detailed output from the Python tests.
