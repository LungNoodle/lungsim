===================
How to add a module
===================

First using the correct fortran style, begin to write your module:

.. toctree::
   :maxdepth: 1

   fortranstyles 

Next we need to include it as a module for CMake. From your lungsim root directory, open CMakeLists.txt in the ./src/lib directory and insert the module into the list of source files by adding the following line in the appropriate place (source files are listed in alphabetical order under the *LIB_SRCS* variable, CMake will sort out the compilation order)::

    module_name.f90

Now we need to set up bindings. Create the following files, and fill in their contents appropriately (the details of which can be found in :doc:`wrapping`)::

    ./src/bindings/c/module_name.c

    ./src/bindings/c/module_name.f90

    ./src/bindings/c/module_name.h

    ./src/bindings/interface/module_name.i

and add the c bindings sources named after your module into the appropriate CMake variable (one of *C_FORTRAN_LIB_SRCS*, *C_C_LIB_SRCS*, *C_LIB_HDRS*) in::

    ./src/bindings/c/CMakeLists.txt

and lastly add the interface source named after your module into the CMake variable *INTERFACE_SRCS* in::

    ./src/bindings/python/CMakeLists.txt
