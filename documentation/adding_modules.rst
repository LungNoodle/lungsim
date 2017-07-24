===================
How to add a module
===================

First using the correct fortran style, begin to write your module:

.. toctree::
   :maxdepth: 1

   fortranstyles 

First we need to include it as a module for CMake. From your lungsim root directory, open CMakeLists.txt in the ./src/lib directory and insert the module into your list of source files by adding the following line in the appropriate place (source files are listed in alphabetical order, CMake will sort out the compilation order):

*module_name.f90*

Now we need to set up bindings. Create the following files, and fill in their contents appropriately (instructions to come in future)

./src/bindings/c/src/module_name.c

./src/bindings/c/src/module_name.f90

./src/bindings/c/src/module_name.h

./src/bindings/interface/module_name.i

and add interface sources named after your module in:

./src/bindings/python/CMakeLists.txt



