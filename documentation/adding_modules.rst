===================
How to add a module
===================

First using the correct fortran style, begin to write your module:

.. toctree::
   :maxdepth: 1

   fortranstyles 

First we need to include it as a module for cmake. In your lungsim root directory, open CMakeLists.txt and insert the module into your list of source files by adding the following line in the appropriate place:

*/src/module_name.f90*

Now we need to set up bindings. Create the following files, and fill in their contents appropriately (instructions to come in future)

./bindings/c/src/module_name.c

./bindings/c/src/module_name.f90

./bindings/c/src/module_name.h

./bindings/interface/module_name.i

and add interface sources named after your module in:

./bindings/python/CMakeLists.txt



