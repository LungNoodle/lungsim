
========
Building
========

There are two ways to build the pulmonary simulation library:

#. Using CMake (recommended)
#. Using a makefile which calls the libraries and control files that you need

The CMake builds the library outside of the source directory and allows for configuration via the command line or GUI.  CMake will configure the build files required for the current environment, on GNU/Linux and OS X this could be a Makefile or Xcode project and on Windows a Visual Studio solution file.  The supplied makefile shows an example of building your library. The library will build in the source directory with the flags set in the makefile itself.  

-----
CMake
-----

CMake is designed for out-of-source builds this enables us to have different builds with different configurations available.  Typically we create sibling directories of the source directory to build the application within, this is not necessary though the build directory can be anywhere.  To simply build the library we would run the following commands in the terminal (starting from the parent directory of *lungsim*)::

  mkdir lungsim-build
  cd lungsim-build
  cmake ../lungsim
  make

This will build a **Release** version of the application by default.  To build a debug version we would run the following commands::

  mkdir lungsim-build-debug
  cd lungsim-build-debug
  cmake -DBUILD_TYPE=Debug ../lungsim
  make

Here we use the **-D** to set a configuration option, in this case *BUILD_TYPE*, to the value **Debug**.  For the library we can configure three different build types; **Release**, **Debug**, and **Pedantic**.  The **Release** build type creates an optimized application, the **Debug** build type creates an application with debugging symbols present and the **Pedantic** build type turns on more warnings and tests to help create reliable software.  The **Pedantic** option is only available with the GNU Fortran compiler at this time.

The build can also be configured with a CMake GUI application, for instance you could use the ncurses based CMake configuration application called *ccmake* to configure a build.  When configuring the build with CMake on Windows and OS X there are easily installable binaries provided for these platforms that will install a GUI.  When using the GUI you must specify the source and build directory and the type of generator to generate the build files for.  With these requirements set options for setting the build like build type become available.

Targets
=======

Below is a list of the more important targets that can be built.  Each target can be built either from the command line on make based scripts or through a project for IDE build scripts.

aether
------

The *aether* target builds the aether fortran libary.

cbindings
---------

The *cbindings* target builds the aether C library.  This target is synonymous with aether_c.

pybindings
----------

The *pybindings* target builds the aether Python package and associated modules.

.. note:: The *pybindings* target is only available if both Python and SWIG are available.

docs
----

The *docs* target builds the documentation from the restructured text into html which can be viewed with a webbrowser from the build directory (for example some_path/lungsim-build/html/index.html).

.. note::  This target is only available if Sphinx is available.

clean
-----

The *clean* target removes all generated files.

-----------------
Supplied makefile
-----------------

From the terminal change into the 'lungsim' directory, then run the **make** command.  Edit the compiler flags by editing the makefile in this directory.

.. note:: Not recently checked to see if this is still working.
