
Building
========

There are two ways to build executables supplied with this application:

#. Using a makefile which calls the libraries and control files that you need
#. Using CMake

The supplied makefile shows an example of building your executable, this will use the libraries that are required to run the ventilation model, and an exposure of the ventilation model. The executable will build in the source directory with the flags set in the makefile itself.  The CMake alternative builds the application outside or the source directory and allows for configuration via the command line or GUI.

Supplied makefile
-----------------

From the terminal change into the 'NEONcode' directory, then run the **make** command.  Edit the compiler flags by editing the makefile in this directory.

CMake
-----

CMake is designed for out-of-source builds this enables us to have different builds with different configurations available.  Typically we create sibling directories of the source directory to build the application within, this is not necessary though the build directory can be anywhere.  To simply build the application we would run the following commands in the terminal (starting from the parent directory of *NEONcode*)::

  mkdir NEONcode-build
  cd NEONcode-build
  cmake ../NEONcode
  make

This will build a **Release** version of the application by default.  To build a debug version we would run the following commands::

  mkdir NEONcode-build-debug
  cd NEONcode-build-debug
  cmake -DBUILD_TYPE=Debug ../NEONcode
  make

Here we use the **-D** to set a configuration option, in this case *BUILD_TYPE*, to the value **Debug**.  For the application we can configure three different build types; **Release**, **Debug**, and **Pedantic**.  The **Release** build type creates an optimized application, the **Debug** build type creates an application with debugging symbols present and the **Pedantic** build type turns on more warnings and tests to help create reliable software.  The **Pedantic** option is only available with the GNU Fortran compiler at this time.

The build can also be configured with a CMake GUI application, for instance you could use the ncurses based CMake configuration application called *ccmake* to configure a build.  

