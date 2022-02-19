
========
Building
========

There are two ways to build the pulmonary simulation library:

#. Using CMake (recommended)
#. Using a makefile which calls the libraries and control files that you need

The CMake builds the library outside of the source directory and allows for configuration via the command line or GUI.
CMake will configure the build files required for the current environment, on GNU/Linux and OS X this could be a Makefile or Xcode project and on Windows a Visual Studio solution file.
The supplied makefile shows an example of building your library.
The library will build in the source directory with the flags set in the makefile itself.  

------------
Requirements
------------

In order to build the Aether library there are some tools that are required:

* Compiler toolchain
* CMake
* SWIG (optional)
* Python (optional)
* Sphinx (optional)

  * sphinx-fortran

If you wish to build the Python bindings for the library then Python and SWIG become necessary requirements.
Sphinx is used to generate nicely formatted output for the documentation, you can still edit and read the documentation without Sphinx.
The 'docs' target in the build will generate the html version of the documentation using Sphinx, without Sphinx this target will not be available.

Virtual environment
===================

We recommend that you create a virtual environment to use for building the library.
The following sections will assume that you have created a virtual environment and activated it.

To create a virtual environment for building the library execute the following commands, these commands assume a bash shell on GNU/Linux or macOS::

  python -m venv venv-aether
  source venv-aether/bin/activate

For Windows *cmd* execute the following::

  python -m venv venv-aether
  venv-aether\Scripts\activate

The Python bindings require *numpy* to build, we install *numpy* with::

  pip install numpy

For installing the packages required for building the documentation with Sphinx::

  pip install -r lungsim/documentation/requirements.txt

Windows
=======

On Windows, Visual Studio is the recommended toolchain with the Intel fortran compiler.
CMake is readily available and a binary is supplied on the CMake `download page <CMakeDownload_>`_.
SWIG 4.0 is available from the SWIG `download page <SWIGDownload_>`_.
The latest release at this time is version 4.0.2.
For Sphinx you will first need to have Python installed, which is required when creating the Python bindings anyway.
Python 3.9 works well with Visual Studio 2017 when building Python extension libraries (which is what the bindings are when used from Python).
Python 3.9 is availble from the Python `download page <PythonDownload>`_.

macOS
=====

Use brew to install gcc, which includes gfortran::

  brew install gcc

If you don't have brew, install it by following the instructions from `brew.sh <http://brew.sh/>`_.
CMake is readily available and a binary is supplied on the CMake `download page <CMakeDownload_>`_.
SWIG can be installed through brew::

  brew install swig
  
GNU/Linux
=========

The package manager for the distro will (most likely) have the required packages to install.
Before installing check to see if any of the requirements are already available::

  gfortran --version
  cmake --version
  python --version
  swig -version
  
In the case of the python package we require the *development* package for python this must be installed for the python bindings to become available.
For the Ubuntu distribution you can get the missing packages with the following commands::

  sudo apt-get install gfortran
  sudo apt-get install cmake
  sudo apt-get install pythonX.Y-dev # Where X and Y are the major and minor version numbers of the Python you want to install, any version above 3.7 will work
  sudo apt-get install swig

-----
CMake
-----

CMake is designed for out-of-source builds this enables us to have different builds with different configurations available from the same source.
Typically we create sibling directories of the source directory to build the application within, this is not necessary though the build directory can be anywhere.
To simply build the library we would run the following commands in the terminal (starting from the parent directory of *lungsim*)::

  cmake -S lungsim -B build-lungsim
  cd build-lungsim
  make

This will build a **Release** version of the library by default.
To build a debug version we would run the following commands::

  cmake -S lungsim -B build-lungsim -D BUILD_TYPE=Debug
  cd lungsim-build-debug
  make

Here we use the **-D** to set a configuration option, in this case *BUILD_TYPE*, to the value **Debug**.
For the library we can configure three different build types; **Release**, **Debug**, and **Pedantic**.
The **Release** build type creates an optimized application, the **Debug** build type creates an application with debugging symbols present and the **Pedantic** build type turns on more warnings and tests to help create reliable software.
The **Pedantic** option is only available with the GNU Fortran compiler at this time.

The build can also be configured with a CMake GUI application, for instance you could use the ncurses based CMake configuration application called *ccmake* to configure a build.
When configuring the build with CMake on Windows and OS X there are easily installable binaries provided for these platforms that will install a GUI.
When using the GUI you must specify the source and build directory and the type of generator to generate the build files for.  With these requirements set options for setting the build like build type become available.

Targets
=======

Below is a list of the more important targets that can be built.
Each target can be built either from the command line on make based scripts or through a project for IDE build scripts.

aether
------

The *aether* target builds the aether fortran libary.

cbindings
---------

The *cbindings* target builds the aether C library.
This target is synonymous with aether_c.

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


.. _CMakeDownload: https://cmake.org/download

.. _SWIGDownload: http://www.swig.org/download.html

.. _PythonDownload: https://www.python.org/downloads/
