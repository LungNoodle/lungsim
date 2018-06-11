
========
Building
========

There are two ways to build the pulmonary simulation library:

#. Using CMake (recommended)
#. Using a makefile which calls the libraries and control files that you need

The CMake builds the library outside of the source directory and allows for configuration via the command line or GUI.  CMake will configure the build files required for the current environment, on GNU/Linux and OS X this could be a Makefile or Xcode project and on Windows a Visual Studio solution file.  The supplied makefile shows an example of building your library. The library will build in the source directory with the flags set in the makefile itself.  

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

If you wish to build the Python bindings for the library then Python and SWIG become necessary requirements.  Sphinx is used to generate nicely formatted output for the documentation, you can still edit and read the documentation without Sphinx.  The 'docs' target in the build will generate the html version of the documentation using Sphinx, without Sphinx this target will not be available.

Windows
=======

On Windows, Visual Studio is the recommended toolchain with intel fortran compiler.  CMake is readily available and a binary is supplied on the CMake `download page <CMakeDownload_>`_.  SWIG 3.0 is available from the SWIG `download page <SWIGDownload_>`_.  The latest release at this time is version 3.0.8.  For Sphinx you will first need to have Python installed, which is required when creating the Python bindings anyway.  Python 3.5 works well with Visual Studio 2015 when building Python extension libraries (which is what the bindings are when used from Python).  Python 3.5 is availble from the Python `download page <PythonDownload>`_.  Once Python is installed Sphinx can be installed with pip::

  pip install sphinx
  pip install sphinx-fortran
  
Also install the sphinx-fortran extension for Sphinx for dealing with Fortran code.

It is also useful to make use of python virtual environments, this is not critical on windows but very useful if you want to have easy access to multiple branches of compiled code through python.  Install virtualenv with pip::

  pip install virtualenv

OS X
====

Use brew to install gcc, which includes gfortran.  If you don't have brew install it by following the instructions from `brew.sh <http://brew.sh/>`_.  CMake is readily available and a binary is supplied on the CMake `download page <CMakeDownload_>`_.  SWIG can be installed through brew::

  brew install swig
  
To install Sphinx you will need pip, if you don't already have it you can use the `get-pip.py script <https://bootstrap.pypa.io/get-pip.py>`_::

  curl https://bootstrap.pypa.io/get-pip.py > get-pip.py
  sudo -H python get-pip.py
  
Once pip is installed you can get Sphinx and sphinx-fortran as above::

  sudo -H pip install sphinx
  sudo -H pip install sphinx-fortran
  
It is also useful to make use of python virtual environments, this is very useful if you want to have easy access to multiple branches of compiled code through python.  Install virtualenv with pip::
  
  sudo -H pip install virtualenv

GNU/Linux
=========

The package manager for the distro will (most likely) have the required packages to install.  Before installing check to see if any of the requirements are already available::

  gfortran --version
  cmake --version
  python --version
  swig -version
  sphinx-build --version
  virtualenv --version
  
In the case of the python package we require the *development* package for python this must be installed for the python bindings to become available.  For the Ubuntu distribution you can get the missing packages with the following commands::

  sudo apt-get install gfortran
  sudo apt-get install cmake
  sudo apt-get install pythonX.Y-dev # Where X and Y are the major and minor version numbers of the Python you want to install, any version above 2.6 will work
  sudo apt-get install swig
  sudo apt-get install python-sphinx
  sudo apt-get install python-virtualenv
  
Install sphinx-fortran extension with pip::

  sudo pip install sphinx-fortran

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


.. _CMakeDownload: https://cmake.org/download

.. _SWIGDownload: http://www.swig.org/download.html

.. _PythonDownload: https://www.python.org/downloads/
