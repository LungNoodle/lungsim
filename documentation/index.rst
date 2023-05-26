.. Auckland Bioengineering Lung Simulator master file, created by
   Alys Clark on Fri 18 March 2016.

========================================
Welcome to the documentation for Aether!
========================================

Hello everyone! Welcome to the documentation for our lung code!

The first thing we need to do to get going using this code is to :doc:`acquire </acquiring>` it then we can :doc:`build </building>` the aether library.
Once we have built the library we can :doc:`run </running>` it.
Before you start to develop the code make sure you read the section on :doc:`testing </testing>`.
Lastly, before you ask to get your code merged make sure you have looked at the :doc:`modules </modules>` and :doc:`fortran styles </fortranstyles>` sections to ensure that you have followed the conventions for this codebase.

Thanks and enjoy!

Quick Start
===========

If you have all the requirements for building the aether libraries and bindings then the following commands are a quick description on how to get setup and using the library, for further information or clarification read the full documentation.
Quick start instructions for GNU/Linux and macOS (bash)::

  git clone https://github.com/LungNoodle/lungsim.git
  python -m venv venv-aether
  source venv-aether/bin/activate
  pip install --upgrade pip
  pip install numpy
  here=$(pwd)
  cmake -S lungsim -B build-lungsim -D Python_EXECUTABLE=$here/venv-aether/bin/python
  cd lungsim-build
  make
  pip install -e src/bindings/python
  
Quick start instructions for Windows (cmd)::

  git clone https://github.com/LungNoodle/lungsim.git
  python -m venv venv-aether
  venv-aether/Scripts/activate
  pip install --upgrade pip
  pip install numpy
  cmake -S lungsim -B build-lungsim -D Python_EXECUTABLE=%cd%/venv-aether/Scripts/python
  cd lungsim-build
  cmake --build . --config Release
  pip install -e src/bindings/python/Release

Please note: On Windows using the Intel oneAPI Fortran compilers we have found that we are only able to build the Aether library with Visual Studio 16 2019.
At this time, we have not been able to successfully compile the Aether library with Visual Studio 17 2022 and Intel oneAPI.
The most recent attempt was performed with CMake version 3.26.4.

Contents
========

.. toctree::
   :maxdepth: 2

   acquiring
   building
   testing
   running
   modules
   fortranstyles
