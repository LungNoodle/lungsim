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
Install the requirements for :ref:`building:GNU/Linux` or :ref:`building:Windows`

If you have all the requirements for building the aether libraries and bindings then the following commands are a quick description on how to get setup and using the library, for further information or clarification read the full documentation.

Quick start instructions for GNU/Linux and macOS

(bash)::

  git clone https://github.com/LungNoodle/lungsim.git
  python -m venv venv-aether
  source venv-aether/bin/activate
  pip install --upgrade pip
  pip install numpy
  here=$(pwd)
  cmake -S lungsim -B build-lungsim -D Python_EXECUTABLE=$here/venv-aether/bin/python
  cd build-lungsim
  make
  pip install -e src/bindings/python
  python ../lungsim/.github/scripts/diagnostics_test.py
  
Quick start instructions for Windows

(Note: replace paths for -D SWIG_EXECUTABLE and -D SWIG_DIR with the relevant paths for your SWIG install)

(cmd)::

  git clone https://github.com/LungNoodle/lungsim.git
  python -m venv venv-aether
  python -m pip install --upgrade pip
  venv-aether\Scripts\activate
  pip install numpy
  cmake -S lungsim -B build-lungsim -D Python_EXECUTABLE=%cd%/venv-aether/Scripts/python -D SWIG_EXECUTABLE="C:/Program Files (x86)/SWIG/swigwin-4.1.1/swig.exe" -D SWIG_DIR="C:/Program Files (x86)/SWIG/swigwin-4.1.1/lib"
  cd build-lungsim
  cmake --build . --config Release
  pip install -e src\bindings\python\Release
  python ..\lungsim\.github\scripts\diagnostics_test.py



If all is successful, the final command should run a test python script with the output "test succeeded".