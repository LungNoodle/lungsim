.. Auckland Bioengineering Lung Simulator master file, created by
   Alys Clark on Fri 18 March 2016.

========================================
Welcome to the documentation for Aether!
========================================

Hello everyone! Welcome to the documentation for our lung code!

The first thing we need to do to get going using this code is to :doc:`acquire </acquiring>` it then we can :doc:`build </building>` the aether library.  Once we have built the library we can :doc:`run </running>` it.  Before you start to develop the code make sure you read the section on :doc:`testing </testing>`.  Lastly, before you ask to get your code merged make sure you have looked at the :doc:`modules </modules>` and :doc:`fortran styles </fortranstyles>` sections to ensure that you have followed the conventions for this codebase.

Thanks and enjoy!

Quick Start
===========

If you have all the requirements for building the aether libraries and bindings then the following commands are a quick description on how to get setup and using the library, for further information or clarification read the full documentation::

  git clone https://github.com/LungNoodle/lungsim.git
  mkdir lungsim-build
  cd lungsim-build
  cmake ../lungsim
  make
  cd ..
  virtualenv --system-site-packages venv_aether
  source venv_aether/bin/activate
  pip install -e lungsim-build/src/bindings/python

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
