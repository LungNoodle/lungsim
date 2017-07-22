
=======
Running
=======

The pulmonary simulation library, as it's name suggests, is a library and as such the concept of running the library doesn't really exist.  However you can make use of the applications that rely on this library.  A small repository of applications is availble from `here <https://github.com/LungNoodle/lungapps>`_.  

Making Use of Virtual Environments
==================================

Possibly the best way to make use of the Aether library is through the python bindings made available through virtual environments.  This allows a correspondance of library configuration to virtual environment that can be easily moved between.

Follow these steps for creating a python virtual environment from which the aether library will be available.

Create a home for all virtual environments
------------------------------------------

The first task is to create a directory to hold the virtual environment installations::

  mkdir virtual_environments
  
This directory can be created anywhere.

Create a virtual environment
----------------------------

The second task is to create a python virtual environment to install the aether python modules into::

  cd virtual_environments # change directory to where the virtual environment should be created
  virtualenv --system-site-packages develop
  
The *--system-site-packages* flag allows the virtual environment to access all the packages that the system python has installed.  This is useful for big packages which may be required, for example; numpy or scipy.  The name of the virtual environment (in this case *develop*) is determined from the branch of the aether library that is going to be available.

Activate virtual environment
----------------------------

The third task is to activate the python environment.  This can be done by executing a shell script made available in the installation, for POSIX systems execute the command::

  source /path/to/env/bin/activate
  
for Windows the equivalent command is::

  \path\to\env\Scripts\activate
  
The activate script may alter the command prompt to indicate the active virtual environment.  This script will also make changes to your path variables.  To undo these changes execute the *deactivate* script::

  deactivate
  
Install Aether into virtual environment
---------------------------------------

With an active python virtual environment change directory into the lungsim build directory::

  cd /path/to/lungsim-build/
  
From this directory change into the *bindings/python* directory::

  cd bindings/python
  
in this directory a python file named *setup.py* should exist.  To make the aether library available via the active virtual environment execute the following command::

  python setup.py develop
  
This will create a link from the active virtual environment to the aether library.  Thus making the aether python library available from the the currently active python environment.

Test Aether in virtual environment
----------------------------------

With the virtual environment active that aether is linked to, run python to get a command prompt::
  
  python
  
and at the command prompt enter the following::

  >>> from aether.diagnostics import set_diagnostics_on
  
if all has gone correctly ... nothing should happen! Another command prompt should appear::

  >>>

If the above command was successful then the python applications given above will also run successfully.

Finally
-------

This procedure of making the aether library available through a python virtual environment can be repeated for different builds of the aether library.  The virtual environment is lightweight and provides great encapsulation for development of the library.
