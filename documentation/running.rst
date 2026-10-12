
=======
Running
=======

The pulmonary simulation library, as it's name suggests, is a library and as such the concept of running the library doesn't really exist.  However you can make use of the applications that rely on this library.  A small repository of applications is availble from `here <https://github.com/LungNoodle/lungapps>`_.  

Making use of virtual environments
==================================

Possibly the best way to make use of the Aether library is through the python bindings made available through virtual environments.  This allows a correspondance of library configuration to virtual environment that can be easily moved between i.e. for each build configuration Debug, Release we create a corresponding virtual environment.

Follow these steps for creating a Python virtual environment from which the Aether library will be available.

Create a home for all virtual environments
------------------------------------------

The first task is to create a directory to hold the virtual environment installations::

  mkdir virtual_environments
  
This directory can be created anywhere.

Create a virtual environment
----------------------------

The second task is to create a Python virtual environment to install the Aether python modules into::

  cd virtual_environments # change directory to where the virtual environment should be created
  python -m venv venv-develop-release
  
The name of the virtual environment (in this case *venv-develop-release*) is determined from the branch of the Aether library and the configuration that is going to be available.

Activate virtual environment
----------------------------

The third task is to activate the Python environment.  This can be done by executing a shell script made available in the installation, for POSIX systems execute the command::

  source /path/to/virtual_environments/venv-develop-release/bin/activate
  
for Windows the equivalent command is::

  \path\to\virtual_environments\venv-develop-release\Scripts\activate
  
The activate script may alter the command prompt to indicate the active virtual environment.
This script will also make changes to your path variables.
To undo these changes execute the *deactivate* script::

  deactivate
  
Install Aether into virtual environment
---------------------------------------

With the virtual environment activated from above, change directory into the lungsim build directory::

  cd /path/to/lungsim-build/

From this directory install the aether package with *pip* (POSIX systems)::

  pip install -e src/bindings/python

For Windows the command is slightly different::

  pip install -e src/bindings/python/Release

Here the *Release* directory is the configuration that was build previously, we may also build a *Debug* configuration and we would have to change the location if we were wanting to install that configuration to our virtual environment.
  
This will create a link from the active virtual environment to the Aether library.
Thus making the Aether python library available from the the currently active Python environment.

Test Aether in virtual environment
----------------------------------

With the virtual environment active that Aether is linked to, run python to get a command prompt::
  
  python
  
and at the command prompt enter the following::

  >>> from aether.diagnostics import set_diagnostics_on
  
if all has gone correctly ... nothing should happen! Another command prompt should appear::

  >>>

If the above command was successful then the Python applications given above will also run successfully.

Prescribing selected element flows
----------------------------------

The pressure-resistance-flow solver retains its original behaviour unless prescribed
perfusion flows are explicitly configured.  From Python, pass matching lists of
one-based element identifiers and flow values before calling ``evaluate_prq``::

  from aether.pressure_resistance_flow import (
      clear_perfusion_flows,
      evaluate_prq,
      set_perfusion_flows,
  )

  set_perfusion_flows([12, 19], [125.0, 75.0])
  evaluate_prq(
      "simple_tree", "rigid", 3, 1.0, "flow", 200.0, 0.0, 0.0
  )

Each selected element flow is fixed to the supplied value while the solver computes
the remaining flows and pressures.  Nodal flow conservation is retained, so fixing
the flow in the unique element entering a subtree fixes the total flow delivered to
that subtree.  The prescribed-flow condition replaces that element's pressure-drop
equation, allowing the pressure jump required to impose the flow to be determined by
the rest of the network.  The solver keeps its topology-aware equation ordering so
the reduced matrix has the structural diagonal required by the ILU preconditioner.

The setting persists for later PRQ solves in the same process.  Clear it to restore
the original solver behaviour::

  clear_perfusion_flows()

Calling ``set_perfusion_flows([], [])`` also clears the setting.  Element identifiers
must be unique and valid for the current mesh, and a selected element cannot also be
the inlet element fixed by the standard ``flow`` boundary condition.

Exporting time-dependent ventilation flows
-------------------------------------------

``evaluate_vent`` always exports the flow and distal pressure for every airway
element at every solver time step.  Pass the two required output filenames::

  from aether.ventilation import evaluate_vent

  evaluate_vent(
      "output/flows",
      "output/pressures",
  )

``.txt`` is appended to each filename when needed.  The example above therefore
writes ``flows.txt`` and ``pressures.txt``.  The established whole-lung ``.exflow``
summary is also retained and uses the flow filename as its stem.

Both text exports contain the final simulated breath, since earlier breaths are
used to converge the ventilation solution.  The first column contains global
element numbers.  Each remaining column contains the signed element flow or the
pressure at the element's distal node at the time shown in the header.  The number
and spacing of the time columns follow the solver time step; they are not fixed to
a particular number of samples.  Comment lines at the top of each file state the
time and value units.  Flows use the model's reported flow unit and pressures are
reported in ``cmH2O``.

Exporting time-dependent perfusion waveforms
---------------------------------------------

``evaluate_wave_transmission`` writes the total flow and distal pressure for every
arterial element to ``flows.txt`` and ``pressures.txt``.  Both files cover one
cardiac cycle and use the same column layout as the ventilation exports: the first
column contains element numbers and each remaining column contains the value at the
time shown in the header.  Perfusion flows are reported in ``mm^3/s`` and pressures
in ``Pa``.

Finally
-------

This procedure of making the Aether library available through a python virtual environment can be repeated for different builds of the Aether library.
The virtual environment is lightweight and provides great encapsulation for development of the library.
