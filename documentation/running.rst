
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
one-based element identifiers and flow values before calling ``evaluate_prq``.  Flow
values use the same units as ``inlet_bc`` and the element flow field (normally
``mm^3/s``)::

  from aether.pressure_resistance_flow import (
      clear_perfusion_flows,
      evaluate_prq,
      set_perfusion_flows,
  )

  set_perfusion_flows([12, 19], [125.0, 75.0])
  evaluate_prq(
      "simple_tree", "rigid", 3, 1.0, "flow", 200.0, 0.0, 0.0, "arterial"
  )

Each selected flow is achieved by adjusting vessel resistances; Poiseuille's law
(plus the hydrostatic term) still holds in every element, so no artificial pressure
jump is introduced.  A controlled element multiplies the resistance of its whole
perfusion territory -- every vessel whose flow passes entirely through it -- by a
factor that the solver finds iteratively.  Factors are nested, so any number of flows can be prescribed at
different depths of the tree, and they combine with elastic vessels, gravity and
remodelling (remodelling changes the radii; the factors are applied on top of the
resulting resistances).  The converged factors are printed after the solve.

Flows implied by conservation are worked out first, starting from the inlet flow.
At a junction where the parent flow is known and every child flow is known
(prescribed or implied), all children are controlled together, so the answer does
not depend on which children were named.  For example, with a fixed MPA flow,
prescribing 70% to the RPA or 30% to the LPA gives identical results.  Such a
junction has one more resistance factor than it has constraints; the remaining
choice is made by keeping the geometric mean of the children's factors at one,
i.e. the smallest symmetric change in resistance.

A prescribed element whose junction is not fully determined (for example a single
branch deep in the tree whose siblings are free) is controlled on its own.  Prescribed
flows must be positive and must not be on an inlet or capillary element.  Before
solving, the prescribed flows (with the inlet flow) are propagated by conservation
up and down the tree: redundant but consistent values are accepted (for example
both the LPA and RPA when they add up to the inlet flow, or two elements of the same
unbranched vessel with the same flow), while values that cannot all hold, or that
would force a zero or negative flow elsewhere, stop the run with an error naming
the node or element.  Flows can be prescribed with ``flow`` or ``pressure``
boundary conditions.

Prescribed flows work with both ``simple_tree`` and ``full_plus_ladder`` meshes.
In a ``simple_tree`` the territory of an element is its subtree.  In a
``full_plus_ladder`` network (built with ``append_units`` and ``add_matching_mesh``)
the territory of an artery is its arterial subtree, the capillaries it feeds and
the veins that drain only those capillaries; for a vein it is the mirror image.
Scaling the whole territory matters: the capillaries and veins carry most of the
resistance along a flow path, so scaling the arteries alone cannot raise a
branch's flow much above its natural value.  Because each arterial territory drains
through its own veins, prescribing the flow in an artery also fixes the flow in
the matching vein, so prescribe a given split on one side only.  In a ladder
network the factor updates are under-relaxed, because the capillary resistances
are recomputed from the solution in the same iteration.

The resistance controller considers the flows matched when the maximum relative
error over all controlled elements is at most ``1e-5``.  If it cannot meet that
tolerance within 500 outer iterations, ``evaluate_prq`` terminates with a non-zero
status instead of returning an approximate solution.  Requested flows must therefore
be achievable under the selected pressure or flow boundary conditions.

The setting persists for later PRQ solves in the same process.  Clear it to restore
the original solver behaviour::

  clear_perfusion_flows()

Calling ``set_perfusion_flows([], [])`` also clears the setting.  Element identifiers
must be unique and valid for the current mesh, and a selected element cannot also be
the inlet element fixed by the standard ``flow`` boundary condition.

Exporting time-dependent ventilation flows
-------------------------------------------

``evaluate_vent`` retains the flow and distal pressure for every airway element
at every solver time step. Export those results explicitly after the solve::

  from aether.ventilation import (
      evaluate_vent,
      export_time_dependent_ventilation,
  )

  evaluate_vent("output/ventilation")
  export_time_dependent_ventilation(
      "output/flows",
      "output/pressures",
  )

``.txt`` is appended to each export filename when needed. The example above
therefore writes ``flows.txt`` and ``pressures.txt``. The established whole-lung
``.exflow`` summary is written by ``evaluate_vent`` using its filename argument.

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

``evaluate_wave_transmission`` retains the total flow and upstream pressure for
every arterial and venous element. Write the waveforms after the solve with::

  from aether.wave_transmission import export_time_dependent_perfusion

  export_time_dependent_perfusion("output/flows", "output/pressures")

This matches the complete vascular waveform export used by the ``peter`` model.
Both files cover one cardiac cycle and use the same column layout as the
ventilation exports:
the first column contains element numbers and each remaining column contains the
value at the time shown in the header.  Perfusion flows are reported in ``mm^3/s``
and pressures in ``Pa``.

Exporting terminal seed mappings
--------------------------------

After ``grow_tree`` completes, the association between each seed data point and
the terminal element generated for it can be exported explicitly. This is useful
when values defined at seed points need to be transferred to the completed tree::

  from aether.growtree import export_terminal_seed_mapping, grow_tree

  grow_tree(
      [22],
      60.0,
      40.0,
      0.4,
      1.0,
      1.2,
      180.0,
      grouping="split",
  )
  export_terminal_seed_mapping("output/RUL_seed_terminal_mapping")

LungSim appends ``.txt`` to the filename. The first column of the resulting file
is the data-point number and the second is the global terminal element number.
The mapping retained in memory is replaced by each subsequent ``grow_tree`` call,
so export it before growing another lobe when separate mapping files are required.

Finally
-------

This procedure of making the Aether library available through a python virtual environment can be repeated for different builds of the Aether library.
The virtual environment is lightweight and provides great encapsulation for development of the library.
