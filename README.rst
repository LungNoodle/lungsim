
============================
Pulmonary Simulation Library
============================

The code in this repository is a library written in Fortran with bindings for C and Python for simulating pulmonary function.

Documentation
=============

Documentation for this library is available on `readthedocs <http://lungsim.readthedocs.io/>`_, and the status of the documentation build is: |docs_build_badge|.

.. |docs_build_badge| image:: https://readthedocs.org/projects/lungsim/badge/?version=latest
   :target: http://lungsim.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

License
-------

LungSim is available under a dual-licensing model. Non-commercial use is
available under the Apache License 2.0 terms described in `LICENCE.txt
<LICENCE.txt>`_. Commercial use requires a commercial license from the
copyright holder.

This dual-licensing model applies to repository revisions from the commit
that introduced it onward. Earlier revisions remain available under the
license terms that accompanied those revisions.

Citation
--------

If you use LungSim in research, please cite it using the metadata in
`CITATION.cff <CITATION.cff>`_.
Microvascular remodeling
------------------------

For a full_plus_ladder perfusion mesh, use the Python binding
``evaluate_prq(mesh_type, vessel_type, grav_dirn, grav_factor,
bc_type, inlet_bc, outlet_bc, remodeling_grade, remodeling_side)``
from ``aether.pressure_resistance_flow``. ``remodeling_side`` is optional
and defaults to ``"none"``. Grade zero is a healthy run: leave the side out
(any side given is ignored). For a grade above zero, set
``remodeling_side`` to ``"arterial"``, ``"venous"``, or ``"both"``; the run
stops with an error if it is missing or ``"none"``.

For wave simulations, ``evaluate_wave_transmission`` has the same optional
final ``remodeling_side`` argument.
Use the same grade and side in both calls.

Selection uses the healthy nominal vessel diameter within each ladder unit.
Hypertrophy affects 100-1000 micrometres; narrowing affects 30-300
micrometres; fibrosis affects 30-500 micrometres; and effective pruning
affects diameters below 300 micrometres (below 500 from grade 50). The symmetric ladder represents
pruning as a fractional loss of parallel vessels, rather than individual
occlusions. The ninth generation contains capillary sheets only, so tubular
arteriole and venule remodeling applies to generations 1-8. The capillary
sheets themselves are unchanged.
