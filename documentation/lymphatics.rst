==============================
Pulmonary lymphatic transport
==============================

The ``lymphatics`` module implements the pulmonary fluid and lymphatic
transport equations published by Ashworth et al. (2023).  It was restored from
the historical EdeTede implementation and integrated with the current
``develop`` data structures without including the later surfactant model.

Model inputs
============

The workflow combines two independently calculated inputs:

* perfusion supplies mean capillary pressure, capillary sheet surface area,
  and blood transit time; and
* ventilation supplies the minimum and maximum terminal elastic recoil
  pressure through ``export_terminal_lymphatic_inputs`` and
  ``import_terminal``.

After perfusion geometry and results have been established, call
``define_problem_type('lymphatic_transport')``.  This preserves the existing
perfusion unit fields while allocating the additional lymphatic fields.  Then
import any archived inputs and call ``lymphatic_transport``.

Model parameters
================

The published human values remain the defaults, so existing workflows do not
need to set any parameters.  A value can be changed from Python before calling
``lymphatic_transport`` with::

  from aether.parameter_types import update_lymphatics

  update_lymphatics('lung_mass_g', 1.5)

Calling ``update_lymphatics('help', 0.0)`` prints the current values.
Parameter names and human defaults are:

.. list-table::
   :header-rows: 1

   * - Name
     - Default
     - Meaning
   * - ``lung_mass_g``
     - 639
     - Total lung mass in g
   * - ``breathing_rate_bpm``
     - 15
     - Respiratory rate in breaths/min
   * - ``capillary_hydraulic_conductivity``
     - 4.41335e-8
     - Historical capillary conductivity value used by the published model
   * - ``interstitial_capacity_ml_per_100g``
     - 30
     - Maximum interstitial fluid volume per 100 g of lung
   * - ``initial_interstitial_saturation``
     - 0.48
     - Initial interstitial volume as a fraction of capacity
   * - ``interstitial_compartment_a_fraction``
     - 0.005
     - Fraction of interstitial capacity assigned to compartment A
   * - ``interstitial_pressure_min_mmhg``
     - -8
     - Minimum interstitial pressure in mmHg
   * - ``interstitial_pressure_max_mmhg``
     - -1
     - Maximum interstitial pressure in mmHg
   * - ``lymphatic_pressure_min_mmhg``
     - -8
     - Minimum initial-lymphatic pressure in mmHg
   * - ``lymphatic_pressure_max_mmhg``
     - 1
     - Maximum initial-lymphatic pressure in mmHg
   * - ``lymphatic_density``
     - 1
     - Lymphatic-to-capillary exchange-area ratio
   * - ``lymphatic_saturation_threshold``
     - 0.3
     - Saturation at which the conductivity polynomial becomes active
   * - ``lymphatic_baseline_conductivity_ratio``
     - 1.48
     - Baseline lymphatic-to-capillary conductivity ratio
   * - ``lymphatic_conductivity_coefficient_1`` to ``_6``
     - 845.87, -2416.7, 2388.5, -922.24, 125.85, -0.0067
     - Fifth-order saturation/conductivity polynomial coefficients
   * - ``pressure_phase_offset_radians``
     - 1.570796326794895
     - Lymphatic pressure phase offset relative to breathing
   * - ``integration_steps_per_transit``
     - 96
     - Number of integration steps per capillary transit time
   * - ``convergence_tolerance``
     - 5e-6
     - Interstitial-saturation convergence tolerance

``lymphatic_surface_area_ratio`` is accepted as a descriptive alias for
``lymphatic_density``.  The retained ``lymphatic_integrity``, ``test_time``,
and ``reflection_coefficient`` names do not affect the active published
hydrostatic model: integrity has no implemented valve/backflow equation,
``test_time`` is superseded by convergence control, and the osmotic pathway is
inactive with its default reflection coefficient of zero.

Output fields
=============

``get_ne_lymph_flux`` and ``get_ne_lymph_intsat`` return the element-field
indices for average flux and interstitial saturation.  These are dedicated
fields: lymphatic transport does not overwrite the unstrained vascular radii.
Terminal results are available from ``export_terminal_lymphatic``.

Reproducibility and provenance
==============================

The equations in ``src/lib/lymphatics.f90`` retain the published
implementation.  Integration fixes are limited to elapsed-time reporting,
safe unit-field resizing, explicit terminal-to-unit matching, and dedicated
lymphatic export fields.  Runtime can vary with processor load; numerical
outputs, rather than the elapsed-time message, should be compared when checking
reproducibility.

The complete Python workflow is maintained in the matching
``lymphatics_Ashworth2023`` example in the lung-group-examples repository.

The source implementation contains a historical note that output units may
need further verification for a possible factor-of-1000 conversion.  That
scientific calibration has not been altered by this restoration and should be
resolved in a separately validated change.
