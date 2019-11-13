.. _configuration_solver:

solver
================================================================================

Parameters for the core LETKF solver, diagnostic output, and covariance inflation.


.. rubric:: Parameters

:save_diag:

  **type:** *boolean*, **default:** *true*
  
  Whether a diagnostic file is saved at the end of the UMD-LETKF run.
  
  This file contains diagnostic information such as the number of observations
  used per grid cell, maximum horizontal localization radius, etc. See
  :ref:`diagnostics_diag.solver`

:diag_file:

  **type:** *string*,  **default:** *"diag.solver.nc"*
  
  The file name of the diagnostic output.

  Only used if ``save_diag`` is set to true.


.. rubric:: Example
.. code-block:: yaml

  solver:
    save_diag: true
    diag_file: "diag.solver.nc"


.. _configuration_solver_inflation:

inflation
--------------------------------------------------------------------------------

Covariance inflation methods for increasing the ensemble spread. This section is
optional. If not provided the default of "no inflation" will be used.

.. note::
   The `RTPS` and `RTPP` methods cannot be enabled simultaneously.


.. rubric:: Parameters

:mul:

   **type:** *float*, **default:** *1.0*

   The amount of :ref:`mul_infl` inflation to apply.
  
   Valid parameters are greater than or equal to 1.0. Value of 1.0 indicates
   multiplicative inflation is off.


:rtpp:

   **type:** *float*, **default:** *0.0*

   The percentage of :ref:`RTPP` to apply.

   Valid parameters are between 0.0 and 1.0. Value of 0.0 indicates RTPS is off,
   1.0 indicates the analysis ensemble perturbations are relaxed 100% back toward
   the background perturbation. Unless you know what you are doing, you are better
   off using :ref:`RTPS`. Cannot be used if RTPS is enabled

:rtps:

   **type:** *float*, **default:** *0.0*

   The percentage of :ref:`RTPS` to apply.

   Valid parameters are between 0.0 and 1.0. A value of 0.0 indicates RTPS is off,
   1.0 indicates analysis spread is relaxed 100% back toward the background spread.
   Values between 0.5 and 0.8 are often good choices. You could in theory use values
   greater than 1.0 to result in analysis spread that is larger than background
   spread, but I have no idea why you would want to do this. Cannot be used if RTPP
   is enabled.


.. rubric:: Example
.. code-block:: yaml

  solver:
    inflation:
      rtps: 0.6
      rtpp: 0.0
      mul:  1.0

