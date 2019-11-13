.. _diagnostics:

Diagnostic Output
========================


.. _diagnostics_diag.solver:

LETKF Solver Diagnostics
-------------------------

The core LETKF solver will output several fields of diagnostic information at the end of the program if enabled in the yaml configuration file. See
:ref:`solver.save_diag<configuration_solver>` and
:ref:`solver.diag_file<configuration_solver>` configuration sections.

These help get a sense of how many observations are being used by each grid-point.

* **col_maxhz** - the maximum horizontal search radius for the grid column
* **col_obscount** - the number of observations, for each grid column, that were found within the given col_maxhz radius
* **lg_obscount** - the subset of observations that were allowed to be used for each localization group
* **lg_obsloc** - the sum of the localization values for observations used by each localization group. A single observation can have a localization value between `0.0` and `1.0`. This gives a general sense of the amount of impact observations have.


.. _diagnostics_diag.loc_ocean:

loc_ocean Diagnostics
----------------------
     
If the :ref:`configuration_localization_loc_ocean` is used,  additional fields will be saved to diagnose how vertical localization is performed.

* **vtloc_surf_lvl** - the depth, in number of levels, of the surface localization group
* **vtloc_surf_depth** - The depth, in meters, of the surface localization group
