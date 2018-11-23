Diagnostic Output
========================


diag.solver
-----------------

The core LETKF solver will output several fields of diagnostic information to `diag.solver.nc` at the end of the program if `solver.save_diag=true` in the json configuration file. These help get a sense of how many observations are being used by each gridpoint.

TODO: see other documentation for info about localization groups

* **col_maxhz** - the maximum horizontal seach radius for the grid column
* **col_obscount** - the number of observations, for each grid column, that were found within the given col_maxhz radius
* **lg_obscount** - the subset of observations that were allowed to be used for each localization group
* **lg_obsloc** - the sum of the localization values for observations used by each localization group. A single observation can have a localization value between `0.0` and `1.0`.


	    
diag.loc_ocean
--------------------
     
If the :ref:`ocean localizer <loc_ocean>` is used, diagnostic additional fields will be saved to `diag.loc_ocean.nc`.

* **vtloc_surf_lvl** - the depth, in number of levels, of the surface localization group
* **vtloc_surf_depth** - The depth, in meters, of the surface localization group
