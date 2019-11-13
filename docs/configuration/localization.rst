.. _configuration_localization:

localization
================================================================================

The localization class determines how spatial and temporal localization is performed,
a crucial aspect of how an LETKF operates. This include localization for horizontal,
vertical, temporal, and state variable components.

.. warning::
   The localization configuration is perhaps the least finished part of UMD-LETKF.
   Things here will likely change quite a bit as localization methods are added
   for other domains, and/or a way to generically and flexibly specify localization
   is added.
   
.. rubric:: Parameters

The following parameters are available regardless of which localization class is
selected. Additional parameters, described in subsequent sections here, will be
required depending on which localization class is selected.

:class:
   
  **type:** *string*, ``required``

  The name of the localization class to use.

  Currently, UMD-LETKF has two built-in classes that can be used, additional
  localization classes can be implemented by the user. The following options are
  available:

  * | :ref:`configuration_localization_loc_novrt` -
      A generic class that has no localization in the vertical, only
      in the horizontal.

  * | :ref:`configuration_localization_loc_ocean` -
      Localization specific to the ocean.

.. note::
   All localization radii defined below are given as a standard deviation
   of a Gaussian. (Even though they are implemented as a compact Gaspari-Cohn
   function)


   
.. _configuration_localization_loc_novrt:

loc_novrt
--------------------------------------------------------------------------------

The ``loc_novrt`` localization class implements a basic horizontal-only localization.
The bare minimum needed to have a working LETKF. No vertical localization is performed.

.. rubric:: Parameters

:hzloc:
   **type:** :ref:`configuration_localization_hzloc`, ``required``

   Horizontal localization specification used for all observation types.

   
.. rubric:: Example

.. code-block:: yaml	    

  localization:
    class: loc_novrt
    hzloc:
      type: linearinterp_lat
      value:
      - {lat: 0.0,  radius: 500.0e3}
      - {lat: 90.0, radius:  50.0e3}


.. _configuration_localization_loc_ocean:

loc_ocean
--------------------------------------------------------------------------------

The ``loc_ocean`` localization class implements a localization strategy specific to
the ocean. Namely, satellite and insitu observations can be given a different
horizontal localization radius (given the abundance of satellite observations
compared to insitu, satellite observations should be given a smaller horizontal
localization radius). Also, vertical localization of the satellite observations to
just the ocean mixed layer, can be enabled

.. rubric:: Parameters

:save_diag:
   **type:** *logical*, **default:** *true*

   If true, diagnostic information specific to the ocean localization will
   be saved. See :ref:`diagnostics_diag.loc_ocean` for more information on the
   fields that are saved.
   
:diag_file:
   **type:** *string*, **default:** *diag.loc_ocean.nc*

   The file to which ocean localization diagnostics are saved, if ``save_diag``
   is set to true. See :ref:`diagnostics_diag.loc_ocean` for more information
   
:hzloc_prof:
   **type:** :ref:`configuration_localization_hzloc`, ``required``

   The horizontal localization specification for insitu profiles.

   Insitu profiles are determined to be the observations and platform types
   that are NOT included in the following ``sat_obs`` or ``sat_plats``
   parameters.

:hzloc_sat:
   **type:** :ref:`configuration_localization_hzloc`, ``required``

   The horizontal localization specification for satellite observations.

   Satellite observations are determined to be the observations and platform
   types that are included in the following ``sat_obs`` or ``sat_plats``
   parameters. For each observation, if its type matches one listed in ``sat_obs``,
   or its platform type matches one listed in ``sat_plats``, it is considered a
   satellite observation (it does not have to match both).
   
:tloc_prof:
   **type:** *float*, **default:** *-1.0*

   Temporal localization for insitu profiles (in hours).
   If < 0, temporal localization is disabled.
   
:tloc_sat:
   **type:** *float*, **default:** *-1.0*

   Temporal localization for satellite observations (in hours).
   If < 0, temporal localization is disabled.
   
:vtloc_surf:
  **type:** :ref:`configuration_localization_vtloc`, **default:** *type=none*

  The vertical localization specification for satellite observations.

  Insitu profiles do not have any vertical localization.
  
:sat_obs:
  **type:** *array of strings*, **optional**

  An array of observation names that are to be treated as satellite observation
  for localization purposes. See :ref:`configuration_observation_names`.
  
:sat_plats:
  **type:** *array of strings*, **optional**

  An array of platform names that are to be treated as satellite observations
  for localization purposed. See :ref:`configuration_observation_names`
  

.. _configuration_localization_vtloc:

.. rubric:: ``vtloc`` Parameters

Specification of the vertical localization.

:type:
  **type:** *string*, **default:** *none*

  The type of vertical localization to use for the ocean. Currently two options
  are available:
  
  - **none** - vertical localization is off, observations impact the
    entire vertical column.

  - **bkg_t** - surface observations are localized to the surface mixed layer,
    as calculated from a change in background temperature criteria.

:bkg_t_delta:
  **type:** *float*, ``required``

  The change in background temperature (Celsius) from the surface to some depth,
  used for calculating the depth of the ocean mixed layer.
  
:bkg_t_var:
  **type:** *string*

  The name of the background temperature variable used for calculating the mixed
  layer depth. This state variable name must be one of those given in 
  :ref:`state.statedef<configuration_state_statedef>`.

       
.. rubric:: Example

.. code-block:: yaml
		
  localization:
    class: loc_ocean
    save_diag: true
    hzloc_prof:
      type: linearinterp_lat
      value:
      - {lat: 0.0,  radius: 720.0e3}
      - {lat: 90.0, radius: 200.0e3}
    hzloc_sat:
      type: linearinterp_lat
      value:
      - {lat: 0.0,  radius: 500.0e3}
      - {lat: 90.0, radius: 50.0e3}
    sat_plats:
     - ocn_sat
    vtloc_surf:
      type: bkg_t
      bkg_t_delta: 0.2
      bkg_t_var: ocn_t



Common Types
--------------------------------------------------------------------------------

.. _configuration_localization_hzloc:

.. rubric:: ``hzloc`` Parameters

This parameter type is used to specify the characteristics of the horizontal
localization.

:type:
  **type:** *string*, ``required``

  The type of horizontal localization to use. Currently, the only valid option
  is ``linearinterp_lat``. This type gives a horizontal localization radius
  that changes with latitude. Several latitudes are specified, along with the
  desired radius, and linear interpolation is used to calculate the radius
  for any valid latitude.

:value:
  **type:** *array of lat/radius values*

  An array of ``lat`` / ``radius`` pairs. See the example below for clarification.
  
  - **lat:** absolute value of latitude in degrees
  - **radius:** horizontal localization radius, meters. Given as the standard
    deviation of a Gaussian.

Note that all latitude values are positive. Currently, different values cannot
be given for southern/northern hemisphere. If 0.0 and 90.0 are not included
in the list of latitudes, they are implicitly added using the radius of the
nearest given latitude.
    
.. rubric:: Example

Note that in this example a latitudes between 0.0 degrees and 5.0 degrees have a
localization radius of 500 km, all latitudes above 50.0 degrees have a radius of
100 km. In between they are appropriately linearly interpolated.

.. code-block:: yaml

  hzloc:
    type: linearinterp_lat
    value:
      - {lat:  5.0, radius: 500.0e3}
      - {lat: 10.0, radius: 300.0e3}
      - {lat: 50.0, radius: 100.0e3}
