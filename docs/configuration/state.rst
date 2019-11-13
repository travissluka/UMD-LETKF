.. _configuration_state:

state
================================================================================

This section defines what the model state looks like, both in terms of the state
variables (:ref:`configuration_state_statedef` section) and the horizontal/vertical
grid (:ref:`configuration_state_hzgrid` / :ref:`configuration_state_vtgrid`
sections). The specifics of what is in each subsections may include additional
parameters, depending on the stateio class being used (but everything is currently
identical for the two builtin state I/O classes provided by default.)

.. _configuration_state_filespec:

.. rubric:: File Specification Format
  
Wherever the ``filevar`` type is used in the sections below, the following format
is used to define which file and variable name the data is pulled from.:

.. code-block:: yaml

  {file: file_name, variable: variable_name}


.. rubric:: Parameters

The following parameters are available regardless of which state I/O class is
selected. Additional parameters, described in subsequent sections here, will
be required depending on which state I/O class is selected.

:class:

   **type:** *string*,  ``required``
   
   The name of the I/O class used to handle the state.
   
   Currently, UMD-LETKF has two built-in classes that can be used, additional stateio
   can be implemented by the user. The following options are available:

   *  **stateio_grib** -
      Handles state I/O through grib2 formatted files. Word of caution: this has
      not been well tested yet! But it should work, I think.
      Only available if built with the :ref:`-DLETKF_ENABLE_GRIB=ON<cmake-options>`
      option.

   *  **stateio_nc** -
      Handles state I/O through NetCDF formatted files.
   
:hzgrid:
   **type:** :ref:`configuration_state_hzgrid` ``required``

   The definitions of the horizontal grid(s).

:statedef:
   **type:** :ref:`configuration_state_statedef` ``required``

   The definitions of the state variables.
   
:verbose:

   **type:** *boolean*, **default:** *false*

   Sets if diagnostic information is printed to the console.
   
   If true, diagnostic information indicating which processor is responsible for
   reading or writing which file will be displayed.

:vtgrid:
   **type:** :ref:`configuration_state_vtgrid` ``required``

   The definitions of the vertical grid(s).

.. rubric:: Example
.. code-block:: yaml
		
  state:
    class: stateio_nc
    verbose: true


.. _configuration_state_hzgrid:

hzgrid
________________________________________________________________________________

One or more horizontal grid specification(s) are given by defining how to read the
latitude, longitude, and optional mask. Each grid will have the following parameters.

.. note::
   Currently only **one** horizontal grid can be specified. Although models
   often produce some variables on a staggered grid, you can still use the
   lat/lon of the grid center, and assuming the localization radius is not
   too small, this limitation should make very little difference.

   
.. rubric:: Parameters

In the following parameters for the horizontal grid specification, at least one
(or both) of the ``1d`` and ``2d`` set of parameters needs to  be defined for
latitude and longitude.
   
:name:

   **type:** *string*, ``required``

   A unique name for the horizontal grid.

   The exact name doesn't really matter, but it is referenced in the subsequent
   :ref:`configuration_state_statedef` sections for assigning a horizontal grid
   to each state variable.
   
:lat2d:

   **type:** :ref:`filevar<configuration_state_filespec>`

   Specifies the 2D latitude grid, in degrees.

:lat1d:

   **type:** :ref:`filevar<configuration_state_filespec>`

   Specifies the 1D latitude grid, in degrees.

   The latitude for each row of the grid will be identical. If ``lat2d`` is
   specified as well, ``lat1d`` will only be used as the nominal latitude for the
   output files. It will not be used to determine lat/lon for each grid-point
   in the LETKF algorithm.

:lon2d:

   **type:** :ref:`filevar<configuration_state_filespec>`

   Specifies the 2D longitude grid, in degrees.

:lon1d:
   
   **type:** :ref:`filevar<configuration_state_filespec>`

   Specifies the 1D longitude grid, in degrees.

   The longitude for each column of the grid will be identical. If ``lon2d`` is
   specified as well, ``lon1d`` will only be used as the nominal longitude for the
   output files. It will not be used to determine lat/lon for each grid-point in
   the LETKF algorithm. 
      
:mask:

   **type:** :ref:`filevar<configuration_state_filespec>`,  **(optional)**

   Specifies the optional mask.

   The mask is optional, but can increase the UMD-LETKF speed in domains such as
   the ocean where land points should skipped over. For the input data, grid-points
   with values of 0.0 are masked out and not used.

.. rubric:: Example

In the following example, a single horizontal grid named ``hz1`` is specified, the
latitude, longitude, and mask of the grid are obtained from the appropriate variables
of the ``grid/ocean.hgrid.nc`` file.

.. code-block:: yaml

  state:
    hzgrid:
    - name: hz1
      lat2d: {file: grid/ocean.hgrid.nc, variable: geolat}
      lon2d: {file: grid/ocean.hgrid.nc, variable: geolon}
      lat1d: {file: grid/ocean.hgrid.nc, variable: lath}
      lon1d: {file: grid/ocean.hgrid.nc, variable: lonh}
      mask:  {file: grid/ocean.hgrid.nc, variable: wet}

      
.. _configuration_state_vtgrid:		

vtgrid
________________________________________________________________________________

Definitions for depth/height information of the vertical grid(s) are specified
here. One or more sets of vertical grids can be defined. 

.. rubric:: Parameters

:name:
   **type:** *string*, ``required``

   A unique name for the vertical grid.

   The exact name doesn't really matter, but it is referenced in the subsequent
   :ref:`configuration_state_statedef` sections for assigning a vertical grid to each
   state variable.
   
:vert0d:
   ``not yet implemented``
   
:vert1d:
   **type:** :ref:`filevar<configuration_state_filespec>`

   Vertical coordinates for a column that don't vary in the horizontal direction.
   
:vart2d:
   ``not yet implemented``

:vert3d:
   ``not yet implemented``

.. note::
   ``vert1d`` can also use a constant value specification for now, see the following
   example. This is needed for surface fields, and is a temporary work around
   until the ``vert0d`` parameter is implemented.
 

.. rubric:: Example  

In the following example one vertical grid named ``vt1`` is specified for the 3D
variables, and another ``vt_surf`` is specified with a constant value (surface)
for the surface only variables

.. code-block:: yaml

  state:
    vtgrid:
    - name: vt1
      vert1d: {file: Vertical_coordinate.nc, variable: Layer}
    - name: vt_surf
      vert1d: {constant: 0.0}



.. _configuration_state_statedef:

statedef
________________________________________________________________________________

This section defines one or more state variables. It defines what the state
variables are that should be read and written by UMD-LETKF, which grid specification
they use, and if there are any optional bounds checking on the final state value
or the analysis increment that is applied
to the background.


.. _configuration_state_placeholders:

.. rubric:: Filename String Placeholders

The ``input`` and ``output`` parameters below can use special placeholders in the
filename string that get replaced at run-time.

- ``#ENSX#``
  This placeholder is replaced with the ensemble number (starting at 1), padded
  with zeros to  ensure the number is of length ``X``. This can also be replaced
  with ``mean`` or ``sprd`` for the ensemble mean and spread output files.
  
- ``#TYPE#``
  This placeholder is replaced with either ``ana`` or ``bkg`` if the output
  file is for the analysis or background.

As an example, the specification string ``ocn.#TYPE#.#ENS4#.nc`` will be used
to generate the following files ``ocn.bkg.mean.nc``, ``ocn.bkg.sprd.nc``,
``ocn.ana.mean.nc``, ``ocn.ana.sprd.nc``, ``ocn.ana.0001.nc``, ``ocn.ana.0002.nc``,
``ocn.ana.0003.nc``, ``...``

.. rubric:: Parameters 

:name:
   **type:** *string*, ``required``

   A unique name for the state variable.

   The exact name doesn't really matter, but it may be referenced in other sections
   of the configuration (such as localization)

:hzgrid:
   **type:** *string*, ``required``

   The name of the horizontal grid to use from the :ref:`configuration_state_hzgrid`
   section of the configuration file.

   The x/y dimensions of the input data given below must match the dimensions
   given for the specified horizontal grid.
   
:vtgrid:
   **type:** *string*, ``required``

   The name of the vertical grid to use from the :ref:`configuration_state_vtgrid`
   section of the configuration file.

   The z dimension of the input data given below must match the dimensions
   given for the specified vertical grid.
   
:input:
   **type:** :ref:`filevar<configuration_state_filespec>`, ``required``

   The file and variable name of per-ensemble background data.

   Each ensemble member is assumed to be in a separate file, and so the input
   filename should use the #ENSX# placeholder.
   (See :ref:`configuration_state_placeholders`)
   
:output:
   **type:** :ref:`filevar<configuration_state_filespec>`, ``required``

   The file and variable name of the per-ensemble, and mean/spread data.

   Each ensemble member is assumed to be in a separate file, and so the
   output filename should used the #ENSX# and #TYPE# placeholders
   (see :ref:`configuration_state_placeholders`). In addition to the analysis
   per-ensemble output, this handles the mean and spread output files for
   the analysis and background.

:ana_bounds:
   **type:** *float[2]*, **(optional)**

   The bounds to which the final analysis should be clamped.
   
:ana_inc_max:
   **type:** *float*, **(optional)**

   The maximum absolute value allowed for the analysis increment.

   Any increment with an absolute value greater than this will be clamped
   (respecting the original sign of the increment).

   
.. rubric:: Example
.. code-block:: yaml

  state:
    statedef:
    - name: ocn_s
      hzgrid: hz1
      vtgrid: vt1
      ana_bounds: [0, 50.0]
      ana_inc_max: 2
      input:  {variable: salt, file: "ocn.bkg.#ENS4#.nc"}
      output: {variable: salt, file: "ocn.#TYPE#.#ENS4#.nc"}
    - name: ocn_ssh
      hzgrid: hz1
      vtgrid: vt_surf
      input:  {variable: ssh, file: "ocn.bkg.#ENS4#.nc"}
      output: {variable: ssh, file: "ocn.#TYPE#.#ENS4#.nc"}
