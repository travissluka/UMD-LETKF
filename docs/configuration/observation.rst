.. _configuration_observation:

observation
================================================================================

Parameters for specifying how the observations are read in, or how synthetic
test observations are generated, are controlled under this section. The exact
contents of this configuration section depends on which I/O ``class`` is used.


.. rubric:: Parameters

The following parameters are available regardless of which observation reader
class is selected. Additional parameters, described in subsequent sections here,
will be required depending on which localization class is selected.

:class:
   **type:** *string*, ``required``

   The name of the observation I/O class to use.

   Currently, the UMD-LETKF has three built-in classes. Additional classes
   may be implemented by the user. The following options are available, and
   their specific configuration requirements are described in the following
   sections.

   *  :ref:`configuration_observation_obsio_ioda` -
      Provides the ability to read observation files in the JEDI IODA format
      from the JEDI hofx application. 

   * :ref:`configuration_observation_obsio_nc` -
     A generic NetCDF4 file reader.

   * :ref:`configuration_observation_obsio_test` -
     Generates synthetic observations from a specified increment value.
     
.. note::
   If you're hooking up your model to the UMD-LETKF for the first time, you're
   best bet is to use the :ref:`configuration_observation_obsio_test` reader
   first (to make sure everything else is hooked up correctly), before trying the
   :ref:`configuration_observation_obsio_nc` or
   :ref:`configuration_observation_obsio_ioda` readers with real observations.
   
.. _configuration_observation_placeholder:

.. rubric:: Filename Ensemble Placeholder

Some of the classes below require filenames for the per-ensemble member observation
input. In these cases the ``#ENSX#`` placeholder can be used within the string of
the filename. It is replaced with the ensemble member number (starting at 1),
padded with zeros to ensure the number is ``X`` digits long. For example
``sst_obs.#ENS4#.nc`` will be substituted as
``sst_obs.0001.nc``, ``sst_obs.0002.nc``, ...


.. _configuration_observation_names:

.. rubric:: Observation and Platform Names

The observation I/O classes require that names are given for different observations
and platforms. These can be set to whatever the user wants, and their use can
be considered optional. The exact name is not important, but may be referenced
by other sections of the configuration (such as :ref:`configuration_localization`).
As a general reference, the observation type should reflect which variable is observed
(e.g. ocn_sst, ocn_t) and the platform type can reflect either specific platforms
(e.g viirs, avhrr) or a general satellite vs. insitu.  Hopefully this will make
more sense when seen how it is used in the :ref:`configuration_localization` section.

.. note::
   The observation and platform names should be short, with a 10 character max.

The following documentation describes the observation reader classes that are
available for use.

.. _configuration_observation_obsio_ioda:

obsio_ioda
--------------------------------------------------------------------------------

This observation I/O class can be used to read observation operator output files
from the Joint Effort for Data Assimilation Integration (JEDI) based applications.
Files are in the IODA NetCDF format.
(More of an explanation about this will likely be added once the JEDI repositories
are made public.) Currently only the `hofx` or `hofx3d` applications are supported,
not the `enshofx`. (Odds are you will be wanting to use the `hofx3d` application
only anyway).

.. warning::
  Efficient distribution of the read operations across PEs has not been implemented
  for this class. Large operational size datasets might be a little slow until this
  is fixed.
   
.. rubric:: Parameters

:ioda_files:

   **type:** :ref:`ioda_file<configuration_observation_ioda_file>` list, ``required``
   
   This section contains a list of files that should be loaded, each with the
   following parameters:

   
.. _configuration_observation_ioda_file:

.. rubric:: ``ioda_file`` Parameters

Each set of ioda files to be read requires the following parameters:

:file:

  **type:** *string*, ``required``

  The base name of the file to read.

  The notation of the :ref:`configuration_observation_placeholder`
  should be used since there should be separate files for each individual
  ensemble member. Also, JEDI applications currently produce output files
  for each PE of the application, so the filename given will automatically
  try appending the appropriate ``_0001.nc``, ``_0002.nc``, ``...`` suffixes.
      
:vars:

   **type:** *list of array(s)*, ``required``

   For each desired variable in the input file, an array is given with three
   values that have the following meaning:

   1. observation name - see :ref:`configuration_observation_names`
   2. platform name - see :ref:`configuration_observation_names`
   3. variable name as given in the IODA observation file

   
.. rubric:: Example

.. code-block:: yaml

  observation:
    class: obsio_ioda
    ioda_files:
    - file: mem#ENS1#/sst.out
      vars:
      - [ocn_sst, sst_viirs, sea_surface_temperature]
    - file: mem#ENS1#/insitu.out
      vars:
      - [ocn_s, insitu, sea_water_salinity]
      - [ocn_t, insitu, sea_water_temperature]

	
.. _configuration_observation_obsio_nc:

obsio_nc
--------------------------------------------------------------------------------
The NetCDF reader will read in two types of files. The first is the main observation
file given by the ``filename_obs`` parameter below and the format of which is described
by :ref:`configuration_observation_obsfile`. This file provides each observation
type, location, and value. The second set of files are the per-ensemble member observation
operator files, given by the ``filename_obshx`` parameter below and the format of which is
described by :ref:`configuration_observation_obshxfile`.

.. note::
   Although the configuration here allows for observation data that is common across
   all ensemble members to be specified in a separate ``filename_obs`` file, they
   do not have to be. All observation data could be in the per-ensemble member
   ``filename_obshx`` files. In this case, observation files should contain all
   the data required by both the :ref:`configuration_observation_obshxfile` and
   :ref:`configuration_observation_obsfile` specs, and the ``filename_obs`` should
   simply point to one of the ensemble files.

   
.. rubric:: Parameters

:filename_obs:

   **type:** *string*, ``required``

   The name of the observation file to read in.

   The expected contents of this NetCDF file are specified by
   :ref:`configuration_observation_obsfile`. 

:filename_obshx:

   **type:** *string*, ``required``

   The name of the per-ensemble observation operator file.
   
   The :ref:`configuration_observation_placeholder` should be used to read in
   each individual ensemble member file. The expected contents of this file are
   specified by :ref:`configuration_observation_obshxfile`.
   
   
:obsdef:
   **type:** list of :ref:`obsplat_def<configuration_observation_obsplatdef>`,
   ``required``

   Provides a mapping from the integer values of the observation type in the
   NetCDF file with a human readable name. See also :ref:`configuration_observation_names`.
   
:platdef:
   **type:** list of :ref:`obsplat_def<configuration_observation_obsplatdef>`,
   ``required``

   Provides a mapping from the integer values of the platform type in the
   NetCDF file with a human readable name.
   See also :ref:`configuration_observation_names`.

:read_inc:

   **type:** *boolean*, ``required``

   If true, the values given in the per-ensemble member files are given as
   increments, :math:`y^o - h(x)`, otherwise they are taken as the direct
   output of an observation operator, :math:`h(x)`.

      
.. _configuration_observation_obsplatdef:

.. rubric:: ``obsplat_def`` Parameters

These parameters are required for the ``obsdef`` and ``platdef`` sections of
:ref:`configuration_observation_obsio_nc` and are used to associate a human
readable ``name`` with the integer ``id`` that is stored in the NetCDF file

:name:
   **type:** *string*

   Name of the observation or platform.
   Note the advice of :ref:`configuration_observation_names`
   
:id:

   **type:** *integer*

   The integer value in the NetCDF file.

:description:

   **type:** *string*

   Optional description of the observation or platform type.
   Not needed by UMD-LETKF other than for the sanity of the user.
   
   
.. rubric:: Example

.. code-block:: yaml

  observation:
    class: obsio_nc
    obsdef:
    - name: ocn_t
      id:   2210
      description: "ocean insitu temperature (C)"
    - name: ocn_s
      id:   2220
      description: "ocean salinity (PSU)"
    platdef:
    - name: ocn_prf
      id:   1
      description: "all insitu obs"
    - name: ocn_sat
      id:   1000
      description: "all satellite based obs"
    filename_obs:   obs.nc
    filename_obshx: "obs.#ENS4#.nc"
    read_inc: false


.. _configuration_observation_obsfile:

Observation File Format
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The NetCDF file containing observation data needs to contain the following
dimensions and variables of the same name. An example file can be found in the
test data for UMD-LETKF.

.. note::
   I realize the variable "depth" is required and that that "height" is not
   a valid option. Since UMD-LETKF was started for ocean DA, this will
   be addressed once non-ocean localization classes are implemented.


.. rubric:: dimensions
	    
:obs:	    
  Number of observations in the file
   
.. rubric:: variables
	    
All variables here are of size ``obs``

:depth:
   **type:** *float*

   The depth of the observation in meters.

   
:err:
   **type:** *float*

   The standard deviation of the observation error

:hr:
   **type:** *float*

   The time offset (in hours) from the analysis time. Only actually used if temporal
   localization is used.
   
:lat:
   **type:** *float*

   Latitude in degrees
   
:lon:
   **type:** *float*

   Longitude in degrees
   
:obid:
   **type:** *integer*

   The observation id. See :ref:`configuration_observation_obsplatdef`
   
:plat:
   **type:** *integer*

   The platform id. See :ref:`configuration_observation_obsplatdef`

:qc:
   **type:** *integer*

   Quality control flag. Observation is used by UMD-LETKF only if ``qc`` is zero.

:val: 
   **type:** *float*

   The value of the observation.

   
.. _configuration_observation_obshxfile:
	    
Observation H(x) File Format
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The NetCDF file containing per-ensemble member observation operator data needs
to contain the following dimensions and variables of the same name. An example
file can be found in the test data for UMD-LETKF.

.. rubric:: dimensions

:obs:

   Number of observations in the file.
   
.. rubric:: variables	    

All variables here are of size ``obs``

:hx:
   **type:** *float*

   The value of the observation operator from a single ensemble member background.
   This can either contain the value (:math:`h(x)`), or the observation increment
   (:math:`y^o-h(x)`), depending on the value of ``read_inc`` in
   :ref:`configuration_observation_obsio_nc`
   
.. _configuration_observation_obsio_test:
   
obsio_test
--------------------------------------------------------------------------------

This observation I/O class can be used to generate synthetic observations from
the state background mean using a specified increment. This method can be useful
when wanting to perform a quick single-obs test, bypassing the need to generate
observation files. Test observations can only be generated directly from the state
background (i.e. the identity observation operator is used.)

.. rubric:: Parameters

:synthetic_obs:

   This section contains an array of arrays (see the example below if that doesn't
   make sense). Each observation specification contains an array of nine values,
   in the following order
   
   1. **observation_id** - A string reflecting the type of observation
      (see :ref:`configuration_observation_names`).
   2. **platform_id** - A string reflecting the type of platform
      (see :ref:`configuration_observation_names`).
   3. **state_variable** - The state variable that this observation is generated
      from. The value given must be one of the name of one of the state variables
      given in the :ref:`state.statedef<configuration_state_statedef>` section.
   4. **latitude** - in degrees
   5. **longitude** - in degrees
   6. **depth/height** - The value in the vertical coordinate. If this observation
      is being generated from a 2D surface state field then the depth/height
      here is ignored.
   7. **time** - the time offset (in hours) from the analysis time. This value
      is only used if temporal localization is enabled.
   8. **increment** - The value of this observation will be generated as
      the increment plus the background 
   9. **error** - standard deviation of the observation error


.. rubric:: Example

This example generates two observations from the background temperature, both
with an observation increment of 1 degree and observation error of 0.2 degree.

.. code-block:: yaml

  observations:
    class: obsio_test
    synthetic_obs:
    - [ocn_sst, satellite, ocn_t, 20.0, -140.0,  0.0, 0.0, 1.0, 0.2]
    - [ocn_t,   insitu,    ocn_t, 25.0, -162.0, 10.0, 0.0, 1.0, 0.2]

