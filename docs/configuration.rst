Configuration
===================

All configuration of UMD-LETKF is done through a single `YAML <https://yaml.org/spec/1.1/>`_ configuration file.

The UMD-LETKF, when run, will by default look for a configuration file in the same directory named `letkf.yaml`. Or, a different file can be specfied on the command line::

  ./letkfdriver <somefile.yaml>


An introduction to the YAML format can be found `here <https://yaml.org/spec/1.1/>`_. Be warned: the UMD-LETKF does not yet perform extensive error checking on the configuration file, meaning an improperly defined configuration file might result in cryptic messages and crash (if you come across this, please let `me know <mailto:tsluka@umd.edu>`_). 


Each of the major sections of the configuration file are described below. Each main section is required, though within a section specific parameters may be optional. There is a :ref:`config_ex` available at the end to help clarify the format of these configuration files.



mpi
---------

These parameters directly affect how UMD-LETKF scatters the model state across the processors.

* **mpi.ens_size=** *<integer>*

  | The number of ensemble members.  **(required)**

* **mpi.ppn=** *<integer>*

  | The number of processors per node. **(default: 1)**
  | Optional, but helps improve I/O performance by evenly distributing the PEs which perform simulttaneous I/O across nodes.

  |
  | example:

.. code-block :: yaml

  mpi:
    ens_size: 20
    ppn: 10



solver
---------

Parameters for the core LETKF solver, mainly the parameters for chosing the type of covariance inflation.

* **solver.save_diag=** *<logical>*

  | If true, the `diag.solver.nc` file is saved at the end of the UMD-LETKF run. **(default: True)**
  | This file contains diagnostic information such as the number of observations used per grid cell.


* **solver.inflation.rtps=** *<real>*

  | The percentage of :ref:`RTPS` to apply. **(default: 0.0)**
  | Valid parameters are betwen 0.0 and 1.0. Value of 0.0 indicates RTPS is off, 1.0 indicates analysis spread is relaxed 100% back toward the background spread. Values between 0.5 and 0.8 are often good choices. You could in theory use values greater than 1.0 to result in analysis spread that is larger than background spread, but I have no idea why you would want to do this. Cannot be used if RTPP is enabled.

* **solver.inflation.rtpp=** *<real>*

  | The percentage of :ref:`RTPP` to apply. **(default: 0.0)**
  | Valid parameters are between 0.0 and 1.0. Value of 0.0 indicates RTPS is off, 1.0 indicates the analysis ensemble perturbations are relaxed 100% back toward the background perturbation. Unless you know what you are doing, you are better off using :ref:`RTPS`. Cannot be used if RTPS is enabled

* **solver.inflation.mul=** *<real>*

  | The amount of :ref:`mul_infl` inflation to apply. **(default: 1.0)**
  | Valid parameters are greater than or equal to 1.0. Value of 1.0 indicates multiplicative inflation is off.

  |
  | example:

.. code-block:: yaml

  solver:
    save_diag: true
    inflation:
      rtps: 0.6
      rtpp: 0.0
      mul:  1.0


state
--------

This section defines what the model state looks like, and how it should read and written. The actual contet of the subsections here may include additional paramters, depending on the stateio class being used.


* **state.verbose=** *<logical>*

  | If true, diagnostic information indicating which processor is responsible for reasing or writing which file. **(default: false)**


* **ioclass=** *<string>*

  | The name of the stateio class to use for reading and writing the model state variables and the grid **(required)**
  | Currently, UMD-LETKF has built in classes
  
  * :ref:`stateio_nc` for handling NetCDF files
  * :ref:`stateio_grib` for GRIB2 files.

    
  | Additional stateio class can be implemented by the user.


* **hzgrid**

  | Definitions for the latitude/longitude/mask of the horizontal grid(s) are specified here. *Currently only ONE horizontal grid can be specified.*

  | Each horizontal grid can specify either a 2D latitude (`lat2d`) or a 1D latitude (`lat1d`). If `lat2d` is specified, `lat1d` can be specified as well. In this case, `lat2d` is used for the actual latitude point of a grid cell in the computation o the LETKF solver, whereas `lat1d` is used as the nominal latitude when saving diagnostic and state output files. This is useful when dealing with curvilinear grids. The same applies for longitude via `lon2d` and `lon1d`.  Additionally, the `mask` field is optional, and is used for masking out gridpoints that should not be evaluated in UMD-LETKF (useful in the ocean.)

  | In the following example, a single horizontal grid named `hz1` is specified, the latitude, longitude, and mask of the grid are obtained from the appropriate variables of the `grid/ocean.hgrid.nc` file.

.. code-block:: yaml
		
  state:
  
    hzgrid:
    - name: hz1
      lat2d: ["geolat", "grid/ocean.hgrid.nc"]
      lon2d: ["geolon", "grid/ocean.hgrid.nc"]
      lat1d: ["lath",   "grid/ocean.hgrid.nc"]
      lon1d: ["lonh",   "grid/ocean.hgrid.nc"]
      mask:  ["wet",    "grid/ocean.hgrid.nc"]


		
  
* **vtgrid**

  | Definitions for depth/height information of the vertical grid(s) are specified here. *Currently only ONE vertical grid can be specified*

  | The vertical grid can be specified as one of the following:

  * **vert0d** - *not yet implemented*
  * **vert1d** - Depth/height values are given, and do not vary with the horizontal grid. If a higher order vertical dimension is specified (**vert2d** or **vert3d**) then this variable will instead serve as the nominal vertical coordinate used when saving the state and diagnostic output files.
  * **vert2d** - *not yet implemented*
  * **vert3d** - *not yet implemented*

  | In the following example, a single vertical grid named `vt1` is specified, the values of the grid are obtained from the appropriate variables of the `grid/ocean.vgrid.nc` file. 

.. code-block:: yaml

  state:
  
    vtgrid:
    -  name: vt1
       vert1d: ["Layer", "grid/ocean.vgrid.nc"]

       
* **statedef**

  | TODO describe the filename templates #ENS?# and #TYPE#

  | This section defines what the state variables that are to be read and written by UMD-LETKF, what their grid specification is, and if there is any optional bounds checking on the final state value or the analysis incrememnt that is applied to the background.

  | A list of variables should be given, which each variable containing the following parameters:

  * **hzgrid** - One of the  horizontal grids as speficifed in the **hzgrid** section.
  * **vtgrid** - One of the vertical grids as specified in the **vtgrid** section.
  * **inpput** - Variable name and file path of the input background ensemble member.
  * **output** - Variable name and file path of the output analysis ensemble member.
  * **ana_bounds** - **(optional)** TODO describe
  * **ana_inc_max** - **(optional)** TODO describe
  
  | 
  | Example:

.. code-block:: yaml

  state:
  
    statedef:
    - name: ocn_temp
      hzgrid: hz1
      vtgrid: vt1
      input:  ["Temp", "bkg/#ENS4#.nc"]
      output: ["Temp", "#TYPE#.#ENS4#.nc"]

    - name: ocn_salt
      hzgrid: hz1
      vtgrid: vt1
      input:  ["Salt", "bkg/#ENS4#.nc"]
      output: ["Salt", "#TYPE#.#ENS4#.nc"]
      ana_bounds: [0, 40]
      ana_inc_max: 1.0

  
localization
-------------------

* **class** - current built in classes are

  * **loc_novrt** -
  * **loc_ocean** -
    
* everything else depends on the specific localization class being used

observation
---------------

* **class**
* obsdef
* platdef




.. _config_ex:

example configuratin file
---------------------------

TODO explain
TODO explain include statement

letkf.yaml

.. code-block:: yaml

  ---
  "#import": "letkf.mpi.yaml"

  solver:
    inflation:
      rtps: 0.6
      rtpp: 0.0
      mul:  1.0


  state:
    class: stateio_nc

    hzgrid:
    - name: hz1
      lat2d: ["geolat", "grid/ocean.hgrid.nc"]
      lon2d: ["geolon", "grid/ocean.hgrid.nc"]
      lat1d: ["lath",   "grid/ocean.hgrid.nc"]
      lon1d: ["lonh",   "grid/ocean.hgrid.nc"]
      mask:  ["wet",    "grid/ocean.hgrid.nc"]

    vtgrid:
    - name: vt1
      vert1d: ["Layer", "grid/ocean.vgrid.nc"]

    statedef:
    - name: ocn_t
      hzgrid: hz1
      vtgrid: vt1
      input:  ["Temp", "bkg/bkg.#ENS4#.TS.nc"]
      output: ["Temp", "#TYPE#.#ENS4#.nc"]

    - name: ocn_s
      hzgrid: hz1
      vtgrid: vt1
      input:  ["Salt", "bkg/bkg.#ENS4#.TS.nc"]
      output: ["Salt", "#TYPE#.#ENS4#.nc"]

    - name: ocn_u
      hzgrid: hz1
      vtgrid: vt1
      input:  ["u", "bkg/bkg.#ENS4#.UV.nc"]
      output: ["u", "#TYPE#.#ENS4#.nc"]

    - name: ocn_v
      hzgrid: hz1
      vtgrid: vt1
      input:  ["v", "bkg/bkg.#ENS4#.UV.nc"]
      output: ["v", "#TYPE#.#ENS4#.nc"]


  localization:
    class: loc_ocean
    hzloc_prof: "0.0 720e3 / 90.0 200.0e3"



letkf.mpi.yaml

.. code-block:: yaml

   ---
   mpi:
     ens_size: 20
     ppn: 10
