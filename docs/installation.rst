================================================================================
Installation
================================================================================


Dependencies
--------------------------------------------------------------------------------

The following dependencies are required in order to compile UMD-LETKF.
*(These should all be available in package managers for standard Linux
installations)*

* CMake
* Fortran Compiler *( gfortran >= 5.0, Intel >= 16.0, tested)*
* NetCDF4_ w/ API for Fortran
* MPI *(openmpi and intelmpi tested)*
* BLAS
* LAPACK

The following dependencies are **optional**, depending on the features the
user wants enabled for UMD-LETKF at compile time, and can be installed separately,
or built as part of UMD-LETKF with the appropriate :ref:`cmake-options`

* `WGRIB2 API`_ from NCEP *(optional)*

.. _NetCDF4: https://www.unidata.ucar.edu/downloads/netcdf/index.jsp
.. _WGRIB2 API: http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/



Compiling
--------------------------------------------------------------------------------

First, download the source code, which includes external repositories
(libyaml_ and geoKdTree_)
::

   git clone https://github.com/travissluka/UMD-LETKF.git
   cd UMD-LETKF
   git submodule update --init


Then, create a directory in which UMD-LETKF will be built
::

   mkdir build
   cd build

Configure with cmake, pointing it to the location of the source directory
(the parent directory in this example), and build
::

   cmake ../
   make

CMake might complain about certain libraries not being found, such as NetCDF4_. If
this happens, you might need to specify the path to these libraries with the
`-DCMAKE_PREFIX_PATH=` option. (see :ref:`cmake-options`)

.. _libyaml: https://github.com/yaml/libyaml
.. _geoKdTree: https://github.com/travissluka/geoKdTree


Running Tests
--------------------------------------------------------------------------------

To run the test cases, specify the `-DLETKF_ENABLE_TESTS=ON` when running `cmake`
::

   cmake ../ -DLETKF_ENABLE_TESTS=ON
   make
   ctest

This will download test data and the correct reference solutions, and test to
make sure the answers are within margin of error of being identical to the reference
solutions. Note that the test data and reference solutions are updated every now
and then on Dropbox, so old versions of the GitHub are not guaranteed to pass the
tests. It is therefore recommended that you are using the latest version of UMD-LETKF
before running the test cases.

The above commands will simply tell you whether or not the test cases pass or fail.
To see the actual output of the UMD-LETKF, run with the `-VV` flag. You can also
run specific subsets of tests. Run with the `-N` flag to see a list of the available
tests, and run with `-R <testname>` to run a specific set of tests.

.. warning::
   The reference answers for these tests were generated with GCC7 and so may not match
   Intel compilers that have I have not extensively tested (version 17.0+ ). Hopefully
   this will be updated in the near future. So, don't be alarmed if your ctests fail
   if you're using Intel compiler.


.. _cmake-options:



CMake options
--------------------------------------------------------------------------------

There are several command line options that can be passed to CMake when configuring
the build. The following options are a small subset of the options most relevant
to UMD-LETKF

* **-DCMAKE_BUILD_TYPE=Debug**

  |  compile in debug mode, the code runs slower, but more likely to produce a useful
     error message if there is a unexpected run-time error. Without this flag UMD-LETKF
     will by default compile in **Release** mode

* **-DCMAKE_PREFIX_PATH=...**

  |  used to specify one or more directories in which to search for the required
     libraries (NetCDF, MPI, WGRIB2, etc...) If more than one directory is specified,
     they should be separated by a semicolon and surrounded by quotation marks. E.g.
  
  .. code-block:: bash

     cmake ../ -DCMAKE_PREFIX_PATH="$NETCDF_DIR;$WGRIB2_API_DIR"

* **-DLETKF_ENABLE_GRIB=ON**

  |  Builds the optional grib :ref:`configuration_state` I/O module.
     If enabled, either the path to the wgrib2 api needs to be specified
     in `CMAKE_PREFIX_PATH`, or `LETKF_BUILD_GRIB` needs to be enabled.

* **-DLETKF_BUILD_GRIB=ON**

  |  If `LETKF_ENABLE_GRIB` is on, the wgrib2 api will be downloaded from the NCEP
     server and built before the rest of the UMD-LETKF library is built

* **-DLETKF_ENABLE_TESTS=ON**

  |  The test cases and reference solutions are downloaded. After compiling the tests
     can be run by using `ctest`
