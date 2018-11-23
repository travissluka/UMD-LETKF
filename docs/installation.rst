==============
Installation
==============

Dependencies
----------------
* Fortran ( gfortran >= 5.0, intel 16.0, tested)
* NetCDF4_ w/ API for Fortran
* MPI (openmpi, intelmpi, ...)
* `WGRIB2 API`_ from NCEP (optional)

.. _NetCDF4: https://www.unidata.ucar.edu/downloads/netcdf/index.jsp
.. _WGRIB2 API: http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/ 

Compiling
---------------

First, download the source code, including the required external repositories (fson_ and geoKdTree_)
::
   
   git clone https://github.com/travissluka/UMD-LETKF.git
   cd UMD-LETKF
   git submodule update --init


Then, create a directory in which UMD-LETKF will be built
::
   
   mkdir build
   cd build

Configure with cmake, pointing it ot the location of the source directory (the parent directory in this example), and build
::
   
   cmake ../
   make

cmake might complain about certain libraries not being found, such as NetCDF4_. If this happens, you might need to specify the path to these libraries with the `-DCMAKE_PREFIX_PATH=` option.

.. _fson: https://github.com/josephalevin/fson
.. _geoKdTree: https://github.com/travissluka/geoKdTree


Running Tests
---------------

To run the test cases, specify the `-DLETKF_ENABLE_TESTS=ON` when running `cmake`
::
   
   cmake ../ -DLETKF_ENABLE_TESTS=ON
   make
   ctest

This will download test data and the correct reference solutions, and test to make sure the answers are within margin of error of being identical to the reference solutions. Note that the test data and reference solutions are updated every now and then on Dropbox, so old versions of the github are not guaranteed to pass the tests. It is therefore recommended that you are using the latest version of UMD-LETKF before running the test cases.

The above commands will simply tell you whether or not the test cases pass or fail. To see the actual output of the UMD-LETKF, run with the `-VV` flag. You can also run specific subsests of tests. Run with the `-N` flag to see a list of the available tests, and run with `-R <testname>` to run a specific set of tests.


CMake options
------------------

There are several command line options that can be passed to CMake when configuring the build. The following options are a small subset of the options most relevant to UMD-LETKF

* **-DCMAKE_BUILD_TYPE=Debug**

  |  compile in debug mode, the code runs slower, but more likely to produce a useful error message if there is a unexpected runtime error

* **-DCMAKE_PREFIX_PATH=...**

  |  used to specify one or more directories in which to search for the required libraries (NetCDF, MPI, WGRIB2, etc...) If more than one directory is specified, they should be separatated by a comma and surrounded by quotation marks.
  |  E.g. `cmake ../ -DCMAKE_PREFIX_PATH="$NETCDF_DIR;$WGRIB2_API_DIR"`

* **-DLETKF_ENABLE_GRIB=ON**

  |  Builds the optional grib state I/O module. If enabled either the path to the wgrib2 api needs to be specified in `CMAKE_PREFIX_PATH`, or `LETKF_BUILD_GRIB` needs to be enabled.

* **-DLETKF_BUILD_GRIB=ON**

  |  If `LETKF_ENABLE_GRIB` is on, the wgrib2 api will be downloaded from the NCEP server and built before the rest of the UMD-LETKF library is built

* **-DLETKF_ENABLE_TESTS=ON**
  
  |  The test cases and reference solutions are downloaded. The after compiling, the test can be run by using `ctest`
 
