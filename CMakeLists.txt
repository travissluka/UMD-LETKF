# Copyright 2016-2019 Travis Sluka
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#-------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(letkf Fortran C CXX)

# set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
  "Build type, options are: Debug Release RelWithDebInfo" FORCE)
endif()


#------------------------------------------------------------
# determine the repository version
#------------------------------------------------------------
execute_process(COMMAND git describe --dirty --tags --always
  OUTPUT_VARIABLE GIT_VERSION
  RESULT_VARIABLE RES
  OUTPUT_STRIP_TRAILING_WHITESPACE
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

if(NOT ${RES} EQUAL "0")
  set(GIT_VERSION "Unknown")
endif()
add_definitions(-DCVERSION=\"${GIT_VERSION}\")
message(STATUS "Git repository version: ${GIT_VERSION}")



#------------------------------------------------------------
# Testing and code coverage
#------------------------------------------------------------
option( LETKF_ENABLE_CODECOVERAGE "enable code coverage testing support")
if(LETKF_ENABLE_CODECOVERAGE)
  # TODO, this is only valid for GCC
  add_definitions( -fprofile-arcs -ftest-coverage )
  set( CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} --coverage )
endif()

option( LETKF_ENABLE_TESTS "enable testing suite, including downloading of test data")
if(LETKF_ENABLE_TESTS)
  enable_testing()
endif()



#------------------------------------------------------------
#  dependencies
#------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# MPI
find_package(MPI REQUIRED)

# NetCDF
set(NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)

# LAPACK / BLAS
find_package(LAPACK REQUIRED)

# grib2
option( LETKF_ENABLE_GRIB "enable grib2 support for state I/O" OFF)
option( LETKF_BUILD_GRIB "If the wgrib2 api is not found, it is downloaded and built" OFF)
IF(LETKF_ENABLE_GRIB)
  # need OpenMP
  find_package(OpenMP REQUIRED)

  IF(LETKF_BUILD_GRIB)
    # If we are going to build grib2 api ourselves in cmake
    set(GRIB2_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/dep/wgrib2-prefix/src/wgrib2/lib)
    set(GRIB2_INC_DIR ${GRIB2_BUILD_DIR})
    set(GRIB2_LIBRARY ${GRIB2_BUILD_DIR}/libwgrib2.a)
    set(GRIB2_API_LIBRARY ${GRIB2_BUILD_DIR}/libwgrib2_api.a)
  ELSE()
    # otherwise find the library
    find_path(GRIB2_INC_DIR wgrib2api.mod)
    find_library(GRIB2_API_LIBRARY libwgrib2_api.a)
    find_library(GRIB2_LIBRARY libwgrib2.a)
    IF(NOT GRIB2_INC_DIR OR NOT GRIB2_API_LIBRARY OR NOT GRIB2_LIBRARY)
      MESSAGE(FATAL_ERROR "cannot find wgrib2 API")
    ENDIF()
  ENDIF()
  mark_as_advanced(GRIB2_API_LIBRARY GRIB2_LIBRARY GRIB2_INC_DIR)
  set(GRIB2_LIBRARIES   ${GRIB2_API_LIBRARY} ${GRIB2_LIBRARY})
ENDIF()


#------------------------------------------------------------
# Set output directories
#------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)


#------------------------------------------------------------
# subdirectories
#------------------------------------------------------------

# included source of external dependencies
add_subdirectory(dep)

# main source code for library and default driver
add_subdirectory(src)

# tests
if(LETKF_ENABLE_TESTS)
  add_subdirectory(test)
endif()
