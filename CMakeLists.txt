cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(letkf Fortran)


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


# external dependencies
#------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# MPI
find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})

# NetCDF
set(NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIRS})

# LAPACK / BLAS
#TODO, cmake not finding lapack with intel's MKL, but still compiles fine. fix this at some point
# so that blas and lapack are REQUIRED.
find_package(LAPACK)


# subdirectories
#------------------------------------------------------------

# external dependency source
add_subdirectory(dep)

# main source code for library and default driver
add_subdirectory(src)

# documentation
option ( LETKF_ENABLE_DOCS "Generate documentation")
if(LETKF_ENABLE_DOCS)
  add_subdirectory(doc)
endif()

# tests
if(LETKF_ENABLE_TESTS)
  add_subdirectory(test)
endif()
