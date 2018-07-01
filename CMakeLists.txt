cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(letkf Fortran)


# Testing and code coverage
#------------------------------------------------------------
option( ENABLE_TESTS "enable testing suite, including downloading of test data")
option( ENABLE_CODECOVERAGE "enable code coverage testing support")
if(ENABLE_CODECOVERAGE)
  add_definitions( -fprofile-arcs -ftest-coverage )
  set( CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} --coverage )
endif(ENABLE_CODECOVERAGE)


# external dependencies
#------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

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



# subdirectories
#------------------------------------------------------------
add_subdirectory(pkg)
add_subdirectory(src)
option ( ENABLE_DOCS "Generate documentation")
if(${ENABLE_DOCS})
  add_subdirectory(doc)
endif(${ENABLE_DOCS})


if( ${ENABLE_TESTS} )
  enable_testing()
  add_subdirectory(test)
endif(${ENABLE_TESTS})
