cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

## TODO, shouldn't need this!!
#include_directories( "/usr/lib64/gfortran/modules")


# global compiler options
# ------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
       set(CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build, options are: Debug Release" FORCE)
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_definitions(-Wall)
  add_definitions(-fcheck=all)
  add_definitions(-g)
else()
  add_definitions(-O3)
endif()


# other global variables
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)
enable_language( Fortran )
include ( "../cmake/pickFortranCompilerFlags.cmake" )
include ( "../cmake/checkOutOfSource.cmake" )


## setup project files, dependencies, install locations
## ------------------------------------------------------------
project(letkf)

## output file directory
set ( MODULE_DIR "${CMAKE_BINARY_DIR}/include")
set ( CMAKE_Fortran_MODULE_DIRECTORY ${MODULE_DIR} )
set ( INSTALL_LIB_DIR ${INSTALL_MOD_DIR} )
set ( INSTALL_MOD_DIR ${PACKAGE_NAME}-foo/lib )


## required dependencies
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)

find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})

find_package(LAPACK REQUIRED)

set(NETCDF_F77 "YES")
find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDES})
#find_package(OpenMP_Fortran)


## build shared library
## ------------------------------------------------------------
file( GLOB LETKF_LIB_SRCS  "module/*.f90" )
set( LIB_NAME ${CMAKE_PROJECT_NAME})
add_library( ${LIB_NAME} SHARED ${LETKF_LIB_SRCS} )
set_target_properties( ${LIB_NAME}
  PROPERTIES
  PREFIX lib
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)


## build driver
## ------------------------------------------------------------
file( GLOB LETKF_DRIV_SRCS "driver/*.f90" )
set( DRIVER_NAME letkfdriver)
add_executable( ${DRIVER_NAME} ${LETKF_DRIV_SRCS})
target_link_libraries( ${DRIVER_NAME} ${LIB_NAME} )
target_link_libraries( ${DRIVER_NAME} ${NETCDF_LIBRARIES_F77})
target_link_libraries( ${DRIVER_NAME} ${MPI_Fortran_LIBRARIES} )
target_link_libraries( ${DRIVER_NAME} ${LAPACK_LIBRARIES})
