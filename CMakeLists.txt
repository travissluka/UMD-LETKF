cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

## TODO, shouldn't need this!!
include_directories( "/usr/lib64/gfortran/modules")


## set the configuration of build to preform
## ------------------------------------------------------------
#set (CMAKE_CONFIGURATION_TYPES "Debug" "Release")
#set (CMAKE_BUILD_TYPE "Release"
#  CACHE STRING "Select which configuration to build.")
#set_property ( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_CONFIGURATION_TYPES} )

enable_language( Fortran )
include ( "cmake/pickFortranCompilerFlags.cmake" )
include ( "cmake/checkOutOfSource.cmake" )



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
find_package(MPI)
include_directories(${MPI_INCLUDE})
find_package(LAPACK REQUIRED)
set(NETCDF_F77 "YES")
find_package(NetCDF REQUIRED)
find_package(OpenMP_Fortran)


## build shared library
file( GLOB LETKF_LIB_SRCS  "src/*.f90" )
set( LIB_NAME ${CMAKE_PROJECT_NAME})
add_library( ${LIB_NAME} SHARED ${LETKF_LIB_SRCS} )
set_target_properties( ${LIB_NAME}
  PROPERTIES
  PREFIX lib
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)


## build test executables
enable_testing()
file( GLOB LETKF_TEST_SRCS "src/test/*.f90" )
foreach( test_source_file ${LETKF_TEST_SRCS} )
  string( REPLACE ".f90" "" test_name ${test_source_file} )
  string( FIND  ${test_name} "/" pos REVERSE )
  math( EXPR pos "${pos}+1" )
  string( SUBSTRING ${test_name} ${pos} -1 test_name )
  add_executable( ${test_name} ${test_source_file} )
  target_link_libraries( ${test_name} ${LIB_NAME} )
  target_link_libraries( ${test_name} ${LAPACK_LIBRARIES} )
  add_test(${test_name} ${test_name})
endforeach( test_source_file ${LETKF_TEST_SRCS} )
