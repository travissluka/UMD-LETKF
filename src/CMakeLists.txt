cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(letkf Fortran)


# global compiler options
# ------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
       set(CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: Debug Release" FORCE)
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_definitions(-Wall)
  add_definitions(-fcheck=all)
  add_definitions(-g)
else()
  add_definitions(-O3)
endif()


# determine the repository version and time 
execute_process(COMMAND git describe --dirty --tags --always OUTPUT_VARIABLE GIT_VERSION RESULT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT ${RES} EQUAL "0")
  set(GIT_VERSION "Unknown")
endif()
add_definitions(-DCVERSION=\"${GIT_VERSION}\")
string(TIMESTAMP TS %Y-%m-%d )
add_definitions(-DCTIME=\"${TS}\")


# other global variables
# ------------------------------------------------------------
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
include ( "cmake/pickFortranCompilerFlags.cmake" )
include ( "cmake/checkOutOfSource.cmake" )


# required packages
# ------------------------------------------------------------
find_package(LAPACK)
find_package(MPI REQUIRED)

set(NETCDF_INCLUDES $ENV{NETCDF_INCS})
set(NETCDF_LIBRARIES_F77 $ENV{NETCDF_LIBS})

# subprojects
# ------------------------------------------------------------
add_subdirectory(module)
add_subdirectory(driver)

