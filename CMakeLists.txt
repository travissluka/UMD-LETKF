cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(letkf Fortran)

#------------------------------------------------------------
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


#------------------------------------------------------------
# determine the repository version and time
#------------------------------------------------------------
execute_process(COMMAND git describe --dirty --tags --always OUTPUT_VARIABLE
			GIT_VERSION RESULT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT ${RES} EQUAL "0")
  set(GIT_VERSION "Unknown")
endif()
add_definitions(-DCVERSION=\"${GIT_VERSION}\")
#string(TIMESTAMP TS %Y-%m-%d )
#add_definitions(-DCTIME=\"${TS}\")


# subdirectories
#------------------------------------------------------------
enable_testing()
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(doc)
