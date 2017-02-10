cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(letkf Fortran)


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



# ## TODO, shouldn't need this!!
# #include_directories( "/usr/lib64/gfortran/modules")



# other global variables
# ------------------------------------------------------------
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
include ( "../cmake/pickFortranCompilerFlags.cmake" )
include ( "../cmake/checkOutOfSource.cmake" )


# required packages
# ------------------------------------------------------------
set(NETCDF_F77 "YES")
find_package(NetCDF REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)


# subprojects
# ------------------------------------------------------------
add_subdirectory(module)
add_subdirectory(driver)

