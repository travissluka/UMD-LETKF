cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(letkf Fortran)


# set ouptut directory locations
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)



#------------------------------------------------------------
# external libraries
#------------------------------------------------------------
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# netcdf
set(NETCDF_F90 "YES")
set(NETCDF_USE_DEFAULT_PATHS "YES")
find_package(NetCDF REQUIRED)

# geoKdTree
set ( KDTREE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/geoKdTree CACHE FILEPATH "source of kdtree submodule" )
file ( GLOB KDTREE_SRC "${KDTREE_DIR}/*.[f,F]90")

# MPI
find_package(MPI REQUIRED)

# LAPACK / BLAS
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


# ------------------------------------------------------------
# global compiler options
# ------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
  "Choose the type of build, options are: Debug Release" FORCE)
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
#  add_definitions(-g)
  if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )
    add_definitions(-traceback -check all -ftrapuv -debug all)    
  elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )
    add_definitions(-Wall -fcheck=all)
  endif()
elseif(CMAKE_BUILD_TYPE MATCHES "Release")
#   add_definitions(-O3)
endif()



#------------------------------------------------------------
# determine the repository version and time
#------------------------------------------------------------
execute_process(COMMAND git describe --dirty --tags --always OUTPUT_VARIABLE
			GIT_VERSION RESULT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT ${RES} EQUAL "0")
  set(GIT_VERSION "Unknown")
endif()
add_definitions(-DCVERSION=\"${GIT_VERSION}\")
string(TIMESTAMP TS %Y-%m-%d )
add_definitions(-DCTIME=\"${TS}\")



# ------------------------------------------------------------
# subprojects
# ------------------------------------------------------------
add_subdirectory(module)
add_subdirectory(driver)
