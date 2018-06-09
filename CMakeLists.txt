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
find_package(LAPACK REQUIRED)


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


#------------------------------------------------------------
# main letkf library and driver program
#------------------------------------------------------------
add_subdirectory(src)


#------------------------------------------------------------
# testing
#------------------------------------------------------------
enable_testing()
add_subdirectory(test)


#------------------------------------------------------------
# documentation
#------------------------------------------------------------
# find_package(Doxygen)
# if (DOXYGEN_FOUND)
#   # set input and output files
#   set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
#   set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

#   # request to configure the file
#   configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
#   message("Doxygen build started")

#   # note the option ALL which allows to build the docs together with the application
#   add_custom_target( doc_doxygen ALL
#     COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
#     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#     COMMENT "Generating API documentation with Doxygen"
#     VERBATIM )
# else (DOXYGEN_FOUND)
#   message("Doxygen need to be installed to generate the doxygen documentation")
# endif (DOXYGEN_FOUND)
