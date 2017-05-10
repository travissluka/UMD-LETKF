project(letkf_module Fortran)

set( LIB_NAME letkf)

# external variables to export
# ------------------------------------------------------------
set(LETKF_LIBRARIES_SHARED ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${LIB_NAME}.so)
set(LETKF_LIBRARIES_STATIC ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/lib${LIB_NAME}.a)

set(LETKF_LIBRARIES_SHARED ${LETKF_LIBRARIES_SHARED})
set(LETKF_LIBRARIES_STATIC ${LETKF_LIBRARIES_STATIC})
set(LETKF_LIBRARIES ${LETKF_LIBRARIES_STATIC} PARENT_SCOPE)


# required packages
# ------------------------------------------------------------
include_directories(${NETCDF_INCLUDES})
include_directories(${MPI_Fortran_INCLUDE_PATH})


# build target
#------------------------------------------------------------
file( GLOB SRCS "*.f90")
#add_library( ${LIB_NAME}_shared SHARED ${SRCS})
#set_target_properties( ${LIB_NAME}_shared
#  PROPERTIES
#  OUTPUT_NAME ${LIB_NAME}
#)

add_library( ${LIB_NAME}_static STATIC ${SRCS})
set_target_properties( ${LIB_NAME}_static
  PROPERTIES
  OUTPUT_NAME ${LIB_NAME}
)
