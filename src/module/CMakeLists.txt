project(letkf_module Fortran)

set( LIB_NAME letkf)


# source files for this module
#------------------------------------------------------------
file( GLOB SRCS "*.[f,F]90")


# required libraries
#------------------------------------------------------------

# netcdf
include_directories(${NETCDF_INCLUDE_DIRS})

# mpi
#include_directories(${MPI_Fortran_INCLUDE_PATH})



# static linked library
#------------------------------------------------------------
add_library( ${LIB_NAME}_static STATIC ${SRCS} ${KDTREE_SRC})
set_target_properties( ${LIB_NAME}_static
  PROPERTIES
  OUTPUT_NAME ${LIB_NAME}
)

# # shared library
# #------------------------------------------------------------
# add_library( ${LIB_NAME}_shared SHARED ${SRCS} ${KDTREE_SRC})
# set_target_properties( ${LIB_NAME}_shared
#   PROPERTIES
#   OUTPUT_NAME ${LIB_NAME}
# )



# external variables to export
# ------------------------------------------------------------
# set(LETKF_LIBRARIES_SHARED ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${LIB_NAME}.so)
set(LETKF_LIBRARIES_STATIC ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/lib${LIB_NAME}.a)

set(LETKF_LIBRARIES ${LETKF_LIBRARIES_STATIC} PARENT_SCOPE)




