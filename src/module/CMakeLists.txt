project(letkf_module Fortran)

set( LIB_NAME letkf)

# required packages
# ------------------------------------------------------------
include_directories(${NETCDF_INCLUDES})
include_directories(${MPI_Fortran_INCLUDE_PATH})
message(${MPI_Fortran_INCLUDE_PATH})


# build target
#------------------------------------------------------------
file( GLOB SRCS "*.f90")
add_library( ${LIB_NAME}_shared SHARED ${SRCS})
set_target_properties( ${LIB_NAME}_shared
  PROPERTIES
  OUTPUT_NAME ${LIB_NAME}
)

add_library( ${LIB_NAME}_static STATIC ${SRCS})
set_target_properties( ${LIB_NAME}_static
  PROPERTIES
  OUTPUT_NAME ${LIB_NAME}
)