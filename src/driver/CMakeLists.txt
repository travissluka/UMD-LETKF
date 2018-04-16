project(letkf_driver Fortran)

set( EXE_NAME letkfdriver)

# # todo, allow selection of static vs shared


#------------------------------------------------------------
# required packages
#------------------------------------------------------------

# netcdf
include_directories(${NETCDF_INCLUDE_DIRS})

# mpi
include_directories(${MPI_Fortran_INCLUDE_PATH})


# ------------------------------------------------------------
# build target
# ------------------------------------------------------------
file( GLOB SRCS "*.f90")
add_executable(${EXE_NAME} ${SRCS})
target_link_libraries(${EXE_NAME} ${LETKF_LIBRARIES})
target_link_libraries(${EXE_NAME} ${NETCDF_LIBRARIES})
target_link_libraries(${EXE_NAME} ${LAPACK_LIBRARIES})
# target_link_libraries(${EXE_NAME} ${MPI_Fortran_LIBRARIES})
