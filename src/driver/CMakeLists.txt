project(letkf_driver Fortran)

set( EXE_NAME letkfdriver)

# required packages
#------------------------------------------------------------
include_directories(${NETCDF_INCLUDES})
include_directories(${MPI_Fortran_INCLUDE_PATH})

# todo, allow selection of static vs shared
set(LETKF_LIB_STATIC ${CMAKE_BINARY_DIR}/lib/libletkf.a)
set(LETKF_LIB_SHARED ${CMAKE_BINARY_DIR}/lib/libletkf.so)


# build target
# ------------------------------------------------------------
file( GLOB SRCS "*.f90")
add_executable(${EXE_NAME} ${SRCS})
target_link_libraries(${EXE_NAME} ${LETKF_LIB_STATIC})
target_link_libraries(${EXE_NAME} ${NETCDF_LIBRARIES_F77})
target_link_libraries(${EXE_NAME} ${LAPACK_LIBRARIES})
target_link_libraries(${EXE_NAME} ${MPI_Fortran_LIBRARIES})