set( letkf_lib_srcs
  getmem.F90
  timing.F90
  running_stats.F90
  geoKdTree/kdtree.f90
  
  
  letkf.F90
  letkf_core.F90
  letkf_loc.F90
  letkf_loc_novrt.F90
  letkf_mpi.F90
  letkf_obs.F90
  letkf_obs_nc.F90
  letkf_obs_test.F90
  letkf_solver.F90
  letkf_state.F90
  letkf_state_nc.F90
  )

add_library( letkf STATIC ${letkf_lib_srcs} )



set (letkf_driver_srcs
  letkf_driver.F90
  )

add_executable( letkfdriver ${letkf_driver_srcs})
target_link_libraries( letkfdriver letkf )
target_link_libraries( letkfdriver ${MPI_Fortran_LIBRARIES})
target_link_libraries( letkfdriver ${NETCDF_LIBRARIES})
target_link_libraries( letkfdriver ${LAPACK_LIBRARIES})
