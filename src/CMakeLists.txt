# All of the source files
set( letkf_lib_srcs
  getmem.F90
  timing.F90
  running_stats.F90
   
  letkf.F90
  letkf_config.F90
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

  ${JSON_SRC_FILES}
  ${KDTREE_SRC_FILES}
  )

# The modules that should be exported in the install directory
set( export_mods
  letkf.mod
  letkf_loc.mod
  letkf_mpi.mod
  letkf_obs.mod
  letkf_state.mod
  )


# static linked library
add_library( letkf STATIC ${letkf_lib_srcs} )


# main driver program
set (letkf_driver_srcs
  letkf_driver.F90
  )
add_executable( letkfdriver ${letkf_driver_srcs})
target_link_libraries( letkfdriver letkf )
target_link_libraries( letkfdriver ${MPI_Fortran_LIBRARIES})
target_link_libraries( letkfdriver ${NETCDF_LIBRARIES})
target_link_libraries( letkfdriver ${LAPACK_LIBRARIES})


# install directory
install(TARGETS letkfdriver letkf
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib)

foreach (m ${export_mods})
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${m}
    DESTINATION include)
endforeach(m)
