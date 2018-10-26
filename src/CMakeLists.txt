#--------------------------------------------------------------------------------
# LETKF library
#--------------------------------------------------------------------------------
include(letkf_pluginloader.cmake)

# the core source files
set( letkf_lib_srcs
  getmem.F90
  timing.F90
  running_stats.F90
  ${JSON_SRC_FILES}
  ${KDTREE_SRC_FILES}
  letkf_config.F90
  letkf_core.F90
  letkf_loc.F90
  letkf_mpi.F90
  letkf_obs.F90
  letkf_solver.F90
  letkf_state.F90
  letkf_state_helper.F90
  letkf.F90
  )


# define the LETKF static linked library
add_library( letkf STATIC ${letkf_lib_srcs} )
target_link_libraries(letkf
  PRIVATE
    ${MPI_Fortran_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${LAPACK_LIBRARIES}
  )


# installation directory for the library
install(TARGETS letkf ARCHIVE DESTINATION lib)


# add the default localization, state, observation classes source code to the library
add_subdirectory(subclass.letkf_localizer)
add_subdirectory(subclass.letkf_stateio)
add_subdirectory(subclass.letkf_obsio)


# The modules that should be exported to the install directory
set( letkf_export_mods
  letkf.mod
  letkf_loc.mod
  letkf_mpi.mod
  letkf_obs.mod
  letkf_state.mod
  )
foreach (m ${letkf_export_mods})
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${m}
    DESTINATION include)
endforeach(m)



#--------------------------------------------------------------------------------
# Default LETKF driver program
#--------------------------------------------------------------------------------
add_executable( letkfdriver letkf_driver.F90)
target_link_libraries( letkfdriver letkf )
install(TARGETS letkfdriver RUNTIME DESTINATION bin)
