#--------------------------------------------------------------------------------
# LETKF library
#--------------------------------------------------------------------------------
include(letkf_pluginloader.cmake)


# the core source files
set( letkf_lib_srcs
  getmem.F90
  timing.F90
  running_stats.F90
  letkf_config.F90
  letkf_config_c.cpp
  letkf_core.F90
  letkf_loc.F90
  letkf_mpi.F90
  letkf_obs.F90
  letkf_solver.F90
  letkf_state.F90
  letkf_state_helper.F90
  letkf.F90
  ${KDTREE_SRC_FILES}
  )
add_library( letkf STATIC ${letkf_lib_srcs} )


# add the default plugins source code to the library
# (localization, state, observation classes)
add_subdirectory(subclass.letkf_localizer)
add_subdirectory(subclass.letkf_stateio)
add_subdirectory(subclass.letkf_obsio)


# external dependencies
#-------------------------------------------------------------------------------
# MPI
include_directories(${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(letkf  PRIVATE  ${MPI_Fortran_LIBRARIES})

#NetCDF
include_directories(${NETCDF_INCLUDE_DIRS})
target_link_libraries(letkf  PRIVATE  ${NETCDF_LIBRARIES})

# LAPACK
target_link_libraries(letkf  PRIVATE ${LAPACK_LIBRARIES})

# GRIB2
IF(LETKF_ENABLE_GRIB)
  include_directories(${GRIB2_INC_DIR})
  target_link_libraries(letkf  PRIVATE  ${GRIB2_LIBRARIES})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_Fortran_FLAGS}")
ENDIF()

# yaml, TODO: set the include directory correctly
include_directories(${CMAKE_SOURCE_DIR}/dep/libyaml/include)
target_link_libraries(letkf PRIVATE yaml)


# installation directory for the library
#-------------------------------------------------------------------------------
install(TARGETS letkf ARCHIVE DESTINATION lib)

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
set_target_properties(letkfdriver PROPERTIES LINKER_LANGUAGE Fortran) 
target_link_libraries( letkfdriver letkf )
install(TARGETS letkfdriver RUNTIME DESTINATION bin)
