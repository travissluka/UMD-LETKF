# Copyright 2016-2019 Travis Sluka
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
  letkf_util.F90
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

# yaml
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
