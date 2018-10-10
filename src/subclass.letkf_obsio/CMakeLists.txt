set(srcs
  letkf_obs_nc.F90
  letkf_obs_dat.F90
  letkf_obs_test.F90
  )


foreach(s ${srcs})
  target_sources(letkf PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${s})
endforeach()
