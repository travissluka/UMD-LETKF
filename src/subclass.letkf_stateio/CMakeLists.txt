set(srcs
  letkf_state_nc.F90
  letkf_state_grib.F90
  )

foreach(s ${srcs})
  target_sources(letkf PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${s})
endforeach()
