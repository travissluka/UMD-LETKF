set(srcs
  letkf_loc_novrt.F90
  letkf_loc_ocean.F90
  )

foreach(s ${srcs})
  target_sources(letkf PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${s})
endforeach()
