file(GLOB input FOLLOW_SYMLINKS config/*)
foreach(F ${input})
  get_filename_component(F2 ${F} NAME)
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${F} ${CMAKE_CURRENT_BINARY_DIR}/${F2})
endforeach(F)

add_test(
  NAME ocean_mom6_1deg_obstest
  COMMAND letkfdriver letkf.json
)

#add_test(
#  NAME ocean_mom6_1deg_obsnc
#  COMMAND letkfdriver letkf.json
#)

#add_test(
#  NAME ocean_mom6_1deg
#  COMMAND letkfdriver full.namelist
#)

