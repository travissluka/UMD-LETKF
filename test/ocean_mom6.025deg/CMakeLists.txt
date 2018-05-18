file(GLOB_RECURSE input FOLLOW_SYMLINKS config/*)
foreach(F ${input})
  get_filename_component(F2 ${F} NAME)
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${F} ${CMAKE_CURRENT_BINARY_DIR}/${F2})
endforeach(F)

# add_test(
#   NAME test_ocean_mom6_025deg
#   COMMAND letkfdriver
# )
