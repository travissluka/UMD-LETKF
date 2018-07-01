file(GLOB input FOLLOW_SYMLINKS *)
foreach(F ${input})
  get_filename_component(F2 ${F} NAME)
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${F} ${CMAKE_CURRENT_BINARY_DIR}/${F2})
endforeach(F)

add_test(
  NAME simple3d
  COMMAND letkfdriver default.json
  )

add_test(
  NAME simple3d.inflation.rtps
  COMMAND letkfdriver inflation_rtps.json
)

add_test(
  NAME simple3d.inflation.rtpp
  COMMAND letkfdriver inflation_rtpp.json
)

add_test(
  NAME simple3d.inflation.mul
  COMMAND letkfdriver inflation_mul.json
)
