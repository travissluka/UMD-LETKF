set(PLUGINLOADER_TMPL ${CMAKE_CURRENT_SOURCE_DIR}/letkf_pluginloader.F90.tmpl)

#-------------------------------------------------------------------------------
# This automatically generates a source code file
# "letkf_pluginloader_${PLUGIN_TYPE}.F90" that will be added to the build with
#  the purpose of registering each plugin with the main program.
#-------------------------------------------------------------------------------
function( add_letkf_plugins PLUGIN_TYPE PLUGIN_MOD PLUGIN_SRCS)

  set(PLUGIN_MODS "")
  set(PLUGIN_ALLOCS "")
  foreach(p ${PLUGIN_SRCS})
      target_sources(letkf PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${p}.F90)
      STRING(CONCAT PLUGIN_MODS ${PLUGIN_MODS}
        "USE ${p}_mod\n")
      STRING(CONCAT PLUGIN_ALLOCS ${PLUGIN_ALLOCS}
       "ALLOCATE(${p}::ptr)\n"
       "CALL ${PLUGIN_MOD}_register(ptr)\n")
  endforeach()

  # generate the source file, and add it to the build
  configure_file(${PLUGINLOADER_TMPL} letkf_pluginloader_${PLUGIN_TYPE}.F90)
  target_sources(letkf PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/letkf_pluginloader_${PLUGIN_TYPE}.F90)

endfunction()
