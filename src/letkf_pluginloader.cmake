# Copyright 2018-2019 Travis Sluka
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
#-------------------------------------------------------------------------------

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
