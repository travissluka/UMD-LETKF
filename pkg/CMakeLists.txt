set(JSON_SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/json-fortran/src/json_file_module.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/json-fortran/src/json_kinds.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/json-fortran/src/json_module.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/json-fortran/src/json_parameters.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/json-fortran/src/json_string_utilities.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/json-fortran/src/json_value_module.F90
  PARENT_SCOPE)

set(KDTREE_SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/geoKdTree/kdtree.f90
  PARENT_SCOPE)
  
