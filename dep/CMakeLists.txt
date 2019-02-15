set(KDTREE_SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/geoKdTree/kdtree.f90
  PARENT_SCOPE)

IF(LETKF_ENABLE_GRIB)
  include(wgrib2.cmake)
ENDIF()

option(BUILD_TESTING "" OFF)
add_subdirectory(libyaml)
