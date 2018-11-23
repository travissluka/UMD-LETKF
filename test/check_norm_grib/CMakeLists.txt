# source files
set(srcs check_norm_grib.F90)
set(exe check_norm_grib)

# executable program
add_executable(${exe} ${srcs})

# dependencies
IF(LETKF_ENABLE_GRIB)
  include_directories(${GRIB2_INC_DIR})
  target_link_libraries(${exe}  PRIVATE  ${GRIB2_LIBRARIES})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
ENDIF()
