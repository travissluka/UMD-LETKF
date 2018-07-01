set(srcs check_norm_nc.F90)
set(exe check_norm_nc)

add_executable(${exe} ${srcs})
target_link_libraries( ${exe} ${NETCDF_LIBRARIES})

