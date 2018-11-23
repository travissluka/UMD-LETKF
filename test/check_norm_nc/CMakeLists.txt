# source files
set(srcs check_norm_nc.F90)
set(exe check_norm_nc)

# executable program
add_executable(${exe} ${srcs})

# dependencies
include_directories(${NETCDF_INCLUDE_DIRS})
target_link_libraries( ${exe} ${NETCDF_LIBRARIES})
