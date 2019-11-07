# function to check and download test datasets
function(letkf_data localdir remotefile)
  if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${localdir}")
    message(STATUS "Downloading test data file for '${localdir}'...")
    file(
      DOWNLOAD ${remotefile}
      ${CMAKE_CURRENT_BINARY_DIR}/${localdir}.tgz )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E tar zxf ${localdir}.tgz
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endfunction()


set(LETKF_TEST_SH ${CMAKE_CURRENT_SOURCE_DIR}/letkf_test.sh)

# function that actually generates the tests
function(letkf_test name config)
  
  IF(TARGET check_norm_grib)
    SET(CHECK_GRIB_FILE $<TARGET_FILE:check_norm_grib>)
  ELSE()
    SET(CHECK_GRIB_FILE none)
  ENDIF()
    
  add_test(
    NAME ${name}.${config}
    COMMAND ${LETKF_TEST_SH}
    $<TARGET_FILE:check_norm_nc>
    ${CHECK_GRIB_FILE}
    $<TARGET_FILE:letkfdriver>
    ${name}
    ${config}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    )
endfunction()


# support executables for the tests
add_subdirectory(check_norm_nc)
IF(LETKF_ENABLE_GRIB)
  add_subdirectory(check_norm_grib)
ENDIF()

#--------------------------------------------------------------------------------

# test cases
add_subdirectory(simple3d)
add_subdirectory(ocean)
add_subdirectory(wave)
add_subdirectory(ioda)
