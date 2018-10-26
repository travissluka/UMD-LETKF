
# stateio plugins that are always built
set(LETKF_PLUGINS
 letkf_stateio_nc
)



# GRIB stateio plugin is optional
#-------------------------------------------------------------------------------
option( LETKF_ENABLE_GRIB "enable grib2 support for state I/O" 1)
#TODO this is temporary for working with grib API, need to clean up
set(GRIB_INC_DIR  /home/tsluka/proj/UMD-LETKF/grib2/ftn_api)
set(GRIB_LIB /home/tsluka/proj/UMD-LETKF/grib2/ftn_api/wgrib2api.o)
IF(LETKF_ENABLE_GRIB)
  LIST(APPEND LETKF_PLUGINS letkf_stateio_grib)
  include_directories(${GRIB_INC_DIR})
  target_sources(letkf PRIVATE ${GRIB_LIB})
ENDIF()


add_letkf_plugins("letkf_stateio" "letkf_state" "${LETKF_PLUGINS}")
