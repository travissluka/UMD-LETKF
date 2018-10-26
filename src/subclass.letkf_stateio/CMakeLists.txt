
# stateio plugins that are always built
set(LETKF_PLUGINS
 letkf_stateio_nc
)

# GRIB stateio plugin is optional
IF(LETKF_ENABLE_GRIB)
  LIST(APPEND LETKF_PLUGINS letkf_stateio_grib)
ENDIF()


add_letkf_plugins("letkf_stateio" "letkf_state" "${LETKF_PLUGINS}")
