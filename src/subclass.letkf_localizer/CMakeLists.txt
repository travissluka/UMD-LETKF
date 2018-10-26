set(LETKF_PLUGINS
  letkf_loc_novrt
  letkf_loc_ocean
  )

add_letkf_plugins("letkf_localizer" "letkf_loc" "${LETKF_PLUGINS}")
