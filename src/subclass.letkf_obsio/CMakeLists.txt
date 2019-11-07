# stateio plugins that are always built
set(LETKF_PLUGINS
 letkf_obsio_nc
 letkf_obsio_dat
 letkf_obsio_test
 letkf_obsio_ioda
)

add_letkf_plugins("letkf_obsio" "letkf_obs" "${LETKF_PLUGINS}")
