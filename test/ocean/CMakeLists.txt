# make sure the data is downloaded
letkf_data(
  ocean.testdata
  https://www.dropbox.com/s/4m5bpvb53vtsvlp/ocean.testdata.tgz )
letkf_data(
  ocean.ref_solutions
  https://www.dropbox.com/s/v95au84np83svas/ocean.ref_solutions.tgz )


# run a test for each of these .json configuration files
letkf_test(ocean obsio_nc)
letkf_test(ocean obsio_test)
letkf_test(ocean vtloc)
