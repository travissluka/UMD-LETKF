# make sure the data is downloaded
letkf_data(
  mom6_1deg.testdata
  https://www.dropbox.com/s/4m5bpvb53vtsvlp/mom6_1deg.testdata.tgz )
letkf_data(
  mom6_1deg.ref_solutions
  https://www.dropbox.com/s/v95au84np83svas/mom6_1deg.ref_solutions.tgz )


# run a test for each of these .json configuration files
letkf_test(mom6_1deg obsio_nc)
letkf_test(mom6_1deg obsio_test)
letkf_test(mom6_1deg vtloc)
