# make sure the data is downloaded
letkf_data(
  ioda.testdata
  https://www.dropbox.com/s/lozi4karbjvsmxx/ioda.testdata.tgz )
letkf_data(
  ioda.ref_solutions
  https://www.dropbox.com/s/oamigw8o558mu54/ioda.ref_solutions.tgz )



# run a test for each of these .json configuration files
letkf_test(ioda default)
