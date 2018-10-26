# make sure test data is downloaded
#letkf_data(
#  simple3d.testdata
#  https://www.dropbox.com/s/8t9thu4om0iichv/simple3d.testdata.tgz)
#letkf_data(
  #simple3d.ref_solutions
#https://www.dropbox.com/s/v6uier1xzjbnsav/simple3d.ref_solutions.tgz)


# run a test for each of these .json configuration files
letkf_test(wave default)
