# wave DA tests need grib IO
IF(LETKF_ENABLE_GRIB)

  # make sure test data is downloaded
  letkf_data(
    wave.testdata
    https://www.dropbox.com/s/1l4rplpnl49tp6d/wave.testdata.tgz)
  letkf_data(
    wave.ref_solutions
    https://www.dropbox.com/s/lzb5kmorpwhkcqa/wave.ref_solutions.tgz)

  # run a test for each of these .json configuration files
  letkf_test(wave default)

ENDIF()
