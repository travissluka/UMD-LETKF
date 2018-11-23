set(WGRIB2_URL http://www.ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/wgrib2.tgz)
set(WGRIB2_TGZ ${CMAKE_CURRENT_BINARY_DIR}/wgrib2.tgz)

IF(LETKF_BUILD_GRIB)
  include(ExternalProject)

  # make sure the tgz file is downloaded
  IF(NOT EXISTS ${WGRIB2_TGZ})
    message(STATUS "Downloading NCEP wgrib2 source code...")
    FILE (DOWNLOAD ${WGRIB2_URL} ${WGRIB2_TGZ})
  ENDIF()

  # details for compiling the wgrib2 library
  ExternalProject_Add(wgrib2
    URL ${WGRIB2_TGZ}
    URL_MD5 fd08ac8b988c310c125b50550aebd9f2
    CONFIGURE_COMMAND ""
    PATCH_COMMAND patch -p1 < ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2.patch
    BUILD_COMMAND make lib FC=gfortran CC=gcc
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
  )

ENDIF()
