! Copyright 2019-2019 Travis Sluka
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!================================================================================
!>
!================================================================================
MODULE letkf_obsio_ioda_mod
  USE netcdf
  USE letkf_config
  USE letkf_obs
  USE letkf_mpi

  IMPLICIT NONE
  PRIVATE

  !================================================================================
  !================================================================================
  ! Public module components
  !================================================================================
  !================================================================================

  INTEGER, PARAMETER :: MAX_FILENAME_LEN = 1024

  !================================================================================
  !> observation file I/O class for handling JEDI IODA NetCDF files
  !--------------------------------------------------------------------------------
  TYPE, EXTENDS(letkf_obsio), PUBLIC :: letkf_obsio_ioda
     TYPE(configuration) :: obsfiles_config

   CONTAINS
     PROCEDURE, NOPASS :: name => obsio_ioda_get_name
     PROCEDURE, NOPASS :: desc => obsio_ioda_get_desc
     PROCEDURE         :: init => obsio_ioda_init
     PROCEDURE         :: read_obs => obsio_ioda_read_obs
     PROCEDURE         :: read_hx  => obsio_ioda_read_hx
  END TYPE letkf_obsio_ioda
  !================================================================================





CONTAINS



  !================================================================================
  !> Get the unique name of the observation I/O class
  !! @param name short (less than ~10 character) class name
  !--------------------------------------------------------------------------------
  FUNCTION obsio_ioda_get_name() RESULT(name)
    CHARACTER(:), ALLOCATABLE :: name
    name = "OBSIO_IODA"
  END FUNCTION obsio_ioda_get_name
  !================================================================================



  !================================================================================
  !> Get a human readable description of this observation I/O class
  !! @param desc description of the class
  !--------------------------------------------------------------------------------
  FUNCTION obsio_ioda_get_desc() RESULT(desc)
    CHARACTER(:), ALLOCATABLE :: desc
    desc = "NetCDF formatted JEDI-IODA observation I/O"
  END FUNCTION obsio_ioda_get_desc
  !================================================================================



  !================================================================================
  !> Initialize the netcdf obsio reader/writer.
  !! This loads in a section of the namelist (&letkf_obs_nc) to determine what the
  !! input filenames should be.
  !! @param nml_filename path to the namelist we are to open
  !--------------------------------------------------------------------------------
  SUBROUTINE obsio_ioda_init(self, config, obsdef, platdef)
    CLASS(letkf_obsio_ioda) :: self
    TYPE(configuration), INTENT(in) :: config
    TYPE(letkf_obsplatdef_list), INTENT(out) :: obsdef
    TYPE(letkf_obsplatdef_list), INTENT(out) :: platdef

    TYPE(configuration) :: config_tmp1, config_tmp2, config_tmp3
    INTEGER :: n, i, nvars
    CHARACTER(len=:), ALLOCATABLE :: str

    TYPE(letkf_obsplatdef) :: obsdef_item

    ! get the files configurations
    CALL config%get("ioda_files", self%obsfiles_config)

    ! generate the "obs_def section"
    nvars=0
    DO n=1, self%obsfiles_config%count()
      CALL self%obsfiles_config%get(n, config_tmp1)
      CALL config_tmp1%get("vars", config_tmp2)
      DO i=1, config_tmp2%count()
        CALL config_tmp2%get(i, config_tmp3)
        CALL config_tmp3%get(1, str)
        IF (.NOT. obsdef%has(str)) nvars = nvars+1
        obsdef_item%id = nvars
        obsdef_item%name = TRIM(str)
        CALL config_tmp3%get(3, str)
        obsdef_item%name_long = TRIM(str)
        CALL obsdef%add(obsdef_item, allow_duplicate=.TRUE.)
      END DO
    END DO

    ! generate the platdef section
    nvars=0
    DO n=1, self%obsfiles_config%count()
      CALL self%obsfiles_config%get(n, config_tmp1)
      CALL config_tmp1%get("vars", config_tmp2)
      DO i=1, config_tmp2%count()
        CALL config_tmp2%get(i, config_tmp3)
        CALL config_tmp3%get(2, str)
        IF (.NOT. platdef%has(str)) nvars = nvars+1
        obsdef_item%id = nvars
        obsdef_item%name = TRIM(str)
        obsdef_item%name_long = ""
        CALL platdef%add(obsdef_item, allow_duplicate=.TRUE.)
      END DO
    END DO

  END SUBROUTINE obsio_ioda_init
  !================================================================================



  !================================================================================
  !> Read the observations from a netcdf file.
  !! Note: this is the observation only (location, value, etc.) and NOT the
  !! per-ensemble observation departure
  !! @param obs the list of observations as they are loaded in.
  !--------------------------------------------------------------------------------
  SUBROUTINE obsio_ioda_read_obs(self, obs)
    CLASS(letkf_obsio_ioda) :: self
    TYPE(letkf_observation), ALLOCATABLE, INTENT(out) :: obs(:)

    INTEGER :: count_local, count, file_num

    INTEGER :: obsid, platid
    INTEGER :: ncid, dimid
    INTEGER :: nlocs, offset, nvars
    INTEGER :: i, n, num_files
    LOGICAL :: ex
    TYPE(configuration) :: local_config, vars_config, var_config
    TYPE(letkf_observation), ALLOCATABLE :: obs_tmp(:)
    TYPE(letkf_obsplatdef) :: obsdef

    CHARACTER(len=:), ALLOCATABLE :: str1, var_short, var_long
    CHARACTER(len=MAX_FILENAME_LEN) :: str2

    ! Determine the number of obs that need to be allocated
    count = 0
    num_files = self%obsfiles_config%count()
    do n = 1, num_files
      count_local = 0

      ! get filename base to read
      CALL self%obsfiles_config%get(n, local_config)
      CALL local_config%get("file", str1)
      CALL str_pattern_ens(str1, str1, 1)

      ! get the number of variables (according to the yaml file)
      CALL local_config%get("vars", var_config)
      nvars = var_config%count()

      ! there are likely multiple files that end with a 4 digit number
      ! Check to see which ones exist and parse them
      file_num = 0
      DO
        WRITE(str2, '(A,I0.4,A)') str1//"_", file_num, ".nc"

        ! make sure the file exists
        INQUIRE(file=str2, exist=ex)
        IF (.NOT. ex) EXIT

        ! file exists, get the number of observations in the file
        CALL check( nf90_open(str2, nf90_nowrite, ncid))
        CALL check( nf90_inq_dimid(ncid, "nlocs", dimid), &
            'Reading dimension "nlocs"')
        CALL check( nf90_inquire_dimension(ncid, dimid, len=nlocs) )
        CALL check( nf90_close(ncid))
        count_local = count_local + (nlocs * nvars)

        ! find next file...
        file_num = file_num + 1
      END DO

      count = count + count_local
    end do
    PRINT *, "total obs count: ", count


    ! allocate space depending on global number of observations
    ALLOCATE(obs(count))

    ! read in each file, and append to the observation list
    offset = 0
    do n = 1, num_files
      ! get filename base to read
      CALL self%obsfiles_config%get(n, local_config)
      CALL local_config%get("file", str1)
      CALL str_pattern_ens(str1, str1, 1)

      ! there are likely multiple files that end with a 4 digit number
      ! Check to see which ones exist and parse them
      file_num = 0
      DO
        WRITE(str2, '(A,I0.4,A)') str1//"_", file_num, ".nc"

        ! make sure the file exists
        INQUIRE(file=str2, exist=ex)
        IF (.NOT. ex) EXIT
        !IF (pe_isroot) PRINT *, "  reading: ", TRIM(str2)

        ! for each variable listed in the yaml for this file
        CALL local_config%get("vars", vars_config)
        DO i = 1, vars_config%count()
          CALL vars_config%get(i, var_config)
          CALL var_config%get(1, var_short)
          obsdef=letkf_obs_getdef('O', var_short)
          obsid=obsdef%id
          CALL var_config%get(2, var_short)
          obsdef=letkf_obs_getdef('P', var_short)
          platid=obsdef%id
          CALL var_config%get(3, var_long)

          CALL read_file(str2, obs_tmp, obsid, platid, var_long)
          IF (.NOT. ALLOCATED(obs_tmp)) &
            CALL letkf_mpi_abort("ERROR loading file")

          ! cat array
          obs(offset+1:offset+1+SIZE(obs_tmp)) = obs_tmp
          offset = offset + SIZE(obs_tmp)
        END DO

        ! find next file...
        file_num = file_num + 1
      END DO

    end do

  END SUBROUTINE obsio_ioda_read_obs
  !================================================================================



  !================================================================================
  !> Read the observation operator output for a single ensemble member.
  !! @param ensmem ensemble member number to load in
  !! @param hx the ensemble member state mapped to observation space
  !--------------------------------------------------------------------------------
  SUBROUTINE obsio_ioda_read_hx(self, ensmem, hx)
    CLASS(letkf_obsio_ioda) :: self
    INTEGER, INTENT(in) :: ensmem
    REAL, ALLOCATABLE, INTENT(inout) :: hx(:)

    INTEGER :: n,i
    LOGICAL :: ex
    CHARACTER(len=:), ALLOCATABLE :: str1, var_long
    INTEGER :: offset, file_num
    TYPE(configuration) :: local_config, vars_config, var_config
    CHARACTER(len=MAX_FILENAME_LEN) :: str2

    REAL, ALLOCATABLE :: hx_tmp(:)

    ! read in each file, and append to the hofx list
    offset = 0
    DO n=1, self%obsfiles_config%count()
       ! get filename base to read
       CALL self%obsfiles_config%get(n, local_config)
       CALL local_config%get("file", str1)
       CALL str_pattern_ens(str1, str1, ensmem)

       ! there are likely multiple files that end with a 4 digit number
       ! check to see which ones exist and parse them
       file_num = 0
       DO
          WRITE(str2, '(A,I0.4,A)') str1 // "_", file_num, ".nc"

          ! make sure file exists
          ! TODO check verbose
          INQUIRE(file=str2, exist=ex)
          IF (.NOT. ex) EXIT
          !PRINT *, "reading hofx: ", TRIM(str2)

          ! for each variable in the yaml for this file
          CALL local_config%get("vars", vars_config)
          DO i=1,vars_config%count()
            CALL vars_config%get(1, var_config)
            CALL var_config%get(3, var_long)
            CALL read_file_hx(str2, hx_tmp, var_long)

            IF (.NOT. ALLOCATED(hx_tmp)) CALL letkf_mpi_abort( &
               "ERROR loading hx file")

            ! safety check to make sure we won't go out of bounds
            IF (offset + SIZE(hx_tmp) > SIZE(hx)) CALL letkf_mpi_abort( &
              "ERROR in ioda_read_hx(), hx array out of bounds")

            ! cat array
            hx(offset+1:offset+SIZE(hx_tmp)) = hx_tmp
            offset = offset + SIZE(hx_tmp)
          END DO

          ! find next file
          file_num = file_num + 1
       END DO
    END DO
  END SUBROUTINE obsio_ioda_read_hx
  !================================================================================


  !================================================================================
  SUBROUTINE read_file_hx(filename, hx, varname)
    CHARACTER(*), INTENT(in) :: filename, varname
    REAL, ALLOCATABLE, INTENT(OUT) :: hx(:)

    INTEGER :: ncid, dimid, nlocs, varid
    LOGICAL :: ex
    CHARACTER(len=:), ALLOCATABLE :: varname2

    ! make sure the file exists
    INQUIRE(file=filename, exist=ex)
    IF( .NOT. ex) CALL letkf_mpi_abort( &
         "observation hx file does not exist: " &
         // TRIM(filename))

    ! read the dimension sizes
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid), &
         "Openning "//TRIM(filename) )
    CALL check( nf90_inq_dimid(ncid, "nlocs", dimid), &
         '"nlocs" in '//TRIM(filename))
    CALL check( nf90_inquire_dimension(ncid, dimid, len=nlocs))
    ALLOCATE( hx(nlocs))

    ! read in the data
    varname2=varname//"@hofx"
    CALL check( nf90_inq_varid(ncid, varname2, varid), &
         varname2 // ' in '//filename )
    CALL check( nf90_get_var(ncid, varid, hx) )

    ! done, close
    CALL check( nf90_close(ncid))
  END SUBROUTINE read_file_hx
  !================================================================================


  !================================================================================
  !--------------------------------------------------------------------------------
  SUBROUTINE read_file(filename, obs, obsid, platid, varname)
    CHARACTER(*), INTENT(in) :: filename, varname
    INTEGER, INTENT(in) :: obsid, platid
    TYPE(letkf_observation), ALLOCATABLE, INTENT(out) :: obs(:)

    INTEGER :: n
    LOGICAL :: ex
    INTEGER :: dimid, ncid, nlocs, vid

    INTEGER, ALLOCATABLE :: tmp_i(:)
    REAL(4), ALLOCATABLE :: tmp_r(:)

    CHARACTER(len=:), ALLOCATABLE :: varname2

    ! need to read in, for each variable
    ! variable_names@VarMetaData
    ! @ObsValue, @ObsError, @PreQC
    ! latitude@MetaData, longitude@MetaData, depth@MetaData

    ! make sure the file exists
    INQUIRE(file=filename, exist=ex)
    IF (.NOT. ex) CALL letkf_mpi_abort(&
        "File not found: "//filename)

    ! open the file, get number of variables and size of observations
    CALL check( nf90_open(filename, nf90_nowrite, ncid))
    CALL check( nf90_inq_dimid(ncid, "nlocs", dimid), &
         'Reading dimension "nlocs"')
    CALL check( nf90_inquire_dimension(ncid, dimid, len=nlocs) )
    !CALL check( nf90_inq_dimid(ncid, "nvars", dimid), &
    !     'Reading dimension "nvars"')
    !CALL check( nf90_inquire_dimension(ncid, dimid, len=nvars) )

    ALLOCATE(obs(nlocs))
    ALLOCATE(tmp_r(nlocs))
    ALLOCATE(tmp_i(nlocs))
    DO n=1,nlocs
       obs(n)%obsid = obsid
       obs(n)%platid = platid
    END DO

    ! read in the variables
    CALL check( nf90_inq_varid(ncid, "latitude@MetaData", vid),&
         '"lat" in '//filename)
    CALL check( nf90_get_var(ncid, vid, tmp_r) )
    DO n=1,nlocs
        obs(n)%lat = tmp_r(n)
    END DO

    CALL check( nf90_inq_varid(ncid, "longitude@MetaData", vid),&
         '"lon" in '//filename )
    CALL check( nf90_get_var(ncid, vid, tmp_r) )
    DO n=1,nlocs
       obs(n)%lon = tmp_r(n)
    END DO

    ! TODO, check if z coord is given
    ! CALL check( nf90_inq_varid(ncid, "depth", varid),&
    DO n=1,nlocs
       obs(n)%zdim = 0.0
    END DO

    ! TODO check hour
    ! CALL check( nf90_inq_varid(ncid, "hr", varid),&
    DO n=1,nlocs
       obs(n)%zdim = 0.0
    END DO

    varname2=varname//"@ObsValue"
    CALL check( nf90_inq_varid(ncid, varname2, vid), &
         varname2// ' in '//filename )
    CALL check( nf90_get_var(ncid, vid, tmp_r) )
    DO n=1,nlocs
       obs(n)%val = tmp_r(n)
    END DO

    varname2=varname//"@ObsError"
    CALL check( nf90_inq_varid(ncid, varname2, vid), &
         varname2// ' in '//filename )
    CALL check( nf90_get_var(ncid, vid, tmp_r) )
    DO n=1,nlocs
       obs(n)%err = tmp_r(n)
    END DO

    varname2=varname//"@PreQC"
    ! TODO, use the actual QC from hofx files?
    CALL check( nf90_inq_varid(ncid, varname2, vid), &
         varname2// ' in '//filename )
    CALL check( nf90_get_var(ncid, vid, tmp_i) )
    DO n=1,nlocs
       obs(n)%qc = tmp_i(n)
    END DO

    ! cleanup
    CALL check( nf90_close(ncid) )
    DEALLOCATE(tmp_i, tmp_r)

  END SUBROUTINE read_file
  !================================================================================


  !================================================================================
  !> Helper function to check status of netcdf calls
  !--------------------------------------------------------------------------------
  SUBROUTINE check(status, str)
    INTEGER, INTENT(in) :: status
    CHARACTER(*), OPTIONAL, INTENT(in) :: str

    IF(status /= nf90_noerr) THEN
       IF (PRESENT(str)) THEN
          WRITE (*,*) TRIM(nf90_strerror(status)), ": ",str
       ELSE
          WRITE (*,*) TRIM(nf90_strerror(status))
       END IF
       CALL letkf_mpi_abort("NetCDF error")
    END IF
  END SUBROUTINE check
  !================================================================================


  !================================================================================
  !--------------------------------------------------------------------------------
  SUBROUTINE str_pattern_ens (str_in, str_out, ensmem)
    CHARACTER(len=*), INTENT(in) :: str_in
    CHARACTER(len=:), ALLOCATABLE, INTENT(inout) :: str_out
    INTEGER,  INTENT(in) :: ensmem

    INTEGER :: n, i
    CHARACTER(len=6)  :: pattern
    CHARACTER(len=10) :: fmt
    CHARACTER(len=1024) :: str_tmp

    DO n=1,9
       ! we are looking for "#ENSx#" strings, where "x" is the number of digits.
       ! For example, "#ENS4#" in the filename will be replaced with "0001" for
       ! ensemble member number 1.
       WRITE (pattern, "(A,I0,A)") "#ENS",n, "#"
       i = INDEX(str_in, pattern)

       IF(i > 0) THEN
          WRITE (fmt, '(A,I0,A)') '(A,I0.', n, ',A)'
          WRITE (str_tmp, fmt) str_in(1:i-1), ensmem, &
               str_in(i+LEN(pattern):LEN(str_in))
          EXIT
       END IF
    END DO
    str_out = TRIM(str_tmp)
  END SUBROUTINE str_pattern_ens
  !================================================================================

END MODULE letkf_obsio_ioda_mod
