MODULE letkf_state_nc
  USE netcdf
  USE letkf_state
  USE letkf_mpi

  IMPLICIT NONE
  PRIVATE


  ! Public types
  !------------------------------------------------------------

  !> state file I/O class for handling NetCDF files
  PUBLIC :: stateio_nc
  TYPE, EXTENDS(letkf_stateio) :: stateio_nc
     INTEGER :: compression !< nc4 compression (0-9)
   CONTAINS
     PROCEDURE, NOPASS :: name => stateio_nc_get_name
     PROCEDURE, NOPASS :: desc => stateio_nc_get_desc
     PROCEDURE         :: init => stateio_nc_init
     PROCEDURE         :: read_specs  => stateio_nc_read_specs
     PROCEDURE         :: read_state  => stateio_nc_read_state
     PROCEDURE         :: write_init  => stateio_nc_write_init
     PROCEDURE         :: write_state => stateio_nc_write_state
  END TYPE stateio_nc



  TYPE data_entry
     CHARACTER(len=20)   :: domain
     CHARACTER(len=20)   :: name
     INTEGER             :: dims
     CHARACTER(len=20)   :: zcoord
     CHARACTER(len=20)   :: input_var
     CHARACTER(len=1024) :: input_file
  END TYPE data_entry



  INTEGER, PARAMETER :: fields_max = 100
  INTEGER            :: fields_num = 0
  TYPE(data_entry),TARGET   :: data_entries(fields_max)


CONTAINS



  !--------------------------------------------------------------------------------

  !> Get the unique name of the state I/O class
  FUNCTION stateio_nc_get_name() RESULT(name)
    CHARACTER(:), ALLOCATABLE :: name
    name = "STATEIO_NC"
  END FUNCTION stateio_nc_get_name



  !--------------------------------------------------------------------------------

  !> Get the description of the state I/O class
  FUNCTION stateio_nc_get_desc() RESULT(desc)
    CHARACTER(:), ALLOCATABLE :: desc
    desc = "NetCDF formatted model state I/O"
  END FUNCTION stateio_nc_get_desc



  !--------------------------------------------------------------------------------

  !>
  SUBROUTINE stateio_nc_init(self, nml_filename)
    CLASS(stateio_nc) :: self
    CHARACTER(:), ALLOCATABLE, INTENT(in) :: nml_filename

    CHARACTER(:), ALLOCATABLE :: statedef_file
    CHARACTER(len=1024) :: line, strs(6)
    INTEGER :: unit, iostat, i
    LOGICAL :: ex
    TYPE(data_entry) :: newEntry
    INTEGER :: compression = 0

    NAMELIST /stateio_nc/ statedef_file, compression

    IF (pe_isroot) THEN
       PRINT '(/A)', ""
       PRINT *, " letkf_stateio_nc_init() : "
       PRINT *, "------------------------------------------------------------"
    END IF

    ! load in our section of the namelist
    ALLOCATE(CHARACTER(1024) :: statedef_file);   WRITE(statedef_file, *)  "UNDEFINED"
    OPEN(newunit=unit, file=nml_filename, status='OLD')
    READ(unit, nml=stateio_nc)
    CLOSE(unit)
    statedef_file=TRIM(statedef_file)
    IF (pe_isroot) PRINT stateio_nc
    self%compression = compression

    ! make sure statedef file exists
    INQUIRE(file=statedef_file, exist=ex)
    IF (.NOT. ex) THEN
       CALL letkf_mpi_abort("File does not exist: " // TRIM(statedef_file))
    END IF

    ! read in the statedef file
    OPEN(newunit=unit, file=statedef_file, action='read')
    DO WHILE(.TRUE.)
       ! read a new line
       READ(unit, "(A)", iostat=iostat) line
       IF (iostat<0) EXIT
       IF (iostat>0) CALL letkf_mpi_abort("Error reading file.")

       ! convert tabs to spaces
       DO i=1, LEN(line)
          IF(line(i:i) == CHAR(9)) line(i:i) = ' '
       END DO

       ! ignore comments / empty lines
       line = ADJUSTL(line)
       IF(line(1:1) == '#') CYCLE
       IF(LEN(TRIM(line)) == 0) CYCLE

       ! parse into substrings
       READ(line,*) strs
       newEntry%domain = TRIM(strs(1))
       newEntry%name   = TRIM(strs(2))
       READ(strs(3), *) newEntry%dims
       newEntry%zcoord = TRIM(strs(4))
       newEntry%input_var  = TRIM(strs(5))
       newEntry%input_file = TRIM(strs(6))

       ! add latest entry
       fields_num = fields_num + 1
       data_entries(fields_num) = newEntry
    END DO
  END SUBROUTINE stateio_nc_init



  !--------------------------------------------------------------------------------
  !>
  SUBROUTINE stateio_nc_read_specs(self, hzgrids, vtgrids, statevars)
    CLASS(stateio_nc) :: self
    TYPE(letkf_hzgrid_spec),   ALLOCATABLE, INTENT(out) :: hzgrids(:)
    TYPE(letkf_vtgrid_spec),   ALLOCATABLE, INTENT(out) :: vtgrids(:)
    TYPE(letkf_statevar_spec), ALLOCATABLE, INTENT(out) :: statevars(:)

    REAL, ALLOCATABLE :: tmp_r_1d(:)
    REAL, ALLOCATABLE :: tmp_r_2d(:,:)

    INTEGER :: nx, ny, i, s

    nx = -1
    ny = -1


    ! Read horizontal grid
    ! ------------------------------------------------------------
    !TODO, only handling single hzgrid spec, modify to handle
    ! an arbitrary number
    ALLOCATE(hzgrids(1))
    hzgrids(1)%name = "hzgrd1"

    PRINT *, ""
    PRINT *, "Loading horizontal grid specs..."
    ! read latitud1e
    CALL read_nc_d2("grid_hz", "lat", hzgrids(1)%lat)
    CALL read_nc_d1("grid_hz", "lat", hzgrids(1)%lat_nom)
    IF (ALLOCATED(hzgrids(1)%lat)) THEN
       ! 2d fields was read in
       nx = SIZE(hzgrids(1)%lat, 1)
       ny = SIZE(hzgrids(1)%lat, 2)
    END IF
    IF (ny < 0 .AND. ALLOCATED(hzgrids(1)%lat_nom))  ny = SIZE(hzgrids(1)%lat_nom)

    ! read longitude
    CALL read_nc_d2("grid_hz", "lon", hzgrids(1)%lon)
    CALL read_nc_d1("grid_hz", "lon", hzgrids(1)%lon_nom)
    IF (nx < 0 .AND. ALLOCATED(hzgrids(1)%lon_nom))  nx = SIZE(hzgrids(1)%lon_nom)

    ! if we need to generate the 2d lat/lon from the 1d lat/lon
    IF( .NOT. ALLOCATED(hzgrids(1)%lat) .OR. .NOT. ALLOCATED(hzgrids(1)%lon)) THEN
       IF (.NOT. ALLOCATED(hzgrids(1)%lat_nom) .OR. .NOT. ALLOCATED(hzgrids(1)%lon_nom) ) &
            CALL letkf_mpi_abort("neither a 2D lat/lon nor 1D lat/lon was specified")
       ! TODO test to make sure this is actually correct
       ALLOCATE(hzgrids(1)%lat(nx,ny))
       ALLOCATE(hzgrids(1)%lon(nx,ny))
       DO i=1,nx
          hzgrids(1)%lon(i,:) = hzgrids(1)%lon_nom(i)
       END DO
       DO i=1,ny
          hzgrids(1)%lat(i,:) = hzgrids(1)%lat_nom(i)
       END DO
    END IF

    ! if we need to generate the 1d nominal lat/lon from the 2d lat/lon
    IF (.NOT. ALLOCATED(hzgrids(1)%lat_nom) .OR. .NOT. ALLOCATED(hzgrids(1)%lon_nom)) THEN
       PRINT *, ""
       PRINT *, "WARNING: a 2D lat/lon grid was specified, but a nominal ",&
            "1D lat/lon was not. Estimating a viable nominal 1D lat/lon, but ",&
            "cordinates of resulting files might not be ideal."
       PRINT *, ""
       hzgrids(1)%lat_nom = SUM(hzgrids(1)%lat, dim=1)/nx
       hzgrids(1)%lon_nom = SUM(hzgrids(1)%lon, dim=2)/ny
    END IF

    ! read mask
    CALL read_nc_d2("grid_hz", "mask", tmp_r_2d)
    ALLOCATE(hzgrids(1)%mask(nx,ny))
    hzgrids(1)%mask = .FALSE.
    IF (ALLOCATED(tmp_r_2d))   WHERE (tmp_r_2d <= 0) hzgrids(1)%mask = .TRUE.


    ! read vertical grid
    !------------------------------------------------------------
    ! TODO allow for arbitrary number of vertical grids
    ! TODO, interpret depth/height vs thickness
    ALLOCATE(vtgrids(1))
    PRINT *, ""
    PRINT *, "Loading vertical grids specs..."
    DO i=1, fields_num
       IF (data_entries(i)%domain == 'grid_vt') THEN
          vtgrids(1)%name = data_entries(i)%name
          
          IF (data_entries(i)%dims == 1) THEN
             ! 1 dimensional case
             vtgrids(1)%dims=1
             IF (ALLOCATED(tmp_r_1d)) DEALLOCATE(tmp_r_1d)
             CALL read_nc_d1("grid_vt", data_entries(i)%name, tmp_r_1d)
             ALLOCATE(vtgrids(1)%vert(SIZE(tmp_r_1d),1,1))
             ALLOCATE(vtgrids(1)%vert_nom(size(tmp_r_1d)))
             vtgrids(1)%vert(:,1,1) = tmp_r_1d
             vtgrids(1)%vert_nom = tmp_r_1d
             
          ELSE IF(data_entries(i)%dims == 2) THEN
             ! 2 dimensional case
             !TODO finish this
             CALL letkf_mpi_abort("not yet implemented")
             IF (ALLOCATED(tmp_r_2d)) DEALLOCATE(tmp_r_2d)
             CALL read_nc_d2("grid_vt", data_entries(i)%name, tmp_r_2d)
             
          ELSE
             ! TODO handle 3D case
             CALL letkf_mpi_abort("variable '"//data_entries(i)%name//&
                  "' defined with illegal number of dimensions, must be 1 or 2")
          END IF
       END IF
    END DO


    ! determine the state variables that will be loaded
    !TODO  save the configuration
    PRINT *, ""
    PRINT *, "Loading model state specs..."
    s=0
    DO i=1, fields_num
       IF (data_entries(i)%domain == 'state_in') THEN
          IF (data_entries(i)%dims /= 3) CALL letkf_mpi_abort ("NYI, state with dims /= 3")
          s = s + 1
       END IF
    END DO
    ALLOCATE(statevars(s))
    s=0
    DO i=1, fields_num
       IF (data_entries(i)%domain == 'state_in') THEN
          IF (data_entries(i)%dims /= 3) CALL letkf_mpi_abort ("NYI, state with dims /= 3")
          s = s + 1
          statevars(s)%name=data_entries(i)%name
          statevars(s)%hzgrid = 'hzgrd1' ! TODO, fill in with correct hz grid names
          statevars(s)%vtgrid = data_entries(i)%zcoord
       END IF
    END DO
  END SUBROUTINE stateio_nc_read_specs



  !-----------------------------------------------------------------------------
  !>
  SUBROUTINE stateio_nc_write_init(self, ftype, ensmem)
    CLASS(stateio_nc) :: self
    CHARACTER(len=*), INTENT(in) :: ftype
    INTEGER, INTENT(in) :: ensmem

    CHARACTER(:), ALLOCATABLE :: filename, arg1, arg2
    INTEGER :: ncid, d_t, d_x, d_y, d_z(100), varid ! TODO, remove hardcoded dz size
    INTEGER :: i, j

    ! TODO, use a template from the namelist
    ! determine the filename to write out to

    ! TODO, i'm getting weird trailing whitespace
    filename="#TYPE#.#ENS4#.nc"
    arg1="#TYPE#"
    arg2=TRIM(ftype)
    filename=str_pattern(filename, arg1, arg2)
    filename=str_ens_pattern(filename, ensmem)

    if (self%verbose) &
         PRINT '(X,A,I0.3,2A)', " PE ", pe_rank," initializing output for ", filename

    ! start creating output file
    CALL check(nf90_create(filename, nf90_hdf5, ncid))

    CALL check(nf90_def_dim(ncid, "time", 0, d_t))
    call check(nf90_def_var(ncid, "time", nf90_real, (/d_t/), varid))
    call check(nf90_put_att(ncid, varid, "axis", "T"))
    ! TODO: fill in time data    
    !call chekc(nf90_put_var(

    ! TODO: add support for multiple horizontal grids
    CALL check(nf90_def_dim(ncid, "lon", grid_nx, d_x))
    call check(nf90_def_var(ncid, "lon", nf90_real, (/d_x/), varid))
    call check(nf90_put_att(ncid, varid, "long_name", "longitude"))
    call check(nf90_put_att(ncid, varid, "units", "degrees_east"))
    call check(nf90_put_att(ncid, varid, "axis", "X"))
    
    CALL check(nf90_def_dim(ncid, "lat", grid_ny, d_y))
    call check(nf90_def_var(ncid, "lat", nf90_real, (/d_y/), varid))
    call check(nf90_put_att(ncid, varid, "long_name", "latitude"))    
    call check(nf90_put_att(ncid, varid, "units", "degrees_north"))
    call check(nf90_put_att(ncid, varid, "axis", "Y"))

    do i=1,size(vtgrids)
       !TODO, handle 2D, 3D, and constant vertical grids
       if (vtgrids(i)%dims > 1)&
            call letkf_mpi_abort ("LETKF_STATE_NC: 2D/3D vertical grids not yet supported.")
       CALL check(nf90_def_dim(ncid, vtgrids(i)%name, size(vtgrids(i)%vert, 1), d_z(i)))
       call check(nf90_def_var(ncid, vtgrids(i)%name, nf90_real, (/d_z(i)/), varid))
       call check(nf90_put_att(ncid, varid, "axis", "Z"))
    end do

    ! for each "state_out" defined in the config file, create the output variable
    !TODO, handle file name output templating where variables are split out to multiple files
    !TODO, remove the "state_out" section, and just add another set of variable name, file, fields?
    DO i = 1, SIZE(data_entries)
       IF (data_entries(i)%domain=="state_out") THEN
          ! determine which vertical coordinate is the matching one
          do j=1, size(vtgrids)
             if (vtgrids(j)%name == data_entries(i)%zcoord) exit
          end do
          if(j > size(vtgrids))&
               call letkf_mpi_abort("cannot find vertical grid in stateio_nc_write_init")
          
          CALL check(nf90_def_var(ncid, data_entries(i)%input_var,&
               nf90_real, (/d_x,d_y,d_z(j),d_t/),varid))
          CALL check(nf90_def_var_deflate(ncid, varid, 1, 1, self%compression))
       END IF
    END DO

    ! write out nominal lat/lon/vert
    CALL check(nf90_enddef(ncid))
    call check(nf90_inq_varid(ncid, "lon", varid))    
    call check(nf90_put_var(ncid, varid, hzgrids(1)%lon_nom))
    call check(nf90_inq_varid(ncid, "lat", varid))    
    call check(nf90_put_var(ncid, varid, hzgrids(1)%lat_nom))
    do i =1, size(vtgrids)
       call check(nf90_inq_varid(ncid, vtgrids(i)%name, varid))    
       call check(nf90_put_var(ncid, varid, vtgrids(i)%vert_nom))
    end do

    ! all done
    CALL check(nf90_close(ncid))
  END SUBROUTINE stateio_nc_write_init



  !-----------------------------------------------------------------------------
  !>
  SUBROUTINE stateio_nc_write_state(self, ftype, ensmem, state_var, state_val)
    CLASS(stateio_nc) :: self
    CHARACTER(len=*), INTENT(in) :: ftype
    INTEGER,        INTENT(in) :: ensmem
    CHARACTER(*),   INTENT(in) :: state_var
    REAL,           INTENT(in) :: state_val(:,:,:)

    CHARACTER(:), ALLOCATABLE :: filename, varname, arg, arg2

    INTEGER :: i
    INTEGER :: ncid, varid

    ! TODO, use a template from the namelist to determine the filename to write out to
    ! TODO, i'm getting weird trailing whitespace
    filename="#TYPE#.#ENS4#.nc"
    arg="#TYPE#"
    arg2=TRIM(ftype)
    filename=str_pattern(filename,arg,arg2)
    filename=str_ens_pattern(filename, ensmem)

    !PRINT *, "Writing out to file "    // TRIM(filename)//"..."
    CALL check(nf90_open(filename, nf90_write, ncid))

    ! find which variable this data should be written to
    varname = ""
    DO i=1,SIZE(data_entries)
       IF (data_entries(i)%domain =="state_out" .AND. &
            data_entries(i)%name == state_var) THEN
          varname = data_entries(i)%input_var
          EXIT
       END IF
    END DO

    ! if line in the config file was not found, abort
    IF (varname == "") CALL letkf_mpi_abort("Variable "//TRIM(state_var)//&
         " was not " // "found in the configuration file")

    ! write out the variable to the file
    CALL check(nf90_inq_varid(ncid, varname, varid))
    CALL check(nf90_put_var(ncid, varid, state_val))
    CALL check(nf90_sync(ncid))

    ! all done
    CALL check(nf90_close(ncid))

  END SUBROUTINE stateio_nc_write_state



  !--------------------------------------------------------------------------------
  !>
  !! TODO, modify so that we can read/scatter the state 1 variable and/or 1 level
  !! at a time in order to reduce the peak memory per compute node
  SUBROUTINE stateio_nc_read_state(self, ensmem, state_var, state_val)
    CLASS(stateio_nc)    :: self
    INTEGER,      INTENT(in)  :: ensmem
    CHARACTER(*), INTENT(in)  :: state_var
    REAL, ALLOCATABLE, INTENT(out) :: state_val(:,:,:)

    INTEGER :: i
    TYPE(data_entry), POINTER :: def
    CHARACTER(:), ALLOCATABLE :: filename

    INTEGER :: nx, ny, nz
    INTEGER :: ncid, dimids(4), varid

    ! find matching data entry from config file;
    ! TODO, use the toupper function
    NULLIFY(def)
    DO i = 1, fields_num
       IF ( (data_entries(i)%domain == "state_in") .AND. &
            (data_entries(i)%name == state_var)) THEN
          def => data_entries(i)
          EXIT
       END IF
    END DO
    IF (.NOT. ASSOCIATED(def)) CALL letkf_mpi_abort("state variable "//&
         state_var//" not found in configuration file.")

    ! determine the filename to load in. Subsituting #ENSx# with the ensemble num
    filename=def%input_file
    filename=str_ens_pattern(filename, ensmem)

    ! make sure this file exists
    ! TODO
    if (self%verbose) &
         PRINT '(X, A,I0.3,A,I0.3,6A)',' PE ',pe_rank,' loading ens ', ensmem, ' "',&
         TRIM(state_var),'" from "',TRIM(def%input_var),'" of "', TRIM(filename)

    ! load field
    ! TODO, handle 2D vars
    ! TODO, handle empty time variables
    CALL check(nf90_open(filename, nf90_nowrite, ncid))
    CALL check(nf90_inq_varid(ncid, def%input_var, varid))
    CALL check(nf90_inquire_variable(ncid, varid, ndims=i))
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids))
    CALL check(nf90_inquire_dimension(ncid, dimids(1), len=nx))
    CALL check(nf90_inquire_dimension(ncid, dimids(2), len=ny))
    CALL check(nf90_inquire_dimension(ncid, dimids(3), len=nz))
    ALLOCATE(state_val(nx,ny,nz))
    CALL check(nf90_get_var(ncid, varid, state_val))
    CALL check(nf90_close(ncid))
  END SUBROUTINE stateio_nc_read_state



  !-----------------------------------------------------------------------------
  !>TODO why am I tracking 'bkg' and 'ana' separately??
  FUNCTION str_ens_pattern(str_in, ensmem) RESULT(str_out)
    CHARACTER(:), ALLOCATABLE, INTENT(in) :: str_in
    INTEGER, INTENT(in) :: ensmem
    CHARACTER(:), ALLOCATABLE :: str_out

    INTEGER :: i, n
    CHARACTER(len=6)  :: pattern
    CHARACTER(len=10) :: fmt

    !TODO handle mean and spread
    str_out = str_in
    DO n=1,9
       WRITE (pattern, "(A,I0,A)") "#ENS",n,"#"
       i = INDEX(str_out, pattern)
       IF(i > 0) THEN
          IF (ensmem >= 0) THEN
             WRITE (fmt, '(A,I0,A)') '(A,I0.',n,',A)'
             WRITE (str_out, fmt) str_out(1:i-1), ensmem, &
                  str_out(i+LEN(pattern):LEN(str_out))
          ELSE IF ( ensmem == ENS_BKG_MEAN .OR. ensmem == ENS_ANA_MEAN) THEN
             WRITE (str_out, '(3A)') str_out(1:i-1), "mean", &
                  str_out(i+LEN(pattern):LEN(str_out))
          ELSE IF ( ensmem == ENS_BKG_SPRD .OR. ensmem == ENS_ANA_SPRD) THEN
             WRITE (str_out, '(3A)') str_out(1:i-1), "sprd", &
                  str_out(i+LEN(pattern):LEN(str_out))
          ELSE
             CALL letkf_mpi_abort("Illegal 'ensmem' given to str_ens_pattern")
          END IF
          EXIT
       END IF
    END DO
  END FUNCTION str_ens_pattern



  !-----------------------------------------------------------------------------
  FUNCTION str_pattern(str_in, key, val) RESULT(str_out)
    CHARACTER(:), ALLOCATABLE, INTENT(in) :: str_in
    CHARACTER(:), ALLOCATABLE, INTENT(in) :: key
    CHARACTER(:), ALLOCATABLE, INTENT(in) :: val
    CHARACTER(:), ALLOCATABLE :: str_out

    INTEGER :: i

    ALLOCATE(CHARACTER(len=LEN(str_in)+LEN(val)) :: str_out)
    i = INDEX(str_in, key)

    IF(i>0) THEN
       str_out = str_in(1:i-1) // TRIM(val) // str_in(i+LEN(key):LEN(str_in))
    ELSE
       str_out = str_in
    END IF

  END FUNCTION str_pattern




  !--------------------------------------------------------------------------------
  !>
  SUBROUTINE read_nc_d2(domain, name, array)
    CHARACTER(*), INTENT(in)  :: domain
    CHARACTER(*), INTENT(in)  :: name
    REAL, ALLOCATABLE, INTENT(out) :: array(:,:)

    INTEGER :: i, dimids(2), nx, ny
    TYPE(data_entry), POINTER :: def
    INTEGER :: ncid, varid

    ! find matching data entry from config file;
    ! TODO, use the toupper function
    NULLIFY(def)
    DO i = 1, fields_num
       IF ( (data_entries(i)%domain == domain) .AND. &
            (data_entries(i)%name == name) .AND. &
            (data_entries(i)%dims == 2)) THEN
          def => data_entries(i)
          EXIT
       END IF
    END DO

    ! if we didn't find a definition:
    IF (.NOT. ASSOCIATED(def)) RETURN

    PRINT *,' loading 2D "', TRIM(name),'" from "', TRIM(def%input_var),&
         '" of "',TRIM(def%input_file)

    ! open file, find dimensions
    ! TODO, handle a time dimension if stored that way
    CALL check(nf90_open(def%input_file, nf90_nowrite, ncid))
    CALL check(nf90_inq_varid(ncid, def%input_var, varid))
    CALL check(nf90_inquire_variable(ncid, varid, ndims=i))
    IF ( i /= 2) CALL letkf_mpi_abort("variable dimensions /= 2")
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids))
    CALL check(nf90_inquire_dimension(ncid, dimids(1), len=nx))
    CALL check(nf90_inquire_dimension(ncid, dimids(2), len=ny))

    ! read in the variable
    ALLOCATE(array(nx,ny))
    CALL check(nf90_get_var(ncid, varid, array))

    ! cleanup
    CALL check(nf90_close(ncid))

  END SUBROUTINE read_nc_d2



  !--------------------------------------------------------------------------------
  !>
  SUBROUTINE read_nc_d1(domain, name, array)
    CHARACTER(*), INTENT(in)  :: domain
    CHARACTER(*), INTENT(in)  :: name
    REAL, ALLOCATABLE, INTENT(out) :: array(:)

    INTEGER :: i, dimids(1), ni
    TYPE(data_entry), POINTER :: def
    INTEGER :: ncid, varid

    NULLIFY(def)
    ! find matching data entry from config file;
    ! TODO, use the toupper function
    DO i = 1, fields_num
       IF ( (data_entries(i)%domain == domain) .AND. &
            (data_entries(i)%name == name) .AND. &
            (data_entries(i)%dims == 1)) THEN
          def => data_entries(i)
          EXIT
       END IF
    END DO

    ! if we didn't find a definition:
    IF (.NOT. ASSOCIATED(def)) RETURN

    PRINT *,' loading 1D "', TRIM(name), '" from "', TRIM(def%input_var), &
         '" of "', TRIM(def%input_file)

    ! open file, find dimensions
    ! TODO, handle a time dimension if stored that way
    CALL check(nf90_open(def%input_file, nf90_nowrite, ncid))
    CALL check(nf90_inq_varid(ncid, def%input_var, varid))
    CALL check(nf90_inquire_variable(ncid, varid, ndims=i))
    IF ( i /= 1) CALL letkf_mpi_abort("variable dimensions /= 1")
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids))
    CALL check(nf90_inquire_dimension(ncid, dimids(1), len=ni))

    ! read in the variable
    ALLOCATE(array(ni))
    CALL check(nf90_get_var(ncid, varid, array))

    ! cleanup
    CALL check(nf90_close(ncid))

  END SUBROUTINE read_nc_d1



  !-----------------------------------------------------------------------------
  !> Helper function to check status of netcdf calls
  SUBROUTINE check(status, str)
    INTEGER, INTENT(in) :: status
    CHARACTER(*), OPTIONAL, INTENT(in) :: str

    IF(status /= nf90_noerr) THEN
       IF(PRESENT(str)) THEN
          WRITE (*,*) TRIM(nf90_strerror(status)), ": ", str
       ELSE
          WRITE (*,*) TRIM(nf90_strerror(status))
       END IF
       CALL letkf_mpi_abort("NetCDF error")
    END IF
  END SUBROUTINE check

END MODULE letkf_state_nc
