MODULE letkf_state_nc
  USE netcdf
  use letkf_config
  USE letkf_state
  USE letkf_mpi

  IMPLICIT NONE
  PRIVATE

  !TODO the configuration files are a bit of a mess
  ! TODO, replace the fixed sized arrayes with linked lists


  ! Public types
  !------------------------------------------------------------

  INTEGER, PARAMETER :: fields_max = 100


  !> state file I/O class for handling NetCDF files
  PUBLIC :: stateio_nc  
  TYPE, EXTENDS(letkf_stateio) :: stateio_nc
     INTEGER :: compression !< nc4 compression (0-9)
     type(configuration) :: hzgriddef
     type(configuration) :: vtgriddef
     type(configuration) :: statedef
     
   CONTAINS

     PROCEDURE, NOPASS :: name => stateio_nc_get_name
     PROCEDURE, NOPASS :: desc => stateio_nc_get_desc
     PROCEDURE         :: init => stateio_nc_init
     PROCEDURE         :: read_specs  => stateio_nc_read_specs
     PROCEDURE         :: read_state  => stateio_nc_read_state
     PROCEDURE         :: write_init  => stateio_nc_write_init
     PROCEDURE         :: write_state => stateio_nc_write_state
  END TYPE stateio_nc



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
  SUBROUTINE stateio_nc_init(self, config)
    CLASS(stateio_nc) :: self
    type(configuration), intent(in) :: config

    CHARACTER(len=1024) :: line, strs(6)
    INTEGER :: unit, iostat, i
    LOGICAL :: ex
    INTEGER :: compression = 0


    IF (pe_isroot) THEN
       PRINT '(/A)', ""
       PRINT *, " letkf_stateio_nc_init() : "
       PRINT *, "------------------------------------------------------------"
    END IF

    ! load in our section of the configuration
    call config%get("compression", self%compression)
    if (pe_isroot) print *, " state.compression=", self%compression
    call config%get("hzgrid", self%hzgriddef)
    call config%get("vtgrid", self%vtgriddef)
    call config%get("statedef", self%statedef)

  END SUBROUTINE stateio_nc_init



  !--------------------------------------------------------------------------------
  !>
  SUBROUTINE stateio_nc_read_specs(self, hzgrids, vtgrids, statevars)
    CLASS(stateio_nc) :: self
    TYPE(letkf_hzgrid_spec),   ALLOCATABLE, INTENT(out) :: hzgrids(:)
    TYPE(letkf_vtgrid_spec),   ALLOCATABLE, INTENT(out) :: vtgrids(:)
    TYPE(letkf_statevar_spec), ALLOCATABLE, INTENT(out) :: statevars(:)

    type(configuration) :: config, cfg2
    
    character(len=:), allocatable :: str, filename, varname
    REAL, ALLOCATABLE :: tmp_r_1d(:)
    REAL, ALLOCATABLE :: tmp_r_2d(:,:)

    INTEGER :: nx, ny, i, s, c, grd_cnt, j
    character(len=100) :: strs(100) ! TODO, make a linked list
    nx = -1
    ny = -1

    ! Read horizontal grids
    ! ------------------------------------------------------------
    PRINT *, ""
    PRINT *, "Loading horizontal grid specs..."
    grd_cnt = self%hzgriddef%count()
    allocate(hzgrids(grd_cnt))
    print *, "Found ",grd_cnt," horizontal grids definitions."
    print *, ""
    
    !  for each horizontal grid:
    do i=1,grd_cnt
       
       ! grid name
       call self%hzgriddef%get(i, config, name=str)
       print *, "loading horizontal grid: ", str
       hzgrids(i)%name = str
       
       ! latitude
       if(config%found("lat2d")) then
          call config%get("lat2d", cfg2)
          call cfg2%get(1, varname)
          call cfg2%get(2, filename)
          PRINT *,' loading 2D "lat2d" from "', varname, '" of "',filename,'"'          
          call read_nc_d2(varname, filename, hzgrids(i)%lat)
       end if
       if(config%found("lat1d")) then
          call config%get("lat1d", cfg2)
          call cfg2%get(1, varname)
          call cfg2%get(2, filename)
          PRINT *,' loading 1D "lat1d" from "', varname, '" of "',filename,'"'          
          call read_nc_d1(varname, filename, hzgrids(i)%lat_nom)
       end if
       IF (ALLOCATED(hzgrids(i)%lat)) THEN
          ! 2d fields was read in
          nx = size(hzgrids(i)%lat, 1)
          ny = size(hzgrids(i)%lat, 2)
       end if
       IF (ny < 0 .AND. ALLOCATED(hzgrids(i)%lat_nom)) &
            ny = SIZE(hzgrids(i)%lat_nom)


       ! longitude
       if(config%found("lon2d")) then
          call config%get("lon2d", cfg2)
          call cfg2%get(1, varname)
          call cfg2%get(2, filename)
          PRINT *,' loading 2D "lon2d" from "', varname, '" of "',filename,'"'          
          call read_nc_d2(varname, filename, hzgrids(i)%lon)
       end if
       if(config%found("lon1d")) then
          call config%get("lon1d", cfg2)
          call cfg2%get(1, varname)
          call cfg2%get(2, filename)
          PRINT *,' loading 1D "lon1d" from "', varname, '" of "',filename,'"'          
          call read_nc_d1(varname, filename, hzgrids(i)%lon_nom)
       end if
       IF (nx < 0 .AND. ALLOCATED(hzgrids(i)%lon_nom)) &
            nx = SIZE(hzgrids(i)%lon_nom)

       ! if we need to generate the 2d lat/lon from the 1d lat/lon
       IF( .NOT. ALLOCATED(hzgrids(i)%lat) .OR. .NOT. ALLOCATED(hzgrids(i)%lon)) THEN
          IF (.NOT. ALLOCATED(hzgrids(i)%lat_nom) .OR. .NOT. ALLOCATED(hzgrids(i)%lon_nom) ) &
               CALL letkf_mpi_abort("neither a 2D lat/lon nor 1D lat/lon was specified")
          ALLOCATE(hzgrids(i)%lat(nx,ny))
          ALLOCATE(hzgrids(i)%lon(nx,ny))
          DO j=1,nx
             hzgrids(i)%lon(j,:) = hzgrids(i)%lon_nom(j)
          END DO
          DO j=1,ny
             hzgrids(i)%lat(j,:) = hzgrids(i)%lat_nom(j)
          END DO
       END IF

       ! if we need to generate the 1d nominal lat/lon from the 2d lat/lon
       IF (.NOT. ALLOCATED(hzgrids(i)%lat_nom) .OR. .NOT. ALLOCATED(hzgrids(i)%lon_nom)) THEN
          PRINT *, ""
          PRINT *, "WARNING: a 2D lat/lon grid was specified, but a nominal ",&
               "1D lat/lon was not. Estimating a viable nominal 1D lat/lon, but ",&
               "cordinates of resulting files might not be ideal."
          PRINT *, ""
          hzgrids(i)%lat_nom = SUM(hzgrids(i)%lat, dim=1)/nx
          hzgrids(i)%lon_nom = SUM(hzgrids(i)%lon, dim=2)/ny
       END IF

       ! read mask
       allocate(hzgrids(i)%mask(nx,ny))
       hzgrids(i)%mask = .FALSE.       
       if(config%found("mask")) then
          call config%get("mask", cfg2)
          call cfg2%get(1, varname)
          call cfg2%get(2, filename)
          PRINT *,' loading 2D "mask"  from "', varname, '" of "',filename,'"'          
          CALL read_nc_d2(varname, filename,tmp_r_2d)
          WHERE (tmp_r_2d <= 0) hzgrids(i)%mask = .TRUE.          
       end if

       print *, ""
    end do


    ! read vertical grid
    !------------------------------------------------------------
    ! TODO, interpret depth/height vs thickness
    PRINT *, ""
    PRINT *, "Loading vertical grid specs..."
    grd_cnt = self%vtgriddef%count()
    allocate(vtgrids(grd_cnt))
    print *, "Found ",grd_cnt," vertical grids definitions."
    print *, ""
    
    do i=1,grd_cnt
       ! grid name
       call self%vtgriddef%get(i, config, name=str)
       print *, "loading vertical grid: ", str
       vtgrids(i)%name = str
       
       ! levels
       IF (ALLOCATED(tmp_r_1d)) DEALLOCATE(tmp_r_1d)       
       if(config%found("vert1d")) then
          call config%get("vert1d", cfg2)
          call cfg2%get(1, varname)
          call cfg2%get(2, filename)
          PRINT *,' loading 1D "vert1d" from "', varname, '" of "',filename,'"'          
          call read_nc_d1(varname, filename, tmp_r_1d)
          vtgrids(i)%dims=1
          ALLOCATE(vtgrids(i)%vert(SIZE(tmp_r_1d),1,1))
          ALLOCATE(vtgrids(i)%vert_nom(size(tmp_r_1d)))
          vtgrids(i)%vert(:,1,1) = tmp_r_1d
          vtgrids(i)%vert_nom = tmp_r_1d          
       else
          call letkf_mpi_abort("vert1d is missing for "//trim(vtgrids(i)%name))          
       end if
       
       print *, ""
    END DO


    ! Read state variable config
    ! ------------------------------------------------------------
    PRINT *, ""
    PRINT *, "Loading model state specs..."
    grd_cnt = self%statedef%count()
    allocate(statevars(grd_cnt))
    print *, "Found ",grd_cnt," state definitions."
    print *, ""
    do i=1, grd_cnt
       ! grid name
       call self%statedef%get(i, config, name=str)
       statevars(i)%name = str
       call config%get("hzgrid", str)
       statevars(i)%hzgrid = str
       call config%get("vtgrid", str)
       statevars(i)%vtgrid = str
    end do

  END SUBROUTINE stateio_nc_read_specs



  !-----------------------------------------------------------------------------
  !>
  SUBROUTINE stateio_nc_write_init(self, ftype, ensmem)
    CLASS(stateio_nc) :: self
    CHARACTER(len=*), INTENT(in) :: ftype
    INTEGER, INTENT(in) :: ensmem
    type(configuration) :: config, config2
    
    CHARACTER(:), ALLOCATABLE :: filename, arg1, arg2, name, str
    INTEGER :: ncid, d_t, d_x, d_y, d_z(100), varid ! TODO, remove hardcoded dz size
    INTEGER :: i, j

    ! TODO, use a template from the namelist
    ! determine the filename to write out to

    ! Determine the filename to save to
    filename = str_ens_pattern(str_pattern("#TYPE#.#ENS4#.nc", "#TYPE#", TRIM(ftype)),ensmem)
    if (self%verbose) &
         PRINT '(X,A,I0.3,2A)', " PE ", pe_rank,' initializing output for "'//filename//'"'

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
    DO i = 1, self%statedef%count()
       call self%statedef%get(i, config, name)
       
       ! determine which vertical coordinate is the matching one
       call config%get("vtgrid", str)
       do j=1, size(vtgrids)
          if (vtgrids(j)%name == str) exit
       end do
       if(j > size(vtgrids))&
            call letkf_mpi_abort("cannot find vertical grid in stateio_nc_write_init")

       call config%get("output", config2)
       call config2%get(1, str)
       CALL check(nf90_def_var(ncid, str,&
            nf90_real, (/d_x,d_y,d_z(j),d_t/),varid))
       CALL check(nf90_def_var_deflate(ncid, varid, 1, 1, self%compression))
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

    type(configuration) :: config, config2
    CHARACTER(:), ALLOCATABLE :: filename, varname, arg, arg2, name

    INTEGER :: i
    INTEGER :: ncid, varid

    ! TODO, use a template from the namelist to determine the filename to write out to
    filename="#TYPE#.#ENS4#.nc"
    filename=str_pattern(filename,"#TYPE#",trim(ftype))
    filename=trim(str_ens_pattern(filename, ensmem))

    if (self%verbose) &
         PRINT '(X,A,I0.3,A)', " PE ", pe_rank,' writing "'//state_var//'" out to file "'//filename//'"'

    CALL check(nf90_open(filename, nf90_write, ncid))

    ! find which variable this data should be written to
    varname = ""
    DO i=1, self%statedef%count()
       call self%statedef%get(i, config, name)
       IF (name == state_var) THEN
          call config%get("output", config2)
          call config2%get(1, varname)
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

    type(configuration) :: config, config2
    INTEGER :: i, j
    CHARACTER(:), ALLOCATABLE :: filename, invar, str

    INTEGER :: nx, ny, nz
    INTEGER :: ncid, dimids(4), varid
    type(letkf_hzgrid_spec) :: hgrid
    type(letkf_vtgrid_spec) :: vgrid

    ! TODO, use the toupper function
    
    ! find matching data entry from config file
    call self%statedef%get(state_var, config)

    ! determine the filename to load in. Subsituting #ENSx# with the ensemble num
    call config%get("input", config2)
    call config2%get(1, invar)
    call config2%get(2, filename)
    filename=str_ens_pattern(filename, ensmem)

    ! make sure this file exists
    ! TODO
    if (self%verbose) &
         PRINT '(X, A,I0.3,A,I0.3,6A)',' PE ',pe_rank,' loading ens ', ensmem, ' "',&
         TRIM(state_var),'" from "',TRIM(invar),'" of "', TRIM(filename)

    ! get the x,y,z dimensions based on the specified grid
    call config%get("hzgrid", str)
    hgrid = letkf_state_hzgrid_getdef(str)
    call config%get("vtgrid", str)
    vgrid = letkf_state_vtgrid_getdef(str)

    ! load field
    ! TODO, handle 2D vars
    ! TODO, handle empty time variables
    CALL check(nf90_open(filename, nf90_nowrite, ncid))
    CALL check(nf90_inq_varid(ncid, invar, varid))
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
    CHARACTER(*), INTENT(in) :: str_in
    INTEGER, INTENT(in) :: ensmem
    CHARACTER(:), ALLOCATABLE :: str_out

    character(:), allocatable :: str
    INTEGER :: i, n
    CHARACTER(len=6)  :: pattern
    CHARACTER(len=10) :: fmt

    !TODO handle mean and spread
    str = str_in
    DO n=1,9
       WRITE (pattern, "(A,I0,A)") "#ENS",n,"#"
       i = INDEX(str, pattern)
       IF(i > 0) THEN
          IF (ensmem >= 0) THEN
             WRITE (fmt, '(A,I0,A)') '(A,I0.',n,',A)'
             WRITE (str, fmt) str(1:i-1), ensmem, &
                  str(i+LEN(pattern):LEN(str))
          ELSE IF ( ensmem == ENS_BKG_MEAN .OR. ensmem == ENS_ANA_MEAN) THEN
             WRITE (str, '(3A)') str(1:i-1), "mean", &
                  str(i+LEN(pattern):LEN(str))
          ELSE IF ( ensmem == ENS_BKG_SPRD .OR. ensmem == ENS_ANA_SPRD) THEN
             WRITE (str, '(3A)') str(1:i-1), "sprd", &
                  str(i+LEN(pattern):LEN(str))
          ELSE
             CALL letkf_mpi_abort("Illegal 'ensmem' given to str_ens_pattern")
          END IF
          EXIT
       END IF
    END DO
    str_out=trim(str)
  END FUNCTION str_ens_pattern



  !-----------------------------------------------------------------------------
  FUNCTION str_pattern(str_in, key, val) RESULT(str_out)
    CHARACTER(*), INTENT(in) :: str_in
    CHARACTER(*), INTENT(in) :: key
    CHARACTER(*), INTENT(in) :: val
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
  SUBROUTINE read_nc_d2(varname, filename, array)
    CHARACTER(*), INTENT(in)  :: varname, filename
    REAL, ALLOCATABLE, INTENT(out) :: array(:,:)

    INTEGER :: i, dimids(2), nx, ny
    integer :: def
    INTEGER :: ncid, varid


    ! open file, find dimensions
    ! TODO, handle a time dimension if stored that way
    CALL check(nf90_open(filename, nf90_nowrite, ncid))
    CALL check(nf90_inq_varid(ncid, varname, varid))
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
  SUBROUTINE read_nc_d1(varname, filename,  array)
    CHARACTER(*), INTENT(in)  :: varname, filename
    REAL, ALLOCATABLE, INTENT(out) :: array(:)

    INTEGER :: i, dimids(1), ni, def
    INTEGER :: ncid, varid

    ! open file, find dimensions
    ! TODO, handle a time dimension if stored that way
    CALL check(nf90_open(filename, nf90_nowrite, ncid))
    call check(nf90_inq_varid(ncid, varname, varid))
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
