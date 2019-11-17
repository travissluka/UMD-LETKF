! Copyright 2016-2019 Travis Sluka
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
!> NetCDF based state I/O
!================================================================================
MODULE letkf_stateio_nc_mod
  USE netcdf
  USE letkf_config
  USE letkf_state
  USE letkf_state_helper
  USE letkf_mpi

  IMPLICIT NONE
  PRIVATE



  !================================================================================
  !================================================================================
  ! public module components
  !================================================================================
  !================================================================================



  !================================================================================
  !> state file I/O class for handling NetCDF files
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC, EXTENDS(letkf_stateio) :: letkf_stateio_nc
     INTEGER :: compression !< nc4 compression (0-9)

     TYPE(statevar_files), ALLOCATABLE :: statevars_files(:)

   CONTAINS

     PROCEDURE, NOPASS :: name => stateio_nc_get_name
     PROCEDURE, NOPASS :: desc => stateio_nc_get_desc
     PROCEDURE         :: init => stateio_nc_init
     PROCEDURE         :: read_state  => stateio_nc_read_state
     PROCEDURE         :: write_state => stateio_nc_write_state
  END TYPE letkf_stateio_nc
  !================================================================================




CONTAINS



  !================================================================================
  !> Get the unique name of the state I/O class
  !--------------------------------------------------------------------------------
  FUNCTION stateio_nc_get_name() RESULT(name)
    CHARACTER(:), ALLOCATABLE :: name
    name = "STATEIO_NC"
  END FUNCTION stateio_nc_get_name
  !================================================================================




  !================================================================================
  !> Get the description of the state I/O class
  !--------------------------------------------------------------------------------
  FUNCTION stateio_nc_get_desc() RESULT(desc)
    CHARACTER(:), ALLOCATABLE :: desc
    desc = "NetCDF formatted model state I/O"
  END FUNCTION stateio_nc_get_desc
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE stateio_nc_init(self, config, hzgrids, vtgrids, statevars)
    CLASS(letkf_stateio_nc) :: self
    TYPE(configuration), INTENT(in) :: config
    TYPE(letkf_hzgrid_spec),   ALLOCATABLE, INTENT(out) :: hzgrids(:)
    TYPE(letkf_vtgrid_spec),   ALLOCATABLE, INTENT(out) :: vtgrids(:)
    TYPE(letkf_statevar_spec), ALLOCATABLE, INTENT(out) :: statevars(:)

    TYPE(hzgrid_files),        ALLOCATABLE :: hzgrids_files(:)
    TYPE(vtgrid_files),        ALLOCATABLE :: vtgrids_files(:)

    REAL, ALLOCATABLE :: tmp_r_2d(:,:)
    INTEGER :: nx, ny, i


    IF (pe_isroot) THEN
       PRINT '(/A)', ""
       PRINT *, " letkf_stateio_nc_init() : "
       PRINT *, "------------------------------------------------------------"
    END IF


    ! load in our section of the configuration
    CALL config%get("compression", self%compression, 0)
    IF (pe_isroot) PRINT *, " state.compression=", self%compression


    ! Read horizontal configuration
    CALL parse_hzgrids(config, hzgrids, hzgrids_files)

    ! onlyt the root PE needs to read in the actual horizontal grid
    ! letkf_state will handle scattering across the PEs
    IF (pe_isroot) THEN
       DO i = 1, SIZE(hzgrids)
          ! load lats
          IF (hzgrids_files(i)%lat%var /= "") CALL read_nc_d2( &
               hzgrids_files(i)%lat%var, hzgrids_files(i)%lat%file, &
               hzgrids(i)%lat)
          IF (hzgrids_files(i)%nomlat%var /= "") CALL read_nc_d1( &
               hzgrids_files(i)%nomlat%var, hzgrids_files(i)%nomlat%file, &
               hzgrids(i)%lat_nom)

          ! load lons
          IF (hzgrids_files(i)%lon%var /= "") CALL read_nc_d2( &
               hzgrids_files(i)%lon%var, hzgrids_files(i)%lon%file, &
               hzgrids(i)%lon)
          IF (hzgrids_files(i)%nomlon%var /= "") CALL read_nc_d1( &
               hzgrids_files(i)%nomlon%var, hzgrids_files(i)%nomlon%file, &
               hzgrids(i)%lon_nom)

          ! generate nominal 1D lat/lon, or fill in 2D lat/lon from 1D lat/lon
          CALL check_hzgrid(hzgrids(i))
          nx = SIZE(hzgrids(i)%lon_nom)
          ny = SIZE(hzgrids(i)%lat_nom)

          ! read in mask
          ALLOCATE(hzgrids(i)%mask(nx,ny))
          hzgrids(i)%mask = .FALSE.
          IF(hzgrids_files(i)%mask%var /= "") THEN
             CALL read_nc_d2( &
                  hzgrids_files(i)%mask%var, hzgrids_files(i)%mask%file, &
                  tmp_r_2d)
             WHERE (tmp_r_2d <= 0) hzgrids(i)%mask =.TRUE.
             DEALLOCATE(tmp_r_2d)
          END IF
       END DO
    END IF


    ! read vertical grid configuration
    CALL parse_vtgrids(config, vtgrids, vtgrids_files)

    ! only the root PE needs to read in the actual vertical grid
    ! letkf_state will handle broadcasting to the other PEs
    IF (pe_isroot) THEN
       DO i = 1, SIZE(vtgrids)
          IF (vtgrids_files(i)%vt1d%file == "#CONST#") THEN
             ! special case: if a 1 level constant vertical coord is desired
             ALLOCATE(vtgrids(i)%vert_nom(1))
             READ(vtgrids_files(i)%vt1d%var, *) vtgrids(i)%vert_nom(1)
          ELSE
             ! read in a 1 dimensional vertical coordinate
             CALL read_nc_d1( &
               vtgrids_files(i)%vt1d%var, vtgrids_files(i)%vt1d%file, &
               vtgrids(i)%vert_nom)
          END IF

          ! generate a 3d vertical coordinate field from the 1D field
          CALL check_vtgrid(vtgrids(i))
       END DO
    END IF


    ! Read state variable config
    ! The actual values are read by other subroutines in this module
    ! statevars_files needs to be explicitly saved by this module because
    ! it is used later by the ens read/write subroutines
    CALL parse_statedef(config, statevars, self%statevars_files)

  END SUBROUTINE stateio_nc_init
  !================================================================================




  !================================================================================
  !>
  !! TODO, modify so that we can read/scatter the state 1 variable and/or 1 level
  !! at a time in order to reduce the peak memory per compute node
  !--------------------------------------------------------------------------------
  SUBROUTINE stateio_nc_read_state(self, ensmem, state_var, state_val)
    CLASS(letkf_stateio_nc) :: self
    INTEGER,      INTENT(in)  :: ensmem     !< ensemble member number
    CHARACTER(*), INTENT(in)  :: state_var  !< name of the state variable
    REAL, ALLOCATABLE, INTENT(out) :: state_val(:,:,:)

    CHARACTER(:), ALLOCATABLE :: filename, invar
    INTEGER :: nx, ny, nz, i, idx
    INTEGER :: ncid, dimids(4), varid
    TYPE(letkf_hzgrid_spec) :: hgrid
    TYPE(letkf_vtgrid_spec) :: vgrid

    ! find matching state definition
    DO idx = 1, SIZE(statevars)
       IF (state_var == statevars(idx)%name) EXIT
    END DO

    ! determine the filename to load in. Subsituting #ENSx# with the ensemble num
    filename=parse_ens_filename(self%statevars_files(idx)%input%file, ensmem)
    invar = self%statevars_files(idx)%input%var

    ! TODO make sure this file exists
    IF (self%verbose) &
         PRINT '(X, A,I0.3,A,I0.3,6A)',' PE ',pe_rank,' loading ens ', ensmem, ' "',&
         TRIM(state_var),'" from "',TRIM(invar),'" of "', TRIM(filename)

    ! get the x,y,z dimensions based on the specified grid
    hgrid = letkf_state_hzgrid_getdef(statevars(idx)%hzgrid)
    vgrid = letkf_state_vtgrid_getdef(statevars(idx)%vtgrid)

    ! load field
    ! TODO, handle 2D vars
    ! TODO, handle empty time variables
    CALL check(nf90_open(filename, nf90_nowrite, ncid), &
         TRIM(filename))
    CALL check(nf90_inq_varid(ncid, invar, varid), &
         TRIM(filename)//" "//TRIM(invar))
    CALL check(nf90_inquire_variable(ncid, varid, ndims=i))
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids))
    CALL check(nf90_inquire_dimension(ncid, dimids(1), len=nx))
    CALL check(nf90_inquire_dimension(ncid, dimids(2), len=ny))
    CALL check(nf90_inquire_dimension(ncid, dimids(3), len=nz))
    ALLOCATE(state_val(nx,ny,nz))
    CALL check(nf90_get_var(ncid, varid, state_val))
    CALL check(nf90_close(ncid))

  END SUBROUTINE stateio_nc_read_state
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE stateio_nc_write_state(self, ensmem, state_vals)
    CLASS(letkf_stateio_nc) :: self
    INTEGER,        INTENT(in) :: ensmem
    REAL,           INTENT(in) :: state_vals(:,:,:)

    INTEGER :: i, j, idx1, idx2
    CHARACTER(:), ALLOCATABLE :: str
    INTEGER :: ncid, d_t, d_x, d_y, d_z(100), varid

    LOGICAL :: var_unsaved(SIZE(statevars))
    LOGICAL :: var_tosave(SIZE(statevars))
    CHARACTER(:), ALLOCATABLE :: curr_filename, var_filename


    var_unsaved=.true.
    ! while there are variables that have not been saved yet
    !------------------------------------------------------------
    DO WHILE(ANY(var_unsaved))
       var_tosave=.false.
       curr_filename = ""

       ! look at all the variables, and determine the next batch to save
       ! (i.e. variables with the same output filename)
       DO i=1, SIZE(statevars)
          ! skip this variable if already saved
          IF(.NOT. var_unsaved(i)) CYCLE

          ! determine the filename for thie variable
          var_filename=parse_ens_filename(self%statevars_files(i)%output%file, ensmem)

          ! if current filename empty, init
          IF(curr_filename=="") curr_filename = var_filename

          ! if this variable's filename does equals current filename, add to list to save
          IF(curr_filename == var_filename) var_tosave(i) = .true.
       END DO


       ! initialize the file
       ! TODO: add support for multiple horizontal grids
       ! TODO handle 2D, 3D, and constant vertical grids
       IF (self%verbose) &
         PRINT '(X,A,I0.3,2A)', " PE ", pe_rank,' writing output for "'//curr_filename//'"'
       CALL check(nf90_create(curr_filename, nf90_hdf5, ncid))
       CALL check(nf90_def_dim(ncid, "time", 0, d_t))
       CALL check(nf90_def_var(ncid, "time", nf90_real, (/d_t/), varid))
       CALL check(nf90_put_att(ncid, varid, "axis", "T"))
       CALL check(nf90_def_dim(ncid, "lon", grid_nx, d_x))
       CALL check(nf90_def_var(ncid, "lon", nf90_real, (/d_x/), varid))
       CALL check(nf90_put_att(ncid, varid, "long_name", "longitude"))
       CALL check(nf90_put_att(ncid, varid, "units", "degrees_east"))
       CALL check(nf90_put_att(ncid, varid, "axis", "X"))
       CALL check(nf90_def_dim(ncid, "lat", grid_ny, d_y))
       CALL check(nf90_def_var(ncid, "lat", nf90_real, (/d_y/), varid))
       CALL check(nf90_put_att(ncid, varid, "long_name", "latitude"))
       CALL check(nf90_put_att(ncid, varid, "units", "degrees_north"))
       CALL check(nf90_put_att(ncid, varid, "axis", "Y"))
       DO i=1,SIZE(vtgrids)
          IF (vtgrids(i)%dims > 1)&
            CALL letkf_mpi_abort ("LETKF_STATE_NC: 2D/3D vertical grids not yet supported.")
          CALL check(nf90_def_dim(ncid, vtgrids(i)%name, SIZE(vtgrids(i)%vert, 1), d_z(i)))
          CALL check(nf90_def_var(ncid, vtgrids(i)%name, nf90_real, (/d_z(i)/), varid))
          CALL check(nf90_put_att(ncid, varid, "axis", "Z"))
       END DO


       ! for each variable to save in this file, add its definition to file
       DO i=1, SIZE(statevars)
          ! skip if not being saved in this file
          IF( .NOT. var_tosave(i)) CYCLE

          ! determine which vertical coordinate is the matching one
          str = statevars(i)%vtgrid
          DO j=1, SIZE(vtgrids)
             IF (vtgrids(j)%name == str) EXIT
          END DO

          ! create variable
          IF(j > SIZE(vtgrids))&
               CALL letkf_mpi_abort("cannot find vertical grid in stateio_nc_write_init")
          CALL check(nf90_def_var(ncid, self%statevars_files(i)%output%var,&
               nf90_real, (/d_x,d_y,d_z(j),d_t/),varid))
          CALL check(nf90_def_var_deflate(ncid, varid, 1, 1, self%compression))
       END DO

       ! done with file definitions
       CALL check(nf90_enddef(ncid))

       ! write out horizontal/vertical grids
       CALL check(nf90_inq_varid(ncid, "lon", varid))
       CALL check(nf90_put_var(ncid, varid, hzgrids(1)%lon_nom))
       CALL check(nf90_inq_varid(ncid, "lat", varid))
       CALL check(nf90_put_var(ncid, varid, hzgrids(1)%lat_nom))
       DO i =1, SIZE(vtgrids)
          CALL check(nf90_inq_varid(ncid, vtgrids(i)%name, varid))
          CALL check(nf90_put_var(ncid, varid, vtgrids(i)%vert_nom))
       END DO

       !TODO fill in the time variable

       ! write actual values of the variables, mark variable as saved
       DO i=1, SIZE(statevars)
          IF(.NOT. var_tosave(i)) CYCLE

          idx1 = statevars(i)%grid_s_idx
          idx2 = idx1 + statevars(i)%levels - 1

          ! TODO, check to make sure strings are valid
          CALL check(nf90_inq_varid(ncid, self%statevars_files(i)%output%var, varid))
          CALL check(nf90_put_var(ncid, varid, state_vals(:,:,idx1:idx2)))

          var_unsaved(i)=.false.
       END DO

       ! close the file
       CALL check(nf90_close(ncid))

    END DO


  END SUBROUTINE stateio_nc_write_state
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE read_nc_d2(varname, filename, array)
    CHARACTER(*), INTENT(in)  :: varname, filename
    REAL, ALLOCATABLE, INTENT(out) :: array(:,:)

    INTEGER :: i, j, l, nx, ny
    INTEGER, ALLOCATABLE :: dimids(:)
    INTEGER :: ncid, varid

    ! open file, find dimensions
    CALL check(nf90_open(filename, nf90_nowrite, ncid), &
         '"'//TRIM(filename)//'"')

    CALL check(nf90_inq_varid(ncid, varname, varid), &
         '"'//TRIM(varname)//'" in "'//TRIM(filename)//'"')

    CALL check(nf90_inquire_variable(ncid, varid, ndims=i))
    ALLOCATE( dimids(i) )

    ! if more than 2 dimensions, make sure all but the first 2 are of size 1
    ! (since a lot of files tend to include time dimension of size 1)
    IF ( i < 2) CALL letkf_mpi_abort("variable dimensions /= 2")
    IF ( i > 2 ) THEN
      CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids))
      DO j=i,3,-1
        CALL check(nf90_inquire_dimension(ncid, dimids(j), len=l))
        IF ( l > 1 ) CALL letkf_mpi_abort("variable dimensions /= 2")
      ENDDO
    ENDIF

    ! get the dimension sizes
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids))
    CALL check(nf90_inquire_dimension(ncid, dimids(1), len=nx))
    CALL check(nf90_inquire_dimension(ncid, dimids(2), len=ny))

    ! read in the variable
    ALLOCATE(array(nx,ny))
    CALL check(nf90_get_var(ncid, varid, array))

    ! cleanup
    CALL check(nf90_close(ncid))

  END SUBROUTINE read_nc_d2
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE read_nc_d1(varname, filename,  array)
    CHARACTER(*), INTENT(in)  :: varname, filename
    REAL, ALLOCATABLE, INTENT(out) :: array(:)

    INTEGER :: i, dimids(1), ni
    INTEGER :: ncid, varid

    ! open file, find dimensions
    ! TODO, handle a time dimension if stored that way
    CALL check(nf90_open(filename, nf90_nowrite, ncid), &
         '"'//TRIM(filename)//'"')

    CALL check(nf90_inq_varid(ncid, varname, varid), &
         '"'//TRIM(varname)//'" in "'//TRIM(filename)//'"')

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
  !================================================================================



  !================================================================================
  !> Helper function to check status of netcdf calls
  !--------------------------------------------------------------------------------
  SUBROUTINE check(status, str)
    INTEGER, INTENT(in) :: status
    CHARACTER(*), OPTIONAL, INTENT(in) :: str

    CHARACTER(:), ALLOCATABLE :: str2

    IF(status /= nf90_noerr) THEN
       str2=TRIM(nf90_strerror(status))

       IF(PRESENT(str)) THEN
          str2=str2//": "//TRIM(str)
       END IF

       CALL letkf_mpi_abort(str2)
    END IF
  END SUBROUTINE check
  !================================================================================


END MODULE letkf_stateio_nc_mod
