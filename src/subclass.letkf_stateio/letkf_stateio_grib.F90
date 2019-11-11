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
!> grib based state I/O
!!
!================================================================================
MODULE letkf_stateio_grib_mod
  USE letkf_config
  USE letkf_state
  USE letkf_state_helper
  USE letkf_mpi
  USE wgrib2api

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
  TYPE, PUBLIC, EXTENDS(letkf_stateio) :: letkf_stateio_grib
     TYPE(statevar_files), ALLOCATABLE :: statevars_files(:)
     INTEGER :: grib2_regex
   CONTAINS
     PROCEDURE, NOPASS :: name => stateio_grib_get_name
     PROCEDURE, NOPASS :: desc => stateio_grib_get_desc
     PROCEDURE         :: init => stateio_grib_init
     PROCEDURE         :: read_state  => stateio_grib_read_state
     PROCEDURE         :: write_state => stateio_grib_write_state
  END TYPE letkf_stateio_grib
  !================================================================================


  CHARACTER(len=1024) :: grib2_inv ='@mem:'
  INTEGER(kind=8), PARAMETER :: GRIB2_BUFFER_SIZE = 1024*1024*10 ! 10Mb is enough, right?
  CHARACTER(len=GRIB2_BUFFER_SIZE) :: grib2_buffer

CONTAINS



  !================================================================================
  !> Get the unique name of the state I/O class
  !--------------------------------------------------------------------------------
  FUNCTION stateio_grib_get_name() RESULT(name)
    CHARACTER(:), ALLOCATABLE :: name
    name = "STATEIO_GRIB"
  END FUNCTION stateio_grib_get_name
  !================================================================================




  !================================================================================
  !> Get the description of the state I/O class
  !--------------------------------------------------------------------------------
  FUNCTION stateio_grib_get_desc() RESULT(desc)
    CHARACTER(:), ALLOCATABLE :: desc
    desc = "grib2 formatted model state I/O"
  END FUNCTION stateio_grib_get_desc
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE stateio_grib_init(self, config, hzgrids, vtgrids, statevars)
    CLASS(letkf_stateio_grib) :: self
    TYPE(configuration), INTENT(in) :: config
    TYPE(letkf_hzgrid_spec),   ALLOCATABLE, INTENT(out) :: hzgrids(:)
    TYPE(letkf_vtgrid_spec),   ALLOCATABLE, INTENT(out) :: vtgrids(:)
    TYPE(letkf_statevar_spec), ALLOCATABLE, INTENT(out) :: statevars(:)

    TYPE(hzgrid_files), ALLOCATABLE :: hzgrids_files(:)
    TYPE(vtgrid_files), ALLOCATABLE :: vtgrids_files(:)

    INTEGER :: ret, i, nx, ny

    IF (pe_isroot) THEN
       PRINT '(/A)', ""
       PRINT *, " letkf_stateio_grib_init() : "
       PRINT *, "------------------------------------------------------------"
    END IF

    ! load in our section of the configuration
    CALL config%get("grib_regex", self%grib2_regex, 0)

    ! create memory files in wgrib api for later use
    ret = wgrib2_set_mem_buffer(grib2_buffer, GRIB2_BUFFER_SIZE, 1)
    WRITE (grib2_inv, '(A,I0)') '@mem:', 1

    ! read the horizontal configuration
    CALL parse_hzgrids(config, hzgrids, hzgrids_files)

    ! only the root PE needs to read in the actual horizontal grid
    ! letkf_state will handle scattering across the PEs
    IF (pe_isroot) THEN
      DO i = 1, SIZE(hzgrids)

          ! load lats, lons
          ret = grb2_mk_inv(hzgrids_files(i)%lat%file, grib2_inv)
          if (ret/=0) call letkf_mpi_abort(&
            "Cannot make grib2 inventory for "//hzgrids_files(i)%lat%file)

          ret = grb2_inq(hzgrids_files(i)%lat%file, grib2_inv, &
            hzgrids_files(i)%lat%var, lat=hzgrids(i)%lat, lon=hzgrids(i)%lon, &
            regex=self%grib2_regex)
          if (ret <= 0) call letkf_mpi_abort(&
            "Cannot find variable in grib file")

          ! check the grids
          CALL check_hzgrid(hzgrids(i))
          nx = SIZE(hzgrids(i)%lon_nom)
          ny = SIZE(hzgrids(i)%lat_nom)

          ! TODO read in a mask
          ALLOCATE(hzgrids(i)%mask(nx,ny))
          hzgrids(i)%mask = .FALSE.

       END DO
    END IF


    ! read vertical grid configuration
    CALL parse_vtgrids(config, vtgrids, vtgrids_files)

    ! only the root PE needs to read in the actual vertical grid
    ! letkf_state will handle broadcasting to the other PEs
    IF (pe_isroot) THEN
       DO i = 1, SIZE(vtgrids)

          ! TODO do this correctly
          ! do a temporary 1 level vertical coord for now
          ALLOCATE(vtgrids(i)%vert_nom(1))
          vtgrids(i)%vert_nom = 1.0

          ! generate a 3d vertical coordinate field from the 1D field
          CALL check_vtgrid(vtgrids(i))
       END DO
    END IF

    ! Read state variable config
    ! The actual values are read by other subroutines in this module
    ! statevars_files needs to be explicitly saved by this module because
    ! it is used later by the ens read/write subroutines
    CALL parse_statedef(config, statevars, self%statevars_files)

  END SUBROUTINE stateio_grib_init
  !================================================================================




  !================================================================================
  !>
  !! TODO, modify so that we can read/scatter the state 1 variable and/or 1 level
  !! at a time in order to reduce the peak memory per compute node
  !--------------------------------------------------------------------------------
  SUBROUTINE stateio_grib_read_state(self, ensmem, state_var, state_val)
    CLASS(letkf_stateio_grib) :: self
    INTEGER,      INTENT(in)  :: ensmem     !< ensemble member number
    CHARACTER(*), INTENT(in)  :: state_var  !< name of the state variable
    REAL, ALLOCATABLE, INTENT(out) :: state_val(:,:,:)

    CHARACTER(:), ALLOCATABLE :: filename, invar
    INTEGER :: nx, ny, nz, i, idx, ret, z
    TYPE(letkf_hzgrid_spec) :: hgrid
    TYPE(letkf_vtgrid_spec) :: vgrid

    REAL, ALLOCATABLE :: tmp2d(:,:)

    ! find matching state definition
    DO idx = 1, SIZE(statevars)
       IF (state_var == statevars(idx)%name) EXIT
    END DO

    ! determine the filename to load in. Subsituting #ENSx# with the ensemble num
    ! TODO replace "var name" with list of match strings
    filename=parse_ens_filename(self%statevars_files(idx)%input%file, ensmem)
    invar = trim(self%statevars_files(idx)%input%var)

    ! TODO make sure this file exists
    IF (self%verbose) &
         PRINT '(X, A,I0.3,A,I0.3,6A)',' PE ',pe_rank,' loading ens ', ensmem, ' "',&
         TRIM(state_var),'" from "',TRIM(invar),'" of "', TRIM(filename)

    ! get the x,y,z dimensions based on the specified grid
    hgrid = letkf_state_hzgrid_getdef(statevars(idx)%hzgrid)
    vgrid = letkf_state_vtgrid_getdef(statevars(idx)%vtgrid)


    ! load field
    ! TODO, modify to handle 3d variables
    ret = grb2_mk_inv(filename, grib2_inv)
    IF (ret /= 0) CALL letkf_mpi_abort(&
      "unable to create grib2 inventory for file "//filename)

    ret = grb2_inq(filename, grib2_inv, invar, nx=nx, ny=ny, &
      regex=self%grib2_regex)
    IF (ret <= 0) CALL letkf_mpi_abort(&
      " variable does not exist "//invar)

    ret = grb2_inq(filename, grib2_inv, invar, nx=nx, ny=ny, data2=tmp2d, &
      regex=self%grib2_regex)
    IF (ret <= 0) CALL letkf_mpi_abort(&
      "unable to read data from file "//filename)

    ALLOCATE(state_val(nx,ny,1))
    state_val(:,:,1) = tmp2d
    DEALLOCATE(tmp2d)

    ret = grb2_free_file(filename)
    IF (ret /= 0) CALL letkf_mpi_abort( &
      "unable to release file "//filename)

  END SUBROUTINE stateio_grib_read_state
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE stateio_grib_write_state(self, ensmem, state_vals)
    CLASS(letkf_stateio_grib) :: self
    INTEGER,        INTENT(in) :: ensmem
    REAL,           INTENT(in) :: state_vals(:,:,:)

    CHARACTER(:), ALLOCATABLE :: filename_out, filename_prev
    CHARACTER(:), ALLOCATABLE :: filename_tmpl
    CHARACTER(:), ALLOCATABLE :: meta_str

    REAL, ALLOCATABLE :: tmp2d(:,:)
    INTEGER :: ret, idx1, idx2, i, z

    ALLOCATE(tmp2d(grid_nx, grid_ny))
    filename_prev = ""

    ! For each state variable
    DO i = 1, SIZE(statevars)
       idx1 = statevars(i)%grid_s_idx
       idx2 = idx1 + statevars(i)%levels - 1

       ! determine the filename to use as a template (from the input file for now)
       ! TODO, add option for specifying a different template file
       filename_tmpl = parse_ens_filename(self%statevars_files(i)%input%file, 1)

       ! determine the output filename
       filename_out = parse_ens_filename(self%statevars_files(i)%output%file, ensmem)
       IF (self%verbose .AND. filename_out /= filename_prev) THEN
         ! If this is a different filename than the one for the previous variable,
         ! print out the filename
         PRINT '(X,A,I0.3,2A)', " PE ", pe_rank,' writing output for "'//filename_out//'"'
         filename_prev = filename_out
       END if

       ! determine the meta data string of the output variable, based on the template
       ret = grb2_mk_inv(filename_tmpl, grib2_inv)
       IF (ret/=0) CALL letkf_mpi_abort(&
         "Cannot make grib2 inventory for "//filename_tmpl)
       ALLOCATE(CHARACTER(1024) :: meta_str)
       ret = grb2_inq(filename_tmpl, grib2_inv, &
         self%statevars_files(i)%output%var, desc=meta_str, &
         regex=self%grib2_regex)
       IF (ret <= 0) CALL letkf_mpi_abort( &
         "Cannot get metadata for "//TRIM(self%statevars_files(i)%output%var))

       ! write out the levels of this variable
       ! TODO make sure this works for 3d variables
       DO z=1,statevars(i)%levels
         tmp2d = state_vals(:,:,idx1+z-1)
         ret = grb2_wrt(filename_out, filename_tmpl, z, &
                        meta=meta_str, &
                        data2=tmp2d)
         IF (ret /= 0) call letkf_mpi_abort( &
            "unable to write to file "//filename_out)
       END DO

    END DO
    DEALLOCATE(tmp2d)

  END SUBROUTINE stateio_grib_write_state
  !================================================================================


END MODULE letkf_stateio_grib_mod
