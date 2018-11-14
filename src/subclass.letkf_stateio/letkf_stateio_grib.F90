!================================================================================
!> grib based state I/O
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
   CONTAINS

     PROCEDURE, NOPASS :: name => stateio_grib_get_name
     PROCEDURE, NOPASS :: desc => stateio_grib_get_desc
     PROCEDURE         :: init => stateio_grib_init
     PROCEDURE         :: read_state  => stateio_grib_read_state
     PROCEDURE         :: write_state => stateio_grib_write_state
  END TYPE letkf_stateio_grib
  !================================================================================


  CHARACTER(len=1024) :: grib2_inv ='@mem:'
  INTEGER(kind=8), PARAMETER :: GRIB2_BUFFER_SIZE = 1024*1024*10
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
    desc = "grib formatted model state I/O"
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

    ! create memory files for later use
    ret = wgrib2_set_mem_buffer(grib2_buffer, GRIB2_BUFFER_SIZE, 1)
    WRITE (grib2_inv, '(A,I0)') '@mem:', 1


    ! read the horizontal configuration
    CALL parse_hzgrids(config, hzgrids, hzgrids_files)

    ! onlyt the root PE needs to read in the actual horizontal grid
    ! letkf_state will handle scattering across the PEs
    IF (pe_isroot) THEN
      DO i = 1, SIZE(hzgrids)

          ! TODO, clean this up, use the correct lat/lon/mask files
          ! TODO, use the correct mask file

          ! load lats, lons
          ret = grb2_mk_inv(hzgrids_files(i)%lat%file, grib2_inv)
          if (ret/=0) call letkf_mpi_abort(&
            "Cannot make grib2 inventory for "//hzgrids_files(i)%lat%file)

          ret = grb2_inq(hzgrids_files(i)%lat%file, grib2_inv, "HTSGW", lat=&
            hzgrids(i)%lat, lon=hzgrids(i)%lon)
          if (ret == 0) call letkf_mpi_abort(&
            "Cannot find variable in grib file")

          CALL check_hzgrid(hzgrids(i))
          nx = SIZE(hzgrids(i)%lon_nom)
          ny = SIZE(hzgrids(i)%lat_nom)

          ! TODO read in mask
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
    invar = self%statevars_files(idx)%input%var

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

    ret = grb2_inq(filename, grib2_inv, invar, nx=nx, ny=ny, data2=tmp2d)
    IF (ret == 0) CALL letkf_mpi_abort(&
      "unable to read data from file "//filename)

    ALLOCATE(state_val(nx,ny,1))
    state_val(:,:,1) = tmp2d
    DEALLOCATE(tmp2d)

    ret = grb2_free_file(filename)

  END SUBROUTINE stateio_grib_read_state
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE stateio_grib_write_state(self, ensmem, state_vals)
    CLASS(letkf_stateio_grib) :: self
    INTEGER,        INTENT(in) :: ensmem
    REAL,           INTENT(in) :: state_vals(:,:,:)

    CHARACTER(:), ALLOCATABLE :: filename
    CHARACTER(:), ALLOCATABLE :: filename_tmpl

    REAL, ALLOCATABLE :: tmp2d(:,:)
    INTEGER :: ret, idx1, idx2, i, z, c

    ! determine the filename to save
    ! TODO use the correct file template
    filename = parse_ens_filename("#TYPE#.#ENS4#.grib2", ensmem)
    IF (self%verbose) &
         PRINT '(X,A,I0.3,2A)', " PE ", pe_rank,' writing output for "'//filename//'"'

    ! TODO write the state to file
    ! write out the state variables
    filename_tmpl = 'bkg/gwes01.glo_30m.t12z.grib2'
    ALLOCATE(tmp2d(grid_nx, grid_ny))
    DO i = 1, SIZE(statevars)
       idx1 = statevars(i)%grid_s_idx
       idx2 = idx1 + statevars(i)%levels - 1

       DO z=1,statevars(i)%levels
         tmp2d = state_vals(:,:,idx1+z-1)
         ret = grb2_wrt(filename, filename_tmpl, z, &
                        var=self%statevars_files(i)%output%var, &
                        data2=tmp2d)
       END DO
    END DO
    DEALLOCATE(tmp2d)

    ret = grb2_free_file(filename)


  END SUBROUTINE stateio_grib_write_state
  !================================================================================



  !================================================================================
  !> @brief Subroutine checking the return status from the grb2_mk_inv
  !! and grb2_inq.
  !! @param[in] iret return status
  !<
  SUBROUTINE check (iret)
    INTEGER, INTENT(in) :: iret
    CHARACTER (len=25),PARAMETER ::myname='check_iret'
    IF (iret/=1) THEN
       IF (iret<=0) THEN
          PRINT*,TRIM(myname),': No match! Check input to grb2'
       ELSE IF (iret>1)THEN
          PRINT*,TRIM(myname),'There are ', iret,' messages that matched'
       END IF
       PRINT*, 'Instant Termination!'
       STOP
    END IF
  END SUBROUTINE check
  !================================================================================


END MODULE letkf_stateio_grib_mod
