!================================================================================
!> grib based state I/O
!================================================================================
MODULE letkf_state_grib
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
  TYPE, PUBLIC, EXTENDS(letkf_stateio) :: stateio_grib

   CONTAINS

     PROCEDURE, NOPASS :: name => stateio_grib_get_name
     PROCEDURE, NOPASS :: desc => stateio_grib_get_desc
     PROCEDURE         :: init => stateio_grib_init
     PROCEDURE         :: read_state  => stateio_grib_read_state
     PROCEDURE         :: write_state => stateio_grib_write_state
  END TYPE stateio_grib
  !================================================================================




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
    CLASS(stateio_grib) :: self
    TYPE(configuration), INTENT(in) :: config
    TYPE(letkf_hzgrid_spec),   ALLOCATABLE, INTENT(out) :: hzgrids(:)
    TYPE(letkf_vtgrid_spec),   ALLOCATABLE, INTENT(out) :: vtgrids(:)
    TYPE(letkf_statevar_spec), ALLOCATABLE, INTENT(out) :: statevars(:)

    IF (pe_isroot) THEN
       PRINT '(/A)', ""
       PRINT *, " letkf_stateio_grib_init() : "
       PRINT *, "------------------------------------------------------------"
    END IF


  END SUBROUTINE stateio_grib_init
  !================================================================================



  
  !================================================================================
  !>
  !! TODO, modify so that we can read/scatter the state 1 variable and/or 1 level
  !! at a time in order to reduce the peak memory per compute node
  !--------------------------------------------------------------------------------
  SUBROUTINE stateio_grib_read_state(self, ensmem, state_var, state_val)
    CLASS(stateio_grib) :: self
    INTEGER,      INTENT(in)  :: ensmem     !< ensemble member number
    CHARACTER(*), INTENT(in)  :: state_var  !< name of the state variable
    REAL, ALLOCATABLE, INTENT(out) :: state_val(:,:,:)

    CHARACTER(:), ALLOCATABLE :: filename, invar
    INTEGER :: nx, ny, nz, i
    TYPE(letkf_hzgrid_spec) :: hgrid
    TYPE(letkf_vtgrid_spec) :: vgrid
    TYPE(letkf_statevar_spec) :: statevar

    ! ! find matching state definition
    ! statevar = letkf_state_var_getdef(state_var)

    ! ! determine the filename to load in. Subsituting #ENSx# with the ensemble num
    ! filename=parse_ens_filename(statevar%input_file, ensmem)
    ! invar = statevar%input_var

    ! ! TODO make sure this file exists
    ! IF (self%verbose) &
    !      PRINT '(X, A,I0.3,A,I0.3,6A)',' PE ',pe_rank,' loading ens ', ensmem, ' "',&
    !      TRIM(state_var),'" from "',TRIM(invar),'" of "', TRIM(filename)

    ! ! get the x,y,z dimensions based on the specified grid
    ! hgrid = letkf_state_hzgrid_getdef(statevar%hzgrid)
    ! vgrid = letkf_state_vtgrid_getdef(statevar%vtgrid)

    ! ! load field
    ! ! TODO Load the state variable  "invar" from file "filename"
 
  END SUBROUTINE stateio_grib_read_state
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE stateio_grib_write_state(self, ensmem, state_vals)
    CLASS(stateio_grib) :: self
    INTEGER,        INTENT(in) :: ensmem
    REAL,           INTENT(in) :: state_vals(:,:,:)

    CHARACTER(:), ALLOCATABLE :: filename


    ! determine the filename to save
    ! TODO use the correct file template
    filename = parse_ens_filename("#TYPE#.#ENS4#.nc", ensmem)
    IF (self%verbose) &
         PRINT '(X,A,I0.3,2A)', " PE ", pe_rank,' writing output for "'//filename//'"'

    ! TODO write the state to file

  END SUBROUTINE stateio_grib_write_state
  !================================================================================

END MODULE letkf_state_grib
