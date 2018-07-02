!================================================================================
!> The main LETKF library module
!================================================================================
MODULE letkf
  
  use letkf_config
  USE letkf_diag
  USE letkf_loc
  USE letkf_loc_novrt
  USE letkf_mpi
  USE letkf_obs
  USE letkf_obs_nc
  USE letkf_obs_test
  USE letkf_solver
  USE letkf_state
  USE letkf_state_nc
  USE timing
  USE getmem

  IMPLICIT NONE
  PRIVATE

  
  ! use this to get the repository version at compile time
#ifndef CVERSION
#define CVERSION "Unknown"
#endif

  
  !================================================================================
  !================================================================================
  ! Public module components
  !================================================================================
  !================================================================================

  PUBLIC :: letkf_init
  PUBLIC :: letkf_run
  PUBLIC :: letkf_register_hook



  !================================================================================
  !================================================================================
  ! Private module components
  !================================================================================
  !================================================================================

  type(configuration) :: config


  
CONTAINS



  !================================================================================
  !> Initialize the LETKF library.
  !! This needs to be done before any other user-called functions are performed.
  !! This initializes the MPI backend, and loads the default stateio, obsio, and
  !! localizer classes
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_init(config_filename)
    character(len=*) :: config_filename !< json configuration file to load

    ! temprary pointers for initializing and registering the default classes
    CLASS(letkf_obsio),     POINTER :: obsio_ptr
    CLASS(letkf_stateio),   POINTER :: stateio_ptr
    CLASS(letkf_localizer), POINTER :: localizer_ptr

    ! initialize the mpi backend
    CALL letkf_mpi_preinit()
    CALL getmem_init(pe_root, letkf_mpi_comm)

    ! initialize global timers
    CALL timing_init(letkf_mpi_comm, 0)
    CALL timing_start("pre-init")


    ! print out the welcome screen
    IF (pe_isroot) THEN
       PRINT *, "======================================================================"
       PRINT *, "======================================================================"
       PRINT *, "Universal Multi-Domain Local Ensemble Transform Kalman Filter"
       PRINT *, " (UMD-LETKF)"
       print *, " version:  ", CVERSION
       PRINT *, "======================================================================"
       PRINT *, "======================================================================"
    END IF

    ! load in the configuration file
    if (pe_isroot) then
       print *, "Using configuration file: " // trim(config_filename)
       print *, ""
    end if
    letkf_config_log = pe_isroot
    call letkf_config_loadfile(config_filename, config)


    ! initialize the rest of MPI
    ! (determines the processor distribution for I/O)
    ! TODO

    ! setup the default observation I/O classes
    ALLOCATE(obsio_nc :: obsio_ptr)
    CALL letkf_obs_register(obsio_ptr)
    ALLOCATE(obsio_test :: obsio_ptr)
    CALL letkf_obs_register(obsio_ptr)

    ! setup the default state I/O classes
    ALLOCATE(stateio_nc :: stateio_ptr)
    CALL letkf_state_register(stateio_ptr)

    ! setup the default localizer classes
    ALLOCATE(loc_novrt :: localizer_ptr)
    CALL letkf_loc_register(localizer_ptr)

    CALL timing_stop("pre-init")

  END SUBROUTINE letkf_init
  !================================================================================



  !================================================================================
  !> Does all the fun stuff of actually running the LETKF
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_run()

    ! LETKF initialization
    ! module initialization order is somewhat important
    ! "obs" must be initialized before "state" because, in the case of the test obs
    ! reader, the "obs" module might create hooks to get information from the state
    ! before "obs_read" is called.
    CALL timing_start("init", TIMER_SYNC)

    CALL letkf_mpi_init(config%get_child("mpi"))
    call letkf_mpi_barrier()

    CALL letkf_state_init(config%get_child("state"))
    call letkf_mpi_barrier()

    CALL letkf_obs_init(config%get_child("observation"))
    call letkf_mpi_barrier()

    CALL letkf_obs_read()
    call letkf_mpi_barrier()

    CALL letkf_loc_init(config%get_child("localization"))
    call letkf_mpi_barrier()

    call letkf_solver_init(config%get_child("solver"))
    call letkf_mpi_barrier()

    CALL timing_stop("init")

    ! run the LETKF solver
    CALL letkf_solver_run()
    CALL letkf_mpi_barrier()


    ! LETKF final output
    ! ------------------------------------------------------------
    CALL timing_start("output", TIMER_SYNC)
    IF(pe_isroot) THEN
       PRINT *, ""
       PRINT *, "Saving output..."
    END IF

    ! miscellaneous diagnostics from the LETKF solver
    call letkf_diag_write()

    ! ensemble mean, spread
    call letkf_state_write_meansprd("ana")

    ! save ensemble members
    call letkf_state_write_ens()

    CALL timing_stop("output")



    ! all done
    call letkf_mpi_barrier()
    CALL timing_print()
    CALL getmem_print()
    CALL letkf_mpi_final()

  END SUBROUTINE letkf_run
  !================================================================================


  !================================================================================
  !> Register a callback function for one of the hooks where
  !! users can implement custom functionality
  !! @TODO implement this
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_register_hook
    PRINT *, "NOT yet implemented (letkf_register_hook)"
    STOP 1
  END SUBROUTINE letkf_register_hook
  !================================================================================


  
END MODULE letkf
