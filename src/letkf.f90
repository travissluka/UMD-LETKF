module letkf
  use iso_fortran_env
  use timing
  use letkf_obs
  implicit none

  type :: letkf_solver
   contains
     procedure :: run => letkf_solver_run
     procedure :: read_config => letkf_solver_read_config
     procedure :: read_obs => letkf_solver_read_obs
  end type letkf_solver

  
contains

  
  subroutine letkf_solver_run(self)
    class(letkf_solver) :: self
    print *, "============================================================"
    print *, " Unified Multi-Domain Local Ensemble Transform Kalman Filter"
    print *, " (UMD-LETKF)"
    print *, " version 0.0.0"
    print *, "============================================================"
    print *, ""

    call self%read_config
    call self%read_obs
    call timing_print()
  end subroutine letkf_solver_run


  
  subroutine letkf_solver_read_config(self)
    class(letkf_solver) :: self
    logical :: ex

    
    call timing_start("config")

    print *, ""
    print *, "Reading LETKF configuration"
    print *, "==========================================================="

    call letkf_obs_init("obsdef.cfg", "platdef.cfg")
    ! inquire(file=config_file, exist=ex)
    ! if (.not. ex) then
    !    write (error_unit, *) "ERROR: LETKF configuration file does not exist: ", trim(config_file)
    !    stop 1
    ! end if

!    open(99, file=config_file, status="old")
!    read(99, nml=master)
!    close(99)
!    print master
    
    call timing_end("config")
  end subroutine letkf_solver_read_config

  
  
  subroutine letkf_solver_read_obs(self)
    class(letkf_solver) :: self
    call timing_start("obs")

    print *, "Reading Observations"
    print *, "==========================================================="
    
    call timing_end("obs")    
  end subroutine letkf_solver_read_obs
  
end module letkf
