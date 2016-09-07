module letkf
  use timing
  use mpi
  use letkf_obs
  use letkf_mpi
  
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

    call letkf_mpi_init()

    if (mpi_rank == 0) then
       print *, "============================================================"
       print *, " Unified Multi-Domain Local Ensemble Transform Kalman Filter"
       print *, " (UMD-LETKF)"
       print *, " version 0.0.0"
       print *, "============================================================"
       print *, ""
    end if

    call self%read_config
    call self%read_obs
    if (mpi_rank == 0) call timing_print()
    
    call letkf_mpi_end()
  end subroutine letkf_solver_run


  
  subroutine letkf_solver_read_config(self)
    class(letkf_solver) :: self
    logical :: ex
    integer :: ierr

    
    call timing_start("config")
    if (mpi_rank == 0) then
       print *, ""
       print *, "Reading LETKF configuration"
       print *, "============================================================"
       print *, ""       
    end if
    call mpi_barrier(mpi_comm_world,ierr)

    call letkf_obs_init("obsdef.cfg", "platdef.cfg")
    call timing_end("config")
  end subroutine letkf_solver_read_config

  
  
  subroutine letkf_solver_read_obs(self)
    class(letkf_solver) :: self
    integer :: ierr
    
    call timing_start("obs")

    if (mpi_rank == 0) then
       print *, ""
       print *, "Reading Observations"
       print *, "============================================================"
       print *, ""       
    end if
    call mpi_barrier(mpi_comm_world,ierr)
    
    call timing_end("obs")    
  end subroutine letkf_solver_read_obs
  
end module letkf
