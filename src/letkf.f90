module letkf
  use timing
  use mpi
  use letkf_obs
  use letkf_mpi
  use global
  use letkf_obs_dat
  
  implicit none

  type :: letkf_solver
   contains
     procedure :: run => letkf_solver_run
     procedure :: read_config => letkf_solver_read_config
     procedure :: read_obs => letkf_solver_read_obs
  end type letkf_solver


  
contains


  !############################################################
  subroutine letkf_solver_run(self)
    class(letkf_solver) :: self
    integer :: t_total, t_init, t_letkf


    namelist /letkf_settings/ mem
    
    t_total = timer_init("Total Runtime")
    call timer_start(t_total)

    ! Initialize
    t_init = timer_init("Initialization")
    call timer_start(t_init)

    open(90, file="letkf.nml")
    read(90, nml=letkf_settings)
    close(90)
    
    call letkf_mpi_init()
    if (pe_isroot) then
       print *, "============================================================"
       print *, " Unified Multi-Domain Local Ensemble Transform Kalman Filter"
       print *, " (UMD-LETKF)"
       print *, " version 0.0.0"
       print *, "============================================================"
       print *, ""
       print letkf_settings
    end if
    call letkf_mpi_init2(mem)
    
    call self%read_config
    call self%read_obs
    call timer_stop(t_init)
    
    ! run LETKF
    t_letkf = timer_init("core", TIMER_SYNC)
    call timer_start(t_letkf)
    call timer_stop(t_letkf)
    
    ! all done
    
    call timer_stop(t_total)

    call timer_print()
    
    call letkf_mpi_end()

    
  end subroutine letkf_solver_run


  
  !############################################################
  subroutine letkf_solver_read_config(self)
    class(letkf_solver) :: self
    logical :: ex
    integer :: ierr

    integer :: timer

    timer = timer_init("read config", TIMER_SYNC)
    call timer_start(timer)
    

    call letkf_obs_init("obsdef.cfg", "platdef.cfg")
    call timer_stop(timer)
    
  end subroutine letkf_solver_read_config

  

  !############################################################  
  subroutine letkf_solver_read_obs(self)
    class(letkf_solver) :: self
    type(obsio_dat) :: ioclass

    call letkf_obs_read(ioclass)
  end subroutine letkf_solver_read_obs
  
end module letkf
