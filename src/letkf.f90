module letkf
  use timing
  use mpi
  use letkf_obs
  use letkf_mpi

  use letkf_obs_dat
  
  implicit none

  type :: letkf_solver
   contains
     procedure :: run => letkf_solver_run
     procedure :: read_config => letkf_solver_read_config
     procedure :: read_obs => letkf_solver_read_obs
  end type letkf_solver


  integer :: mem
  integer, parameter :: dp = kind(0.0d0)
  
contains

  
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

  
  
  subroutine letkf_solver_read_obs(self)
    class(letkf_solver) :: self
    integer :: ierr

    integer :: i, j
    integer :: timer, timer2, timer3

    character(len=1024) :: filename    
    type(obsio_dat) :: reader
    
    type(observation), allocatable :: obs(:)    
    real(dp), allocatable :: obs_inov(:,:)
    integer, allocatable :: obs_qc(:,:)

    type(observation), allocatable :: obs_t(:)
    real(dp), allocatable :: obs_inov_t(:)
    integer, allocatable :: obs_qc_t(:)

    integer, allocatable :: obstat_count(:,:)
    integer, allocatable :: obplat_count(:,:)

    timer  = timer_init("obs read", TIMER_SYNC)
    timer2 = timer_init("obs read I/O")
    timer3 = timer_init("obs read MPI", timer_sync)
    call timer_start(timer)

    if (pe_isroot) then
       print *, ""
       print *, "Reading Observations"
       print *, "============================================================"
    end if

    ! parallel read of the observation innovation files for each ensemble member
    do i=1,size(ens_list)
       write (filename, '(A,I0.3,A)') 'obsop/2005100100_atm',ens_list(i),'.dat'

       ! read the file
       !TODO, not the most efficient, general observation information is not
       ! needed with every single ensemble member, this should be in a separate
       ! file, with each ens member file only containing the obs space value and qc
       call timer_start(timer2)       
       call reader%read(filename, obs_t, obs_inov_t, obs_qc_t)
       call timer_stop(timer2)

       ! create the obs storage if it hasn't already been done
       if (.not. allocated(obs)) then
          allocate(obs(size(obs_t)))
          obs = obs_t
          allocate(obs_inov(mem, size(obs_inov_t)))
          allocate(obs_qc  (mem, size(obs_qc_t)))
       end if

       ! copy the per ensemble member innovation and qc
       obs_inov(ens_list(i), :) = obs_inov_t(:)
       obs_qc  (ens_list(i), :) = obs_qc_t(:)

       !cleanup
       deallocate(obs_t)
       deallocate(obs_inov_t)
       deallocate(obs_qc_t)
    end do

    
    ! distribute the qc and innovation values
    ! TODO, using MPI_SUM is likely inefficient, do this another way
    call timer_start(timer3)    
    call mpi_allreduce(mpi_in_place, obs_inov, mem*size(obs), mpi_real8, mpi_sum, mpi_comm_letkf, ierr)
    call mpi_allreduce(mpi_in_place, obs_qc  , mem*size(obs), mpi_integer, mpi_sum, mpi_comm_letkf, ierr)
    call timer_stop(timer3)


    ! print statistics about the observations
    if (pe_isroot) then
       print '(I9,A)', size(obs), " observations loaded"

       allocate(obstat_count (size(obsdef_list)+1,  3))
       allocate(obplat_count (size(platdef_list)+1, 3))
       
       obstat_count = 0
       obplat_count = 0

       ! count
       do i=1,size(obs)
          do j=1,size(obsdef_list)
             if (obs(i)%id == obsdef_list(j)%id) then
                obstat_count(j,1) = obstat_count(j,1) + 1
                if (all( obs_qc(:,i) == 1)) &
                     obstat_count(j,2) = obstat_count(j,2) + 1
                exit
             end if
          end do
          if (j > size(obsdef_list)) &
               obstat_count(size(obstat_count),1) = obstat_count(size(obstat_count),1) + 1
          do j=1,size(platdef_list)
             if (obs(i)%plat == platdef_list(j)%id) then
                obplat_count(j,1) = obplat_count(j,1) + 1
                if (all( obs_qc(:,i) == 1)) &
                     obplat_count(j,2) = obplat_count(j,2) + 1                
                exit
             end if
          end do
          if (j > size(platdef_list)) &
               obplat_count(size(obplat_count),1) = obplat_count(size(obplat_count),1) + 1
       end do

       ! print stats
       print *, ""
       print '(A7,A8,A8)', '','total','good'
       do i=1,size(obstat_count,1)
          if (obstat_count(i,1) == 0) cycle
          if (i < size(obstat_count)) then
             print '(A7,I8,I8,A2,F5.1,A)',&
                  obsdef_list(i)%name_short, obstat_count(i,1), obstat_count(i,2), &
                  '(',real(obstat_count(i,2))/obstat_count(i,1)*100, ')%'
          else
             print '(A7,I8)', 'unknown', obstat_count(i,1)
          end if
       end do

       print *, ""
       print '(A7,A8,A8)', '','total','good'
       do i=1,size(obplat_count,1)
          if (obplat_count(i,1) == 0) cycle
          if (i < size(obplat_count)) then
             print '(A7,I8,I8,A2,F5.1,A)',&
                  platdef_list(i)%name_short, obplat_count(i,1), obplat_count(i,2), &
                  '(',real(obplat_count(i,2))/obplat_count(i,1)*100, ')%'
          else
             print '(A7,I8)', 'unknown', obplat_count(i,1)
          end if
       end do
       
       !cleanup
       deallocate(obstat_count)
       deallocate(obplat_count)
    end if
    

    call timer_stop(timer)
  end subroutine letkf_solver_read_obs
  
end module letkf
