module letkf
  use timing
  use letkf_common
  use letkf_obs
  use letkf_mpi
  use letkf_core
  use letkf_state
  use letkf_loc

  implicit none
  private

  public :: letkf_driver_init
  public :: letkf_driver_run



contains



  subroutine letkf_driver_init()
    !! Initialize the LETKF module. This must be called before anything else.
    call letkf_mpi_preinit()
    call timing_init(mpi_comm_letkf, pe_root)
    if (isroot) then
       print "(A)", "============================================================"
       print "(A)", " Universal Multi-Domain Local Ensemble Transform Kalman Filter"
       print "(A)", " (UMD-LETKF)"
       print "(A)", " version 0.0.0"
       print "(A)", " Travis Sluka (tsluka@umd.edu, travis.sluka@noaa.gov)"
       print "(A)", "============================================================"
       print "(A)", ""
    end if
  end subroutine letkf_driver_init



  subroutine letkf_driver_run()
    real, allocatable :: wrk(:,:,:)
    real, allocatable :: wrk2(:,:,:,:)

    character(len=1024) :: filename
    integer :: t_total, t_init, t_letkf
    integer :: unit
    integer :: m, i

    namelist /letkf_settings/ mem, grid_nx, grid_ny


    ! Initialize timers
    t_total = timer_init("Total")
    t_init = timer_init("Initialization")
    call timer_start(t_total)
    call timer_start(t_init)


    ! read in main section of the  namelist
    open(newunit=unit, file=nml_filename)
    read(unit, nml=letkf_settings)
    close(unit)
    if (isroot) then
       print letkf_settings
    end if


    ! initialize individual modules
    call letkf_mpi_init(mem, grid_nx*grid_ny)
    call letkf_obs_init("obsdef.cfg", "platdef.cfg")
    call letkf_state_init()
    call letkf_core_init(mem)
    call timer_stop(t_init)


    ! ensure output directory is setup right
    call system('mkdir -p OUTPUT')


    ! run LETKF core
    if(isroot) then
       print *, new_line('a'),&
            new_line('a'), '============================================================',&
            new_line('a'), '    Running LETKF core',&
            new_line('a'), '============================================================'
    end if
    t_letkf = timer_init("loop", TIMER_SYNC)
    call timer_start(t_letkf)
    call letkf_do_letkf()
    call timer_stop(t_letkf)


    ! gather the analysis mean/sprd and write out
    allocate(wrk(grid_nx, grid_ny, grid_ns))
    do i=1,grid_ns
       call letkf_mpi_ij2grd(ana_mean_ij(:,i), wrk(:,:,i))
    end do
    if (isroot) call stateio_class%write('OUTPUT/ana_mean.nc', wrk)
    do i=1,grid_ns
       call letkf_mpi_ij2grd(ana_sprd_ij(:,i), wrk(:,:,i))
    end do
    if (isroot) call stateio_class%write('OUTPUT/ana_mean.nc', wrk)
    deallocate(wrk)


    ! write the analysis ensemble members
    allocate(wrk2(grid_nx, grid_ny, grid_ns, size(ens_list)))
    if(isroot) print *, "Writing analysis ensemble members..."
    call letkf_mpi_ij2ens(ana_ij, wrk2)
    do m=1,size(ens_list)
       write (filename, '(A,I0.4,A)') 'OUTPUT/',ens_list(m),'.nc'
       call stateio_class%write(filename, wrk2(:,:,:,m))
    end do
    deallocate(wrk2)


    ! all done, cleanup
    call timer_stop(t_total)
    call timer_print()
    call letkf_mpi_end()

  end subroutine letkf_driver_run





  !============================================================
  subroutine letkf_do_letkf()
    !TODO move this to another module
    integer :: ij
    integer, parameter :: maxpt = 100000
    integer :: rpoints(maxpt)
    real :: rdistance(maxpt)
    integer :: rnum, ob_cnt

    real :: hdxb(maxpt,mem), rdiag(maxpt), rloc(maxpt), dep(maxpt)
    real :: trans(mem,mem)
    integer :: timer1, timer2, n, timer3
    integer :: i

    real, allocatable:: wrk1(:,:), wrk2(:,:)
    real :: loc_h


    allocate(wrk1( grid_ns, mem))
    allocate(wrk2( grid_ns, mem))

    timer1 = timer_init("obs search")
    timer2 = timer_init("letkf_core_solve")
    timer3 = timer_init("letkf_core trans")

    ana_ij = 0

    do ij=1,ij_count
       ! search for all observations in a given radius of this gridpoint
       call timer_start(timer1)
       call letkf_obs_get(lat_ij(ij), lon_ij(ij), 1000.0e3, rpoints, rdistance, rnum)
       call timer_stop(timer1)

       ! if there are observations found, process them
       ob_cnt = 0
       do i=1,rnum
          n = rpoints(i)

          !TODO, do this earlier in program
          ! get rid of obs with bad QC values
          if (obs_qc(n) /= 0) cycle

          ! calculate observation localization, and
          ! get rid of obs outside of localization radius
          loc_h = loc_gc(rdistance(i), 100.0e3)
          if (loc_h <= 0 ) cycle

          ! use this observation
          ob_cnt = ob_cnt + 1
          hdxb(ob_cnt,:) = obs_ohx(:,n)  !TODO: should hdxb be transposed for efficiency?
          rdiag(ob_cnt)  = obs_list(n)%err
          rloc(ob_cnt) = loc_h
          dep(ob_cnt) = obs_list(n)%val - obs_ohx_mean(n)
       end do


       ! if there are still good quality observations to assimilate, do so
       if (ob_cnt > 0) then

          ! main LETKF equations
          call timer_start(timer2)
          call letkf_core_solve(&
               hdxb(:ob_cnt,:), rdiag(:ob_cnt), rloc(:ob_cnt),&
               dep(:ob_cnt), 1.0e0, trans)
          call timer_stop(timer2)

          ! calculate the ensemble increments
          call timer_start(timer3)
          !!@todo is there anyway to call sgemm without needing to copy in/out of
          !! the wrk1/wrk2 arrays (anal/gues need to be copied anyway, even if by compiler
          !! because they are being sliced from the front, somewhat faster on gfortran if
          !! i copy myself in a preallocated array)
          !! this might not be necessary with better intel compiler, which should do everythinh on the stack
          wrk1(:,:) = bkg_ij(ij,:,:)
          call sgemm('n','n', grid_ns, mem, mem, &
               1.0e0, wrk1(:,:), grid_ns, &
               trans, mem, 0.0e0, wrk2(:,:), grid_ns)
          ana_ij(ij,:,:) = wrk2(:,:)

          call timer_stop(timer3)
       else
          ana_ij(ij,:,:) = bkg_ij(ij,:,:)
       end if

    end do

    ! add the mean back to the analysis
    do i=1,mem
       ana_ij(:,:,i) = ana_ij(:,:,i) + bkg_mean_ij
    end do
    ana_mean_ij = sum(ana_ij,3) / mem

    !calculate spread
    ana_sprd_ij = 0
    do i=1,mem
       ana_sprd_ij = ana_sprd_ij + &
            (ana_ij(:,:,i) - ana_mean_ij) * &
            (ana_ij(:,:,i) - ana_mean_ij)
    end do
    ana_sprd_ij = sqrt(ana_sprd_ij/mem)

    ! done, cleanup
    deallocate(wrk1, wrk2)

  end subroutine letkf_do_letkf



end module letkf
