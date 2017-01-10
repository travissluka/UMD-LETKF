module letkf
  use timing
  use mpi
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


  !! The local variables
  real, allocatable :: lon_ij(:)
  real, allocatable :: lat_ij(:)

  real, allocatable :: gues_ij(:,:,:)
  real, allocatable :: gues_mean_ij(:,:)
  real, allocatable :: gues_sprd_ij(:,:)

  real, allocatable :: anal_ij(:,:,:)
  real, allocatable :: anal_mean_ij(:,:)
  real, allocatable :: anal_sprd_ij(:,:)

  logical :: initialized = .false.
  integer :: t_total




contains



  subroutine letkf_driver_init()
    call letkf_mpi_init()
    if (isroot) then
       print "(A)", "============================================================"
       print "(A)", " Universal Multi-Domain Local Ensemble Transform Kalman Filter"
       print "(A)", " (UMD-LETKF)"
       print "(A)", " version 0.0.0"
       print "(A)", " Travis Sluka (tsluka@umd.edu, travis.sluka@noaa.gov)"
       print "(A)", "============================================================"
       print "(A)", ""
       print "(A)", "*** beginning LETKF pre-initialization ***"
       print "(A)", ""
    end if
  end subroutine letkf_driver_init




  ! subroutine letkf_set_obsio(obsio_class_i)
  !   class(obsio), target, intent(in) :: obsio_class_i

  !   if(isroot) print *, 'Setting obs I/O class to: "','','"'
  !   obsio_class => obsio_class_i

  !   if(isroot) print *, 'Initializing obs I/O class...'
  !   call obsio_class%init()
  !   if(isroot) print *, ''
  ! end subroutine letkf_set_obsio



  ! subroutine letkf_set_stateio(stateio_class_i)
  !   class(stateio), target, intent(in) :: stateio_class_i
  !   stateio_class => stateio_class_i
  !   call stateio_class%init()
  !   if(isroot) print *, 'Setting state I/O class to: "',&
  !        trim(stateio_class%description),'"'
  !   if(isroot) print *, ""
  ! end subroutine letkf_set_stateio



  subroutine letkf_driver_run()
    integer ::  t_letkf, ierr

    integer :: i

    real, allocatable :: anal_mean(:,:,:)
    real, allocatable :: anal_sprd(:,:,:)
    integer :: t_init
    namelist /letkf_settings/ mem, obsqc_maxstd



    initialized = .true.
    t_total = timer_init("Total Runtime")
    call timer_start(t_total)

    ! Initialize
    t_init = timer_init("Initialization")
    call timer_start(t_init)

    !TODO
    ! ensure required preinitialized components are set

    ! read in namelist
    open(90, file=nml_filename)
    read(90, nml=letkf_settings)
    close(90)

    if (isroot) then
       print letkf_settings
    end if

    call letkf_state_init()
    call letkf_mpi_init2(mem, grid_nx*grid_ny)
    call letkf_obs_init("obsdef.cfg", "platdef.cfg")
    call letkf_core_init(mem)
    call timer_stop(t_init)


    ! pre-run setup
    ! ------------------------------

    ! read in observation departure ensemble
    call letkf_obs_read(obsio_class)

    ! read in background ensemble
    call letkf_read_gues()

    ! ensure output directory is setup right
    call system('mkdir -p OUTPUT')


    ! run LETKF
    t_letkf = timer_init("loop", TIMER_SYNC)
    call timer_start(t_letkf)
    call letkf_do_letkf()
    call timer_stop(t_letkf)
    call mpi_barrier(mpi_comm_letkf, ierr)

    ! gather the mean/sprd and write out
    if (isroot) then
       allocate(anal_mean(grid_nx, grid_ny, grid_ns))
       allocate(anal_sprd(grid_nx, grid_ny, grid_ns))
    end if
    do i=1,grid_ns
       call mpi_gatherv(anal_mean_ij(:,i), ij_count, mpi_real, &
            anal_mean(:,:,i), scatterv_count, scatterv_displ, mpi_real,&
            pe_root, mpi_comm_letkf, ierr)
       call mpi_gatherv(anal_sprd_ij(:,i), ij_count, mpi_real, &
            anal_sprd(:,:,i), scatterv_count, scatterv_displ, mpi_real,&
            pe_root, mpi_comm_letkf, ierr)
    end do
    if (isroot) then
       call stateio_class%write('OUTPUT/anal_mean.nc', anal_mean)
       call stateio_class%write('OUTPUT/anal_sprd.nc', anal_sprd)
       deallocate(anal_mean)
       deallocate(anal_sprd)
    end if


    call letkf_write_ana()
    ! all done

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

    if (isroot) then
       print *, ""
       print *, ""
       print *, "Main LETKF loop"
       print *, "============================================================"
    end if

    timer1 = timer_init("obs search")
    timer2 = timer_init("letkf_core_solve")
    timer3 = timer_init("letkf_core trans")

    anal_ij = 0

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
          wrk1(:,:) = gues_ij(ij,:,:)
          call sgemm('n','n', grid_ns, mem, mem, &
               1.0e0, wrk1(:,:), grid_ns, &
               trans, mem, 0.0e0, wrk2(:,:), grid_ns)
          anal_ij(ij,:,:) = wrk2(:,:)

          call timer_stop(timer3)
       else
          anal_ij(ij,:,:) = gues_ij(ij,:,:)
       end if

    end do

    ! add the mean back to the analysis
    do i=1,mem
       anal_ij(:,:,i) = anal_ij(:,:,i) + gues_mean_ij
    end do
    anal_mean_ij = sum(anal_ij,3) / mem

    !calculate spread
    anal_sprd_ij = 0
    do i=1,mem
       anal_sprd_ij = anal_sprd_ij + &
            (anal_ij(:,:,i) - anal_mean_ij) * &
            (anal_ij(:,:,i) - anal_mean_ij)
    end do
    anal_sprd_ij = sqrt(anal_sprd_ij/mem)


    deallocate(wrk1, wrk2)

    ! print *,"done",pe_rank


  end subroutine letkf_do_letkf





  !============================================================
  subroutine letkf_write_ana()
    real, allocatable :: ana(:,:,:,:)
    character(len=1024) :: filename
    integer :: i, m, ierr
    allocate(ana(grid_nx, grid_ny, grid_ns, size(ens_list)))

    if(isroot) then
       print *, "writing analysis ensemble members..."
    end if

    ! gather the ensemble members
    do m=1,mem
       !distribute variables
       do i=1,grid_ns
          call mpi_gatherv(&
               anal_ij(:,i,m), ij_count, mpi_real, &
               ana(:,:,i,ens_idx(m)), scatterv_count, scatterv_displ, mpi_real, &
               ens_map(m), mpi_comm_letkf, ierr)
       end do
    end do

    ! write out the members this proc is responsible for
    do m=1,size(ens_list)
       write (filename, '(A,I0.4,A)') 'OUTPUT/',ens_list(m),'.nc'
       call stateio_class%write(filename, ana(:,:,:,m))
    end do

    deallocate(ana)

  end subroutine letkf_write_ana



  !============================================================
  subroutine letkf_read_gues()
    !! @todo move this out of this module
    integer :: timer, timer2, timer3
    integer :: m, i
    character(len=1024) :: filename

    integer :: ierr, revcount

    real, allocatable :: gues(:,:,:,:)
    real, allocatable :: gues_mean(:,:,:)
    real, allocatable :: gues_sprd(:,:,:)

    timer = timer_init("read bg", TIMER_SYNC)
    timer2 = timer_init("read bg I/O", TIMER_SYNC)
    timer3 = timer_init("read bg MPI", TIMER_SYNC)
    call timer_start(timer)

    if (isroot) then
       print *, ""
       print *, ""
       print *, "Ensemble background"
       print *, "============================================================"
       print *,"Reading ensemble background ..."
       print *, "slabs = ", grid_ns
       print "(A,I6,A,I6,A,I6,A)","  shape = (",grid_nx," x ",grid_ny," x ",grid_ns,")"
    end if

    allocate(gues(grid_nx, grid_ny, grid_ns, size(ens_list)))
    allocate(gues_ij(ij_count, grid_ns, mem))
    allocate(gues_mean_ij(ij_count, grid_ns))
    allocate(gues_sprd_ij(ij_count, grid_ns))
    allocate(lon_ij(ij_count))
    allocate(lat_ij(ij_count))

    allocate(anal_ij(ij_count, grid_ns, mem))
    allocate(anal_mean_ij(ij_count,grid_ns))
    allocate(anal_sprd_ij(ij_count,grid_ns))


    ! for each ensemble member we are responsible for loading, read in
    ! the ensemble background file
    call timer_start(timer2)
    !! @todo un-hardcode this
    do m=1,size(ens_list)
       write (filename, '(A,I0.4,A)') 'INPUT/gues/',ens_list(m),'.nc'
       call stateio_class%read(filename, gues(:,:,:,m))
    end do
    call timer_stop(timer2)


    !distribute segments of each ensemble member to the appropriate proc
    call timer_start(timer3)
    revcount=ij_count
    if (isroot) print*, "distributing via MPI..."
    do m=1,mem
       !distribute variables
       do i=1,grid_ns
          call mpi_scatterv(&
               gues(:,:,i,ens_idx(m)), scatterv_count, scatterv_displ, mpi_real, &
               gues_ij(:,i,m), revcount, mpi_real, &
               ens_map(m), mpi_comm_letkf, ierr)
       end do
    end do
    call timer_stop(timer3)


    !determine lat/lon allocation
    if (isroot) print *,"distributing lat/lon grid..."
!    call getlatlon()
    call mpi_scatterv(&
         lon, scatterv_count, scatterv_displ, mpi_real, &
         lon_ij, revcount, mpi_real, &
         pe_root, mpi_comm_letkf, ierr)
    call mpi_scatterv(&
         lat, scatterv_count, scatterv_displ, mpi_real, &
         lat_ij, revcount, mpi_real, &
         pe_root, mpi_comm_letkf, ierr)


    ! calculate the mean and spread of the background
    ! after this, guesxd_ij will contain perturbations only
    if(isroot) print *, "calculating background mean/spread..."
    gues_mean_ij  = sum(gues_ij, 3) / mem
    gues_sprd_ij = 0
    do m=1,mem
       gues_ij(:,:,m) = gues_ij(:,:,m) - gues_mean_ij
       gues_sprd_ij   = gues_sprd_ij   + gues_ij(:,:,m)*gues_ij(:,:,m)
    end do
    gues_sprd_ij = sqrt(gues_sprd_ij/mem)


    ! collect and save the combined mean/spread
    ! todo, write separate mean/spread calculating function
    if (isroot) then
       allocate(gues_mean(grid_nx, grid_ny, grid_ns))
       allocate(gues_sprd(grid_nx, grid_ny, grid_ns))
    end if
    do i=1,grid_ns
       call mpi_gatherv(gues_mean_ij(:,i), ij_count, mpi_real, &
            gues_mean(:,:,i), scatterv_count, scatterv_displ, mpi_real,&
            pe_root, mpi_comm_letkf, ierr)
       call mpi_gatherv(gues_sprd_ij(:,i), ij_count, mpi_real, &
            gues_sprd(:,:,i), scatterv_count, scatterv_displ, mpi_real,&
            pe_root, mpi_comm_letkf, ierr)
    end do
    if (isroot) then
       call stateio_class%write('OUTPUT/gues_mean.nc', gues_mean)
       call stateio_class%write('OUTPUT/gues_sprd.nc', gues_sprd)

       !done cleanup
       deallocate(gues_mean)
       deallocate(gues_sprd)
    end if


    !cleanup
    deallocate(gues)

    call timer_stop(timer)
  end subroutine letkf_read_gues


  ! subroutine letkf_set_grid(lon_grd_in, lat_grd_in)
  !   real, allocatable, intent(in) :: lon_grd_in(:,:)
  !   real, allocatable, intent(in) :: lat_grd_in(:,:)

  !   !! TODO, do some error checking on these
  !   lat = lat_grd_in
  !   lon = lon_grd_in
  ! end subroutine letkf_set_grid

end module letkf
