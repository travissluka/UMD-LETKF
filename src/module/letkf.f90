module letkf
  use timing
  use mpi
  use letkf_common
  use letkf_obs
  use letkf_mpi
  use letkf_core
  use letkf_state
  use letkf_loc
  use kdtree

  implicit none


  !TODO: move these elsewhere
  real, allocatable :: lat(:,:), lon(:,:)
  real, allocatable :: gues3d_ij(:,:,:,:)
  real, allocatable :: gues2d_ij(:,:,:)
  real, allocatable :: anal3d_ij(:,:,:,:)
  real, allocatable :: anal2d_ij(:,:,:)

  real, allocatable :: gues3d_mean_ij(:,:,:)
  real, allocatable :: gues3d_sprd_ij(:,:,:)
  real, allocatable :: gues2d_mean_ij(:,:)
  real, allocatable :: gues2d_sprd_ij(:,:)
  real, allocatable :: lon_ij(:)
  real, allocatable :: lat_ij(:)

  real, allocatable :: anal3d_mean_ij(:,:,:)
  real, allocatable :: anal2d_mean_ij(:,:)
  real, allocatable :: anal3d_sprd_ij(:,:,:)
  real, allocatable :: anal2d_sprd_ij(:,:)

  logical :: initialized = .false.
  integer :: t_total
  class(obsio), pointer :: usr_obsio_class
contains



  !############################################################
  subroutine letkf_driver_init()
    integer :: t_init
    namelist /letkf_settings/ mem, obsqc_maxstd

    if (initialized) then
       print *, "ERROR: letkf_driver_init() called more than once"
       stop 1
    end if
    initialized = .true.
    t_total = timer_init("Total Runtime")
    call timer_start(t_total)

    ! Initialize
    t_init = timer_init("Initialization")
    call timer_start(t_init)

    open(90, file=nml_filename)
    read(90, nml=letkf_settings)
    close(90)

    call letkf_mpi_init()
    if (pe_isroot) then
       print "(A)", "============================================================"
       print "(A)", " Universal Multi-Domain Local Ensemble Transform Kalman Filter"
       print "(A)", " (UMD-LETKF)"
       print "(A)", " version 0.0.0"
       print "(A)", " Travis Sluka (tsluka@umd.edu, travis.sluka@noaa.gov)"
       print "(A)", "============================================================"
       print "(A)", ""
       print letkf_settings
    end if

    call letkf_state_init()
    call letkf_mpi_init2(mem, grid_x*grid_y)

    call letkf_read_config
    call letkf_core_init(mem)
    call timer_stop(t_init)

  end subroutine letkf_driver_init



  subroutine letkf_set_obsio(obsio_class)
    class(obsio), target, intent(in) :: obsio_class
    call obsio_class%init()
    usr_obsio_class => obsio_class
  end subroutine letkf_set_obsio



  subroutine letkf_set_stateio(stateio_class)
    class(stateio), intent(in) :: stateio_class
  end subroutine letkf_set_stateio


  !############################################################
  subroutine letkf_driver_run()
    integer ::  t_letkf, ierr

    integer :: i, k
    real, allocatable :: anal3d_mean(:,:,:,:)
    real, allocatable :: anal2d_mean(:,:,:)
    real, allocatable :: anal3d_sprd(:,:,:,:)
    real, allocatable :: anal2d_sprd(:,:,:)





    call letkf_obs_read(usr_obsio_class)
    call letkf_read_gues()

    call system('mkdir -p OUTPUT')

    ! run LETKF
    t_letkf = timer_init("loop", TIMER_SYNC)
    call timer_start(t_letkf)
    call letkf_do_letkf()
    call timer_stop(t_letkf)
    call mpi_barrier(mpi_comm_letkf, ierr)

    ! gather the mean/sprd and write out
    if (pe_isroot) then
       allocate(anal3d_mean(grid_x, grid_y, grid_z, grid_3d))
       allocate(anal2d_mean(grid_x, grid_y, grid_2d))
       allocate(anal3d_sprd(grid_x, grid_y, grid_z, grid_3d))
       allocate(anal2d_sprd(grid_x, grid_y, grid_2d))
    end if
    do i=1,grid_3d
       do k=1,grid_z
          call mpi_gatherv(anal3d_mean_ij(:,k,i), ij_count, mpi_real, &
               anal3d_mean(:,:,k,i), scatterv_count, scatterv_displ, mpi_real,&
               pe_root, mpi_comm_letkf, ierr)
          call mpi_gatherv(anal3d_sprd_ij(:,k,i), ij_count, mpi_real, &
               anal3d_sprd(:,:,k,i), scatterv_count, scatterv_displ, mpi_real,&
               pe_root, mpi_comm_letkf, ierr)
       end do
    end do
    do i=1,grid_2d
       call mpi_gatherv(anal2d_mean_ij(:,i), ij_count, mpi_real, &
            anal2d_mean(:,:,i), scatterv_count, scatterv_displ, mpi_real,&
            pe_root, mpi_comm_letkf, ierr)
       call mpi_gatherv(anal2d_sprd_ij(:,i), ij_count, mpi_real, &
            anal2d_sprd(:,:,i), scatterv_count, scatterv_displ, mpi_real,&
            pe_root, mpi_comm_letkf, ierr)
    end do
    if (pe_isroot) then
       call letkf_state_write('OUTPUT/anal_mean.nc', anal3d_mean, anal2d_mean)
       call letkf_state_write('OUTPUT/anal_sprd.nc', anal3d_sprd, anal2d_sprd)
       deallocate(anal3d_mean)
       deallocate(anal2d_mean)
       deallocate(anal3d_sprd)
       deallocate(anal2d_sprd)
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
    integer,parameter :: maxpt = 100000
    integer :: rpoints(maxpt)
    real :: rdistance(maxpt)
    integer :: rnum, ob_cnt

    real :: hdxb(maxpt,mem), rdiag(maxpt), rloc(maxpt), dep(maxpt)
    real :: trans(mem,mem)
    integer :: timer1, timer2, n, timer3
    integer :: i

    real, allocatable:: wrk3d1(:,:,:), wrk3d2(:,:,:)
    real, allocatable:: wrk2d1(:,:), wrk2d2(:,:)
    real :: loc_h


    allocate(wrk3d1( grid_z, grid_3d, mem))
    allocate(wrk3d2( grid_z, grid_3d, mem))
    allocate(wrk2d1( grid_2d, mem))
    allocate(wrk2d2( grid_2d, mem))

    if (pe_isroot) then
       print *, ""
       print *, ""
       print *, "Main LETKF loop"
       print *, "============================================================"
    end if

    timer1 = timer_init("obs search")
    timer2 = timer_init("letkf_core_solve")
    timer3 = timer_init("letkf_core trans")

    anal3d_ij = 0
    anal2d_ij = 0

    do ij=1,ij_count

       ! search for all observations in a given radius of this gridpoint
       call timer_start(timer1)
       call kd_search_radius(obs_tree, &
            (/lon_ij(ij)*1.0e0, lat_ij(ij)*1.0e0 /), 1000.0e3, &
            rpoints, rdistance, rnum, .false.)
       call timer_stop(timer1)

       ! if there are observations found, process them
       ob_cnt = 0
       do i=1,rnum
          n = rpoints(i)

          ! get rid of obs with bad QC values
          if (obs_qc(n) /= 0) cycle

          ! calculate observation localization, and
          ! get rid of obs outside of localization radius
          loc_h = loc_gc(rdistance(i), 500.0e3)
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
          wrk3d1 = gues3d_ij(ij,:,:,:)
          call sgemm('n','n', grid_3d*grid_z, mem, mem, &
               1.0e0, wrk3d1, grid_3d*grid_z, &
               trans, mem, 0.0e0, wrk3d2, grid_3d*grid_z)
          anal3d_ij(ij,:,:,:) = wrk3d2

          wrk2d1(:,:) = gues2d_ij(ij,:,:)
          call sgemm('n','n', grid_2d, mem, mem, &
               1.0e0, wrk2d1(:,:), grid_2d, &
               trans, mem, 0.0e0, wrk2d2(:,:), grid_2d)
          anal2d_ij(ij,:,:) = wrk2d2(:,:)

          call timer_stop(timer3)
       else
          anal3d_ij(ij,:,:,:) = gues3d_ij(ij,:,:,:)
          anal2d_ij(ij,:,:) = gues2d_ij(ij,:,:)
       end if

    end do

    ! add the mean back to the analysis
    do i=1,mem
       anal3d_ij(:,:,:,i) = anal3d_ij(:,:,:,i) + gues3d_mean_ij
       anal2d_ij(:,:,i) = anal2d_ij(:,:,i) + gues2d_mean_ij
    end do
    anal3d_mean_ij = sum(anal3d_ij,4) / mem
    anal2d_mean_ij = sum(anal2d_ij,3) / mem

    !calculate spread
    anal3d_sprd_ij = 0
    anal2d_sprd_ij = 0
    do i=1,mem
       anal3d_sprd_ij = anal3d_sprd_ij + &
                        (anal3d_ij(:,:,:,i) - anal3d_mean_ij) * &
                        (anal3d_ij(:,:,:,i) - anal3d_mean_ij)
       anal2d_sprd_ij = anal2d_sprd_ij + &
                        (anal2d_ij(:,:,i) - anal2d_mean_ij) * &
                        (anal2d_ij(:,:,i) - anal2d_mean_ij)
    end do
    anal3d_sprd_ij = sqrt(anal3d_sprd_ij/mem)
    anal2d_sprd_ij = sqrt(anal2d_sprd_ij/mem)


    deallocate(wrk3d1, wrk3d2)
    deallocate(wrk2d1, wrk2d2)

    print *,"done",pe_rank


  end subroutine letkf_do_letkf


  !############################################################
  subroutine letkf_read_config()
    integer :: timer

    timer = timer_init("read config", TIMER_SYNC)
    call timer_start(timer)


    call letkf_obs_init("obsdef.cfg", "platdef.cfg")
    call timer_stop(timer)

  end subroutine letkf_read_config



  !############################################################
!  subroutine letkf_read_obs()
!  end subroutine letkf_read_obs



  !============================================================
  subroutine letkf_write_ana()
    real, allocatable :: ana3d(:,:,:,:,:)
    real, allocatable :: ana2d(:,:,:,:)
    character(len=1024) :: filename
    integer :: i,k,m, ierr
    allocate(ana3d(grid_x, grid_y, grid_z, grid_3d, size(ens_list)))
    allocate(ana2d(grid_x, grid_y, grid_2d, size(ens_list)))

    if(pe_isroot) then
       print *, "writing analysis ensemble members..."
    end if

    ! gather the ensemble members
    do m=1,mem
       ! distribute 3d variables
       do i=1,grid_3d
          do k=1,grid_z
             call mpi_gatherv(&
                  anal3d_ij(:,k,i,m), ij_count, mpi_real, &
                  ana3d(:,:,k,i,ens_idx(m)), scatterv_count, scatterv_displ, mpi_real, &
                  ens_map(m), mpi_comm_letkf, ierr)
          end do
       end do

       !distribute 2d variables
       do i=1,grid_2d
          call mpi_gatherv(&
               anal2d_ij(:,i,m), ij_count, mpi_real, &
               ana2d(:,:,i,ens_idx(m)), scatterv_count, scatterv_displ, mpi_real, &
               ens_map(m), mpi_comm_letkf, ierr)
       end do
    end do

    ! write out the members this proc is responsible for
    do m=1,size(ens_list)
       write (filename, '(A,I0.4,A)') 'OUTPUT/',ens_list(m),'.nc'
       call letkf_state_write(filename, ana3d(:,:,:,:,m), ana2d(:,:,:,m))
    end do

    deallocate(ana3d)
    deallocate(ana2d)

  end subroutine letkf_write_ana



  !============================================================
  subroutine letkf_read_gues()
    !! @todo move this out of this module
    integer :: timer, timer2, timer3
    integer :: m, i, k
    character(len=1024) :: filename

    integer :: ierr, revcount

    !! @TODO, can i collapse these down into a single variable?
    !! have 2d just be a single level 3d?
    real, allocatable :: gues3d(:,:,:,:,:)
    real, allocatable :: gues2d(:,:,:,:)


    real, allocatable :: gues3d_mean(:,:,:,:)
    real, allocatable :: gues2d_mean(:,:,:)
    real, allocatable :: gues3d_sprd(:,:,:,:)
    real, allocatable :: gues2d_sprd(:,:,:)

    timer = timer_init("read bg", TIMER_SYNC)
    timer2 = timer_init("read bg I/O", TIMER_SYNC)
    timer3 = timer_init("read bg MPI", TIMER_SYNC)
    call timer_start(timer)

    if (pe_isroot) then
       print *, ""
       print *, ""
       print *, "Ensemble background"
       print *, "============================================================"
       print *,"Reading ensemble background ..."
       print *,"3D vars = ", grid_3d
       print *,"2D vars = ", grid_2d
       print "(A,I6,A,I6,A,I6,A)","  shape = (",grid_x," x ",grid_y," x ",grid_z,")"
    end if

    allocate(gues3d(grid_x, grid_y, grid_z, grid_3d, size(ens_list)))
    allocate(gues2d(grid_x, grid_y, grid_2d, size(ens_list)))

    allocate(gues3d_ij(ij_count, grid_z, grid_3d, mem))
    allocate(gues3d_mean_ij(ij_count, grid_z, grid_3d))
    allocate(gues3d_sprd_ij(ij_count, grid_z, grid_3d))
    allocate(gues2d_ij(ij_count, grid_2d, mem))
    allocate(gues2d_mean_ij(ij_count, grid_2d))
    allocate(gues2d_sprd_ij(ij_count, grid_2d))
    allocate(lon_ij(ij_count))
    allocate(lat_ij(ij_count))

    allocate(anal3d_ij(ij_count, grid_z, grid_3d, mem))
    allocate(anal2d_ij(ij_count, grid_2d, mem))
    allocate(anal3d_mean_ij(ij_count, grid_z, grid_3d))
    allocate(anal2d_mean_ij(ij_count,grid_3d))
    allocate(anal3d_sprd_ij(ij_count, grid_z, grid_3d))
    allocate(anal2d_sprd_ij(ij_count,grid_2d))


    ! for each ensemble member we are responsible for loading, read in
    ! the ensemble background file
    call timer_start(timer2)
    !! @todo un-hardcode this
    do m=1,size(ens_list)
       write (filename, '(A,I0.4,A)') 'INPUT/gues/',ens_list(m),'.nc'
       call letkf_state_read(filename, gues3d(:,:,:,:,m), gues2d(:,:,:,m))
    end do
    call timer_stop(timer2)


    !distribute segments of each ensemble member to the appropriate proc
    call timer_start(timer3)
    revcount=ij_count
    if (pe_isroot) print*, "distributing via MPI..."
    do m=1,mem
       ! distribute 3d variables
       do i=1,grid_3d
          do k=1,grid_z
             call mpi_scatterv(&
                  gues3d(:,:,k,i,ens_idx(m)), scatterv_count, scatterv_displ, mpi_real, &
                  gues3d_ij(:,k,i,m), revcount, mpi_real,&
                  ens_map(m), mpi_comm_letkf, ierr)
          end do
       end do

       !distribute 2d variables
       do i=1,grid_2d
          call mpi_scatterv(&
               gues2d(:,:,i,ens_idx(m)), scatterv_count, scatterv_displ, mpi_real, &
               gues2d_ij(:,i,m), revcount, mpi_real, &
               ens_map(m), mpi_comm_letkf, ierr)
       end do
    end do
    call timer_stop(timer3)


    !determine lat/lon allocation
    if (pe_isroot) print *,"distributing lat/lon grid..."
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
    if(pe_isroot) print *, "calculating background mean/spread..."
    gues3d_mean_ij  = sum(gues3d_ij, 4) / mem
    gues3d_sprd_ij = 0
    gues2d_mean_ij  = sum(gues2d_ij, 3) / mem
    gues2d_sprd_ij = 0
    do m=1,mem
       gues3d_ij(:,:,:,m) = gues3d_ij(:,:,:,m) - gues3d_mean_ij
       gues2d_ij(:,:,m) = gues2d_ij(:,:,m) - gues2d_mean_ij

       gues3d_sprd_ij = gues3d_sprd_ij + gues3d_ij(:,:,:,m)*gues3d_ij(:,:,:,m)
       gues2d_sprd_ij = gues2d_sprd_ij + gues2d_ij(:,:,m)*gues2d_ij(:,:,m)
    end do
    gues3d_sprd_ij = sqrt(gues3d_sprd_ij/mem)
    gues2d_sprd_ij = sqrt(gues2d_sprd_ij/mem)


    ! collect and save the combined mean/spread
    ! todo, write separate mean/spread calculating function
    if (pe_isroot) then
       allocate(gues3d_mean(grid_x, grid_y, grid_z, grid_3d))
       allocate(gues2d_mean(grid_x, grid_y, grid_2d))
       allocate(gues3d_sprd(grid_x, grid_y, grid_z, grid_3d))
       allocate(gues2d_sprd(grid_x, grid_y, grid_2d))
    end if
    do i=1,grid_3d
       do k=1,grid_z
          call mpi_gatherv(gues3d_mean_ij(:,k,i), ij_count, mpi_real, &
               gues3d_mean(:,:,k,i), scatterv_count, scatterv_displ, mpi_real,&
               pe_root, mpi_comm_letkf, ierr)
          call mpi_gatherv(gues3d_sprd_ij(:,k,i), ij_count, mpi_real, &
               gues3d_sprd(:,:,k,i), scatterv_count, scatterv_displ, mpi_real,&
               pe_root, mpi_comm_letkf, ierr)
       end do
    end do
    do i=1,grid_2d
       call mpi_gatherv(gues2d_mean_ij(:,i), ij_count, mpi_real, &
            gues2d_mean(:,:,i), scatterv_count, scatterv_displ, mpi_real,&
            pe_root, mpi_comm_letkf, ierr)
       call mpi_gatherv(gues2d_sprd_ij(:,i), ij_count, mpi_real, &
            gues2d_sprd(:,:,i), scatterv_count, scatterv_displ, mpi_real,&
            pe_root, mpi_comm_letkf, ierr)
    end do
    if (pe_isroot) then
       call letkf_state_write('OUTPUT/gues_mean.nc', gues3d_mean, gues2d_mean)
       call letkf_state_write('OUTPUT/gues_sprd.nc', gues3d_sprd, gues2d_sprd)

       !done cleanup
       deallocate(gues3d_mean)
       deallocate(gues2d_mean)
       deallocate(gues3d_sprd)
       deallocate(gues2d_sprd)
    end if


    !cleanup
    deallocate(gues3d)
    deallocate(gues2d)

    call timer_stop(timer)
  end subroutine letkf_read_gues


  subroutine letkf_set_grid(lon_grd_in, lat_grd_in)
    real, allocatable, intent(in) :: lon_grd_in(:,:)
    real, allocatable, intent(in) :: lat_grd_in(:,:)

    !! TODO, do some error checking on these
    lat = lat_grd_in
    lon = lon_grd_in
  end subroutine letkf_set_grid

end module letkf