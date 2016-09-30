module letkf
  use timing
  use mpi
  use letkf_obs
  use letkf_mpi
  use global
  use letkf_obs_dat
  use netcdf
  use letkf_core
  
  implicit none

  type :: letkf_solver
   contains
     procedure :: run => letkf_solver_run
     procedure :: read_config => letkf_solver_read_config
     procedure :: read_obs => letkf_solver_read_obs
  end type letkf_solver


  !TODO: move these elsewhere
  real, allocatable :: lat(:,:), lon(:,:)
  real, allocatable :: gues3d_ij(:,:,:,:)
  real, allocatable :: gues2d_ij(:,:,:)
  real, allocatable :: lon_ij(:)
  real, allocatable :: lat_ij(:)

contains


  !############################################################
  subroutine letkf_solver_run(self)
    class(letkf_solver) :: self
    integer :: t_total, t_init, t_letkf, ierr


    namelist /letkf_settings/ mem, obsqc_maxstd
    
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
       print "(A)", "============================================================"
       print "(A)", " Universal Multi-Domain Local Ensemble Transform Kalman Filter"
       print "(A)", " (UMD-LETKF)"
       print "(A)", " version 0.0.0"
       print "(A)", " tsluka@umd.edu"
       print "(A)", "============================================================"
       print "(A)", ""
       print letkf_settings
    end if
    call letkf_mpi_init2(mem)
    
    call self%read_config
    call letkf_core_init(mem)
    call timer_stop(t_init)

    call self%read_obs
    call letkf_read_gues()


    
    ! run LETKF
    t_letkf = timer_init("core", TIMER_SYNC)
    call timer_start(t_letkf)
    call letkf_do_letkf()    
    call timer_stop(t_letkf)
    call mpi_barrier(mpi_comm_letkf, ierr)
    
    
    ! all done
    
    call timer_stop(t_total)

    call timer_print()
    
    call letkf_mpi_end()

    
  end subroutine letkf_solver_run



  subroutine letkf_do_letkf()
    !TODO move this to another module    
    integer :: ij
    integer,parameter :: maxpt = 100000
    integer :: rpoints(maxpt)
    real(dp) :: rdistance(maxpt)
    integer :: rnum

    real(dp) :: hdxb(maxpt,mem), rdiag(maxpt), rloc(maxpt), dep(maxpt)
    real(dp) :: trans(mem,mem)
    integer :: timer1, timer2, i, n, ierr

    if (pe_isroot) then
       print *, ""
       print *, ""       
       print *, "Main LETKF loop"
       print *, "============================================================"
    end if
    
    timer1 = timer_init("kd_search")
    timer2 = timer_init("letkf_core_solve")
    do ij=1,ij_count

       call timer_start(timer1)
       call kd_search_radius(obs_tree, &
            (/lon_ij(ij)*1.0d0, lat_ij(ij)*1.0d0 /), 3000.0d3, &
            rpoints, rdistance, rnum, .false.)
        ! if (rnum > 1000) then
        !    call kd_search_nnearest(obs_tree, &
        !         (/lon_ij(ij)*1.0d0, lat_ij(ij)*1.0d0 /), 500, &
        !         rpoints,rdistance,rnum, .false.)
        ! end if

       ! if(pe_isroot) then
       !    print *, ij, rnum
       ! end if
       call timer_stop(timer1)

       ! only if there are observations to assimilate
       if (rnum > 0) then
          do i=1,rnum
             n = rpoints(i)
             !TODO: should hdxb be transposed for efficiency?
             hdxb(i,:) = obs_ohx(:,n)
             rdiag(i)  = obs_list(n)%err
             rloc(i) = 1.0
             dep(i) = obs_ohx_mean(n) - obs_list(n)%val
          end do
       
          call timer_start(timer2)
!          if  (pe_rank == 1) then
!             print *, ij, rnum
!          end if
          call letkf_core_solve(&
               hdxb(:rnum,:), rdiag(:rnum), rloc(:rnum),&
               dep(:rnum), 1.0d0, trans)
          call timer_stop(timer2)
       end if     
    end do
    print *,"done",pe_rank


  end subroutine letkf_do_letkf

  
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



  !============================================================
  subroutine letkf_read_gues()
    !! @todo move this out of this module
    integer :: timer, timer2, timer3
    integer :: m, i, j, k, l
    character(len=1024) :: filename

    integer :: unit, nrec
    logical :: ex
    integer :: ierr, revcount
    
    !! @TODO, can i colapse these down into a single variable?
    !! have 2d just be a single level 3d?
    real, allocatable :: gues3d(:,:,:,:,:)
    real, allocatable :: gues2d(:,:,:,:)

    


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
    allocate(gues2d_ij(ij_count, grid_2d, mem))
    allocate(lon_ij(ij_count))
    allocate(lat_ij(ij_count))    
    

    ! for each ensemble member we are responsible for loading, read in
    ! the ensemble background file
    call timer_start(timer2)
    !! @todo un-hardcode this
    do m=1,size(ens_list)
       write (filename, '(A,I0.3,A)') 'data/gues/gues',ens_list(m),'.grd'

       ! check to make sure file exists
       inquire(file=filename, exist=ex)
       if (.not. ex) then
          print *, "ERROR: file does not exist ",filename
          stop 1
       end if

       !read in the file
       open(newunit=unit, file=filename, form='unformatted', access='direct', recl=grid_x*grid_y*4)
       nrec=1
       do l=1,grid_3d
          do k=1,grid_z
             read(unit,rec=nrec) gues3d(:,:,k,l,m)
             nrec = nrec + 1
          end do
       end do
       do l=1,grid_2d
          read(unit, rec=nrec) gues2d(:,:,l,m)
          nrec = nrec + 1
       end do      
       close(unit)
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
    call getlatlon()
    call mpi_scatterv(&
         lon, scatterv_count, scatterv_displ, mpi_real, &
         lon_ij, revcount, mpi_real, &
         pe_root, mpi_comm_letkf, ierr)
    call mpi_scatterv(&
         lat, scatterv_count, scatterv_displ, mpi_real, &
         lat_ij, revcount, mpi_real, &
         pe_root, mpi_comm_letkf, ierr)
    
    !cleanup
    deallocate(gues3d)
    deallocate(gues2d)
    
    call timer_stop(timer)
  end subroutine letkf_read_gues


  
  subroutine getlatlon()
    integer :: ncid, ierr, varid, n

    allocate(lat(grid_x, grid_y))
    allocate(lon(grid_x, grid_y))    
    
    ierr = nf90_open('data/grid_spec.nc', nf90_nowrite, ncid)
    ierr = nf90_inq_varid(ncid, 'xta', varid)
    ierr = nf90_get_var(ncid, varid, lon(:,1))
    ierr = nf90_inq_varid(ncid, 'yta', varid)
    ierr = nf90_get_var(ncid, varid, lat(1,:))
    ierr = nf90_close(ncid)

    do n=2,grid_x
       lat(n,:) = lat(1,:)
    end do
    do n=2, grid_y
       lon(:,n) = lon(:,1)
    end do
  end subroutine getlatlon
end module letkf
