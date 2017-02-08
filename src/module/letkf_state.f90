module letkf_state
  !! performs I/O for the background and analysis states
  !!
  ! library modules
  use timing
  use letkf_mpi
  use letkf_state_I
  use letkf_state_generic
  use letkf_obs
  use kdtree

  implicit none
  private

  public :: letkf_state_init

  ! public module variables
  ! ------------------------------------------------------------
  integer, public, protected :: grid_nz
    !! number of vertical levels for 3D variables
  real,    public, protected, allocatable :: lat_ij(:)
    !! ***Size is ( [[letkf_mpi:ij_count]] ) ***
  real,    public, protected, allocatable :: lon_ij(:)
    !! ***Size is ( [[letkf_mpi:ij_count]] ) ***
  real,    public, protected, allocatable :: mask_ij(:)

  real,    public, protected, allocatable :: bkg_ij(:,:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_mpi:grid_ns]], [[letkf_mpi:mem]] ) ***
  real,    public, protected, allocatable :: bkg_mean_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_mpi:grid_ns]] ) ***
  real,    public, protected, allocatable :: bkg_sprd_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_mpi:grid_ns]] ) ***

  real,    public,            allocatable :: ana_ij(:,:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_mpi:grid_ns]], [[letkf_mpi:mem]] ) ***
  real,    public,            allocatable :: ana_mean_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_mpi:grid_ns]] ) ***
  real,    public,            allocatable :: ana_sprd_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_mpi:grid_ns]] ) ***


  ! private module variables
  class(stateio), public, pointer :: stateio_class



contains



  !============================================================
  subroutine letkf_state_init(nml_filename)
    !! Initialize the state vector by reading in the background ensemble
    !! and distributing across cores

    character(len=*), intent(in) :: nml_filename
      !! filename of namelist to read in for the ***grid_def*** section

    real, allocatable :: gues(:,:,:,:)
    real, allocatable :: wrk(:,:,:)
    real, allocatable :: wrk2(:)
    real, allocatable :: lat(:,:), lon(:,:), mask(:,:)

    integer :: unit, tbgscatter, tbgread, tbgms
    character(len=1024) :: filename
    integer :: m, i, j

    ! lat/lon gridpoint kdtree, if needed for test observations
    real :: r
    integer :: x,y
    type(kd_root) :: llkd_root
    integer, allocatable :: llkd_x(:), llkd_y(:)
    real,    allocatable :: llkd_lons(:), llkd_lats(:)
    integer :: r_points(1), r_num
    real    :: r_distance(1)
    

    if(pe_isroot) then
       print *, new_line('a'), &
            new_line('a'), '============================================================',&
            new_line('a'), ' letkf_state_init() : ',&
            new_line('a'), '============================================================'
    end if

    ! read in our section of the namelist
!    open(newunit=unit, file=nml_filename)
!    read(unit, nml=grid_def)
!    close(unit)
!    if (pe_isroot) then
!       print grid_def
!    end if


    !create the state io class, and initialize
    allocate(stateio_generic :: stateio_class)
    call stateio_class%init(grid_nx, grid_ny, grid_nz)
    allocate(lat(grid_nx, grid_ny))
    allocate(lon(grid_nx, grid_ny))
    call stateio_class%latlon(lat,lon)
    allocate(mask(grid_nx, grid_ny))
    call stateio_class%mask(mask)


    ! read in the background members that our process is responsible for
     if (pe_isroot) then
        print *, ""
        print *, "Reading ensemble background..."
     end if


     ! console output synchronization
     call letkf_mpi_barrier(syncio=.true.)


     ! read in the files
     allocate(gues(grid_nx, grid_ny, grid_ns, size(ens_list)))
     tbgread = timer_init("    bkg_read", TIMER_SYNC)
     call timer_start(tbgread)
     do m=1,size(ens_list)
        write (filename, '(A,I0.4,A)') 'INPUT/gues/',ens_list(m)
        print '(A,I5,3A)', " PROC ",pe_rank," is READING file: ",&
             trim(filename),trim(stateio_class%extension)
        call stateio_class%read(filename, gues(:,:,:,m))
     end do
     call timer_stop(tbgread)


     ! console output synchronization
     call letkf_mpi_barrier(syncio=.true.)


     ! before scattering, see if the obs module needs the full
     ! background in order to generate test observations
     if (obs_test) then
        if(pe_isroot) then
           print *,""
           print *, "Generating test observations from background..."
           print *, size(obs_list), "observations"
        end if

        ! generate the kd tree for obs location lookup
        ! TODO: filter out masked points
        allocate(llkd_lons(grid_nx*grid_ny))
        allocate(llkd_lats(grid_nx*grid_ny))
        allocate(llkd_x(grid_nx*grid_ny))
        allocate(llkd_y(grid_nx*grid_ny))
        do x=1,grid_nx
           do y=1,grid_ny
              llkd_lons((y-1)*grid_nx + x) = lon(x,y)
              llkd_lats((y-1)*grid_nx + x) = lat(x,y)
              llkd_x((y-1)*grid_nx + x) = x
              llkd_y((y-1)*grid_nx + x) = y
           end do
        end do
        call kd_init(llkd_root, llkd_lons, llkd_lats)

        ! for each test ob, get the closest grid point and use that as its value
        !TODO use the state config to get the right slab
        do i =1,size(obs_list)
           call kd_search_nnearest(llkd_root, obs_list(i)%lon, obs_list(i)%lat,&
                1, r_points, r_distance, r_num, .false.)
           do m=1,size(ens_list)
              obs_ohx(ens_list(m), i) = gues(&
                   llkd_x(r_points(1)), llkd_y(r_points(1)), &
                   merge(1,41,obs_list(i)%id==2210), m)
           end do
        end do
        call letkf_mpi_obs(obs_ohx)        

        ! update the observation ensemble mean/departures using
        ! the desired increment currently stored in obs_list%val
        obs_ohx_mean = sum(obs_ohx,1)/mem
        do i = 1,size(obs_list)
           obs_ohx(:,i) = obs_ohx(:,i) - obs_ohx_mean(i)
           r =  obs_list(i)%val
           obs_list(i)%val= obs_ohx_mean(i)
           obs_ohx_mean(i) = obs_ohx_mean(i) + r
        end do
        
        deallocate(llkd_lons)
        deallocate(llkd_lats)
        deallocate(llkd_x)
        deallocate(llkd_y)
     end if
        

     ! scatter grids to mpi procs
     if (pe_isroot) then
        print *,""
        print *, "Scattering across procs..."
     end if
     allocate(bkg_ij(mem, grid_ns,ij_count))
     allocate(bkg_mean_ij(grid_ns, ij_count))
     allocate(bkg_sprd_ij(grid_ns, ij_count))
     allocate(ana_ij(mem, grid_ns, ij_count))
     allocate(ana_mean_ij(grid_ns, ij_count))
     allocate(ana_sprd_ij(grid_ns, ij_count))
     allocate(lon_ij(ij_count))
     allocate(lat_ij(ij_count))
     allocate(mask_ij(ij_count))
     tbgscatter = timer_init("    bkg_scatter", TIMER_SYNC)
     call timer_start(tbgscatter)
     call letkf_mpi_ens2ij(gues, bkg_ij)
     call letkf_mpi_grd2ij(lon,  lon_ij)
     call letkf_mpi_grd2ij(lat,  lat_ij)
     call letkf_mpi_grd2ij(mask, mask_ij)
     call timer_stop(tbgscatter)
     deallocate(gues)
     deallocate(lon)
     deallocate(lat)


     ! caculate mean / spread
     tbgms = timer_init("    bkg_meansprd")
     call timer_start(tbgms)
     if (pe_isroot) then
        print *,""
        print *, "Calculating background mean / spread..."
     end if

     bkg_mean_ij = sum(bkg_ij, 1) / mem
     bkg_sprd_ij = 0
     do i=1,grid_ns
        do j=1,ij_count
           bkg_ij(:,i,j) = bkg_ij(:,i,j) - bkg_mean_ij(i,j)
           bkg_sprd_ij(i,j) =  sqrt(dot_product(bkg_ij(:,i,j),bkg_ij(:,i,j)) / mem)
        end do
     end do


     ! save the mean / spread out to files
     allocate(wrk(grid_nx, grid_ny, grid_ns))
     allocate(wrk2(ij_count))
     do i = 1,grid_ns
        wrk2 = bkg_mean_ij(i,:)
        call letkf_mpi_ij2grd(wrk2, wrk(:,:,i))
     end do
     if (pe_isroot) then
        write (filename, '(A)') 'OUTPUT/bkg_mean'
        print '(A,I5,3A)', " PROC ",pe_rank," is WRITING file: ",&
             trim(filename),trim(stateio_class%extension)
        call stateio_class%write(filename, wrk)
     end if
     do i = 1,grid_ns
        wrk2 = bkg_sprd_ij(i,:)
        call letkf_mpi_ij2grd(wrk2, wrk(:,:,i))
     end do
     if (pe_isroot) then
        write (filename, '(A)') 'OUTPUT/bkg_sprd'
        print '(A,I5,3A)', " PROC ",pe_rank," is WRITING file: ",&
             trim(filename),trim(stateio_class%extension)
        call stateio_class%write(filename, wrk)
     end if
     deallocate(wrk)
     deallocate(wrk2)
     call timer_stop(tbgms)

  end subroutine letkf_state_init



end module letkf_state
