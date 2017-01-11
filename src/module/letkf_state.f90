module letkf_state
  !! performs I/O for the background and analysis states
  !!
  ! library modules
  use letkf_mpi
  use letkf_state_I
  use letkf_state_generic

  implicit none
  private

  public :: letkf_state_init

  ! public module variables
  integer, public  :: grid_nx
    !! number of grid points in the x direction
  integer, public  :: grid_ny
    !! number of grid points in the y direction
  integer, public, protected :: grid_nz
    !! number of vertical levels for 3D variables
  integer, public, protected :: grid_ns
    !! total number of 2D grid slices. This is equal to grid_nz*num_3D_vars + num_2D_vars
  real,    public, protected, allocatable :: lat_ij(:), lon_ij(:)
  real,    public, protected, allocatable :: bkg_ij(:,:,:)
  real,    public, protected, allocatable :: bkg_mean_ij(:,:)
  real,    public, protected, allocatable :: bkg_sprd_ij(:,:)
  real,    public,            allocatable :: ana_ij(:,:,:)
  real,    public,            allocatable :: ana_mean_ij(:,:)
  real,    public,            allocatable :: ana_sprd_ij(:,:)


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
    real, allocatable :: lat(:,:), lon(:,:)

    integer :: unit
    character(len=1024) :: filename
    integer :: m, i

    namelist /grid_def/ grid_nx, grid_ny, grid_nz, grid_ns

    if(pe_isroot) then
       print *, new_line('a'), &
            new_line('a'), '============================================================',&
            new_line('a'), 'letkf_state_init() : ',&
            new_line('a'), '============================================================'
    end if

    ! read in our section of the namelist
    open(newunit=unit, file=nml_filename)
    read(unit, nml=grid_def)
    close(unit)
    if (pe_isroot) then
       print grid_def
    end if


    !create the state io class, and initialize
    allocate(stateio_generic :: stateio_class)
    call stateio_class%init(grid_nx, grid_ny, grid_nz)
    allocate(lat(grid_nx,grid_ny))
    allocate(lon(grid_nx, grid_ny))
    call stateio_class%latlon(lat,lon)


    ! read in the background members that our process is responsible for
     if (pe_isroot) then
        print *, ""
        print *, "Reading ensemble background"
        print *, "------------------------------------------------------------"
        print *,"Reading ensemble background ..."
        print *, "slabs = ", grid_ns
        print "(A,I6,A,I6,A,I6,A)"," shape = (",grid_nx," x ",grid_ny," x ",grid_nz,")"
     end if
     allocate(gues(grid_nx, grid_ny, grid_ns, size(ens_list)))
     do m=1,size(ens_list)
        write (filename, '(A,I0.4,A)') 'INPUT/gues/',ens_list(m),'.nc'
        print *, "reading",trim(filename)
        call stateio_class%read(filename, gues(:,:,:,m))
     end do


     ! scatter grids to mpi procs
     if (pe_isroot) print *, "Scattering across procs..."
     allocate(bkg_ij(ij_count, grid_ns, mem))
     allocate(bkg_mean_ij(ij_count, grid_ns))
     allocate(bkg_sprd_ij(ij_count, grid_ns))
     allocate(ana_ij(ij_count, grid_ns, mem))
     allocate(ana_mean_ij(ij_count, grid_ns))
     allocate(ana_sprd_ij(ij_count, grid_ns))
     allocate(lon_ij(ij_count))
     allocate(lat_ij(ij_count))
     call letkf_mpi_ens2ij(gues, bkg_ij)
     call letkf_mpi_grd2ij(lon,  lon_ij)
     call letkf_mpi_grd2ij(lat,  lat_ij)
     deallocate(gues)
     deallocate(lon)
     deallocate(lat)


     ! caculate mean / spread
     if (pe_isroot) print *, "Calculating background mean / spread..."
     bkg_mean_ij = sum(bkg_ij, 3) / mem
     bkg_sprd_ij = 0
     do m = 1,mem
        bkg_ij(:,:,m) = bkg_ij(:,:,m) - bkg_mean_ij
        bkg_sprd_ij = bkg_sprd_ij + bkg_ij(:,:,m)*bkg_ij(:,:,m)
     end do
     bkg_sprd_ij = sqrt(bkg_sprd_ij/mem)


     ! save the mean / spread out to files
     allocate(wrk(grid_nx, grid_ny, grid_ns))
     do i = 1,grid_ns
        call letkf_mpi_ij2grd(bkg_mean_ij(:,i), wrk(:,:,i))
     end do
     if (pe_isroot) call stateio_class%write('OUTPUT/bkg_mean.nc', wrk)
     do i = 1,grid_ns
        call letkf_mpi_ij2grd(bkg_sprd_ij(:,i), wrk(:,:,i))
     end do
     if (pe_isroot) call stateio_class%write('OUTPUT/bkg_sprd.nc', wrk)
     deallocate(wrk)


  end subroutine letkf_state_init



end module letkf_state
