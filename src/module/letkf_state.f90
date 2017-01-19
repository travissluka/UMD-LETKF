module letkf_state
  !! performs I/O for the background and analysis states
  !!
  ! library modules
  use timing
  use letkf_mpi
  use letkf_state_I
  use letkf_state_generic

  use iso_fortran_env
  implicit none
  private

  public :: letkf_state_init

  ! public module variables
  integer, public  :: grid_nx
    !! number of grid points in the x direction
  integer, public  :: grid_ny
    !! number of grid points in the y direction
  integer, public  :: grid_ns
    !! total number of 2D grid slices. This is equal to grid_nz*num_3D_vars + num_2D_vars
  integer, public, protected :: grid_nz
    !! number of vertical levels for 3D variables

  real,    public, protected, allocatable :: lat_ij(:)
    !! ***Size is ( [[letkf_mpi:ij_count]] ) ***
  real,    public, protected, allocatable :: lon_ij(:)
    !! ***Size is ( [[letkf_mpi:ij_count]] ) ***

  real,    public, protected, allocatable :: bkg_ij(:,:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_state:grid_ns]], [[letkf_mpi:mem]] ) ***
  real,    public, protected, allocatable :: bkg_mean_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_state:grid_ns]] ) ***
  real,    public, protected, allocatable :: bkg_sprd_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_state:grid_ns]] ) ***

  real,    public,            allocatable :: ana_ij(:,:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_state:grid_ns]], [[letkf_mpi:mem]] ) ***
  real,    public,            allocatable :: ana_mean_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_state:grid_ns]] ) ***
  real,    public,            allocatable :: ana_sprd_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:ij_count]], [[letkf_state:grid_ns]] ) ***


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
    real, allocatable :: lat(:,:), lon(:,:)

    integer :: unit, tbgscatter, tbgread, tbgms
    character(len=1024) :: filename
    integer :: m, i, j

    namelist /grid_def/ grid_nx, grid_ny, grid_nz, grid_ns

    if(pe_isroot) then
       print *, new_line('a'), &
            new_line('a'), '============================================================',&
            new_line('a'), ' letkf_state_init() : ',&
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
        print *, "Reading ensemble background..."
        ! print *, "------------------------------------------------------------"
        ! print *, "slabs = ", grid_ns
        ! print "(A,I6,A,I6,A,I6,A)"," shape = (",grid_nx," x ",grid_ny," x ",grid_nz,")"
        ! print *,""
     end if

     flush(output_unit)
     call system('sleep 0')
     call letkf_mpi_barrier()

     ! read in the files
     allocate(gues(grid_nx, grid_ny, grid_ns, size(ens_list)))
     tbgread = timer_init("    bkg_read", TIMER_SYNC)
     call timer_start(tbgread)
     do m=1,size(ens_list)
        write (filename, '(A,I0.4,A)') 'INPUT/gues/',ens_list(m)
        print '(A,I5,3A)', " PROC ",pe_rank," is READING file: ",trim(filename),trim(stateio_class%extension)
        call stateio_class%read(filename, gues(:,:,:,m))
     end do
     call timer_stop(tbgread)


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
     tbgscatter = timer_init("    bkg_scatter", TIMER_SYNC)
     call timer_start(tbgscatter)
     call letkf_mpi_ens2ij(gues, bkg_ij)
     call letkf_mpi_grd2ij(lon,  lon_ij)
     call letkf_mpi_grd2ij(lat,  lat_ij)
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
        print '(A,I5,3A)', " PROC ",pe_rank," is WRITING file: ",trim(filename),trim(stateio_class%extension)
        call stateio_class%write(filename, wrk)
     end if
     do i = 1,grid_ns
        wrk2 = bkg_sprd_ij(i,:)
        call letkf_mpi_ij2grd(wrk2, wrk(:,:,i))
     end do
     if (pe_isroot) then
        write (filename, '(A)') 'OUTPUT/bkg_sprd'
        print '(A,I5,3A)', " PROC ",pe_rank," is WRITING file: ",trim(filename),trim(stateio_class%extension)
        call stateio_class%write(filename, wrk)
     end if
     deallocate(wrk)
     deallocate(wrk2)
     call timer_stop(tbgms)

  end subroutine letkf_state_init



end module letkf_state
