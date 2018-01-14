module letkf_state
  !! performs I/O for the background and analysis states

  use timing
  use letkf_mpi
  use letkf_obs
  use kdtree

  implicit none
  private

  ! public module methods
  !------------------------------------------------------------
  public :: letkf_state_init
  public :: letkf_state_read
  public :: letkf_state_write
  public :: letkf_state_register
  public :: stateio, stateioptr


  ! public module variables
  ! ------------------------------------------------------------
  integer, public, protected :: grid_nx
  integer, public, protected :: grid_ny
  integer, public, protected :: grid_ns

  real,    public, protected, allocatable :: lat_ij(:)
    !! ***Size is ( [[letkf_mpi:ij_count]] ) ***
  real,    public, protected, allocatable :: lon_ij(:)
    !! ***Size is ( [[letkf_mpi:ij_count]] ) ***
  logical,  public, protected, allocatable :: mask_ij(:)

  real,    public, protected, allocatable :: bkg_ij(:,:,:)
    !! ***Shape is ( [[letkf_mpi:mem]], [[letkf_mpi:grid_ns]], [[letkf_mpi:ij_count]] ) ***
  real,    public, protected, allocatable :: bkg_mean_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:grid_ns]], [[letkf_mpi:ij_count]] ) ***
  real,    public, protected, allocatable :: bkg_sprd_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:grid_ns]], [[letkf_mpi:ij_count]] ) ***

  real,    public,            allocatable :: ana_ij(:,:,:)
    !! ***Shape is ( [[letkf_mpi:mem]], [[letkf_mpi:grid_ns]], [[letkf_mpi:ij_count]] ) ***
  real,    public,            allocatable :: ana_mean_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:grid_ns]], [[letkf_mpi:ij_count]] ) ***
  real,    public,            allocatable :: ana_sprd_ij(:,:)
    !! ***Shape is ( [[letkf_mpi:grid_ns]], [[letkf_mpi:ij_count]] ) ***



  !! ------------------------------------------------------------



  type, abstract :: stateio
     !! abstract base class for reading and writing of state files
     integer :: i !don't worry about this, trust me
   contains
     procedure(I_stateio_getstr), deferred :: get_name
     procedure(I_stateio_getstr), deferred :: get_desc
     procedure(I_stateio_init),   deferred :: init
     procedure(I_stateio_read),   deferred :: read
     procedure(I_stateio_write),  deferred :: write
!     procedure(I_stateio_latlon), deferred :: latlon
!     procedure(I_stateio_mask),   deferred :: mask
  end type stateio

  abstract interface
     function I_stateio_getstr(self) result(str)
       import stateio
       class(stateio) :: self
       character(:), allocatable :: str
     end function I_stateio_getstr

     subroutine I_stateio_init(self, lat, lon, mask)
       import stateio
       class(stateio) :: self
       real,    intent(inout) :: lat(:,:), lon(:,:)
       logical, intent(inout) :: mask(:,:)
     end subroutine I_stateio_init

     subroutine I_stateio_read(self, ens, state)
       import stateio
       class(stateio) :: self
!       character(len=*),  intent(in)  :: filename
       integer, intent(in) :: ens
       real, intent(out) :: state(:,:,:)
     end subroutine I_stateio_read

     subroutine I_stateio_write(self, filename, state)
       import stateio
       class(stateio) :: self
       character(len=*),  intent(in)  :: filename
       real, intent(in)  :: state(:,:,:)
     end subroutine I_stateio_write

     subroutine I_stateio_latlon(self, lat, lon)
       import stateio
       class(stateio) :: self
       real, intent(inout) :: lat(:,:), lon(:,:)
     end subroutine I_stateio_latlon

     subroutine I_stateio_mask(self, mask)
       import stateio
       class(stateio) :: self
       logical, intent(inout) :: mask(:,:)
     end subroutine I_stateio_mask

  end interface


  !! ------------------------------------------------------------


  type stateioptr
     class(stateio), pointer :: p
  end type stateioptr


  ! private module variables
  !------------------------------------------------------------
  !TODO, remove public from this
  public :: stateio_class

  ! registration of built-in and user-defined stateio classes
  ! TODO allow for 2D/3D mask, 1D/2D lat/lon
  integer, parameter      :: stateio_reg_max = 100
  integer                 :: stateio_reg_num = 0
  type(stateioptr)        :: stateio_reg(stateio_reg_max)
  class(stateio), pointer :: stateio_class
  real,    allocatable    :: lat(:,:), lon(:,:)
  logical, allocatable    :: mask(:,:)

  !TODO, make these private, and/or protected
  character(len=12), public, allocatable :: state_var(:)
  !! name of a state variable
  integer,          public, allocatable :: slab_var(:)
  !! state variable for each slab, gives an index in the [[letkf_state::state_var]] array
  integer,          public, allocatable :: slab_lvl(:)
  !! level number for each slab




contains




  !================================================================================
  !================================================================================
  subroutine letkf_state_init(nml_filename)
    !! Initialize the state vector by reading in the background ensemble
    !! and distributing across cores

    character(len=*), intent(in) :: nml_filename
      !! filename of namelist to read in for the ***grid_def*** section
    integer :: unit, i   
    character(len=:), allocatable :: ioclass


    namelist /letkf_state/ ioclass, grid_nx, grid_ny, grid_ns


    if(pe_isroot) then
       print "(A)", ""
       print "(A)", ""
       print "(A)", '============================================================'
       print "(A)", ' letkf_state_init() : '
       print "(A)", '============================================================'
    end if


    ! read in our section of the namelist
    allocate(character(1024) :: ioclass)
    open(newunit=unit, file=nml_filename)
    read(unit, nml=letkf_state)
    close(unit)
    ioclass = trim(ioclass)
    if (pe_isroot) then
       print letkf_state
    end if
    

    ! print a list of all stateio classes that have been registered
    if(pe_isroot) then
       print *, ""
       print *, "List of stateio classes registered:"
       do i=1,stateio_reg_num
          print "(A,A,3A)", " * ", toupper(stateio_reg(i)%p%get_name()), &
               " (",stateio_reg(i)%p%get_desc(), ")"
       end do
       print *, ""
    end if


    ! determine the stateio class to use
    nullify(stateio_class)
    do i=1,stateio_reg_num
       if(trim(toupper(stateio_reg(i)%p%get_name())) == trim(toupper(ioclass))) then
          stateio_class => stateio_reg(i)%p
          exit
       end if
    end do
    if(.not. associated(stateio_class)) then
       if(pe_isroot) &
            print *, 'ERROR: stateio class "',toupper(trim(ioclass)), &
            '" not found. Check that the name is in the list of registered classes'
       stop 1
    end if
    if(pe_isroot) print *, 'Using "', trim(toupper(stateio_class%get_name())),'"'


    ! initialize the user-defined grid
    !TODO, should this be done inside the stateio instance (if we want stateio to determine
    ! the grid_ns/nx/ny parameters)?
    allocate(slab_var(grid_ns))
    allocate(slab_lvl(grid_ns))
    slab_var = -1
    slab_lvl = -1
    ! TODO, only one process needs to read in lat/lon/mask and scatter it
    allocate(lat (grid_nx, grid_ny))
    allocate(lon (grid_nx, grid_ny))
    allocate(mask(grid_nx, grid_ny))
    call stateio_class%init(lat, lon, mask) 

    ! tell the mpi module about the grid layout
    call letkf_mpi_setgrid(grid_nx, grid_ny, grid_ns)


    ! scatter grids to mpi procs
    if (pe_isroot) then
       print *,""
       print *, "Scattering grid definition across procs..."
    end if
    allocate(lon_ij(ij_count))
    allocate(lat_ij(ij_count))
    allocate(mask_ij(ij_count))
    !TODO, get the generic mpi_grd2ij interface working
    call letkf_mpi_grd2ij_real(lat,  lat_ij)
    call letkf_mpi_grd2ij_real(lon,  lon_ij)
    call letkf_mpi_grd2ij_logical(mask, mask_ij)
    ! TODO, deallocate lat/lon here,
    ! Err, can't right now because it is needed when the state bg is read in 
    ! if doing test observations
!    deallocate(lon)
!    deallocate(lat)

    ! print out statistics about the grid layout
    if (pe_isroot) then
       print *,""
       print *, "grid definition:"
       print *, "------------------------------------------------------------"
       print '(A,I6,A,I6)', " horizontal grid size", grid_nx," X",grid_ny
       print '(A,F6.1,A,F6.1)', "  longitude range: ", minval(lon), " to ",maxval(lon)
       print '(A,F6.1,A,F6.1)', "  latitude  range: ", minval(lat), " to ",maxval(lat)       
       print '(A,I6)', " vertical/variable slabs:", grid_ns
!       do i =1, grid_ns
!          print '(A,I4,A,A10,A,I6)', "  slab", i, "       var=",&
!               letkf_state_getvarname(slab_var(i)), "      lvl=",slab_lvl(i)
!       end do
    end if

  end subroutine letkf_state_init
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_state_read
    real, allocatable :: gues(:,:,:,:)
    integer :: tbgscatter, tbgread
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
       print "(A)", ""
       print "(A)", ""
       print "(A)", '============================================================'
       print "(A)", ' letkf_state_read() :'
       print "(A)", '============================================================'
    end if


    if (pe_isroot) then
       print *, ""
       print *, "Reading ensemble background..."
    end if


    ! console output synchronization
    call letkf_mpi_barrier(syncio=.true.)


    ! read in the background members that our process is responsible for
    allocate(gues(grid_nx, grid_ny, grid_ns, size(ens_list)))
    tbgread = timer_init("    bkg_read", TIMER_SYNC)
    call timer_start(tbgread)
    do m=1,size(ens_list)
!       write (filename, '(A,I0.4,A)') 'INPUT/gues/',ens_list(m)
       print '(A,I5,A,I5)', " PROC ",pe_rank," is READING ens member: ", ens_list(m)

       call stateio_class%read(ens_list(m), gues(:,:,:, m))
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

       ! all done with the kd-tree
       deallocate(llkd_lons)
       deallocate(llkd_lats)
       deallocate(llkd_x)
       deallocate(llkd_y)
       call kd_free(llkd_root)
    end if


    ! scatter background, 
    if(pe_isroot) then
       print *,""
       print *, "scattering background..."
    end if
    tbgscatter = timer_init("    bkg_scatter", TIMER_SYNC)

    call timer_start(tbgscatter)
    allocate(bkg_ij(mem, grid_ns, ij_count))
    allocate(bkg_mean_ij(grid_ns, ij_count))
    allocate(bkg_sprd_ij(grid_ns, ij_count))
    allocate(ana_ij(mem, grid_ns, ij_count))
    allocate(ana_mean_ij(grid_ns, ij_count))
    allocate(ana_sprd_ij(grid_ns, ij_count))
    call letkf_mpi_ens2ij(gues, bkg_ij)
    call timer_stop(tbgscatter)


    ! cleanup no longer needed arrays
    ! TODO, is lat/lon needed all the way to this point?
    deallocate(gues)
    deallocate(lon)
    deallocate(lat)


    ! caculate mean / spread
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

  end subroutine letkf_state_read
  !================================================================================



  !================================================================================
  !================================================================================
  subroutine letkf_state_write()
    integer :: i, m 
    integer :: rank_bm, rank_bs, rank_am, rank_as

    real, allocatable :: wrkij(:)
    real, allocatable :: wrk_bm(:,:,:)
    real, allocatable :: wrk_bs(:,:,:)
    real, allocatable :: wrk_am(:,:,:)
    real, allocatable :: wrk_as(:,:,:)
    real, allocatable :: wrk4(:,:,:,:)
    character(len=1024) :: filename

    integer :: t_mpi, t_write

    t_mpi   = timer_init("  output_scatter")
    t_write = timer_init("  output_write")

    if(pe_isroot) then
       print *, ""
       print *, ""
       print "(A)", '============================================================'
       print "(A)", ' letkf_state_write() :'
       print "(A)", '============================================================'
    end if

    ! determine which processes are going to output which of the bkg/ana mean/spread
    ! so that we can write them in parallel
    rank_bm = mod(pe_size  , min(4,pe_size))
    rank_bs = mod(pe_size+1, min(4,pe_size))
    rank_am = mod(pe_size+2, min(4,pe_size))
    rank_as = mod(pe_size+3, min(4,pe_size))
    if (pe_isroot) then
       print *, ""
       print *, "Saving ana/bkg mean/sprd"
       print '(A,I5,3A)', " PROC ",rank_bm, " is WRITING file: OUTPUT/bkg_mean"
       print '(A,I5,3A)', " PROC ",rank_bs, " is WRITING file: OUTPUT/bkg_sprd"
       print '(A,I5,3A)', " PROC ",rank_am, " is WRITING file: OUTPUT/ana_mean"
       print '(A,I5,3A)', " PROC ",rank_as, " is WRITING file: OUTPUT/ana_sprd"
    end if

    ! gather the mean/spread
    if(pe_isroot) print *, "scattering..."
    allocate(wrkij(ij_count))
    call timer_start(t_mpi)
    if (pe_rank == rank_bm) allocate(wrk_bm(grid_nx, grid_ny, grid_ns))
    do i=1,grid_ns
       wrkij = bkg_mean_ij(i,:)
       call letkf_mpi_ij2grd(wrkij, wrk_bm(:,:,i), rank_bm)
    end do    
    if (pe_rank == rank_bs) allocate(wrk_bs(grid_nx, grid_ny, grid_ns))
    do i=1,grid_ns
       wrkij = bkg_sprd_ij(i,:)
       call letkf_mpi_ij2grd(wrkij, wrk_bs(:,:,i), rank_bs)
    end do    
    if (pe_rank == rank_am) allocate(wrk_am(grid_nx, grid_ny, grid_ns))
    do i=1,grid_ns
       wrkij = ana_mean_ij(i,:)
       call letkf_mpi_ij2grd(wrkij, wrk_am(:,:,i), rank_am)
    end do    
    if (pe_rank == rank_as) allocate(wrk_as(grid_nx, grid_ny, grid_ns))
    do i=1,grid_ns
       wrkij = ana_sprd_ij(i,:)
       call letkf_mpi_ij2grd(wrkij, wrk_as(:,:,i), rank_as)
    end do    
    call timer_stop(t_mpi)
    deallocate(wrkij)


    ! write out the mean/spread
    if(pe_isroot) print *, "writing..."
    call timer_start(t_write)
    if (pe_rank == rank_bm) then
       call stateio_class%write('OUTPUT/bkg_mean', wrk_bm)       
       deallocate(wrk_bm)
    end if
    if (pe_rank == rank_bs) then
       call stateio_class%write('OUTPUT/bkg_sprd', wrk_bs)
       deallocate(wrk_bs)
    end if
    if (pe_rank == rank_am) then
       call stateio_class%write('OUTPUT/ana_mean', wrk_am)
       deallocate(wrk_am)
    end if
    if (pe_rank == rank_as) then
       call stateio_class%write('OUTPUT/ana_sprd', wrk_as)
       deallocate(wrk_as)
    end if
    call timer_stop(t_write)
    call letkf_mpi_barrier(.true.)
    if(pe_isroot) print *, "Done."


    ! analysis ensemble members
    if(pe_isroot) then
       print *,""
       print *, "Saving analyais ensemble members"
       print *, "mpi collect..."
    end if

    call timer_start(t_mpi)
    allocate(wrk4(grid_nx,grid_ny, grid_ns, size(ens_list)))
    call letkf_mpi_ij2ens(ana_ij, wrk4)
    call timer_stop(t_mpi)

    call timer_start(t_write)
    do m=1,size(ens_list)
       write (filename, '(A,I0.4,A)') 'OUTPUT/', ens_list(m)
       print '(A,I5,3A)', " PROC ",pe_rank, " is WRITING file: ",trim(filename)
       call stateio_class%write(filename, wrk4(:,:,:,m))
    end do
    call timer_stop(t_write)

    call letkf_mpi_barrier()
    deallocate(wrk4)


    
  end subroutine letkf_state_write
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_state_register(cls)
    class(stateio), pointer :: cls
    integer :: i

    if(stateio_reg_num == stateio_reg_max) then
       print *, "ERROR: too many stateio classes registered"
       stop 1
    end if

    ! ensure a class of that same name doesn't already exist
    do i=1, stateio_reg_num
       if (toupper(stateio_reg(i)%p%get_name()) == toupper(cls%get_name())) then
          print *, "ERROR: can't regsiter stateio class '",toupper(cls%get_name()),&
               "', a class by that name has already been registered"
          stop 1
       end if
    end do

    stateio_reg_num = stateio_reg_num + 1
    stateio_reg(stateio_reg_num)%p => cls
  end subroutine letkf_state_register
  !================================================================================




  !================================================================================
  !================================================================================
  function letkf_state_getvarid(str) result(id)
    character(*), intent(in) :: str
    integer :: id
    id = -1
  end function letkf_state_getvarid
  !================================================================================




  !================================================================================
  !================================================================================
  function letkf_state_getvarname(id) result(str)
    integer, intent(in) :: id
    character(len=:), allocatable :: str
    if (id <= 0 .or. id > size(state_var)) then
       str = "ERR"
       print *, "ERROR, trying to get non-existant state variable name in slot ", id
       stop 1
    end if   
    str = state_var(id)
  end function letkf_state_getvarname
  !================================================================================



  !================================================================================
  !================================================================================
  function toupper(in_str) result(out_str)
    character(*), intent(in) :: in_str
    character(len(in_str)) :: out_str
    integer :: i
    integer, parameter :: offset = 32

    out_str = in_str
    do i=1, len(out_str)
       if(out_str(i:i) >= "a" .and. out_str(i:i) <= "z") then
          out_str(i:i) = achar(iachar(out_str(i:i))-offset)
       end if
    end do
  end function toupper
  !================================================================================

end module letkf_state
