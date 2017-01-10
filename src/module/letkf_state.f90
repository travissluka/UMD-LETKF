module letkf_state
  !! performs I/O for the background and analysis states
  !!
  ! library modules
  use letkf_common
  use letkf_state_I
  use letkf_state_generic

  implicit none
  private

  public :: letkf_state_init

  ! public module variables
  integer, public, protected :: grid_nx
  integer, public, protected :: grid_ny
  integer, public, protected :: grid_nz
  integer, public, protected :: grid_ns

  real, public, protected,  allocatable :: lat(:,:), lon(:,:)

!  real,    public, protected, allocatable :: state_2d(:,:,:)
    !! a series of 2D slabs that represent the 2d AND 3d state
    !! the level and variable type of each slab is given by
    !! state_lvl and state_var respectively (x,y,lvl)
!  integer, public, protected, allocatable :: state_lvl(:)
    !! the level of the corresponding slab of state_2d, or 0
    !! if a 2D variable
!  integer, public, protected, allocatable :: state_var(:)
    !! the variable ID of

  ! private module variables
!  type(statedef), allocatable :: statedef_list(:)


  class(stateio), public, pointer :: stateio_class


contains



  subroutine letkf_state_init()
    integer :: unit

    namelist /grid_def/ grid_nx, grid_ny, grid_nz, grid_ns

    ! read in our section of the namelist
    open(newunit=unit, file=nml_filename)
    read(unit, nml=grid_def)
    close(unit)


    allocate(stateio_generic :: stateio_class)
    call stateio_class%init(grid_nx, grid_ny, grid_nz)
    allocate(lat(grid_nx,grid_ny))
    allocate(lon(grid_nx, grid_ny))
    call stateio_class%latlon(lat,lon)
  end subroutine letkf_state_init


end module letkf_state
