module letkf_state_n
  !! performs I/O for the background and analysis states
  !!
  ! library modules
  use letkf_common
  use letkf_state_types
  use letkf_state_generic
  use str

  implicit none
  private

  ! public module variables
  integer, public, protected :: grid_nx
  integer, public, protected :: grid_ny
  integer, public, protected :: grid_nz
  integer, public, protected :: grid_nvar
  real,    public, protected, allocatable :: state_2d(:,:,:)
  integer, public, protected, allocatable :: state_lvl(:)
  integer, public, protected, allocatable :: state_var(:)

  ! private module variables
  type(statedef), allocatable :: statedef_list(:)


contains



  subroutine letkf_state_init()
    integer :: unit

    ! read in our section of the namelist
  end subroutine letkf_state_init


end module letkf_state_n
