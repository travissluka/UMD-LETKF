module global


  integer, parameter :: dp = kind(0.0d0)

  ! the following are read in from namelists
  integer :: mem
  real :: obsqc_maxstd

  !! @todo make the following configurable
  integer :: grid_x = 192
  integer :: grid_y = 94
  integer :: grid_z = 64
  integer :: grid_3d = 5
  integer :: grid_2d = 1
end module global
