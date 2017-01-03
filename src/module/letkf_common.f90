module letkf_common
  implicit none

  real, parameter :: pi = 4*atan(1.0)
  real, parameter :: re = 6371d3
  integer, parameter :: dp = kind(0.00)

  character(len=1024) :: nml_filename = "namelist.letkf"

  ! set by mpi module
  logical :: isroot = .true.

  ! the following are read in from namelists
  integer :: mem
  real :: obsqc_maxstd

end module letkf_common
