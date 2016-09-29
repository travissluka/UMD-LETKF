module global

  integer :: mem

  integer, parameter :: dp = kind(0.0d0)

  !TODO: read the following in from namelist
  real :: obsqc_maxstd = 3.0

end module global
