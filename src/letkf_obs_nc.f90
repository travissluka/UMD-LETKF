module letkf_obs_nc
  use netcdf
  use letkf_obs

  implicit none
  private

  public :: obsio_nc
  
  !------------------------------------------------------------

  integer, parameter :: dp=kind(0.0d0)

  type, extends(obsio) :: obsio_nc
   contains
     procedure :: write => obsio_nc_write
     procedure :: read => obsio_nc_read
  end type obsio_nc


contains

  subroutine obsio_nc_write(self, file, obs, iostat)
    class(obsio_nc) :: self
    character(len=*), intent(in) :: file
    class(observation), intent(in) :: obs(:)
    integer, optional, intent(out) :: iostat
  end subroutine obsio_nc_write



  subroutine obsio_nc_read(self, file, obs, iostat)
    class(obsio_nc) :: self
    character(len=*), intent(in) :: file
    class(observation), intent(out) :: obs(:)
    integer, optional, intent(out) :: iostat
  end subroutine obsio_nc_read

  
end module letkf_obs_nc
