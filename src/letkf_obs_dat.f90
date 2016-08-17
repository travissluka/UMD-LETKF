module letkf_obs_dat
  use letkf_obs

  implicit none
  private

  public :: obsio_dat
  
  !------------------------------------------------------------

  integer, parameter :: dp=kind(0.0d0)

  type, extends(obsio) :: obsio_dat
   contains
     procedure :: write => obsio_dat_write
     procedure :: read => obsio_dat_read
  end type obsio_dat


contains

  subroutine obsio_dat_write(self, file, obs, iostat)
    class(obsio_dat) :: self
    character(len=*), intent(in) :: file
    class(observation), intent(in) :: obs(:)
    integer, optional, intent(out) :: iostat
  end subroutine obsio_dat_write



  subroutine obsio_dat_read(self, file, obs, iostat)
    class(obsio_dat) :: self
    character(len=*), intent(in) :: file
    class(observation), intent(out) :: obs(:)
    integer, optional, intent(out) :: iostat
  end subroutine obsio_dat_read

  
end module letkf_obs_dat
