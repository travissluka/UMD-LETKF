module letkf_obs_I
  implicit none
  private

  integer, parameter :: dp = kind(0.0)

  ! public types
  public :: observation
  public :: obsio, obsioptr



  type observation
     !! information for a single observation
     integer :: id = -1
     !! observation type id, as specified in the `obsdef` configuration file.
     !! A value of -1 indicates no observation data has been read in.
     integer :: plat
     !! Platform type id, as specified in the `platdef` configuration file.
     real(dp) :: lat
     !! Latitude (degrees)
     real(dp) :: lon
     !! Longitude (degrees)
     real(dp) :: depth
     !! Depth, in whatever units is appropriate for the domain.
     !! If this is a 2D observation, this value is ignored.
     real(dp) :: time
     !! Time (hours) from base time
     real(dp) :: val
     !! Observation value
     real(dp) :: err
     !! standard deviation of observation error
!     integer :: qc
     !! Quality control flag. ( 0 = no errors, > 0 removed by obs op, < 0
     !! removed by LETKF qc
  end type observation



  type, abstract :: obsio
     !! Abstract base class for observation file reading and writing.
     !! All user-defined, and built-in, obs file I/O modules for
     !! specific file types should be built from a class extending this
     integer :: i ! don't worry about this
   contains
     procedure(I_obsio_getstr),  deferred :: get_name
     procedure(I_obsio_getstr),  deferred :: get_desc
     procedure(I_obsio_write),   deferred :: write
     !! write a list of observatiosn to the given file
     procedure(I_obsio_read),    deferred :: read
     !! read a list of observatiosn from the given file
  end type obsio

  abstract interface

     
     function I_obsio_getstr(self)
       import obsio
       class(obsio) :: self
       character(:), allocatable :: I_obsio_getstr
     end function I_obsio_getstr

     
     subroutine I_obsio_init(self)
       import obsio
       class(obsio) :: self
     end subroutine I_obsio_init

     
     subroutine I_obsio_write(self, file, obs, iostat)
       !! interface for procedures to write observation data
       import obsio
       import observation
       class(obsio) :: self
       character(len=*), intent(in)   :: file
       !! filename to write observations to
       type(observation), intent(in) :: obs(:)
       !! list of 1 or more observatiosn to write
       integer, optional, intent(out) :: iostat
     end subroutine I_obsio_write

     
     subroutine I_obsio_read(self, file, obs, obs_innov, obs_qc, iostat)
       !! interface for procedures to load observation data
       import observation
       import obsio
       import dp
       class(obsio) :: self
       character(len=*), intent(in) :: file

       type(observation),allocatable, intent(out) :: obs(:)
       real(dp), allocatable, intent(out) :: obs_innov(:)
       integer,  allocatable, intent(out) :: obs_qc(:)

       integer, optional, intent(out) :: iostat
     end subroutine I_obsio_read
  end interface


  type obsioptr
     class(obsio), pointer :: p
  end type obsioptr


end module letkf_obs_I
