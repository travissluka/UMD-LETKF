module letkf_obs
  implicit none
  private

  public :: letkf_obs_init
  public :: observation
  public :: obsdef,   obsdef_read,  obsdef_getbyname,  obsdef_getbyid
  public :: platdef, platdef_read, platdef_getbyname, platdef_getbyid
  public :: obsio, I_obsio_write, I_obsio_read


  !------------------------------------------------------------
  
  integer, parameter :: dp=kind(0.0d0)


  
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
     integer :: qc
     !! Quality control flag. ( 0 = no errors, > 0 observation type
     !! dependant error)
  end type observation


  
  type obsdef
     !! The specification for a single user-defined observation type
     integer :: id
     !! unique identifier
     character(len=:), allocatable :: name_short
     !! short (3-4 characters) name of observation type
     character(len=:), allocatable :: name_long
     !! longer more descriptive name
     character(len=:), allocatable :: units     
  end type obsdef


  

  type platdef
     !!specification for a single user-defined observation platform type
     integer :: id
     !! unique id for the platform type
     character(len=:), allocatable :: name_short
     !! short (3-4 character) name of platform type
     character(len=:), allocatable :: name_long
     !! longer more descriptive name of platform type
  end type platdef

  

  type, abstract :: obsio
     !! Abstract base class for observation file reading and writing.
     !! All user-defined, and built-in, obs file I/O modules for
     !! specific file types should be built from a class extending this
   contains
     procedure(I_obsio_write), deferred :: write
     !! write a list of observatiosn to the given file     
     procedure(I_obsio_read),  deferred :: read
     !! read a list of observatiosn from the given file
  end type obsio
  
  abstract interface
     subroutine I_obsio_write(self, file, obs, iostat)
       !! interface for procedures to write observation data
       import obsio
       import observation
       class(obsio) :: self
       character(len=*), intent(in)   :: file
       !! filename to write observations to
       class(observation), intent(in) :: obs(:)
       !! list of 1 or more observatiosn to write
       integer, optional, intent(out) :: iostat
     end subroutine I_obsio_write

     subroutine I_obsio_read(self, file, obs, iostat)
       !! interface for procedures to load observation data
       import observation
       import obsio
       class(obsio) :: self
       character(len=*), intent(in) :: file
       class(observation), intent(out) :: obs(:)
       integer, optional, intent(out) :: iostat
     end subroutine I_obsio_read
  end interface


  !------------------------------------------------------------
  type(obsdef),  allocatable ::  obsdef_list(:)
  type(platdef), allocatable :: platdef_list(:)

  
  !------------------------------------------------------------

  
contains

  subroutine letkf_obs_init(reader, obsdef_file, platdef_file)
    !! @warning this has not been implemented
    class(obsio), intent(in) :: reader
    character(len=*), optional, intent(in) :: obsdef_file
    !! observation definition file to read in. By default `letkf.obsdef`
    !! will be used.
    character(len=*), optional, intent(in) :: platdef_file
    !! observation platform definition file to read in. By default
    !! `letkf.platdef` will be used.

    print *, ''
    print *, 'LETKF observation configuration'
    print *, '----------------------------------------'
    print *, '  observation definition file: "', trim(obsdef_file),'"'
    print *, '  platform    definition file: "', trim(platdef_file),'"'
    print *, '  I/O format: ???'

    call obsdef_read(obsdef_file)
    call platdef_read(platdef_file)
    
  end subroutine letkf_obs_init


  
  subroutine obsdef_read(file)    
    !! @warning this has not been implemented        
    character(len=*), intent(in) :: file
  end subroutine obsdef_read


  
  subroutine platdef_read(file)
    !! @warning this has not been implemented        
    character(len=*), intent(in) :: file
  end subroutine platdef_read


  
  function obsdef_getbyid(id) result(res)
    !! returns information about an observation type given its id number.
    !! An error is thrown if the id is not found
    !! @warning this has not been implemented    
    integer, intent(in) :: id
    type(obsdef) :: res

    integer :: i
  end function obsdef_getbyid

  

  function obsdef_getbyname(name) result(res)
    !! @warning this has not been implemented    
    character(len=*), intent(in) :: name
    type(obsdef) :: res

    integer :: i
  end function obsdef_getbyname


  
  function platdef_getbyid(id) result(res)
    !! @warning this has not been implemented    
    integer, intent(in) :: id
    type(platdef) :: res

    integer :: i
  end function platdef_getbyid


  
  function platdef_getbyname(name) result(res)
    !! @warning this has not been implemented
    
    character(len=*), intent(in) :: name
    type(platdef) :: res

    integer :: i
  end function platdef_getbyname

end module letkf_obs
