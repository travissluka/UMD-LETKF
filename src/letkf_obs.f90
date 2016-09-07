module letkf_obs
  implicit none
  private

  public :: letkf_obs_init
  public :: observation
  public :: obsio, I_obsio_write, I_obsio_read

  ! observation definitions
  public :: obsdef
  public :: obsdef_read,  obsdef_getbyname,  obsdef_getbyid

  !platform definitions
  public :: platdef
  public :: platdef_read, platdef_getbyname, platdef_getbyid
 

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
   contains
     procedure :: print => obsdef_print
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

  subroutine letkf_obs_init(obsdef_file, platdef_file)
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
    character(len=*), intent(in) :: file

    integer :: unit, pos, iostat
    character(len=1024) :: line, s_tmp
    logical :: ex
    type(obsdef) :: new_ob
    integer, parameter :: MAX_OBSDEF = 1024
    type(obsdef) :: obsdef_list_tmp(MAX_OBSDEF)
    integer :: obsdef_list_tmp_len
    integer :: i,j

    print *, ""
    print *, "Observation Definition File"

    ! make sure the file exists
    inquire(file=file, exist=ex)
    if (.not. ex) then
       print *, 'ERROR: file does not exists: "',trim(file),'"'
       stop 1
    end if

    obsdef_list_tmp_len = 0
    
    ! open it up for reading
    print *, 'reading file: "',trim(file),'"'
    open(newunit=unit, file=file, action='read')
    do while (1==1)
       ! read a new line
       read(unit, '(A)', iostat=iostat) line
       if (iostat < 0) exit
       if (iostat > 0) then
          print *, 'ERROR: there was a problem reading "', &
            trim(file), '" error code: ', iostat
            stop 1
       end if

       ! ignore comments
       s_tmp = adjustl(line)
       if (s_tmp(1:1) == '#') cycle

       ! ignore empty lines
       if (len(trim(adjustl(line))) == 0) cycle

       ! read ID
       line = adjustl(line)
       pos = findspace(line)
       read(line(1:pos), *) new_ob%id

       ! read short name
       line = adjustl(line(pos+1:))
       pos = findspace(line)
       new_ob%name_short = toupper(trim(line(:pos-1)))

       ! read units
       line = adjustl(line(pos+1:))
       pos = findspace(line)
       new_ob%units = trim(line(:pos-1))

       ! read long name
       line = adjustl(line(pos+1:))
       new_ob%name_long = trim(line)

       ! add this one to the list
       obsdef_list_tmp_len = obsdef_list_tmp_len + 1
       if (obsdef_list_tmp_len > MAX_OBSDEF) then
          print *, 'ERROR, there are more observation definitions ',&
               'than there is room for. Check the "',trim(file),&
               '" file, and/or increase MAX_OBSDEF.'
          print *, 'MAX_OBSDEF = ',MAX_OBSDEF
          stop 1
       end if
       obsdef_list_tmp(obsdef_list_tmp_len) = new_ob         
    end do
    
    close (unit)
    allocate(obsdef_list(obsdef_list_tmp_len))
    obsdef_list = obsdef_list_tmp(1:obsdef_list_tmp_len)

    ! write summary
    print *, 'obs defined = ', size(obsdef_list)
    print "(A6,A10,A10,A5,A)", "ID", "NAME", "UNITS", "","FULL NAME"
    do pos=1, size(obsdef_list)
       call obsdef_list(pos)%print()
    end do

    ! check for duplicate ID / short name
    i = 1
    do while(i < size(obsdef_list))
       j = i + 1
       do while(j < size(obsdef_list))
          if ( obsdef_list(i)%id == obsdef_list(j)%id) then
             print *, "ERROR: multiple definitions for ID = ",&
                  obsdef_list(j)%id
             stop 1
          else if( obsdef_list(i)%name_short==obsdef_list(j)%name_short) then
             print *, "ERROR: multiple definitions for NAME = ",&
                  obsdef_list(j)%name_short
             stop 1
          end if          
          j = j + 1
       end do
       i = i + 1
    end do

    ! all done
    print *, ""
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


  subroutine obsdef_print(ob)
    class(obsdef), intent(in) :: ob
    print "(I6,A10,A10,A5,A)", ob%id, ob%name_short, &
         ob%units, " ", ob%name_long
  end subroutine obsdef_print


  function tolower(in_str) result(out_str)
    character(*), intent(in) :: in_str
    character(len(in_str)) :: out_str
    integer :: i,n
    integer, parameter :: offset = 32
    out_str = in_str

    do i = 1, len(out_str)
       if (out_str(i:i) >= "A" .and. out_str(i:i) <= "Z") then
          out_str(i:i) = achar(iachar(out_str(i:i)) + offset)
       end if
    end do
  end function tolower


  
  function toupper(in_str) result(out_str)
    character(*), intent(in) :: in_str
    character(len(in_str)) :: out_str
    integer :: i,n
    integer, parameter :: offset = 32
    out_str = in_str

    do i = 1, len(out_str)
       if (out_str(i:i) >= "a" .and. out_str(i:i) <= "z") then
          out_str(i:i) = achar(iachar(out_str(i:i)) - offset)
       end if
    end do
  end function toupper


  
  function findspace(string)
    integer :: findspace
    character(len=*), intent(in) :: string
    integer :: i, start = 1

    i = start
    do while(string(i:i) /= ' ' .and. i < len(string))
       i = i + 1
    end do
    if (i >= len(string)) i = -1
    findspace = i    
  end function findspace
end module letkf_obs
