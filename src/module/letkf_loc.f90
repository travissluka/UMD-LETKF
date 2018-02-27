module letkf_loc
  !! observation localization routines
  !! TODO:
  !!  * still need to generalize the horizontal / vertical / variable localization routines
  
  use letkf_mpi
  use letkf_obs
  use letkf_state

  implicit none
  private

  
  ! public module methods
  !------------------------------------------------------------
  public :: letkf_loc_init
  public :: letkf_loc_register
  public :: loc_gc, loc_gaus

  !! TODO: this is temporary for the ocean development, need to make horizontal localization specification generalized
  real, private :: loc_hz(2) = (/720e3, 200e3/)

  !================================================================================
  !================================================================================

  
  type, public ::  localizer_group
     !! contains the definition for a single localization group,
     !! consisting of the list of model state slabs,
     !! and the localization parameters
     integer, allocatable :: slab(:)
     real :: hz_loc = -1
     real :: vt_loc = -1
  end type localizer_group

  !--------------------------------------------------------------------------------
  
  type, public, abstract :: localizer
     !! the abstract class that all localization specification
     !! classes are to be derived from. Responsible for
     !! determining the logic of how LETKF should perform
     !! horizontal / vertical / variable localization
   contains

     procedure(I_letkf_loc_getstr),     deferred :: get_name
     procedure(I_letkf_loc_getstr),     deferred :: get_desc
     procedure(I_letkf_loc_getgroups),  deferred :: get_groups
     procedure(I_letkf_loc_maxhz),      deferred :: get_maxhz
     procedure(I_letkf_loc_localize),   deferred :: localize

  end type localizer

  !--------------------------------------------------------------------------------  
  
  abstract interface

     function I_letkf_loc_getstr(self)
       import localizer
       class(localizer) :: self
       character(:), allocatable :: I_letkf_loc_getstr
     end function I_letkf_loc_getstr
     
     pure subroutine I_letkf_loc_getgroups(self, ij, grps)
       import localizer
       import localizer_group
       class(localizer), intent(in) :: self
       integer,          intent(in) :: ij
       type(localizer_group), intent(out),  allocatable :: grps(:)
     end subroutine I_letkf_loc_getgroups

     pure function I_letkf_loc_maxhz(self, ij)
       import localizer
       class(localizer), intent(in) :: self
       integer,          intent(in) :: ij
       real :: I_letkf_loc_maxhz
     end function I_letkf_loc_maxhz
     
     pure subroutine I_letkf_loc_localize(self, ij, lg, ob_num, ob_dist, ob_idx, rloc)
       import localizer
       import localizer_group
       class(localizer),       intent(in) :: self
       integer,                intent(in) :: ij
       class(localizer_group), intent(in) :: lg
       integer,                intent(in) :: ob_num
       real,                   intent(in) :: ob_dist(:)
       integer,                intent(in) :: ob_idx(:)
       real,                   intent(out):: rloc(:)
     end subroutine I_letkf_loc_localize

  end interface

  !--------------------------------------------------------------------------------
  
  type localizer_ptr
     class(localizer), pointer :: p
  end type localizer_ptr

  
  !================================================================================
  !================================================================================


  ! registration of built-in and user-defined localization classes
  integer, parameter        :: localizer_reg_max = 100
  integer                   :: localizer_reg_num = 0
  type(localizer_ptr)       :: localizer_reg(localizer_reg_max)
  class(localizer),public, protected, pointer :: localizer_class



  !================================================================================
  !================================================================================
  ! built-in localizer classes

  type, public, extends(localizer) :: localizer_novrt
   contains
     procedure :: get_name   => localizer_novrt_name
     procedure :: get_desc   => localizer_novrt_desc
     procedure :: get_groups => localizer_novrt_groups
     procedure :: get_maxhz  => localizer_novrt_maxhz
     procedure :: localize   => localizer_novrt_localize     
  end type localizer_novrt


  

contains



  
  !================================================================================
  !================================================================================
  subroutine letkf_loc_init(nml_filename)
    character(len=*), intent(in) :: nml_filename

    character(len=:), allocatable :: locclass
    integer :: unit, i
    
    namelist /letkf_loc/ locclass, loc_hz
    
    if(pe_isroot)then
       print "(A)", ""
       print "(A)", ""
       print "(A)", '============================================================'
       print "(A)", ' letkf_loc_init() : '
       print "(A)", '============================================================'
    end if

    ! read in our section of the namelist
    allocate(character(1024) :: locclass)
    open(newunit=unit, file=nml_filename)
    read(unit, nml=letkf_loc)
    close(unit)
    locclass = trim(locclass)
    if (pe_isroot) then
       print letkf_loc
    end if

    ! print a list of all loc classes that have been registered
    if(pe_isroot) then
       print *, ""
       print *, "List of localization classes registered:"
       do i =1, localizer_reg_num
          print "(A,A,3A)", " * ", toupper(localizer_reg(i)%p%get_name()), &
               " (",localizer_reg(i)%p%get_desc(), ")"
       end do
       print *, ""
    end if

    ! determine the localizer class to use
    nullify (localizer_class)
    do i =1,localizer_reg_num
       if(trim(toupper(localizer_reg(i)%p%get_name())) == trim(toupper(locclass))) then
          localizer_class => localizer_reg(i)%p
          exit
       end if
    end do    
    if (.not. associated(localizer_class)) then
       if(pe_isroot) &
            print *, 'ERROR: localizer class "', toupper(trim(locclass)), &
            '" not found. Check that the name is in the list of registered classes.'
       stop 1
    end if
    if (pe_isroot) print *, 'Using "', trim(localizer_class%get_name()), '"'
    
  end subroutine letkf_loc_init
  !================================================================================


  
  
  !================================================================================
  !================================================================================
  pure function loc_gc(z, L)
    !! Gaspari-Cohn localization function
    !! Possibly faster than the Gaussian function, depending on computer architecture.
    !! Similar shape to Gaussian, except it is compact, goes to 0 at (\ 2L sqrt( 0.3) \ )

    real, intent(in) :: z
    real, intent(in) :: L
    !! (\ e^(0.5) \)
    real :: loc_gc
    real :: c
    real :: abs_z, z_c

    c = L / sqrt(0.3)
    abs_z = abs(z)
    z_c = abs_z / c

    if (abs_z >= 2*c) then
       loc_gc = 0.0
    elseif (abs_z < 2*c .and. abs_z > c) then
       loc_gc = &
            (1.0/12.0)*z_c**5 - 0.5*z_c**4 + &
            (5.0/8.0)*z_c**3 + (5.0/3.0)*z_c**2 &
            - 5.0*z_c + 4 - (2.0/3.0)*c/abs_z
    else
       loc_gc = &
            -0.25*z_c**5 + 0.5*z_c**4 + &
            (5.0/8.0)*z_c**3 - (5.0/3.0)*z_c**2 + 1
    end if
  end function loc_gc
  !================================================================================




  !================================================================================
  !================================================================================
  pure function loc_gaus(z, L)
    real, intent(in) :: z, L
    real :: loc_gaus
    loc_gaus = exp( -0.5 *  z*z / L*L)
  end function loc_gaus
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_loc_register(cls)
    class(localizer), pointer :: cls
    integer :: i

    if(localizer_reg_num == localizer_reg_max) then
       print *, "ERROR: too many localizer classes registered"
       stop 1
    end if

    do i=1, localizer_reg_num
       if(toupper(localizer_reg(i)%p%get_name()) == toupper(cls%get_name())) then
          print *, "ERROR: can't register localizer class '", toupper(cls%get_name()), &
               "', a class by that name has already been registered"
          stop 1
       end if
    end do

    localizer_reg_num = localizer_reg_num + 1
    localizer_reg(localizer_reg_num)%p => cls
  end subroutine letkf_loc_register
  !================================================================================




  !================================================================================
  !================================================================================ 
  function toupper(in_str) result(out_str)
    character(*), intent(in) :: in_str
    character(len(in_str)) :: out_str
    integer :: i
    integer, parameter :: offset = 32

    out_str = in_str
    do i = 1, len(out_str)
       if (out_str(i:i) >= "a" .and. out_str(i:i) <= "z") then
          out_str(i:i) = achar(iachar(out_str(i:i)) - offset)
       end if
    end do
  end function toupper
  !================================================================================




  !================================================================================
  !================================================================================
  ! no vertical localization, the simlest case.
  ! only horizontal localization needs to be specified
  
  function localizer_novrt_name(self)
    class(localizer_novrt) :: self
    character(:), allocatable :: localizer_novrt_name
    localizer_novrt_name  = "LOC_NOVT"
  end function localizer_novrt_name

  
  function localizer_novrt_desc(self)
    class(localizer_novrt) :: self
    character(:), allocatable :: localizer_novrt_desc
    localizer_novrt_desc = "Single column horizontal localization only scheme"
  end function localizer_novrt_desc

  ! TODO: not a clean way of doing this, make it so that i'm not returning
  ! an allocatable array.
  pure subroutine localizer_novrt_groups(self, ij, grps)
    class(localizer_novrt), intent(in) :: self
    integer, intent(in) :: ij
    type(localizer_group), intent(out), allocatable :: grps(:)
    integer:: i

    if (allocated(grps)) then
       do i=1,size(grps)
          deallocate(grps(i)%slab)
       end do
       deallocate(grps)
    end if

    allocate(grps(1))
    allocate(grps(1)%slab(grid_ns))
    do i=1,grid_ns
       grps(1)%slab(i) = i
    end do
  end subroutine localizer_novrt_groups

  
  pure function localizer_novrt_maxhz(self, ij) result(d)
    class(localizer_novrt), intent(in) :: self
    integer,                intent(in) :: ij
    real :: d, r
    r=abs(lat_ij(ij))/90.0
    d=r*loc_hz(1) + (1.0-r)*loc_hz(0)
  end function localizer_novrt_maxhz

  
  pure subroutine localizer_novrt_localize(self, ij, lg, ob_num, ob_dist, ob_idx, rloc)
    class(localizer_novrt), intent(in) :: self
    integer,                intent(in) :: ij
    class(localizer_group), intent(in) :: lg
    integer,                intent(in) :: ob_num
    real,                   intent(in) :: ob_dist(:)
    integer,                intent(in) :: ob_idx(:)
    real,                   intent(out):: rloc(:)

    integer :: i
    real :: r,l

    do i = 1, ob_num
       ! Horizontal localization 
       !TODO: hardcoded to lat based hz loc distance, generalize this
       r=abs(lat_ij(ij))/90.0
       l=r*loc_hz(1) + (1.0-r)*loc_hz(0)
       rloc(i) = loc_gc(ob_dist(i), l)

    end do

  end subroutine localizer_novrt_localize


  !================================================================================
  !================================================================================

  
end module letkf_loc
