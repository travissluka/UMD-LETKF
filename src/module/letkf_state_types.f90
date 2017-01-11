module letkf_state_I
  implicit none
  private

!  public :: statedef
  public :: stateio
  public :: slab


  !! ------------------------------------------------------------


  ! type statedef
  !    character(len=:), allocatable :: name_short
  !    character(len=:), allocatable :: name_long
  !    character(len=:), allocatable :: file_field
  !    character(len=:), allocatable :: units
  !    integer                       :: levels
  ! end type statedef


  !! ------------------------------------------------------------


  type slab
     integer :: lvl
     integer :: var
     real, allocatable :: val(:,:)
  end type slab


  !! ------------------------------------------------------------


  type, abstract :: stateio
     !! abstract base class for reading and writing of state files
     character(len=1024)  :: description
   contains
     procedure(I_stateio_init),   deferred :: init
     procedure(I_stateio_latlon), deferred :: latlon
     procedure(I_stateio_read),   deferred :: read
     procedure(I_stateio_write),  deferred :: write
  end type stateio

  abstract interface
     subroutine I_stateio_init(self,x,y,z)
       import stateio
       class(stateio) :: self
       integer, intent(in) :: x,y,z
     end subroutine I_stateio_init


     subroutine I_stateio_latlon(self, lat, lon)
       import stateio
       class(stateio) :: self
       real, intent(inout) :: lat(:,:), lon(:,:)
     end subroutine I_stateio_latlon

     subroutine I_stateio_read(self, filename, state)
       import stateio
       class(stateio) :: self
       character(len=*),  intent(in)  :: filename
       real, intent(out) :: state(:,:,:)
     end subroutine I_stateio_read

     subroutine I_stateio_write(self, filename, state)
       import stateio
       class(stateio) :: self
       character(len=*),  intent(in)  :: filename
       real, intent(in)  :: state(:,:,:)
     end subroutine I_stateio_write

  end interface


  !! ------------------------------------------------------------

end module letkf_state_I
