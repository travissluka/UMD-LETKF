module letkf_state_types
  use letkf_common


  implicit none
  private

  public :: statedef
  public :: stateio

  type statedef
     character(len=:), allocatable :: name_short
     character(len=:), allocatable :: name_long
     character(len=:), allocatable :: file_field
     character(len=:), allocatable :: units
     integer                       :: levels
!   contains
!     procedure :: print => statedef_print
  end type statedef



  type, abstract :: stateio
     !! abstract base class for reading and writing of state files
   contains
     procedure(I_stateio_init),  deferred :: init
     procedure(I_stateio_read),  deferred :: read
     procedure(I_stateio_write), deferred :: write
  end type stateio

  abstract interface

     subroutine I_stateio_init(self)
       import stateio
       class(stateio) :: self
     end subroutine I_stateio_init

     subroutine I_stateio_read(self, filename, state)
       import stateio
       class(stateio) :: self
       character(len=*),  intent(in)  :: filename
       real, allocatable, intent(out) :: state(:,:,:)
     end subroutine I_stateio_read

     subroutine I_stateio_write(self, filename, state)
       import stateio
       class(stateio) :: self
       character(len=*),  intent(in)  :: filename
       real, allocatable, intent(in)  :: state(:,:,:)
     end subroutine I_stateio_write

  end interface


end module letkf_state_types
