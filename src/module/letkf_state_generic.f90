module letkf_state_generic
  use letkf_state_types

  implicit none
  private

  public :: stateio_generic

  type, extends(stateio) :: stateio_generic
   contains
     procedure :: init  => stateio_generic_init
     procedure :: read  => stateio_generic_read
     procedure :: write => stateio_generic_write
  end type stateio_generic



contains


  subroutine stateio_generic_init(self)
    class(stateio_generic) :: self
  end subroutine stateio_generic_init


  subroutine stateio_generic_read(self, filename, state)
    class(stateio_generic) :: self
    character(len=*), intent(in)  :: filename
    real, allocatable, intent(out) :: state(:,:,:)
  end subroutine stateio_generic_read

  subroutine stateio_generic_write(self, filename, state)
    class(stateio_generic) :: self
    character(len=*), intent(in)  :: filename
    real, allocatable, intent(in) :: state(:,:,:)
  end subroutine stateio_generic_write

end module letkf_state_generic
