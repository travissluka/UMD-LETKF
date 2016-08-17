module timing
  implicit none

  public :: timing_start, timing_end
  public :: timing_print

  integer :: timers = 0
  
contains

  function gettimer(timer) result(id)
    character(len=*) :: timer
    integer :: id

    id = -1
  end function gettimer

  
  subroutine timing_start(timer)
    character(len=*) :: timer    
  end subroutine timing_start


  subroutine timing_end(timer)
    character(len=*) :: timer    
  end subroutine timing_end


  subroutine timing_print
  end subroutine timing_print
  
end module timing
