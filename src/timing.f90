module timing
  implicit none
  private
  
  public :: timing_start, timing_end
  public :: timing_print


  
  type timer_obj
     integer(kind=8):: start, end
     character(len=1024) :: name
  end type timer_obj

  integer, parameter :: max_timers = 1024
  integer, save :: active_timers = 0
  type(timer_obj) :: timer_objs(max_timers)
  
contains

  function gettimer(timer) result(id)
    character(len=*) :: timer
    integer :: id
    integer :: i

    if (active_timers == 0) then
       active_timers = 1
       timer_objs(1)%name='all'
       call system_clock(timer_objs(1)%start)
    end if
    
    ! see if the timer already exists
    i = 1
    do while (i <= active_timers)
       if (trim(timer_objs(i)%name) == trim(timer)) exit
       i = i + 1
    end do

    ! if not, create one
    if ( i > active_timers) then
       active_timers = active_timers + 1
       i = active_timers
       if (active_timers > max_timers) then
          print *, "ERROR, too many timers have been created, increase max_timers"
          stop 1
       end if
       timer_objs(i)%name = timer
    end if

    id  = i
  end function gettimer

  
  subroutine timing_start(timer)
    character(len=*), intent(in) :: timer
    type(timer_obj) :: t
    integer :: id    
    id = gettimer(timer)
    call system_clock(timer_objs(id)%start)
  end subroutine timing_start


  subroutine timing_end(timer)
    character(len=*) :: timer
    integer :: i
    integer :: id
    id = gettimer(timer)
    call system_clock(timer_objs(id)%end)
  end subroutine timing_end


  subroutine timing_print
    integer :: i
    real :: t

    integer(kind=8) :: timer_rate
    real :: total

    print *,""    
    print *,"Timing:"
    print *,"============================================================"    
    print *,""
    

    call system_clock(timer_objs(1)%end, timer_rate)
    total = (timer_objs(1)%end-timer_objs(1)%start)/timer_rate
    
    i = 2
    do while (i <= active_timers)
       t = (timer_objs(i)%end-timer_objs(i)%start)/timer_rate       
       print '(A15,F10.3,F5.1,A1)', trim(timer_objs(i)%name), t, t/total*100,"%"
       i = i +1
    end do
    print *,"      -------------------------"
    print '(A15,F10.3)', "total", total
  end subroutine timing_print
  
end module timing
