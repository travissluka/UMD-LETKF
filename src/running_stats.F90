!================================================================================
!>
!================================================================================
module running_stats_mod
  implicit none
  private

  
  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  type, public :: running_stats
     integer :: count = 0
     real :: vmin(2) = huge(0.0)
     real :: vmax(2) = -huge(0.0)
     real :: M(2) = 0.0
     real :: S = 1.0

   contains

     procedure :: add  => running_stats_add
     procedure :: min  => running_stats_min
     procedure :: max  => running_stats_max
     procedure :: mean => running_stats_mean 
     procedure :: variance => running_stats_variance 
  end type running_stats
  !================================================================================


  
contains

  

  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  subroutine running_stats_add(self, val)
    class(running_stats) :: self
    real, intent(in) :: val

    real :: d(2)

    self%count = self%count + 1
    
    self%vmin(1) = min(self%vmin(1), val)
    self%vmax(1) = max(self%vmax(1), val)
    self%vmin(2) = min(self%vmin(2), val*val)
    self%vmax(2) = max(self%vmax(2), val*val)

    d(1) = val-self%M(1)
    d(2) = (val*val)-self%M(2)
    
    self%M(1) = self%M(1) + d(1)/self%count
    self%M(2) = self%M(2) + d(2)/self%count
    
    self%S  = self%S + (val-self%M(1))*d(1)
  end subroutine running_stats_add
  !================================================================================


  
  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  function running_stats_min(self, m) result(val)
    class(running_stats) :: self
    integer, optional :: m
    integer :: m0
    real ::val

    m0 = merge(m, 1, present(m))
    val = self%vmin(m0)
  end function running_stats_min
  !================================================================================


  
  !================================================================================
  !>
  !--------------------------------------------------------------------------------  
  function running_stats_max(self, m) result(val)
    class(running_stats) :: self
    integer, optional :: m
    integer :: m0
    real ::val

    m0 = merge(m, 1, present(m))
    val = self%vmax(m0)
  end function running_stats_max
  !================================================================================


  
  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  function running_stats_mean(self, m) result(val)
    class(running_stats) :: self
    real ::val
    integer, optional :: m
    integer :: m0

    m0 = merge(m, 1, present(m))
    val = self%M(m0)
  end function running_stats_mean
  !================================================================================


  
  !================================================================================
  !>
  !-------------------------------------------------------------------------------- 
  function running_stats_variance(self) result(val)
    class(running_stats) :: self
    real ::val
    if (self%count < 2) then
       val = 0.0
    else
       val = self%S/(self%count - 1)
    end if
  end function running_stats_variance
  !================================================================================

end module running_stats_mod
