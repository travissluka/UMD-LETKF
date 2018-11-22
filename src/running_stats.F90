!================================================================================
!>
!================================================================================
MODULE running_stats_mod
  IMPLICIT NONE
  PRIVATE


  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC :: running_stats
     INTEGER :: count = 0
     REAL :: vmin(2) = HUGE(0.0)
     REAL :: vmax(2) = -HUGE(0.0)
     REAL :: M(2) = 0.0
     REAL :: S = 0.0

   CONTAINS

     PROCEDURE :: add  => running_stats_add
     PROCEDURE :: min  => running_stats_min
     PROCEDURE :: max  => running_stats_max
     PROCEDURE :: mean => running_stats_mean
     PROCEDURE :: variance => running_stats_variance
  END TYPE running_stats
  !================================================================================



CONTAINS



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE running_stats_add(self, val)
    CLASS(running_stats) :: self
    REAL, INTENT(in) :: val

    REAL :: d(2)

    self%count = self%count + 1

    self%vmin(1) = MIN(self%vmin(1), val)
    self%vmax(1) = MAX(self%vmax(1), val)
    self%vmin(2) = MIN(self%vmin(2), val*val)
    self%vmax(2) = MAX(self%vmax(2), val*val)

    d(1) = val-self%M(1)
    d(2) = (val*val)-self%M(2)

    self%M(1) = self%M(1) + d(1)/self%count
    self%M(2) = self%M(2) + d(2)/self%count

    self%S  = self%S + (val-self%M(1))*d(1)
  END SUBROUTINE running_stats_add
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION running_stats_min(self, m) RESULT(val)
    CLASS(running_stats) :: self
    INTEGER, OPTIONAL :: m
    INTEGER :: m0
    REAL ::val

    m0 = MERGE(m, 1, PRESENT(m))
    val = self%vmin(m0)
  END FUNCTION running_stats_min
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION running_stats_max(self, m) RESULT(val)
    CLASS(running_stats) :: self
    INTEGER, OPTIONAL :: m
    INTEGER :: m0
    REAL ::val

    m0 = MERGE(m, 1, PRESENT(m))
    val = self%vmax(m0)
  END FUNCTION running_stats_max
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION running_stats_mean(self, m) RESULT(val)
    CLASS(running_stats) :: self
    REAL ::val
    INTEGER, OPTIONAL :: m
    INTEGER :: m0

    m0 = MERGE(m, 1, PRESENT(m))
    val = self%M(m0)
  END FUNCTION running_stats_mean
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION running_stats_variance(self) RESULT(val)
    CLASS(running_stats) :: self
    REAL ::val
    IF (self%count < 2) THEN
       val = 0.0
    ELSE
       val = self%S/(self%count - 1)
    END IF
  END FUNCTION running_stats_variance
  !================================================================================

END MODULE running_stats_mod
