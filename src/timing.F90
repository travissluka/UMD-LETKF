! Copyright 2016-2019 Travis Sluka
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!================================================================================
!>
!================================================================================
MODULE timing
  USE mpi

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: timing_init
  PUBLIC :: timing_print
  PUBLIC :: timing_start
  PUBLIC :: timing_stop

  INTEGER, PUBLIC, PARAMETER :: TIMER_SYNC = 1

  TYPE timer
     CHARACTER(len=20) :: name
     INTEGER(kind=8)   :: total_ticks
     INTEGER(kind=8)   :: start_tick
     INTEGER           :: level

     TYPE(timer), POINTER :: parent
     TYPE(timer), POINTER :: child
     TYPE(timer), POINTER :: sibling
  END TYPE timer


  !================================================================================
  !================================================================================
  ! private module components
  !================================================================================
  !================================================================================


  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  TYPE timer_ptr
     TYPE(timer), POINTER :: p
  END TYPE timer_ptr
  !================================================================================


  INTEGER :: pe_comm, pe_root, pe_rank, pe_size
  TYPE(timer), POINTER :: root_timer
  TYPE(timer), POINTER :: cur_timer



CONTAINS



  !================================================================================
  !> initialize the timing module.
  !! Must be called before any other functions
  !--------------------------------------------------------------------------------
  SUBROUTINE timing_init(comm, root)
    INTEGER, INTENT(in) :: comm
    INTEGER, INTENT(in) :: root
    INTEGER :: ierr

    pe_comm = comm
    pe_root = root

    CALL mpi_comm_rank(pe_comm, pe_rank, ierr)
    CALL mpi_comm_size(pe_comm, pe_size, ierr)

    ALLOCATE(root_timer)
    NULLIFY(root_timer%parent)
    NULLIFY(root_timer%child)
    NULLIFY(root_timer%sibling)
    root_timer%name="Total"
    root_timer%total_ticks = 0
    root_timer%level = 0
    CALL system_CLOCK(root_timer%start_tick)

    cur_timer => root_timer
  END SUBROUTINE timing_init
  !================================================================================



  !================================================================================
  !> Start a timer
  !--------------------------------------------------------------------------------
  SUBROUTINE timing_start(timer_name, flags)
    CHARACTER(*), INTENT(in) :: timer_name
    TYPE(timer), POINTER :: new_timer
    TYPE(timer), POINTER :: ct
    INTEGER, OPTIONAL, INTENT(in) :: flags

    INTEGER :: ierr
    INTEGER :: flags0

    ! set default flags if none given
    flags0 = MERGE(flags, 0, PRESENT(flags))

    ! make sure the timer doesn't already exist as the current timer or any of its parents
    ! TODO

    ! check to see if timer is an already existing child
    NULLIFY(new_timer)
    IF ( ASSOCIATED(cur_timer%child) ) THEN
       ct=>cur_timer%child
       DO WHILE(ASSOCIATED(ct))
          IF (ct%name == timer_name) THEN
             new_timer => ct
             EXIT
          END IF
          ct=> ct%sibling
       END DO
    END IF


    ! add timer as a child
    IF (.NOT. ASSOCIATED(new_timer) ) THEN
       ! create a new timer
       ALLOCATE(new_timer)
       new_timer%parent=>cur_timer
       NULLIFY(new_timer%child)
       new_timer%sibling => cur_timer%child
       new_timer%name=timer_name
       new_timer%total_ticks = 0
       new_timer%level = cur_timer%level + 1

       cur_timer%child => new_timer
    END IF

    ! wait for synchronization, if desired
    IF (IAND(flags0, TIMER_SYNC) > 0) THEN
       CALL mpi_barrier(pe_comm, ierr)
    END IF

    ! set the clock start time
    CALL system_CLOCK(new_timer%start_tick)
    cur_timer => new_timer

  END SUBROUTINE timing_start
  !================================================================================



  !================================================================================
  !> Stop a timer
  !--------------------------------------------------------------------------------
  SUBROUTINE timing_stop(timer_name)
    CHARACTER(*), INTENT(in) :: timer_name
    INTEGER(kind=8) :: ticks


    ! make sure the names match
    IF (cur_timer%name /= timer_name) THEN
       PRINT *, "ERROR: trying to stop non-current timer: ", timer_name
       STOP 1
    END IF

    ! make sure not root
    IF ( ASSOCIATED(cur_timer,root_timer) ) THEN
       PRINT *, "ERROR: trying to stop the main root timer."
       STOP 1
    END IF

    ! stop this timer
    CALL system_CLOCK(ticks)
    cur_timer%total_ticks = cur_timer%total_ticks + (ticks - cur_timer%start_tick)
    cur_timer%start_tick = 0

    ! set current timer to the parent
    cur_timer => cur_timer%parent

  END SUBROUTINE timing_stop
  !================================================================================



  !================================================================================
  !> Print out the across-PE timing statistics
  !! TODO: check for common multiple timers in the tree to sum up
  !--------------------------------------------------------------------------------
  SUBROUTINE timing_print()
    TYPE(timer_ptr) :: stack(100)
    INTEGER :: pos
    TYPE(timer), POINTER :: p, c
    REAL :: t, tmin, tmax, tave, tdif

    INTEGER(kind=8) :: timer_rate, ticks
    INTEGER :: ierr
    CHARACTER(len=500) :: fmt

    CALL system_CLOCK(count_rate=timer_rate)
    IF (pe_root == pe_rank) THEN
       PRINT *, ""
       PRINT *, ""
       PRINT *, "================================================================================"
       PRINT *, " Timing"
       PRINT *, "================================================================================"
       PRINT *, "calculating timer statistics across ", pe_size, " cores..."
       PRINT *, ""
       PRINT '(A18, 4A10)', "","ave","min","max","std"
    END IF

    CALL system_CLOCK(ticks)

    pos = 1
    stack(pos)%p => root_timer
    DO WHILE( pos > 0)
       p=>stack(pos)%p
       pos = pos -1

       ! get mpi stats
       t = REAL(p%total_ticks) / REAL(timer_rate)
       IF ( p%start_tick > 0 ) THEN
          t = t + REAL(ticks - p%start_tick) / REAL(timer_rate)
       END IF
       CALL mpi_reduce   (t, tmin, 1, mpi_real, mpi_min, pe_root, pe_comm, ierr)
       CALL mpi_reduce   (t, tmax, 1, mpi_real, mpi_max, pe_root, pe_comm, ierr)
       CALL mpi_allreduce(t, tave, 1, mpi_real, mpi_sum, pe_comm, ierr)
       tave = tave / pe_size
       t = (t-tave)*(t-tave)
       CALL mpi_reduce   (t, tdif, 1, mpi_real, mpi_sum, pe_root, pe_comm, ierr)
       tdif = SQRT(tdif/pe_size)

       ! print
       IF (pe_root == pe_rank) THEN
          WRITE (fmt, '(A,I0,A,I0,A)')  "(A",2*(p%level+1),",A",18-2*(p%level+1),",4F10.2)"
          PRINT fmt, "",p%name, tave, tmin, tmax, tdif
       END IF

       ! add children to processing stack
       c=>p%child
       DO WHILE(ASSOCIATED(c))
          pos = pos+1
          stack(pos)%p => c
          c => c%sibling
       END DO

    END DO

  END SUBROUTINE timing_print
  !================================================================================


END MODULE timing
