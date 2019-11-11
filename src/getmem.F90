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
!> \cond INTERNAL
!> Basic memory profiler
!--------------------------------------------------------------------------------
MODULE getmem

  USE mpi
  USE iso_fortran_env
#ifdef __INTEL_COMPILER
  USE ifport, ONLY : getpid
#endif

  IMPLICIT NONE
  PRIVATE


  !================================================================================
  !================================================================================
  ! Public module components
  !================================================================================
  !================================================================================

  PUBLIC :: getmem_init
  PUBLIC :: getmem_print



  !================================================================================
  !================================================================================
  ! Private module components
  !================================================================================
  !================================================================================

  INTEGER :: pe_rank
  INTEGER :: pe_root
  INTEGER :: pe_comm
  INTEGER :: pe_size
  LOGICAL :: pe_isroot
  REAL    :: all_min
  REAL    :: all_max
  REAL    :: all_avg



CONTAINS



  !================================================================================
  !> initailize the memory reporting module.
  !--------------------------------------------------------------------------------
  SUBROUTINE getmem_init(p_pe_root, p_pe_comm)
    INTEGER, INTENT(in) :: p_pe_root !< PE number of the root (usually 0)
    INTEGER, INTENT(in) :: p_pe_comm !< the MPI communicator

    INTEGER :: ierr

    pe_root = p_pe_root
    pe_comm = p_pe_comm

    CALL mpi_comm_rank(pe_comm, pe_rank, ierr)
    CALL mpi_comm_size(pe_comm, pe_size, ierr)
    pe_isroot = pe_rank == pe_root

    all_min = 0
    all_max = 0
    all_avg = 0
  END SUBROUTINE getmem_init
  !================================================================================



  !================================================================================
  !> print out the current memory usage statistics
  !--------------------------------------------------------------------------------
  SUBROUTINE getmem_print
    REAL :: cur, max
    REAL :: mem_min, mem_max, mem_avg
    INTEGER :: ierr

    CALL readmem(cur, max)
    CALL mpi_reduce(cur, mem_min, 1, mpi_real, mpi_min, pe_root, pe_comm, ierr)
    CALL mpi_reduce(cur, mem_max, 1, mpi_real, mpi_max, pe_root, pe_comm, ierr)
    CALL mpi_reduce(cur, mem_avg, 1, mpi_real, mpi_sum, pe_root, pe_comm, ierr)
    mem_avg = mem_avg / pe_size
    IF(pe_isroot) THEN
       PRINT *, ""
       PRINT '(A,4F8.2)', "Current Mem in GB (min/max/avg/tot): ",&
            mem_min, mem_max, mem_avg, mem_avg*pe_size
    END IF

    !TODO: there's a way with some linux kernels to clear the previous max memory
    ! try this out? otherwise the process total max is going to be wrong near
    ! the end of the program
    CALL mpi_reduce(max, mem_min, 1, mpi_real, mpi_min, pe_root, pe_comm, ierr)
    CALL mpi_reduce(max, mem_max, 1, mpi_real, mpi_max, pe_root, pe_comm, ierr)
    CALL mpi_reduce(max, mem_avg, 1, mpi_real, mpi_sum, pe_root, pe_comm, ierr)
    mem_avg = mem_avg / pe_size
    IF(pe_isroot) THEN
       PRINT '(A,4F8.2)', "Maximum Mem in GB (min/max/avg/tot): ",&
            mem_min, mem_max, mem_avg, mem_avg * pe_size
    END IF

    IF(mem_min > all_min) all_min = mem_min
    IF(mem_max > all_max) all_max = mem_max
    IF(mem_avg > all_avg) all_avg = mem_avg

  END SUBROUTINE getmem_print
  !================================================================================



  !================================================================================
  !> read the current memory usage for the PE from the linux system
  !--------------------------------------------------------------------------------
  SUBROUTINE readmem(cur, max)
    !! return memory hiwater mark, in gigabytes.
    !! this relies on the VmHWM field set in the /proc/(pid)/status file
    !! if this file or field can't be found, -1 is returned
    REAL, INTENT(out)   :: cur, max
    CHARACTER(len=30)   :: pid_char
    CHARACTER(len=1024) :: proc_file, line
    INTEGER             :: unit, iostat, i
    LOGICAL             :: ex

    ! get process id
    WRITE(pid_char,'(I10)') getpid()
    pid_char = ADJUSTL(pid_char)

    cur = -1
    max = -1

    ! make sure the file exists that contains memory hi usage info in linux
    proc_file="/proc/"//TRIM(pid_char)//"/status"
    INQUIRE(file=proc_file, exist=ex)
    IF (.NOT. ex) THEN
       PRINT *, "can't open ", TRIM(proc_file)
       RETURN
    END IF

    ! read in the VmHWM line of the file
    OPEN(newunit=unit, file=proc_file, action='read')
    iostat = 0
    DO WHILE (iostat==0)
       READ(unit,'(A)',iostat=iostat) line
       i = SCAN(line, ':')
       IF(i<=0) CYCLE
       IF(line(:i) == "VmHWM:") THEN
          ! the VmHWM line has been found, parse out the number of kilobytes
          line = line(i+1:)
          DO i=1,LEN(line)
             IF(line(i:i) == CHAR(9)) line(i:i) = ' '
          END DO
          line = ADJUSTL(line)
          READ(line,*) max
          ! return the number of gigabytes
          max = max/1024/1024
       ELSE  IF(line(:i) == "VmRSS:") THEN
          ! the VmHWM line has been found, parse out the number of kilobytes
          line = line(i+1:)
          DO i=1,LEN(line)
             IF(line(i:i) == CHAR(9)) line(i:i) = ' '
          END DO
          line = ADJUSTL(line)
          READ(line,*) cur
          ! return the number of gigabytes
          cur = cur/1024/1024
       END IF
    END DO
    CLOSE(unit)
  END SUBROUTINE readmem
  !================================================================================



END MODULE getmem
!>\endcond
