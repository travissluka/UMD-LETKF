! Copyright 2018-2019 Travis Sluka
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

PROGRAM check_norm_nc
  USE netcdf
  IMPLICIT NONE

  CHARACTER(len=1024) :: file1, file2, str
  REAL :: max_norm

  INTEGER :: i, j, s(4)
  INTEGER :: nc1, nc2
  INTEGER :: ndim1, ndim2, nvar1, nvar2
  INTEGER, ALLOCATABLE :: dimlen1(:), dimlen2(:)
  INTEGER :: dimids1(NF90_MAX_VAR_DIMS), dimids2(NF90_MAX_VAR_DIMS)
  REAL, ALLOCATABLE :: tmpr1(:,:,:,:), tmpr2(:,:,:,:)
  REAL :: n


  CALL get_command_ARGUMENT(1, file1)
  CALL get_command_ARGUMENT(2, file2)
  CALL get_command_ARGUMENT(3, str)
  READ (str, *) max_norm


  CALL check(nf90_open(file1, nf90_nowrite, nc1))
  CALL check(nf90_open(file2, nf90_nowrite, nc2))

  ! check the dimensions
  CALL check(nf90_inquire(nc1, ndim1, nvar1))
  CALL check(nf90_inquire(nc2, ndim2, nvar2))
  IF(ndim1 /= ndim2) THEN
     PRINT *, "ndim mismatch"
     STOP 1
  END IF
  ALLOCATE( dimlen1(ndim1), dimlen2(ndim2) )
  DO i =1, ndim1
     ! check lengths
     CALL check(nf90_inquire_dimension(nc1, i, len=dimlen1(i)))
     CALL check(nf90_inquire_dimension(nc2, i, len=dimlen2(i)))
     IF(dimlen1(i) /= dimlen2(i)) THEN
        PRINT *, "dim length mismatch"
        STOP 1
     END IF
     ! TODO: check name
     ! TODO: check attributes
  END DO


  ! check the variables
  IF(nvar1 /= nvar2) THEN
     PRINT *, "nvar mismatch"
     STOP 1
  END IF
  DO i = 1, nvar1
     ! dimension lengths
     CALL check(nf90_inquire_variable(nc1, i, ndims=ndim1))
     CALL check(nf90_inquire_variable(nc2, i, ndims=ndim2))
     IF(ndim1 /= ndim2) THEN
        PRINT *, "dimension length of variable do not match"
        STOP 1
     END IF

     ! dimension ids
     CALL check(nf90_inquire_variable(nc1, i, dimids=dimids1))
     CALL check(nf90_inquire_variable(nc2, i, dimids=dimids2))

     ! variable size
     s=1
     DO j=1,ndim1
        s(j) = dimlen1(dimids1(j))
     END DO

     ! variable value
     ALLOCATE(tmpr1(s(1),s(2),s(3),s(4)))
     ALLOCATE(tmpr2(s(1),s(2),s(3),s(4)))
     CALL check(nf90_get_var(nc1,i, tmpr1))
     CALL check(nf90_get_var(nc2,i, tmpr2))
     n = SQRT(SUM((tmpr1-tmpr2)*(tmpr1-tmpr2)))
     !n = MAXVAL((tmpr1-tmpr2)*(tmpr1-tmpr2))
     IF(n>max_norm) THEN
        PRINT *, "max difference of ",n
        STOP 1
     END IF
     DEALLOCATE(tmpr1, tmpr2)
  END DO


CONTAINS


  SUBROUTINE check(status)
    INTEGER, INTENT(in) :: status

    IF(status /= nf90_noerr) THEN
       WRITE (*,*) TRIM(nf90_strerror(status))
       STOP 100
    END IF
  END SUBROUTINE check

END PROGRAM check_norm_nc
