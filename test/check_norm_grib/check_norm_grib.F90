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

PROGRAM check_norm_grib
  USE wgrib2api
  IMPLICIT NONE

  CHARACTER(len=1024) :: file1, file2, str1, str2, str3
  REAL :: max_norm

  INTEGER(kind=8), PARAMETER :: BUFFER_SIZE = 1024*1024*1
  CHARACTER(len=BUFFER_SIZE) :: buffer1, buffer2
  REAL, ALLOCATABLE :: data1(:,:), data2(:,:)

  REAL :: n
  INTEGER :: ret, ret2, nrec, i

  CALL get_command_ARGUMENT(1, file1)
  CALL get_command_ARGUMENT(2, file2)
  CALL get_command_ARGUMENT(3, str1)
  READ (str1, *) max_norm

  ret = wgrib2_set_mem_buffer(buffer1, BUFFER_SIZE, 1)
  ret = wgrib2_set_mem_buffer(buffer2, BUFFER_SIZE, 2)
  ret = grb2_mk_inv(file1, '@mem:1')
  ret = grb2_mk_inv(file2, '@mem:2')

  ! make sure the two files have the same number of records
  nrec = grb2_inq(file1, '@mem:1', '.*', regex=1)
  i = grb2_inq(file2, '@mem:2', '.*', regex=1)
  IF (i/=nrec) THEN
     PRINT *, "ERROR: different number of records"
     STOP 1
  END IF

  ! check each record
  DO i = 1, nrec
     WRITE (str1, '(A,I0,A)' ) '^', i, '\:.*'
     ret  = grb2_inq(file1, '@mem:1', str1, regex=1, desc=str2, data2=data1)
     ret2 = grb2_inq(file2, '@mem:2', str1, regex=1, desc=str3, data2=data2)
     IF(ret/=ret2) STOP 1
     IF(TRIM(str2)/=TRIM(str3)) THEN
        PRINT *, i, "record description is different"
        STOP 1
     END IF
     n = SQRT(SUM( (data1-data2)*(data1-data2)))
     IF(n>max_norm) THEN
        PRINT *, "max difference of ",n
        STOP 1
     END IF
  END DO
END PROGRAM check_norm_grib
