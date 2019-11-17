! Copyright 2019-2019 Travis Sluka
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


MODULE letkf_util

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: str_tolower


CONTAINS

  !================================================================================
  !> Convert a string to lowercase
  !--------------------------------------------------------------------------------
  FUNCTION str_tolower(in_str) RESULT(out_str)
    CHARACTER(*), INTENT(in) :: in_str !< input string
    CHARACTER(LEN(in_str)) :: out_str  !< output string

    INTEGER :: i
    INTEGER, PARAMETER :: offset = 32

    out_str = in_str
    DO i = 1, LEN(out_str)
       IF (out_str(i:i) >= "A" .AND. out_str(i:i) <= "Z") THEN
          out_str(i:i) = ACHAR(IACHAR(out_str(i:i)) + offset)
       END IF
    END DO

  END FUNCTION str_tolower
  !================================================================================

END MODULE letkf_util
