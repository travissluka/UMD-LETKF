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
