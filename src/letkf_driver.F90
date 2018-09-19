!================================================================================
!> Basic driver program that calls the LETKF initialize and run routines
!! with the default classes.
!================================================================================
PROGRAM letkf_driver
  USE letkf
  IMPLICIT NONE

  INTEGER :: i
  CHARACTER(len=:),ALLOCATABLE :: config_filename

  ! Get the command line argument pointing to the namelist to use
  i=command_argument_COUNT()
  IF (i == 1) THEN
     ALLOCATE(CHARACTER(1024) :: config_filename)
     CALL get_command_ARGUMENT(1, config_filename)
     config_filename=TRIM(config_filename)
  ELSE
     config_filename="letkf.json"
  END IF

  ! initialize and run
  CALL letkf_init(config_filename)
  CALL letkf_run()

END PROGRAM letkf_driver
