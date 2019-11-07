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
     config_filename="letkf.yaml"
  END IF

  ! initialize and run
  CALL letkf_init(config_filename)
  CALL letkf_run()

END PROGRAM letkf_driver
