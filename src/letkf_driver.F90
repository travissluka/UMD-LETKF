!================================================================================
!> Basic driver program that calls the LETKF initialize and run routines
!! with the default classes.
!================================================================================
program letkf_driver
  use letkf
  implicit none

  integer :: i
  character(len=:),allocatable :: config_filename

  ! Get the command line argument pointing to the namelist to use
  i=command_argument_count()
  if (i == 1) then
    allocate(character(1024) :: config_filename)
    call get_command_argument(1, config_filename)
    config_filename=trim(config_filename)
  else
     config_filename="letkf.json"
  end if

  ! initialize and run
  call letkf_init(config_filename)
  call letkf_run()

end program
