program letkf_driver
  use letkf
  implicit none

  integer :: i
  character(len=1024) ::nml_filename

  ! Get the command line argument pointing to the namelist to use
  i=command_argument_count()
  if (i == 1) then
     call get_command_argument(1, nml_filename)
  else
     nml_filename="namelist.letkf"
  end if

  ! initialize and run
  call letkf_init(nml_filename)
  call letkf_run()
  
end program
