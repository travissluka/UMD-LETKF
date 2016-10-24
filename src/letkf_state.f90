module letkf_state
  !! performs I/O for the background and analysis states
  use global
  use letkf_mpi
  use str
  implicit none


  integer :: grid_x
  integer :: grid_y
  integer :: grid_z
  !! number of vertical levels for the 3D state variables
  integer :: grid_3d = 5
  !! number of 3D state variables
  integer :: grid_2d = 1
  !! number of 2D state variables


  type statedef
     character(len=:), allocatable :: name_short
     character(len=:), allocatable :: name_long     
     character(len=:), allocatable :: units
     integer :: levels
   contains
     procedure :: print => statedef_print
  end type statedef


  type(statedef), allocatable :: statedef_list(:)

contains

  subroutine letkf_state_init()
    character(len=:), allocatable :: statedef_file
    integer :: unit

    namelist /grid_def/  grid_x, grid_y, grid_z, statedef_file

    if (pe_isroot) then
       print *, ''
       print *, 'model state definition configuration file'
       print *, '------------------------------------------------------------'
    end if

    ! read in namelist
    allocate(character(len=1024) :: statedef_file)
    open(newunit=unit, file='letkf.nml')
    read(unit, nml=grid_def)
    close(unit)

    statedef_file = trim(statedef_file)
    if (pe_isroot) then
       print grid_def
       print *,''
    end if
    
    ! do some basic checks    
    if (pe_isroot) then
       if (grid_x <= 0) stop 1
       if (grid_y <= 0) stop 1
       if (grid_z <= 0) stop 1
    end if
          
    

    ! read in statedef config file
    call statedef_read(statedef_file)


  end subroutine letkf_state_init



  
  subroutine statedef_read(file)
    character(len=*), intent(in) :: file

    logical :: ex
    integer :: unit, iostat, i, pos
    character(len=1024) :: line, s_tmp
    type(statedef) :: newdef
    integer, parameter :: MAX_STATEDEF = 1024
    type(statedef) :: statedef_list_tmp(MAX_STATEDEF)
    integer :: statedef_list_tmp_len

    if (pe_isroot) then
       print *, 'Reading file "',trim(file),'" ...'
    end if

    ! make sure the file exists
    inquire(file=file, exist=ex)
    if (.not. ex) then
       print *, 'ERROR: file does not exist: "',trim(file),'"'
       stop 1
    end if

    statedef_list_tmp_len = 0
    
    ! open it up for reading
    open(newunit=unit, file=file, action='read')
    do while (.true.)
       ! read a new line
       read(unit, '(A)', iostat=iostat) line
       if (iostat < 0) exit
       if (iostat > 0) then
          print *, 'ERROR: there was a problem reading "', &
               trim(file), '" error code: ',iostat
          stop 1
       end if

       ! convert tabs to spaces
       do i=1, len(line)
          if (line(i:i) == char(9)) line(i:i) = ' '
       end do

       !ignore comments
       s_tmp = adjustl(line)
       if (s_tmp(1:1) == '#') cycle

       !ignore empty lines
       if (len(trim(adjustl(line))) == 0) cycle

       line = adjustl(line)
       pos = findspace(line)
       newdef%name_short = toupper(trim(line(:pos-1)))

       line= adjustl(line(pos+1:))
       pos = findspace(line)      
       s_tmp = tolower(trim(line(:pos-1)))
       if(s_tmp == 'nlev') then
          newdef%levels = grid_z
       else
          read(s_tmp,*) newdef%levels
       end if

       line= adjustl(line(pos+1:))
       pos = findspace(line)      
       newdef%units = toupper(trim(line(:pos-1)))

       line = adjustl(line(pos+1:))
       newdef%name_long = trim(line)

       ! add this one to the list
       statedef_list_tmp_len = statedef_list_tmp_len + 1
       if (statedef_list_tmp_len > MAX_STATEDEF) then
          print *, "ERROR: more states than expecting. check the ",&
               trim(file),'" file, and/or increase MAX_STATEDEF.'
          stop 1
       end if
       statedef_list_tmp(statedef_list_tmp_len) = newdef          
    end do

    close (unit)
    allocate(statedef_list(statedef_list_tmp_len))
    statedef_list = statedef_list_tmp(1:statedef_list_tmp_len)


    ! write a summary
    if (pe_isroot) then
       print *, 'state variables = ', size(statedef_list)
       do pos=1, size(statedef_list)
          call statedef_list(pos)%print()
       end do
    end if

    ! check for duplicate name
    ! TODO

    ! all done
    if(pe_isroot) print *, ""
  end subroutine statedef_read


  subroutine statedef_print(st)
    class(statedef), intent(in) :: st
    print "(A10, I6, A,A)", st%name_short, st%levels, " ",st%name_long
  end subroutine statedef_print
  
  !============================================================
  subroutine letkf_state_read(filename, state_3d, state_2d)
    !! read a complete model state from a file.
    !! Terminates program if the file is not found.

    character(len=*), intent(in) :: filename
    !! name of file from which state variables will be read
    real, intent(inout) :: state_3d(:,:,:,:)
    !! 3D state variables of shape \( (x,y,z,var) \)
    real, intent(inout) :: state_2d(:,:,:)
    !! 2D state variables of shape \( (x,y,var) \)

    integer :: unit, nrec
    integer :: l,k
    real :: r
    logical :: ex
    
    ! check to make sure the file exists
    inquire(file=filename, exist=ex)
    if (.not. ex) then
       print *, "ERROR: file does not exists ",trim(filename)
       stop 1
    end if

    ! read in the file
    open(newunit=unit, file=filename, form='unformatted', access='direct',&
         recl=grid_x*grid_y*sizeof(r))
    nrec = 1
    do l=1,grid_3d
       do k=1,grid_z
          read(unit, rec=nrec) state_3d(:,:,k,l)
          nrec = nrec + 1
       end do
    end do
    do l=1,grid_2d
       read(unit, rec=nrec) state_2d(:,:,l)
    end do
    close(unit)
  end subroutine letkf_state_read

  
  !============================================================
  subroutine letkf_state_write(filename, state_3d, state_2d)
    !! writes a given complete model state out to a file.

    character(len=*), intent(in) :: filename
    !! name of file to which state variables will be written
    real, intent(in) :: state_3d(:,:,:,:)
    !! 3D state variables of shape \( (x,y,z,var) \)
    real, intent(in) :: state_2d(:,:,:)
    !! 2D state variables of shape \( (x,y,var) \)

    integer :: unit, nrec
    real :: r
    integer :: l,k
    open(newunit=unit, file=filename, form='unformatted',&
         access='direct', recl=grid_x*grid_y*sizeof(r))
    nrec = 1
    do l=1,grid_3d
       do k=1, grid_z
          write(unit, rec=nrec) state_3d(:,:,k,l)
          nrec = nrec+1
       end do
    end do
    do l=1,grid_2d
       write(unit, rec=nrec) state_2d(:,:,l)
       nrec = nrec + 1
    end do

    close(unit)
       
  end subroutine letkf_state_write

  
end module letkf_state
