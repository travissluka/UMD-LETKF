module letkf_state_nc
  ! TODO, make an actual generic one... this is for the MOM4 grid
  use letkf_mpi ! TODO, remove this
  use letkf_state
  use netcdf

  implicit none
  private

  public :: stateio_nc

  type, extends(stateio) :: stateio_nc
   contains
     procedure :: get_name => stateio_get_name
     procedure :: get_desc => stateio_get_desc
     procedure :: init     => stateio_init
     procedure :: read     => stateio_read
     procedure :: write    => stateio_write
!     procedure :: latlon   => stateio_latlon
!     procedure :: mask     => stateio_mask
  end type stateio_nc

!  integer, parameter :: grid_nz = 40

  real, allocatable :: nom_lon(:), nom_lat(:), depths(:)


  
  character(len=*), parameter :: statedef_file = "statedef.cfg"


  type data_entry
     character(len=30)  :: w1
     character(len=30)  :: w2
     character(len=512) :: w3
     integer :: i
  end type data_entry

  integer, parameter :: max_vars = 100
  integer            :: num_vars = 0
  type(data_entry)   :: data_entries(max_vars)


  
contains



  !================================================================================
  !================================================================================
  function stateio_get_name(self) result(str)
    class(stateio_nc) :: self
    character(:), allocatable :: str
    str = "LETKF_NC"
    self%i = self%i
  end function stateio_get_name
  !================================================================================




  !================================================================================
  !================================================================================
  function stateio_get_desc(self) result(str)
    class(stateio_nc) :: self
    character(:), allocatable :: str
    str = "generic netCDF state I/O"
    self%i = self%i
  end function stateio_get_desc
  !================================================================================    




  !================================================================================
  !================================================================================
  function getline(str) result(line)
    character(len=*), intent(in) :: str
    integer :: line, i
    line = -1
    do i=1,num_vars       
       if(toupper(trim(str)) == toupper(trim(data_entries(i)%w1))) then
          line = i
          return
       end if
    end do
  end function getline
  !================================================================================



  !================================================================================
  !================================================================================
  subroutine stateio_init(self, lat, lon, mask)
    class(stateio_nc)   :: self
    real, intent(inout) :: lat(:,:)
    real, intent(inout) :: lon(:,:)
    logical, intent(inout) :: mask(:,:)

    real :: tmp(grid_nx, grid_ny)

    integer :: ncid, varid
    integer :: i, unit, iostat, i1,i2,i3
    logical :: ex
    character(len=1024) :: line
    type(data_entry) :: newline  
     
    ! read in the statedef configuration file
    ! ------------------------------------------------------------
    if(pe_isroot) then
       print *, ''
       print *, ""
       print *, "state_generic initialization"
       print *, "------------------------------------------------------------"
       print *, "Reading file '", trim(statedef_file),"' ..."
    end if

    ! make sure the file exists
    inquire(file=statedef_file, exist=ex)
    if (.not. ex) then
       print *, 'ERROR: file does not exist: "', trim(statedef_file), '"'
       stop 1
    end if

    ! open it up for reading
    open(newunit=unit, file=statedef_file, action='read')
    do while(.true.)
       ! read a new line
       read(unit, '(A)', iostat=iostat) line
       if (iostat < 0) exit
       if (iostat > 0) then
          print *, 'ERROR: there was a problem reading the file.'
          print *, '  error code: ',iostat
          stop 1
       end if

       ! convert tabs to spaces
       do i=1, len(line)
          if(line(i:i) == char(9)) line(i:i) = ' '
       end do

       ! ignore comments / empty lines
       line = adjustl(line)
       if(line(1:1) == '#') cycle
       if(len(trim(adjustl(line))) == 0) cycle

       ! read in the line
       read(line, *, iostat=iostat) newline
       newline%w1 = toupper(newline%w1)
       if(num_vars == max_vars) then
          print *, "ERROR: too many lines in ",trim(statedef_file)
          stop 1
       end if
       if(getline(newline%w1) >0) then
          print *, "ERROR: duplicate entry '",trim(newline%w1),"'"
          stop 1
       end if
       num_vars = num_vars + 1
       data_entries(num_vars) = newline
    end do

    ! ------------------------------------------------------------
    ! read in grid definitions
    if(pe_isroot) then
       print *, ""
       print *, "grid defnition..."
    end if
    ! TODO clean this up
    ! TODO give choice between 1d or 2d
    ! TODO allow specification in cfg file of type (real, logical, int) and inversion (for mask)

    ! read in the lat
    i1 = getline('grid_lat_1d')
    i2 = getline('grid_lat_2d')
    if( .not. xor(i1 > 0, i2 >0) ) then
       print *, "ERROR: one, and only one of the following must be specified: ",&
            "'grid_lat_2d' 'grid_lat_1d'"
       stop 1
    end if
    if(pe_isroot) then
       print *, 'reading ',trim(data_entries(i2)%w1),' from "', trim(data_entries(i2)%w2),&
            '" of file "',trim(data_entries(i2)%w3),'"'
    end if

    call check(nf90_open(data_entries(i2)%w3, nf90_nowrite, ncid),&
         "File: "//trim(data_entries(i2)%w3))
    call check(nf90_inq_varid(ncid, data_entries(i2)%w2, varid),&
         "File: "//trim(data_entries(i2)%w3)//" Var: "//trim(data_entries(i2)%w2))
    call check(nf90_get_var(ncid, varid, lat),&
         "File: "//trim(data_entries(i2)%w3)//" Var: "//trim(data_entries(i2)%w2))
    call check(nf90_close(ncid),&
         "File: "//trim(data_entries(i2)%w3))


    ! read in the lon
    i1 = getline('grid_lon_1d')
    i2 = getline('grid_lon_2d')
    if( .not. xor(i1 > 0, i2 >0) ) then
       print *, "ERROR: one, and only one of the following must be specified: ",&
            "'grid_lon_2d' 'grid_lon_1d'"
       stop 1
    end if
    if(pe_isroot) then
       print *, 'reading ',trim(data_entries(i2)%w1),' from "', trim(data_entries(i2)%w2),&
            '" of file "',trim(data_entries(i2)%w3),'"'
    end if
    call check(nf90_open(data_entries(i2)%w3, nf90_nowrite, ncid))
    call check(nf90_inq_varid(ncid, data_entries(i2)%w2, varid))
    call check(nf90_get_var(ncid, varid, lon))
    call check(nf90_close(ncid))

    ! read in the nom lat 
    allocate(nom_lon(grid_nx))
    allocate(nom_lat(grid_ny))
    allocate(depths(grid_nz))

    i2 = getline('grid_lat_nom')
    if(pe_isroot) then
       print *, 'reading ',trim(data_entries(i2)%w1),' from "', trim(data_entries(i2)%w2),&
            '" of file "',trim(data_entries(i2)%w3),'"'
    end if
    call check(nf90_open(data_entries(i2)%w3, nf90_nowrite, ncid))
    call check(nf90_inq_varid(ncid, data_entries(i2)%w2, varid))
    call check(nf90_get_var(ncid, varid, nom_lat))
    call check(nf90_close(ncid))

    i2 = getline('grid_lon_nom')
    if(pe_isroot) then
       print *, 'reading ',trim(data_entries(i2)%w1),' from "', trim(data_entries(i2)%w2),&
            '" of file "',trim(data_entries(i2)%w3),'"'
    end if
    call check(nf90_open(data_entries(i2)%w3, nf90_nowrite, ncid))
    call check(nf90_inq_varid(ncid, data_entries(i2)%w2, varid))
    call check(nf90_get_var(ncid, varid, nom_lon))
    call check(nf90_close(ncid))

    i2 = getline('grid_z_1d')
    if(pe_isroot) then
       print *, 'reading ',trim(data_entries(i2)%w1),' from "', trim(data_entries(i2)%w2),&
            '" of file "',trim(data_entries(i2)%w3),'"'
    end if
    call check(nf90_open(data_entries(i2)%w3, nf90_nowrite, ncid))
    call check(nf90_inq_varid(ncid, data_entries(i2)%w2, varid))
    call check(nf90_get_var(ncid, varid, depths))
    call check(nf90_close(ncid))
       
    
    ! read in the mask
    i2 = getline('grid_mask_2d')
    if(pe_isroot) then
       print *, 'reading ',trim(data_entries(i2)%w1),' from "', trim(data_entries(i2)%w2),&
            '" of file "',trim(data_entries(i2)%w3),'"'
    end if
    call check(nf90_open(data_entries(i2)%w3, nf90_nowrite, ncid))
    call check(nf90_inq_varid(ncid, data_entries(i2)%w2, varid))
    call check(nf90_get_var(ncid, varid, tmp))
    call check(nf90_close(ncid))
    where (tmp < tiny(0.0))
       mask = .true.
    elsewhere
       mask = .false.
    end where



    !------------------------------------------------------------   
    ! read in the state variable configurations
    if(pe_isroot) then
       print *, ""
       print *, "state variables defnition..."
    end if
    
    ! get the number of state variables
    i1 = 0
    do i=1,num_vars
       if( data_entries(i)%w1(1:5) == 'GRID_') cycle
       i1 = i1 +1
    end do
    if(pe_isroot) then
       print *, "processing ",i1,"state variables"
    end if
    allocate(state_var(i1))
    i1 = 0
    i2 = 1
    i3 = 1
    do i=1,num_vars       
       if( data_entries(i)%w1(1:5) == 'GRID_') cycle
       i1 = i1 + 1
       i3 = max(data_entries(i)%i,1)-1+i2
       if(pe_isroot) then
          print '(A,I5,A,I5,A,A10,5A)', "  slab",i2,' to ',i3,' is ',trim(data_entries(i)%w1), ' from "',&
               trim(data_entries(i)%w2),'" of "', trim(data_entries(i)%w3),'"'
       end if

       state_var(i1) = trim(data_entries(i)%w1)
       slab_var(i2:i3) = i1       

       i2 = i3+1
    end do

    if(i3 /= grid_ns) then
       if(pe_isroot) then
          print *, "ERROR: slabs defined only add up to",i3,". But grid_ns is set to ",grid_ns
       end if
       
       call letkf_mpi_barrier(.true.)
       stop 1

    end if
  end subroutine stateio_init
  !================================================================================



  subroutine stateio_read(self, ens, state)
    class(stateio_nc) :: self
    integer, intent(in) :: ens
    real, intent(out) :: state(:,:,:)

    character(len=1024)  :: filename
    integer :: ncid, varid, stat
    integer :: i, j, i1,i2, i3, i4, k

    do i=1,size(state_var)
        j =  getline(state_var(i))
!       print *, state_var(i), trim(data_entries(j)%w2), "  ",trim(data_entries(j)%w3)


        ! TODO, add logic to use ENS0X, and/or ENSX
        i1=index(toupper(data_entries(j)%w3), '{ENS04}')
        i2 = i1+7
        write (filename, '(A,I0.4,A)') data_entries(j)%w3(1:i1-1), ens, trim(data_entries(j)%w3(i2:))

        ! determine start/stop slabs for this var
        ! TODO, handle differently if slabs are arranged with entire LEVELS grouped together
        do i3=1,grid_ns
           if(slab_var(i3) == i) exit
        end do
        do k=1,grid_ns
           if(slab_var(k) == i) i4 = k
        end do

        call check(nf90_open(trim(filename), nf90_nowrite, ncid), &
             "File: "//trim(filename))
        call check(nf90_inq_varid(ncid, trim(data_entries(j)%w2), varid),&
             "File: "//trim(filename)//" Var:"//trim(data_entries(j)%w2))
        call check(nf90_get_var(ncid, varid, state(:,:,i3:i4)),&
             "File: "//trim(filename)//" Var:"//trim(data_entries(j)%w2))
        call check(nf90_close(ncid), &
             "File: "//trim(filename))
    end do

  end subroutine stateio_read



  subroutine stateio_write(self, filename, state)
    class(stateio_nc) :: self
    character(len=*), intent(in)  :: filename
    real, intent(in) :: state(:,:,:)

    integer :: ncid
    integer :: d_x, d_y, d_z, d_t
    integer :: v_x, v_y, v_z, v_t
    integer :: v_v(size(state_var))
    integer :: i3, i4, i ,j, k

    call check(nf90_create(trim(filename)//'.nc', nf90_write, ncid))

    call check(nf90_def_dim(ncid, "time", nf90_unlimited,  d_t))
    call check(nf90_def_var(ncid, "time", nf90_real, (/d_t/), v_t))
    call check(nf90_put_att(ncid, v_t, "axis", "T"))

    call check(nf90_def_dim(ncid, "grid_x", grid_nx, d_x))
    call check(nf90_def_var(ncid, "grid_x", nf90_real, (/d_x/), v_x))
    call check(nf90_put_att(ncid, v_x, "units", "degrees_east"))
    call check(nf90_put_att(ncid, v_x, "axis", "X"))

    call check(nf90_def_dim(ncid, "grid_y", grid_ny, d_y))
    call check(nf90_def_var(ncid, "grid_y", nf90_real, (/d_y/), v_y))
    call check(nf90_put_att(ncid, v_y, "units", "degrees_north"))
    call check(nf90_put_att(ncid, v_y, "axis", "Y"))

    call check(nf90_def_dim(ncid, "grid_z", grid_nz, d_z))
    call check(nf90_def_var(ncid, "grid_z", nf90_real, (/d_z/), v_z))
    call check(nf90_put_att(ncid, v_z, "units", "meters"))
    call check(nf90_put_att(ncid, v_z, "axis", "Z"))

    
    do i=1,size(state_var)
        j =  getline(state_var(i))
        do i3=1,grid_ns
           if(slab_var(i3) == i) exit
        end do
        do k=1,grid_ns
           if(slab_var(k) == i) i4 = k
        end do
        
        !TODO, allow for 2d variables as well
        call check(nf90_def_var(ncid, trim(data_entries(j)%w2),  nf90_real, (/d_x,d_y,d_z,d_t/), v_v(i)))
!        call check(nf90_def_var_deflate(ncid, v_v(i), 1, 1, 1))
    end do
    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, v_x, nom_lon))
    call check(nf90_put_var(ncid, v_y, nom_lat))
    call check(nf90_put_var(ncid, v_z, depths))

    do i=1,size(state_var)
        j =  getline(state_var(i))
        do i3=1,grid_ns
           if(slab_var(i3) == i) exit
        end do
        do k=1,grid_ns
           if(slab_var(k) == i) i4 = k
        end do
        
        !TODO, allow for 2d variables as well
        call check(nf90_put_var(ncid, v_v(i), state(:,:,i3:i4)))              
    end do


    call check(nf90_close(ncid))

  end subroutine stateio_write



  function toupper(in_str) result(out_str)
    character(*), intent(in) :: in_str
    character(len(in_str))   :: out_str
    integer :: i
    integer, parameter :: offset = 32

    out_str = in_str
    do i =1, len(out_str)
       if (out_str(i:i) >= "a" .and. out_str(i:i) <= "z") then
          out_str(i:i) = achar(iachar(out_str(i:i)) - offset)
       end if
    end do
  end function toupper


  subroutine check(status, message)
    integer, intent(in) :: status
    character(*), optional, intent(in) :: message
    if(status /= nf90_noerr) then
       if(present(message)) then
          write (*,*) trim(nf90_strerror(status)),":  ", trim(message)
       else
          write (*,*) trim(nf90_strerror(status))
       end if
       stop 1
    end if
  end subroutine check
end module letkf_state_nc
