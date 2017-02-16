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
     procedure :: latlon   => stateio_latlon
     procedure :: mask     => stateio_mask
  end type stateio_nc

  integer, parameter :: grid_nz = 40

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
  subroutine stateio_init(self)
    class(stateio_nc) :: self
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

    ! read in the lat/lon
    i1 = getline('grid_lon_1d')
    i2 = getline('grid_lon_2d')
    if( .not. xor(i1 > 0, i2 >0) ) then
       print *, "ERROR: one, and only one of the following must be specified: 'grid_lon_2d' 'grid_lon_1d'"
       stop 1
    end if
    i1 = getline('grid_lat_1d')
    i2 = getline('grid_lat_2d')
    if( .not. xor(i1 > 0, i2 >0) ) then
       print *, "ERROR: only one of the following can be specified: 'grid_lat_2d' 'grid_lat_1d'"
       stop 1
    end if
    
    ! read the z
    i1 = getline('grid_z_1d')
    i2 = getline('grid_z_2d')
    if( .not. xor(i1 > 0, i2 >0) ) then
       print *, "ERROR: only one of the following can be specified: 'grid_z_2d' 'grid_z_1d'"
       stop 1
    end if

   
    !------------------------------------------------------------
    
    ! every other variable remaining is a state variable
    
    
    ! set the variable name properties
    !TODO read this from a configuration file
    allocate(state_var(2))
    state_var(1) = "OCN_T"
    state_var(2) = "OCN_S"
    slab_var(1:40) = 1
    slab_var(41:80) = 2
    do i=1,40
       slab_lvl(i)=i
       slab_lvl(i+40) = i
    end do
    

    ! read in some things from the grid spec file
    ! TODO, move this to stateio_latlon ?
    call check(nf90_open('INPUT/grid_spec.nc', nf90_write, ncid))
    allocate(nom_lon(grid_nx))
    allocate(nom_lat(grid_ny))
    allocate(depths(grid_nz))
    call check(nf90_inq_varid(ncid, "grid_x_T", varid))
    call check(nf90_get_var(ncid, varid, nom_lon))
    call check(nf90_inq_varid(ncid, "grid_y_T", varid))
    call check(nf90_get_var(ncid, varid, nom_lat))
    call check(nf90_inq_varid(ncid, "zt", varid))
    call check(nf90_get_var(ncid, varid, depths))
    call check(nf90_close(ncid))
  end subroutine stateio_init
  !================================================================================


  subroutine stateio_latlon(self, lat, lon)
    class(stateio_nc) :: self
    real, intent(inout) :: lat(:,:)
    real, intent(inout) :: lon(:,:)
    integer :: ncid, varid

    ! pointless statement to prevent "self" not used warnings
    self%i = self%i

    call check(nf90_open('INPUT/grid_spec.nc', nf90_nowrite, ncid))
    call check(nf90_inq_varid(ncid, 'x_T', varid))
    call check(nf90_get_var(ncid, varid, lon))
    call check(nf90_inq_varid(ncid, 'y_T', varid))
    call check(nf90_get_var(ncid, varid, lat))
    call check(nf90_close(ncid))
  end subroutine stateio_latlon



  subroutine stateio_mask(self, mask)
    class(stateio_nc) :: self
    logical, intent(inout) :: mask(:,:)
    integer :: ncid, varid

    real :: wrk(grid_nx,grid_ny)

    ! pointless statement to prevent "self" not used warnings
    self%i = self%i

    call check(nf90_open('INPUT/grid_spec.nc', nf90_nowrite, ncid))
    call check(nf90_inq_varid(ncid, 'wet', varid))
    call check(nf90_get_var(ncid, varid, wrk))
    where(wrk <= tiny(0.0))
       mask = .true.
    elsewhere
       mask = .false.
    end where
    call check(nf90_close(ncid))
  end subroutine stateio_mask






  subroutine stateio_read(self, filename, state)
    class(stateio_nc) :: self
    character(len=*), intent(in)  :: filename
    real, intent(out) :: state(:,:,:)

    integer :: ncid, varid

    call check(nf90_open(trim(filename)//'.nc', nf90_nowrite, ncid))
    call check(nf90_inq_varid(ncid, "temp", varid))
    call check(nf90_get_var(ncid, varid, state(:,:,1:40)))
    call check(nf90_inq_varid(ncid, "salt", varid))
    call check(nf90_get_var(ncid, varid, state(:,:,41:80)))
    call check(nf90_close(ncid))

  end subroutine stateio_read



  subroutine stateio_write(self, filename, state)
    class(stateio_nc) :: self
    character(len=*), intent(in)  :: filename
    real, intent(in) :: state(:,:,:)

    integer :: ncid
    integer :: d_x, d_y, d_z, v_t, v_s, v_x, v_y, v_z

    call check(nf90_create(trim(filename)//'.nc', nf90_write, ncid))
    call check(nf90_def_dim(ncid, "grid_x", grid_nx, d_x))
    call check(nf90_def_var(ncid, "grid_x", nf90_real, (/d_x/), v_x))
    call check(nf90_put_att(ncid, v_x, "units", "degrees_east"))

    call check(nf90_def_dim(ncid, "grid_y", grid_ny, d_y))
    call check(nf90_def_var(ncid, "grid_y", nf90_real, (/d_y/), v_y))
    call check(nf90_put_att(ncid, v_y, "units", "degrees_north"))

    call check(nf90_def_dim(ncid, "grid_z", grid_nz, d_z))
    call check(nf90_def_var(ncid, "grid_z", nf90_real, (/d_z/), v_z))
    call check(nf90_put_att(ncid, v_z, "units", "meters"))

    call check(nf90_def_var(ncid, "temp",  nf90_real, (/d_x,d_y,d_z/), v_t))
    call check(nf90_def_var(ncid, "salt",  nf90_real, (/d_x,d_y,d_z/), v_s))
    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, v_t, state(:,:,1:40)))
    call check(nf90_put_var(ncid, v_s, state(:,:,41:80)))
    call check(nf90_put_var(ncid, v_x, nom_lon))
    call check(nf90_put_var(ncid, v_y, nom_lat))
    call check(nf90_put_var(ncid, v_z, depths))
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


  subroutine check(status)
    integer, intent(in) :: status
    if(status /= nf90_noerr) then
       write (*,*) trim(nf90_strerror(status))
       stop 1
    end if
  end subroutine check
end module letkf_state_nc
