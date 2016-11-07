module gfs_mod
  use letkf
  use netcdf 
  implicit none

  type, extends(stateio) :: gfs_stateio
   contains
     procedure :: read => gfs_obsio_read
     procedure :: write => gfs_obsio_write
  end type gfs_stateio

contains
  subroutine vcoord(i,j,z)
    integer, intent(in) :: i, j
    real, intent(inout) :: z(:)
  end subroutine vcoord

  subroutine gfs_obsio_read(self, file, state)
    class(gfs_stateio) :: self
    character(len=*), intent(in) :: file
    real, allocatable, intent(out) :: state(:,:,:)
  end subroutine gfs_obsio_read

  subroutine gfs_obsio_write(self, file, state)
    class(gfs_stateio) :: self
    character(len=*), intent(in) :: file
    real, allocatable, intent(in) :: state(:,:,:)
  end subroutine gfs_obsio_write
  
  subroutine gfs_read_latlon(lon_grid, lat_grid)
    implicit none
    real, allocatable, intent(out) :: lon_grid(:,:)
    real, allocatable, intent(out) :: lat_grid(:,:)

    integer :: ncid, ierr, varid, n
  
    allocate(lon_grid(grid_x, grid_y))
    allocate(lat_grid(grid_x, grid_y))

    ierr = nf90_open('data/grid_spec.nc', nf90_nowrite, ncid)                                           
    ierr = nf90_inq_varid(ncid, 'xta', varid)
    ierr = nf90_get_var(ncid, varid, lon_grid(:,1))
    ierr = nf90_inq_varid(ncid, 'yta', varid)
    ierr = nf90_get_var(ncid, varid, lat_grid(1,:))
    ierr = nf90_close(ncid)

    do n=2,grid_x
       lat_grid(n,:) = lat_grid(1,:)
    end do
    do n=2, grid_y
       lon_grid(:,n) = lon_grid(:,1)
    end do
  end subroutine gfs_read_latlon
  
end module gfs_mod

! ------------------------------------------------------------
! ------------------------------------------------------------

program test
  use gfs_mod
  use letkf
  use letkf_obs_dat
  implicit none
 
  type(gfs_stateio) ::gfsio_class
  type(obsio_dat) :: obsio_class
  real, allocatable :: lon_grid(:,:)
  real, allocatable :: lat_grid(:,:)  
  call letkf_driver_init()
  
  call gfs_read_latlon(lon_grid, lat_grid)
  call letkf_set_grid(lon_grid, lat_grid)
  call letkf_set_vcoord(vcoord)
  
  call letkf_set_stateio(gfsio_class)
  call letkf_set_obsio(obsio_class)
  
  call letkf_driver_run()

end program test

