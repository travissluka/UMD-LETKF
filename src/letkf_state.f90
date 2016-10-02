module letkf_state
  use global

  implicit none


  !! @todo make these configurable
  integer :: grid_x = 192
  integer :: grid_y = 94
  integer :: grid_z = 64
  integer :: grid_3d = 5
  integer :: grid_2d = 1
  
contains


  !============================================================
  subroutine letkf_state_read(filename, state_3d, state_2d)
    character(len=*), intent(in) :: filename
    real, intent(inout) :: state_3d(:,:,:,:)
    real, intent(inout) :: state_2d(:,:,:)
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
    character(len=*), intent(in) :: filename
    real, intent(in) :: state_3d(:,:,:,:)
    real, intent(in) :: state_2d(:,:,:)

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
