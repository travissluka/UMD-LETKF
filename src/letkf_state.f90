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

  subroutine letkf_state_write(filename, state_3d, state_2d)
    character(len=*), intent(in) :: filename
    real, allocatable, intent(in) :: state_3d(:,:,:,:)
    real, allocatable, intent(in) :: state_2d(:,:,:)

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
