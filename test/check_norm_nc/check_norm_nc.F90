program check_norm_nc
  use netcdf
  implicit none
  
  character(len=1024) :: file1, file2, str
  real :: max_norm

  integer :: i, j, s(4)
  integer :: nc1, nc2
  integer :: ndim1, ndim2, nvar1, nvar2
  integer, allocatable :: dimlen1(:), dimlen2(:)
  integer :: dimids1(NF90_MAX_VAR_DIMS), dimids2(NF90_MAX_VAR_DIMS)
  real, allocatable :: tmpr1(:,:,:,:), tmpr2(:,:,:,:)
  real :: n

  
  call get_command_argument(1, file1)
  call get_command_argument(2, file2)
  call get_command_argument(3, str)
  read (str, *) max_norm

  
  call check(nf90_open(file1, nf90_nowrite, nc1))
  call check(nf90_open(file2, nf90_nowrite, nc2))

  ! check the dimensions
  call check(nf90_inquire(nc1, ndim1, nvar1))
  call check(nf90_inquire(nc2, ndim2, nvar2))
  if(ndim1 /= ndim2) then
     print *, "ndim mismatch"
     stop 1
  end if
  allocate( dimlen1(ndim1), dimlen2(ndim2) )
  do i =1, ndim1
     ! check lengths
     call check(nf90_inquire_dimension(nc1, i, len=dimlen1(i)))
     call check(nf90_inquire_dimension(nc2, i, len=dimlen2(i)))
     if(dimlen1(i) /= dimlen2(i)) then
        print *, "dim length mismatch"
        stop 1
     end if
     ! TODO: check name
     ! TODO: check attributes
  end do

  
  ! check the variables
  if(nvar1 /= nvar2) then
     print *, "nvar mismatch"
     stop 1
  end if
  do i = 1, nvar1
     ! dimension lengths
     call check(nf90_inquire_variable(nc1, i, ndims=ndim1))
     call check(nf90_inquire_variable(nc2, i, ndims=ndim2))
     if(ndim1 /= ndim2) then
        print *, "dimension length of variable do not match"
        stop 1
     end if

     ! dimension ids
     call check(nf90_inquire_variable(nc1, i, dimids=dimids1))
     call check(nf90_inquire_variable(nc2, i, dimids=dimids2))

     ! variable size
     s=1     
     do j=1,ndim1
        s(j) = dimlen1(dimids1(j))
     end do

     ! variable value
     allocate(tmpr1(s(1),s(2),s(3),s(4)))
     allocate(tmpr2(s(1),s(2),s(3),s(4)))     
     call check(nf90_get_var(nc1,i, tmpr1))
     call check(nf90_get_var(nc2,i, tmpr2))
     n = maxval((tmpr1-tmpr2)*(tmpr1-tmpr2))
     if(n>max_norm) then
        print *, "max difference of ",n
        stop 1
     end if
     deallocate(tmpr1, tmpr2)
  end do

  
contains

  
  subroutine check(status)
    integer, intent(in) :: status
    
    if(status /= nf90_noerr) then
       write (*,*) trim(nf90_strerror(status))
       stop 100
    end if
  end subroutine check
  
end program check_norm_nc
