module letkf_obs_dat
  use letkf_common
  use letkf_obs
  use letkf_mpi

  implicit none
  private

  public :: obsio_dat

  !------------------------------------------------------------

  type, extends(obsio) :: obsio_dat
   contains
     procedure :: init => obsio_dat_init
     procedure :: write => obsio_dat_write
     procedure :: read => obsio_dat_read
  end type obsio_dat


contains


  subroutine obsio_dat_init(self)
    class(obsio_dat) :: self

    self%description = "raw observation I/O"
    self%extension   = "dat"
  end subroutine obsio_dat_init


  subroutine obsio_dat_write(self, file, obs, iostat)
    class(obsio_dat) :: self
    character(len=*), intent(in) :: file
    type(observation), intent(in) :: obs(:)
    integer, optional, intent(out) :: iostat

    if (present(iostat)) iostat = 0
    print *, "ERROR: did not write ",trim(file), &
         " this method not yet implemented"
    stop 1
  end subroutine obsio_dat_write



  subroutine obsio_dat_read(self, file, obs, obs_innov, obs_qc, iostat)
    class(obsio_dat) :: self
    character(len=*), intent(in) :: file
    type(observation), allocatable, intent(out) :: obs(:)
    real(dp), allocatable, intent(out) :: obs_innov(:)
    integer,  allocatable, intent(out) :: obs_qc(:)
    integer, optional, intent(out) :: iostat

    logical :: ex
    integer :: filesize, unit, i
    real(kind=4) :: record(10)

    ! determine the number of observations that will be read in
    inquire(file=file, size=filesize)
    allocate( obs(filesize/4/12) ) !TODO, do this better
    allocate( obs_innov(filesize/4/12) )
    allocate( obs_qc(filesize/4/12) )

    ! make sure th at the desired file exists
    inquire(file=file, exist=ex)
    if(.not. ex) then
       print *, "ERROR: Unable to open observation file ",trim(file)
       stop 1
    end if

    ! open the file and read in observations
    open(newunit=unit, file=file, form='unformatted', access='sequential', action='read')
    do i = 1, size(obs)
       read(unit) record
       obs(i)%id    = int(record(1))
       obs(i)%lon   = record(2)
       obs(i)%lat   = record(3)
       obs(i)%depth = record(4)
       obs(i)%val   = record(5)
       obs(i)%err   = record(6)
       obs(i)%plat  = int(record(7))
       obs(i)%time  = record(8)
       obs_innov(i) = record(9)
       obs_qc(i)    = int(record(10))

       !TODO: temporary, remove this
       if (obs(i)%id == 1100) obs_innov(i) = obs_innov(i) / 100.0
!       if (obs_qc(i) == 0) then
!          obs_qc(i) = 1
!       else
!          obs_qc(i) = 0
!       end if
    end do
    close(unit)

    !!@todo do something actually useful with iostat, or get rid of it
    if (present(iostat)) iostat = 1
  end subroutine obsio_dat_read


end module letkf_obs_dat
