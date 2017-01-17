module letkf_mpi
  !! MPI handler, handles scattering and gathering of state grids and ensemble members
  use mpi

  implicit none
  private

  public :: letkf_mpi_preinit, letkf_mpi_init, letkf_mpi_final
  public :: letkf_mpi_obs, letkf_mpi_grd2ij, letkf_mpi_ij2grd
  public :: letkf_mpi_ens2ij, letkf_mpi_ij2ens
  public :: ens_list
  public :: letkf_mpi_barrier

  integer, parameter :: dp = kind(0.0)


  integer, public :: mem
    !! Number of ensemble members

  ! public variables
  ! ------------------------------
  integer,public, allocatable :: ens_map(:)
    !! for each index MEM, the number of the PE
    !! respsonsible for its I/O

  integer, public :: ij_count
    !! Number of grid points the current process is responsible for handling

  integer, public, allocatable :: scatterv_count(:)
  integer, public, allocatable :: scatterv_displ(:)


  ! private variables
  ! ------------------------------
  integer, allocatable :: ens_list(:)
    !! lists which process is responsible for the I/O for
    !! each ensemble member number.
    !! ***Size is ([[letkf_mpi:mem]])***

  integer :: ens_count
    !! number of ensemble members for which is process is responsible for the I/O

  integer, allocatable :: ij_list(:)
    !! A list of grid point locations that this process
    !! is responsible for doing core LETKF algorithm on

  real, allocatable :: load_weights(:)


  integer, public :: pe_rank, pe_size, pe_root
  logical, public :: pe_isroot
  integer, public :: mpi_comm_letkf


  integer :: grid_nx, grid_ny, grid_ns

contains


  subroutine  letkf_mpi_barrier()
    integer :: ierr
    call mpi_barrier(mpi_comm_letkf, ierr)
    if(ierr /= 0) then
       print *, "ERROR: with letkf_mpi_barrier()"
       stop 1
    end if
  end subroutine letkf_mpi_barrier


  subroutine letkf_mpi_ens2ij(ens, ij)
    !! takes ensemble members read in by 1 or more process and distributes
    !! to all processes. Afterward each process will have all ensemble members for a
    !! subset of grid points.

    real, intent(in)    :: ens(grid_nx, grid_ny, grid_ns, mem)
      !! The ensemble members this given process is responsible for loading in.
      !! ***Shape is ([[letkf_state:grid_nx]], [[letkf_state:grid_ny]],
      !!    [[letkf_state:grid_ny]], [[letkf_mpi:mem]])***

    real, intent(inout) :: ij(mem, grid_ns, ij_count)
      !! The gridpoints this given process is responsible for computing.
      !! ***Shape is ([[letkf_mpi:ij_count]], [[letkf_state:grid_ns]], [[letkf_mpi:mem]])***

    integer :: ierr, i ,m
    real :: wrk(ij_count,grid_ns,mem)

    do m=1,mem
       do i=1,size(ens,3)
!          call mpi_scatterv(ens(:,:,i,ens_idx(m)), scatterv_count, scatterv_displ, mpi_real, &
!               ij(:,i,m), ij_count, mpi_real, ens_map(m), mpi_comm_letkf, ierr)
          call mpi_scatterv(ens(:,:,i,ens_idx(m)), scatterv_count, scatterv_displ, mpi_real, &
               wrk(:,i,m), ij_count, mpi_real, ens_map(m), mpi_comm_letkf, ierr)
          if(ierr /= 0) then
             print *, "ERROR: with letkf_mpi_ens2ij",ierr
             stop 1
          end if
       end do
    end do
!    ij = reshape(wrk, (/ij_count, grid_ns, mem/),order=(/1,2,3/))
    ij = reshape(wrk, (/mem, grid_ns, ij_count/), order=(/3,2,1/))

  end subroutine letkf_mpi_ens2ij



  subroutine letkf_mpi_ij2ens(ij, ens)
    !! Takes gridpoints that are scattered across the processors and
    !! combines them and sends individual whole ensemble members to the processors
    !! responsible for saving them

    real, intent(in)    :: ij(:,:,:)
      !! The gridpoints this given process is responsible for computing.
      !! *** Shape is ([[letkf_mpi:ij_count]], [[letkf_state:grid_ns]], [[letkf_mpi:mem]])***

    real, intent(inout) :: ens(:,:,:,:)
      !! The ensemble members this given process is responsible for loading in.
      !! *** Shape is ([[letkf_state:grid_nx]], [[letkf_state:grid_ny]],
      !!  [[letkf_state:grid_ns]], [[letkf_mpi:mem]])***

    integer :: ierr, i, m

    do m=1,mem
       do i=1,size(ens,3)
          call mpi_gatherv(ij(:,i,m), ij_count, mpi_real,&
               ens(:,:,i,ens_idx(m)), scatterv_count, scatterv_displ, mpi_real,&
               ens_map(m), mpi_comm_letkf, ierr)
          if(ierr /= 0) then
             print *, "ERROR: with letkf_mpi_ij2ens",ierr
             stop 1
          end if
       end do
    end do
  end subroutine letkf_mpi_ij2ens



  subroutine letkf_mpi_grd2ij(grd, ij)
    !! Takes a single grid on the root process and distributes
    !! portions of it to the worker processes

    real, intent(in) :: grd(grid_nx, grid_ny)
      !! The 2D grid to send.
      !! ***Shape is ([[letkf_state:grid_nx]], [[letkf_state:grid_ny]])***

    real, intent(inout) :: ij(ij_count)
      !! The gridpoints this process is responsible for using.
      !! ***Shape is ([[letkf_mpi:ij_count]])***

    integer :: ierr

    call mpi_scatterv(grd, scatterv_count, scatterv_displ, mpi_real, &
         ij, ij_count, mpi_real, pe_root, mpi_comm_letkf, ierr)
    if(ierr /= 0) then
       print *, "ERROR: with letkf_mpi_grd2ij",ierr
       stop 1
    end if
  end subroutine letkf_mpi_grd2ij



  subroutine letkf_mpi_ij2grd(ij, grd)
    real, intent(in)    :: ij(ij_count)
    real, intent(inout) :: grd(grid_nx, grid_ny)
    integer :: ierr
    call mpi_gatherv(ij, ij_count, mpi_real, &
         grd, scatterv_count, scatterv_displ, mpi_real, pe_root, mpi_comm_letkf, ierr)
    if(ierr /= 0) then
       print *, "ERROR: with letkf_mpi_ij2grd", ierr
       stop 1
    end if
  end subroutine letkf_mpi_ij2grd



  subroutine letkf_mpi_obs(ohx, qc)
    real(dp), intent(inout) :: ohx(:,:)
    integer,  intent(inout) :: qc(:,:)
    integer :: ierr, i

    !TODO, this is inefficient
    call mpi_allreduce(mpi_in_place, ohx, mem*size(ohx,2), mpi_real, mpi_sum, mpi_comm_letkf, ierr)
    call mpi_allreduce(mpi_in_place, qc, mem*size(ohx,2), mpi_integer, mpi_sum, mpi_comm_letkf, ierr)
  end subroutine letkf_mpi_obs


  function ens_idx(m)
    !! Get the global ensemble number of the m-th member this process is responsible for.
    !! Returns -1 if m is out of bounds, or this process does not have any members to load
    integer, intent(in) :: m
      !! the \(m^{th}\) ensemble member this process is responsible for handling for I/O.
      !! ***Range is [1, [[letkf_mpi:ens_count]] ]***
    integer ::  ens_idx
    integer :: i
    do i=1,ens_count
       if (ens_list(i) == m) then
          ens_idx = i
          return
       end if
    end do
    ens_idx = -1
  end function ens_idx



  !############################################################
  subroutine letkf_mpi_preinit()
    integer :: ierr

    call mpi_init(ierr)

    ! rightn now all processes are used for the LETKF,
    ! but this leave the option open for someone needing to define
    ! a comm subgroup
    mpi_comm_letkf = mpi_comm_world
    pe_root = 0

    call mpi_comm_size(mpi_comm_letkf, pe_size, ierr)
    call mpi_comm_rank(mpi_comm_letkf, pe_rank, ierr)

    pe_isroot = pe_root == pe_rank


  end subroutine letkf_mpi_preinit



  !############################################################
  subroutine letkf_mpi_init(mem, nx, ny, ns)
    integer, intent(in) :: mem
    integer, intent(in) :: nx, ny, ns
    integer :: i, j
    integer :: count, prev

    grid_nx=nx
    grid_ny=ny
    grid_ns=ns

    if(pe_isroot) then
       print *, new_line('a'), &
            new_line('a'), '============================================================', &
            new_line('a'), ' letkf_mpi_init() : ', &
            new_line('a'), '============================================================'
    end if
    allocate (ens_map(mem))

    if(pe_isroot) then
       print *, ""
       print '(A,I5)', " MPI processors: ", pe_size
       if (mpi_comm_letkf == mpi_comm_world) then
          print *,"MPI communicator: MPI_COMM_WORLD"
       else
          print *, "MPI communicator: ",mpi_comm_letkf
       end if

    end if

    ! determine which ensemble members number this PE should deal with
    ! ----------------------------------------
    prev = 0
    if (pe_isroot) then
       print *, ""
       print *, "ensemble member I/O list:"
    end if
    do i=0, pe_size-1
       ! for each proc, determine the number of members required
       count = mem / pe_size
       if (i < mod(mem, pe_size)) count = count + 1
       if (i == pe_rank) then
          ! if we are the mpi rank we are currently calculating for... save this
          ens_count = count
          allocate(ens_list(count))
          do j=1,count
             ens_list(j) = prev+j
          end do
       end if
       do j=1,count
          ens_map(prev+j) = i
       end do

       ! if root proc, print out info about each PE
       if (pe_isroot) then
          if (count > 1) then
             print '(A5,I4,A,I4,A,I4,A,I4)', "proc", i,' is I/O for ',count,' ensemble member(s): ', prev+1, ' to ',count+prev
          else if (count == 1) then
             print '(A5,I4,A,I4)', " proc ", i,' is I/O for    1 ensemble member(s): ', prev+1
          else
             print '(A5,I4,A)', " proc ", i,' is I/O for    0 ensemble member(s)'
          end if
       end if

       prev = prev + count
    end do



    ! determine how many gridpoints this PE should deal with
    ! ----------------------------------------
    if (pe_isroot) print *, ""

    ! determine the process load weights
    allocate(load_weights(pe_size))
    load_weights = 1.0
    load_weights = load_weights / sum(load_weights)


    ! calculate actual numer of gridpoints to use
    if (pe_isroot) then
       print *, ""
       print *, "gridpoint assignment list:"
    end if
    allocate(scatterv_count(pe_size))
    allocate(scatterv_displ(pe_size))
    prev = 0
    do i=0, pe_size-1
       count = nint(grid_nx*grid_ny * load_weights(i+1))
       if (i == pe_size-1) count = grid_nx*grid_ny - prev

       if (i == pe_rank) then
          ij_count = count
          allocate(ij_list(count))
          do j = 1,count
             ij_list(j) = prev + j
          end do
       end if
       scatterv_count(i+1) = count
       scatterv_displ(i+1) = prev

       if (pe_isroot) then
          print '(A5,I4,A,I7,A,I7,A,I7)', "proc",i,' calcs for', count,' grid points: ',prev+1,' to ',prev+count
       end if
       prev = prev+count
    end do

  end subroutine letkf_mpi_init



  !############################################################
  subroutine letkf_mpi_final
    integer :: ierr
    real :: mem, mem_avg, mem_min, mem_max, mem_sum
    ! calculate memore hiwater mark
    call getMaxMem(mem)
    call mpi_reduce(mem, mem_min, 1, mpi_real, mpi_min, pe_root, mpi_comm_letkf, ierr)
    call mpi_reduce(mem, mem_max, 1, mpi_real, mpi_max, pe_root, mpi_comm_letkf, ierr)
    call mpi_reduce(mem, mem_sum, 1, mpi_real, mpi_sum, pe_root, mpi_comm_letkf, ierr)
    mem_avg = mem_sum / pe_size

    if(pe_isroot) then
       print *, new_line('a'), &
            new_line('a'), "============================================================",&
            new_line('a'), " Max memory usage",&
            new_line('a'), "============================================================"
       print *, "per core:"
       print '(A,F5.2,A)', "   avg: ",mem_avg," GB"
       print '(A,F5.2,A)', "   min: ",mem_min," GB"
       print '(A,F5.2,A)', "   max: ",mem_max," GB"
       print *,"total:"
       print '(A,F5.2,A)', "    ",mem_sum," GB"
    end if
    call mpi_finalize(ierr)
  end subroutine letkf_mpi_final



  subroutine getMaxMem(mem)
    !! return memory hiwater mark, in gigabytes.
    !! this relies on the VmHWM field set in the /proc/(pid)/status file
    !! if this file or field can't be found, -1 is returned
    real, intent(out)   :: mem
    character(len=30)   :: pid_char
    character(len=1024) :: proc_file, line
    integer             :: unit, iostat, i
    logical             :: ex

    ! get process id
    write(pid_char,'(I10)') getpid()
    pid_char = adjustl(pid_char)

    ! make sure the file exists that contains memory hi usage info in linux
    proc_file="/proc/"//trim(pid_char)//"/status"
    inquire(file=proc_file, exist=ex)
    if (.not. ex) then
       print *, "can't open ",trim(proc_file)
       mem = -1
       return
    end if

    ! read in the VmHWM line of the file
    open(newunit=unit, file=proc_file, action='read')
    iostat = 0
    do while (iostat==0)
       read(unit,'(A)',iostat=iostat) line
       i = scan(line, ':')
       if(i<=0) cycle
       if(line(:i) == "VmHWM:") then
          ! the VmHWM line has been found, parse out the number of kilobytes
          line = line(i+1:)
          do i=1,len(line)
             if(line(i:i) == char(9)) line(i:i) = ' '
          end do
          line = adjustl(line)
          read(line,*) mem
          ! return the number of gigabytes
          mem = mem/1024/1024
          close(unit)
          return
       end if
    end do
    mem = -1
    close(unit)
  end subroutine getMaxMem
end module letkf_mpi
