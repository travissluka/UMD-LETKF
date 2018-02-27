module getmem

  use mpi
  use iso_fortran_env
#ifdef __INTEL_COMPILER
  use ifport, only : getpid
#endif

  implicit none
  private

  public :: getmem_init
  public :: getmem_print


  integer :: pe_rank
  integer :: pe_root
  integer :: pe_comm
  integer :: pe_size
  logical :: pe_isroot

  real    :: all_min, all_max, all_avg

contains


  !================================================================================
  !================================================================================
  subroutine getmem_init(p_pe_root, p_pe_comm)
    integer ::  p_pe_root, p_pe_comm, ierr

    pe_root = p_pe_root
    pe_comm = p_pe_comm

    call mpi_comm_rank(pe_comm, pe_rank, ierr)
    call mpi_comm_size(pe_comm, pe_size, ierr)
    pe_isroot = pe_rank == pe_root

    all_min = 0
    all_max = 0
    all_avg = 0
  end subroutine getmem_init
  !================================================================================



  !================================================================================
  !================================================================================
  subroutine getmem_print
    real :: cur, max
    real :: mem_min, mem_max, mem_avg
    integer :: ierr

    call readmem(cur, max)
    call mpi_reduce(cur, mem_min, 1, mpi_real, mpi_min, pe_root, pe_comm, ierr)
    call mpi_reduce(cur, mem_max, 1, mpi_real, mpi_max, pe_root, pe_comm, ierr)
    call mpi_reduce(cur, mem_avg, 1, mpi_real, mpi_sum, pe_root, pe_comm, ierr)
    mem_avg = mem_avg / pe_size
    if(pe_isroot) then
       print *, ""
       print '(A,4F8.2)', "Current Mem in GB (min/max/avg/tot): ",&
            mem_min, mem_max, mem_avg, mem_avg*pe_size
    end if

    !TODO: there's a way with some linux kernels to clear the previous max memory
    ! try this out? otherwise the process total max is going to be wrong near
    ! the end of the program
    call mpi_reduce(max, mem_min, 1, mpi_real, mpi_min, pe_root, pe_comm, ierr)
    call mpi_reduce(max, mem_max, 1, mpi_real, mpi_max, pe_root, pe_comm, ierr)
    call mpi_reduce(max, mem_avg, 1, mpi_real, mpi_sum, pe_root, pe_comm, ierr)
    mem_avg = mem_avg / pe_size
    if(pe_isroot) then
       print '(A,4F8.2)', "Maximum Mem in GB (min/max/avg/tot): ",&
            mem_min, mem_max, mem_avg, mem_avg * pe_size
    end if

    if(mem_min > all_min) all_min = mem_min
    if(mem_max > all_max) all_max = mem_max
    if(mem_avg > all_avg) all_avg = mem_avg

!    call memReset()
  end subroutine getmem_print
  !================================================================================


  subroutine getmem_max
  end subroutine getmem_max

  ! !================================================================================
  ! !================================================================================
  ! subroutine letkf_mpi_final
  !   integer :: ierr
  !   real :: cur, mem, mem_avg, mem_min, mem_max, mem_sum
  !   ! calculate memore hiwater mark
  !   call getMem(cur, mem)
  !   call mpi_reduce(mem, mem_min, 1, mpi_real, mpi_min, pe_root, pe_comm, ierr)
  !   call mpi_reduce(mem, mem_max, 1, mpi_real, mpi_max, pe_root, pe_comm, ierr)
  !   call mpi_reduce(mem, mem_sum, 1, mpi_real, mpi_sum, pe_root, pe_comm, ierr)
  !   mem_avg = mem_sum / pe_size

  !   if(pe_isroot) then
  !      print *, new_line('a'), &
  !           new_line('a'), "============================================================",&
  !           new_line('a'), " Max memory usage",&
  !           new_line('a'), "============================================================"
  !      print *, "per core:"
  !      print '(A,F5.2,A)', "   avg: ",mem_avg," GB"
  !      print '(A,F5.2,A)', "   min: ",mem_min," GB"
  !      print '(A,F5.2,A)', "   max: ",mem_max," GB"
  !      print *,"total:"
  !      print '(A,F5.2,A)', "    ",mem_sum," GB"
  !   end if
  !   call mpi_finalize(ierr)
  ! end subroutine letkf_mpi_final
  ! !================================================================================



  ! subroutine memReset()
  !   character(len=1024) :: proc_file
  !   character(len=30)   :: pid_char

  !   write(pid_char, '(I10)') getpid()
  !   pid_char = adjustl(pid_char)

  !   proc_file="/proc/"//trim(pid_char)//"/clear_refs"
  !   call system("echo 1 > "//trim(proc_file))
  ! end subroutine memReset



  !================================================================================
  !================================================================================
  subroutine readmem(cur, max)
    !! return memory hiwater mark, in gigabytes.
    !! this relies on the VmHWM field set in the /proc/(pid)/status file
    !! if this file or field can't be found, -1 is returned
    real, intent(out)   :: cur, max
    character(len=30)   :: pid_char
    character(len=1024) :: proc_file, line
    integer             :: unit, iostat, i
    logical             :: ex

    ! get process id
    write(pid_char,'(I10)') getpid()
    pid_char = adjustl(pid_char)

    cur = -1
    max = -1

    ! make sure the file exists that contains memory hi usage info in linux
    proc_file="/proc/"//trim(pid_char)//"/status"
    inquire(file=proc_file, exist=ex)
    if (.not. ex) then
       print *, "can't open ", trim(proc_file)
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
          read(line,*) max
          ! return the number of gigabytes
          max = max/1024/1024
       else  if(line(:i) == "VmRSS:") then
          ! the VmHWM line has been found, parse out the number of kilobytes
          line = line(i+1:)
          do i=1,len(line)
             if(line(i:i) == char(9)) line(i:i) = ' '
          end do
          line = adjustl(line)
          read(line,*) cur
          ! return the number of gigabytes
          cur = cur/1024/1024
       end if
    end do
    close(unit)
  end subroutine readmem
  !================================================================================


end module getmem
