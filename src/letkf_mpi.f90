module letkf_mpi
  use mpi
  implicit none
  
  integer, save :: pe_root, pe_rank, pe_size
  integer, save :: mpi_comm_letkf

  logical, save :: pe_isroot

  integer, save, allocatable :: ens_list(:)
  integer, save :: ens_count


  
contains


    
  !############################################################  
  subroutine letkf_mpi_init()
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
    
  end subroutine letkf_mpi_init



  !############################################################
  subroutine letkf_mpi_init2(mem)
    integer, intent(in) :: mem    
    integer :: ierr, i, j
    integer :: count, prev
    !TODO: allow for more complicated layouts
        
    call mpi_barrier(mpi_comm_letkf, ierr)

    ! print header info
    if (pe_isroot) then
       print '(A)', ""
       print '(A)', "MPI configuration"
       print '(A)', "------------------------------------------------------------"
       print '(A,I4)', "Using MPI nproc =", pe_size
       print '(A,I4)', "Using MPI root  =", pe_root
       print '(A)', "ensemble member I/O list:"
    end if


    prev = 0
    do i=0, pe_size-1      
       !determine how many ensemble members number this PE should deal with
       count = mem / pe_size
       if (i < mod(mem, pe_size)) count = count + 1

       !
       if (i == pe_rank) then
          ens_count = count
          allocate(ens_list(count))
          do j=1,count
             ens_list(j) = prev+j
          end do
       end if

       ! if root proc, print out info about each PE
       if (pe_isroot) then
          if (count > 1) then
             print '(A,I4,A,I4,A,I4,A,I4)', " proc ", i,' is I/O for ',count,' ensemble member(s): ', prev+1, ' to ',count+prev
          else if (count == 1) then
             print '(A,I4,A,I4)', " proc ", i,' is I/O for    1 ensemble member(s): ', prev+1
          else
             print '(A,I4,A)', " proc ", i,' is I/O for    0 ensemble member(s)'
          end if
       end if
       
       prev = prev + count
    end do
    
  end subroutine letkf_mpi_init2
  

  !############################################################
  subroutine letkf_mpi_end
    integer :: ierr
    call mpi_finalize(ierr)
  end subroutine letkf_mpi_end


  
end module letkf_mpi
