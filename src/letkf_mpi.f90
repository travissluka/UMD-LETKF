module letkf_mpi
  use mpi

  integer, save :: pe_root, pe_rank, pe_size
  integer, save :: mpi_comm_letkf

  logical, save :: pe_isroot
  
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
    

    if (pe_isroot) then       
       print *, "Using MPI nproc =", pe_size
       print *, "Using MPI root  =", pe_root
       print *,""
    end if

  end subroutine letkf_mpi_init

  

  !############################################################
  subroutine letkf_mpi_end
    integer :: ierr
    call mpi_finalize(ierr)
  end subroutine letkf_mpi_end


  
end module letkf_mpi
