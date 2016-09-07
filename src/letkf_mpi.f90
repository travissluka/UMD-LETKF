module letkf_mpi
  use mpi
  use timing

  integer :: mpi_rank, mpi_size
  
  
contains
  
  subroutine letkf_mpi_init()
    integer :: ierr

    call timing_start('mpi_init')

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, mpi_size, ierr)
    call mpi_comm_rank(mpi_comm_world, mpi_rank, ierr)

    if (mpi_rank == 0) then       
       print *, "Using MPI nproc =",mpi_size
       print *,""
    end if
    
    call timing_end('mpi_init')
  end subroutine letkf_mpi_init

  subroutine letkf_mpi_end
    integer :: ierr
    call mpi_finalize(ierr)
  end subroutine letkf_mpi_end
end module letkf_mpi
