module letkf_mpi
  use mpi
  use global
  
  implicit none
  
  integer, save :: pe_root, pe_rank, pe_size
  integer, save :: mpi_comm_letkf

  logical, save :: pe_isroot

  integer, save, allocatable :: ens_list(:)
  !! A list of ensemble numbers (from 0 to MEM) that this
  !! process is responsible for the I/O
  integer, save :: ens_count

  integer, save, allocatable :: ens_map(:)
  !! for each index MEM, the number of the PE
  !! respsonsible for its I/O
  
  
  integer, save, allocatable :: ij_list(:)
  !! A list of grid point locations that this process
  !! is responsible for doing core LETKF algorithm on
  integer, save :: ij_count


  real, allocatable, save :: load_weights(:)
  
  integer, allocatable, save :: scatterv_count(:)
  integer, allocatable, save :: scatterv_displ(:)
  
contains


  function ens_idx(m)
    integer, intent(in) :: m
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
  subroutine letkf_mpi_init2(mem, grid_ij)
    integer, intent(in) :: mem
    integer, intent(in) :: grid_ij
    integer :: ierr, i, j
    integer :: count, prev
    !TODO: allow for more complicated layouts
        
    call mpi_barrier(mpi_comm_letkf, ierr)

    !TODO load weights
    allocate(load_weights(pe_size))
    load_weights = 1
    load_weights(1) = 2.325
    load_weights(2) = 0.7
    load_weights(3) = 0.375
    load_weights(4) = 0.75
    load_weights = load_weights / sum(load_weights)
    
    ! print header info
    if (pe_isroot) then
       print '(A)', ""
       print '(A)', "MPI configuration"
       print '(A)', "------------------------------------------------------------"
       print '(A,I4)', "Using MPI nproc =", pe_size
       print '(A,I4)', "Using MPI root  =", pe_root
    end if

    allocate (ens_map(mem))
    
    ! determine which ensemble members number this PE should deal with
    ! ----------------------------------------    
    prev = 0   
    if (pe_isroot) then
       print *, ""
       print '(A)', "ensemble member I/O list:"
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
    allocate(scatterv_count(pe_size))
    allocate(scatterv_displ(pe_size))    
    if (pe_isroot) then
       print *, ""
       print "(A)", "gridpoint assignment list:"
    end if    
    prev = 0
    do i=0, pe_size-1      
       count = nint(grid_ij * load_weights(i+1))
       if (i == pe_size-1) count = grid_ij - prev
       
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

  end subroutine letkf_mpi_init2
  

  
  !############################################################
  subroutine letkf_mpi_end
    integer :: ierr
    call mpi_finalize(ierr)
  end subroutine letkf_mpi_end


  
end module letkf_mpi
