module letkf_mpi
  !! MPI handler, handles scattering and gathering of state grids and ensemble members

  use mpi
  use iso_fortran_env

  implicit none
  private

#define SCATTER_CUSTOMPACK


  ! public module methods
  !------------------------------------------------------------
  public :: letkf_mpi_preinit, letkf_mpi_init, letkf_mpi_final
  public :: letkf_mpi_setgrid
  public :: letkf_mpi_obs, letkf_mpi_ij2grd !, letkf_mpi_grd2ij
  public :: letkf_mpi_ens2ij, letkf_mpi_ij2ens
  public :: letkf_mpi_barrier

  public :: letkf_mpi_grd2ij_real
  public :: letkf_mpi_grd2ij_logical
  !TODO, why does this interface not work
!   interface letkf_mpi_grd2ij
!      module procedure letkf_mpi_grd2ij_real
!      module procedure letkf_mpi_grd2ij_logical
!   end interface

  ! public module variables
  ! ------------------------------------------------------------
  integer, public, protected :: mem
  !! Number of ensemble members

  integer, public, protected, allocatable :: ens_list(:)
    !! lists which process is responsible for the I/O for
    !! each ensemble member number.
    !! ***Size is ([[letkf_mpi:mem]])***

  integer,public, protected, allocatable :: ens_map(:)
    !! for each index MEM, the number of the PE
    !! respsonsible for its I/O

  integer, public, protected :: ij_count
    !! Number of grid points the current process is responsible for handling

  integer, public, protected :: pe_rank, pe_size, pe_root
  logical, public, protected :: pe_isroot
  integer, public, protected :: mpi_comm_letkf

  

  ! private variables
  ! ------------------------------
  integer :: ens_count
    !! number of ensemble members for which is process is responsible for the I/O

  integer, allocatable :: ij_list(:)
    !! A list of grid point locations that this process
    !! is responsible for doing core LETKF algorithm on

  integer, allocatable :: scatterv_count(:)
  integer, allocatable :: scatterv_displ(:)

  logical :: interleave = .true.
  integer :: grid_nx
    !! number of grid points in the x direction (longitude usually)  

  integer :: grid_ny
    !! number of grid points in the y direction (latitude usually)  

  integer :: grid_ns
    !! total number of 2d slices. This is equal to grid_nz *num_3d_vars + num_2d_vars

  ! custom MPI object types
  !------------------------------------------------------------
  integer :: state_mpi_type
  
contains

  


  !================================================================================
  !================================================================================
  subroutine letkf_mpi_setgrid(nx,ny,ns)
    integer, intent(in) :: nx, ny, ns
    integer :: i, j
    integer :: count, prev

    real, allocatable :: load_weights(:)

    grid_nx = nx
    grid_ny = ny
    grid_ns = ns

    ! determine how many gridpoints this PE should deal with
    ! ----------------------------------------
    if (pe_isroot) print *, ""

    ! determine the process load weights
    ! TODO: remove load weights... this was leftover from a previous failed idea
    allocate(load_weights(pe_size))
    load_weights = 1.0
    load_weights = load_weights / sum(load_weights)


    ! calculate actual numer of gridpoints to use
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
       prev = prev+count
    end do
    deallocate(load_weights)
   
    if (pe_isroot) then
       print *, ""
       print *, "gridpoint assignment balancing:"
       i = minval(scatterv_count)
       j = maxval(scatterv_count)
       if (i==j) then
          print *, " all procs assigned ",i,"gridpoints"
       else
          print *, " all procs assigned between ",i,"and",j,"gridpoints"
       end if
    end if

    ! initialize the custom MPI types
    call init_state_mpi_type()

  end subroutine letkf_mpi_setgrid
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine  letkf_mpi_barrier(syncio)
    integer :: ierr
    logical, optional :: syncio

    call mpi_barrier(mpi_comm_letkf, ierr)    

    if(present(syncio) .and. syncio) then
       flush(output_unit)
       call system('usleep 1')
       call mpi_barrier(mpi_comm_letkf, ierr)
    end if
    
  end subroutine letkf_mpi_barrier
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_ens2ij(ens, ij)
    !! takes ensemble members read in by 1 or more process and distributes
    !! to all processes. Afterward each process will have all ensemble members for a
    !! subset of grid points.

    real, intent(in)    :: ens(grid_nx*grid_ny, grid_ns, size(ens_list))
      !! The ensemble members this given process is responsible for loading in.
      !! ***Shape is ([[letkf_mpi:grid_nx]], [[letkf_mpi:grid_ny]],
      !!    [[letkf_mpi:grid_ny]], [[letkf_mpi:mem]])***

    real, intent(inout) :: ij(mem, grid_ns, ij_count)
      !! The gridpoints this given process is responsible for computing.
      !! ***Shape is ([[letkf_mpi:mem]], [[letkf_mpi:grid_ns]],[[letkf_mpi:ij_count]] )***

#ifdef SCATTER_CUSTOMPACK
    !custom buffer / packing method
    !------------------------------------------------------------
       integer :: ierr, s, i, j, m
       real :: wrk(ij_count)
       real :: wrk2(grid_nx*grid_ny)
    ! TODO, is there any way to do this with fewer mpi_scatter calls?
    do m=1,mem
       do s=1,size(ens,2)
          if(.not. interleave) then
             call mpi_scatterv(ens(:,s,ens_idx(m)), scatterv_count, scatterv_displ, mpi_real, &
                  wrk, ij_count, mpi_real, ens_map(m), mpi_comm_letkf, ierr)
          else
             if(ens_map(m) == pe_rank) then
                do i=1,pe_size
                   do j=1,scatterv_count(i)
                      wrk2(scatterv_displ(i)+j) = ens((j-1)*pe_size+i,s,ens_idx(m))
                   end do
                end do
             end if
             call mpi_scatterv(wrk2, scatterv_count, scatterv_displ, mpi_real, &
                  wrk, ij_count, mpi_real, ens_map(m), mpi_comm_letkf, ierr)
          end if

          ij(m,s,:) = wrk

          if(ierr /= 0) then
             print *, "ERROR: with letkf_mpi_ens2ij",ierr
             stop 1
          end if
       end do
    end do


#else

    ! using MPI custom types method
    !------------------------------------------------------------
    ! TODO, possible cache miss performance issues by specifying data by ns length then nij?
    !  trying toing by nij then x ns?
    integer :: mpitype_col_send
    integer :: mpitype_col_recv
    integer(kind=mpi_address_kind) :: lb, ex, ex_real
    integer :: j, m, p, ierr
    integer :: status(mpi_status_size), recv_req(mem), send_req(pe_size*ens_count)


    ! how big is a real?
    call mpi_type_get_extent(mpi_real, lb, ex_real, ierr)

    ! define the outgoing structure for a single grid column
    call mpi_type_vector(grid_ns, 1, grid_nx*grid_ny, mpi_real, mpitype_col_send, ierr)
    lb=0; ex=ex_real*pe_size
    call mpi_type_create_resized(mpitype_col_send, lb, ex, mpitype_col_send, ierr)
    call mpi_type_commit(mpitype_col_send, ierr)

    ! define the incoming structure for a single grid column
    call mpi_type_vector(grid_ns, 1, mem, mpi_real, mpitype_col_recv, ierr)
    lb=0; ex=ex_real*grid_ns*mem
    call mpi_type_create_resized(mpitype_col_recv, lb, ex, mpitype_col_recv, ierr)
    call mpi_type_commit(mpitype_col_recv, ierr)

    do m=1,mem      
       ! initialize receive
       call mpi_Irecv(ij(m,1,1), scatterv_count(pe_rank+1), mpitype_col_recv, ens_map(m), 0, mpi_comm_letkf, recv_req(m), ierr)
    end do

    ! FAST 
    do m=1,ens_count
       do p=0,pe_size-1
          call mpi_Isend(ens(p+1,1,m), scatterv_count(p+1), mpitype_col_send, p, 0, mpi_comm_letkf, &
              send_req((m-1)*pe_size+p+1), ierr)
       end do
    end do
    ! do m=1,ens_count
    !    do p=0,pe_size-1
    !       call mpi_send(ens(p+1,1,m), scatterv_count(p+1), mpitype_col_send, p, 0, mpi_comm_letkf, ierr)
    !    end do
    ! end do



    ! wait for receive to finish
    call mpi_waitall(mem, recv_req, status, ierr)
    if(j>0) then
       call mpi_waitall(j, send_req, status, ierr)
    end if

    call mpi_type_free(mpitype_col_send, ierr)
    call mpi_type_free(mpitype_col_recv, ierr)

#endif

  end subroutine letkf_mpi_ens2ij
  !================================================================================



  subroutine init_state_mpi_type
!    mpi_type_Create_indexed_block(grid_nx*grid_ny, 1, 
!    state_mpi_type
!    integer :: ierr
!    integer :: state_col_snd_mpi_type

!    call mpi_type_vector(grid_ns, 1, grid_nx*grid_ny, mpi_real, state_col_snd_mpi_type, ierr)
!    call mpi_type_commit(state_col_mpi_type, ierr)

  end subroutine init_state_mpi_type


  !================================================================================
  !================================================================================
  subroutine letkf_mpi_ij2ens(ij, ens)
    !! Takes gridpoints that are scattered across the processors and
    !! combines them and sends individual whole ensemble members to the processors
    !! responsible for saving them

    real, intent(in) :: ij(mem, grid_ns, ij_count)
      !! The gridpoints this given process is responsible for computing.
      !! ***Shape is ([[letkf_mpi:mem]], [[letkf_mpi:grid_ns]],[[letkf_mpi:ij_count]] )***

    real, intent(inout)    :: ens(grid_nx*grid_ny, grid_ns, size(ens_list))
      !! The ensemble members this given process is responsible for loading in.
      !! ***Shape is ([[letkf_mpi:grid_nx]], [[letkf_mpi:grid_ny]],
      !!    [[letkf_mpi:grid_ny]], [[letkf_mpi:mem]])***

    integer :: ierr, s, m, i, j
    real :: wrk(ij_count)
    real :: wrk2(grid_nx*grid_ny)

    ! TODO, is there any way to do this with fewer mpi_gather calls?
    do m=1,mem
       do s=1,size(ens,2)
          wrk = ij(m,s,:)
          if(.not. interleave) then
             call mpi_gatherv(wrk, ij_count, mpi_real,&
                  ens(:,s,ens_idx(m)), scatterv_count, scatterv_displ, mpi_real,&
                  ens_map(m), mpi_comm_letkf, ierr)
          else
             call mpi_gatherv(wrk, ij_count, mpi_real,&
                  wrk2, scatterv_count, scatterv_displ, mpi_real,&
                  ens_map(m), mpi_comm_letkf, ierr)
             if(ens_map(m) == pe_rank) then
                do i=1,pe_size
                   do j=1,scatterv_count(i)
                      ens((j-1)*pe_size+i,s,ens_idx(m)) = wrk2(scatterv_displ(i)+j)
                   end do
                end do
             end if
          end if
          if(ierr /= 0) then
             print *, "ERROR: with letkf_mpi_ij2ens",ierr
             stop 1
          end if
       end do
    end do

  end subroutine letkf_mpi_ij2ens
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_grd2ij_real(grd, ij)
    !! Takes a single grid on the root process and distributes
    !! portions of it to the worker processes

    real, intent(in) :: grd(grid_nx*grid_ny)
      !! The 2D grid to send.
      !! ***Shape is ([[letkf_mpi:grid_nx]], [[letkf_mpi:grid_ny]])***

    real, intent(inout) :: ij(ij_count)
      !! The gridpoints this process is responsible for using.
      !! ***Shape is ([[letkf_mpi:ij_count]])***

    integer :: ierr,i,j
    real :: wrk(grid_nx*grid_ny)

    if(.not.interleave) then
       call mpi_scatterv(grd, scatterv_count, scatterv_displ, mpi_real, &
            ij, ij_count, mpi_real, pe_root, mpi_comm_letkf, ierr)
    else
       if (pe_isroot) then
          do i=1,pe_size
             do j=1,scatterv_count(i)
                wrk(scatterv_displ(i)+j) = grd((j-1)*pe_size+i)
             end do
          end do
       end if
       call mpi_scatterv(wrk, scatterv_count, scatterv_displ, mpi_real, &
            ij, ij_count, mpi_real, pe_root, mpi_comm_letkf, ierr)
    end if

    if(ierr /= 0) then
       print *, "ERROR: with letkf_mpi_grd2ij",ierr
       stop 1
    end if
  end subroutine letkf_mpi_grd2ij_real
  !================================================================================



  !================================================================================
  !================================================================================
  subroutine letkf_mpi_grd2ij_logical(grd, ij)
    !! Takes a single grid on the root process and distributes
    !! portions of it to the worker processes

    logical, intent(in) :: grd(grid_nx*grid_ny)
      !! The 2D grid to send.
      !! ***Shape is ([[letkf_mpi:grid_nx]], [[letkf_mpi:grid_ny]])***

    logical, intent(inout) :: ij(ij_count)
      !! The gridpoints this process is responsible for using.
      !! ***Shape is ([[letkf_mpi:ij_count]])***

    integer :: ierr,i,j
    logical :: wrk(grid_nx*grid_ny)

    if(.not.interleave) then
       call mpi_scatterv(grd, scatterv_count, scatterv_displ, mpi_logical, &
            ij, ij_count, mpi_logical, pe_root, mpi_comm_letkf, ierr)
    else
       if (pe_isroot) then
          do i=1,pe_size
             do j=1,scatterv_count(i)
                wrk(scatterv_displ(i)+j) = grd((j-1)*pe_size+i)
             end do
          end do
       end if
       call mpi_scatterv(wrk, scatterv_count, scatterv_displ, mpi_logical, &
            ij, ij_count, mpi_logical, pe_root, mpi_comm_letkf, ierr)
    end if

    if(ierr /= 0) then
       print *, "ERROR: with letkf_mpi_grd2ij",ierr
       stop 1
    end if
  end subroutine letkf_mpi_grd2ij_logical
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_ij2grd(ij, grd, rank)
    real, intent(in)    :: ij(ij_count)
    real, intent(inout) :: grd(grid_nx* grid_ny)
    integer, intent(in), optional :: rank

    integer :: rank0
    integer :: ierr,i,j
    real :: wrk(grid_nx*grid_ny)

    if(present(rank)) then
       rank0 = rank
    else
       rank0 = pe_root
    end if

    if(.not.interleave) then
       call mpi_gatherv(ij, ij_count, mpi_real, &
            grd, scatterv_count, scatterv_displ, mpi_real, rank0, mpi_comm_letkf, ierr)
    else
       call mpi_gatherv(ij, ij_count, mpi_real, &
            wrk, scatterv_count, scatterv_displ, mpi_real, rank0, mpi_comm_letkf, ierr)
       if(rank0 == pe_rank) then
          do i=1,pe_size
             do j=1,scatterv_count(i)
                grd((j-1)*pe_size+i) = wrk(scatterv_displ(i)+j)
             end do
          end do
       end if
    end if

    if(ierr /= 0) then
       print *, "ERROR: with letkf_mpi_ij2grd", ierr
       stop 1
    end if
  end subroutine letkf_mpi_ij2grd
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_obs(ohx, qc)
    real,     intent(inout)           :: ohx(:,:)
    integer,  intent(inout), optional :: qc(:,:)
    integer :: ierr
    !TODO, this is inefficient
    call mpi_allreduce(mpi_in_place, ohx, mem*size(ohx,2), mpi_real, mpi_sum, mpi_comm_letkf, ierr)
    if(present(qc)) then
       call mpi_allreduce(mpi_in_place, qc, mem*size(ohx,2), mpi_integer, mpi_sum, mpi_comm_letkf, ierr)
    end if
  end subroutine letkf_mpi_obs
  !================================================================================




  !================================================================================
  !================================================================================
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
  !================================================================================




  !================================================================================
  !================================================================================
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
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_init(filename)
    character(len=*), intent(in) :: filename
    integer :: i, j
    integer :: count, prev, unit

    namelist /letkf_mpi/ mem, interleave

    if(pe_isroot) then
       print "(A)", ""
       print "(A)", ""
       print "(A)", '============================================================'
       print "(A)", ' letkf_mpi_init() : '
       print "(A)", '============================================================'
    end if

    ! read in namelist
    open(newunit=unit, file=filename)
    read(unit, nml=letkf_mpi)
    close(unit)
    if (pe_isroot) then
       print letkf_mpi
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
             print '(A5,I4,A,I4,A,I4,A,I4)', "proc", i,' is I/O for ',count,&
                  ' ensemble member(s): ', prev+1, ' to ',count+prev
          else if (count == 1) then
             print '(A5,I4,A,I4)', " proc ", i,' is I/O for    1 ensemble member(s): ', prev+1
          else
             print '(A5,I4,A)', " proc ", i,' is I/O for    0 ensemble member(s)'
          end if
       end if

       prev = prev + count
    end do

#ifdef SCATTER_CUSTOMPACK
    if (pe_isroot) then
       print *,""
       print *,"NOTE: using secondary buffering for scatter /gether calls"
    end if
#endif

  end subroutine letkf_mpi_init
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_final
    integer :: ierr
    ! real :: cur, mem, mem_avg, mem_min, mem_max, mem_sum
    ! ! calculate memore hiwater mark
    ! call getMem(cur, mem)
    ! call mpi_reduce(mem, mem_min, 1, mpi_real, mpi_min, pe_root, mpi_comm_letkf, ierr)
    ! call mpi_reduce(mem, mem_max, 1, mpi_real, mpi_max, pe_root, mpi_comm_letkf, ierr)
    ! call mpi_reduce(mem, mem_sum, 1, mpi_real, mpi_sum, pe_root, mpi_comm_letkf, ierr)
    ! mem_avg = mem_sum / pe_size

    ! if(pe_isroot) then
    !    print *, new_line('a'), &
    !         new_line('a'), "============================================================",&
    !         new_line('a'), " Max memory usage",&
    !         new_line('a'), "============================================================"
    !    print *, "per core:"
    !    print '(A,F5.2,A)', "   avg: ",mem_avg," GB"
    !    print '(A,F5.2,A)', "   min: ",mem_min," GB"
    !    print '(A,F5.2,A)', "   max: ",mem_max," GB"
    !    print *,"total:"
    !    print '(A,F5.2,A)', "    ",mem_sum," GB"
    ! end if
    call mpi_finalize(ierr)
  end subroutine letkf_mpi_final
  !================================================================================





end module letkf_mpi
