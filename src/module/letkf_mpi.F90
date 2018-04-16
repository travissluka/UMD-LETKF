module letkf_mpi
  !! MPI handler, handles scattering and gathering of state grids and ensemble members

  use mpi
  use iso_fortran_env

  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: letkf_mpi_preinit, letkf_mpi_init, letkf_mpi_final
  public :: letkf_mpi_setgrid
  public :: letkf_mpi_obs, letkf_mpi_ij2grd
  public :: letkf_mpi_ens2ij, letkf_mpi_ij2ens
  public :: letkf_mpi_barrier

  public :: letkf_mpi_wait, letkf_mpi_wait_recv

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

  integer, public, protected :: ens_count
    !! number of ensemble members for which this process is responsible for the I/O

  integer, public, protected, allocatable :: ens_list(:)
    !! lists which I/O members this proc is responsible for
    !! ***Size is ([[letkf_mpi:ens_count]])***

  integer,public, protected, allocatable :: ens_map(:)
    !! for each index MEM, the number of the PE
    !! respsonsible for its I/O

  integer, public, protected :: ij_count
    !! Number of grid points the current process is responsible for handling

  integer, public, protected :: pe_rank, pe_size, pe_root
  logical, public, protected :: pe_isroot
  integer, public, protected :: mpi_comm_letkf

  integer, public, protected, allocatable :: io_order(:)


  ! private variables
  ! ------------------------------

  logical :: skip_masked = .true.

  integer, allocatable :: ij_list(:)
    !! A list of grid point locations that this process
    !! is responsible for doing core LETKF algorithm on

  integer, allocatable :: scatterv_count(:)
  integer, allocatable :: scatterv_displ(:)

  integer :: grid_nx
    !! number of grid points in the x direction (longitude usually)  

  integer :: grid_ny
    !! number of grid points in the y direction (latitude usually)  

  integer :: grid_ns
    !! total number of 2d slices. This is equal to grid_nz *num_3d_vars + num_2d_vars

  ! custom MPI object types
  integer :: mpitype_grd_nxy_ns
  integer :: mpitype_grd_nxy
  integer :: mpitype_grd_nm_ns_nij
  integer :: mpitype_grd_ns_nij

  ! queued mpi send/recv wait requests
  integer, parameter :: max_wait_size = 1024 ! TODO, remove this hardcoding
  integer :: rcv_wait(max_wait_size)  
  integer :: snd_wait(max_wait_size)
  integer :: rcv_wait_count = 0
  integer :: snd_wait_count = 0
  


contains
  



  !================================================================================
  !================================================================================
  subroutine letkf_mpi_setgrid(nx,ny,ns, mask)
    integer, intent(in) :: nx, ny, ns
    logical, intent(in) :: mask(:,:)
    integer :: i, j
    integer :: gridpoint_count, cnt, prev

    grid_nx = nx
    grid_ny = ny
    grid_ns = ns

    
    ! determine how many gridpoints this PE should deal with
    ! ----------------------------------------
    if (pe_isroot) then
       print *,""
       print *, "Gridpoint/proc distributions..."
    end if


    ! determine the actual number of gridpoints DA is to be done for
    ! if we are using a mask
    gridpoint_count = count(.not. mask)
    if (pe_isroot)  then
       print *, " total    gridpoints: ", grid_nx*grid_ny
       print *, " unmasked gridpoints: ", gridpoint_count
       print *, ""
       if (skip_masked) then
          print *, " will scatter only UNMAKSED gridpoints"
          print *, "NOT YET IMPLEMENTED"
          stop 1
       else
          print *, " will scatter ALL gridpoints"
       end if
    end if
    if (.not. skip_masked) gridpoint_count = grid_nx*grid_ny

    cnt = nint(1.0*gridpoint_count / pe_size)
    
    ! calculate actual numer of gridpoints to use for each proc
    allocate(scatterv_count(pe_size))
    allocate(scatterv_displ(pe_size))
    prev = 0
    do i=0, pe_size-1
       if (i == pe_size-1) cnt = gridpoint_count - prev

       if (i == pe_rank) then
          ij_count = cnt
          allocate(ij_list(cnt))
          do j = 1,cnt
             ij_list(j) = prev + j
          end do
       end if
       scatterv_count(i+1) = cnt
       scatterv_displ(i+1) = prev
       prev = prev+cnt
    end do
   
    if (pe_isroot) then
       print *, ""
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
  subroutine letkf_mpi_wait()
    call letkf_mpi_wait_recv()
    call letkf_mpi_wait_send()
  end subroutine letkf_mpi_wait
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_wait_send()
    integer :: status(mpi_status_size, max_wait_size)
    integer :: ierr
    if (snd_wait_count > 0) then
       call mpi_waitall(snd_wait_count, snd_wait, status, ierr)
       snd_wait_count = 0
    end if
  end subroutine letkf_mpi_wait_send
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_wait_recv()
    integer :: status(mpi_status_size, max_wait_size)
    integer :: ierr
    if (rcv_wait_count > 0) then
       call mpi_waitall(rcv_wait_count, rcv_wait, status, ierr)
       rcv_wait_count = 0
    end if
  end subroutine letkf_mpi_wait_recv
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

    integer :: m, p, ierr
    integer :: status(mpi_status_size, max_wait_size), recv_req(mem), send_req(pe_size*ens_count)

    !------------------------------------------------------------
    ! TODO, possible cache miss performance issues by specifying data by ns length then nij?
    !  trying toing by nij then x ns?

    ! initialize asynchronous receives for our segment of  each ensemble member
    do m=1,mem      
       call mpi_Irecv(ij(m,1,1), scatterv_count(pe_rank+1), mpitype_grd_nm_ns_nij,&
            ens_map(m), 0, mpi_comm_letkf, recv_req(m), ierr)
    end do

    ! initialize asynchronous scatters of each whole ensemble we have stored
    do m=1,ens_count
       do p=0,pe_size-1
          call mpi_Isend(ens(p+1,1,m), scatterv_count(p+1), mpitype_grd_nxy_ns,&
               p, 0, mpi_comm_letkf, send_req((m-1)*pe_size+p+1), ierr)
       end do
    end do

    ! wait for receives to finish
    call mpi_waitall(mem, recv_req, status, ierr)

    ! wait for the sends to finish
    if(ens_count>0) call mpi_waitall(pe_size*ens_count, send_req, status, ierr)

  end subroutine letkf_mpi_ens2ij
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine init_state_mpi_type
    integer(kind=mpi_address_kind) :: lb, ex, ex_real
    integer :: ierr

    ! TODO, possible cache miss performance issues by specifying data by ns length then nij?

    ! how big is a real?
    call mpi_type_get_extent(mpi_real, lb, ex_real, ierr)


    call mpi_type_vector(grid_ns, 1, grid_nx*grid_ny, mpi_real, mpitype_grd_nxy_ns, ierr)
    lb=0; ex=ex_real*pe_size
    call mpi_type_create_resized(mpitype_grd_nxy_ns, lb, ex, mpitype_grd_nxy_ns, ierr)
    call mpi_type_commit(mpitype_grd_nxy_ns, ierr)


    call mpi_type_vector(grid_ns, 1, mem, mpi_real, mpitype_grd_nm_ns_nij, ierr)
    lb=0; ex=ex_real*grid_ns*mem
    call mpi_type_create_resized(mpitype_grd_nm_ns_nij, lb, ex, mpitype_grd_nm_ns_nij, ierr)
    call mpi_type_commit(mpitype_grd_nm_ns_nij, ierr)


    lb=0; ex=ex_real*pe_size
    call mpi_type_create_resized(mpi_real, lb, ex, mpitype_grd_nxy, ierr)
    call mpi_type_commit(mpitype_grd_nxy, ierr)


    call mpi_type_vector(grid_ns, 1, 1, mpi_real, mpitype_grd_ns_nij, ierr)
    call mpi_type_commit(mpitype_grd_ns_nij, ierr)

  end subroutine init_state_mpi_type
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_ij2ens(ij, ens, async)
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

    logical, intent(in), optional :: async

    logical :: async0
    integer ::m, p, ierr
    integer :: status(mpi_status_size, max_wait_size), send_req(mem), recv_req(pe_size*ens_count)

    async0 = merge(async, .false., async)

    ! initialize the asynchronous receives for our whole ensembme members
    do m=1,ens_count
       do p=0,pe_size-1
          call mpi_Irecv(ens(p+1,1,m), scatterv_count(p+1), mpitype_grd_nxy_ns, &
               p, 0, mpi_comm_letkf, recv_req((m-1)*pe_size+p+1), ierr)
       end do
    end do

    !initialize the asynchronous sends
    do m=1,mem
       call mpi_Isend(ij(m,1,1), scatterv_count(pe_rank+1), mpitype_grd_nm_ns_nij, &
            ens_map(m), 0, mpi_comm_letkf, send_req(m), ierr)
    end do

    if(async0) then
       ! we will wait somewhere OUTSIDE this subroutine for the send/recv to finish
       if(ens_count>0) then
          rcv_wait(rcv_wait_count+1:rcv_wait_count+(pe_size*ens_count)) = recv_req(:)
          rcv_wait_count = rcv_wait_count + pe_size*ens_count
       end if
       snd_wait(snd_wait_count+1:snd_wait_count+mem) = send_req(:)
       snd_wait_count = snd_wait_count + mem
    else
       ! wait here for the send/recv to finish
       if(ens_count>0) call mpi_waitall(pe_size*ens_count, recv_req, status, ierr)
       call mpi_waitall(mem, send_req, status, ierr)
    end if

  end subroutine letkf_mpi_ij2ens
  !================================================================================



  
  !================================================================================
  !================================================================================
  subroutine letkf_mpi_grd2ij_real(grd, ij, root)
    !! Takes a single grid on the root process and distributes
    !! portions of it to the worker processes

    real, intent(in) :: grd(grid_nx*grid_ny)
      !! The 2D grid to send.
      !! ***Shape is ([[letkf_mpi:grid_nx]], [[letkf_mpi:grid_ny]])***

    real, intent(inout) :: ij(ij_count)
      !! The gridpoints this process is responsible for using.
      !! ***Shape is ([[letkf_mpi:ij_count]])***

    integer, intent(in), optional :: root
    
    integer :: recv_req, send_req(pe_size)
    integer :: status(mpi_status_size, pe_size)

    integer :: root0, p, ierr
    
    root0 = merge(root, pe_root, present(root))

    ! initialize asynchronous receives
    call mpi_Irecv(ij, scatterv_count(pe_rank+1), mpi_real,&
         root0, 0, mpi_comm_letkf, recv_req, ierr)
         
    ! initialize the asynchronous send
    if(root0 == pe_rank) then
       do p=0,pe_size-1
          call mpi_Isend(grd(p+1), scatterv_count(p+1), mpitype_grd_nxy,&
               p, 0, mpi_comm_letkf, send_req(p+1), ierr)
       end do
    end if

    ! wait for receives to finish
    call mpi_wait(recv_req, status(:,1), ierr)
    
    ! wait for sends to finish
    if( root0 ==pe_rank) call mpi_waitall(pe_size, send_req, status, ierr)


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

    if (pe_isroot) then
       do i=1,pe_size
          do j=1,scatterv_count(i)
             wrk(scatterv_displ(i)+j) = grd((j-1)*pe_size+i)
          end do
       end do
    end if
    call mpi_scatterv(wrk, scatterv_count, scatterv_displ, mpi_logical, &
         ij, ij_count, mpi_logical, pe_root, mpi_comm_letkf, ierr)

    if(ierr /= 0) then
       print *, "ERROR: with letkf_mpi_grd2ij",ierr
       stop 1
    end if
  end subroutine letkf_mpi_grd2ij_logical
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_ij2grd(ij, grd, rank, async)
    real, intent(in)    :: ij(grid_ns, ij_count)
    real, intent(out)   :: grd(grid_nx* grid_ny, grid_ns)
    integer, intent(in), optional :: rank
    logical, intent(in), optional :: async

    logical :: async0
    integer :: rank0
    integer :: recv_req(pe_size), send_req
    integer :: status(mpi_status_size, max_wait_size)
    integer :: p, ierr

    async0=merge(async, .false., present(async))
    rank0=merge(rank, pe_root,present(rank))

    ! initialize asynchronous receives
    if(rank0 == pe_rank) then
       grd=0
       do p=0, pe_size-1
          call mpi_Irecv(grd(p+1,1), scatterv_count(p+1), mpitype_grd_nxy_ns, &
               p, 0, mpi_comm_letkf, recv_req(p+1), ierr)
       end do
    end if

    ! initialize the asynchronous send
    call mpi_Isend(ij, scatterv_count(pe_rank+1), mpitype_grd_ns_nij, &
         rank0, 0, mpi_comm_letkf, send_req, ierr)


    if(async0) then
       ! we will wait somewhere OUTSIDE this subroutine for the send/recv to finish
       if(rank0 == pe_rank) then
          rcv_wait(rcv_wait_count+1:rcv_wait_count+pe_size) = recv_req(:)
          rcv_wait_count = rcv_wait_count+pe_size
       end if
       snd_wait_count = snd_wait_count + 1
       snd_wait(snd_wait_count) = send_req
    else
       ! wait here for the send /recv to finish
       if(rank0 == pe_rank) call mpi_waitall(pe_size, recv_req, status, ierr)
       call mpi_wait(send_req, status(:,1), ierr)
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

    rcv_wait_count = 0
    snd_wait_count = 0
  end subroutine letkf_mpi_preinit
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_init(filename)
    character(len=*), intent(in) :: filename
    integer :: i, j
    integer :: count, prev, unit
    integer :: counts(pe_size)
    character(len=1024) :: c
    integer :: ppn, n

    namelist /letkf_mpi/ mem,  ppn, skip_masked

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

    ! determine the order that we should choose PE numbers for IO
    allocate(io_order(pe_size))
    n=((pe_size-1)/ppn)+1
    do i=0,pe_size
       j=mod(i,n)
       io_order(i+1) = j*ppn+(i-j)/n
    end do


    ! determine which ensemble members number this PE should deal with
    ! ----------------------------------------
    counts(:)=0   
    do i=0, mem-1
       j=io_order(mod(i, pe_size)+1)+1
       counts(j)=counts(j)+1
    end do

    count = 0
    prev = 0
    if (pe_isroot) then
       print *, ""
       print *, "ensemble member I/O list:"
    end if
    do i=0, pe_size-1
       count = counts(i+1)

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
       ! print out info about each PE
       ! TODO, this is temporary
       if (pe_rank == i) then
          call hostnm(c)

          if (count > 1) then
             print '(A5,I4,3A,I4,A,I4,A,I4)', "proc", i,' on ',trim(c),' is I/O for ',count,&
                   ' ensemble member(s): ', prev+1, ' to ',count+prev
          else if (count == 1) then
             print '(A5,I4,3A,I4)', " proc ", i,' on ',trim(c),' is I/O for    1 ensemble member:    ', prev+1
          end if
       end if
       call letkf_mpi_barrier(.true.)

       prev = prev + count
    end do
    
  end subroutine letkf_mpi_init
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_mpi_final
    integer :: ierr
  
    call mpi_finalize(ierr)
  end subroutine letkf_mpi_final
  !================================================================================





end module letkf_mpi
