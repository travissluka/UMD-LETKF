! Copyright 2016-2019 Travis Sluka
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!================================================================================
!> mpi methods for grid decomposition and io scattering
!================================================================================
MODULE letkf_mpi
  USE mpi
  USE letkf_config

  IMPLICIT NONE
  PRIVATE



  !================================================================================
  !================================================================================
  ! Public module components
  !================================================================================
  !================================================================================

  PUBLIC :: letkf_mpi_preinit
  PUBLIC :: letkf_mpi_init
  PUBLIC :: letkf_mpi_final
  PUBLIC :: letkf_mpi_abort

  PUBLIC :: letkf_mpi_setgrid
  PUBLIC :: letkf_mpi_barrier
  PUBLIC :: letkf_mpi_nextio

  PUBLIC :: letkf_mpi_grd2ij
  INTERFACE letkf_mpi_grd2ij
     MODULE PROCEDURE letkf_mpi_grd2ij_real
     MODULE PROCEDURE letkf_mpi_grd2ij_logical
  END INTERFACE letkf_mpi_grd2ij

  PUBLIC :: letkf_mpi_ij2grd
  INTERFACE letkf_mpi_ij2grd
     MODULE PROCEDURE letkf_mpi_ij2grd_real
  END INTERFACE letkf_mpi_ij2grd


  ! MPI parameters
  !--------------------------------------------------------------------------------
  INTEGER, PUBLIC, PROTECTED :: pe_rank   !< PE of this MPI rank
  INTEGER, PUBLIC, PROTECTED :: pe_size   !< number of MPI PEs
  INTEGER, PUBLIC, PROTECTED :: pe_root   !< the root MPI PE
  LOGICAL, PUBLIC, PROTECTED :: pe_isroot !< true of this PE is root
  INTEGER, PUBLIC, PROTECTED :: letkf_mpi_comm !< the MPI communicator LETKF should use
  INTEGER, PUBLIC, PROTECTED :: ij_count  !< number of grid points this PE is responsible for
  INTEGER, PUBLIC, PROTECTED, ALLOCATABLE :: ij_count_pe(:)  !< number of gridpoints, for each PE
  INTEGER, PUBLIC :: mpitype_grid_nxy_real
  INTEGER, PUBLIC :: mpitype_grid_nk_ns
  INTEGER, PUBLIC :: mpitype_grid_ns

  ! ensemble parameters
  !--------------------------------------------------------------------------------
  INTEGER, PUBLIC, PROTECTED :: ens_size  !< ensemble size




  !================================================================================
  !================================================================================
  ! Private module components
  !================================================================================
  !================================================================================

  ! grid definition
  !--------------------------------------------------------------------------------
  INTEGER :: grid_nx  !< size of the global x dimension
  INTEGER :: grid_ny  !< size of the global y dimension
  INTEGER :: grid_ns  !< size of the vertical/variable dimension (slabs)

  ! MPI scatter/gather parameters
  !--------------------------------------------------------------------------------
  INTEGER :: mpitype_grid_nxy_logical
  INTEGER :: nextio
  INTEGER :: ppn



CONTAINS




  !================================================================================
  !> Initialize MPI (so that all console output is done through the root PE).
  !! Should be done before pretty much every other module
  !! initialization routine. Including before letkf_mpi_init().
  !! TODO - modify this to allow specification of an MPI communicator
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_mpi_preinit()
    INTEGER :: ierr

    ! initialize MPI
    CALL mpi_init(ierr)

    ! determine the communicator
    letkf_mpi_comm = mpi_comm_world
    pe_root = 0

    nextio = 0

    ! other initialization
    CALL mpi_comm_size(letkf_mpi_comm, pe_size, ierr)
    CALL mpi_comm_rank(letkf_mpi_comm, pe_rank, ierr)
    pe_isroot = pe_root == pe_rank

  END SUBROUTINE letkf_mpi_preinit
  !================================================================================



  !================================================================================
  !> initialize the MPI module.
  !! Read in a namelist, determine grid distribution.
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_mpi_init(config)
    TYPE(configuration), INTENT(in) :: config


    ! read in our section of the configuration
    CALL config%get("ens_size", ens_size)
    CALL config%get("ppn", ppn, 1)

    ! output basic MPI info
    IF(pe_isroot) THEN
       PRINT *, ""
       PRINT *, "MPI processors: ", pe_size
       IF (letkf_mpi_comm == mpi_comm_world) THEN
          PRINT *, "MPI communicator: MPI_COMM_WORLD"
       ELSE
          PRINT *, "MPI communicaotr: ", letkf_mpi_comm
       END IF
       PRINT *, ""
       PRINT '(A,I0)', " mpi.ens_size=", ens_size
       PRINT '(A,I0)', " mpi.ppn=", ppn
    ENDIF

  END SUBROUTINE letkf_mpi_init
  !================================================================================



  !================================================================================
  !> Abort the program in a clean way
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_mpi_abort(str)
    CHARACTER(*), INTENT(in) :: str
    INTEGER :: ierr

    PRINT *, ""
    PRINT *, "============================================================"
    PRINT *, "FATAL ERROR: ", str
    PRINT *, "============================================================"
    print *, ""
    CALL mpi_abort(letkf_mpi_comm, 1, ierr)
  END SUBROUTINE letkf_mpi_abort
  !================================================================================



  !================================================================================
  !> Shutdown MPI.
  !! This should be the last thing called by the LETKF library
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_mpi_final
    INTEGER :: ierr

    CALL mpi_finalize(ierr)
  END SUBROUTINE letkf_mpi_final
  !================================================================================



  !================================================================================
  !> Given the grid dimension, determine how to split it up amongst the PEs
  !! for any future scatter/gather calls. This should only be called once, at
  !! initialization time, by the letkf_state module.
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_mpi_setgrid(nx, ny, ns)
    INTEGER, INTENT(in) :: nx, ny, ns
    INTEGER :: prev, cnt, i

    grid_nx = nx
    grid_ny = ny
    grid_ns = ns

    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "setting MPI scatter parameters using..."
       PRINT *, " grid_nx:",grid_nx
       PRINT *, " grid_ny:",grid_ny
       PRINT *, " grid_ns:",grid_ns
       PRINT *, ""
    ENDIF

    ! TODO implement masking

    ! TODO allow for a tiled distribution, not just round robin

    ! Calculate the number of gridpoints to be handled for each proc.
    ! afterward, ij_count and ij_count_pe will be correctly set
    ALLOCATE(ij_count_pe(0:pe_size-1))
    prev = 0
    cnt=NINT(1.0*grid_nx*grid_ny/pe_size)
    DO i=0, pe_size-1
       IF (i== pe_size-1) cnt = grid_nx*grid_ny - prev
       IF( i==pe_rank) THEN
          ij_count = cnt
       END IF
       ij_count_pe(i) = cnt
       prev = prev + cnt
    END DO

    ! initialize the custom MPI types
    CALL mpitypes()

  CONTAINS

    SUBROUTINE mpitypes()
      INTEGER(kind=mpi_address_kind) :: lb, ex, ex_real, ex_logical
      INTEGER :: ierr
      ! TODO, possible cache miss performance issues caused by specifying data by ns
      ! length and then nij?

      ! how big is a real, logical
      CALL mpi_type_get_extent(mpi_real, lb, ex_real, ierr)
      CALL mpi_type_get_extent(mpi_logical, lb, ex_logical, ierr)
      lb=0;

      ! scatter of 2D grid, skipping across PE counts
      ex=ex_real*pe_size
      CALL mpi_type_create_resized(mpi_real, lb, ex, mpitype_grid_nxy_real, ierr)
      CALL mpi_type_commit(mpitype_grid_nxy_real, ierr)

      ex=ex_logical*pe_size
      CALL mpi_type_create_resized(mpi_logical, lb, ex, mpitype_grid_nxy_logical, ierr)
      CALL mpi_type_commit(mpitype_grid_nxy_logical, ierr)

      ! scatter of 2D grid, skipping across PE*
      ex=ex_real*ens_size*grid_ns
      !call mpi_type_vector(grid_ns, 1, grid_nx*grid_ny, mpi_real,&
      !   mpitype_grid_nxy_ns_real, ierr)
      CALL mpi_type_create_resized(mpi_real, lb, ex, mpitype_grid_nk_ns, ierr)
      CALL mpi_type_commit(mpitype_grid_nk_ns, ierr)

      ex=ex_real*grid_ns
      CALL mpi_type_create_resized(mpi_real, lb, ex, mpitype_grid_ns, ierr)
      CALL mpi_type_commit(mpitype_grid_ns, ierr)

    END SUBROUTINE mpitypes
  END SUBROUTINE letkf_mpi_setgrid
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_mpi_ij2grd_real(root, ij, grd)
    INTEGER, INTENT(in) :: root
    REAL, INTENT(in) :: ij(ij_count)
    REAL, INTENT(out) :: grd(grid_nx,grid_ny)

    INTEGER :: p, ierr
    INTEGER :: recvs(pe_size)


    ! initialize asynchronous receives
    IF (pe_rank == root) THEN
       DO p=0,pe_size-1
          CALL mpi_irecv(grd(p+1,1), ij_count_pe(p), mpitype_grid_nxy_real,&
               p, 0, letkf_mpi_comm, recvs(p+1), ierr)
       END DO
    END IF

    ! do the send
    CALL mpi_send(ij, ij_count, mpi_real, root, 0, letkf_mpi_comm, ierr)

    ! wait for the receives to finish
    IF (pe_rank == root) &
         CALL mpi_waitall(pe_size, recvs, MPI_STATUSES_IGNORE, ierr)
  END SUBROUTINE letkf_mpi_ij2grd_real
  !================================================================================



  !================================================================================
  !> Decompose a 2D grid from a single PE and scatter it across the PEs
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_mpi_grd2ij_real(root, grd, ij)
    INTEGER, INTENT(in) :: root
    REAL, ALLOCATABLE, INTENT(in) :: grd(:,:)
    REAL, INTENT(out) :: ij(ij_count)

    INTEGER :: recv_req, send_req(0:pe_size-1), ierr, p

    IF (pe_isroot .AND. .NOT. ALLOCATED(grd)) &
         CALL letkf_mpi_abort("2D grid not allocated on root node")

    ! TODO, change this so that we're using mpi_scatter call instead?

    ! initialize asynchronous receives
    CALL mpi_irecv(ij, ij_count_pe(pe_rank), mpi_real, &
         root, 0, letkf_mpi_comm, recv_req, ierr)

    ! initialize asynchronous sends
    IF ( pe_rank == root) THEN
       DO p=0,pe_size-1
          CALL mpi_Isend(grd(p+1,1), ij_count_pe(p), mpitype_grid_nxy_real, &
               p, 0, letkf_mpi_comm, send_req(p), ierr)
       END DO
    END IF

    ! wait for receves and sends to finish
    CALL mpi_wait(recv_req, MPI_STATUS_IGNORE, ierr)
    IF (pe_rank == root) &
         CALL mpi_waitall(pe_size, send_req, MPI_STATUSES_IGNORE, ierr)

  END SUBROUTINE letkf_mpi_grd2ij_real
  !================================================================================



  !================================================================================
  !> Decompose a 2D grid from a single PE and scatter it across the PEs
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_mpi_grd2ij_logical(root, grd, ij)
    INTEGER, INTENT(in) :: root
    LOGICAL, INTENT(in) :: grd(grid_nx, grid_ny)
    LOGICAL, INTENT(out) :: ij(ij_count)

    INTEGER :: recv_req, send_req(0:pe_size-1), ierr, p

    ! TODO, change this so that we're using mpi_scatter call instead?

    ! initialize asynchronous receives
    CALL mpi_irecv(ij, ij_count_pe(pe_rank), mpi_logical, &
         root, 0, letkf_mpi_comm, recv_req, ierr)

    ! initialize asynchronous sends
    IF ( pe_rank == root) THEN
       DO p=0,pe_size-1
          CALL mpi_Isend(grd(p+1,1), ij_count_pe(p), mpitype_grid_nxy_logical, &
               p, 0, letkf_mpi_comm, send_req(p), ierr)
       END DO
    END IF

    ! wait for receves and sends to finish
    CALL mpi_wait(recv_req, MPI_STATUS_IGNORE, ierr)
    IF (pe_rank == root)&
         CALL mpi_waitall(pe_size, send_req, MPI_STATUSES_IGNORE, ierr)

  END SUBROUTINE letkf_mpi_grd2ij_logical
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_mpi_barrier
    INTEGER :: ierr
    CALL mpi_barrier(letkf_mpi_comm, ierr)

  END SUBROUTINE letkf_mpi_barrier
  !================================================================================



  !================================================================================
  !> Return the next PE that should handle IO. This function distributes the
  !! PEs so that work should be spread out evenly across nodes
  !--------------------------------------------------------------------------------
  FUNCTION letkf_mpi_nextio(forced_pe) RESULT(pe)
    INTEGER, OPTIONAL, INTENT(in) :: forced_pe
    ! TODO IO proc choice should be strided across nodes
    INTEGER :: pe

    pe = MERGE(forced_pe, nextio, PRESENT(forced_pe))
    ! TODO allow option to place an upper bound on this, in case we are
    ! severely memory limited
    nextio = pe+ppn
    IF (nextio >= pe_size) THEN
       nextio = MOD(MOD(nextio, pe_size) + 1, ppn)
    END IF
  END FUNCTION letkf_mpi_nextio
  !================================================================================

END MODULE letkf_mpi
