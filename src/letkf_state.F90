MODULE letkf_state
  USE timing
  USE mpi
  USE letkf_mpi

  IMPLICIT NONE
  PRIVATE


  !-----------------------------------------------------------------------------
  ! Public methods
  !-----------------------------------------------------------------------------
  PUBLIC :: letkf_state_register
  PUBLIC :: letkf_state_init
  PUBLIC :: letkf_state_write_ens
  PUBLIC :: letkf_state_write_meansprd
  public :: letkf_state_var_getdef

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------

  !> horizontal grid spec
  !! TODO, remove the nominal lat/lon and keep only in child IO class?
  TYPE, PUBLIC :: letkf_hzgrid_spec
     CHARACTER(len=20)         :: name
     REAL,         ALLOCATABLE :: lat(:,:)
     REAL,         ALLOCATABLE :: lon(:,:)
     real,         ALLOCATABLE :: lat_nom(:) 
     real,         ALLOCATABLE :: lon_nom(:)
     LOGICAL,      ALLOCATABLE :: mask(:,:)
  END TYPE letkf_hzgrid_spec


  !> veritical grid specifications
  !> TODO remove nominal vertical grid and keep only in child IO class?
  ! TODO need to scatter the vertical grids the same as the horizontal
  TYPE, PUBLIC :: letkf_vtgrid_spec
     CHARACTER(len=20)         :: name
     INTEGER                   :: dims !< 1 or 3
     REAL,         ALLOCATABLE :: vert(:,:,:)
     real,         ALLOCATABLE :: vert_nom(:)
  END TYPE letkf_vtgrid_spec


  !> state variable specifications
  TYPE, PUBLIC :: letkf_statevar_spec
     INTEGER           :: grid_s_idx
     INTEGER           :: levels
     CHARACTER(len=20) :: name
     CHARACTER(len=20) :: hzgrid
     CHARACTER(len=20) :: vtgrid
  END TYPE letkf_statevar_spec


  !> Abstract base class for model state reading and writing.
  !! All user-defined and built-in state file I/O classes for
  !! specific file types should inherit this class
  TYPE, ABSTRACT, PUBLIC :: letkf_stateio
     logical :: verbose = .false.
   CONTAINS
     PROCEDURE(I_letkf_stateio_getstr), NOPASS, DEFERRED :: name
     PROCEDURE(I_letkf_stateio_getstr), NOPASS, DEFERRED :: desc
     PROCEDURE(I_letkf_stateio_init),           DEFERRED :: init
     PROCEDURE(I_letkf_stateio_read_specs),     DEFERRED :: read_specs
     PROCEDURE(I_letkf_stateio_read_state),     DEFERRED :: read_state
     PROCEDURE(I_letkf_stateio_write_state),    DEFERRED :: write_state
     PROCEDURE(I_letkf_stateio_write_init),     DEFERRED :: write_init
  END TYPE letkf_stateio

  INTEGER, PUBLIC, PARAMETER :: ENS_BKG_MEAN = -1
  INTEGER, PUBLIC, PARAMETER :: ENS_ANA_MEAN = -2
  INTEGER, PUBLIC, PARAMETER :: ENS_BKG_SPRD = -3
  INTEGER, PUBLIC, PARAMETER :: ENS_ANA_SPRD = -4

  ABSTRACT INTERFACE
     FUNCTION I_letkf_stateio_getstr()
       CHARACTER(:), ALLOCATABLE :: I_letkf_stateio_getstr
     END FUNCTION I_letkf_stateio_getstr

     SUBROUTINE I_letkf_stateio_init(self, nml_filename)
       IMPORT letkf_stateio
       CLASS(letkf_stateio)                  :: self
       CHARACTER(:), ALLOCATABLE, INTENT(in) :: nml_filename
     END SUBROUTINE I_letkf_stateio_init

     SUBROUTINE I_letkf_stateio_read_specs(self, hzgrids, vtgrids, statevars)
       IMPORT letkf_stateio, letkf_hzgrid_spec
       IMPORT letkf_vtgrid_spec, letkf_statevar_spec
       CLASS(letkf_stateio)                               :: self
       TYPE(letkf_hzgrid_spec),  ALLOCATABLE, INTENT(out) :: hzgrids(:)
       TYPE(letkf_vtgrid_spec),  ALLOCATABLE, INTENT(out) :: vtgrids(:)
       TYPE(letkf_statevar_spec),ALLOCATABLE, INTENT(out) :: statevars(:)
     END SUBROUTINE I_letkf_stateio_read_specs

     SUBROUTINE I_letkf_stateio_read_state(self, ensmem, state_var, state_val)
       IMPORT letkf_stateio
       CLASS(letkf_stateio)           :: self
       INTEGER,           INTENT(in)  :: ensmem
       CHARACTER(*),      INTENT(in)  :: state_var
       REAL, ALLOCATABLE, INTENT(out) :: state_val(:,:,:)
     END SUBROUTINE I_letkf_stateio_read_state

     SUBROUTINE I_letkf_stateio_write_init(self, ftype, ensmem)
       IMPORT letkf_stateio
       CLASS(letkf_stateio) :: self
       CHARACTER(len=*), INTENT(in) :: ftype
       INTEGER, INTENT(in) :: ensmem
     END SUBROUTINE I_letkf_stateio_write_init

     SUBROUTINE I_letkf_stateio_write_state(self, ftype, ensmem, state_var, state_val)
       IMPORT letkf_stateio
       CLASS(letkf_stateio)      :: self
       CHARACTER(len=*), INTENT(in) :: ftype
       INTEGER,       INTENT(in) :: ensmem
       CHARACTER(*),  INTENT(in) :: state_var
       REAL,          INTENT(in) :: state_val(:,:,:)
     END SUBROUTINE I_letkf_stateio_write_state

  END INTERFACE


  !-----------------------------------------------------------------------------
  ! Public variables
  !-----------------------------------------------------------------------------
  ! grid definition, after being distributed across PEs
  REAL,    PUBLIC, PROTECTED, ALLOCATABLE :: lat_ij(:)
  REAL,    PUBLIC, PROTECTED, ALLOCATABLE :: lon_ij(:)
  LOGICAL, PUBLIC, PROTECTED, ALLOCATABLE :: mask_ij(:)

  ! state variables after being distributed across PEs
  REAL, PUBLIC, ALLOCATABLE :: state_ij(:,:,:)
  REAL, PUBLIC, ALLOCATABLE :: state_mean_ij(:,:)
  REAL, PUBLIC, ALLOCATABLE :: state_sprd_ij(:,:)

  INTEGER, PUBLIC, PROTECTED :: grid_nx
  INTEGER, PUBLIC, PROTECTED :: grid_ny
  INTEGER, PUBLIC, PROTECTED :: grid_ns

  ! master grid and variable definitions
  TYPE(letkf_hzgrid_spec),  PUBLIC, PROTECTED, ALLOCATABLE :: hzgrids(:)
  TYPE(letkf_vtgrid_spec),  PUBLIC, PROTECTED, ALLOCATABLE :: vtgrids(:)
  TYPE(letkf_statevar_spec),PUBLIC, PROTECTED, ALLOCATABLE :: statevars(:)

  !-----------------------------------------------------------------------------
  ! Private methods
  !-----------------------------------------------------------------------------
  ! only needed here temporarily for my own sanity



  !-----------------------------------------------------------------------------
  ! Private types
  !-----------------------------------------------------------------------------

  !> simple wrapper for letkf_stateio so that we can have an array of pointers
  !! to abstract classes
  TYPE stateio_ptr
     CLASS(letkf_stateio), POINTER :: p
  END TYPE stateio_ptr


  !-----------------------------------------------------------------------------
  ! Private variables
  !-----------------------------------------------------------------------------


  ! registration of built-in and user-defined stateio classes
  INTEGER, PARAMETER :: stateio_reg_max = 100         !< max classes to handle
  INTEGER            :: stateio_reg_num = 0           !< current number of classes registered
  TYPE(stateio_ptr)  :: stateio_reg(stateio_reg_max)  !< list of registered classes
  CLASS(letkf_stateio), POINTER :: stateio_class      !< The actual I/O module selected



  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------

CONTAINS



  !--------------------------------------------------------------------------------
  !> Initialize the letkf state module.
  !! The namefile is read, stateio class is initialized (reading in the grid specs),
  !! ensemble background state read in in parallel, mean and spread are calculated
  !! and saved to a file.
  SUBROUTINE letkf_state_init(nml_filename)
    CHARACTER(:), ALLOCATABLE, INTENT(in) :: nml_filename
    INTEGER :: unit, i,  s, n
    CHARACTER(:), ALLOCATABLE :: ioclass
    LOGICAL :: write_bkg_meansprd = .TRUE.
    logical :: verbose = .false.
    
    NAMELIST /letkf_state/ ioclass, write_bkg_meansprd, verbose

    ! print header
    IF (pe_isroot) THEN
       PRINT '(//A)', ""
       PRINT *, "============================================================"
       PRINT *, " letkf_state_init() : Model state I/O module initialization"
       PRINT *, "============================================================"
       PRINT *, ""
    END IF

    ! read our section of the namelist
    ALLOCATE(CHARACTER(1024) :: ioclass);   WRITE (ioclass, *) "UNDEFINED"
    OPEN(newunit=unit, file=nml_filename, status='OLD')
    READ(unit, nml=letkf_state)
    CLOSE(unit)
    ioclass=toupper(TRIM(ioclass))
    IF (pe_isroot) PRINT letkf_state

    ! print a list of all stateio classes that have been registered
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "List of stateio classes registered:"
       DO i=1, stateio_reg_num
          PRINT *, " * ", toupper(stateio_reg(i)%p%name()), &
               "  (", stateio_reg(i)%p%desc(), ")"
       END DO
       PRINT *, ""
    END IF

    ! determine the stateio class to use
    NULLIFY(stateio_class)
    DO i=1, stateio_reg_num
       IF (stateio_reg(i)%p%name() == ioclass) THEN
          stateio_class => stateio_reg(i)%p
          EXIT
       END IF
    END DO
    IF (.NOT. ASSOCIATED(stateio_class)) THEN
       CALL letkf_mpi_abort("stateio class "//ioclass//" not found.")
    END IF

    ! initialize the stateio class
    stateio_class%verbose = verbose
    CALL stateio_class%init(nml_filename)

    ! Read in the grid /variable specs
    CALL letkf_state_init_spec()
    
    ! allocate the memory needed by this pe for the ensemble model state
    ALLOCATE( state_ij( ens_size, grid_ns, ij_count))
    ALLOCATE( state_mean_ij(      grid_ns, ij_count))
    ALLOCATE( state_sprd_ij(      grid_ns, ij_count))

    ! load in the ensemble and scatter to the PEs
    CALL letkf_state_read_ens()


    CALL timing_start("bkg_mean_sprd")

    ! calculate state mean and spread, and remove the mean
    ! from the perturbations. state_ij will then be the ensemble background
    ! perturbations
    CALL timing_start("calc")
    state_mean_ij = SUM(state_ij, 1) / ens_size
    state_sprd_ij = 0
    DO s=1,grid_ns
       DO n=1,ij_count
          state_ij(:,s,n) = state_ij(:,s,n) - state_mean_ij(s,n)
          state_sprd_ij(s,n) = SQRT(&
               dot_PRODUCT(state_ij(:,s,n),state_ij(:,s,n))/(ens_size-1))
       END DO
    END DO
    CALL timing_stop("calc")

    ! Writ the mean and spread out to files
    IF (write_bkg_meansprd) &
         CALL letkf_state_write_meansprd("bkg")

    CALL timing_stop("bkg_mean_sprd")

  END SUBROUTINE letkf_state_init



  !-----------------------------------------------------------------------------
  !> Writes the ensemble mean and spread out to a set of files.
  !! This subroutine handles both the analysis and background cases, depending
  !! on the value of the argument. The current contents of "state_mean_ij" and
  !! "state_sprd_ij" are used for the output.
  SUBROUTINE letkf_state_write_meansprd(mode)
    CHARACTER(*), INTENT(in) :: mode !< either "ana" or "bkg"
    INTEGER, ALLOCATABLE :: sends(:)
    INTEGER, ALLOCATABLE :: recvs(:)
    INTEGER :: sends_num, recvs_num

    REAL,    ALLOCATABLE :: tmp_r_3d(:,:,:)
    INTEGER :: dest, ierr,  i, j, k, p, s, tag
    INTEGER :: mode_mean, mode_sprd

    ! check input mode argument
    IF (mode == "bkg" ) THEN
       mode_mean = ENS_BKG_MEAN
       mode_sprd = ENS_BKG_SPRD
    ELSE IF (mode == "ana") THEN
       mode_mean = ENS_ANA_MEAN
       mode_sprd = ENS_ANA_SPRD
    ELSE
       IF(pe_isroot) CALL letkf_mpi_abort("Invalid mode: "//mode)
    END IF


    CALL timing_start("write_meansprd")

    ! initialize the background mean/spread output files
    IF (pe_isroot) PRINT '(/,X,A)', "Writing background mean/spread..."
    CALL timing_start("io_init", TIMER_SYNC)
    IF (letkf_mpi_nextio() == pe_rank)&
         CALL stateio_class%write_init(mode, mode_mean)
    IF (letkf_mpi_nextio() == pe_rank)&
         CALL stateio_class%write_init(mode, mode_sprd)
    CALL timing_stop("io_init")
    CALL letkf_mpi_barrier()


    ! initialize all the sends
    ALLOCATE(sends(grid_ns*2))
    sends_num=0
    dest=letkf_mpi_nextio(0)
    DO j=1,2 ! one loop for "mean", one for "sprd"
       s=0
       dest = letkf_mpi_nextio()
       DO i=1,SIZE(statevars) ! for each state variable
          DO k=1,statevars(i)%levels
             s=s+1
             tag=k-1+statevars(i)%grid_s_idx + (j-1)*grid_ns
             sends_num = sends_num + 1
             IF (j==1) THEN
                CALL MPI_ISend(state_mean_ij(s,1), ij_count, mpitype_grid_ns, &
                     dest, tag, letkf_mpi_comm, sends(sends_num), ierr)
             ELSE
                CALL MPI_ISend(state_sprd_ij(s,1), ij_count, mpitype_grid_ns, &
                     dest, tag, letkf_mpi_comm, sends(sends_num), ierr)
             END IF
          END DO
       END DO
    END DO


    ! initialize all the receives
    ALLOCATE(recvs(grid_ns*pe_size))
    dest = letkf_mpi_nextio(0)
    DO j=1,2 ! one loop for "mean", one for "sprd"
       dest = letkf_mpi_nextio()
       DO i=1,SIZE(statevars) ! for each state variable

          CALL timing_start("io_write")
          CALL timing_stop("io_write")
          CALL timing_start("mpi_recv")
          CALL timing_stop("mpi_recv")

          IF (dest == pe_rank) THEN
             CALL timing_start("mpi_recv")
             IF(ALLOCATED(tmp_r_3d)) DEALLOCATE(tmp_r_3d)
             ALLOCATE(tmp_r_3d(grid_nx,grid_ny,statevars(i)%levels))
             recvs_num=0
             ! initialize the receives
             DO p=0,pe_size-1
                DO k=1,statevars(i)%levels
                   recvs_num = recvs_num + 1
                   tag=k-1+statevars(i)%grid_s_idx + (j-1)*grid_ns
                   CALL mpi_Irecv(tmp_r_3d(p+1,1,k), ij_count_pe(p),&
                        mpitype_grid_nxy_real, p, tag, letkf_mpi_comm,&
                        recvs(recvs_num),ierr)
                END DO
             END DO

             ! wait for the receives to finish
             CALL mpi_waitall(recvs_num, recvs, MPI_STATUSES_IGNORE, ierr)
             CALL timing_stop("mpi_recv")

             ! write out to file
             CALL timing_start("io_write")
             if (stateio_class%verbose) &
                  PRINT '(X,A,I0.3,4A)', " PE ", pe_rank, " writing ", &
                  TRIM(statevars(i)%name), " for ",TRIM(mode)// &
                  MERGE(" mean  "," spread", j==1)
             CALL stateio_class%write_state(mode, &
                  MERGE(ENS_BKG_MEAN, ENS_BKG_SPRD,j==1),&
                  statevars(i)%name, tmp_r_3d(:,:,:))
             CALL timing_stop("io_write")

          END IF
       END DO
    END DO

    ! wait for the sends to finish
    CALL mpi_waitall(sends_num, sends, MPI_STATUSES_IGNORE, ierr)
    CALL timing_stop("write_meansprd")

  END SUBROUTINE letkf_state_write_meansprd




  !-----------------------------------------------------------------------------
  !> Reads in the grid specification (horizontal and vertical grids) and the
  !! specification of what state variables are involved with the state IO
  SUBROUTINE letkf_state_init_spec()
    INTEGER :: i, j, ierr

    CALL timing_start("read_state_specs")


    ! Load the grid/state specification
    IF (pe_isroot) THEN

       ! TODO currently only handling the first hzgrid, and first vtgrid,
       ! modify to handle an arbitrary number
       !TODO, no point keeping the separate lat,lon,mask variables around
       CALL stateio_class%read_specs(hzgrids, vtgrids, statevars)

       ! sanity checks to make sure the data are okay
       IF (.NOT. ALLOCATED(hzgrids(1)%lat) ) &
            CALL letkf_mpi_abort("LAT uninitiallized from stateio_class%read_specs()")
       IF (.NOT. ALLOCATED(hzgrids(1)%lon) ) &
            CALL letkf_mpi_abort("LON uninitiallized from stateio_class%read_specs()")
       IF (.NOT. ALLOCATED(hzgrids(1)%mask) ) &
            CALL letkf_mpi_abort("MASK uninitiallized from stateio_class%read_specs()")
       DO i=1,2
          IF ((SIZE(hzgrids(1)%lat, dim=i)/=SIZE(hzgrids(1)%lon,dim=i)) .OR. &
               (SIZE(hzgrids(1)%lat,dim=i)/=SIZE(hzgrids(1)%mask, dim=i)))  THEN
             CALL letkf_mpi_abort("LAT,LON,MASK must have same shapes "//&
                  "from stateio_class%read_specs()")
          END IF
       END DO

       ! TODO make sure all the names are uppercase
       
       ! determine the number of slabs required, and the slab start/stop
       ! point of each variable
       grid_ns = 0
       DO i=1,SIZE(statevars)
          ! determine the starting and stopping slab index
          statevars(i)%grid_s_idx=grid_ns + 1
          DO j = 1, SIZE(vtgrids)
             IF (vtgrids(j)%name == statevars(i)%vtgrid) THEN
                statevars(i)%levels = SIZE(vtgrids(j)%vert,1)
                EXIT
             END IF
          END DO
          IF (j > SIZE(vtgrids)) THEN
             call letkf_mpi_abort('vertical grid "'//trim(statevars(i)%vtgrid)//&
                  '" not found. Needed for "'//trim(statevars(i)%name)//'".')
             STOP 1
          END IF
          grid_ns = grid_ns + statevars(i)%levels
       END DO

       ! get the overall horizontal grid dimensions
       ! TODO use the actually specified hz and vt grids
       grid_nx = SIZE(hzgrids(1)%lat,dim=1)
       grid_ny = SIZE(hzgrids(1)%lat,dim=2)
    END IF


    ! broadcast general grid information to all PEs
    CALL mpi_bcast(grid_nx, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
    CALL mpi_bcast(grid_ny, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
    CALL mpi_bcast(grid_ns, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)

    ! tell the mpi module about the grid layout
    CALL letkf_mpi_setgrid(grid_nx, grid_ny, grid_ns)

    
    ! broadcast horizontal grid spec to all mpi procs
    IF(pe_isroot) i=SIZE(hzgrids)
    CALL mpi_bcast(i, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
    IF (.NOT. pe_isroot) ALLOCATE(hzgrids(i))
    DO i=1,SIZE(hzgrids)
       CALL mpi_bcast(hzgrids(i)%name, 20, mpi_character, pe_root, &
            letkf_mpi_comm, ierr)
       
       ! send nominal latitude
       if(pe_isroot) j=size(hzgrids(i)%lat_nom)
       call mpi_bcast(j, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
       if(.not. pe_isroot) allocate(hzgrids(i)%lat_nom(j))
       call mpi_bcast(hzgrids(i)%lat_nom, j, mpi_real, pe_root, letkf_mpi_comm, ierr)
       
       ! send nominal longitude
       if(pe_isroot) j=size(hzgrids(i)%lon_nom)
       call mpi_bcast(j, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
       if(.not. pe_isroot) allocate(hzgrids(i)%lon_nom(j))
       call mpi_bcast(hzgrids(i)%lon_nom, j, mpi_real, pe_root, letkf_mpi_comm, ierr)
       
       ! TODO, initializing with 0 size on non root nodes to silence gfortran
       ! runtime debug errors. Do I ever need the full 2d grid on other PEs?
       IF (.NOT. pe_isroot) THEN
          ALLOCATE(hzgrids(i)%lat(0,0))
          ALLOCATE(hzgrids(i)%lon(0,0))
          ALLOCATE(hzgrids(i)%mask(0,0))
       END IF       
    END DO

    
    ! scatter horizontal grid parameters
    ALLOCATE(lat_ij(ij_count))
    ALLOCATE(lon_ij(ij_count))
    ALLOCATE(mask_ij(ij_count))
    CALL letkf_mpi_grd2ij(pe_root, hzgrids(1)%lat, lat_ij)
    CALL letkf_mpi_grd2ij(pe_root, hzgrids(1)%lon, lon_ij)
    CALL letkf_mpi_grd2ij(pe_root, hzgrids(1)%mask, mask_ij)

    
    ! print summary of horizontal grids
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "Horizontal grids summary"
       PRINT *, "----------------------------------------"
       DO i=1,SIZE(hzgrids)
          PRINT *, "Name: ", hzgrids(i)%name
          PRINT "(X,A,I12)", " grid_nx:    ", grid_nx
          PRINT "(X,A,I12)", " grid_ny:    ", grid_ny
          PRINT "(X,A,I12)", " grid total: ", grid_nx*grid_ny
          PRINT "(X,A,I12,A,F4.1,A)", " grid masked:", COUNT(hzgrids(i)%mask), &
               " (",REAL(COUNT(hzgrids(i)%mask)*100) / REAL(grid_nx*grid_ny),"%)"
          PRINT *, " lat range: ", MINVAL(hzgrids(i)%lat), MAXVAL(hzgrids(i)%lat)
          PRINT *, " lon range: ", MINVAL(hzgrids(i)%lon), MAXVAL(hzgrids(i)%lon)
          PRINT *, ""
       END DO
    END IF


    ! broadcast vertical grid specs to all mpi PEs
    IF (pe_isroot) i=SIZE(vtgrids)
    CALL mpi_bcast(i, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
    IF (.NOT. pe_isroot) ALLOCATE(vtgrids(i))
    DO i=1, SIZE(vtgrids)
       CALL mpi_bcast(vtgrids(i)%name, 20, mpi_character, pe_root, letkf_mpi_comm, ierr)
       CALL mpi_bcast(vtgrids(i)%dims, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
       !TODO handle 3D vertical grid definitions
       ! send full 1D/2D/3D vertical grid
       IF (vtgrids(i)%dims /= 1) CALL letkf_mpi_abort("only 1D vertical grids handled now")
       IF (pe_isroot) j=SIZE(vtgrids(i)%vert)
       CALL mpi_bcast(j, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
       IF (.NOT. pe_isroot) ALLOCATE(vtgrids(i)%vert(j,1,1))
       CALL mpi_bcast(vtgrids(i)%vert, j, mpi_real, pe_root, letkf_mpi_comm, ierr)

       ! send nominal vertical grid
       ! send nominal longitude
       if(pe_isroot) j=size(vtgrids(i)%vert_nom)
       call mpi_bcast(j, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
       if(.not. pe_isroot) allocate(vtgrids(i)%vert_nom(j))
       call mpi_bcast(vtgrids(i)%vert_nom, j, mpi_real, pe_root, letkf_mpi_comm, ierr)
       
    END DO


    ! print summary of grids
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "Vertical grids summary"
       PRINT *, "----------------------------------------"
       DO i=1,SIZE(vtgrids)
          PRINT *, "Name: ", vtgrids(i)%name
          PRINT *, " dims:", vtgrids(i)%dims
          IF (vtgrids(i)%dims == 1) THEN
             PRINT *, " size:", SIZE(vtgrids(i)%vert, dim=1)
          ELSE
             CALL letkf_mpi_abort("Not yet implemented")
          END IF
          PRINT *, " range: ", MINVAL(vtgrids(i)%vert), MAXVAL(vtgrids(i)%vert)
          PRINT *, ""
       END DO
    END IF



    ! Load the state variable specifications
    !---------------------------------------------------------------------------
    ! broadcast the state variable definitions
    IF (pe_isroot) i=SIZE(statevars)
    CALL mpi_bcast(i, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
    IF (.NOT. pe_isroot)  ALLOCATE(statevars(i))
    DO i=1,SIZE(statevars)
       CALL mpi_bcast(statevars(i)%name, 20, mpi_character,&
            pe_root, letkf_mpi_comm, ierr)
       CALL mpi_bcast(statevars(i)%hzgrid, 20, mpi_character,&
            pe_root, letkf_mpi_comm, ierr)
       CALL mpi_bcast(statevars(i)%vtgrid, 20, mpi_character,&
            pe_root, letkf_mpi_comm, ierr)
       CALL mpi_bcast(statevars(i)%levels, 1, mpi_integer, &
            pe_root, letkf_mpi_comm, ierr)
       CALL mpi_bcast(statevars(i)%grid_s_idx, 1, mpi_integer, &
            pe_root, letkf_mpi_comm, ierr)
    END DO

    ! print out stats
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "State variables summary"
       PRINT *, "----------------------------------------"
       DO i=1,SIZE(statevars)
          PRINT *, 'Name: "', TRIM(statevars(i)%name),'"'
          PRINT *, ' hz grid: "', TRIM(statevars(i)%hzgrid),'"'
          PRINT *, ' vt grid: "', TRIM(statevars(i)%vtgrid),'"'
          PRINT '(X,A,I0,A,I0)', ' slab:    ',statevars(i)%grid_s_idx,' to ',&
               statevars(i)%grid_s_idx+statevars(i)%levels-1
          PRINT *, ""
       END DO
    END IF

    CALL timing_stop("read_state_specs")

  END SUBROUTINE letkf_state_init_spec




  !--------------------------------------------------------------------------------
  ! > Writes the state ensemble in parallel
  SUBROUTINE letkf_state_write_ens()
    integer :: i, j, k, p, s, tag, sends_num, recvs_num, ierr
    integer, allocatable :: recvs(:), sends(:)
    real, allocatable :: tmp_r_3d(:,:,:)
    integer :: pe_ens_io(ens_size)
    
    if(pe_isroot) print '(/,X,A)', "Writing analysis ensemble members..."
    call timing_start("write_ens", TIMER_SYNC)

    !initialize output files
    do i=1,ens_size
       if (letkf_mpi_nextio() == pe_rank) &
            call stateio_class%write_init('ana', i)
    end do
    call letkf_mpi_barrier()

    ! determine who does the IO
    do i=1,ens_size
       pe_ens_io(i) = letkf_mpi_nextio()
    end do

    ! initialize all the sends
    allocate(sends(grid_ns*ens_size))
    sends_num=0
    do j=1,ens_size
       s=0
       do i=1,size(statevars)
          do k=1,statevars(i)%levels
             s = s + 1
             tag = k-1+statevars(i)%grid_s_idx + (j-1)*grid_ns
             sends_num = sends_num + 1
             call mpi_isend(state_ij(j,s,1), ij_count, mpitype_grid_nk_ns, &
                  pe_ens_io(j), tag, letkf_mpi_comm, sends(sends_num), ierr)
          end do
       end do
    end do

    ! initialize all the receives
    allocate(recvs(grid_ns*pe_size))
    do j=1,ens_size
       do i=1,size(statevars)
          call timing_start("io_write")
          call timing_stop("io_write")
          call timing_start("mpi_gather")
          call timing_stop("mpi_gather")

          if (pe_ens_io(i) == pe_rank) then
             call timing_start("mpi_gather")
             if (allocated(tmp_r_3d)) deallocate(tmp_r_3d)
             allocate(tmp_r_3d(grid_nx, grid_ny, statevars(i)%levels))
             recvs_num =0
             
             ! initialize the receives for this variable
             do p=0, pe_size-1
                do k=1, statevars(i)%levels
                   recvs_num = recvs_num +1
                   tag=k-1+statevars(i)%grid_s_idx + (j-1)*grid_ns
                   call mpi_irecv(tmp_r_3d(p+1,1,k), ij_count_pe(p), &
                        mpitype_grid_nxy_real, p, tag, letkf_mpi_comm, &
                        recvs(recvs_num), ierr)
                end do
             end do

             ! wait for the receives to finish
             call mpi_waitall(recvs_num, recvs, MPI_STATUSES_IGNORE, ierr)
             call timing_stop("mpi_gather")
             

             ! write out to file
             call timing_start("io_write")
             ! TODO add verbose output
             call stateio_class%write_state('ana', j, statevars(i)%name, tmp_r_3d)
             call timing_stop("io_write")
             
          end if
       end do
    end do

    ! wait for the sends to finish
    call mpi_waitall(sends_num, sends, MPI_STATUSES_IGNORE, ierr)
    call timing_stop("write_ens")

    
!    IF (pe_isroot) call letkf_mpi_abort("WARNING: state ens write not yet implemented!!")
  END SUBROUTINE letkf_state_write_ens



  !--------------------------------------------------------------------------------
  !> Reads the state ensemble in parallel.
  !! Responsibility for reading in each variable of each ensemble member is spread
  !! across the PEs. After each variable is read, it is immediately scattered to
  !! all participating PEs
  SUBROUTINE letkf_state_read_ens()
    INTEGER, ALLOCATABLE :: requests(:)
    INTEGER, ALLOCATABLE :: sends(:)
    REAL,    ALLOCATABLE :: tmp_r_3d(:,:,:)
    INTEGER :: sends_cnt
    LOGICAL :: do_send
    INTEGER :: i, j, k, p, s, ierr, tag

    CALL timing_start("read_state")

    ! Load and scatter the ensemble info
    ! ---------------------------------------------
    ! TODO more efficient to do entire variable in one call?
    ALLOCATE(requests(grid_ns*ens_size))
    ALLOCATE(sends(grid_ns*pe_size))

    IF (pe_isroot) THEN
       PRINT *, "Reading and distributing ensemble state..."
    END IF

    ! initialize all the receives
    !------------------------------
    DO j=1,ens_size  ! for each ensemble member
       s=0
       DO i=1,SIZE(statevars)  ! for each state variable
          DO k=1,statevars(i)%levels
             s=s+1
             tag=s + (j-1)*grid_ns
             CALL mpi_Irecv(state_ij(j,s,1), ij_count, mpitype_grid_nk_ns, &
                  MPI_ANY_SOURCE, tag, letkf_mpi_comm, requests(tag), ierr)
          END DO
       END DO
    END DO

    ! initialize all the sends
    DO j=1,ens_size  ! for each ensemble member
       s=0
       DO i=1,SIZE(statevars)  ! for each state variable

          !determine who should do the file loading for this state var
          do_send = letkf_mpi_nextio() == pe_rank

          ! read the file
          CALL timing_start("io_read")
          IF(do_send) CALL stateio_class%read_state(j, statevars(i)%name, tmp_r_3d)
          CALL timing_stop("io_read")

          ! scatter the variable
          CALL timing_start("mpi_send")
          sends_cnt=0
          DO k=1,statevars(i)%levels
             s=s+1
             tag=s + (j-1)*grid_ns
             IF(do_send) THEN
                DO p=0,pe_size-1
                   sends_cnt = sends_cnt+1
                   CALL mpi_Isend(tmp_r_3d(p+1,1,k), ij_count_pe(p),&
                        mpitype_grid_nxy_real, p, tag, letkf_mpi_comm,&
                        sends(sends_cnt), ierr)
                END DO
             END IF
          END DO

          ! wait for the send to finish
          IF(do_send) THEN
             CALL mpi_waitall(sends_cnt, sends, MPI_STATUSES_IGNORE, ierr)
             DEALLOCATE(tmp_r_3d)
          END IF
          CALL timing_stop("mpi_send")
       END DO
    END DO

    ! wait for all the receives to finish
    CALL mpi_waitall(grid_ns*ens_size, requests, MPI_STATUSES_IGNORE, ierr)

    CALL timing_stop("read_state")

  END SUBROUTINE letkf_state_read_ens



  !--------------------------------------------------------------------------------
  !> register user-defined and built-in state I/O classes to be available
  !! for use by LETKF. Actual state I/O class to be used is specified in the
  !! namelist
  SUBROUTINE letkf_state_register(ioclass)
    CLASS(letkf_stateio), POINTER :: ioclass
    INTEGER :: i

    ! make sure we han't reached our max number of classes
    IF ( pe_isroot ) THEN
       IF (stateio_reg_num == stateio_reg_max) THEN
          CALL letkf_mpi_abort("too many stateio classes have been registered.")
       END IF
    END IF

    ! make sure a class of this name hasn't already been registered
    IF ( pe_isroot ) THEN
       DO i=1, stateio_reg_num
          IF (toupper(stateio_reg(i)%p%name()) == toupper(ioclass%name())) THEN
             CALL letkf_mpi_abort("can't register stateio class '"// &
                  toupper(ioclass%name())// &
                  "', a class by that name already has been registered.")
          END IF
       END DO
    END IF

    ! add in the class to the list
    stateio_reg_num = stateio_reg_num + 1
    stateio_reg(stateio_reg_num)%p => ioclass

  END SUBROUTINE letkf_state_register


  
  !--------------------------------------------------------------------------------
  function letkf_state_var_getdef(name) result(res)
    character(len=*), intent(in) :: name
    type(letkf_statevar_spec) :: res

    character(:), allocatable :: name0
    integer :: i

    ! TODO, convert to upper, but make sure all the other instances
    ! are already in uppercase
    !    name0 = trim(toupper(name))
    name0 = trim(name)
    do i=1,size(statevars)
       if(statevars(i)%name == name0) exit
    end do
    if (i > size(statevars)) &
         call letkf_mpi_abort("state definition for variable '"//name0//"' not found")
    res = statevars(i)
  end function letkf_state_var_getdef

  
  !--------------------------------------------------------------------------------

  !< Convert a string to uppercase
  FUNCTION toupper(in_str) RESULT(out_str)
    CHARACTER(*), INTENT(in) :: in_str
    CHARACTER(LEN(in_str)) :: out_str
    INTEGER :: i
    INTEGER, PARAMETER :: offset = 32

    out_str = in_str
    DO i = 1, LEN(out_str)
       IF (out_str(i:i) >= "a" .AND. out_str(i:i) <= "z") THEN
          out_str(i:i) = ACHAR(IACHAR(out_str(i:i)) - offset)
       END IF
    END DO

  END FUNCTION toupper


END MODULE letkf_state
