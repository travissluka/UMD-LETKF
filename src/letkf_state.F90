!================================================================================
!> Module providing model state IO and access to the PE-scattered state.
!!
!--------------------------------------------------------------------------------
MODULE letkf_state
  USE timing
  USE mpi
  USE letkf_config
  USE letkf_mpi

  IMPLICIT NONE
  PRIVATE




  !================================================================================
  !================================================================================
  ! Public module components
  !================================================================================
  !================================================================================


  PUBLIC :: letkf_state_register
  PUBLIC :: letkf_state_init
  PUBLIC :: letkf_state_write_ens
  PUBLIC :: letkf_state_write_meansprd
  PUBLIC :: letkf_state_var_getdef
  PUBLIC :: letkf_state_hzgrid_getdef
  PUBLIC :: letkf_state_vtgrid_getdef



  !================================================================================
  !> Horizontal grid specification
  !!
  !! \todo remove the nominal lat/lon and keep only in child IO class?
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC :: letkf_hzgrid_spec
     CHARACTER(len=20)         :: name      !< unique name of the horizontal grid
     REAL,         ALLOCATABLE :: lat(:,:)  !< 2D lattidue grid (in degrees)
     REAL,         ALLOCATABLE :: lon(:,:)  !< 2D longitude grid (in degrees)
     REAL,         ALLOCATABLE :: lat_nom(:)!< nominal 1D lattitude grid (in degrees)
     REAL,         ALLOCATABLE :: lon_nom(:)!< nominal 1D longitude grid (in degrees)
     LOGICAL,      ALLOCATABLE :: mask(:,:) !< 2D mask
  END TYPE letkf_hzgrid_spec
  !================================================================================



  !================================================================================
  !> Veritical grid specification
  !!
  !! \todo remove nominal vertical grid and keep only in child IO class?
  !! \todo need to scatter the vertical grids the same as the horizontal
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC :: letkf_vtgrid_spec
     CHARACTER(len=20)         :: name !< unique name of the vertical grid
     INTEGER                   :: dims !< dimensions of the grid (0, 1, 2, or 3)
     REAL,         ALLOCATABLE :: vert(:,:,:) !< vertical grid
     REAL,         ALLOCATABLE :: vert_nom(:) !< nominal 1D vertival grid
  END TYPE letkf_vtgrid_spec
  !--------------------------------------------------------------------------------



  !================================================================================
  !> state variable specifications
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC :: letkf_statevar_spec
     !> unique name of the state variable
     CHARACTER(len=20) :: name
     !> name of horizontal grid (letkf_obs::letkf_hzgrid_spec) this variable is on
     CHARACTER(len=20) :: hzgrid
     !> name of vertial grid (letkf_obs::letkf_vtgrid_spec) this variable is on
     CHARACTER(len=20) :: vtgrid
     !> the internal starting index of the "slab" this variable is on
     INTEGER           :: grid_s_idx
     !> number of vertical levels this state variable consists of
     INTEGER           :: levels
  END TYPE letkf_statevar_spec
  !--------------------------------------------------------------------------------



  !================================================================================
  !> Abstract base class for model state reading and writing.
  !! All user-defined and built-in state file I/O classes for
  !! specific file types should inherit this class
  !--------------------------------------------------------------------------------
  TYPE, ABSTRACT, PUBLIC :: letkf_stateio
     LOGICAL :: verbose = .FALSE.
   CONTAINS
     PROCEDURE(I_letkf_stateio_getstr), NOPASS, DEFERRED :: name
     PROCEDURE(I_letkf_stateio_getstr), NOPASS, DEFERRED :: desc
     PROCEDURE(I_letkf_stateio_init),           DEFERRED :: init
     PROCEDURE(I_letkf_stateio_read_specs),     DEFERRED :: read_specs
     PROCEDURE(I_letkf_stateio_read_state),     DEFERRED :: read_state
     PROCEDURE(I_letkf_stateio_write_state),    DEFERRED :: write_state
     PROCEDURE(I_letkf_stateio_write_init),     DEFERRED :: write_init
  END TYPE letkf_stateio

  ABSTRACT INTERFACE
     FUNCTION I_letkf_stateio_getstr()
       CHARACTER(:), ALLOCATABLE :: I_letkf_stateio_getstr
     END FUNCTION I_letkf_stateio_getstr

     SUBROUTINE I_letkf_stateio_init(self, config)
       IMPORT letkf_stateio
       IMPORT configuration
       CLASS(letkf_stateio)             :: self
       TYPE(configuration), INTENT(in)  :: config
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
  !================================================================================


  ! enumerations for the stateio_write_* methods
  !------------------------------------------------------------
  INTEGER, PUBLIC, PARAMETER :: ENS_BKG_MEAN = -1
  INTEGER, PUBLIC, PARAMETER :: ENS_ANA_MEAN = -2
  INTEGER, PUBLIC, PARAMETER :: ENS_BKG_SPRD = -3
  INTEGER, PUBLIC, PARAMETER :: ENS_ANA_SPRD = -4


  ! grid definition, after being distributed across PEs
  !------------------------------------------------------------

  !> Latitude of gridpoints (in degrees) that have been scattered to this PE
  REAL,    PUBLIC, PROTECTED, ALLOCATABLE :: lat_ij(:)

  !> Longitude of gridpoints (in degrees) that have been scattered to this PE
  REAL,    PUBLIC, PROTECTED, ALLOCATABLE :: lon_ij(:)

  !> mask values of gridpoints that have been scattered to this PE
  LOGICAL, PUBLIC, PROTECTED, ALLOCATABLE :: mask_ij(:)


  ! state variables after being distributed across PEs
  !------------------------------------------------------------

  !> ensemble state variables that have been scattered to this PE.
  !! initially this is the background state perturbations, but stores
  !! the final analysis state after the LETKF solver is run
  REAL, PUBLIC, ALLOCATABLE :: state_ij(:,:,:)

  !> state variables mean that have been scattered to this PE.
  !! initially this is the background state mean, but stores the analysis
  !! state mean after the LETKF solver is run
  REAL, PUBLIC, ALLOCATABLE :: state_mean_ij(:,:)

  !> state variables spread that have been scattered to this PE.
  !! Initially this is the background state spread, but stores the analysis
  !! state spread after the LETKF solver is run
  REAL, PUBLIC, ALLOCATABLE :: state_sprd_ij(:,:)

  !> number of gridpoints in the X direction.
  !! \todo generalize this for handling multiple horizontal grids
  INTEGER, PUBLIC, PROTECTED :: grid_nx

  !> number of gridpoints in the Y direction.
  !! \todo generalize this for handling multiple horizontal grids
  INTEGER, PUBLIC, PROTECTED :: grid_ny

  !> number of state slabs (i.e. state variables times vertical levels
  INTEGER, PUBLIC, PROTECTED :: grid_ns

  ! master grid and variable definitions
  TYPE(letkf_hzgrid_spec),  PUBLIC, PROTECTED, ALLOCATABLE :: hzgrids(:)
  TYPE(letkf_vtgrid_spec),  PUBLIC, PROTECTED, ALLOCATABLE :: vtgrids(:)
  TYPE(letkf_statevar_spec),PUBLIC, PROTECTED, ALLOCATABLE :: statevars(:)




  !================================================================================
  !================================================================================
  ! Private module components
  !================================================================================
  !================================================================================



  !================================================================================
  !> simple wrapper for letkf_stateio so that we can have an array of pointers.
  !! to abstract classes
  !--------------------------------------------------------------------------------
  TYPE stateio_ptr
     CLASS(letkf_stateio), POINTER :: p
  END TYPE stateio_ptr
  !================================================================================



  ! registration of built-in and user-defined stateio classes
  !------------------------------------------------------------

  INTEGER, PARAMETER :: stateio_reg_max = 100         !< max classes to handle
  INTEGER            :: stateio_reg_num = 0           !< current number of classes registered
  TYPE(stateio_ptr)  :: stateio_reg(stateio_reg_max)  !< list of registered classes
  CLASS(letkf_stateio), POINTER :: stateio_class      !< The actual I/O module selected




  !================================================================================
  !================================================================================
CONTAINS
  !================================================================================
  !================================================================================



  !================================================================================
  !> \cond INTERNAL
  !> Initialize the letkf state module.
  !!
  !! The namefile is read, stateio class is initialized (reading in the grid specs),
  !! ensemble background state read in in parallel, mean and spread are calculated
  !! and saved to a file.
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_state_init(config)
    TYPE(configuration), INTENT(in) :: config

    TYPE(configuration) :: config_ioclass
    INTEGER ::  i,  s, n
    CHARACTER(:), ALLOCATABLE :: ioclass
    LOGICAL :: write_bkg_meansprd = .TRUE.
    LOGICAL :: verbose = .FALSE.


    ! print header
    IF (pe_isroot) THEN
       PRINT '(//A)', ""
       PRINT *, "============================================================"
       PRINT *, " letkf_state_init() : Model state I/O module initialization"
       PRINT *, "============================================================"
       PRINT *, ""
    END IF

    ! read some global configuration settings
    CALL config%get("verbose", verbose, default=.FALSE.)
    IF (pe_isroot) PRINT *, "state.verbose=",verbose

    ! print a list of all stateio classes that have been registered
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "List of stateio classes registered:"
       DO i=1, stateio_reg_num
          PRINT *, ' * "', tolower(stateio_reg(i)%p%name()), &
               '"  (', stateio_reg(i)%p%desc(), ")"
       END DO
       PRINT *, ""
    END IF

    ! determine the stateio class to use
    CALL config%get("ioclass", ioclass)
    ioclass = tolower(ioclass)
    IF (pe_isroot) PRINT '(A,A)',  " state.ioclass=",ioclass
    NULLIFY(stateio_class)
    DO i=1, stateio_reg_num
       IF (tolower(stateio_reg(i)%p%name()) == ioclass) THEN
          stateio_class => stateio_reg(i)%p
          EXIT
       END IF
    END DO
    IF (.NOT. ASSOCIATED(stateio_class)) THEN
       CALL letkf_mpi_abort('stateio class "'//ioclass//'" not found.')
    END IF

    ! initialize the stateio class
    stateio_class%verbose = verbose
    config_ioclass = config ! TODO, get the correct tree
    CALL stateio_class%init(config_ioclass)

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
  !> \endcond
  !================================================================================



  !================================================================================
  !> \cond INTERNAL
  !> Writes the ensemble mean and spread out to a set of files.
  !!
  !! This subroutine handles both the analysis and background cases, depending
  !! on the value of the argument. The current contents of "state_mean_ij" and
  !! "state_sprd_ij" are used for the output.
  !--------------------------------------------------------------------------------
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
             IF(ALLOCATED(tmp_r_3d)) DEALLOCATE(tmp_r_3d)
             ALLOCATE(tmp_r_3d(grid_nx,grid_ny,statevars(i)%levels))

             ! initialize the receives
             CALL timing_start("mpi_recv")
             recvs_num=0
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
             CALL stateio_class%write_state(mode, &
                  MERGE(ENS_BKG_MEAN, ENS_BKG_SPRD,j==1),&
                  TRIM(statevars(i)%name), tmp_r_3d(:,:,:))
             CALL timing_stop("io_write")

          END IF
       END DO
    END DO

    ! wait for the sends to finish
    CALL mpi_waitall(sends_num, sends, MPI_STATUSES_IGNORE, ierr)
    CALL timing_stop("write_meansprd")

  END SUBROUTINE letkf_state_write_meansprd
  !> \endcond
  !================================================================================



  !================================================================================
  !> \cond INTERNAL
  !> Reads in the grid specification (horizontal and vertical grids) and the
  !! specification of what state variables are involved with the state IO
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_state_init_spec()
    INTEGER :: i, j, ierr

    CALL timing_start("read_state_specs")


    !------------------------------------------------------------------------
    ! Load the grid/state specification
    !------------------------------------------------------------------------
    IF (pe_isroot) THEN

       ! TODO currently only handling the first hzgrid, and first vtgrid,
       ! modify to handle an arbitrary number
       !TODO, no point keeping the separate lat,lon,mask variables around
       CALL stateio_class%read_specs(hzgrids, vtgrids, statevars)


       ! Check the horizontal grid, print summary of grid stats
       !------------------------------------------------------------------------
       IF (.NOT. ALLOCATED(hzgrids) )&
            CALL letkf_mpi_abort('"HZGRIDS" was not allocated by state_io class.')

       DO i=1,SIZE(hzgrids)
          IF (.NOT. ALLOCATED(hzgrids(i)%lat) ) &
               CALL letkf_mpi_abort("LAT uninitiallized from stateio_class%read_specs()")
          IF (.NOT. ALLOCATED(hzgrids(i)%lon) ) &
               CALL letkf_mpi_abort("LON uninitiallized from stateio_class%read_specs()")
          IF (.NOT. ALLOCATED(hzgrids(i)%mask) ) &
               CALL letkf_mpi_abort("MASK uninitiallized from stateio_class%read_specs()")
          DO j=1,2
             IF ((SIZE(hzgrids(i)%lat, dim=j)/=SIZE(hzgrids(i)%lon,dim=j)) .OR. &
                  (SIZE(hzgrids(i)%lat,dim=j)/=SIZE(hzgrids(i)%mask,dim=j)))  THEN
                CALL letkf_mpi_abort("LAT,LON,MASK must have same shapes "//&
                     "from stateio_class%read_specs()")
             END IF
          END DO
          IF (.NOT. ALLOCATED(hzgrids(i)%lat_nom)) &
               CALL letkf_mpi_abort("LAT_NOM uninitialized from stateio_class%read_specs()")
          IF (.NOT. ALLOCATED(hzgrids(i)%lon_nom)) &
               CALL letkf_mpi_abort("LON_NOM uninitialized from stateio_class%read_specs()")
          !TODO, create nominal lat/lon if not already given
       END DO

       PRINT *, ""
       PRINT *, "Horizontal grids summary"
       PRINT *, "----------------------------------------"
       DO i=1,SIZE(hzgrids)
          PRINT *, "Name: ", hzgrids(i)%name
          PRINT "(X,A,I12)", " grid_nx:    ", SIZE(hzgrids(i)%lat, dim=1)
          PRINT "(X,A,I12)", " grid_ny:    ", SIZE(hzgrids(i)%lat, dim=2)
          PRINT "(X,A,I12)", " grid total: ", SIZE(hzgrids(i)%lat)
          PRINT "(X,A,I12,A,F4.1,A)", " grid masked:", COUNT(hzgrids(i)%mask), &
               " (",REAL(COUNT(hzgrids(i)%mask)*100) / REAL(SIZE(hzgrids(i)%lat)),"%)"
          PRINT *, " lat range: ", MINVAL(hzgrids(i)%lat), MAXVAL(hzgrids(i)%lat)
          PRINT *, " lon range: ", MINVAL(hzgrids(i)%lon), MAXVAL(hzgrids(i)%lon)
          PRINT *, ""
       END DO


       ! check vertical grids, print summary of grid stats
       !------------------------------------------------------------------------
       IF (.NOT. ALLOCATED(vtgrids)) &
            CALL letkf_mpi_abort('"vtgrids" not allocated by stateio_class')

       DO i=1,SIZE(vtgrids)
          ! TODO check name
          IF(.NOT. ALLOCATED(vtgrids(i)%vert)) &
               CALL letkf_mpi_abort('VERT not defined for vtgrid "'//&
               TRIM(vtgrids(i)%name)//'"')
          IF(.NOT. ALLOCATED(vtgrids(i)%vert_nom)) &
               CALL letkf_mpi_abort('VERT_NOM not defined for vtgrid "'//&
               TRIM(vtgrids(i)%name)//'"')
       END DO

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


       ! check state vars
       !------------------------------------------------------------------------
       IF (.NOT. ALLOCATED(statevars)) &
            CALL letkf_mpi_abort('"statevars" not allocated by stateio_class')


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
             CALL letkf_mpi_abort('vertical grid "'//TRIM(statevars(i)%vtgrid)//&
                  '" not found. Needed for state variable "'//&
                  TRIM(statevars(i)%name)//'".')
             STOP 1
          END IF
          grid_ns = grid_ns + statevars(i)%levels
       END DO

       PRINT *, ""
       PRINT *, "State variables summary"
       PRINT *, "----------------------------------------"
       DO i=1,SIZE(statevars)
          PRINT *, 'Name: "', TRIM(statevars(i)%name),'"'
          PRINT *, ' hz grid: "', TRIM(statevars(i)%hzgrid),'"'
          PRINT *, ' vt grid: "', TRIM(statevars(i)%vtgrid),'"'
          PRINT '(X,A,I0,A,I0)', ' slabs:    ',statevars(i)%grid_s_idx,' to ',&
               statevars(i)%grid_s_idx+statevars(i)%levels-1
          PRINT *, ""
       END DO


       ! get the overall horizontal grid dimensions
       ! TODO use the actually specified hz and vt grids
       grid_nx = SIZE(hzgrids(1)%lat,dim=1)
       grid_ny = SIZE(hzgrids(1)%lat,dim=2)
    END IF



    !---------------------------------------------------------------------------
    ! Scatter the grid data and variable definitions
    !---------------------------------------------------------------------------

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
       IF(pe_isroot) j=SIZE(hzgrids(i)%lat_nom)
       CALL mpi_bcast(j, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
       IF(.NOT. pe_isroot) ALLOCATE(hzgrids(i)%lat_nom(j))
       CALL mpi_bcast(hzgrids(i)%lat_nom, j, mpi_real, pe_root, letkf_mpi_comm, ierr)

       ! send nominal longitude
       IF(pe_isroot) j=SIZE(hzgrids(i)%lon_nom)
       CALL mpi_bcast(j, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
       IF(.NOT. pe_isroot) ALLOCATE(hzgrids(i)%lon_nom(j))
       CALL mpi_bcast(hzgrids(i)%lon_nom, j, mpi_real, pe_root, letkf_mpi_comm, ierr)

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
       IF(pe_isroot) j=SIZE(vtgrids(i)%vert_nom)
       CALL mpi_bcast(j, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
       IF(.NOT. pe_isroot) ALLOCATE(vtgrids(i)%vert_nom(j))
       CALL mpi_bcast(vtgrids(i)%vert_nom, j, mpi_real, pe_root, letkf_mpi_comm, ierr)

    END DO


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


    CALL timing_stop("read_state_specs")

  END SUBROUTINE letkf_state_init_spec
  !> \endcond
  !================================================================================



  !================================================================================
  !> \cond INTERNAL
  !> Writes the state ensemble in parallel
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_state_write_ens()
    INTEGER :: i, j, k, p, s, tag, sends_num, recvs_num, ierr
    INTEGER, ALLOCATABLE :: recvs(:), sends(:)
    REAL, ALLOCATABLE :: tmp_r_3d(:,:,:)
    INTEGER :: pe_ens_io(ens_size)

    IF(pe_isroot) PRINT '(/,X,A)', "Writing analysis ensemble members..."
    CALL timing_start("write_ens", TIMER_SYNC)


    ! determine who does the IO
    DO i=1,ens_size
       pe_ens_io(i) = letkf_mpi_nextio()
    END DO

    !initialize output files
    DO i=1,ens_size
       IF (pe_ens_io(i) == pe_rank) &
            CALL stateio_class%write_init('ana', i)
    END DO


    ! initialize all the sends
    ALLOCATE(sends(grid_ns*ens_size))
    sends_num=0
    DO j=1,ens_size
       s=0
       DO i=1,SIZE(statevars)
          DO k=1,statevars(i)%levels
             s = s + 1
             tag = k-1+statevars(i)%grid_s_idx + (j-1)*grid_ns
             sends_num = sends_num + 1
             CALL mpi_isend(state_ij(j,s,1), ij_count, mpitype_grid_nk_ns, &
                  pe_ens_io(j), tag, letkf_mpi_comm, sends(sends_num), ierr)
          END DO
       END DO
    END DO

    ! initialize all the receives
    ALLOCATE(recvs(grid_ns*pe_size))
    DO j=1,ens_size
       DO i=1,SIZE(statevars)
          CALL timing_start("io_write")
          CALL timing_stop("io_write")
          CALL timing_start("mpi_gather")
          CALL timing_stop("mpi_gather")

          IF (pe_ens_io(j) == pe_rank) THEN
             CALL timing_start("mpi_gather")
             IF (ALLOCATED(tmp_r_3d)) DEALLOCATE(tmp_r_3d)
             ALLOCATE(tmp_r_3d(grid_nx, grid_ny, statevars(i)%levels))
             recvs_num =0

             ! initialize the receives for this variable
             DO p=0, pe_size-1
                DO k=1, statevars(i)%levels
                   recvs_num = recvs_num +1
                   tag=k-1+statevars(i)%grid_s_idx + (j-1)*grid_ns
                   CALL mpi_irecv(tmp_r_3d(p+1,1,k), ij_count_pe(p), &
                        mpitype_grid_nxy_real, p, tag, letkf_mpi_comm, &
                        recvs(recvs_num), ierr)
                END DO
             END DO

             ! wait for the receives to finish
             CALL mpi_waitall(recvs_num, recvs, MPI_STATUSES_IGNORE, ierr)
             CALL timing_stop("mpi_gather")


             ! write out to file
             CALL timing_start("io_write")
             ! TODO add verbose output
             CALL stateio_class%write_state(&
                  'ana', j, TRIM(statevars(i)%name), tmp_r_3d)
             CALL timing_stop("io_write")

          END IF
       END DO
    END DO

    ! wait for the sends to finish
    CALL mpi_waitall(sends_num, sends, MPI_STATUSES_IGNORE, ierr)
    CALL timing_stop("write_ens")

  END SUBROUTINE letkf_state_write_ens
  !> \endcond
  !================================================================================



  !================================================================================
  !> \cond INTERNAL
  !> Reads the state ensemble in parallel.
  !!
  !! Responsibility for reading in each variable of each ensemble member is spread
  !! across the PEs. After each variable is read, it is immediately scattered to
  !! all participating PEs
  !--------------------------------------------------------------------------------
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
          IF(do_send) THEN
             CALL stateio_class%read_state(j, statevars(i)%name, tmp_r_3d)
             IF(.NOT. ALLOCATED(tmp_r_3d)) &
                  CALL letkf_mpi_abort("no data returned from stateio_class%read_state()")
          END IF
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
  !> \endcond
  !================================================================================



  !================================================================================
  !> register user-defined and built-in state I/O classes to be available
  !! for use by LETKF. Actual state I/O class to be used is specified in the
  !! namelist
  !--------------------------------------------------------------------------------
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
          IF (tolower(stateio_reg(i)%p%name()) == tolower(ioclass%name())) THEN
             CALL letkf_mpi_abort("can't register stateio class "//'"'// &
                  tolower(ioclass%name())// &
                  '", a class by that name already has been registered.')
          END IF
       END DO
    END IF

    ! add in the class to the list
    stateio_reg_num = stateio_reg_num + 1
    stateio_reg(stateio_reg_num)%p => ioclass

  END SUBROUTINE letkf_state_register
  !================================================================================



  !================================================================================
  !> Get horizontal grid definition
  !--------------------------------------------------------------------------------
  FUNCTION letkf_state_hzgrid_getdef(name) RESULT(res)
    CHARACTER(len=*), INTENT(in) :: name
    TYPE(letkf_hzgrid_spec) :: res

    CHARACTER(:), ALLOCATABLE :: name0
    INTEGER :: i

    ! TODO, convert to upper, but make sure all the other instances
    ! are already in uppercase
    !    name0 = trim(toupper(name))
    name0 = TRIM(name)
    DO i=1,SIZE(hzgrids)
       IF(hzgrids(i)%name == name0) EXIT
    END DO
    IF (i > SIZE(hzgrids)) &
         CALL letkf_mpi_abort("hzgrid definition for '"//name0//"' not found")
    res = hzgrids(i)
  END FUNCTION letkf_state_hzgrid_getdef
  !================================================================================



  !================================================================================
  !> get vertical grid definition
  !--------------------------------------------------------------------------------
  FUNCTION letkf_state_vtgrid_getdef(name) RESULT(res)
    CHARACTER(len=*), INTENT(in) :: name
    TYPE(letkf_vtgrid_spec) :: res

    CHARACTER(:), ALLOCATABLE :: name0
    INTEGER :: i

    ! TODO, convert to upper, but make sure all the other instances
    ! are already in uppercase
    !    name0 = trim(toupper(name))
    name0 = TRIM(name)
    DO i=1,SIZE(vtgrids)
       IF(vtgrids(i)%name == name0) EXIT
    END DO
    IF (i > SIZE(vtgrids)) &
         CALL letkf_mpi_abort("vtgrid definition for '"//name0//"' not found")
    res = vtgrids(i)
  END FUNCTION letkf_state_vtgrid_getdef
  !================================================================================



  !================================================================================
  !> Get state variable definition
  !--------------------------------------------------------------------------------
  FUNCTION letkf_state_var_getdef(name) RESULT(res)
    CHARACTER(len=*), INTENT(in) :: name
    TYPE(letkf_statevar_spec) :: res

    CHARACTER(:), ALLOCATABLE :: name0
    INTEGER :: i

    ! TODO, convert to upper, but make sure all the other instances
    ! are already in uppercase
    !    name0 = trim(toupper(name))
    name0 = TRIM(name)
    DO i=1,SIZE(statevars)
       IF(statevars(i)%name == name0) EXIT
    END DO
    IF (i > SIZE(statevars)) &
         CALL letkf_mpi_abort("state definition for variable '"//name0//"' not found")
    res = statevars(i)
  END FUNCTION letkf_state_var_getdef
  !================================================================================



  !================================================================================
  !> Convert a string to lower
  !--------------------------------------------------------------------------------
  FUNCTION tolower(in_str) RESULT(out_str)
    CHARACTER(*), INTENT(in) :: in_str
    CHARACTER(LEN(in_str)) :: out_str
    INTEGER :: i
    INTEGER, PARAMETER :: offset = 32

    out_str = in_str
    DO i = 1, LEN(out_str)
       IF (out_str(i:i) >= "A" .AND. out_str(i:i) <= "Z") THEN
          out_str(i:i) = ACHAR(IACHAR(out_str(i:i)) + offset)
       END IF
    END DO

  END FUNCTION tolower
  !================================================================================

END MODULE letkf_state
