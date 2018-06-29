!================================================================================
!> Module for handling the observation I/O.
!!
!! There are two classes available by default (OBSIO_NC, OBSIO_TEST).
!! If a user wants to implement their own type of observation IO, they can do so
!! by deriving from the letkf_obs::letkf_obsio class, implmenting the required
!! methods, and registering the class with LETKF via letkf_obs::letkf_obs_register
!--------------------------------------------------------------------------------
MODULE letkf_obs
  USE mpi
  USE timing
  USE kdtree
  use running_stats_mod
  use letkf_config
  USE letkf_mpi

  IMPLICIT NONE
  PRIVATE




  !================================================================================
  !================================================================================
  ! Public module components
  !================================================================================
  !================================================================================


  PUBLIC :: letkf_obs_register
  PUBLIC :: letkf_obs_init
  PUBLIC :: letkf_obs_read
  PUBLIC :: letkf_obs_get
  public :: letkf_obs_getdef



  !================================================================================
  !> Container for a single observation definition.
  !!
  !! The per-ensemble observation opererator values are stored instead in the
  !! letkf_obs::obs_hx variable
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC :: letkf_observation
     INTEGER :: obsid   !< observation type id, as specified in 'obsdef' config file
     INTEGER :: platid  !< platform type id, as specified in 'platdef' config file
     REAL    :: lat     !< latitude (in degrees)
     REAL    :: lon     !< longitude (in degrees)
     REAL    :: zdim    !< z coordinate, in whatever unit is appropriate for the domain.
                        !! If this is a 2d observation, this value is ignored
     REAL    :: time    !< time (in hours) from base analysis time
     REAL    :: val     !< observation value
     REAL    :: err     !< standard deviation of observation error
     INTEGER :: qc      !< quality control flag. obs with non-zero values will not be used
  END TYPE letkf_observation
  !================================================================================



  !================================================================================
  !> container for a single user-defined observation or platform type
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC :: letkf_obsplatdef
     CHARACTER(len=10)  :: name      !< unique short name of observation/plat type
     INTEGER            :: id        !< unique integer identifier, used in obs types
     CHARACTER(len=200) :: name_long !< longer name / description
  END TYPE letkf_obsplatdef
  !================================================================================



  !================================================================================
  !> Abstract base class for observation file reading and writing.
  !!
  !! All user-defined, and built-in, observation file I/O classes for
  !! specific file types should inherit this class
  !--------------------------------------------------------------------------------
  TYPE, ABSTRACT, PUBLIC :: letkf_obsio
   CONTAINS
     !> Get the unique name of the reader
     PROCEDURE(I_letkf_obsio_getstr), NOPASS, DEFERRED :: name

     !> Get a descriptive string for the reader
     PROCEDURE(I_letkf_obsio_getstr), NOPASS, DEFERRED :: desc

     !> initialize the class (usually reading its section of the namelist)
     PROCEDURE(I_letkf_obsio_init),           DEFERRED :: init

     !> read the observations
     PROCEDURE(I_letkf_obsio_read_obs),       DEFERRED :: read_obs

     !> read the observation operator values for a given ensemble member
     PROCEDURE(I_letkf_obsio_read_hx),        DEFERRED :: read_hx
  END TYPE letkf_obsio
  !================================================================================



  !================================================================================
  ! Interface methods for the letkf_obsio class
  !--------------------------------------------------------------------------------
  ABSTRACT INTERFACE
     FUNCTION I_letkf_obsio_getstr()
       CHARACTER(:), ALLOCATABLE :: I_letkf_obsio_getstr
     END FUNCTION I_letkf_obsio_getstr

     SUBROUTINE I_letkf_obsio_init(self, config)
       IMPORT letkf_obsio
       import configuration
       CLASS(letkf_obsio) :: self
       type(configuration), intent(in) :: config
     END SUBROUTINE I_letkf_obsio_init

     SUBROUTINE I_letkf_obsio_read_obs(self, obs)
       IMPORT letkf_obsio
       IMPORT letkf_observation
       CLASS(letkf_obsio) :: self
       TYPE(letkf_observation), ALLOCATABLE, INTENT(out) :: obs(:)
     END SUBROUTINE I_letkf_obsio_read_obs

     SUBROUTINE I_letkf_obsio_read_hx(self, ensmem, hx)
       IMPORT letkf_obsio
       CLASS(letkf_obsio) :: self
       INTEGER, INTENT(in) :: ensmem
       REAL, ALLOCATABLE, INTENT(out) :: hx(:)
     END SUBROUTINE I_letkf_obsio_read_hx
  END INTERFACE
  !================================================================================


  !> The list of valid observation definitions
  TYPE(letkf_observation), PUBLIC, PROTECTED, ALLOCATABLE :: obs_def(:)

  !> per-ensemble observation operator perturbations
  REAL, PUBLIC, PROTECTED, ALLOCATABLE :: obs_hx(:,:)

  !> The mean of the per-ensemble observation operators
  REAL, PUBLIC, PROTECTED, ALLOCATABLE :: obs_hx_mean(:)




  !================================================================================
  !================================================================================
  ! Private module components
  !================================================================================
  !================================================================================


  !================================================================================
  !> simple wrapper for letkf_obsio so that we can have an arary of pointers
  !--------------------------------------------------------------------------------
  TYPE obsio_ptr
     CLASS(letkf_obsio), POINTER :: p
  END TYPE obsio_ptr
  !================================================================================


  ! configuration paramters
  !------------------------------------------------------------

  !> list of all observation types
  TYPE(letkf_obsplatdef), ALLOCATABLE :: obsdef_list(:)

  !> list of all platform types
  TYPE(letkf_obsplatdef), ALLOCATABLE :: platdef_list(:)


  ! mpi derived types
  !------------------------------------------------------------
  !> MPI derived type for letkf_obs::observation
  INTEGER :: mpitype_observation
  INTEGER :: mpitype_real_nk


  ! registration of built-in and user-defined obsio classes
  !------------------------------------------------------------

  !> max number of classes to handle
  INTEGER, PARAMETER  :: obsio_reg_max = 100

  !> current number of classes registered
  INTEGER             :: obsio_reg_num = 0

  !> list of registered classes
  TYPE(obsio_ptr)     :: obsio_reg(obsio_reg_max)

  !> The actual I/O module selected
  CLASS(letkf_obsio), POINTER :: obsio_class

  !------------------------------------------------------------

  !> A kdtree of all the loaded observations, for quick searches of nearest obs
  TYPE(kd_root) :: obs_tree




  !================================================================================
  !================================================================================
CONTAINS
  !================================================================================
  !================================================================================



  !================================================================================
  !> \cond INTERNAL
  !> Initialize the letkf_obs module.
  !!
  !! Reads in our section of the namelist, and any observation definition and platorm
  !! definition configuration files.
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_obs_init(config)
    type(configuration), intent(in) :: config

    type(configuration) :: ioconfig, config_def, config_def0
    INTEGER :: i
    CHARACTER(:), ALLOCATABLE :: ioclass, str

    ! print header
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "============================================================"
       PRINT *, " letkf_obs_init() : observation module initialization"
       PRINT *, "============================================================"
    END IF

    ! print a list of all obsio classes that have been registered
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "List of obsio classes registered:"
       DO i=1, obsio_reg_num
          PRINT *, ' * "', tolower(obsio_reg(i)%p%name()), &
               '"  (', obsio_reg(i)%p%desc(), ")"
       END DO
       PRINT *, ""
    END IF

    
    ! determine which io class to use
    call config%get("ioclass", ioclass)
    ioclass = tolower(ioclass)
    NULLIFY(obsio_class)
    DO i=1,obsio_reg_num
       IF (tolower(obsio_reg(i)%p%name()) == ioclass) THEN
          obsio_class => obsio_reg(i)%p
          EXIT
       END IF
    END DO
    IF (.NOT. ASSOCIATED(obsio_class)) THEN
       CALL letkf_mpi_abort("obsio class "//ioclass// " not found.")
    END IF

    
    ! read in the observation definition configuration
    IF(pe_isroot) THEN
       print *, ""
       PRINT *, "observation definitions"
       PRINT *, "------------------------------------------------------------"
       print "(A8,A8,A5,A)", "NAME", "ID", "","DESCRIPTION"
    END IF
    call config%get("obsdef", config_def)
    allocate(obsdef_list(config_def%count()))
    do i=1,config_def%count()
       call config_def%get(i, config_def0)
       call config_def0%get(1, str)
       obsdef_list(i)%name = tolower(str)
       call config_def0%get(2, obsdef_list(i)%id)
       call config_def0%get(3, str)       
       obsdef_list(i)%name_long = str
       if(pe_isroot) &
            print "(A10, I6, A5, A)", trim(obsdef_list(i)%name), obsdef_list(i)%id, &
            "", trim(obsdef_list(i)%name_long)
    end do
   
    
    ! read in the platform definition configuration
    IF(pe_isroot) THEN
       print *, ""
       PRINT *, "platform definitions"
       PRINT *, "------------------------------------------------------------"
       print "(A8,A8,A5,A)", "NAME", "ID", "","DESCRIPTION"       
    END IF
    call config%get("platdef", config_def)
    allocate(platdef_list(config_def%count()))
    do i=1,config_def%count()
       call config_def%get(i, config_def0)
       call config_def0%get(1, str)
       platdef_list(i)%name = tolower(str)
       call config_def0%get(2, platdef_list(i)%id)
       call config_def0%get(3, str)       
       platdef_list(i)%name_long = str
       if(pe_isroot) &
            print "(A10, I6, A5, A)", trim(platdef_list(i)%name), platdef_list(i)%id, &
            "", trim(platdef_list(i)%name_long)
    end do
    if(pe_isroot) print *, ""
    
    ! initialize MPI object for later sending/receving obsservations
    if(pe_isroot) print *, ""    
    CALL init_mpi_observation()

    ! finish initialize  of the obsio class
    IF (pe_isroot) then
       print *, ""
       PRINT *, "Intializing I/O ioclass: ",obsio_class%name()
    end IF
    call config%get(ioclass, ioconfig)
    CALL obsio_class%init(ioconfig)

  END SUBROUTINE letkf_obs_init
  !> \endcond
  !================================================================================



  !================================================================================
  !> \cond INTERNAL
  !> Initialize the mpi derived type for "letkf_observation"
  !--------------------------------------------------------------------------------
  SUBROUTINE init_mpi_observation()
    INTEGER, PARAMETER :: n = 9
    INTEGER :: TYPE(n), blocklen(n)
    INTEGER(kind=mpi_address_kind) :: disp(n), base, lb, ex, ex_real
    INTEGER :: ierr, i
    TYPE(letkf_observation) :: ob

    ! get the addresses of each variable of the observation type
    i = 0
    i=i+1; CALL mpi_get_address(ob%obsid,  disp(i), ierr); TYPE(i) = mpi_integer
    i=i+1; CALL mpi_get_address(ob%platid, disp(i), ierr); TYPE(i) = mpi_integer
    i=i+1; CALL mpi_get_address(ob%lat,    disp(i), ierr); TYPE(i) = mpi_real
    i=i+1; CALL mpi_get_address(ob%lon,    disp(i), ierr); TYPE(i) = mpi_real
    i=i+1; CALL mpi_get_address(ob%zdim,   disp(i), ierr); TYPE(i) = mpi_real
    i=i+1; CALL mpi_get_address(ob%time,   disp(i), ierr); TYPE(i) = mpi_real
    i=i+1; CALL mpi_get_address(ob%val,    disp(i), ierr); TYPE(i) = mpi_real
    i=i+1; CALL mpi_get_address(ob%err,    disp(i), ierr); TYPE(i) = mpi_real
    i=i+1; CALL mpi_get_address(ob%qc,     disp(i), ierr); TYPE(i) = mpi_integer
    base = disp(1)
    disp(:) = disp(:) - base
    blocklen(:) = 1

    ! construct the observation mpi derived type
    CALL mpi_type_create_struct(n, blocklen, disp, TYPE, mpitype_observation, ierr)
    CALL mpi_type_commit(mpitype_observation, ierr)

    ! derived type for scattering obs_hx
    call mpi_type_get_extent(mpi_real, lb, ex_real, ierr)
    lb = 0
    ex = ex_real*ens_size
    call mpi_type_create_resized(mpi_real, lb, ex, mpitype_real_nk, ierr)
    call mpi_type_commit(mpitype_real_nk, ierr)

  END SUBROUTINE init_mpi_observation
  !> \endcond
  !================================================================================



  !================================================================================
  !> \cond INTERNAL
  !> Read the observations.
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_obs_read()
    INTEGER :: i, nobs, ierr
    REAL, ALLOCATABLE :: obs_lats(:), obs_lons(:), tmp_r(:)
    integer :: obshx_pe(ens_size), obs_pe

    CALL timing_start('read_obs')

    obs_pe = 0
    IF(pe_rank == obs_pe) THEN
       ! have the I/O class read the main observations file
       CALL obsio_class%read_obs(obs_def)

       ! check to make sure everything looks good
       IF (.NOT. ALLOCATED(obs_def)) &
            CALL letkf_mpi_abort("obsio class returned unallocated obs array")

       nobs=SIZE(obs_def)
    END IF

    ! give all PEs a copy of the obs
    CALL mpi_bcast(nobs, 1, mpi_integer, pe_root, letkf_mpi_comm, ierr)
    IF (.NOT. pe_isroot)  ALLOCATE(obs_def(nobs))
    CALL mpi_bcast(obs_def, nobs, mpitype_observation, pe_root, letkf_mpi_comm, ierr)

    ! determine who is reading in which obs_hx
    do i = 1, ens_size
       obshx_pe(i) = letkf_mpi_nextio()
    end do

    ! read in the observation operator
    !> \todo do this in parallel with the above obs read
    ALLOCATE(obs_hx(ens_size, nobs))    
    ALLOCATE(tmp_r(nobs))
    obs_hx = 0.0
    DO i=1,ens_size
       if(pe_rank == obshx_pe(i)) then
          CALL obsio_class%read_hx(i, tmp_r)
          obs_hx(i,:) = tmp_r
       end if
    END DO
    ! broadcast the obs_hx values
    do i=1,ens_size       
       call mpi_bcast(obs_hx(i,1), nobs, mpitype_real_nk, obshx_pe(i), letkf_mpi_comm, ierr)
    end do
    DEALLOCATE(tmp_r)


    !> \todo make sure nobs agrees for all files

    !> \todo make sure bad obs are removed, do extra QC checks

    ! calculate hx_mean
    ALLOCATE(obs_hx_mean(nobs))
    obs_hx_mean = SUM(obs_hx, 1)/ens_size

    ! remove mean from hx
    DO i=1,nobs
       obs_hx(:,i) = obs_hx(:,i) - obs_hx_mean(i)
    END DO


    if (nobs > 0) then
       ! print out observation statistics
       if (pe_isroot) CALL obs_print_stats(obs_def)

       !> \todo, make sure no bad qc obs go into the tree

       ! add obs to KD tree
       IF(pe_isroot) PRINT '(//,X,A)', "Adding observations to local KDTree..."
       ALLOCATE(obs_lons(nobs))
       ALLOCATE(obs_lats(nobs))
       DO i=1,nobs
          obs_lons(i) = obs_def(i)%lon
          obs_lats(i) = obs_def(i)%lat
       END DO
       CALL kd_init(obs_tree, obs_lons, obs_lats)
       DEALLOCATE(obs_lons)
       DEALLOCATE(obs_lats)
    else
       if (pe_isroot) print *, "WARNING: there are NO observations to assimilate"
    end if

    CALL timing_stop('read_obs')

  END SUBROUTINE letkf_obs_read
  !> \endcond
  !================================================================================



  !================================================================================
  !> Register user-defined observation I/O classes to be available for use by LETKF.
  !!
  !! Actual observation I/O class to be used is specified in the namelist.
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_obs_register(ioclass)

    !> An allocated instance of the observation I/O class to be registered
    CLASS(letkf_obsio), POINTER :: ioclass

    INTEGER :: i

    ! make sure we han't reached our max number of classes
    IF ( pe_isroot ) THEN
       IF (obsio_reg_num == obsio_reg_max) THEN
          CALL letkf_mpi_abort("too many obsio classes have been registered.")
       END IF
    END IF

    ! make sure a class of this name hasn't already been registered
    IF ( pe_isroot ) THEN
       DO i=1, obsio_reg_num
          IF (tolower(obsio_reg(i)%p%name()) == tolower(ioclass%name())) THEN
             CALL letkf_mpi_abort("can't register obsio class '"// &
                  tolower(ioclass%name())// &
                  "', a class by that name already has been registered.")
          END IF
       END DO
    END IF

    ! add in the class to the list
    obsio_reg_num = obsio_reg_num + 1
    obsio_reg(obsio_reg_num)%p => ioclass

  END SUBROUTINE letkf_obs_register
  !================================================================================



  !================================================================================
  !> For a given search lat/lon point and radius (slat,slon,sradius), return all
  !! the points that are within that radius.
  !!
  !! An array of indexes (robs) are returned which point to the appropriate
  !! observations in the module's obs arrays (letkf_obs::obs_def,
  !! letkf_obs::obs_hx, letkf_obs::obs_hx_mean)
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_obs_get(slat, slon, sradius, robs, rdist, rnum)
    REAL,    INTENT(in)  :: slat     !< Center of search latitude (degrees)
    REAL,    INTENT(in)  :: slon     !< Center of search longitudee (degrees)
    REAL,    INTENT(in)  :: sradius  !< search radius (meters)
    INTEGER, INTENT(out) :: robs(:)  !< indexes to the observations found
    REAL,    INTENT(out) :: rdist(:) !< distance of each observation (meters)
    INTEGER, INTENT(out) :: rnum     !< number of observations found

    rnum = 0
    if (size(obs_def) > 0) then
       CALL kd_search_radius(obs_tree, slon, slat, sradius, robs, rdist, rnum, .FALSE.)
    end if

  END SUBROUTINE letkf_obs_get
  !================================================================================



  !================================================================================
  !> print out statistics about the loaded observations
  !--------------------------------------------------------------------------------
  SUBROUTINE obs_print_stats(obs_t)

    !> The list of observations to print statistics about
    TYPE(letkf_observation) :: obs_t(:)

    INTEGER, ALLOCATABLE :: obst_count(:,:)
    INTEGER, ALLOCATABLE :: plat_count(:,:)

    type(running_stats), allocatable :: odep_stats(:), pdep_stats(:)
    type(running_stats), allocatable :: osprd_stats(:), psprd_stats(:)

    INTEGER :: cnt, i, j, cnt_total

    IF(.NOT. pe_isroot) RETURN

    cnt_total=0
    ALLOCATE(obst_count(SIZE(obsdef_list)+1, 4))
    ALLOCATE(plat_count(SIZE(platdef_list)+1, 4))
    allocate(odep_stats(SIZE(obsdef_list)))
    allocate(pdep_stats(SIZE(platdef_list)))
    allocate(osprd_stats(SIZE(obsdef_list)))
    allocate(psprd_stats(SIZE(platdef_list)))

    obst_count = 0
    plat_count = 0


    ! calculate the statistics
    DO i=1,SIZE(obs_t)
       ! for each ob, is it good or bad?
       IF(obs_t(i)%qc == 0) THEN
          cnt=2
          cnt_total = cnt_total + 1
       ELSE IF (obs_t(i)%qc > 0 ) THEN
          cnt=3
       ELSE
          cnt=4
       END IF

       ! count by obs type
       DO j=1, SIZE(obsdef_list)
          IF(obs_t(i)%obsid == obsdef_list(j)%id) THEN
             obst_count(j,1) = obst_count(j,1) + 1
             obst_count(j,cnt) = obst_count(j,cnt) +1

             ! departure stats if this is a good ob
             if (obs_t(i)%qc == 0) then
                call odep_stats(j)%add( obs_def(i)%val-obs_hx_mean(i))
                call osprd_stats(j)%add(sum(obs_hx(:,i)*obs_hx(:,i))/ens_size)
             end if

             EXIT
          END IF
       END DO
       IF(j > SIZE(obsdef_list)) obst_count(SIZE(obsdef_list)+1,1) = &
            obst_count(SIZE(obsdef_list)+1,1) + 1

       ! count by plat type
       DO j=1, SIZE(platdef_list)
          IF(obs_t(i)%platid == platdef_list(j)%id) THEN
             plat_count(j,1) = plat_count(j,1) + 1
             plat_count(j,cnt) = plat_count(j,cnt) +1

             ! departure stats if this is a good ob
             if (obs_t(i)%qc == 0) then
                call pdep_stats(j)%add( obs_def(i)%val-obs_hx_mean(i))
                call psprd_stats(j)%add(sum(obs_hx(:,i)*obs_hx(:,i))/ens_size)
             end if

             EXIT
          END IF
       END DO
       IF(j > SIZE(platdef_list)) plat_count(SIZE(platdef_list)+1,1) = &
            plat_count(SIZE(platdef_list)+1,1) + 1
    END DO


    ! Print observations counts
    PRINT *, ""
    PRINT *, ""
    PRINT *, "Observation Counts"
    PRINT *, "------------------------------------------------------------"
    PRINT '(I11,A)', SIZE(obs_t), " observations loaded"
    PRINT *, ""
    PRINT '(4A10,A12)', '','total','bad-obsop','bad-qc','good'
    PRINT *, '         ---------------------------------------------'
    DO i=1,SIZE(obst_count,1)
       IF (obst_count(i,1) == 0) CYCLE
       IF (i < SIZE(obst_count,1)) THEN
          PRINT '(A10,I10, 2I10, I10, A2, F5.1,A)', &
               TRIM(obsdef_list(i)%name), obst_count(i,1), &
               obst_count(i,3), obst_count(i,4), &
               obst_count(i,2), '(', REAL(obst_count(i,2))/obst_count(i,1)*100,&
               ')%'
       ELSE
          PRINT '(A10,I10,3A10,A)', 'unknown', obst_count(i,1), 'X','X','0', ' (  0.0)%'
       END IF
    END DO

    PRINT *, ""
    PRINT *, '         ---------------------------------------------'
    DO i=1,SIZE(plat_count,1)
       IF (plat_count(i,1) == 0) CYCLE
       IF (i < SIZE(plat_count,1)) THEN
          PRINT '(A10,I10, 2I10, I10, A2, F5.1,A)', &
               TRIM(platdef_list(i)%name), plat_count(i,1), &
               plat_count(i,3), plat_count(i,4), &
               plat_count(i,2), '(', REAL(plat_count(i,2))/plat_count(i,1)*100,&
               ')%'
       ELSE
          PRINT '(A10,I10,3A10,A)', 'unknown', plat_count(i,1), 'X','X','0', ' (  0.0)%'
       END IF
    END DO


    ! Print observation departure stats
    PRINT *, ""
    PRINT *, ""
    PRINT *, "Observation departures"
    PRINT *, "------------------------------------------------------------"
    PRINT *, ""
    PRINT '(A10,6A12)', '', 'rms', 'bias', 'min', 'max', 'ens-sprd'
    PRINT *, '         -------------------------------------------------------'
    DO i=1,SIZE(obst_count,1)
       IF (obst_count(i,1) == 0) CYCLE
       IF (i < SIZE(obst_count,1)) THEN
         PRINT '(A10,6F12.5)', &
              TRIM(obsdef_list(i)%name), sqrt(odep_stats(i)%mean(2)), &
              odep_stats(i)%mean(1), odep_stats(i)%min(), odep_stats(i)%max(), &
              sqrt(osprd_stats(i)%mean())
       END IF
    END DO


    PRINT *, ""
    PRINT *, '         -------------------------------------------------------'
    DO i=1,SIZE(plat_count,1)
       IF (plat_count(i,1) == 0) CYCLE
       IF (i < SIZE(plat_count,1)) THEN
         PRINT '(A10,5F12.5)', &
              TRIM(platdef_list(i)%name), sqrt(pdep_stats(i)%mean(2)), &
              pdep_stats(i)%mean(1), pdep_stats(i)%min(), pdep_stats(i)%max(), &
              sqrt(psprd_stats(i)%mean())
       END IF
    END DO

  END SUBROUTINE obs_print_stats
  !================================================================================



  !================================================================================
  !> Get and observation or platform definition information.
  !--------------------------------------------------------------------------------
  function letkf_obs_getdef(obs_plat, name) result(res)
    character(len=1), intent(in) :: obs_plat !< set to 'P' or 'O' if a platform or
                                             !! observation definiont is requested
    character(len=*), intent(in) :: name     !< platform or observation name to get
    type(letkf_obsplatdef) :: res            !< the returned definition

    character(:), allocatable :: name0
    integer :: i

    name0 = trim(tolower(name))
    i =1
    if (obs_plat == 'O') then
       do while (i <= size(obsdef_list))
          if(obsdef_list(i)%name == name0) exit
          i = i + 1
       end do
       if ( i > size(obsdef_list)) &
            call letkf_mpi_abort("observation definition for '"//name0//"' not found")
       res = obsdef_list(i)

    else if (obs_plat == 'P') then
       do while (i <= size(platdef_list))
          if(platdef_list(i)%name == name0) exit
          i = i + 1
       end do
       if ( i > size(platdef_list)) &
            call letkf_mpi_abort("platform definition for '"//name0//"' not found")
       res = platdef_list(i)

    else
       call letkf_mpi_abort("letkf_obs_getdef: first argument must be 'P' or 'O'")

    end if

  end function letkf_obs_getdef
  !================================================================================



  !================================================================================
  !> Convert a string to lowercase
  !--------------------------------------------------------------------------------
  FUNCTION tolower(in_str) RESULT(out_str)
    CHARACTER(*), INTENT(in) :: in_str !< input string
    CHARACTER(LEN(in_str)) :: out_str  !< output string

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


END MODULE letkf_obs
