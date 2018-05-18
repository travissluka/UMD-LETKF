MODULE letkf_obs
  USE mpi
  USE timing
  USE kdtree
  USE letkf_mpi
  use running_stats_mod
  
  IMPLICIT NONE
  PRIVATE


  !--------------------------------------------------------------------------------
  ! public subroutines
  !--------------------------------------------------------------------------------
  PUBLIC :: letkf_obs_register
  PUBLIC :: letkf_obs_init
  PUBLIC :: letkf_obs_read
  PUBLIC :: letkf_obs_get
  public :: letkf_obs_getdef


  !--------------------------------------------------------------------------------
  ! Public types
  !--------------------------------------------------------------------------------

  !> container for a single observation
  TYPE, PUBLIC :: letkf_observation
     INTEGER :: obsid   !< observation type id, as specified in 'obsdef' config file
     INTEGER :: platid  !< platform type id, as specified in 'platdef' config file
     REAL    :: lat     !< latitude (in degrees)
     REAL    :: lon     !< longitude (in degrees)
     REAL    :: zdim    !< z coordinate, in whatever unit is appropriate for the domain
     !! if this is a 2d observation, this value is ignored
     REAL    :: time    !< time (in hours) from base analysis time
     REAL    :: val     !< observation value
     REAL    :: err     !< standard deviation of observation error
     INTEGER :: qc      !< quality control flag. obs with non-zero values will not be used
  END TYPE letkf_observation



  !> container for a single user-defined observation or platform type
  TYPE, PUBLIC :: letkf_obsplatdef
     CHARACTER(len=10)  :: name      !< unique short name of observation/plat type
     INTEGER            :: id        !< unique integer identifier, used in obs types
     CHARACTER(len=200) :: name_long !< longer name / description
  END TYPE letkf_obsplatdef



  !> Abstract base class for observation file reading and writing.
  !! All user-defined, and built-in, observation file I/O classes for
  !! specific file types should inherit this class
  TYPE, ABSTRACT, PUBLIC :: letkf_obsio
   CONTAINS
     PROCEDURE(I_letkf_obsio_getstr), NOPASS, DEFERRED :: name  !< Get the unique name of the reader
     PROCEDURE(I_letkf_obsio_getstr), NOPASS, DEFERRED :: desc  !< Get a descriptive string for the reader
     PROCEDURE(I_letkf_obsio_init),           DEFERRED :: init  !< initialize the class (usually reading its section of the namelist
     PROCEDURE(I_letkf_obsio_read_obs),       DEFERRED :: read_obs
     PROCEDURE(I_letkf_obsio_read_hx),        DEFERRED :: read_hx
  END TYPE letkf_obsio

  ABSTRACT INTERFACE
     FUNCTION I_letkf_obsio_getstr()
       IMPORT letkf_obsio
       CHARACTER(:), ALLOCATABLE :: I_letkf_obsio_getstr
     END FUNCTION I_letkf_obsio_getstr

     SUBROUTINE I_letkf_obsio_init(self, nml_filename)
       IMPORT letkf_obsio
       CLASS(letkf_obsio) :: self
       CHARACTER(:), ALLOCATABLE, INTENT(in) :: nml_filename
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





  !--------------------------------------------------------------------------------
  ! Public variables
  !--------------------------------------------------------------------------------



  !--------------------------------------------------------------------------------
  ! Private types
  !--------------------------------------------------------------------------------

  !< simple wrapper for letkf_obsio so that we can have an arary of pointers
  TYPE obsio_ptr
     CLASS(letkf_obsio), POINTER :: p
  END TYPE obsio_ptr



  !--------------------------------------------------------------------------------
  ! Private variables
  !--------------------------------------------------------------------------------
  ! configuration paramters
  TYPE(letkf_obsplatdef), ALLOCATABLE :: obsdef_list(:)  !< list of all observation types
  TYPE(letkf_obsplatdef), ALLOCATABLE :: platdef_list(:) !< list of all platform types

  ! mpi derived types
  INTEGER :: observation_mpi_type   !< MPI derived type instantiation

  ! registration of built-in and user-defined obsio classes
  INTEGER, PARAMETER  :: obsio_reg_max = 100      !< max classes to handle
  INTEGER             :: obsio_reg_num = 0        !< current number of classes registered
  TYPE(obsio_ptr)     :: obsio_reg(obsio_reg_max) !< list of registered classes
  CLASS(letkf_obsio), POINTER :: obsio_class      !< The actual I/O module selected

  ! stored observations
  TYPE(letkf_observation), PUBLIC, PROTECTED, ALLOCATABLE :: obs_def(:)
  REAL,                    PUBLIC, PROTECTED, ALLOCATABLE :: obs_hx(:,:)
  REAL,                    PUBLIC, PROTECTED, ALLOCATABLE :: obs_hx_mean(:)

  TYPE(kd_root) :: obs_tree

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------


CONTAINS



  !--------------------------------------------------------------------------------

  !> Initialize the letkf_obs module.
  !! Reads in our section of the namelist, and any observation definition and platorm
  !! definition configuration files.
  SUBROUTINE letkf_obs_init(nml_filename)
    CHARACTER(:), ALLOCATABLE, INTENT(in) :: nml_filename

    INTEGER :: unit, i
    CHARACTER(:), ALLOCATABLE :: ioclass, obsdef_file, platdef_file

    NAMELIST /letkf_obs/ ioclass, obsdef_file, platdef_file

    ! print header
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "============================================================"
       PRINT *, " letkf_obs_init() : observation module initialization"
       PRINT *, "============================================================"
    END IF

    ! read in our section of the namelist
    ALLOCATE(CHARACTER(1024) :: ioclass);       WRITE (ioclass,*) "UNDEFINED"
    ALLOCATE(CHARACTER(1024) :: obsdef_file);   WRITE (obsdef_file,'(A)') "obsdef.cfg"
    ALLOCATE(CHARACTER(1024) :: platdef_file);  WRITE (platdef_file, '(A)') "platdef.cfg"
    OPEN(newunit=unit, file=nml_filename, status='OLD')
    READ(unit, nml=letkf_obs)
    CLOSE(unit)
    ioclass = toupper(TRIM(ioclass))
    obsdef_file = TRIM(obsdef_file)
    platdef_file = TRIM(platdef_file)
    IF (pe_isroot) PRINT letkf_obs

    ! print a list of all obsio classes that have been registered
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "List of obsio classes registered:"
       DO i=1, obsio_reg_num
          PRINT *, " * ", toupper(obsio_reg(i)%p%name()), &
               "  (", obsio_reg(i)%p%desc(), ")"
       END DO
       PRINT *, ""
    END IF

    ! determine which io class to use
    NULLIFY(obsio_class)
    DO i=1,obsio_reg_num
       IF (obsio_reg(i)%p%name() == ioclass) THEN
          obsio_class => obsio_reg(i)%p
          EXIT
       END IF
    END DO
    IF (.NOT. ASSOCIATED(obsio_class)) THEN
       CALL letkf_mpi_abort("obsio class "//ioclass// " not found.")
    END IF

    ! read in the obs/plat configuration files
    IF(pe_isroot) THEN
       PRINT *, "observation definition file"
       PRINT *, "------------------------------------------------------------"
    END IF
    CALL obsplatdef_read(obsdef_file, obsdef_list)

    IF(pe_isroot) THEN
       PRINT *, "platform definition file"
       PRINT *, "------------------------------------------------------------"
    END IF
    CALL obsplatdef_read(platdef_file, platdef_list)
    IF(pe_isroot) PRINT *, ""

    ! initialize MPI object for later sending/receving obsservations
    CALL init_mpi_observation()

    ! finish initialize  of the obsio class
    IF (pe_isroot) PRINT *, "Intializing I/O class: ",obsio_class%name()
    CALL obsio_class%init(nml_filename)

  END SUBROUTINE letkf_obs_init



  !--------------------------------------------------------------------------------

  !< Initializze the mpi derived type for "letkf_observation"
  SUBROUTINE init_mpi_observation()
    INTEGER, PARAMETER :: n = 9
    INTEGER :: TYPE(n), blocklen(n)
    INTEGER(kind=mpi_address_kind) :: disp(n), base
    INTEGER :: ierr, i
    TYPE(letkf_observation) :: ob

    IF (pe_isroot) THEN
       PRINT *, "Initializing MPI derived type (letkf_observation)."
    END IF

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

    ! construct the mpi derived type
    CALL mpi_type_create_struct(n, blocklen, disp, TYPE, observation_mpi_type, ierr)
    CALL mpi_type_commit(observation_mpi_type, ierr)

  END SUBROUTINE init_mpi_observation



  !--------------------------------------------------------------------------------

  !> read the observations
  SUBROUTINE letkf_obs_read()
    INTEGER :: i, nobs, ierr
    REAL, ALLOCATABLE :: obs_lats(:), obs_lons(:), tmp_r(:)

    CALL timing_start('read_obs')

    IF(pe_isroot) THEN
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
    CALL mpi_bcast(obs_def, nobs, observation_mpi_type, pe_root, letkf_mpi_comm, ierr)

    ! read in the observation operator
    ! TODO, parallelize this
    ALLOCATE(obs_hx(ens_size, nobs))
    IF (pe_isroot) THEN
       ALLOCATE(tmp_r(nobs))
       DO i=1,ens_size
          CALL obsio_class%read_hx(i, tmp_r)
          obs_hx(i,:) = tmp_r
       END DO
       DEALLOCATE(tmp_r)

       ! TODO make sure nobs agrees for all files

    END IF
    CALL mpi_bcast(obs_hx, nobs*ens_size, mpi_real, pe_root, letkf_mpi_comm, ierr)

    ! TODO make sure bad obs are removed, do extra QC checks

    ! calculate hx_mean
    ALLOCATE(obs_hx_mean(nobs))
    obs_hx_mean = SUM(obs_hx, 1)/ens_size

    ! remove mean from hx
    DO i=1,nobs
       obs_hx(:,i) = obs_hx(:,i) - obs_hx_mean(i)
    END DO

    
    ! print out observation statistics
    if (pe_isroot) CALL obs_print_stats(obs_def)
 
    ! TODO, make sure no bad qc obs go into the tree

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

    CALL timing_stop('read_obs')

  END SUBROUTINE letkf_obs_read



  !--------------------------------------------------------------------------------

  !> register user-defined and built-in observation I/O classes to be
  !! available for use by LETKF. Actual observation I/O class to be used
  !! is specified in the namelist.
  SUBROUTINE letkf_obs_register(ioclass)
    CLASS(letkf_obsio), POINTER :: ioclass !< the observation I/O class to be registered
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
          IF (toupper(obsio_reg(i)%p%name()) == toupper(ioclass%name())) THEN
             CALL letkf_mpi_abort("can't register obsio class '"// &
                  toupper(ioclass%name())// &
                  "', a class by that name already has been registered.")
          END IF
       END DO
    END IF

    ! add in the class to the list
    obsio_reg_num = obsio_reg_num + 1
    obsio_reg(obsio_reg_num)%p => ioclass

  END SUBROUTINE letkf_obs_register



  !--------------------------------------------------------------------------------
  !> For a given search lat/lon point and radius (slat,slon,sradius), return all
  !! the points that are within that radius. An array of indexes (robs) are returned
  !! which point to the appropriate observations in the module's obs arrays
  !! (obs_def, obs_hx, obs_hx_mean)
  SUBROUTINE letkf_obs_get(slat, slon, sradius, robs, rdist, rnum)
    REAL,    INTENT(in)  :: slat     !< Center of search latitude (degrees)
    REAL,    INTENT(in)  :: slon     !< Center of search longitudee (degrees)
    REAL,    INTENT(in)  :: sradius  !< search radius (meters)
    INTEGER, INTENT(out) :: robs(:)  !< indexes to the observations found
    REAL,    INTENT(out) :: rdist(:) !< distance of each observation (meters)
    INTEGER, INTENT(out) :: rnum     !< number of observations found

    CALL kd_search_radius(obs_tree, slon, slat, sradius, robs, rdist, rnum, .FALSE.)

  END SUBROUTINE letkf_obs_get



  !--------------------------------------------------------------------------------

  !> read the obsdef or platdef configuration files
  SUBROUTINE obsplatdef_read(file, array)
    CHARACTER(len=*), INTENT(in) :: file
    TYPE(letkf_obsplatdef), ALLOCATABLE, INTENT(out) :: array(:)
    LOGICAL :: ex
    INTEGER, PARAMETER :: MAX_DEF = 1024
    INTEGER :: i, pos, unit, iostat
    CHARACTER(len=1024) :: line
    TYPE(letkf_obsplatdef) :: array_temp(MAX_DEF)

    ! make sure the file exists
    INQUIRE(file=file, exist=ex)
    IF (.NOT. ex) THEN
       CALL letkf_mpi_abort("File does not exist: "//file)
    END IF

    ! open it up for reading
    pos = 0
    OPEN(newunit=unit, file=file, action='read')
    DO WHILE(.TRUE.)
       ! read a new line
       READ(unit, '(A)', iostat=iostat) line
       IF (iostat < 0) EXIT
       IF (iostat > 0) THEN
          CALL letkf_mpi_abort("Problem reading file "//file)
       END IF

       ! convert tabs to spaces
       DO i=1,LEN(line)
          IF(line(i:i) == CHAR(9)) line(i:i) = ' '
       END DO

       ! ignore comments
       line = ADJUSTL(line)
       IF (line(1:1) == '#') CYCLE

       ! ignore empty lines
       IF (LEN(TRIM(ADJUSTL(line))) == 0) CYCLE

       ! read in
       IF(pos >= MAX_DEF) CALL letkf_mpi_abort( "MAX_DEF reached")
       pos = pos + 1
       READ(line, *, iostat=iostat) array_temp(pos)
       array_temp(pos)%name = TRIM(toupper(array_temp(pos)%name))
       array_temp(pos)%name_long = TRIM(array_temp(pos)%name_long)

    END DO
    CLOSE(unit)

    ! shrink to final array
    ALLOCATE(array(pos))
    array(1:pos) = array_temp(1:pos)

    ! write summary
    IF (pe_isroot) THEN
       PRINT "(A8, A8, A5, A)", "NAME","ID","","DESCRIPTION"
       DO i=1,pos
          PRINT "(A10,I6,A5,A)", TRIM(array(i)%name), array(i)%id, &
               "", TRIM(array(i)%name_long)
       END DO
       PRINT *, ""
    END IF
  END SUBROUTINE obsplatdef_read



  !--------------------------------------------------------------------------------

  !< print out statistics about the loaded observations
  SUBROUTINE obs_print_stats(obs_t)
    TYPE(letkf_observation) :: obs_t(:)

    INTEGER, ALLOCATABLE :: obst_count(:,:) 
    INTEGER, ALLOCATABLE :: plat_count(:,:)

    type(running_stats), allocatable :: odep_stats(:), pdep_stats(:)
    type(running_stats), allocatable :: osprd_stats(:), psprd_stats(:)
    
    INTEGER :: cnt, i, j, cnt_total
    real :: d,d2, r(ens_size)

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

  

  !--------------------------------------------------------------------------------
  function letkf_obs_getdef(obs_plat, name) result(res)
    character(len=1), intent(in) :: obs_plat
    character(len=*), intent(in) :: name
    type(letkf_obsplatdef) :: res
    
    character(:), allocatable :: name0
    integer :: i
    
    name0 = trim(toupper(name))
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


END MODULE letkf_obs
