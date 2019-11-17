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
  USE running_stats_mod
  USE letkf_config
  USE letkf_mpi
  USE letkf_util

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
  PUBLIC :: letkf_obs_getdef



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
  !>
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC :: letkf_obsplatdef_list
    PRIVATE
      TYPE(letkf_obsplatdef), ALLOCATABLE :: list(:)
    CONTAINS
      PROCEDURE :: count => obsplatdef_list_count
      PROCEDURE :: add => obsplatdef_list_add
      PROCEDURE :: get => obsplatdef_list_get
      PROCEDURE :: set => obsplatdef_list_set
      PROCEDURE :: has => obsplatdef_list_has
  END TYPE letkf_obsplatdef_list
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

     SUBROUTINE I_letkf_obsio_init(self, config, obsdef, platdef)
       IMPORT letkf_obsio
       IMPORT configuration
       IMPORT letkf_obsplatdef
       IMPORT letkf_obsplatdef_list
       CLASS(letkf_obsio) :: self
       TYPE(configuration), INTENT(in) :: config
       TYPE(letkf_obsplatdef_list), INTENT(out) :: obsdef
       TYPE(letkf_obsplatdef_list), INTENT(out) :: platdef
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
       REAL, ALLOCATABLE, INTENT(inout) :: hx(:)
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
  TYPE(letkf_obsplatdef_list) :: obsdef_list

  !> list of all platform types
  TYPE(letkf_obsplatdef_list) :: platdef_list


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
    TYPE(configuration), INTENT(in) :: config

    INTEGER :: i
    CHARACTER(:), ALLOCATABLE :: ioclass
    TYPE(letkf_obsplatdef) :: obsplatdef_tmp

    ! print header
    IF (pe_isroot) THEN
       PRINT *, ""
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
          PRINT *, ' * "', str_tolower(obsio_reg(i)%p%name()), &
               '"  (', obsio_reg(i)%p%desc(), ")"
       END DO
       PRINT *, ""
    END IF

    ! determine which io class to use
    CALL config%get("class", ioclass)
    ioclass = str_tolower(ioclass)
    IF (pe_isroot) PRINT *, "observation.class= "//ioclass
    NULLIFY(obsio_class)
    DO i=1,obsio_reg_num
       IF (str_tolower(obsio_reg(i)%p%name()) == ioclass) THEN
          obsio_class => obsio_reg(i)%p
          EXIT
       END IF
    END DO
    IF (.NOT. ASSOCIATED(obsio_class)) THEN
       CALL letkf_mpi_abort("obsio class "//ioclass// " not found.")
    END IF

    ! finish initialize  of the obsio class
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "Intializing I/O ioclass: ",obsio_class%name()
    END IF
    CALL obsio_class%init(config, obsdef_list, platdef_list)

    ! initialize MPI object for later sending/receving obsservations
    CALL init_mpi_observation()

    ! print the observation definition configuration
    IF(pe_isroot) THEN
       PRINT *, ""
       PRINT *, "observation definitions"
       PRINT *, "------------------------------------------------------------"
       PRINT "(A8,A8,A5,A)", "NAME", "ID", "","DESCRIPTION"
    END IF
    IF (obsdef_list%count() <= 0) &
      CALL letkf_mpi_abort("no obsdef list returned by the ioclass")
    DO i=1,obsdef_list%count()
      ! make sure the name is lower case
      obsplatdef_tmp = obsdef_list%get(i)
      obsplatdef_tmp%name = str_tolower(obsplatdef_tmp%name)
      CALL obsdef_list%set(i, obsplatdef_tmp)

      ! print out information
      IF(pe_isroot) &
            PRINT "(A10, I6, A5, A)", TRIM(obsplatdef_tmp%name), &
              obsplatdef_tmp%id, "", TRIM(obsplatdef_tmp%name_long)
    END DO


    ! print the platform definition configuration
    IF(pe_isroot) THEN
       PRINT *, ""
       PRINT *, "platform definitions"
       PRINT *, "------------------------------------------------------------"
       PRINT "(A8,A8,A5,A)", "NAME", "ID", "","DESCRIPTION"
    END IF
    IF (platdef_list%count() <= 0) &
      CALL letkf_mpi_abort("no platdef list returned by the ioclass")
    DO i=1,platdef_list%count()
      ! make sure the name is lower case
      obsplatdef_tmp = platdef_list%get(i)
      obsplatdef_tmp%name = str_tolower(obsplatdef_tmp%name)
      CALL platdef_list%set(i, obsplatdef_tmp)

      ! print out information
      IF(pe_isroot) &
           PRINT "(A10, I6, A5, A)", TRIM(obsplatdef_tmp%name), obsplatdef_tmp%id, &
           "", TRIM(obsplatdef_tmp%name_long)
    END DO
    IF(pe_isroot) PRINT *, ""

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
    CALL mpi_type_get_extent(mpi_real, lb, ex_real, ierr)
    lb = 0
    ex = ex_real*ens_size
    CALL mpi_type_create_resized(mpi_real, lb, ex, mpitype_real_nk, ierr)
    CALL mpi_type_commit(mpitype_real_nk, ierr)

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
    INTEGER :: obshx_pe(ens_size), obs_pe

    CALL timing_start('read_obs')

    IF(pe_isroot) THEN
       PRINT *, ""
       PRINT *, "reading observations"
       PRINT *, "------------------------------------------------------------"
    END IF

    obs_pe = 0
    IF(pe_rank == obs_pe) THEN
       ! have the I/O class read the main observations file
       ! TODO, allow ioclass to distribute this?
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
    DO i = 1, ens_size
       obshx_pe(i) = letkf_mpi_nextio()
    END DO

    ! read in the observation operator
    !> \todo do this in parallel with the above obs read
    ALLOCATE(obs_hx(ens_size, nobs))
    ALLOCATE(tmp_r(nobs))
    obs_hx = 0.0
    DO i=1,ens_size
       IF(pe_rank == obshx_pe(i)) THEN
          CALL obsio_class%read_hx(i, tmp_r)
          IF ( .NOT. ALLOCATED(tmp_r)) &
                CALL letkf_mpi_abort("obsio class returned unallocated hx array")
          IF ( size(tmp_r) /= nobs ) &
                CALL letkf_mpi_abort("obsio class returned hx array of wrong size")
          obs_hx(i,:) = tmp_r
       END IF
    END DO
    ! broadcast the obs_hx values
    DO i=1,ens_size
       CALL mpi_bcast(obs_hx(i,1), nobs, mpitype_real_nk, obshx_pe(i), letkf_mpi_comm, ierr)
    END DO


    !> \todo make sure nobs agrees for all files


    !> \todo do extra QC checks


    ! calculate hx_mean
    ALLOCATE(obs_hx_mean(nobs))
    obs_hx_mean = SUM(obs_hx, 1)/ens_size

    ! remove mean from hx
    DO i=1,nobs
       obs_hx(:,i) = obs_hx(:,i) - obs_hx_mean(i)
    END DO


    ! print out observation statistics
    IF (pe_isroot) THEN
       IF (nobs > 0) THEN
          CALL obs_print_stats(obs_def)
       ELSE
          IF (pe_isroot) PRINT *, "WARNING: there are NO observations to assimilate"
       END IF
    END IF


    ! remove all bad observations, by shifting good obs down in the list, and decrementing "nobs"
    nobs=0
    DO i=1,size(obs_def)
       IF (obs_def(i)%qc > 0) CYCLE
       nobs = nobs +1
       IF(nobs == i) CYCLE
       obs_def(nobs)=obs_def(i)
       obs_hx(:,nobs)=obs_hx(:,i)
       obs_hx_mean(nobs)=obs_hx_mean(i)
    END DO


    ! add obs to KD tree
    IF (nobs > 0) THEN
       ALLOCATE(obs_lons(nobs))
       ALLOCATE(obs_lats(nobs))
       DO i=1,nobs
          obs_lons(i) = obs_def(i)%lon
          obs_lats(i) = obs_def(i)%lat
       END DO
       CALL kd_init(obs_tree, obs_lons, obs_lats)
       DEALLOCATE(obs_lons)
       DEALLOCATE(obs_lats)
    ELSE
    END IF

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
          IF (str_tolower(obsio_reg(i)%p%name()) == str_tolower(ioclass%name())) THEN
             CALL letkf_mpi_abort("can't register obsio class '"// &
                  str_tolower(ioclass%name())// &
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
    IF (SIZE(obs_def) > 0) THEN
       CALL kd_search_radius(obs_tree, slon, slat, sradius, robs, rdist, rnum, .FALSE.)
    END IF

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

    TYPE(running_stats), ALLOCATABLE :: odep_stats(:), pdep_stats(:)
    TYPE(running_stats), ALLOCATABLE :: osprd_stats(:), psprd_stats(:)

    INTEGER :: cnt, i, j, cnt_total
    TYPE(letkf_obsplatdef) :: obsplatdef

    IF(.NOT. pe_isroot) RETURN

    cnt_total=0
    ALLOCATE(obst_count(obsdef_list%count()+1, 4))
    ALLOCATE(plat_count(platdef_list%count()+1, 4))
    ALLOCATE(odep_stats(obsdef_list%count()))
    ALLOCATE(pdep_stats(platdef_list%count()))
    ALLOCATE(osprd_stats(obsdef_list%count()))
    ALLOCATE(psprd_stats(platdef_list%count()))

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
       DO j=1, obsdef_list%count()
          obsplatdef = obsdef_list%get(j)
          IF(obs_t(i)%obsid == obsplatdef%id) THEN
             obst_count(j,1) = obst_count(j,1) + 1
             obst_count(j,cnt) = obst_count(j,cnt) +1

             ! departure stats if this is a good ob
             IF (obs_t(i)%qc == 0) THEN
                CALL odep_stats(j)%add( obs_def(i)%val-obs_hx_mean(i))
                CALL osprd_stats(j)%add(SUM(obs_hx(:,i)*obs_hx(:,i))/ens_size)
             END IF

             EXIT
          END IF
       END DO
       IF(j > obsdef_list%count()) obst_count(obsdef_list%count()+1,1) = &
            obst_count(obsdef_list%count()+1,1) + 1

       ! count by plat type
       DO j=1, platdef_list%count()
          obsplatdef = platdef_list%get(j)
          IF(obs_t(i)%platid == obsplatdef%id) THEN
             plat_count(j,1) = plat_count(j,1) + 1
             plat_count(j,cnt) = plat_count(j,cnt) +1

             ! departure stats if this is a good ob
             IF (obs_t(i)%qc == 0) THEN
                CALL pdep_stats(j)%add( obs_def(i)%val-obs_hx_mean(i))
                CALL psprd_stats(j)%add(SUM(obs_hx(:,i)*obs_hx(:,i))/ens_size)
             END IF

             EXIT
          END IF
       END DO
       IF(j > platdef_list%count()) plat_count(platdef_list%count()+1,1) = &
            plat_count(platdef_list%count()+1,1) + 1
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
          obsplatdef = obsdef_list%get(i)
          PRINT '(A10,I10, 2I10, I10, A2, F5.1,A)', &
               TRIM(obsplatdef%name), obst_count(i,1), &
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
          obsplatdef = platdef_list%get(i)
          PRINT '(A10,I10, 2I10, I10, A2, F5.1,A)', &
               TRIM(obsplatdef%name), plat_count(i,1), &
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
          obsplatdef = obsdef_list%get(i)
          PRINT '(A10,6F12.5)', &
               TRIM(obsplatdef%name), SQRT(odep_stats(i)%mean(2)), &
               odep_stats(i)%mean(1), odep_stats(i)%MIN(), odep_stats(i)%MAX(), &
               SQRT(osprd_stats(i)%mean())
       END IF
    END DO


    PRINT *, ""
    PRINT *, '         -------------------------------------------------------'
    DO i=1,SIZE(plat_count,1)
       IF (plat_count(i,1) == 0) CYCLE
       IF (i < SIZE(plat_count,1)) THEN
          obsplatdef = platdef_list%get(i)
          PRINT '(A10,5F12.5)', &
               TRIM(obsplatdef%name), SQRT(pdep_stats(i)%mean(2)), &
               pdep_stats(i)%mean(1), pdep_stats(i)%MIN(), pdep_stats(i)%MAX(), &
               SQRT(psprd_stats(i)%mean())
       END IF
    END DO

  END SUBROUTINE obs_print_stats
  !================================================================================



  !================================================================================
  !> Get and observation or platform definition information.
  !--------------------------------------------------------------------------------
  FUNCTION letkf_obs_getdef(obs_plat, name) RESULT(res)
    CHARACTER(len=1), INTENT(in) :: obs_plat !< set to 'P' or 'O' if a platform or
    !! observation definiont is requested
    CHARACTER(len=*), INTENT(in) :: name     !< platform or observation name to get
    TYPE(letkf_obsplatdef) :: res            !< the returned definition

    CHARACTER(len=:), ALLOCATABLE :: name0
    INTEGER :: i

    name0 = TRIM(str_tolower(name))
    i =1
    IF (obs_plat == 'O') THEN
       DO WHILE (i <= obsdef_list%count())
          res = obsdef_list%get(i)
          IF(res%name == name0) EXIT
          i = i + 1
       END DO
       IF ( i > obsdef_list%count()) &
            CALL letkf_mpi_abort("observation definition for '"//name0//"' not found")

    ELSE IF (obs_plat == 'P') THEN
       DO WHILE (i <= platdef_list%count())
          res = platdef_list%get(i)
          IF(res%name == name0) EXIT
          i = i + 1
       END DO
       IF ( i > platdef_list%count()) &
            CALL letkf_mpi_abort("platform definition for '"//name0//"' not found")

    ELSE
       CALL letkf_mpi_abort("letkf_obs_getdef: first argument must be 'P' or 'O'")

    END IF

  END FUNCTION letkf_obs_getdef
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION obsplatdef_list_count(self) RESULT(ret)
    CLASS(letkf_obsplatdef_list) :: self
    INTEGER :: ret

    IF(.NOT. ALLOCATED(self%list)) THEN
      ret = 0
    ELSE
      ret = SIZE(self%list)
    END IF
  END FUNCTION obsplatdef_list_count
  !================================================================================



  !================================================================================
  !>
  ! NOTE: this is inefficient as the list is always dealloacted and recreated
  ! should be using a linked list or something similar.
  !--------------------------------------------------------------------------------
  SUBROUTINE obsplatdef_list_add(self, item, allow_duplicate)
    CLASS(letkf_obsplatdef_list) :: self
    TYPE(letkf_obsplatdef), INTENT(in) :: item
    LOGICAL, INTENT(in), OPTIONAL :: allow_duplicate

    INTEGER :: n
    TYPE(letkf_obsplatdef), ALLOCATABLE :: new_list(:)
    LOGICAL :: allow_duplicate0

    allow_duplicate0=.FALSE.
    IF(PRESENT(allow_duplicate)) allow_duplicate0 = allow_duplicate

    ! make sure the defintion doesn't already exist, if so, ignore
    IF ( ALLOCATED(self%list)) THEN
      DO n=1, SIZE(self%list)
        IF (self%list(n)%name == item%name) THEN
          ! TODO, check to make sure they are fully the same
          IF ( allow_duplicate0 ) THEN
            RETURN
          ELSE
            CALL letkf_mpi_abort("duplicate in obsplatdef_list%add()")
          END IF
        END IF
      END DO
    END IF

    ! allocate the new list
    IF (.NOT. ALLOCATED(self%list)) THEN
      ALLOCATE(new_list(1))
    ELSE
      n = SIZE(self%list)
      ALLOCATE(new_list(n+1))
      new_list(1:n) = self%list
    END IF

    ! append to list
    new_list(SIZE(new_list)) = item

    ! swap lists
    IF(ALLOCATED(self%list)) DEALLOCATE(self%list)
    self%list = new_list
  END SUBROUTINE obsplatdef_list_add
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION obsplatdef_list_get(self, idx) RESULT(ret)
    CLASS(letkf_obsplatdef_list) :: self
    INTEGER, INTENT(in) :: idx
    TYPE(letkf_obsplatdef) :: ret

    IF ( .NOT. ALLOCATED(self%list) ) CALL letkf_mpi_abort( &
      "ERROR in obsplatdef_list%set(), unallocated list")
    IF ( idx > SIZE(self%list)) CALL letkf_mpi_abort( &
        "ERROR in obsplatdef_list%set(), index out of range")

    ret = self%list(idx)
  END FUNCTION obsplatdef_list_get
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION obsplatdef_list_has(self, name) RESULT(ret)
    CLASS(letkf_obsplatdef_list) :: self
    CHARACTER(len=*), INTENT(in) :: name
    LOGICAL :: ret

    INTEGER :: i

    ret = .FALSE.
    IF ( .NOT. ALLOCATED(self%list) ) RETURN

    DO i=1,self%count()
      ! TODO to lower
      IF( TRIM(self%list(i)%name) == TRIM(name)) THEN
        ret = .TRUE.
        EXIT
      END IF
    END DO
  END FUNCTION obsplatdef_list_has
  !================================================================================


  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE obsplatdef_list_set(self, idx, item)
    CLASS(letkf_obsplatdef_list) :: self
    INTEGER, INTENT(in) :: idx
    TYPE(letkf_obsplatdef), INTENT(in) :: item

    IF ( .NOT. ALLOCATED(self%list) ) CALL letkf_mpi_abort( &
      "ERROR in obsplatdef_list%set(), unallocated list")
    IF ( idx > SIZE(self%list)) CALL letkf_mpi_abort( &
        "ERROR in obsplatdef_list%set(), index out of range")

    self%list(idx) = item
  END SUBROUTINE obsplatdef_list_set
  !================================================================================


END MODULE letkf_obs
