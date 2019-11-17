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
!> Module for the obesrvation localization methods
!!
!================================================================================
MODULE letkf_loc
  USE letkf_config
  USE letkf_mpi
  USE letkf_obs

  IMPLICIT NONE
  PRIVATE


  !================================================================================
  !================================================================================
  ! Public components
  !================================================================================
  !================================================================================

  PUBLIC :: letkf_loc_init
  PUBLIC :: letkf_loc_register
  PUBLIC :: letkf_loc_gc


  !================================================================================
  !> contains the definition for a single localization group,
  !! consisting of the list of model state slabs and localization
  !! parameters
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC :: letkf_localizer_group
     INTEGER :: id
     INTEGER, ALLOCATABLE :: slab(:)
  END TYPE letkf_localizer_group
  !================================================================================


  !================================================================================
  !> The abstract base class that all localization specification
  !! classes are to be derived from. Responsible for determining
  !! horizontal / temporal / vertical / variable localization
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC, ABSTRACT :: letkf_localizer
     INTEGER :: maxgroups = -1
   CONTAINS
     PROCEDURE(I_letkf_loc_str),  NOPASS,   DEFERRED :: name
     PROCEDURE(I_letkf_loc_str),  NOPASS,   DEFERRED :: desc
     PROCEDURE(I_letkf_loc_init),           DEFERRED :: init
     PROCEDURE(I_letkf_loc_final),          DEFERRED :: FINAL
     PROCEDURE(I_letkf_loc_groups),         DEFERRED :: groups
     PROCEDURE(I_letkf_loc_maxhz),          DEFERRED :: maxhz
     PROCEDURE(I_letkf_loc_group_localize), DEFERRED :: localize
  END TYPE letkf_localizer

  ABSTRACT INTERFACE

     FUNCTION I_letkf_loc_str()
       CHARACTER(:), ALLOCATABLE :: I_letkf_loc_str
     END FUNCTION I_letkf_loc_str

     SUBROUTINE I_letkf_loc_init(self, config)
       IMPORT letkf_localizer
       IMPORT configuration
       CLASS(letkf_localizer), INTENT(inout) :: self
       TYPE(configuration), INTENT(in) :: config
     END SUBROUTINE I_letkf_loc_init

     SUBROUTINE I_letkf_loc_final(self)
       IMPORT letkf_localizer
       CLASS(letkf_localizer), INTENT(inout) :: self
     END SUBROUTINE I_letkf_loc_final

     SUBROUTINE I_letkf_loc_groups(self, ij, groups)
       IMPORT letkf_localizer
       IMPORT letkf_localizer_group
       CLASS(letkf_localizer), INTENT(inout)  :: self
       INTEGER, INTENT(in)    :: ij      !< gridpoint index to localize around
       TYPE(letkf_localizer_group), ALLOCATABLE, INTENT(inout) :: groups(:)
     END SUBROUTINE I_letkf_loc_groups

     FUNCTION I_letkf_loc_maxhz(self, ij) RESULT(dist)
       IMPORT letkf_localizer
       CLASS(letkf_localizer), INTENT(inout) :: self
       INTEGER, INTENT(in) :: ij
       REAL :: dist
     END FUNCTION I_letkf_loc_maxhz

     FUNCTION I_letkf_loc_group_localize(self, ij, group, obs, dist) RESULT(loc)
       IMPORT letkf_localizer
       IMPORT letkf_localizer_group
       IMPORT letkf_observation
       CLASS(letkf_localizer),      INTENT(inout) :: self
       INTEGER,                     INTENT(in) :: ij
       TYPE(letkf_localizer_group), INTENT(in) :: group
       TYPE(letkf_observation),     INTENT(in) :: obs
       REAL,                        INTENT(in) :: dist
       REAL :: loc
     END FUNCTION I_letkf_loc_group_localize
  END INTERFACE
  !================================================================================


  !================================================================================
  !> simple wrapper for localizer so that we can have arrays of pointers.
  !--------------------------------------------------------------------------------
  TYPE :: localizer_ptr
     CLASS(letkf_localizer), POINTER :: p
  END TYPE localizer_ptr
  !================================================================================


  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  TYPE, PUBLIC :: linearinterp_lat
     INTEGER :: npoints
     REAL, ALLOCATABLE :: lat(:)
     REAL, ALLOCATABLE :: dist(:)
   CONTAINS
     PROCEDURE :: get_dist => linearinterp_lat_get_dist
     PROCEDURE :: print => linearinterp_lat_print
  END TYPE linearinterp_lat
  INTERFACE linearinterp_lat
     MODULE PROCEDURE linearinterp_lat_init
  END INTERFACE linearinterp_lat
  !================================================================================


  ! public variables
  !--------------------------------------------------------------------------------
  CLASS(letkf_localizer), PUBLIC, PROTECTED, POINTER :: localizer_class



  !================================================================================
  !================================================================================
  ! Private module components
  !================================================================================
  !================================================================================
  INTEGER, PARAMETER   :: localizer_reg_max = 100
  INTEGER              :: localizer_reg_num = 0
  TYPE(localizer_ptr)  :: localizer_reg(localizer_reg_max)



CONTAINS



  !================================================================================
  !> initialize
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_loc_init(config)
    TYPE(configuration), INTENT(in) :: config

    INTEGER :: i
    CHARACTER(:), ALLOCATABLE :: loc_class

    IF(pe_isroot) THEN
       PRINT "(//A)", ""
       PRINT *, "============================================================"
       PRINT *, " letkf_loc_init() : localization module initialization"
       PRINT *, "============================================================"
    END IF

    ! print a list of all localizer classes that have been registered
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "List of localizer classes registered:"
       DO i=1, localizer_reg_num
          PRINT *, " * ", toupper(localizer_reg(i)%p%name()), &
               "  (", localizer_reg(i)%p%desc(), ")"
       END DO
       PRINT *, ""
    END IF

    ! determine the loc class to use
    CALL config%get("class", loc_class)
    loc_class=toupper(loc_class)
    IF (pe_isroot) PRINT *, "localization.class= "//loc_class
    NULLIFY(localizer_class)
    DO i=1, localizer_reg_num
       IF (localizer_reg(i)%p%name() == loc_class) THEN
          localizer_class => localizer_reg(i)%p
          EXIT
       END IF
    END DO
    IF (.NOT. ASSOCIATED(localizer_class)) THEN
       CALL letkf_mpi_abort("localizer class "//loc_class//" not found.")
    END IF

    !initialize the localizer class
    CALL localizer_class%init(config)

    ! do some sanity checks
    IF(localizer_class%maxgroups <= 0) &
         CALL letkf_mpi_abort("localizer class needs to properly define 'maxgroups'")

  END SUBROUTINE letkf_loc_init
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_loc_register(locclass)
    CLASS(letkf_localizer), POINTER :: locclass
    INTEGER :: i

    !make sure we haven't reached our max number of classes
    IF (pe_isroot) THEN
       IF(localizer_reg_num == localizer_reg_max) THEN
          CALL letkf_mpi_abort("Too many localizer classes have been registered.")
       END IF
    END IF

    ! make sure a class of this name hasn't already been registered
    IF ( pe_isroot) THEN
       DO i=1, localizer_reg_num
          IF (toupper(localizer_reg(i)%p%name() ) == toupper(locclass%name())) THEN
             CALL letkf_mpi_abort("can't registere localizer class '"//&
                  toupper(locclass%name()) // &
                  "', a class by that name alread has been registered.")
          END IF
       END DO
    END IF

    ! add in the class to the list
    localizer_reg_num = localizer_reg_num + 1
    localizer_reg(localizer_reg_num)%p => locclass
  END SUBROUTINE letkf_loc_register
  !=================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION linearinterp_lat_init(config) RESULT(l)
    TYPE(configuration), INTENT(in) :: config
    TYPE(linearinterp_lat) :: l

    TYPE(CONFIGURATION) :: config2
    INTEGER :: npoints
    INTEGER :: i

    REAL :: lat, radius

    ! determine the number of interp points
    npoints = config%count()

    ! allocate space (2 extra in case 0 and 90 not explicitly defined)
    ALLOCATE(l%lat(npoints+2))
    ALLOCATE(l%dist(npoints+2))
    l%npoints=0

    ! parse the configuration
    DO i=1, config%count()
      CALL config%get(i, config2)
      CALL config2%get("lat", lat)
      CALL config2%get("radius", radius)

      ! stick a 0.0 latitude at the beginning, if needed
      ! (If the first point starts after 0.0)
      IF( lat > 0.0 .AND. l%npoints == 0) THEN
        l%npoints = l%npoints + 1
        l%lat(1) = 0.0
        l%dist(1) = radius
      ENDIF

      l%npoints = l%npoints + 1
      l%lat(l%npoints) = lat
      l%dist(l%npoints) = radius
    END DO

    ! add a 90.0 latitue at the end, if needed
    IF (lat < 90) THEN
      l%npoints = l%npoints +1
      l%lat(l%npoints) = 90
      l%dist(l%npoints) = radius
    END IF
  END FUNCTION linearinterp_lat_init
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION linearinterp_lat_get_dist(self, lat) RESULT(dist)
    CLASS(linearinterp_lat), INTENT(in) :: self
    REAL, INTENT(in) :: lat
    REAL :: dist
    INTEGER :: i

    i=1
    DO WHILE(ABS(lat) > self%lat(i))
       i = i+1
    END DO

    IF(i == 1) THEN
       dist = self%dist(i)
    ELSE
       dist = (ABS(lat)-self%lat(i-1)) / (self%lat(i)-self%lat(i-1))
       dist = dist*self%dist(i) + (1.0-dist)*self%dist(i-1)
    END IF

  END FUNCTION linearinterp_lat_get_dist
  !================================================================================



  !================================================================================
  SUBROUTINE linearinterp_lat_print(self)
    CLASS(linearinterp_lat), INTENT(in) :: self

    CHARACTER(:), ALLOCATABLE :: str
    INTEGER :: i

    DO i = 1,self%npoints
      str = REPEAT(' ', 1024)
      WRITE (str, '(A,F4.1,A,EN10.1E1,A)') "    ", self%lat(i), ' N/S ', self%dist(i),' m'
      PRINT *, TRIM(str)
    END DO

  END SUBROUTINE linearinterp_lat_print
  !================================================================================



  !================================================================================
  !> Gaspari-Cohn localization function.
  !! Possibly faster than the Gaussian function, depending on computer architecture.
  !! Similar shape to Gaussian, except it is compact, goes to 0 at 2L*sqrt(0.3)
  !--------------------------------------------------------------------------------
  PURE FUNCTION letkf_loc_gc(z, L) RESULT(res)
    REAL, INTENT(in) :: z  !< value to localize
    REAL, INTENT(in) :: L  !< the equivalent to the Gaussian standard deviation
    REAL :: res
    REAL(8) :: c
    REAL(8) :: abs_z, z_c

    c = L / SQRT(0.3d0)
    abs_z = ABS(z)
    z_c = abs_z / c

    IF (abs_z >= 2*c) THEN
       res = 0.0
    ELSEIF (abs_z > c) THEN
       res = &
              0.08333d0 * z_c**5 &
            - 0.50000d0 * z_c**4 &
            + 0.62500d0 * z_c**3 &
            + 1.66667d0 * z_c**2 &
            - 5.00000d0 * z_c &
            + 4d0 &
            - 0.66667d0 * c/abs_z
    ELSE
       res = &
             -0.25000d0 * z_c**5 &
            + 0.50000d0 * z_c**4 &
            + 0.62500d0 * z_c**3 &
            - 1.66667d0 * z_c**2 &
            + 1d0
    END IF
  END FUNCTION letkf_loc_gc
  !================================================================================



  ! !================================================================================
  ! !> Gaussian localization function
  ! !--------------------------------------------------------------------------------
  ! PURE FUNCTION letkf_loc_gaus(z, L) RESULT(res)
  !   REAL, INTENT(in) :: z   !< value to localize
  !   REAL, INTENT(in) :: L   !< the gaussian standard deviation
  !   REAL :: res
  !
  !   res = EXP( -0.5 *  z*z / (L*L))
  ! END FUNCTION letkf_loc_gaus
  ! !================================================================================
  !


  !================================================================================
  FUNCTION toupper(in_str) RESULT(out_str)
    CHARACTER(*), INTENT(in) :: in_str
    CHARACTER(LEN(in_str)  ) :: out_str
    INTEGER :: i
    INTEGER, PARAMETER :: offset = 32

    out_str = in_str
    DO i=1, LEN(out_str)
       IF (out_str(i:i) >= "a" .AND. out_str(i:i) <= "z") THEN
          out_str(i:i) = ACHAR(IACHAR(out_str(i:i)) - offset)
       END IF
    END DO
  END FUNCTION toupper
  !================================================================================


END MODULE letkf_loc
