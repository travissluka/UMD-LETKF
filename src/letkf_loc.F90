MODULE letkf_loc
  USE letkf_mpi
  USE letkf_obs

  IMPLICIT NONE
  PRIVATE


  !------------------------------------------------------------
  ! Public methods
  !------------------------------------------------------------
  PUBLIC :: letkf_loc_init
  PUBLIC :: letkf_loc_register
  PUBLIC :: letkf_loc_gc
  PUBLIC :: letkf_loc_gaus



  !------------------------------------------------------------
  ! Public types
  !------------------------------------------------------------

  !> contains the definitiona for a single localization group,
  !! consisting of the list of model state slabs and localization
  !! parameters
  TYPE, PUBLIC :: letkf_localizer_group
     INTEGER, ALLOCATABLE :: slab(:)
  END TYPE letkf_localizer_group


  !> simple wrapper for localizer_group so that we can have arrays
  !! of pointers.
!  TYPE, PUBLIC :: letkf_localizer_group_ptr
!     CLASS(letkf_localizer_group), POINTER :: p
!  END TYPE letkf_localizer_group_ptr


  !> The abstract base class that all localization specification
  !! classes are to be derived from. Responsible for determining
  !! horizontal / temporal / vertical / variable localization
  TYPE, PUBLIC, ABSTRACT :: letkf_localizer
   CONTAINS
     PROCEDURE(I_letkf_loc_str),  NOPASS,   DEFERRED :: name
     PROCEDURE(I_letkf_loc_str),  NOPASS,   DEFERRED :: desc
     PROCEDURE(I_letkf_loc_init),           DEFERRED :: init
     PROCEDURE(I_letkf_loc_groups),         DEFERRED :: groups
     PROCEDURE(I_letkf_loc_maxhz),          DEFERRED :: maxhz
     PROCEDURE(I_letkf_loc_group_localize), DEFERRED :: localize
  END TYPE letkf_localizer

  ABSTRACT INTERFACE

     FUNCTION I_letkf_loc_str()
       CHARACTER(:), ALLOCATABLE :: I_letkf_loc_str
     END FUNCTION I_letkf_loc_str

     SUBROUTINE I_letkf_loc_init(self, nml_filename)
       IMPORT letkf_localizer
       CLASS(letkf_localizer), INTENT(inout) :: self
       CHARACTER(:), ALLOCATABLE, INTENT(in) :: nml_filename
     END SUBROUTINE I_letkf_loc_init

     SUBROUTINE I_letkf_loc_groups(self, ij, groups)
       IMPORT letkf_localizer
       IMPORT letkf_localizer_group
       CLASS(letkf_localizer), INTENT(inout)  :: self
       INTEGER, INTENT(in)    :: ij      !< gridpoint index to localize around
       type(letkf_localizer_group), ALLOCATABLE, INTENT(inout) :: groups(:)
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
       type(letkf_localizer_group), INTENT(in) :: group
       TYPE(letkf_observation),     INTENT(in) :: obs
       real,                        INTENT(in) :: dist
       REAL :: loc
     END FUNCTION I_letkf_loc_group_localize
  END INTERFACE

  !> simple wrapper for localizer so that we can have arrays
  !! of pointers.
  TYPE :: localizer_ptr
     CLASS(letkf_localizer), POINTER :: p
  END TYPE localizer_ptr


  CLASS(letkf_localizer), PUBLIC, PROTECTED, POINTER :: localizer_class

  !--------------------------------------------------------------------------------
  ! Private variables
  !--------------------------------------------------------------------------------
  INTEGER, PARAMETER   :: localizer_reg_max = 100
  INTEGER              :: localizer_reg_num = 0
  TYPE(localizer_ptr)  :: localizer_reg(localizer_reg_max)



CONTAINS



  !--------------------------------------------------------------------------------
  !> initialize
  SUBROUTINE letkf_loc_init(nml_filename)
    CHARACTER(:), ALLOCATABLE,  INTENT(in) :: nml_filename
    INTEGER :: i, unit
    CHARACTER(:), ALLOCATABLE :: loc_class

    NAMELIST /letkf_loc/ loc_class

    IF(pe_isroot) THEN
       PRINT "(//A)", ""
       PRINT *, "============================================================"
       PRINT *, " letkf_loc_init() : localization module initialization"
       PRINT *, "============================================================"
    END IF

    ! read in our section of the namelist
    ALLOCATE(CHARACTER(1024) :: loc_class); WRITE (loc_class, *) "UNDEFINED"
    OPEN(newunit=unit, file=nml_filename, status='OLD')
    READ(unit, nml=letkf_loc)
    loc_class=toupper(TRIM(loc_class))
    CLOSE(unit)

    ! print a list of all localizer classes that have been registered
    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "List of localizer classes registered:"
       DO i=1, localizer_reg_num
          PRINT *, " * ", toupper(localizer_reg(i)%p%name()), &
               "  (", localizer_reg(i)%p%desc(), ")"
       END DO
    END IF

    ! determine the loc class to use
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
    CALL localizer_class%init(nml_filename)
  END SUBROUTINE letkf_loc_init



  !--------------------------------------------------------------------------------
  !>
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

  !--------------------------------------------------------------------------------

  !> Gaspari-Cohn localization function.
  !! Possibly faster than the Gaussian function, depending on computer architecture.
  !! Similar shape to Gaussian, except it is compact, goes to 0 at 2L*sqrt(0.3)
  PURE FUNCTION letkf_loc_gc(z, L) RESULT(res)
    REAL, INTENT(in) :: z  !< value to localize
    REAL, INTENT(in) :: L  !< the equivalent to the Gaussian standard deviation
    REAL :: res
    REAL :: c
    REAL :: abs_z, z_c

    c = L / SQRT(0.3)
    abs_z = ABS(z)
    z_c = abs_z / c

    IF (abs_z >= 2*c) THEN
       res = 0.0
    ELSEIF (abs_z < 2*c .AND. abs_z > c) THEN
       res = &
            (1.0/12.0)*z_c**5 - 0.5*z_c**4 + &
            (5.0/8.0)*z_c**3 + (5.0/3.0)*z_c**2 &
            - 5.0*z_c + 4 - (2.0/3.0)*c/abs_z
    ELSE
       res = &
            -0.25*z_c**5 + 0.5*z_c**4 + &
            (5.0/8.0)*z_c**3 - (5.0/3.0)*z_c**2 + 1
    END IF
  END FUNCTION letkf_loc_gc



  !--------------------------------------------------------------------------------

  !> Gaussian localization function
  ! TODO umm, why is this not working
  PURE FUNCTION letkf_loc_gaus(z, L) RESULT(res)
    REAL, INTENT(in) :: z   !< value to localize
    REAL, INTENT(in) :: L   !< the gaussian standard deviation
    REAL :: res

    res = EXP( -0.5 *  z*z / (L*L))
  END FUNCTION letkf_loc_gaus


  !--------------------------------------------------------------------------------
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
END MODULE letkf_loc
