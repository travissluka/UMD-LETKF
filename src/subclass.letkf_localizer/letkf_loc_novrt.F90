!================================================================================
!> default, simple, horizontal only  localization
!================================================================================
MODULE letkf_loc_novrt
  USE letkf_config
  USE letkf_mpi
  USE letkf_loc
  USE letkf_obs
  USE letkf_state

  IMPLICIT NONE
  PRIVATE



  !================================================================================
  !================================================================================
  ! Public module components
  !================================================================================
  !================================================================================


  !================================================================================
  !> localizer class for simple horizontal only localization
  !--------------------------------------------------------------------------------
  TYPE, EXTENDS(letkf_localizer), PUBLIC :: loc_novrt
     REAL :: hzloc(2) !< horizontal localization distance(meters)
     !! at the equator, and pole
   CONTAINS
     PROCEDURE, NOPASS :: name => loc_novrt_name
     PROCEDURE, NOPASS :: desc => loc_novrt_desc
     PROCEDURE         :: init => loc_novrt_init
     PROCEDURE         :: groups => loc_novrt_groups
     PROCEDURE         :: localize => loc_novrt_localize
     PROCEDURE         :: maxhz => loc_novrt_maxhz
  END TYPE loc_novrt
  !================================================================================



CONTAINS


  !================================================================================
  !> Get the name of the class
  !--------------------------------------------------------------------------------
  FUNCTION loc_novrt_name() RESULT(str)
    CHARACTER(:), ALLOCATABLE :: str
    str="LOC_NOVRT"
  END FUNCTION loc_novrt_name
  !================================================================================



  !================================================================================
  !> Get the description of the class
  !--------------------------------------------------------------------------------
  FUNCTION loc_novrt_desc() RESULT(str)
    CHARACTER(:), ALLOCATABLE :: str
    str="Generic horizontal only localization"
  END FUNCTION loc_novrt_desc
  !================================================================================



  !================================================================================
  !> Initialize this class, reading in settings from the namelist mostly
  !--------------------------------------------------------------------------------
  SUBROUTINE loc_novrt_init(self, config)
    CLASS(loc_novrt), INTENT(inout) :: self
    TYPE(configuration), INTENT(in) :: config

    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "LOC_NOVRT initialization"
       PRINT *, "------------------------------------------------------------"
    END IF
    self%hzloc = (/500.0d3, 50.0d3/)
  END SUBROUTINE loc_novrt_init
  !================================================================================



  !================================================================================
  !> For a given gridpoint, get the desired localization parameters
  !! (horizontal search radius and level/variable localization groups)
  !--------------------------------------------------------------------------------
  SUBROUTINE loc_novrt_groups(self, ij, groups)
    CLASS(loc_novrt), INTENT(inout) :: self
    INTEGER, INTENT(in)  :: ij
    TYPE(letkf_localizer_group), ALLOCATABLE, INTENT(inout) :: groups(:)

    INTEGER :: i

    ! TODO, need to correctly set the number of slabs
    ALLOCATE(groups(1))
    ALLOCATE(groups(1)%slab(grid_ns))
    DO i=1, grid_ns
       groups(1)%slab(i) = i
    END DO

  END SUBROUTINE loc_novrt_groups
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION loc_novrt_maxhz(self, ij) RESULT(dist)
    CLASS(loc_novrt), INTENT(inout) :: self
    INTEGER, INTENT(in)  :: ij
    REAL :: dist

    dist=ABS(lat_ij(ij))/90.0
    dist=dist*self%hzloc(2) + (1.0-dist)*self%hzloc(1)

  END FUNCTION loc_novrt_maxhz
  !================================================================================




  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION loc_novrt_localize(self, ij, group, obs, dist) RESULT(loc)
    CLASS(loc_novrt), INTENT(inout) :: self
    INTEGER, INTENT(in) :: ij
    TYPE(letkf_localizer_group), INTENT(in) :: group
    TYPE(letkf_observation), INTENT(in) :: obs
    REAL, INTENT(in) :: dist
    REAL :: loc, r

    ! TODO temporal localization

    ! horizontal localization
    r = ABS(lat_ij(ij))/90.0
    r = r*self%hzloc(2)+(1.0-r)*self%hzloc(1)
    loc = letkf_loc_gc(dist, r)

  END FUNCTION loc_novrt_localize
  !================================================================================


END MODULE letkf_loc_novrt
