MODULE letkf_loc_novrt
  use letkf_config
  USE letkf_mpi
  USE letkf_loc
  USE letkf_obs
  USE letkf_state

  IMPLICIT NONE
  PRIVATE

  TYPE, EXTENDS(letkf_localizer), PUBLIC :: loc_novrt
   CONTAINS
     PROCEDURE, NOPASS :: name => loc_novrt_name
     PROCEDURE, NOPASS :: desc => loc_novrt_desc
     PROCEDURE         :: init => loc_novrt_init
     PROCEDURE         :: groups => loc_novrt_groups
     PROCEDURE         :: localize => loc_novrt_localize
     PROCEDURE         :: maxhz => loc_novrt_maxhz
  END TYPE loc_novrt

  ! TODO get rid of this
  REAL :: hzloc(2)  !< horizontal localization distance(meters)
  !! at the equator, and pole


CONTAINS


  !> Get the name of the class
  FUNCTION loc_novrt_name() RESULT(str)
    CHARACTER(:), ALLOCATABLE :: str
    str="LOC_NOVRT"
  END FUNCTION loc_novrt_name


  !> Get the description of the class
  FUNCTION loc_novrt_desc() RESULT(str)
    CHARACTER(:), ALLOCATABLE :: str
    str="Generic horizontal only localization"
  END FUNCTION loc_novrt_desc



  !> Initialize this class, reading in settings from the namelist mostly
  SUBROUTINE loc_novrt_init(self, config)
    CLASS(loc_novrt), INTENT(inout) :: self
    type(configuration), intent(in) :: config

    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "LOC_NOVRT initialization"
       PRINT *, "------------------------------------------------------------"
    END IF
    hzloc = (/500.0d3, 50.0d3/)
  END SUBROUTINE loc_novrt_init



  !--------------------------------------------------------------------------------

  !> For a given gridpoint, get the desired localization parameters
  !! (horizontal search radius and level/variable localization groups)
  SUBROUTINE loc_novrt_groups(self, ij, groups)
    CLASS(loc_novrt), INTENT(inout) :: self
    INTEGER, INTENT(in)  :: ij
    type(letkf_localizer_group), ALLOCATABLE, INTENT(inout) :: groups(:)

    integer :: i

    ! TODO, need to correctly set the number of slabs
    ALLOCATE(groups(1))
    allocate(groups(1)%slab(grid_ns))
    do i=1, grid_ns
       groups(1)%slab(i) = i
    end do

  END SUBROUTINE loc_novrt_groups



  !--------------------------------------------------------------------------------
  !>
  FUNCTION loc_novrt_maxhz(self, ij) RESULT(dist)
    CLASS(loc_novrt), INTENT(inout) :: self
    INTEGER, INTENT(in)  :: ij
    REAL :: dist

    dist=ABS(lat_ij(ij))/90.0
    dist=dist*hzloc(2) + (1.0-dist)*hzloc(1)

  END FUNCTION loc_novrt_maxhz



  !--------------------------------------------------------------------------------
  !>
  FUNCTION loc_novrt_localize(self, ij, group, obs, dist) RESULT(loc)
    CLASS(loc_novrt), INTENT(inout) :: self
    INTEGER, INTENT(in) :: ij
    type(letkf_localizer_group), INTENT(in) :: group
    TYPE(letkf_observation), INTENT(in) :: obs
    real, intent(in) :: dist
    REAL :: loc, r

    ! TODO temporal localization

    ! horizontal localization
    r = abs(lat_ij(ij))/90.0
    r = r*hzloc(2)+(1.0-r)*hzloc(1)
    loc = letkf_loc_gc(dist, r)

  END FUNCTION loc_novrt_localize


END MODULE letkf_loc_novrt
