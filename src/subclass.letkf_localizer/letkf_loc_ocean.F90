!================================================================================
!>
!================================================================================
MODULE letkf_loc_ocean
  USE letkf_config
  USE letkf_mpi
  USE letkf_loc
  USE letkf_obs
  USE letkf_state

  IMPLICIT NONE
  PRIVATE

  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  TYPE, EXTENDS(letkf_localizer), PUBLIC :: loc_ocean
     LOGICAL :: vt_split_ml
     TYPE(linearinterp_lat) :: hzdist_prof
     TYPE(linearinterp_lat) :: hzdist_sat
     real :: tloc_prof, tloc_sat
   CONTAINS
     PROCEDURE, NOPASS :: name => loc_ocean_name
     PROCEDURE, NOPASS :: desc => loc_ocean_desc
     PROCEDURE         :: init => loc_ocean_init
     PROCEDURE         :: groups => loc_ocean_groups
     PROCEDURE         :: localize => loc_ocean_localize
     PROCEDURE         :: maxhz => loc_ocean_maxhz
  END TYPE loc_ocean
  !================================================================================



CONTAINS




  !================================================================================
  !> Get the name of the class
  !--------------------------------------------------------------------------------
  FUNCTION loc_ocean_name() RESULT(str)
    CHARACTER(:), ALLOCATABLE :: str
    str="LOC_OCEAN"
  END FUNCTION loc_ocean_name
  !================================================================================



  !================================================================================
  !> Get the description of the class
  !--------------------------------------------------------------------------------
  FUNCTION loc_ocean_desc() RESULT(str)
    CHARACTER(:), ALLOCATABLE :: str
    str="Ocean specific localization."
  END FUNCTION loc_ocean_desc
  !================================================================================



  !================================================================================
  !> Initialize this class, reading in settings from the namelist mostly
  !--------------------------------------------------------------------------------
  SUBROUTINE loc_ocean_init(self, config)
    CLASS(loc_ocean), INTENT(inout) :: self
    TYPE(configuration), INTENT(in) :: config

    CHARACTER(:), ALLOCATABLE :: str

    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "LOC_OCEAN initialization"
       PRINT *, "------------------------------------------------------------"
    END IF

    CALL config%get("vt_split_ml", self%vt_split_ml,  default=.FALSE.)
    CALL config%get("hzloc_prof", str)
    self%hzdist_prof = linearinterp_lat(str)
    CALL config%get("hzloc_sat", str)
    self%hzdist_sat = linearinterp_lat(str)
    call config%get("tloc_prof", self%tloc_prof, default=-1.0)
    call config%get("tloc_sat",  self%tloc_sat, default=-1.0)
    
    IF(pe_isroot) THEN
!       PRINT *, "loc_ocean.vt_split_ml=", self%vt_split_ml
       print *, "loc_ocean.tloc_prof=", self%tloc_prof
       print *, "loc_ocean.tloc_sat=",  self%tloc_sat
    END IF

  END SUBROUTINE loc_ocean_init
  !================================================================================



  !================================================================================
  !> For a given gridpoint, get the desired localization parameters
  !! (horizontal search radius and level/variable localization groups)
  SUBROUTINE loc_ocean_groups(self, ij, groups)
    CLASS(loc_ocean), INTENT(inout) :: self
    INTEGER, INTENT(in)  :: ij
    TYPE(letkf_localizer_group), ALLOCATABLE, INTENT(inout) :: groups(:)

    INTEGER :: i

    ! TODO, split based on MLD and thermocline level

    ALLOCATE(groups(1))
    ALLOCATE(groups(1)%slab(grid_ns))
    DO i=1, grid_ns
       groups(1)%slab(i) = i
    END DO

  END SUBROUTINE loc_ocean_groups
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION loc_ocean_maxhz(self, ij) RESULT(dist)
    CLASS(loc_ocean), INTENT(inout) :: self
    INTEGER, INTENT(in)  :: ij
    REAL :: dist

    dist = MAX( self%hzdist_prof%get_dist(lat_ij(ij)),&
         self%hzdist_sat%get_dist(lat_ij(ij)))

  END FUNCTION loc_ocean_maxhz
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION loc_ocean_localize(self, ij, group, obs, dist) RESULT(loc)
    CLASS(loc_ocean), INTENT(inout) :: self
    INTEGER, INTENT(in) :: ij
    TYPE(letkf_localizer_group), INTENT(in) :: group
    TYPE(letkf_observation), INTENT(in) :: obs
    REAL, INTENT(in) :: dist
    REAL :: loc


    ! \todo remove hardcoding of plat ids
    IF(obs%platid == 1000) THEN
       loc = letkf_loc_gc(dist, self%hzdist_sat%get_dist(lat_ij(ij)))
       if(self%tloc_sat > 0) &
            loc = loc * letkf_loc_gc(obs%time, self%tloc_sat)
    ELSE IF(obs%platid == 1) THEN
       loc = letkf_loc_gc(dist, self%hzdist_prof%get_dist( lat_ij(ij)))
       if(self%tloc_prof > 0) &
            loc = loc * letkf_loc_gc(obs%time, self%tloc_prof)       
    ELSE
       PRINT *, "PLATID: ", obs%platid
       CALL letkf_mpi_abort("unhandled platid in loc_ocean_localize")
    ENDIF

  END FUNCTION loc_ocean_localize
  !================================================================================


END MODULE letkf_loc_ocean
