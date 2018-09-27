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
     TYPE(linearinterp_lat) :: hzdist_prof
     TYPE(linearinterp_lat) :: hzdist_sat
     REAL :: tloc_prof, tloc_sat

     INTEGER :: vtloc_surf_src = 0

   CONTAINS
     PROCEDURE, NOPASS :: name => loc_ocean_name
     PROCEDURE, NOPASS :: desc => loc_ocean_desc
     PROCEDURE         :: init => loc_ocean_init
     PROCEDURE         :: groups => loc_ocean_groups
     PROCEDURE         :: localize => loc_ocean_localize
     PROCEDURE         :: maxhz => loc_ocean_maxhz
  END TYPE loc_ocean
  !================================================================================


  INTEGER, PARAMETER :: VTLOC_SURF_SRC_NONE = 0 !< no surface ob vt localization
  INTEGER, PARAMETER :: VTLOC_SURF_SRC_BKGT = 1 !< vtloc based on MLD from bakground T
  !INTEGER, PARAMETER :: VTLOC_SURF_SRC_DPTH = 2
  !INTEGER, PARAMETER :: VTLOC_SURF_SRC_LVLS = 3
  !INTEGER, PARAMETER :: VTLOC_SURF_SRC_SPRD = 4

  INTEGER, PARAMETER :: LOC_GROUP_SURF=0  !< surface localization group
  INTEGER, PARAMETER :: LOC_GROUP_SUBML=1 !< below the mixed layer loc group

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

    LOGICAL :: b
    CHARACTER(:), ALLOCATABLE :: str
    TYPE(configuration) :: config2

    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "LOC_OCEAN initialization"
       PRINT *, "------------------------------------------------------------"
    END IF

    CALL config%get("hzloc_prof", str)
    self%hzdist_prof = linearinterp_lat(str)
    CALL config%get("hzloc_sat", str)
    self%hzdist_sat = linearinterp_lat(str)
    CALL config%get("tloc_prof", self%tloc_prof, default=-1.0)
    CALL config%get("tloc_sat",  self%tloc_sat, default=-1.0)

    self%maxgroups = 1
    self%vtloc_surf_src = VTLOC_SURF_SRC_NONE

    IF(config%found("vtloc_surf"))THEN
       CALL config%get("vtloc_surf", config2)
       CALL config2%get("type", str)
       IF(str == "none") THEN
          CONTINUE
       ELSE IF(str == "mld:bkg_t") THEN
          self%maxgroups = 2
          self%vtloc_surf_src = VTLOC_SURF_SRC_BKGT
       ELSE
          CALL letkf_mpi_abort("illegal vtloc_surf 'type' given: "//str)
       END IF
    END IF

    IF(pe_isroot) THEN
       !       PRINT *, "loc_ocean.vt_split_ml=", self%vt_split_ml
       PRINT *, "loc_ocean.tloc_prof=", self%tloc_prof
       PRINT *, "loc_ocean.tloc_sat=",  self%tloc_sat
    END IF

  END SUBROUTINE loc_ocean_init
  !================================================================================



  !================================================================================
  !> For a given gridpoint, get the desired localization parameters
  !! (horizontal search radius and level/variable localization groups)
  !--------------------------------------------------------------------------------
  SUBROUTINE loc_ocean_groups(self, ij, groups)
    CLASS(loc_ocean), INTENT(inout) :: self
    INTEGER, INTENT(in)  :: ij
    TYPE(letkf_localizer_group), ALLOCATABLE, INTENT(inout) :: groups(:)

    INTEGER :: i, j, vtloc_lvl
    INTEGER :: group1_count, group2_count
    INTEGER :: group1_idx(grid_ns), group2_idx(grid_ns)


    ! 2 group vertical localization for surface obs?
    IF(self%vtloc_surf_src /= VTLOC_SURF_SRC_NONE) THEN

       ! determine the surface localization max depth level
       vtloc_lvl=1
       IF(self%vtloc_surf_src == VTLOC_SURF_SRC_BKGT) THEN
          vtloc_lvl = calc_mldlvl_bkgt(self,ij)
       ELSE
          CALL letkf_mpi_abort("no ocean vertical localization other than BKGT implemented yet.")
       END IF

       ! for each variable, split its levels into the two groups
       group1_count = 0
       group2_count = 0
       DO i =1,SIZE(statevars)
          ! single level variables are assumed to be at the surface (for now)
          ! TODO generalize this, use the actual vtgrid data
          DO j=1, statevars(i)%levels
             IF (j<= vtloc_lvl) THEN
                group1_count = group1_count+1
                group1_idx(group1_count) = statevars(i)%grid_s_idx+j-1
             ELSE
                group2_count = group2_count+1
                group2_idx(group2_count) = statevars(i)%grid_s_idx+j-1
             END IF
          END DO
       END DO

       ! copy the values
       ALLOCATE(groups(2))

       groups(1)%id = LOC_GROUP_SURF
       ALLOCATE(groups(1)%slab(group1_count))
       groups(1)%slab = group1_idx(1:group1_count)

       groups(2)%id = LOC_GROUP_SUBML
       ALLOCATE(groups(2)%slab(group2_count))
       groups(2)%slab = group2_idx(1:group2_count)


    ELSE

       ! otherwise, no vertical localization, so only 1 loc group
       ALLOCATE(groups(1))
       ALLOCATE(groups(1)%slab(grid_ns))
       DO i=1, grid_ns
          groups(1)%slab(i) = i
       END DO

    END IF

  END SUBROUTINE loc_ocean_groups
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION calc_mldlvl_bkgt(self, ij) RESULT(mldlvl)
    CLASS(loc_ocean), INTENT(inout) :: self
    INTEGER, INTENT(in)  :: ij
    INTEGER :: mldlvl

    mldlvl = 1

  END FUNCTION calc_mldlvl_bkgt
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


    ! TODO remove hardcoding of plat ids, obs grid_s_idx
    IF(obs%platid == 1000) THEN
       ! satellite observations
       IF(group%id == LOC_GROUP_SUBML) THEN
          loc = 0.0
       ELSE
          loc = letkf_loc_gc(dist, self%hzdist_sat%get_dist(lat_ij(ij)))
          IF(self%tloc_sat > 0) &
               loc = loc * letkf_loc_gc(obs%time, self%tloc_sat)
       END IF

    ELSE IF(obs%platid == 1) THEN
       ! profile T/S observations
       loc = letkf_loc_gc(dist, self%hzdist_prof%get_dist( lat_ij(ij)))
       IF(self%tloc_prof > 0) &
            loc = loc * letkf_loc_gc(obs%time, self%tloc_prof)

    ELSE
       PRINT *, "PLATID: ", obs%platid
       CALL letkf_mpi_abort("unhandled platid in loc_ocean_localize")

    ENDIF

  END FUNCTION loc_ocean_localize
  !================================================================================


END MODULE letkf_loc_ocean
