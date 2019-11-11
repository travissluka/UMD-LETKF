! Copyright 2018-2019 Travis Sluka
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
!>
!================================================================================
MODULE letkf_loc_ocean_mod
  USE netcdf
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
  TYPE, EXTENDS(letkf_localizer), PUBLIC :: letkf_loc_ocean
     LOGICAL :: save_diag = .FALSE.
     CHARACTER(:), ALLOCATABLE :: diag_file

     ! horizontal and temporal localization
     TYPE(linearinterp_lat) :: hzdist_prof
     TYPE(linearinterp_lat) :: hzdist_sat
     REAL :: tloc_prof, tloc_sat

     ! SST vertical localization type, and parameters to go with it
     INTEGER :: vtloc_surf_src = 0
     INTEGER, ALLOCATABLE :: surf_obs(:)
     INTEGER, ALLOCATABLE :: surf_plats(:)
     REAL :: bkg_t_delta
     INTEGER :: bkg_t_var_idx
     INTEGER :: bkg_t_var_len

     ! optional diagnostics to save
     REAL, ALLOCATABLE :: diag_vtloc_surf_lvl(:)
     REAL, ALLOCATABLE :: diag_vtloc_surf_dpth(:)

   CONTAINS
     PROCEDURE, NOPASS :: name => loc_ocean_name
     PROCEDURE, NOPASS :: desc => loc_ocean_desc
     PROCEDURE         :: init => loc_ocean_init
     PROCEDURE         :: FINAL => loc_ocean_final
     PROCEDURE         :: groups => loc_ocean_groups
     PROCEDURE         :: localize => loc_ocean_localize
     PROCEDURE         :: maxhz => loc_ocean_maxhz
  END TYPE letkf_loc_ocean
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
    CLASS(letkf_loc_ocean), INTENT(inout) :: self
    TYPE(configuration), INTENT(in) :: config

    CHARACTER(:), ALLOCATABLE :: str
    TYPE(configuration) :: config2, config3
    TYPE(letkf_statevar_spec) :: state_def
    TYPE(letkf_obsplatdef) :: obsplat_def
    INTEGER :: i

    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "LOC_OCEAN initialization"
       PRINT *, "------------------------------------------------------------"
    END IF

    CALL config%get("save_diag", self%save_diag, .true.)
    CALL config%get("diag_file", self%diag_file, "diag.loc_ocean.nc" )

    ! TODO, this has gotten a bit messy,
    ! generalize this by making a separate hzloc class

    ! hz loc for profiles
    CALL config%get("hzloc_prof", config2)
    CALL config2%get("type", str)
    IF (str == "linearinterp_lat") THEN
      CALL config2%get("value", config3)
      self%hzdist_prof = linearinterp_lat(config3)
    ELSE
      CALL letkf_mpi_abort('unrecognized hzloc_prof type: "'//TRIM(str)//'"')
    ENDIF

    ! hzloc for satellites
    CALL config%get("hzloc_sat", config2)
    CALL config2%get("type", str)
    IF (str == "linearinterp_lat") THEN
      CALL config2%get("value", config3)
      self%hzdist_sat = linearinterp_lat(config3)
    ELSE
      CALL letkf_mpi_abort('unrecognized hzloc_sat type: "'//TRIM(str)//'"')
    ENDIF

    ! temporal localization
    CALL config%get("tloc_prof", self%tloc_prof, default=-1.0)
    CALL config%get("tloc_sat",  self%tloc_sat, default=-1.0)

    IF(pe_isroot) THEN
       PRINT *, "localization.save_diag  = ", self%save_diag
       IF(self%save_diag) PRINT *, "localization.diag_file  = ", self%diag_file
       PRINT *, "localization.tloc_prof  = ", self%tloc_prof
       PRINT *, "localization.tloc_sat   = ", self%tloc_sat
       PRINT *, "localization.hzloc_prof = "
       CALL self%hzdist_prof%print()
       PRINT *, "localization.hzloc_sat  = "
       CALL self%hzdist_sat%print()
    END IF

    ! get type of vertical localization used for SST obs
    self%maxgroups = 1
    self%vtloc_surf_src = VTLOC_SURF_SRC_NONE
    IF(config%found("vtloc_surf"))THEN
       CALL config%get("vtloc_surf", config2)
       CALL config2%get("type", str)

       IF(pe_isroot) PRINT *, 'vtloc_surf.type  = "', str, '"'

       ! determine the parameters of the vertical localization
       IF(str == "none") THEN
          ! No vertical localization
          CONTINUE

       ELSE IF(str == "bkg_t") THEN
          ! localized to MLD based on temperature delta
          self%maxgroups = 2
          self%vtloc_surf_src = VTLOC_SURF_SRC_BKGT
          CALL config2%get("bkg_t_delta", self%bkg_t_delta)
          CALL config2%get("bkg_t_var", str)
          state_def = letkf_state_var_getdef(str)
          self%bkg_t_var_idx = state_def%grid_s_idx
          self%bkg_t_var_len = state_def%levels

          IF(pe_isroot) THEN
             PRINT *, "vtloc_surf.bkg_t_delta = ", self%bkg_t_delta
             PRINT *, 'vtloc_surf.bkg_t_var   = "', str,'"'
          END IF

       ELSE
          ! ERROR, unknown given type
          CALL letkf_mpi_abort("illegal vtloc_surf 'type' given: "//str)
       END IF
    END IF


    ! TODO clean up this code, this was causing problems with hybrid-godas
    ! so i applied a quick fix


    ! determine the list of observations or platform ids that are considered
    ! surface observations that need to be localized. Convert from a specified
    ! string on the config file to the associated integer id.
    IF(config%found("sat_obs")) THEN
       IF(pe_isroot) PRINT *, "satellite observation types:"

       CALL config%get("sat_obs", config3)
       ALLOCATE(self%surf_obs(config3%COUNT()))
       DO i = 1, SIZE(self%surf_obs)
          CALL config3%get(i,str)
          IF (pe_isroot) PRINT *, "  * ", TRIM(str)
          obsplat_def = letkf_obs_getdef('O',str)
          self%surf_obs(i) = obsplat_def%id
       END DO
    ELSE
       ALLOCATE(self%surf_obs(0))
    END IF

    IF(config%found("sat_plats")) THEN
       IF(pe_isroot) PRINT *, "satellite platform types:"

       CALL config%get("sat_plats", config3)
       ALLOCATE(self%surf_plats(config3%COUNT()))
       DO i = 1, SIZE(self%surf_plats)
          CALL config3%get(i,str)
          IF (pe_isroot) PRINT *, "  * ", TRIM(str)
          obsplat_def = letkf_obs_getdef('P',str)
          self%surf_plats(i) = obsplat_def%id
       END DO
    ELSE
       ALLOCATE(self%surf_plats(0))
    END IF

    ! setup optional diagnostics
    IF(self%save_diag) THEN
       ALLOCATE(self%diag_vtloc_surf_lvl(ij_count))
       ALLOCATE(self%diag_vtloc_surf_dpth(ij_count))
       self%diag_vtloc_surf_lvl = 0
       self%diag_vtloc_surf_dpth = 0
    END IF

  END SUBROUTINE loc_ocean_init
  !================================================================================


  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE loc_ocean_final(self)
    CLASS(letkf_loc_ocean), INTENT(inout) :: self

    INTEGER :: ncid, vid
    INTEGER :: d_x, d_y, d_t
    REAL :: tmp2d(grid_nx, grid_ny)

    IF(self%save_diag) THEN
       !setup oiutput netcdf file
       IF(pe_isroot) THEN
          PRINT *, ""
          PRINT *, 'Saving localization diagnostics to "'//self%diag_file//'"...'
          CALL check(nf90_create(self%diag_file, NF90_CLOBBER, ncid))
          CALL check(nf90_def_dim(ncid, "time", 1, d_t))
          CALL check(nf90_def_dim(ncid, "lat", grid_ny, d_y))
          CALL check(nf90_def_dim(ncid, "lon", grid_nx, d_x))

          CALL check(nf90_def_var(ncid, "lat", nf90_real, (/d_y/), vid))
          CALL check(nf90_put_att(ncid, vid, "axis", "Y"))
          CALL check(nf90_put_att(ncid, vid, "units", "degrees_north"))

          CALL check(nf90_def_var(ncid, "lon", nf90_real, (/d_x/), vid))
          CALL check(nf90_put_att(ncid, vid, "axis", "X"))
          CALL check(nf90_put_att(ncid, vid, "units", "degrees_east"))

          CALL check(nf90_def_var(ncid, "vtloc_surf_lvl", nf90_real, &
               (/d_x, d_y, d_t/), vid))
          CALL check(nf90_put_att(ncid, vid, "long_name", &
               "depth, in levels, of the surface localization group"))

          CALL check(nf90_def_var(ncid, "vtloc_surf_depth", nf90_real, &
               (/d_x, d_y, d_t/), vid))
          CALL check(nf90_put_att(ncid, vid, "long_name", &
               "depth, in meters, of the surface localization group"))
          CALL check(nf90_put_att(ncid, vid, "units", "meters"))

          CALL check(nf90_enddef(ncid))
       END IF

       ! write the values
       CALL letkf_mpi_ij2grd(pe_root, self%diag_vtloc_surf_lvl, tmp2d)
       IF(pe_isroot) THEN
          CALL check(nf90_inq_varid(ncid, "vtloc_surf_lvl", vid))
          CALL check(nf90_put_var(ncid, vid, tmp2d))
       END IF

       CALL letkf_mpi_ij2grd(pe_root, self%diag_vtloc_surf_dpth, tmp2d)
       IF(pe_isroot) THEN
          CALL check(nf90_inq_varid(ncid, "vtloc_surf_depth", vid))
          CALL check(nf90_put_var(ncid, vid, tmp2d))
       END IF

       IF(pe_isroot) THEN
          CALL check(nf90_close(ncid))
       END IF
    END IF

  CONTAINS

    SUBROUTINE check(status, str)
      INTEGER, INTENT(in) :: status
      CHARACTER(*), OPTIONAL, INTENT(in) :: str

      IF(status /= nf90_noerr) THEN
         IF(PRESENT(str)) THEN
            WRITE (*,*) TRIM(nf90_strerror(status)), ": ", str
         ELSE
            WRITE (*,*) TRIM(nf90_strerror(status))
         END IF
         CALL letkf_mpi_abort("NetCDF error")
      END IF
    END SUBROUTINE check

  END SUBROUTINE loc_ocean_final
  !================================================================================



  !================================================================================
  !> For a given gridpoint, get the desired localization parameters
  !! (horizontal search radius and level/variable localization groups)
  !--------------------------------------------------------------------------------
  SUBROUTINE loc_ocean_groups(self, ij, groups)
    CLASS(letkf_loc_ocean), INTENT(inout) :: self
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
          ! if based on MLD from background temperature delta
          vtloc_lvl = calc_mldlvl_bkgt(self,ij)
       ELSE
          CALL letkf_mpi_abort("no ocean vertical localization other than BKGT implemented yet.")
       END IF

       ! save diagnostics for output later
       IF(self%save_diag) THEN
          self%diag_vtloc_surf_lvl(ij) = vtloc_lvl
          ! TODO, handle 3D vertical coordinates
          self%diag_vtloc_surf_dpth(ij) = vtgrids(1)%vert_nom(vtloc_lvl)
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
  !> Calculates the vertical localization depth based on the background temperature
  !--------------------------------------------------------------------------------
  FUNCTION calc_mldlvl_bkgt(self, ij) RESULT(mldlvl)
    CLASS(letkf_loc_ocean), INTENT(inout) :: self
    INTEGER, INTENT(in)  :: ij
    INTEGER :: mldlvl
    INTEGER :: i
    REAL :: t0, t1

    ! find the depth at which the mixed layer ends
    t0 = state_mean_ij(self%bkg_t_var_idx, ij)
    mldlvl=self%bkg_t_var_len
    DO i = 1, self%bkg_t_var_len-1
       t1 = state_mean_ij(self%bkg_t_var_idx+i, ij)
       ! IF(pe_isroot)   PRINT *, self%bkg_t_var_idx+i, t1
       IF (ABS(t0-t1) > self%bkg_t_delta) THEN
          mldlvl = i
          EXIT
       END IF
    END DO

  END FUNCTION calc_mldlvl_bkgt
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  FUNCTION loc_ocean_maxhz(self, ij) RESULT(dist)
    CLASS(letkf_loc_ocean), INTENT(inout) :: self
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
    CLASS(letkf_loc_ocean), INTENT(inout) :: self
    INTEGER, INTENT(in) :: ij
    TYPE(letkf_localizer_group), INTENT(in) :: group
    TYPE(letkf_observation), INTENT(in) :: obs
    REAL, INTENT(in) :: dist
    REAL :: loc
    INTEGER :: i
    LOGICAL :: surf

    ! TODO, this logic has caused problems with hybrid-godas.... need to clean it up


    ! Is this a surface/satellite observation?
    ! test by comparing observation and platform ids against a user
    ! defined list.
    surf = .TRUE.
    IF (SIZE(self%surf_obs) > 0) THEN
       DO i = 1, SIZE(self%surf_obs)
          IF (obs%obsid == self%surf_obs(i)) EXIT
       END DO
       IF (i > SIZE(self%surf_obs)) surf = .FALSE.
    END IF
    IF (SIZE(self%surf_plats) > 0) THEN
       DO i = 1, SIZE(self%surf_plats)
          IF (obs%platid == self%surf_plats(i)) EXIT
       END DO
       IF (i > SIZE(self%surf_plats)) surf = .FALSE.
    END IF


    ! calculate the horizontal / vertical / temporal localization
    IF(surf) THEN
       ! if a surface/satellite observation
       IF(group%id == LOC_GROUP_SUBML) THEN
          ! don't use surface observations below the mixed layer
          loc = 0.0
       ELSE
          ! otherwise, in the surface mixed layer, use the satellite ob
          loc = letkf_loc_gc(dist, self%hzdist_sat%get_dist(lat_ij(ij)))
          IF(self%tloc_sat > 0) &
               loc = loc * letkf_loc_gc(obs%time, self%tloc_sat)
       END IF

    ELSE
       ! If this is a T/S profile observation, apply no vertical localization
       loc = letkf_loc_gc(dist, self%hzdist_prof%get_dist( lat_ij(ij)))
       IF(self%tloc_prof > 0) &
            loc = loc * letkf_loc_gc(obs%time, self%tloc_prof)
    ! ELSE
    !    PRINT *, "PLATID: ", obs%platid
    !    CALL letkf_mpi_abort("unhandled platid in loc_ocean_localize")

    ENDIF

  END FUNCTION loc_ocean_localize
  !================================================================================


END MODULE letkf_loc_ocean_mod
