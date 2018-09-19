!================================================================================
!> does the majority of the real "LETKF" work.
!--------------------------------------------------------------------------------
MODULE letkf_solver
  USE timing
  USE letkf_config
  USE letkf_diag
  USE letkf_mpi
  USE letkf_loc
  USE letkf_state
  USE letkf_obs
  USE letkf_core

  IMPLICIT NONE
  PRIVATE


  !================================================================================
  !================================================================================
  ! public module components
  !================================================================================
  !================================================================================

  PUBLIC :: letkf_solver_init
  PUBLIC :: letkf_solver_run


  !================================================================================
  !================================================================================
  ! private module components
  !================================================================================
  !================================================================================

  REAL, PARAMETER :: stdev2max = SQRT(40.0/3.0)

  ! inflation values
  !--------------------------------------------------------------------------------
  REAL :: infl_rtps
  REAL :: infl_rtpp
  REAL :: infl_mul



CONTAINS



  !================================================================================
  !> initialze the module by loading in configuration settings
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_solver_init(config)
    TYPE(configuration), INTENT(in) :: config

    TYPE(configuration) :: infl_config

    IF (pe_isroot) THEN
       PRINT *, ""
       PRINT *, "======================================================================"
       PRINT *, " letkf_solver_init() : "
       PRINT *, "======================================================================"
    END IF


    ! get the inflation parameters
    CALL config%get("inflation", infl_config)
    CALL infl_config%get("mul", infl_mul, default=1.0)
    CALL infl_config%get("rtps", infl_rtps, default=0.0)
    CALL infl_config%get("rtpp", infl_rtpp, default=0.0)

    IF(pe_isroot) THEN
       PRINT *, "Inflation: "
       PRINT *, "solver.inflation.rtps=",infl_rtps
       PRINT *, "solver.inflation.rtpp=",infl_rtpp
       PRINT *, "solver.inflation.mul=",infl_mul
    END IF

    ! make sure the inflation parameters are correct
    IF(pe_isroot) THEN
       IF(infl_rtps > 1.0 .OR. infl_rtps < 0.0) &
            CALL letkf_mpi_abort("illegal value for infl_rtps "//&
            "(should be <1.0 and >0.0).")
       IF(infl_rtpp > 1.0 .OR. infl_rtpp < 0.0) &
            CALL letkf_mpi_abort("illegal value for infl_rtpp "//&
            "(should be <1.0 and >0.0).")
       IF(infl_rtpp /= 0.0 .AND. infl_rtps /= 0.0) &
            CALL letkf_mpi_abort("cannot have both RTPS and RTPP enabled at the same time. ")
       IF(infl_mul < 1.0) &
            CALL letkf_mpi_abort("infl_mul must be >=1.0.")
    END IF

    CALL letkf_core_init(ens_size)

    ! initialize the diagnostic fields
    CALL letkf_diag_reg("maxhz",&
         desc="max observation search radius")
    CALL letkf_diag_reg("obs_count",&
         desc="number of observations returned from KD tree")
    CALL letkf_diag_reg("obs_count_loc", &
         desc="sum of localization values for obs within maxhz of gridpoint")

  END SUBROUTINE letkf_solver_init
  !================================================================================



  !================================================================================
  !>
  !--------------------------------------------------------------------------------
  SUBROUTINE letkf_solver_run()
    ! TODO, a bit messy, get rid of the duplicate "obs_lg" variables
    REAL :: maxhz
    INTEGER :: i, j, k, idx, ij, lg, slab
    REAL :: obs_cnt_loc
    REAL :: r
    REAL :: max_search_dist

    REAL :: bkg_ij(ens_size, grid_ns)
    REAL :: bkg_mean_ij(grid_ns)
    REAL :: bkg_sprd_ij(grid_ns)

    REAL :: trans(ens_size, ens_size)

    ! TODO, remove the hardcoding for this
    INTEGER,PARAMETER :: max_obs = 100000

    INTEGER :: obs_ij_cnt
    INTEGER :: obs_ij_idx(max_obs)
    REAL    :: obs_ij_dist(max_obs)
    REAL    :: obs_ij_hdxb(ens_size, max_obs)
    REAL    :: obs_ij_rdiag(max_obs)
    REAL    :: obs_ij_dep(max_obs)

    INTEGER :: obs_lg_cnt
    REAL    :: obs_lg_dist(max_obs)
    REAL    :: obs_lg_hdxb(ens_size, max_obs)
    REAL    :: obs_lg_rdiag(max_obs)
    REAL    :: obs_lg_dep(max_obs)
    REAL    :: obs_lg_rloc(max_obs)

    TYPE(letkf_localizer_group), ALLOCATABLE :: loc_groups(:)


    CALL timing_start('solver', TIMER_SYNC)
    IF (pe_isroot) THEN
       PRINT '(//,A)', ""
       PRINT *, "======================================================================"
       PRINT *, "letkf_solver_run() : Core LETKF routine"
       PRINT *, "======================================================================"
    END IF

    ! TODO make diag_obs_cnt_loc multidimensional, depending on the max
    ! number of localization groups expected
    !    ALLOCATE(diag_obs_cnt_loc(ij_count))

    ! perform analysis at each grid point
    ij_loop: DO ij=1,ij_count

       ! skip the gridpoint if it is masked
       ! TODO, don't even include these points in the array_temp
       IF(mask_ij(ij)) CYCLE

       ! get the max horizontal search distance
       maxhz = localizer_class%maxhz(ij)
       CALL letkf_diag_set("maxhz", ij, maxhz)
       max_search_dist = maxhz * stdev2max

       ! search for all observations within a given radius of this point
       ! using a fast O(n log n) kd-tree
       CALL timing_start("obs_search")
       CALL letkf_obs_get(lat_ij(ij), lon_ij(ij), max_search_dist, &
            obs_ij_idx, obs_ij_dist, obs_ij_cnt)
       CALL timing_stop("obs_search")
       CALL letkf_diag_set("obs_count", ij, REAL(obs_ij_cnt))


       ! If there are observations found, prepare them
       DO i=1,obs_ij_cnt
          idx = obs_ij_idx(i)
          obs_ij_hdxb(:,i) = obs_hx(:,idx)
          obs_ij_rdiag(i)  = obs_def(idx)%err
          obs_ij_dep(i)    = obs_def(idx)%val - obs_hx_mean(idx)
       END DO

       ! get the number of localization groups
       CALL localizer_class%groups(ij, loc_groups)

       ! save the current state as the "background", we'll need it again later
       bkg_ij=state_ij(:,:,ij)
       bkg_mean_ij=state_mean_ij(:,ij)
       bkg_sprd_ij=state_sprd_ij(:,ij)

       ! for each localization group
       lg_loop: DO lg=1, SIZE(loc_groups)
          ! TEMPORARY
          IF (lg > 1) CALL letkf_mpi_abort(">1 loc group not yet supported")

          ! localize the obsservations, keeping only the obs with localization > 0.0
          obs_lg_cnt = 0
          obs_cnt_loc = 0.0
          obs_loop: DO i=1,obs_ij_cnt
             idx=obs_ij_idx(i)
             r = localizer_class%localize(ij, loc_groups(lg), &
                  obs_def(idx), obs_ij_dist(i))
             obs_cnt_loc = obs_cnt_loc + r

             IF (r > 0.0) THEN
                ! this observation will be kept
                ! TODO, check the QC flag too? or have I implemented
                ! removal of bad obs higher up in the code somewhere
                obs_lg_cnt = obs_lg_cnt + 1
                obs_lg_rloc(obs_lg_cnt) = r
                obs_lg_dep(obs_lg_cnt) = obs_ij_dep(i)
                obs_lg_dist(obs_lg_cnt) = obs_ij_dist(i)
                obs_lg_hdxb(:,obs_lg_cnt) = obs_ij_hdxb(:,i)
                obs_lg_rdiag(obs_lg_cnt) = obs_ij_rdiag(i)
             END IF
          END DO obs_loop
          CALL letkf_diag_set("obs_count_loc", ij, obs_cnt_loc)

          ! if there are still good quality obs to assimilate, do so
          IF (obs_lg_cnt > 0) THEN

             ! TODO, init the timer outside the loop (in case a PE has no instances
             ! in this loop)
             CALL timing_start("letkf_core")
             CALL letkf_core_solve(obs_lg_cnt, obs_lg_hdxb(:,:obs_lg_cnt),&
                  obs_lg_rdiag(:obs_lg_cnt), obs_lg_rloc(:obs_lg_cnt), &
                  obs_lg_dep(:obs_lg_cnt), 1.0e0, trans)
             CALL timing_stop("letkf_core")

             ! apply the trans matrix
             ! TODO, if slabs are contiguous, can we just do a single sgemm call??
             CALL timing_start("letkf_trans")
             DO i=1, SIZE(loc_groups(lg)%slab)
                slab = loc_groups(lg)%slab(i)
                CALL sgemm('n', 'n', 1, ens_size, ens_size, 1.0e0, bkg_ij(:,slab), &
                     1, trans, ens_size, 0.0e0, state_ij(:,slab,ij), 1)
             END DO
             CALL timing_stop("letkf_trans")
          END IF
       END DO lg_loop


       ! TODO, the following should be rearranged for improved performance
       ! things are being unecessarily recalculated repeatedly

       ! recenter the analysis perturbations on the new analysis mean
       DO i=1, grid_ns
          state_ij(:,i,ij) = state_ij(:,i,ij) + bkg_mean_ij(i)
          state_mean_ij(i,ij) = SUM(state_ij(:,i,ij),1)/ens_size
          state_ij(:,i,ij) = state_ij(:,i,ij) - state_mean_ij(i,ij)
       END DO


       ! clip analysis increment based on both
       ! 1) max analysis increment 2) analysis bounds
       DO i=1, SIZE(statevars)
          DO j = statevars(i)%grid_s_idx, statevars(i)%grid_s_idx+statevars(i)%levels-1
             r = state_mean_ij(j,ij)-bkg_mean_ij(j)
             r = MAX(r, -statevars(i)%ana_inc_max)
             r = MIN(r,  statevars(i)%ana_inc_max)
             r = MAX(r, statevars(i)%ana_bounds(1)-bkg_mean_ij(j))
             r = MIN(r, statevars(i)%ana_bounds(2)-bkg_mean_ij(j))
             state_mean_ij(j,ij) = r + bkg_mean_ij(j)
          END DO
          ! TODO set a flag to alert user that AI was clipped
       END DO


       ! apply RTPS inflation, if enabled
       IF(infl_rtps > 0 ) THEN
          DO i=1,grid_ns
             ! calculate spread
             state_sprd_ij(i,ij) = SQRT(&
                  dot_PRODUCT(state_ij(:,i,ij), state_ij(:,i,ij))/(ens_size-1))

             ! no spread, so skip
             IF (state_sprd_ij(i,ij) <= 0) CYCLE

             ! increase the spread
             state_ij(:,i,ij) = state_ij(:,i,ij) * &
                  (infl_rtps * (bkg_sprd_ij(i) - state_sprd_ij(i,ij))/state_sprd_ij(i,ij)+1)
          END DO
       END IF


       ! apply RTPP inflation, if enabled
       IF(infl_rtpp > 0.0 ) THEN
          DO i=1,grid_ns
             state_ij(:,i,ij) = infl_rtpp*bkg_ij(:,i) + (1.0-infl_rtpp)*state_ij(:,i,ij)
          END DO
       END IF


       ! constant multiplicative inflation, if enabled
       IF(infl_mul /= 1.0) THEN
          state_ij(:,:,ij) = state_ij(:,:,ij) * infl_mul
       END IF


       ! TODO apply additive inflation

       ! TODO apply adaptive inflation

       ! bounds check on individual analysis ensemble members
       DO i=1, SIZE(statevars)
          DO j = statevars(i)%grid_s_idx, statevars(i)%grid_s_idx+statevars(i)%levels-1
             DO k = 1, ens_size
                r = state_ij(k,j,ij)
                r = MAX(r, statevars(i)%ana_bounds(1)-state_mean_ij(j,ij))
                r = MIN(r, statevars(i)%ana_bounds(2)-state_mean_ij(j,ij))
                state_ij(k,j,ij) = r
             END DO
             ! adjusting the mean afterward correctly
             r = SUM(state_ij(:,j,ij),1)/ens_size
             state_ij(:,j,ij) = state_ij(:,j,ij) -r
             state_mean_ij(j,ij) = state_mean_ij(j,ij) +r
             !    ! TODO set a flag to alert user that AI was clipped
          END DO
       END DO


       ! calc final analysis members, spread
       DO i = 1, grid_ns
          state_sprd_ij(i,ij) = SQRT(&
               dot_PRODUCT(state_ij(:,i,ij), state_ij(:,i,ij))/(ens_size-1))
          state_ij(:,i,ij) = state_ij(:,i,ij) + state_mean_ij(i,ij)
       END DO

       ! cleanup
       DEALLOCATE(loc_groups)

    END DO ij_loop

    CALL timing_stop('solver')

  END SUBROUTINE letkf_solver_run
  !================================================================================


END MODULE letkf_solver
