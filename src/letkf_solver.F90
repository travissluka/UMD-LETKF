MODULE letkf_solver
  USE timing
  USE letkf_mpi
  USE letkf_loc
  USE letkf_state
  USE letkf_obs
  use letkf_core

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: letkf_solver_init
  PUBLIC :: letkf_solver_run

  ! diagnostics, in ij space
  REAL, PUBLIC, PROTECTED, ALLOCATABLE :: diag_maxhz(:)
  REAL, PUBLIC, PROTECTED, ALLOCATABLE :: diag_obs_cnt(:)
  REAL, PUBLIC, PROTECTED, ALLOCATABLE :: diag_obs_cnt_loc(:)


  REAL, PARAMETER :: stdev2max = SQRT(40.0/3.0)


CONTAINS


  
  SUBROUTINE letkf_solver_init(nml_filename)
    character(:), allocatable, intent(in) :: nml_filename
    !TODO read in our section of the namelist
    
    !TODO double check some of the parameters

    call letkf_core_init(ens_size)
  END SUBROUTINE letkf_solver_init



  !-----------------------------------------------------------------------------
  !>
  SUBROUTINE letkf_solver_run()
    ! TODO, a bit messy, get rid of the duplicate "obs_lg" variables
    INTEGER :: i, idx, ij, lg, slab
    real :: r
    REAL :: max_search_dist

    REAL :: bkg_ij(ens_size, grid_ns)
    REAL :: bkg_mean_ij(grid_ns)
    REAL :: bkg_sprd_ij(grid_ns)

    real :: trans(ens_size, ens_size)

    INTEGER, PARAMETER :: max_obs = 100

    INTEGER :: obs_ij_cnt
    INTEGER :: obs_ij_idx(max_obs)
    REAL    :: obs_ij_dist(max_obs)
    real    :: obs_ij_hdxb(ens_size, max_obs)
    real    :: obs_ij_rdiag(max_obs)
    real    :: obs_ij_dep(max_obs)

    INTEGER :: obs_lg_cnt
    REAL    :: obs_lg_dist(max_obs)
    real    :: obs_lg_hdxb(ens_size, max_obs)
    real    :: obs_lg_rdiag(max_obs)
    real    :: obs_lg_dep(max_obs)
    real    :: obs_lg_rloc(max_obs)

    class(letkf_localizer_group), allocatable :: loc_groups(:)

    
    
    CALL timing_start('solver', TIMER_SYNC)
    IF (pe_isroot) THEN
       print '(//,A)', ""
       PRINT *, "=========================================================================================="
       PRINT *, "letkf_solver_run() : Core LETKF routine"
       PRINT *, "=========================================================================================="
    END IF

    ! TODO make diag_obs_cnt_loc multidimensional, depending on the max
    ! number of localization groups expected
    allocate(diag_maxhz(ij_count))
    ALLOCATE(diag_obs_cnt(ij_count))    
    ALLOCATE(diag_obs_cnt_loc(ij_count))
    diag_maxhz = 0.0
    diag_obs_cnt = 0
    diag_obs_cnt_loc = 0.0

    ! perform analysis at each grid point
    ij_loop: DO ij=1,ij_count

       ! skip the gridpoint if it is masked
       ! TODO, don't even include these points in the array_temp
       IF(mask_ij(ij)) CYCLE

       ! get the max horizontal search distance
       diag_maxhz(ij) = localizer_class%maxhz(ij)
       max_search_dist = diag_maxhz(ij) * stdev2max

       ! search for all observations within a given radius of this point
       ! using a fast O(n log n) kd-tree
       CALL timing_start("obs_search")
       CALL letkf_obs_get(lat_ij(ij), lon_ij(ij), max_search_dist, &
            obs_ij_idx, obs_ij_dist, obs_ij_cnt)
       CALL timing_stop("obs_search")
       diag_obs_cnt(ij) = obs_ij_cnt

       ! If there are observations found, prepare them
       do i=1,obs_ij_cnt
          idx = obs_ij_idx(i)
          obs_ij_hdxb(:,i) = obs_hx(:,idx)
          obs_ij_rdiag(i)  = obs_def(idx)%err
          obs_ij_dep(i)    = obs_def(idx)%val - obs_hx_mean(idx)
       end do
       
       ! get the number of localization groups
       call localizer_class%groups(ij, loc_groups)       

       ! save the current state as the "background", we'll need it again later
       bkg_ij=state_ij(:,:,ij)
       bkg_mean_ij=state_mean_ij(:,ij)
       bkg_sprd_ij=state_sprd_ij(:,ij)

       ! for each localization group
       lg_loop: DO lg=1, SIZE(loc_groups)
          ! TEMPORARY
          if (lg > 1) call letkf_mpi_abort(">1 loc group not yet supported")
          
          ! localize the obsservations, keeping only the obs with localization > 0.0
          obs_lg_cnt = 0
          obs_loop: do i=1,obs_ij_cnt
             idx=obs_ij_idx(i)
             r = localizer_class%localize(ij, loc_groups(lg), &
                  obs_def(idx), obs_ij_dist(i))
             diag_obs_cnt_loc(ij) = diag_obs_cnt_loc(ij) + r

             if (r > 0.0) then
                ! this observation will be kept
                ! TODO, check the QC flag too? or have I implemented
                ! removal of bad obs higher up in the code somewhere
                obs_lg_cnt = obs_lg_cnt + 1
                obs_lg_rloc(obs_lg_cnt) = r
                obs_lg_dep(obs_lg_cnt) = obs_ij_dep(i)
                obs_lg_dist(obs_lg_cnt) = obs_ij_dist(i)
                obs_lg_hdxb(:,obs_lg_cnt) = obs_ij_hdxb(:,i)
                obs_lg_rdiag(obs_lg_cnt) = obs_ij_rdiag(i)
             end if
          end do obs_loop

          ! if there are still good quality obs to assimilate, do so
          IF (obs_lg_cnt > 0) THEN
             
             ! TODO, init the timer outside the loop (in case a PE has no instances
             ! in this loop)
             call timing_start("letkf_core")             
             call letkf_core_solve(obs_lg_cnt, obs_lg_hdxb(:,:obs_lg_cnt),&
                  obs_lg_rdiag(:obs_lg_cnt), obs_lg_rloc(:obs_lg_cnt), &
                  obs_lg_dep(:obs_lg_cnt), 1.0e0, trans)
             call timing_stop("letkf_core")
             
             ! apply the trans matrix
             ! TODO, if slabs are contiguous, can we just do a single sgemm call??
             call timing_start("letkf_trans")
             do i=1, size(loc_groups(lg)%slab)
                slab = loc_groups(lg)%slab(i)
                call sgemm('n', 'n', 1, ens_size, ens_size, 1.0e0, bkg_ij(:,slab), &
                     1, trans, ens_size, 0.0e0, state_ij(:,slab,ij), 1)
             end do
             call timing_stop("letkf_trans")
          END IF
       END DO lg_loop


       ! recenter the analysis perturbations on the new analysis mean
       do i=1, grid_ns
          state_ij(:,i,ij) = state_ij(:,i,ij) + bkg_mean_ij(i)
          state_mean_ij(i,ij) = sum(state_ij(:,i,ij),1)/ens_size
          state_ij(:,i,ij) = state_ij(:,i,ij) - state_mean_ij(i,ij)
       end do

       
       ! TODO apply RTPP inflation, if enabled

       ! TODO apply RTPS inflation, if enabled

       ! TODO constant multiplicative inflation, if enabled

       ! TODO apply additive inflation

       ! TODO apply adaptive inflation
       
       ! calc final analysis members, spread
       do i = 1, grid_ns
          state_sprd_ij(i,ij) = sqrt(&
               dot_product(state_ij(:,i,ij), state_ij(:,i,ij))/ens_size)
          state_ij(:,i,ij) = state_ij(:,i,ij) + state_mean_ij(i,ij)
       end do

    END DO ij_loop

    CALL timing_stop('solver')

  END SUBROUTINE letkf_solver_run

END MODULE letkf_solver
