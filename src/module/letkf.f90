module letkf
  !! main entry point for the LETKF library

  use timing
  use letkf_mpi
  use letkf_obs
  use letkf_obs_nc
  use letkf_obs_dat
  use letkf_core
  use letkf_loc
  use letkf_state
  use letkf_state_generic

  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: letkf_driver_init
  public :: letkf_driver_run

  
  ! private module variables
  !------------------------------------------------------------
  character(len=1024) :: nml_filename = "namelist.letkf"
  real :: infl_mul = 1.0
  real :: infl_rtps = 0.0
  real :: infl_rtpp = 0.0
  real :: loc_hz(2) = (/-1.0,-1.0/)
  real :: loc_prune = 0.0
  real, allocatable :: diag_count_ij(:,:)
  integer :: t_total,  t_init



  
contains




  !================================================================================
  !================================================================================
  subroutine letkf_driver_init()
    !! Initialize the LETKF module. This must be called before anything else.
    class(obsio),   pointer :: obsio_ptr
    class(stateio), pointer :: stateio_ptr

    integer :: timer


    ! initialize the mpi backend
    ! (mostly just calls mpi_init)
    call letkf_mpi_preinit()
    call timing_init(mpi_comm_letkf, pe_root)

    ! Initialize timers
    t_total = timer_init("Total", TIMER_SYNC)
    t_init  = timer_init("(init) ", TIMER_SYNC)
    call timer_start(t_total)
    call timer_start(t_init)   

    if (pe_isroot) then
       print "(A)", "============================================================"
       print "(A)", " Universal Multi-Domain Local Ensemble Transform Kalman Filter"
       print "(A)", " (UMD-LETKF)"
       print "(A)", " version 0.1.0"
       print "(A)", " Travis Sluka (tsluka@umd.edu, travis.sluka@noaa.gov)"
       print "(A)", "============================================================"
       print "(A)", ""
    end if

    ! initialize the rest of mpi
    ! (determines the processor distribution for I/O
    call letkf_mpi_init(nml_filename)

    ! setup the default I/O classes, user can add their own after calling 
    ! letkf_driver_init
    allocate(obsio_dat :: obsio_ptr)
    call letkf_obs_register(obsio_ptr)

    allocate(obsio_nc :: obsio_ptr)
    call letkf_obs_register(obsio_ptr)

    allocate(stateio_generic :: stateio_ptr)
    call letkf_state_register(stateio_ptr)

    ! initialize other modules
    call letkf_obs_init(nml_filename, "obsdef.cfg", "platdef.cfg")

    timer = timer_init("  state_init", TIMER_SYNC)
    call timer_start(timer)
    call letkf_state_init(nml_filename)
    call timer_stop(timer)
    call letkf_core_init(mem)

    call letkf_mpi_barrier(.true.)

  end subroutine letkf_driver_init
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_driver_run()
    real, allocatable :: wrk1(:)
    real, allocatable :: wrk3(:,:,:)
    integer :: t_letkf, t_output, timer
    integer :: unit
    integer :: i

    namelist /letkf_inflation/ infl_mul, infl_rtps, infl_rtpp
    namelist /letkf_localization/ loc_hz, loc_prune

    if(pe_isroot) then
       print *, new_line('a'),&
            new_line('a'), '============================================================',&
            new_line('a'), ' letkf_driver_run()',&
            new_line('a'), '============================================================'
    end if

    ! read in main section of the  namelist
    open(newunit=unit, file=nml_filename)
    read(unit, nml=letkf_inflation)
    rewind(unit)
    read(unit, nml=letkf_localization)
    close(unit)
    if (loc_hz(2) <= 0) loc_hz(2) = loc_hz(1)
    if (pe_isroot) then
       print letkf_inflation
       print *, ""
       print letkf_localization
    end if

    !make sure inflation parameters are correct
    !TODO, these will be move to a separate module that deals with localization
    if(pe_isroot) then
       if(loc_hz(1) <=0 .or. loc_hz(2) <= 0) then
          print *, "ERROR: illegal values for loc_hz. ", loc_hz
          stop 1
       end if
       if(infl_rtps > 1.0 .or. infl_rtps < 0.0) then
          print *, "ERROR: illegal value for infl_rtps ",&
               "(should be < 1.0 and >0.0). Value given: ", infl_rtps
          stop 1
       end if
       if(infl_rtpp > 1.0 .or. infl_rtpp < 0.0) then
          print *, "ERROR: illegal value for infl_rtpp ",&
               "(should be < 1.0 and >0.0). Value given: ", infl_rtpp
          stop 1
       end if
       if(infl_rtpp /= 0.0 .and. infl_rtps /= 0.0) then
          print *, "ERROR: cannot have both RTPS and RTPP enabled at the same time, ",&
               "check infl_rtpp and infl_rtps"
          stop 1
       end if
       if (infl_mul < 1.0) then
          print *, "WARNING, infl_mul is < 1.0. Are you sure this is what is intended??? ",&
               "Value given: ",infl_mul
       end if
       if (infl_mul <= 0.0) then
          print *, "ERROR: illegal value for infl_mul: ",infl_mul
       end if
    end if


    ! read observations
    timer = timer_init("  obs_read", TIMER_SYNC)
    call timer_start(timer)
    call letkf_obs_read()
    call timer_stop(timer)


    ! read background state
    timer = timer_init("  state_read", TIMER_SYNC)
    call timer_start(timer)
    call letkf_state_read()
    call timer_stop(timer)
    call timer_stop(t_init)


    ! ensure output directory is setup right
    call system('mkdir -p OUTPUT')


    !------------------------------------------------------------

    ! run LETKF core
    if(pe_isroot) then
       print *, new_line('a'),&
            new_line('a'), '============================================================',&
            new_line('a'), ' Running LETKF core',&
            new_line('a'), '============================================================'
       print *, "Beginning core solver..."
    end if

    t_letkf = timer_init("(solve)", TIMER_SYNC)
    allocate(diag_count_ij(grid_ns,ij_count))
    diag_count_ij = 0.0
    call timer_start(t_letkf)
    call letkf_do_letkf()
    call timer_stop(t_letkf)

    call letkf_mpi_barrier()
    if(pe_isroot) print *, "LETKF solver completed."

    !------------------------------------------------------------

    ! begin output
    t_output = timer_init("(output)", TIMER_SYNC)
    call timer_start(t_output)

    ! save the background / analysis
    call letkf_state_write()

    ! write out other diagnostics
    if(pe_isroot) then
       print *, ""
       print *, "Writing diagnostics files..."
    end if
    if (pe_isroot) allocate(wrk3(grid_nx, grid_ny, grid_ns))
    allocate(wrk1(ij_count))
    do i=1,grid_ns
       wrk1=diag_count_ij(i,:)
       call letkf_mpi_ij2grd(wrk1, wrk3(:,:,i))
    end do
    if(pe_isroot) then
       print '(A,I5,3A)', " PROC ",pe_rank, " is WRITING file: ",&
            'OUTPUT/diag_count', '.nc'
       call stateio_class%write('OUTPUT/diag_count', wrk3)
    end if
    if(pe_isroot) deallocate(wrk3)
    deallocate(wrk1)

    ! all done, cleanup
    call timer_stop(t_output)
    call timer_stop(t_total)
    call timer_print()
    call letkf_mpi_final()

  end subroutine letkf_driver_run
  !================================================================================




  !================================================================================
  !================================================================================
  subroutine letkf_do_letkf()
    !TODO move this to another module
    integer :: ij
    integer, parameter :: maxpt = 100000
    integer :: rpoints(maxpt)
    real :: rdistance(maxpt)
    integer :: rnum, ob_cnt, ob_cnt_lg

    integer :: obidx(maxpt)
    real :: hdxb(mem,maxpt), rdiag(maxpt), dep(maxpt), rdist(maxpt)
    real :: lhdxb(mem,maxpt), lrdiag(maxpt), lrloc(maxpt), ldep(maxpt)
    real :: trans(mem,mem)

    integer :: timer1, timer2, n, timer3, timerloc
    integer :: i, lg, slab
    integer, parameter :: progress_bar_size = 45
    real :: rloc, r
    real :: loc_hz_max
    real :: loc_hz_ij
    real, parameter :: pi = 4.0*atan(1.0)

    type(letkf_loc_group), allocatable :: loc_groups(:)

    timer1   = timer_init("  obs_search")
    timerloc = timer_init("  obs_loc")
    timer2   = timer_init("  core_solve")
    timer3   = timer_init("  core_trans")

    ana_ij = bkg_ij
    
    ! setup progress bar
    if(pe_isroot) then
       write(*,*) ""
       write(*,"(A)", advance='no')      "           ┌"
       do i=1,progress_bar_size
          write(*,"(A)", advance='no')   "─"
       end do
       write(*,"(A)") "┐"
       write(*,"(A)", advance='no')      " progress: │"
    end if


    ! TODO, only call loc_groups once if the localization groups are constant with respect
    ! to the grid position
!    loc_groups = letkf_loc_getgroups(0)


    ! perform analysis at each grid point
    ! ------------------------------
    ij_loop: do ij=1,ij_count

       ! print out progress 
       if(pe_isroot) then
          if ( mod(ij,ij_count /progress_bar_size) == 0) then
             write(*,"(A)", advance='no') "█"
          end if
          if (ij == ij_count) then
             write(*,"(A)") "│"
             write(*,*) ""
          end if
       end if


       !TODO, don't include these gridpoints at all in the mpi ij allocation if the point is masked out
       if(mask_ij(ij)) cycle

       !TODO, use a precomputed cos(lat) grid
       !TODO move max hz distance into the localization module
       ! to allow for user-defined localization methods
       loc_hz_ij = loc_hz(2) + (loc_hz(1)-loc_hz(2))*cos(lat_ij(ij) * pi/180.0)
       loc_hz_max = loc_hz_ij * sqrt(40.0/3.0)

       ! search for all observations in a given radius of this gridpoint
       !TODO, do the KD tree call once to get all possible obs to be used by a grid BLOCK
       !TODO allow for multiple kd-trees (observation groups) each with a different
       ! loc_hz_max (e.g probably want ocean and atm obs separate, because ocean obs 
       ! use a smaller radius than atm obs)
       call timer_start(timer1)
       call letkf_obs_get(lat_ij(ij), lon_ij(ij), loc_hz_max, rpoints, rdistance, rnum)

       ! if there are observations found, process them
       ob_cnt = 0
       do i=1,rnum
          n = rpoints(i)

          !TODO, do this earlier in program
          ! get rid of obs with bad QC values
          if (obs_qc(n) /= 0) cycle

          ! use this observation
          ob_cnt = ob_cnt + 1

          obidx(ob_cnt) = n
          hdxb(:,ob_cnt) = obs_ohx(:,n)
          rdiag(ob_cnt)  = obs_list(n)%err
          rdist(ob_cnt) = rdistance(i)
          dep(ob_cnt) = obs_list(n)%val - obs_ohx_mean(n)
          diag_count_ij(:,ij) = diag_count_ij(:,ij)+1!loc_h
       end do
       call timer_stop(timer1)


       ! perform the core letkf equations and apply transformation matrix
       ! to the background ensemble for EACH localization group
       ! ------------------------------------------------------------
       loc_groups = letkf_loc_getgroups(ij)
       lg_loop: do lg=1,size(loc_groups)         

          ! determine observation localization
          ! if an observation has localization of 0, move it to the end of the array
          ! so that it can be ignored by the letkf_core.
          ! TODO, see if its faster to instead move and index pointing to hdxb,rdiag,dep
          !  and copy them into a new array
          ! TODO, this needs to be faster
          !  precalculate bounds for the observation?
          call timer_start(timerloc)
          ob_cnt_lg = 0
          loc_loop: do i=1,ob_cnt
             rloc = letkf_loc_localize(ij, lg, rdist(i), obidx(i), loc_hz_ij)
             ! if localization is 0, ignore this
             if (rloc > 0.0) then
                ob_cnt_lg = ob_cnt_lg + 1
                lrloc(ob_cnt_lg) = rloc
                lhdxb(:,ob_cnt_lg) = hdxb(:,i)
                lrdiag(ob_cnt_lg) = rdiag(i)
                ldep(ob_cnt_lg) = dep(i)
             end if
          end do loc_loop
          call timer_stop(timerloc)

          !if pruning obs based on the ratio of rloc to max(rloc)
          if(loc_prune > 0.0 ) then
             r = maxval(lrloc(:ob_cnt_lg)) * loc_prune
             do i=ob_cnt_lg,1,-1
                if(lrloc(i)<r) then
                   if(i<ob_cnt_lg) then
                      lrloc(i)=lrloc(ob_cnt_lg)
                      lhdxb(:,i) = lhdxb(:,ob_cnt_lg)
                      lrdiag(i) = lrdiag(ob_cnt_lg)
                      ldep(i) = ldep(ob_cnt_lg)
                   end if
                   ob_cnt_lg = ob_cnt_lg -1
                end if
             end do             
          end if

          ! if there are still good quality observations to assimilate, do so
          if (ob_cnt_lg > 0) then
             
             ! main LETKF equations
             call timer_start(timer2)
             call letkf_core_solve(ob_cnt_lg,&
                  lhdxb(:,:ob_cnt_lg), lrdiag(:ob_cnt_lg), lrloc(:ob_cnt_lg),&
                  ldep(:ob_cnt_lg), 1.0e0, trans)
             call timer_stop(timer2)

             ! calculate the ensemble increments by applying the trans matrix
             !TODO, if not doing vertical localization, do all grid_ns layers at once
             ! .. if it turns out to be any faster
             call timer_start(timer3)
             do i=1,loc_groups(lg)%count
                slab = loc_groups(lg)%slab(i)
                call sgemm('n','n', 1, mem, mem, 1.0e0, bkg_ij(:,slab,ij),&
                     1, trans, mem, 0.0e0, ana_ij(:,slab,ij), 1)
             end do
             call timer_stop(timer3)
          end if
       end do lg_loop

          
       ! recenter analysis perturbations on new analysis mean
       do i=1,grid_ns
          ana_ij(:,i,ij) = ana_ij(:,i,ij) + bkg_mean_ij(i,ij)
          ana_mean_ij(i,ij) = sum(ana_ij(:,i,ij),1)/mem
          ana_ij(:,i,ij) = ana_ij(:,i,ij) - ana_mean_ij(i,ij)
       end do
       do i = 1, grid_ns
          ana_sprd_ij(i,ij) = sqrt(dot_product(ana_ij(:,i,ij),ana_ij(:,i,ij)) /mem)
       end do


       ! inflation
       ! ------------------------------
       ! RTPS (relaxation to prior spread)
       if(infl_rtps > 0) then
          do i =1,grid_ns
             if (ana_sprd_ij(i,ij) > 0) then
                ana_ij(:,i,ij) = ana_ij(:,i,ij) * &
                     (infl_rtps * (bkg_sprd_ij(i,ij)-ana_sprd_ij(i,ij))/ana_sprd_ij(i,ij) + 1)
             end if
          end do
       end if

       ! RTPP (relaxation to prior perturbations
       if(infl_rtpp > 0) then
          do i=1,grid_ns
             ana_ij(:,i,ij) = bkg_ij(:,i,ij) * infl_rtpp + (1.0-infl_rtpp)*ana_ij(:,i,ij)
          end do
       end if

       ! constant multiplicative
       if(infl_mul /= 1.0) then
          ana_ij(:,:,ij) = ana_ij(:,:,ij) * infl_mul
       end if

       ! ------------------------------

       !add mean to perturbations
       do i=1,grid_ns
          ana_ij(:,i,ij) = ana_ij(:,i,ij) + ana_mean_ij(i,ij)
       end do


       ! recalculate analysis spread
       do i = 1, grid_ns
          ana_sprd_ij(i,ij) = sqrt(&
               dot_product(ana_ij(:,i,ij)-ana_mean_ij(i,ij),ana_ij(:,i,ij)-ana_mean_ij(i,ij)) /mem)
       end do
       
    end do ij_loop

  end subroutine letkf_do_letkf
  !================================================================================

end module letkf
