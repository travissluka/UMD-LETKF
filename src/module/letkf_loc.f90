module letkf_loc
  !! observation localization routines
  use letkf_obs

  implicit none
  private

  ! public module methods
  !------------------------------------------------------------
  public :: loc_gc, loc_gaus
  public :: letkf_loc_getgroups
  public :: letkf_loc_localize

  public :: letkf_loc_group
  type letkf_loc_group
     integer :: count
     integer, allocatable :: slab(:)
  end type letkf_loc_group




contains




  !================================================================================
  !================================================================================
  function letkf_loc_getgroups(ij) result(grps)
    integer, intent(in) :: ij
    type(letkf_loc_group), allocatable :: grps(:)
    integer:: i


    ! TODO, create configurable loc groups (make sure for this
    !  that a check is done to ensure each slab is placed in
    !  exactly 1 lg)
    ! TODO, do the allocation beforehand and save it if we
    !  are not doing dynamic loc groups, should be faster

    ! allocate(grps(40))
    ! do i=1,40
    !    grps(i)%count=2
    !    allocate(grps(i)%slab(2))
    !    grps(i)%slab(1) = i
    !    grps(i)%slab(2) = i+40
    ! end do


    allocate(grps(1))
    grps(1)%count=80
    allocate(grps(1)%slab(80))
    do i=1,80
       grps(1)%slab(i) = i
    end do
    
  end function letkf_loc_getgroups
  !================================================================================




  !================================================================================
  !================================================================================
  pure function letkf_loc_localize(grd_ij, grd_lg, ob_dist, ob_idx, tmp) result(val)
    integer, intent(in) :: grd_ij, grd_lg
    real,    intent(in) :: ob_dist
    integer, intent(in) :: ob_idx
    real,    intent(in) :: tmp

    real :: val, grd_depth
    

    ! TODO, this needs to be user definable
    ! TODO, caching of horizontal localization, for instances where
    ! horizontal scale is the same for each lg

    ! grd_depth = -5.0 + grd_lg*10.0  
    ! val = loc_gc(ob_dist, tmp) * &
    !       loc_gc( abs(obs_list(ob_idx)%depth - grd_depth), 10.0)

    val = loc_gc(ob_dist, tmp)


  end function letkf_loc_localize
  !================================================================================


  !================================================================================
  !================================================================================
  pure function loc_gc(z, L)
    !! Gaspari-Cohn localization function
    !! Possibly faster than the Gaussian function, depending on computer architecture.
    !! Similar shape to Gaussian, except it is compact, goes to 0 at (\ 2L sqrt( 0.3) \ )
    real, intent(in) :: z
    real, intent(in) :: L
    !! (\ e^(0.5) \)
    real :: loc_gc
    real :: c
    real :: abs_z, z_c

    c = L / sqrt(0.3)
    abs_z = abs(z)
    z_c = abs_z / c

    if (abs_z >= 2*c) then
       loc_gc = 0.0
    elseif (abs_z < 2*c .and. abs_z > c) then
       loc_gc = &
            (1.0/12.0)*z_c**5 - 0.5*z_c**4 + &
            (5.0/8.0)*z_c**3 + (5.0/3.0)*z_c**2 &
            - 5.0*z_c + 4 - (2.0/3.0)*c/abs_z
    else
       loc_gc = &
            -0.25*z_c**5 + 0.5*z_c**4 + &
            (5.0/8.0)*z_c**3 - (5.0/3.0)*z_c**2 + 1
    end if
  end function loc_gc
  !================================================================================



  !================================================================================
  !================================================================================
  pure function loc_gaus(z, L)
    real, intent(in) :: z, L
    real :: loc_gaus
    loc_gaus = exp( -0.5 *  z*z / L*L)
  end function loc_gaus
  !================================================================================

end module letkf_loc
