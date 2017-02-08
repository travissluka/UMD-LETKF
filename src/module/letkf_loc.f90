module letkf_loc
  !! observation localization routines
  implicit none
  private

  public :: loc_gc, loc_gaus

contains

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


  pure function loc_gaus(z, L)
    real, intent(in) :: z, L
    real :: loc_gaus
    loc_gaus = exp( -0.5 *  z*z / L*L)
  end function loc_gaus
end module letkf_loc
