!> \file daylight.f90  Procedures that deal with daylight for libTheSky


!  Copyright (c) 2002-2020  Marc van der Sluys - marc.vandersluys.nl
!   
!  This file is part of the libTheSky package, 
!  see: http://libthesky.sf.net/
!   
!  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!  
!  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
!  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License along with this code.  If not, see 
!  <http://www.gnu.org/licenses/>.


!  extinction_sun_airmass:       Compute an approximation for the bolometric atmospheric extinction for the Sun with a given AM
!  extinction_sun:               Compute an approximation for the bolometric atmospheric extinction for the Sun with a given alt
!  solar_radiation:              Compute the normal and horizontal beam (direct) solar radiation for a given altitude of the Sun
!  diffuse_radiation_Perez87:    Compute diffuse radiation on an inclined surface using the Perez 1987 model


!***********************************************************************************************************************************
!> \brief Daylight procedures

module TheSky_daylight
  implicit none
  save
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute an approximation for the bolometric atmospheric extinction factor for the Sun with a given air mass
  !!
  !! \param airmass  Airmass for the position of the Sun
  !!
  !! \note  - extinction_fac = 1: no extinction, extinction_fac > 1 extinction.
  !!        - Hence, the flux, corrected for extinction, should be  f' = f / extinction_fac(alt,ele)
  !!        - Fit of the power extinction computed by the NREL SMARTS code (valid for airmass <= 38.2):
  !!          - Gueymard, C, Professional Paper FSEC-PF-270-95, 1995
  !!          - Gueymard, C, Solar Energy, Vol. 71, No. 5, pp. 325-346, 2001
  
  function extinction_sun_airmass(airmass)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: airmass
    integer, parameter :: ncoef = 11
    integer :: iCoef
    real(double) :: extinction_sun_airmass, coefs(ncoef), AMpow
    
    if(airmass.gt.38.2d0) then
       
       extinction_sun_airmass = sqrt(huge(extinction_sun_airmass)) + airmass  ! Very bad, getting worse for higher AM for solvers
       
    else
       
       coefs = [ 9.1619283d-2, 2.6098406d-1,-3.6487512d-2, 6.4036283d-3,-8.1993861d-4, 6.9994043d-5,-3.8980993d-6, 1.3929599d-7, &
            -3.0685834d-9, 3.7844273d-11,-1.9955057d-13]
       
       AMpow = 1.d0                                                               ! AM^0
       extinction_sun_airmass = coefs(1)                                          ! c_1 * AM^0
       do iCoef=2,ncoef
          AMpow = AMpow * airmass                                                 ! AM^(i-1)
          extinction_sun_airmass = extinction_sun_airmass + coefs(iCoef) * AMpow  ! + c_i * AM^(i-1)
       end do
       
       extinction_sun_airmass = exp(extinction_sun_airmass)
       
    end if
    
  end function extinction_sun_airmass
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute an approximation for the bolometric atmospheric extinction factor for the Sun with a given altitude
  !!
  !! \param alt  TRUE altitude of the Sun (radians)
  !!
  !! \note  - extinction_fac = 1: no extinction, extinction_fac > 1 extinction.
  !!        - Hence, the flux, corrected for extinction, should be  f' = f / extinction_fac(alt,ele)
  
  function extinction_sun(alt)
    use SUFR_kinds, only: double
    use TheSky_visibility, only: airmass
    
    implicit none
    real(double), intent(in) :: alt
    real(double) :: extinction_sun
    
    extinction_sun = extinction_sun_airmass(airmass(alt))
    
  end function extinction_sun
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the normal and horizontal beam (direct) solar radiation for a given altitude of the Sun, 
  !!         assuming a cloudless sky (and sea level)
  !!
  !! \param alt  TRUE altitude of the Sun (radians)
  !!
  !! \param beam_norm   Normal beam radiation, perpendicular to the position vector of the Sun (W/m2) (output)
  !! \param beam_horiz  Beam radiation on a horizontal surface (W/m2) (output)
  !!
  !! \see  function extinction_sun()
  
  subroutine solar_radiation(alt,  beam_norm, beam_horiz)
    use SUFR_kinds, only: double
    use SUFR_constants, only: solConst
    
    implicit none
    real(double), intent(in) :: alt
    real(double), intent(out) :: beam_norm
    real(double), intent(out), optional :: beam_horiz
    
    
    beam_norm  = solConst / extinction_sun(alt)                ! Normal radiation
    if(present(beam_horiz)) beam_horiz = beam_norm * sin(alt)  ! Radiation on a horizontal surface
    
  end subroutine solar_radiation
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  A simplified clear-sky model for direct and diffuse insolation on horizontal surfaces (a.k.a. as the Bird model)
  !!
  !! \param   alt    Sun altitude above the horizon (rad)
  !!
  !! \param   Io     Solar 'constant' (W/m^2 - optional, default: 1361.5)
  !! \param   Rsun   Sun distance (AU - optional, default: 1)
  !! \param   Press  Air pressure at the observer's site, corrected for altitude (hPa - optional, default: 1013)
  !!
  !! \param   Uo     Ozone abundance in a vertical column (cm - optional, default: 0.34)
  !! \param   Uw     Percipitable water-vapor abundance in a vertical column (cm - optional, default: 1.42)
  !!
  !! \param   Ta5    Aerosol optical depth from surface in vertical path at 500 nm (optional, default: 0.2661)
  !! \param   Ta3    Aerosol optical depth from surface in vertical path at 380 nm (optional, default: 0.3538)
  !! \param   Ba     Aerosol forward-scattering ratio  (optional, 0.82-0.86, default: 0.84)
  !! \param   K1     Aerosol-absorptance constant (optional, rural: 0.0933, urban: 0.385, default: 0.1)
  !!
  !! \param   Rg     Ground albedo (optional, fraction - default: 0.2)
  !!
  !!
  !! \param  Itot   Total insolation on a horizontal surface (W/m^2 - optional) (output)
  !! \param  Idir   Direct (beam) insolation on a horizontal surface (W/m^2 - optional) (output)
  !! \param  Idif   Diffuse insolation on a horizontal surface (W/m^2 - optional) (output)
  !! \param  Igr    Ground-reflection insolation from a horizontal surface (W/m^2 - optional) (output)
  !!
  !!
  !! \see  Bird & Hulstrom, A simplified clear-sky model for direct and diffuse insolation on horizontal surfaces, SERI/TR-642-761 (1981)
  !!
  !! \note The value of Taa does not agree with tabulated values from the paper, and hence neither do dependent values (except for AM~1).
  !!       When I substitute their values for Taa, everything matches perfectly.  Error in their formula, or (hopefully!) in their table?
  
  subroutine clearsky_bird(alt, Io,Rsun,Press, Uo,Uw, Ta5,Ta3,Ba,K1, Rg, Itot,Idir,Idif,Igr)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pio2, r2d  !, solConst
    
    implicit none
    real(double), intent(in) :: alt
    real(double), intent(in), optional :: Io,Rsun,Press, Uo,Uw, Ta5,Ta3,Ba,K1, Rg
    real(double), intent(out), optional :: Itot, Idir, Idif, Igr
    
    real(double) :: Z,cosZ, AM,AMp, Tr, Xo,To, Tum, Xw,Tw, Tau,Ta,Taa,Tas, Rs, tmpVar
    real(double) :: RsunL,IoL,PressL, UoL,UwL, Ta5L,Ta3L,BaL,K1L, RgL,   ItotL,IdirL,IdifL
    
    ! Assign optional variables to local variables:
    RsunL = 1.d0      ! Default Sun distance: 1 AU
    if(present(Rsun)) RsunL = Rsun
    IoL = 1353.d0     ! solConst    ! Default solar "constant": 1353 (1361.5) W/m^2
    if(present(Io)) IoL = Io
    PressL = 1013.d0  ! Default air pressure: 1013 hPa
    if(present(Press)) PressL = Press
    
    UoL = 0.34d0      ! Default ozone abundance: 0.34
    if(present(Uo)) UoL = Uo
    UwL = 1.42d0      ! Default water-vapor abundance: 1.42
    if(present(Uw)) UwL = Uw
    
    Ta5L = 0.2661d0   ! Default aerosol vertical optical depth at 500nm: 0.2661 cm
    if(present(Ta5)) Ta5L = Ta5
    Ta3L = 0.3538d0   ! Default aerosol vertical optical depth at 500nm: 0.3538 cm
    if(present(Ta3)) Ta3L = Ta3
    BaL = 0.84d0      ! Default aerosol forward-scattering ratio: 0.84
    if(present(Ba))  BaL = Ba
    K1L = 0.1d0       ! Default aerosol-absorptance constant  (rural: 0.0933, urban: 0.385, adviced: 0.1)
    if(present(K1))  K1L = K1
    
    RgL = 0.2d0       ! Default ground albedo: 0.2
    if(present(Rg)) RgL = Rg
    
    
    
    Z = pio2 - alt  ! Solar zenith angle
    cosZ = cos(Z)   ! Save a few CPU cycles
    
    
    ! Relative air mass for the solar vector:
    AM  = 1.d0/(cosZ + 0.15d0 * (93.885-Z*r2d)**(-1.25d0))  ! Air mass
    AMp = AM * PressL / 1013.d0                             ! Pressure-corrected air mass
    
    ! TRANSMISSION EQUATIONS:
    ! Rayleigh scattering:
    Tr = exp( -0.0903d0 * AMp**0.84d0 * (1.d0 + AMp - AMp**1.01d0) )
    
    ! Ozone:
    Xo = UoL*AM  ! Amount of ozone in the direction of the Sun
    To = 1.d0  -  0.1611d0 * Xo * (1.d0+139.48d0*Xo)**(-0.3035d0)  -  0.002715 * Xo / (1.d0 + 0.044d0*Xo + 0.0003d0*Xo**2)  ! Transmittance of ozone absorptance
    
    ! Uniformly mixed gases (CO2, O2):
    Tum = exp(-0.0127d0 * AMp**0.26d0)  ! Transmittance of mixed-gas absorptance
    
    ! Water vapor:
    Xw = AM*UwL  ! Amount of water vapor in the direction of the Sun
    Tw = 1.d0 - 2.4959d0 * Xw / ((1.d0 + 79.034d0*Xw)**0.6828d0 + 6.385d0*Xw)     ! Transmittance of water-vapor absorptance - Tw = 1-Aw
    
    ! Daily turbidity:
    Tau = 0.2758d0*Ta3L + 0.35d0*Ta5L                                             ! Broadband turbidity: aerosol optical depth from surface in a vertical column
    Ta  = exp( -Tau**0.873d0  *  (1.d0 + Tau - Tau**0.7088d0)  *  AM**0.9108d0 )  ! Transmittance of aerosol absorptance and scattering
    Taa = 1.d0 - K1L * (1.d0 - AM + AM**1.06d0) * (1.d0-Ta)                       ! Transmittance of aerosol absorptance - this does not agree with tabulated values from the paper (except for AM~1).  When I substitute their values for Taa, everything matches perfectly.  Error in their formula, or in their table?
    Tas = Ta/Taa                                                                  ! Transmittance of aerosol scattering
    Rs  = 0.0685d0 + (1.d0-BaL) * (1.d0-Tas)                                      ! Sky albedo
    
    
    ! IRRADIANCE EQUATIONS - Itot, Idir, Idif are always computed, Igr only if desired:
    
    ! Direct radiation on a horizontal surface:
    tmpVar = IoL * cosZ  *  To * Tum * Tw  ! Save a few CPU cycles
    IdirL = 0.9662d0 * tmpVar  *  Tr * Ta  /  RsunL**2
    
    ! Diffuse (scattered) radiation on a horizontal surface:
    IdifL  = 0.79d0 *  tmpVar        * Taa *  (0.5d0*(1.d0-Tr) + BaL*(1.d0-Tas)) / (1.d0 - AM + AM**1.02d0)
    
    ! Total (direct+diffuse) radiation on a horizontal surface:
    ItotL = (IdirL+IdifL) / (1.d0 - RgL*Rs)
    
    
    ! Copy local variables to optional output variables:
    if(present(Itot)) Itot = ItotL
    if(present(Idir)) Idir = IdirL
    if(present(Idif)) Idif = IdifL
    if(present(Igr))  Igr  = ItotL - (IdirL+IdifL)  ! Ground-reflected radiation from a horizontal surface
    
  end subroutine clearsky_bird
  !*********************************************************************************************************************************
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute diffuse radiation on an inclined surface using the older (1987) Perez model
  !!
  !! \param DoY           Day of year (Nday)
  !! \param alt           Altitude of the Sun (radians)
  !! 
  !! \param surfIncl      Surface inclination wrt horizontal (radians) - 0 = horizontal, pi/2 = vertical
  !! \param theta         Angle between surface normal vector and Sun position vector (radians)
  !!
  !! \param Gbeam_n       Beam (direct) normal radiation (W/m2; in the direction of the Sun)
  !! \param Gdif_hor      Diffuse radiation on a horizontal surface (W/m2)
  !!
  !! \param Gdif_inc      Diffuse irradiation on the inclined surface (W/m2) (output)
  !!
  !! \param Gdif_inc_is   Diffuse irradiation on the inclined surface - isotropic part (optional; W/m2) (output)
  !! \param Gdif_inc_cs   Diffuse irradiation on the inclined surface - circumsolar part (optional; W/m2) (output)
  !! \param Gdif_inc_hz   Diffuse irradiation on the inclined surface - horizon-band part (optional; W/m2) (output)
  !!
  !! \see Perez et al. Solar Energy Vol. 39, Nr. 3, p. 221 (1987) - references to equations and tables are to this paper.
  !! Most equations can be found in the Nomenclature section at the end of the paper (p.230).  I use a and c here, not b and d.
  !!
  !! \todo Implement Perez et al. Solar Energy Vol. 44, Nr. 5, p. 271 (1990)
  
  subroutine diffuse_radiation_Perez87(DoY, alt, surfIncl, theta, Gbeam_n,Gdif_hor,  Gdif_inc,   Gdif_inc_is, Gdif_inc_cs, Gdif_inc_hz)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi2,pio2, d2r,r2d
    
    implicit none
    integer, intent(in) :: DoY
    real(double), intent(in) :: alt, surfIncl, theta, Gbeam_n,Gdif_hor
    real(double), intent(out) :: Gdif_inc
    real(double), intent(out), optional :: Gdif_inc_is, Gdif_inc_cs, Gdif_inc_hz
    integer :: f11,f12,f13, f21,f22,f23
    real(double) :: zeta, AM0rad, Mair,Delta,epsilon, alpha, psiC,psiH, chiC,chiH, F(6),F1,F2, A,C
    real(double) :: Gdif_inc_iso, Gdif_inc_csl, Gdif_inc_hzl
    
    
    ! *** Compute the brightness coefficients for the isotropic (F1), circumsolar (F1) and horizon (F2) regions ***
    
    ! 'External' (AM0) radiation:
    AM0rad = 1370.d0 * (1.d0 + 0.00333d0 * cos(pi2/365.d0 * DoY))
    
    ! Air mass:
    if(alt .lt. -3.885d0*d2r) then
       Mair = 99.d0
    else if(alt .lt. 10.d0*d2r) then
       Mair = 1.d0 / ( sin(alt) + 0.15d0 * (alt*r2d + 3.885d0)**(-1.253d0) )
    else
       Mair = 1.d0 / sin(alt)
    end if
    Delta = Gdif_hor * Mair / AM0rad  ! Brightness of overcast sky - par. 2.2.4 (a)
    
    
    ! Cloud cover: epsilon;  epsilon ~ 1: overcast, epsilon -> infinity: clear  (epsilon ~ 1/fraction of covered sky)
    !   Needed for correct row in Table 1
    if(Gdif_hor.le.0.d0) then  ! Division by zero
       if(Gbeam_n.le.0.d0) then  ! No direct light: 0/0
          epsilon = 0.d0      !   -> completely overcast - first row of Table 1
       else                   ! Some direct daylight: x/0 = large
          epsilon = 99.d0     !   -> completely clear, should be >11 for last row of Table 1
       end if
    else
       epsilon = (Gdif_hor + Gbeam_n) / Gdif_hor  ! Overcast: epsilon ~ 1,  clear: epsilon -> infinity
    end if
    
    
    ! Table 1
    f11=1;  f12=2;  f13=3;  f21=4; f22=5; f23=6
    if(epsilon .le. 1.056d0) then
       F = [-0.011d0,  0.748d0, -0.080d0, -0.048d0,  0.073d0, -0.024d0]
    else if(epsilon .le. 1.253d0) then
       F = [-0.038d0,  1.115d0, -0.109d0, -0.023d0,  0.106d0, -0.037d0]
    else if(epsilon .le. 1.586d0) then
       F = [ 0.166d0,  0.909d0, -0.179d0,  0.062d0, -0.021d0, -0.050d0]
    else if(epsilon .le. 2.134d0) then
       F = [ 0.419d0,  0.646d0, -0.262d0,  0.140d0, -0.167d0, -0.042d0]
    else if(epsilon .le. 3.230d0) then
       F = [ 0.710d0,  0.025d0, -0.290d0,  0.243d0, -0.511d0, -0.004d0]
    else if(epsilon .le. 5.980d0) then
       F = [ 0.857d0, -0.370d0, -0.279d0,  0.267d0, -0.792d0,  0.076d0]
    else if(epsilon .le. 10.080d0) then
       F = [ 0.734d0, -0.073d0, -0.228d0,  0.231d0, -1.180d0,  0.199d0]
    else
       F = [ 0.421d0, -0.661d0,  0.097d0,  0.119d0, -2.125d0,  0.446d0]
    end if
    
    zeta = pio2 - alt  ! Zenith angle = pi/2 - alt
    F1 = F(f11)  +  F(f12) * Delta  +  F(f13) * zeta  ! Isotropic, circumsolar brightness coefficient
    F2 = F(f21)  +  F(f22) * Delta  +  F(f23) * zeta  ! Horizon brightness coefficient
    
    
    
    
    ! *** Compute the mean solid angles occupied by the circumsolar region (C and A, needed for Eq.8) ***
    
    alpha = 25.d0*d2r  ! Half angle of the circumsolar region (degrees -> radians; below Eq.9)
    
    
    ! Solid angle of the circumsolar region weighted by incidence on the HORIZONTAL (variable C, subscript H; 
    !   see Nomenclature, under c):
    ! psiH:
    if(zeta .gt. pio2 - alpha) then
       psiH = 0.5d0 * (pio2 - zeta + alpha) / alpha  ! Dimensionless ratio
    else
       psiH = 1.d0
    end if
    
    ! chiH:
    if(zeta .lt. pio2 - alpha) then
       chiH = cos(zeta)  ! = sin(alt)
    else
       chiH = psiH * sin(psiH*alpha)
    end if
    
    C = 2 * (1.d0 - cos(alpha)) * chiH  ! Solid angle of the circumsolar region, weighted by HORIZONTAL incidence
    
    
    ! Solid angle of the circumsolar region weighted by incidence on the SLOPE (variable A, subscript C; 
    !   see Nomenclature, under c):
    ! psiC:
    psiC = 0.5d0 * (pio2 - theta + alpha) / alpha
    
    ! chiC:
    if(theta .lt. pio2 - alpha) then
       chiC = psiH * cos(theta)
    else if(theta .lt. pio2 + alpha) then
       chiC = psiH * psiC * sin(psiC*alpha)
    else
       chiC = 0.d0
    end if
    
    A = 2 * (1.d0 - cos(alpha)) * chiC  ! Solid angle of the circumsolar region, weighted by SLOPE incidence
    
    
    
    ! Diffuse radiation from circumsolar (F1) and horizon (F2) regions on the inclined surface (Eq.8):
    Gdif_inc_iso = Gdif_hor * 0.5d0 * (1.d0 + cos(surfIncl)) * (1.d0 - F1)  ! Isotropic
    Gdif_inc_csl = Gdif_hor * F1 * A/C                                      ! Circumsolar
    Gdif_inc_hzl = Gdif_hor * F2 * sin(surfIncl)                            ! Horizon band
    
    Gdif_inc = max(Gdif_inc_iso + Gdif_inc_csl + Gdif_inc_hzl, 0.d0)  ! Note: components are sometimes negative
    
    ! Assign optional return values:
    if(present(Gdif_inc_is)) Gdif_inc_is = Gdif_inc_iso
    if(present(Gdif_inc_cs)) Gdif_inc_cs = Gdif_inc_csl
    if(present(Gdif_inc_hz)) Gdif_inc_hz = Gdif_inc_hzl
    
  end subroutine diffuse_radiation_Perez87
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the projection of direct solar radiation (DNI) on a (sloped) surface.
  !!
  !! This function returns the projection factor of direct sunlight on a surface with given orientation,
  !! useful for e.g. insolation on solar panels, solar collectors or windows.  The projection factor equals
  !! the cosine of the angle between the normal vector of the surface and the position vector of the Sun.
  !! 
  !! \param beta    Inclination angle of the surface (w.r.t. the horizontal) (radians)
  !! \param gamma   Azimuth angle of the surface (0=south, pi/2=west, ±pi=north, -pi/2=east) (radians)
  !! \param sunAz   Azimuth of the Sun (0=south, pi/2=west, ±pi=north, -pi/2=east) (radians)
  !! \param sunAlt  Apparent altitude of the Sun (radians)
  !!
  !! \note
  !! Note that different definitions for gamma and sunAz are possible, as long as they correspond.
  !!
  !! \see
  !! Celestial mechanics in a nutshell, Sect. 4.3: Insolation on an inclined surface  (http://CMiaNS.sf.net).
  
  function project_sunlight_on_surface(beta,gamma, sunAz,sunAlt)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: beta,gamma, sunAz,sunAlt
    real(double) :: project_sunlight_on_surface
    
    project_sunlight_on_surface = sin(sunAlt)*cos(beta) + cos(sunAlt)*sin(beta)*cos(sunAz - gamma)

    
  end function project_sunlight_on_surface
  !*********************************************************************************************************************************
  
  
  
  
end module TheSky_daylight
!***********************************************************************************************************************************
