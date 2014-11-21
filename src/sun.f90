!> \file sun.f90  Procedures that calculate a low-accuracy position of the Sun for libTheSky


!  Copyright (c) 2002-2014  Marc van der Sluys - marc.vandersluys.nl
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


!***********************************************************************************************************************************
!> \brief Low-accuracy procedures for the Sun

module TheSky_sun
  implicit none
  save
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Low-accuracy solar coordinates (~0.01deg)
  !! 
  !! \param jd    Julian Day of computation
  !! \param calc  Calculate:  1: l,b,r,  2: & ra,dec,  3: & gmst,agst,  4: & az,alt,  5: & mag + p.a.,  99: & topo alt + refraction
  !!
  !! \note
  !! - lat0, lon0 and height must be provided through the module TheSky_local (rad, rad and m)
  !! - results are returned in the array planpos() in the module TheSky_planetdata
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.25
  !!
  !! \todo  odot is off by ~10" in Meeus, Example 25a.  Would need better L0 or C (or M?)
  
  subroutine sunpos_la(jd, calc)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000, au, earthr,rsun
    use SUFR_angles, only: rev
    
    use TheSky_planetdata, only: planpos
    use TheSky_coordinates, only: eq2horiz, refract
    use TheSky_datetime, only: calc_deltat, calc_gmst
    use TheSky_local, only: lon0,lat0
    
    implicit none
    real(double), intent(in) :: jd
    integer, intent(in) :: calc
    real(double) :: jde,deltat,tjc,tjc2, l0,m,e,c,odot,nu,omg,r, aber, lam,b
    real(double) :: ra,dec,lm,eps,eps0,deps,gmst,agst,lst,dpsi,az,alt,hh
    
    deltat = calc_deltat(jd)
    jde    = jd + deltat/86400.d0
    !tjc    = (jde-jd2000)/36525.d0    ! Julian Centuries after 2000.0 in dynamical time, the T in Meeus, p.163, Eq. 25.1
    tjc    = (jd-jd2000)/36525.d0     ! Julian Centuries after 2000.0;  JD compares better than JD (2015-2114)
    tjc2   = tjc*tjc                  ! T^2                             In fact JD - DeltaT compares even better(!)
    
    l0 = 4.895063168d0 + 628.331966786d0 *tjc + 5.291838d-6    *tjc2  ! Mean longitude, Eq. 25.2
    !l0 = 4.895063113086d0 + 628.331965406500135d0 *tjc + 5.29213354358d-6 *tjc2 + 3.4940522d-10 *tjc2*tjc - &
    !     1.1407666d-10 *tjc2*tjc2 - 8.726646d-14 *tjc2*tjc2*tjc + 9.696d-16 * tjc2*tjc2*tjc2  ! Simon et al, 5.9.3
    m  = 6.240060141d0 + 628.301955152d0 *tjc - 2.682571d-6    *tjc2  ! Mean anomaly, Eq. 25.3
    !m  = 6.24006012726d0 + 628.301955167d0*tjc - 2.680826d-6*tjc2  - 7.1267d-10*tjc2*tjc      ! Mean anomaly, Meeus Eq.47.3
    
    
    ! Eccentricity of the Earth's orbit:
    e  = 0.016708634d0 - 0.000042037d0   *tjc - 0.0000001267d0 *tjc2  ! Meeus, Eq. 25.4
    !e  = 0.0167086342d0 - 0.00004203654d0 *tjc - 0.000000126734d0 *tjc2 + 1.444d-10*tjc2*tjc - 2.d-14*tjc2*tjc2 &
    !     + 3.d-15*tjc2*tjc2*tjc ! Simon ea, 5.9.3
    
    ! Sun's equation of the centre:
    c = (3.34161088d-2 - 8.40725d-5*tjc - 2.443d-7*tjc2)*sin(m)  +  (3.489437d-4 - 1.76278d-6*tjc)*sin(2*m) + 5.044d-6*sin(3*m)
    odot = rev(l0 + c)  ! True longitude
    nu = rev(m + c)     ! True anomaly
    r = 1.000001018d0*(1.d0-e*e)/(1.d0 + e*cos(nu))  ! Heliocentric distance of the Earth / geocentric dist. of the Sun, Eq. 25.5
    
    ! Longitude of the Moon's mean ascending node, Eq.47.7:
    omg  = 2.1824390725d0 - 33.7570464271d0 *tjc  + 3.622256d-5 *tjc2 + 3.7337958d-8 *tjc2*tjc - 2.879321d-10 *tjc2*tjc2
    
    ! Mean longitude of the Moon:
    lm   = 3.8103417d0 + 8399.709113d0*tjc      ! Meeus, p.144
    !lm   = 3.81034102259d0 + 8399.709113521569d0*tjc - 2.315703d-5*tjc2 + 3.23904315d-8*tjc2*tjc - 2.6771317d-10*tjc2*tjc2 !Eq.47.1
    
    dpsi = -8.338795d-5*sin(omg) - 6.39954d-6*sin(2*l0) - 1.115d-6*sin(2*lm) + 1.018d-6*sin(2*omg)  ! Nutation in lon, Meeus, p.144
    aber = -9.93087d-5/r            ! Aberration, Meeus Eq.25.10
    lam  = rev(odot + aber + dpsi)  ! Apparent geocentric longitude, referred to the true equinox of date
    b    = 0.d0
    
    
    planpos(1)  = lam  ! Geocentric longitude
    planpos(2)  = b    ! Geocentric latitude
    planpos(3)  = r    ! Geocentric distance
    planpos(4)  = r    ! Geocentric distance
    
    ! Secondary variables, (almost) for free:
    planpos(7)  = 5.77551830441d-3 * r    ! Light time in days
    planpos(11) = 0.d0                    ! Elongation
    planpos(12) = 2*rsun/(r*au)           ! Apparent diameter
    planpos(13) = -26.73d0                ! Apparent magnitude
    planpos(14) = 1.d0                    ! Illuminated fraction
    planpos(15) = 0.d0                    ! Phase angle
    planpos(17) = asin(earthr/(r*au))     ! Horizontal parallax
    
    planpos(40) = jde                     ! JD
    planpos(41:43) = planpos(1:3)         ! Geocentric, "true" L,B,R of the Sun
    planpos(46) = tjc                     ! App. dyn. time in Julian Centuries since 2000.0
    
    if(calc.eq.1) return
    
    
    
    ! Obliquity of the ecliptic and nutation:
    eps0 = 0.409092804222d0 - 2.26965525d-4*tjc - 2.86d-9*tjc2 + 8.78967d-9*tjc2*tjc       ! Mean obliq. o.t. eclip, Meeus, Eq.22.2
    deps = 4.46d-5*cos(omg) + 2.76d-6*cos(2*l0) + 4.848d-7*cos(2*lm) - 4.36d-7*cos(2*omg)  ! Nutation in obliquity, Meeus, p.144
    eps  = eps0 + deps                                                                     ! True obliquity of the ecliptic
    
    ra   = atan2(cos(eps)*sin(lam),cos(lam))  ! Geocentric right ascension, Eq. 25.6
    dec  = asin(sin(eps)*sin(lam))            ! Geocentric declination,     Eq. 25.7
    
    
    planpos(5)  = ra    ! Geocentric right ascension
    planpos(6)  = dec   ! Geocentric declination
    planpos(47) = dpsi  ! Nutation in longitude
    planpos(48) = eps   ! True obliquity of the ecliptic; corrected for nutation
    planpos(50) = eps0  ! Mean obliquity of the ecliptic; without nutation
    
    if(calc.eq.2) return
    
    
    
    gmst = calc_gmst(jd)              ! Greenwich mean siderial time
    agst = rev(gmst + dpsi*cos(eps))  ! Correction for equation of the equinoxes -> Gr. apparent sid. time
    lst  = rev(agst + lon0)           ! Local apparent siderial time, lon0 > 0 for E
    
    planpos(44) = lst                 ! Local APPARENT siderial time
    planpos(45) = rev(agst)           ! Greenwich mean siderial time
    planpos(49) = rev(gmst)           ! Correction for equation of the equinoxes -> Gr. apparent sid. time
    
    if(calc.eq.3) return
    
    
    
    call eq2horiz(ra,dec,agst, hh,az,alt)
    planpos(8)  = rev(hh)  ! Geocentric hour angle
    planpos(9)  = rev(az)  ! Geocentric azimuth
    planpos(10) = alt      ! Geocentric altitude
    
    if(calc.eq.4) return
    
    
    planpos(13) = sunmagn(r)                                            ! Apparent magnitude
    planpos(16) = atan2(sin(hh),tan(lat0)*cos(dec) - sin(dec)*cos(hh))  ! Parallactic angle
    
    if(calc.eq.5) return
    
    ! calc = 99:
    planpos(30) = alt - asin(sin(planpos(17))*cos(alt))  ! alt' = alt - Hp*cos(alt); topocentric altitude
    planpos(31) = planpos(30) + refract(planpos(30))     ! Topocentric altitude, corrected for atmospheric refraction
    
  end subroutine sunpos_la
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate Sun magnitude
  !!
  !! \param  dist     Distance (AU)
  !! \retval sunmagn  Sun magnitude
  
  function sunmagn(dist)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: dist
    real(double) :: sunmagn,ill
    
    ill = 10.d0**(-0.4d0 * (-26.73d0)) / (dist*dist)  ! Dist in AU, assume average dist = 1
    sunmagn = -2.5d0*log10(ill) 
    
  end function sunmagn
  !*********************************************************************************************************************************
  
  
end module TheSky_sun
!***********************************************************************************************************************************
