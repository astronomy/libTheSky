!> \file asteroids.f90  Procedures to compute asteroid positions and more for libTheSky


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


!***********************************************************************************************************************************
!> \brief Procedures for asteroids

module TheSky_asteroids
  implicit none
  save
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate orbital elements for asteroids
  !!
  !! \param tjm  Dynamical time in Julian Millennia after 2000.0
  !! \param as   Asteroid ID
  !!
  !! \param hcr   Heliocentric distance (AU) (output)
  !! \param om1   Argument of perihelion (rad) (output)
  !! \param nu    True anomaly (rad) (output)
  !!
  !! \param cosi  cos(i): cosine of inclination (output)
  !! \param sini  sin(i): sine of inclination (output)
  !! \param coso  cos(Om): cosine of longitude of ascending node (output)
  !! \param sino  sin(Om): sine of longitude of ascending node (output)
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.33
  !!
  !! \note Used by asteroid_lbr() and asteroid_eq()

  
  subroutine asteroid_elements(tjm,as, hcr, om1,nu, cosi,sini,coso,sino)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000
    
    use TheSky_coordinates, only: precess_orb
    use TheSky_planetdata, only: asterelems
    
    real(double), intent(in) :: tjm
    integer, intent(in) :: as
    real(double), intent(out) :: hcr, om1,nu, sini,cosi,sino,coso
    real(double) :: jd,jd0,djd, sma,ecc,in,om2, m0,mea,mem,ea
    
    
    ! Elements are: 1-Epoch (JD), 2-a, 3-e, 4-i, 5-omega, 6-Omega, 7-M, 8-H, 9-G, for J2000.0:
    jd0 = asterelems(as,1)  ! JD of epoch
    sma = asterelems(as,2)  ! Semi-major axis
    ecc = asterelems(as,3)  ! Eccentricity
    in  = asterelems(as,4)  ! Inlination
    om1 = asterelems(as,5)  ! omega - argument of perihelion
    om2 = asterelems(as,6)  ! Omega - longitude of ascending node
    m0  = asterelems(as,7)  ! Mean anomaly at epoch
    
    ! Test Meeus p.232:
    !jd0 = 2448193.04502000d0
    !sma = 2.2091404d0
    !ecc = 0.8502196d0
    !in  = 11.94524d0*d2r
    !om1 = 186.23352d0*d2r
    !om2 = 334.75006d0*d2r
    !m0 = 0.d0
    
    mem = 0.01720209895d0 * sma**(-1.5d0)   ! Mean motion (rad/day)
    jd  = jd2000 + 365250.d0 * tjm
    djd = jd - jd0                          ! Time since epoch in days
    mea = m0 + mem*djd                      ! Mean anomaly at time t
    
    call precess_orb(jd2000,jd,in,om1,om2)  ! Precess elements to equinox of date (for compatibility with planet_position)
    
    call kepler(ecc,mea, ea)                ! Solve Kepler's equation for eccentricic anomaly ea
    
    nu  = 2*atan(sqrt((1.d0+ecc)/(1.d0-ecc)) * tan(0.5d0*ea))  ! True anomaly, Meeus Eq. 30.1
    hcr = sma*(1.d0 - ecc*cos(ea))                             ! Radius vector, Meeus Eq. 30.2
    
    sini = sin(in)
    cosi = cos(in)
    sino = sin(om2)
    coso = cos(om2)
    
  end subroutine asteroid_elements
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate heliocentric ecliptical coordinates l,b,r for asteroid as at time t (Julian Millennia after 2000.0 in DT)
  !!
  !! \param t    Dynamical time in Julian Millennia after 2000.0
  !! \param as   Asteroid ID
  !!
  !! \param l   Heliocentric ecliptical longitude (rad) (output)
  !! \param b   Heliocentric ecliptical latitude (rad) (output)
  !! \param r   Heliocentric distance (AU) (output)
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.33
  
  subroutine asteroid_lbr(t,as, l,b,r)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev
    
    implicit none
    real(double), intent(in) :: t
    integer, intent(in) :: as
    real(double), intent(out) :: l,b,r
    
    real(double) :: om1,nu, sini,cosi,sino,coso,u,sinu,cosu,x,y,z
    
    
    ! Compute orbital elements:
    call asteroid_elements(t,as, r, om1, nu, cosi,sini, coso,sino)
    
    ! To calculate *ecliptical* coordinates:
    ! Meeus, p.233:
    u = om1 + nu
    cosu = cos(u)
    sinu = sin(u)
    
    ! Heliocentric rectangular position:
    x = r * (coso * cosu  -  sino * sinu * cosi)
    y = r * (sino * cosu  +  coso * sinu * cosi)
    z = r *  sini * sinu
    
    l = rev(atan2(y,x))  ! Heliocentric ecliptic longitude
    b = asin(z/r)        ! Heliocentric ecliptic latitude
    
  end subroutine asteroid_lbr
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate geocentric(?) equatorial coordinates ra, dec, delta for asteroid as at time t (Jul. Millennia DT after 2000)
  !!
  !! \param t    Dynamical time in Julian Millennia after 2000.0
  !! \param as   Asteroid ID
  !! \param l0   Heliocentric longitude of Earth
  !! \param b0   Heliocentric latitude of Earth
  !! \param r0   Heliocentric distance of Earth
  !!
  !! \param ra     Geocentric right ascension (rad) (output)
  !! \param dec    Geocentric declination (rad) (output)
  !! \param delta  Geocentric distance (AU) (output)
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.33
  
  subroutine asteroid_eq(t,as, l0,b0,r0,  ra,dec,delta)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev
    use SUFR_constants, only: eps2000
    
    use TheSky_coordinates, only: calcsunxyz
    
    implicit none
    real(double), intent(in) :: t, l0,b0,r0
    integer, intent(in) :: as
    real(double), intent(out) :: ra,dec,delta
    
    real(double) :: om1,nu,  eps,sine,cose,  ff,gg,hh,pp,qq,rr, a1,a2,b1,b2,c1,c2  ! jd, dpsi,eps0,deps
    real(double) :: sini,cosi,sino,coso,  x,y,z,  r, sx,sy,sz
    
    
    ! Compute orbital elements:
    call asteroid_elements(t,as, r, om1, nu, cosi,sini, coso,sino)
    
    ! To calculate *equatorial* coordinates:
    !jd  = jd2000 + 365250.d0 * t
    !call nutation(t, dpsi,eps0,deps)
    !call nutation2000(jd,dpsi,deps)  ! IAU 2000 Nutation model - doesn't provide eps0
    !eps = eps0 + deps
    
    eps  = eps2000   ! For J2000.0  - this is used by JPL - CHECK
    sine = sin(eps)
    cose = cos(eps)
    
    ! Meeus, Eq. 33.7:
    ff =  coso
    gg =  sino*cose
    hh =  sino*sine
    pp = -sino*cosi
    qq =  coso*cosi*cose - sini*sine
    rr =  coso*cosi*sine + sini*cose
    
    ! Meeus, Eq. 33.8:
    a1 = rev( atan2(ff,pp) )
    b1 = rev( atan2(gg,qq) )
    c1 = rev( atan2(hh,rr) )
    a2 = sqrt( ff**2 + pp**2 )
    b2 = sqrt( gg**2 + qq**2 )
    c2 = sqrt( hh**2 + rr**2 )
    
    ! Meeus, Eq. 33.9:
    x = r * a2 * sin(a1 + om1 + nu)  ! Heliocentric *EQUATORIAL* rectangular coordinates: x,y,z
    y = r * b2 * sin(b1 + om1 + nu)
    z = r * c2 * sin(c1 + om1 + nu)
    
    ! Heliocentric -> geocentric:
    call calcsunxyz(t, l0,b0,r0, sx,sy,sz)   ! Geocentric *equatorial* rectangular solar coordinates
    
    ! Meeus, Eq. 33.10:
    x = x + sx  ! heliocentric -> geocentric
    y = y + sy  ! heliocentric -> geocentric
    z = z + sz  ! heliocentric -> geocentric
    
    delta = sqrt(x**2 + y**2 + z**2)  ! Geocentric distance (AU)
    ra    = atan2(y,x)                ! Geocentric right ascenson
    dec   = asin(z/delta)             ! Geocentric declination
    
  end subroutine asteroid_eq
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate asteroid magnitude for asteroid as
  !! 
  !! \param as     Asteroid number
  !! \param delta  Distance to the Earth
  !! \param r      Distance to the Sun
  !! \param pa     Phase angle
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.33
  
  function asteroid_magn(as, delta,r, pa)
    use SUFR_kinds, only: double
    use TheSky_planetdata, only: asterelems
    
    implicit none
    integer, intent(in) :: as
    real(double), intent(in) :: delta,r,pa
    real(double) :: asteroid_magn,h,g,phi1,phi2
    
    h = asterelems(as,8)  ! Mean absolute visual magnitude
    g = asterelems(as,9)  ! 'Slope parameter'
    
    ! Meeus, p.231:
    phi1 = exp(-3.33d0 * tan(0.5d0*pa)**(0.63d0))
    phi2 = exp(-1.87d0 * tan(0.5d0*pa)**(1.22d0))
    
    asteroid_magn = h + 5.d0*log10(r*delta) - 2.5d0*log10((1.d0-g)*phi1 + g*phi2)  ! Meeus, Eq. 33.14
    
  end function asteroid_magn
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Solve Kepler's equation
  !!
  !! \param  ecc  Eccentricity
  !! \param  mea  Mean anomaly
  !! \param eca   Eccentric anomaly (output)
  !!
  !! Solve Kepler's equation for given eccentricity ecc, mean anomaly mea and eccentric anomaly eca
  !! \see Meeus, Astronomical Algorithms, 1998, p.206
  
  subroutine kepler(ecc,mea, eca)
    use SUFR_kinds, only: double, dbl
    use SUFR_constants, only: pi,pi2, pio2,pio4
    use SUFR_system, only: warn
    use SUFR_angles, only: rev
    
    implicit none
    real(double), intent(in) :: ecc,mea
    real(double), intent(out) :: eca
    
    real(double), parameter :: accur = 1.0e-15_dbl  ! Accur: 1e-10: 34 iterations, 1e-15: 51 iter, 1e-50: 167, 0.d0: 1076
    integer :: sig
    real(double) :: ma,d,de,m1
    
    if(ecc.lt.0.d0 .or. ecc.gt.1.d0)  call warn('kepler() can only be used for orbits with eccentricities between 0 and 1', 0)
    
    ma = mea
    !sig = nint(abs(m)/ma)  ! Sign of ma
    sig = nint(sign(1.d0,ma))  ! Sign of ma
    ma = abs(ma)/pi2
    ma = rev( sig * pi2 * (ma-floor(ma)) )  ! Is INT() in basic really floor() ?
    
    sig = 1
    if(ma.gt.pi) then
       sig = -1
       ma = pi2 - ma
    end if
    
    ! Starting values:
    eca = pio2      ! pi/2
    d   = pio4      ! pi/4
    de  = huge(de)
    
    do while(abs(de).gt.accur)  ! Takes 51 iterations
       m1  =  eca - ecc*sin(eca)
       de  =  d * sign(1.d0,ma-m1)
       eca =  eca + de
       d   =  0.5d0*d
    end do
    
    eca = sig * eca  ! Apply the correct sign
    
  end subroutine kepler
  !*********************************************************************************************************************************
  
  
  
  
end module TheSky_asteroids
!***********************************************************************************************************************************
