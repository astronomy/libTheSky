!> \file coordinates.f90  Procedures to perform coordinate transformations, apply precession, and more for libTheSky


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
!> \brief Procedures for coordinates

module TheSky_coordinates
  implicit none
  save
  
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the geocentric rectangular coordinates of a planet,  from its and the Earth's heliocentric spherical position
  !!
  !! \param l   Heliocentric longitude of planet
  !! \param b   Heliocentric latitude of planet
  !! \param r   Heliocentric distance of planet
  !!
  !! \param l0  Heliocentric longitude of Earth
  !! \param b0  Heliocentric latitude of Earth
  !! \param r0  Heliocentric distance of Earth
  !!
  !! \retval x  Geocentric rectangular X of planet
  !! \retval y  Geocentric rectangular Y of planet
  !! \retval z  Geocentric rectangular Z of planet
  
  subroutine hc_spher_2_gc_rect(l,b,r, l0,b0,r0, x,y,z)  
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: l,b,r, l0,b0,r0
    real(double), intent(out) :: x,y,z
    
    x = r * cos(b) * cos(l)  -  r0 * cos(b0) * cos(l0)
    y = r * cos(b) * sin(l)  -  r0 * cos(b0) * sin(l0)
    z = r * sin(b)           -  r0 * sin(b0)
    
  end subroutine hc_spher_2_gc_rect
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert spherical, ecliptical coordinates to rectangular, equatorial coordinates of an object - both geocentric
  !!
  !! \param l    Heliocentric longitude of planet
  !! \param b    Heliocentric latitude of planet
  !! \param r    Heliocentric distance of planet
  !! \param eps  Obliquity of the ecliptic
  !!
  !! \retval x   Geocentric rectangular X of planet
  !! \retval y   Geocentric rectangular Y of planet
  !! \retval z   Geocentric rectangular Z of planet
  
  subroutine ecl_spher_2_eq_rect(l,b,r, eps,  x,y,z)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: l,b,r,eps
    real(double), intent(out) :: x,y,z
    
    x = r *  cos(b) * cos(l)
    y = r * (cos(b) * sin(l) * cos(eps)  -  sin(b) * sin(eps))
    z = r * (cos(b) * sin(l) * sin(eps)  +  sin(b) * cos(eps))
    
  end subroutine ecl_spher_2_eq_rect
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the geocentric equatorial rectangular coordinates of the Sun, from Earth's heliocentric, spherical position
  !!
  !! \param t1  Dynamical time in Julian Millennia after 2000.0
  !! \param l0  Heliocentric longitude of Earth
  !! \param b0  Heliocentric latitude of Earth
  !! \param r0  Heliocentric distance of Earth
  !!
  !! \retval x  Geocentric rectangular X of planet
  !! \retval y  Geocentric rectangular Y of planet
  !! \retval z  Geocentric rectangular Z of planet
  !!
  !! \todo  To be disposed, has been replcaed by ecl_spher_2_eq_rect above 
  
  subroutine calcsunxyz(t1, l0,b0,r0, x,y,z)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi
    use SUFR_angles, only: rev
    
    use TheSky_nutation, only: nutation
    
    implicit none
    real(double), intent(in) :: t1, l0,b0,r0
    real(double), intent(out) :: x,y,z
    real(double) :: t, l,b,r, dpsi,eps,eps0,deps
    
    t = t1
    l = rev(l0+pi)  ! Heliocentric longitude of the Earth -> geocentric longitude of the Sun
    b = -b0         ! Heliocentric latitude of the Earth  -> geocentric latitude of the Sun
    r = r0          ! Heliocentric distance of the Earth  =  geocentric distance of the Sun
    
    call fk5(t, l,b)
    call nutation(t,dpsi,eps0,deps)
    eps = eps0
    
    x = r *  cos(b) * cos(l)
    y = r * (cos(b) * sin(l) * cos(eps)  -  sin(b) * sin(eps))
    z = r * (cos(b) * sin(l) * sin(eps)  +  sin(b) * cos(eps))
    
    !jd1 = t1*365250.d0 + jd2000      !From equinox of date...
    !jd2 = comepoche                      ! ... to equinox of elements
    !call precess_xyz(jd1,jd2,x,y,z)
    
  end subroutine calcsunxyz
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the precession of the equinoxes in rectangular coordinates, from jd1 to jd2
  !!
  !! \param jd1  Original Julian day 
  !! \param jd2  Target Julian day 
  !!
  !! \retval x  Geocentric rectangular X
  !! \retval y  Geocentric rectangular Y
  !! \retval z  Geocentric rectangular Z
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch. 21
  
  subroutine precess_xyz(jd1,jd2, x,y,z)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000
    
    implicit none
    real(double), intent(in) :: jd1,jd2
    real(double), intent(inout) :: x,y,z
    real(double) :: t1,t2,t22,t23, x1,y1,z1
    real(double) :: dz,ze,th,sd,cd,sz,cz,st,ct
    real(double) :: xx,xy,xz,yx,yy,yz,zx,zy,zz
    
    t1 = (jd1 - jd2000)/36525.d0  !  t since 2000.0 in Julian centuries
    t2 = (jd2 - jd1)/36525.d0         ! dt in Julian centuries
    t22 = t2*t2   ! t2^2
    t23 = t22*t2  ! t2^3
    
    ! Meeus, Eq. 21.2:
    dz = (1.11808609d-2 + 6.770714d-6*t1 - 6.739d-10*t1*t1)*t2 + (1.463556d-6 - 1.668d-9*t1)*t22 + 8.725677d-8*t23
    ze = (1.11808609d-2 + 6.770714d-6*t1 - 6.739d-10*t1*t1)*t2 + (5.307158d-6 + 3.20d-10*t1)*t22 + 8.825063d-8*t23
    th = (9.71717346d-3 - 4.136915d-6*t1 - 1.0520d-9*t1*t1)*t2 - (2.068458d-6 + 1.052d-9*t1)*t22 - 2.028121d-7*t23
    
    sd = sin(dz)
    cd = cos(dz)
    sz = sin(ze)
    cz = cos(ze)
    st = sin(th)
    ct = cos(th)
    
    xx =  cd*cz*ct - sd*sz
    xy =  sd*cz    + cd*sz*ct
    xz =  cd*st
    yx = -cd*sz    - sd*cz*ct
    yy =  cd*cz    - sd*sz*ct
    yz = -sd*st
    zx = -cz*st
    zy = -sz*st
    zz =  ct
    
    x1 = xx*x + yx*y + zx*z
    y1 = xy*x + yy*y + zy*z
    z1 = xz*x + yz*y + zz*z
    
    x = x1
    y = y1
    z = z1
    
  end subroutine precess_xyz
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the precession of the equinoxes in equatorial coordinates, from jd1 to jd2
  !!
  !! \param jd1  Original Julian day 
  !! \param jd2  Target Julian day 
  !! \retval a1  Right ascension (IO, rad)
  !! \retval d1  Declination (IO, rad)
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.21, p.134
  
  subroutine precess_eq(jd1,jd2,a1,d1)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000
    
    implicit none
    real(double), intent(in) :: jd1,jd2
    real(double), intent(inout) :: a1,d1
    real(double) :: t1,t2,t22,t23, dz,ze,th,a,b,c
    
    t1 = (jd1 - jd2000)/36525.d0      !  t since 2000.0 in Julian centuries
    t2 = (jd2 - jd1)/36525.d0         ! dt in Julian centuries
    t22 = t2*t2   ! t2^2
    t23 = t22*t2  ! t2^3
    
    ! Meeus, Eq. 21.2:
    dz = (1.11808609d-2 + 6.770714d-6*t1 - 6.739d-10*t1*t1)*t2 + (1.463556d-6 - 1.668d-9*t1)*t22 + 8.725677d-8*t23
    ze = (1.11808609d-2 + 6.770714d-6*t1 - 6.739d-10*t1*t1)*t2 + (5.307158d-6 + 3.20d-10*t1)*t22 + 8.825063d-8*t23
    th = (9.71717346d-3 - 4.136915d-6*t1 - 1.0520d-9*t1*t1)*t2 - (2.068458d-6 + 1.052d-9*t1)*t22 - 2.028121d-7*t23
    
    ! Meeus, Eq. 21.4:
    a = cos(d1)           * sin(a1+dz)
    b = cos(th) * cos(d1) * cos(a1+dz)  -  sin(th) * sin(d1)
    c = sin(th) * cos(d1) * cos(a1+dz)  +  cos(th) * sin(d1)
    
    a1 = atan2(a,b) + ze
    d1 = asin(c)
    
  end subroutine precess_eq
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the precession of the equinoxes in equatorial coordinates, from yr1 to yr2
  !!
  !! \param yr1  Original year
  !! \param yr2  Target year
  !! \retval a1  Right ascension (IO, rad)
  !! \retval d1  Declination (IO, rad)
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.21, p.134
  
  subroutine precess_eq_yr(yr1,yr2,a1,d1)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: yr1,yr2
    real(double), intent(inout) :: a1,d1
    real(double) :: t1,t2,t22,t23, dz,ze,th,a,b,c
    
    t1 = (yr1 - 2000.d0)/100.d0  !  t since 2000.0 in Julian centuries
    t2 = (yr2 - yr1)/100.d0      ! dt in Julian centuries
    t22 = t2*t2   ! t2^2
    t23 = t22*t2  ! t2^3
    
    ! Meeus, Eq. 21.2:
    dz = (1.11808609d-2 + 6.770714d-6*t1 - 6.739d-10*t1*t1)*t2 + (1.463556d-6 - 1.668d-9*t1)*t22 + 8.725677d-8*t23
    ze = (1.11808609d-2 + 6.770714d-6*t1 - 6.739d-10*t1*t1)*t2 + (5.307158d-6 + 3.20d-10*t1)*t22 + 8.825063d-8*t23
    th = (9.71717346d-3 - 4.136915d-6*t1 - 1.0520d-9*t1*t1)*t2 - (2.068458d-6 + 1.052d-9*t1)*t22 - 2.028121d-7*t23
    
    ! Meeus, Eq. 21.4:
    a = cos(d1)           * sin(a1+dz)
    b = cos(th) * cos(d1) * cos(a1+dz)  -  sin(th) * sin(d1)
    c = sin(th) * cos(d1) * cos(a1+dz)  +  cos(th) * sin(d1)
    
    a1 = atan2(a,b) + ze
    d1 = asin(c)
    
  end subroutine precess_eq_yr
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the precession of the equinoxes in geocentric ecliptical coordinates
  !!
  !! \param  jd1  Original Julian day 
  !! \param  jd2  Target Julian day 
  !! \retval l    Ecliptic longitude (I/O, rad)
  !! \retval b    Ecliptic latitude (I/O, rad)
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.21, p.136
  
  subroutine precess_ecl(jd1,jd2, l,b)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000
    
    implicit none
    real(double), intent(in) :: jd1,jd2
    real(double), intent(inout) :: l,b
    real(double) :: t1,t2,t22,t23, eta,pii,p,aa,bb,cc
    
    t1  = (jd1 - jd2000)/36525.d0      !  t since 2000.0 in Julian centuries
    t2  = (jd2 - jd1)/36525.d0         ! dt in Julian centuries
    t22 = t2*t2   ! t2^2
    t23 = t22*t2  ! t2^3
    
    ! Meeus, Eq. 21.5  -  the different powers of t1/t2 have been checked and are ok:
    eta = (2.278765d-4    - 3.2012d-7     *t1  + 2.899d-9*t1*t1)   *t2  + (-1.60085d-7  + 2.899d-9*t1)   *t22  + 2.909d-10*t23
    pii = 3.0521686858d0  + 1.59478437d-2 *t1  + 2.939037d-6*t1*t1      - (4.2169525d-3 + 2.44787d-6*t1) *t2   + 1.7143d-7*t22
    p   = (2.438174835d-2 + 1.077382d-5   *t1  - 2.036d-10*t1*t1)  *t2  + (5.38691d-6   - 2.036d-10*t1)  *t22  - 2.91d-11*t23
    
    ! Meeus, Eq. 21.7:
    aa = cos(eta) * cos(b) * sin(pii-l)  -  sin(eta) * sin(b)
    bb = cos(b)            * cos(pii-l)
    cc = cos(eta) * sin(b)               +  sin(eta) * cos(b) * sin(pii-l)
    
    l = p + pii - atan2(aa,bb)
    b = asin(cc)
    
  end subroutine precess_ecl
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the precession of the equinoxes in orbital elements
  !!
  !! \param  jd1  Original Julian day 
  !! \param  jd2  Target Julian day 
  !!
  !! \retval i    Inclination
  !! \retval o1   Argument of perihelion
  !! \retval o2   Longitude of ascending node
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.24

  
  subroutine precess_orb(jd1,jd2,i,o1,o2)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000
    
    implicit none
    real(double), intent(in) :: jd1,jd2
    real(double), intent(inout) :: i,o1,o2
    real(double) :: t1,t2,t22,t23,  eta,pii,p,psi,aa,bb,cc,dd
    
    t1  = (jd1 - jd2000)/36525.d0      !  t since 2000.0 in Julian centuries
    t2  = (jd2 - jd1)/36525.d0         ! dt in Julian centuries
    t22 = t2*t2   ! t2^2
    t23 = t22*t2  ! t2^3
    
    ! Meeus, Eq. 21.5:
    eta = (2.278765d-4    - 3.2012d-7 * t1     + 2.899d-9*t1*t1) *t2  + (-1.60085d-7  + 2.899d-9*t1) * t22  + 2.909d-10 * t23
    pii =  3.0521686858d0 + 1.59478437d-2 * t1 + 2.939037d-6*t1 *t1   - (4.2169525d-3 + 2.44787d-6*t1) * t2 + 1.7143d-7 * t22
    p   = (2.438174835d-2 + 1.077382d-5 * t1   - 2.036d-10*t1*t1) *t2 + (5.38691d-6   - 2.036d-10*t1) * t22 - 2.91d-11 * t23
    psi = pii + p
    
    ! Meeus, Eq. 24.2:
    aa =  sin(i)   * sin(o2-pii)
    bb = -sin(eta) * cos(i)      + cos(eta) * sin(i)   * cos(o2-pii)
    cc = -sin(eta) * sin(o2-pii)
    dd =  sin(i)   * cos(eta)    - cos(i)   * sin(eta) * cos(o2-pii)
    
    i  = asin(sqrt(aa*aa+bb*bb))  ! Inclination
    o2 = atan2(aa,bb) + psi       ! Longitude of ascending node (Omega)
    o1 = o1 + atan2(cc,dd)        ! Argument of perihelion (omega)
    
  end subroutine precess_orb
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert rectangular coordinates x,y,z to spherical coordinates l,b,r
  !!
  !! \param  x  Rectangular x coordinate (same unit as r)
  !! \param  y  Rectangular y coordinate (same unit as r)
  !! \param  z  Rectangular z coordinate (same unit as r)
  !!
  !! \retval l  Longitude (rad)
  !! \retval b  Latitude  (rad)
  !! \retval r  Distance  (same unit as x,y,z)
  
  subroutine rect_2_spher(x,y,z, l,b,r) 
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev
    use SUFR_numerics, only: deq
    
    implicit none
    real(double), intent(in) :: x,y,z
    real(double), intent(out) :: l,b,r
    real(double) :: x2,y2
    
    if(deq(x,0.d0) .and. deq(y,0.d0) .and. deq(z,0.d0)) then
       l = 0.d0
       b = 0.d0
       r = 0.d0
    else
       x2 = x*x
       y2 = y*y
       
       l = rev( atan2(y,x) )        ! Longitude
       b = atan2(z, sqrt(x2 + y2))  ! Latitude
       r = sqrt(x2 + y2 + z*z)      ! Distance
    end if
    
  end subroutine rect_2_spher
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Correct ecliptic longitude and latitiude for annual aberration
  !!
  !! \param  t   Dynamical time in Julian Millennia since 2000.0
  !! \param  l0  Earth longitude (rad)
  !!
  !! \retval l   Longitude of the object (rad, I/O)
  !! \retval b   Latitude of the object (rad, I/O)
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.23
  
  subroutine aberration_ecl(t,l0, l,b)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi
    use SUFR_angles, only: rev
    
    implicit none
    real(double), intent(in) :: t,l0
    real(double), intent(inout) :: l,b
    real(double) :: tt,tt2, k,odot,e,pii,dl,db
    
    tt   = t*10   ! Time in Julian Centuries
    tt2  = tt*tt
    
    k    = 9.93650849745d-5  ! Constant of aberration, radians
    odot = rev(l0+pi)        ! Longitude of the Sun
    
    e    = 0.016708634d0 - 4.2037d-5 * tt   - 1.267d-7 * tt2
    pii  = 1.79659568d0  + 3.001024d-2 * tt + 8.0285d-6 * tt2
    
    dl   = k * ( -cos(odot-l)  +  e * cos(pii-l)) / cos(b)  ! Meeus, Eq. 23.2
    db   = -k * sin(b) * (sin(odot-l)  -  e * sin(pii-l))
    
    l    = l + dl
    b    = b + db
    
  end subroutine aberration_ecl
  !*********************************************************************************************************************************
  
  !*********************************************************************************************************************************
  !> \brief  Correct equatorial coordinates for annual aberration - moderate accuracy, use for stars
  !!
  !! \param  jd    Julian day of epoch
  !!
  !! \param  ra    Right ascension of the object (rad)
  !! \param  dec   Declination of the object (rad)
  !!
  !! \retval dra   Correction in right ascension due to aberration (rad)
  !! \retval ddec  Correction in declination due to aberration(rad)
  !!
  !! \param  eps0  The mean obliquity of the ecliptic
  !!
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.23
  
  subroutine aberration_eq(jd, ra,dec, dra, ddec, eps0)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000
    use SUFR_angles, only: rev
    use SUFR_dummy, only: dumdbl1,dumdbl2
    use SUFR_numerics, only: dne
    use TheSky_nutation, only: nutation
    
    implicit none
    real(double), intent(in) :: jd, ra,dec
    real(double), intent(out) :: dra,ddec
    real(double), intent(in), optional :: eps0
    real(double) :: jdold, tjc,tjc2, l0,mas,sec,odot, k,ee,pii, eps
    save jdold, tjc,odot,k,ee,pii, eps
    
    
    ! Compute these constants (independent of ra,dec) only if jd has changed since the last call
    if(dne(jd,jdold)) then
       
       ! Meeus, Ch. 25:
       tjc   =  (jd-jd2000)/36525.0d0  ! Julian centuries since 2000.0, Meeus Eq. 25.1
       tjc2  =  tjc**2
       l0    =  4.8950631684d0 + 628.331966786d0 * tjc + 5.291838d-6 *tjc2  ! Mean longitude of the Sun, Meeus Eq. 25.2
       
       ! Mean anomaly of the Sun, Meeus Eq. 25.3 -> p.144 (nutation):
       mas   =  6.240035881d0 + 628.301956024d0 * tjc - 2.79776d-6 * tjc2 - 5.817764d-8 * tjc**3
       
       ! Sun's equation of the centre, Meeus p.164:
       sec   =  (3.3416109d-2 - 8.40725d-5*tjc - 2.443d-7*tjc2)*sin(mas) + (3.489437d-4-1.763d-6*tjc)*sin(2*mas) &
            + 5.044d-6*sin(3*mas)
       odot  =  rev(l0 + sec)  ! Sun's true longitude
       
       ! Meeus, p.151:
       k     =  9.93650849745d-5  ! Constant of aberration, in radians (kappa in Meeus)
       ee    =  0.016708634d0 - 4.2037d-5 * tjc   - 1.267d-7 * tjc2
       pii   =  1.79659568d0  + 3.001024d-2 * tjc + 8.0285d-6 * tjc2
       
       
       ! If eps0 is provided, use it - otherwise compute it:
       if(present(eps0)) then
          eps = eps0
       else
          call nutation(tjc/10.d0, dumdbl1, eps, dumdbl2)  ! Note: eps is actually eps0
       end if
       
    end if
    jdold = jd
    
    
    ! Note that eps in Meeus is actually eps0, the *mean* obliquity of the ecliptic
    ! Meeus, Eq. 23.3:
    dra =      k * ( -(cos(ra) * cos(odot) * cos(eps)  +  sin(ra) * sin(odot))  / cos(dec)  &
         +    ee *    (cos(ra) * cos(pii)  * cos(eps)  +  sin(ra) * sin(pii))   / cos(dec) )
    
    ddec =      k * ( -(cos(odot) * cos(eps) * (tan(eps) * cos(dec)  -  sin(ra) * sin(dec))  +  cos(ra) * sin(dec) * sin(odot))  &
         +     ee *    (cos(pii)  * cos(eps) * (tan(eps) * cos(dec)  -  sin(ra) * sin(dec))  +  cos(ra) * sin(dec) * sin(pii)) )
    
  end subroutine aberration_eq
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert coordinates to the FK5 system
  !!
  !! \param  t  Dynamical time in Julian Millennia since 2000.0
  !! \retval l  Longitude (rad, I/O)
  !! \retval b  Latitude (rad, I/O)
  
  subroutine fk5(t, l,b)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: t
    real(double), intent(inout) :: l,b
    real(double) :: tt,l2,dl,db
    
    tt = t*10  ! Julian Centuries
    
    l2 = l - 0.02438225d0*tt - 5.41052d-6*tt*tt
    dl = -4.379322d-7 + 1.89853d-7*(cos(l2) + sin(l2))*tan(b)
    db = 1.89853d-7*(cos(l2) - sin(l2))
    
    l  = l + dl
    b  = b + db
    
  end subroutine fk5
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert (geocentric) spherical ecliptical coordinates l,b (and eps) to spherical equatorial coordinates RA, Dec
  !!
  !! \param  l    Longitude (rad)
  !! \param  b    Latitude (rad)
  !! \param  eps  Obliquity of the ecliptic
  !!
  !! \retval ra   Right ascension (rad)
  !! \retval dec  Declination (rad)
  
  subroutine ecl_2_eq(l,b,eps, ra,dec)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev
    
    implicit none
    real(double), intent(in) :: l,b,eps
    real(double), intent(out) :: ra,dec
    
    ra  = rev(atan2( sin(l) * cos(eps)  -  tan(b) * sin(eps),  cos(l) ))
    dec =      asin( sin(b) * cos(eps)  +  cos(b) * sin(eps) * sin(l) )
    
  end subroutine ecl_2_eq
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert (geocentric) spherical equatorial coordinates RA, Dec (and eps) to spherical ecliptical coordinates l,b
  !!
  !! \param  ra   Right ascension (rad)
  !! \param  dec  Declination (rad)
  !! \param  eps  Obliquity of the ecliptic
  !!
  !! \retval l    Longitude (rad)
  !! \retval b    Latitude (rad)
  
  subroutine eq_2_ecl(ra,dec,eps, l,b)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev
    
    implicit none
    real(double), intent(in) :: ra,dec,eps
    real(double), intent(out) :: l,b
    
    l = rev(atan2( sin(ra)  * cos(eps) + tan(dec) * sin(eps),  cos(ra) ))
    b =      asin( sin(dec) * cos(eps) - cos(dec) * sin(eps) * sin(ra) )
    
  end subroutine eq_2_ecl
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert spherical equatorial coordinates (RA, dec, agst) to spherical horizontal coordinates (hh, az, alt)
  !!
  !! \param  ra    Right ascension
  !! \param  dec   Declination
  !! \param  agst  Greenwich siderial time
  !!
  !! \retval hh    Local hour angle
  !! \retval az    Azimuth
  !! \retval alt   Altitude
  
  subroutine eq2horiz(ra,dec,agst, hh,az,alt)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev,rev2
    use TheSky_local, only: lat0,lon0
    
    implicit none
    real(double), intent(in) :: ra,dec,agst
    real(double), intent(out) :: hh,az,alt
    real(double) :: sinhh,coshh, sinlat0,coslat0, sindec,cosdec,tandec
    
    hh  = rev( agst + lon0 - ra )                                 ! Local Hour Angle (agst since ra is also corrected for nutation?)
    
    ! Some preparation, saves ~29%:
    sinhh   = sin(hh)
    coshh   = cos(hh)
    sinlat0 = sin(lat0)
    coslat0 = sqrt(1.d0-sinlat0**2)   ! Cosine of a latitude is always positive
    sindec  = sin(dec)
    cosdec  = sqrt(1.d0-sindec**2)    ! Cosine of a declination is always positive
    tandec  = sindec/cosdec
    
    az  = rev( atan2( sinhh,    coshh  * sinlat0 - tandec * coslat0 ))  ! Azimuth
    alt = rev2( asin( sinlat0 * sindec + coslat0 * cosdec * coshh ))    ! Altitude
    
  end subroutine eq2horiz
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert spherical horizontal coordinates (az, alt, agst) to spherical equatorial coordinates (hh, RA, dec)
  !!
  !! \param  az    Azimuth
  !! \param  alt   Altitude
  !! \param  agst  Greenwich siderial time
  !!
  !! \retval hh    Local hour angle
  !! \retval ra    Right ascension
  !! \retval dec   Declination
  
  subroutine horiz2eq(az,alt,agst, hh,ra,dec)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev,rev2
    use TheSky_local, only: lat0,lon0
    
    implicit none
    real(double), intent(in) :: az,alt,agst
    real(double), intent(out) :: ra,dec,hh
    
    hh  = rev(  atan2( sin(az),    cos(az)  * sin(lat0) + tan(alt) * cos(lat0) ))   ! Local Hour Angle
    dec = rev2( asin(  sin(lat0) * sin(alt) - cos(lat0) * cos(alt) * cos(az)   ))   ! Declination
    ra  = rev( agst + lon0 - hh )                                                   ! Right ascension
    
  end subroutine horiz2eq
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert spherical equatorial coordinates (RA, dec) to spherical galactic coordinates (l,b), for J2000.0!!!
  !!
  !! \param  ra    Right ascension
  !! \param  dec   Declination
  !!
  !! \retval l     Longitude
  !! \retval b     Latitude
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.13
  
  subroutine eq2gal(ra,dec, l,b)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev
    
    implicit none
    real(double), intent(in) :: ra,dec
    real(double), intent(out) :: l,b
    real(double) :: ra0,dec0,l0
    
    ra0  = 3.36603472d0      ! RA of GNP, in rad, J2000.0 (12h51m26.3s ~ 192.8596deg)
    dec0 = 0.473478737d0     ! Decl. of gal. NP in rad, J2000.0 (27d07'42"=27.1283deg)
    l0   = 5.287194d0        ! J2000.0?, = 302.9339deg (the 303deg in Meeus, p.94) = l0+180deg
    
    l = rev( l0 - atan2( sin(ra0-ra), cos(ra0-ra) * sin(dec0)  -  tan(dec) * cos(dec0) ))
    b =            asin( sin(dec)                 * sin(dec0)  +  cos(dec) * cos(dec0) * cos(ra0-ra) )
    
  end subroutine eq2gal
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert spherical galactic coordinates (l,b) to spherical equatorial coordinates (RA, dec), for J2000.0!!!
  !!
  !! \param  l     Longitude
  !! \param  b     Latitude
  !!
  !! \retval ra    Right ascension
  !! \retval dec   Declination
  
  subroutine gal2eq(l,b, ra,dec)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev
    
    real(double), intent(in) :: l,b
    real(double), intent(out) :: ra,dec
    real(double) :: ra0,dec0,l0
    
    ra0  = 0.22444207d0     ! RA of GNP - 12h, in rad, J2000.0 (12h51m26.3s - 12h ~ 12.8596deg)
    dec0 = 0.473478737d0    ! Decl. of gal. NP in rad, J2000.0 (from Sterrengids?: 27d07'42"=27.1283deg)
    l0   = 2.145601d0       ! J2000.0?, = 123.9339deg (the 123deg in Meeus, p.94)
    
    ra  = rev(atan2( sin(l-l0), cos(l-l0) * sin(dec0)  -  tan(b) * cos(dec0)) + ra0)
    dec =      asin( sin(b)   * sin(dec0)              +  cos(b) * cos(dec0) * cos(l-l0))
    
  end subroutine gal2eq
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert spherical ecliptical coordinates from the geocentric to the topocentric system
  !!
  !! \param  gcl   Geocentric longitude
  !! \param  gcb   Geocentric latitude
  !! \param  gcr   Geocentric distance
  !! \param  gcs   Geocentric semi-diameter
  !!
  !! \param  eps  Obliquity of the ecliptic
  !! \param  lst  Local siderial time
  !!
  !! \retval tcl   Topocentric longitude
  !! \retval tcb   Topocentric latitude
  !! \retval tcs   Topocentric semi-diameter
  !!
  !! \note lat0 and height are provided by the module TheSky_local
  !!
  !! \see  Meeus, Astronomical Algorithms, 1998, Ch. 11 and 40
  
  subroutine geoc2topoc_ecl(gcl,gcb,gcr,gcs,eps,lst, tcl,tcb,tcs)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev
    use TheSky_local, only: lat0, height
    
    implicit none
    real(double), intent(in) :: gcl,gcb,gcr,gcs,eps,lst
    real(double), intent(out) :: tcl,tcb,tcs
    real(double) :: ba,re,u,rs,rc,shp,n
    
    ! Meeus, Ch.11, p.82:
    ba = 0.99664710d0   ! b/a = 1-f
    re = 6378140.d0     ! Earth rad in m
    
    u  = atan(ba*tan(lat0))
    rs = ba*sin(u) + height/re*sin(lat0)
    rc = cos(u)    + height/re*cos(lat0)
    
    
    shp = sin(4.26345d-5)/gcr  ! Sine of the horizontal parallax, Meeus, Eq. 40.1
    
    ! Meeus, Ch.40, p.282:
    n  = cos(gcl)*cos(gcb) - rc*shp*cos(lst)
    
    tcl = rev( atan2( sin(gcl)*cos(gcb) - shp*(rs*sin(eps) + rc*cos(eps)*sin(lst)) , n ) )  ! Topocentric longitude
    tcb = atan((cos(tcl)*(sin(gcb) - shp*(rs*cos(eps) - rc*sin(eps)*sin(lst))))/n)          ! Topocentric latitude
    tcs = asin(cos(tcl)*cos(tcb)*sin(gcs)/n)                                                ! Topocentric semi-diameter
    
  end subroutine geoc2topoc_ecl
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert geocentric equatorial coordinates to topocentric
  !!
  !! \param  gcra  Geocentric right ascension
  !! \param  gcd   Geocentric declination
  !! \param  gcr   Geocentric distance
  !! \param  gch   Geocentric hour angle
  !!
  !! \retval tcra  Topocentric right ascension
  !! \retval tcd   Topocentric declination
  !!
  !! \see  Meeus, Astronomical Algorithms, 1998, Ch. 11 and 40
  
  subroutine geoc2topoc_eq(gcra,gcd,gcr,gch, tcra,tcd)
    use SUFR_kinds, only: double
    use TheSky_local, only: lat0, height
    
    implicit none
    real(double), intent(in) :: gcra,gcd,gcr,gch
    real(double), intent(out) :: tcra,tcd
    real(double) :: ba,re,u,rs,rc,shp,dra
    
    ! Meeus, Ch.11, p.82:
    ba = 0.99664710d0   ! 1-f
    re = 6378140.d0     ! Earth radius in m
    
    u  = atan(ba*tan(lat0))
    rs = ba*sin(u) + height/re*sin(lat0)
    rc = cos(u) + height/re*cos(lat0)
    
    ! Meeus Ch.40:
    shp   = sin(4.26345d-5)/gcr                                            ! Sine of the horizontal parallax, Meeus, Eq. 40.1
    dra  = atan2( -rc*shp*sin(gch) , cos(gcd)-rc*shp*cos(gch) )            ! Meeus, Eq. 40.2
    tcra = gcra + dra                                                      ! Topocentric right ascension
    tcd  = atan2( (sin(gcd)-rs*shp)*cos(dra) , cos(gcd)-rc*shp*cos(gch) )  ! Topocentric declination - Meeus, Eq. 40.3
    
  end subroutine geoc2topoc_eq
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the atmospheric refraction for a given altitude.  You should add the result to the uncorrected altitude.
  !!
  !! \param   alt      Altitude (rad)
  !! \param   press    Air pressure (hPa; optional)
  !! \param   temp     Air temperature (degrees Celcius; optional)
  !!
  !! \retval  refract  Refraction in altitude (rad).  You should add the result to the uncorrected altitude.
  !!
  !! \see Meeus (1998), Eq. 16.4 ff, based on Samundsson, Sky & Telescope vol.72, p.70 (1986), converted to radians
  
  function refract(alt, press,temp)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pio2, d2r
    
    implicit none
    real(double), intent(in) :: alt
    real(double), intent(in), optional :: press,temp
    real(double) :: refract
    
    if(abs(alt).ge.pio2) then  ! |alt| >= 90 deg; refraction is meaningless
       refract = 0.d0
    else
       refract = 2.97d-4/tan(alt + 3.14d-3/(alt + 8.92d-2))
       if(present(press)) refract = refract * press/1010.d0              ! Correct for pressure
       if(present(temp))  refract = refract * 283.d0/(273.d0 + temp)     ! Correct for temperature
    end if
    
    if(alt+refract.lt.-0.3d0*d2r) refract = 0.d0  ! No refraction if the object is more than 0.3 deg below the horizon 
    !                                               (maximum apparent Moon radius ~ 0.29deg)
    
  end function refract
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the atmospheric refraction of light for a given true altitude.
  !!         This is a wrapper for aref(), which does the opposite (compute refraction for an observed zenithal angle).
  !!         This is an expensive way to go about(!)
  !!
  !! \param  alt0  The true (theoretical, computed) altitude of the object in radians
  !!
  !! \param  h0    The height of the observer above sea level in metres
  !! \param  lat0  The latitude of the observer in radians
  !!
  !! \param  t0    The temperature at the observer in degrees Celcius
  !! \param  p0    The pressure at the observer in millibars
  !! \param  rh    The relative humidity at the observer (%)
  !!
  !! \param  lam   The wavelength of the light at the observer in nanometres
  !! \param  dTdh  The temperature lapse rate dT/dh in Kelvin/metre in the troposphere (only the absolute value is used)
  !!
  !! \param  eps   The desired precision in arcseconds
  !!
  !! \retval atmospheric_refraction  The refraction at the observer in radians
  !!
  !!
  !! \todo  Adapt aref() to compute the integral the other way around for a direct method(?)
  !!
  
  function atmospheric_refraction(alt0, h0,lat0, t0,p0,rh, lam,dTdh, eps)
    use SUFR_kinds, only: double
    use SUFR_constants, only: d2r,r2d, pio2
    
    implicit none
    real(double), intent(in) :: alt0, h0,lat0, t0,p0, rh,lam, dTdh, eps
    real(double) :: atmospheric_refraction, z0,zi,zio,dz, ph,lamm,t0K,rhfr
    
    atmospheric_refraction = 0.d0
    if(abs(alt0).ge.pio2) return  ! No refraction if |alt| >= 90d
    
    ! Convert some variables/units:
    z0   = 90.d0 - alt0*r2d  ! Altitude (rad) -> zenith angle (deg)
    ph   = lat0*r2d          ! Latitude: rad -> deg
    lamm = lam / 1000.d0     ! Wavelength: nanometre -> micrometre
    t0K  = t0 + 273.15d0     ! Temperature: degC -> K
    rhfr = rh/100.d0         ! Relative humidity: % -> fraction
    
    zi  = z0
    zio = huge(zio)
    dz  = huge(dz)
    
    do while(dz.gt.eps)
       atmospheric_refraction = aref(zi, h0,ph, t0K,p0,rhfr, lamm,dTdh, eps)
       zi = z0 - atmospheric_refraction
       dz = abs(zi-zio)
       zio = zi
    end do
    
    atmospheric_refraction = atmospheric_refraction * d2r
    
  end function atmospheric_refraction
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the atmospheric refraction of light for a given observed zenith angle.
  !!         The method is based on N.A.O Technical Notes 59 and 63 and a paper by Standish and Auer 'Astronomical Refraction: 
  !!         Computational Method for all Zenith Angles'.
  !!
  !! \param  z0    The observed zenith distance of the object in degrees
  !!
  !! \param  h0    The height of the observer above sea level in metres
  !! \param  ph    The latitude of the observer in degrees
  !!
  !! \param  t0    The temperature at the observer in Kelvin
  !! \param  p0    The pressure at the observer in millibars
  !! \param  rh    The relative humidity at the observer
  !!
  !! \param  lam   The wavelength of the light at the observer in micrometres
  !! \param  dTdh  The temperature lapse rate dT/dh in Kelvin/metre in the troposphere, the absolute value is used
  !!
  !! \param  eps   The desired precision in arcseconds
  !!
  !! \retval aref  The refraction at the observer in degrees
  !!
  !! \see  Hohenkerk & Sinclair, HMNAO technical note 63 (1985)
  !!
  
  function aref(z0, h0,ph, t0,p0,rh, lam,dTdh,eps)
    use SUFR_kinds, only: double
    implicit none
    real(double), intent(in) :: z0, h0,ph, t0,p0, rh,lam, dTdh, eps
    
    real(double), parameter :: gcr=8314.36d0, md=28.966d0, mw=18.016d0, s=6378120.d0, gamma=18.36d0
    real(double), parameter :: ht=11000.d0, hs=80000.d0, dgr=0.01745329252d0, z2=11.2684d-06
    
    integer :: i,in,is,istart, j,k
    real(double) :: aref, n,n0,nt,nts,ns,a(10), dndr,dndr0,dndrs,dndrt,dndrts, f,f0, fb,fe,ff,fo,fs,ft,fts,gb, h
    real(double) :: pw0,r,r0,refo,refp,reft,rg,rs,rt,  sk0,step,t,t0o,tg,tt,z,z1,zs,zt,zts
    
    ! Always defined:
    z = 0.d0; reft=0.d0
    
    ! Set up parameters defined at the observer for the atmosphere:
    gb = 9.784d0*(1.d0 - 0.0026d0*cos(2.d0*ph*dgr) - 0.00000028d0*h0)
    z1 = (287.604d0 + 1.6288d0/(lam**2) + 0.0136d0/(lam**4)) * (273.15d0/1013.25d0)*1.d-6
    
    a(1) = abs(dTdh)
    a(2) = (gb*md)/gcr
    a(3) = a(2)/a(1)
    a(4) = gamma
    pw0  = rh*(t0/247.1d0)**a(4)
    a(5) = pw0*(1.d0 - mw/md)*a(3)/(a(4)-a(3))
    a(6) = p0 + a(5)
    a(7) = z1*a(6)/t0
    a(8) = ( z1*a(5) + z2*pw0)/t0
    a(9) = (a(3) - 1.d0)*a(1)*a(7)/t0
    a(10) = (a(4) - 1.d0)*a(1)*a(8)/t0
    
    ! At the Observer:
    r0 = s + h0
    call troposphere_model(r0,t0,a,r0,t0o,n0,dndr0)
    sk0 = n0 * r0 * sin(z0*dgr)
    
    
    f0 = refi(r0,n0,dndr0)
    
    ! At the Tropopause in the Troposphere:
    rt = s + ht
    call troposphere_model(r0,t0,a,rt,tt,nt,dndrt)
    zt = asin(sk0/(rt*nt))/dgr
    ft = refi(rt,nt,dndrt)
    
    ! At the Tropopause in the Stratosphere:
    call stratosfeer_model(rt,tt,nt,a(2),rt,nts,dndrts)
    zts = asin(sk0/(rt*nts))/dgr
    fts = refi(rt,nts,dndrts)
    
    ! At the stratosphere limit:
    rs = s + hs
    call stratosfeer_model(rt,tt,nt,a(2),rs,ns,dndrs)
    zs = asin(sk0/(rs*ns))/dgr
    fs = refi(rs,ns,dndrs)
    
    ! Integrate the refraction integral in the troposhere and stratosphere
    ! ie Ref = Ref troposhere + Ref stratopshere
    
    ! Initial step lengths etc:
    refo = -huge(refo)
    is = 16
    do k = 1,2
       istart = 0
       fe = 0.d0
       fo = 0.d0
       
       if(k.eq.1) then
          h = (zt - z0)/dble(is)
          fb = f0
          ff = ft
       else if(k.eq.2) then
          h = (zs - zts )/dble(is)
          fb = fts
          ff = fs
       end if
       
       in = is - 1
       is = is/2
       step = h
       
200    continue
       
       do i = 1,in
          if(i.eq.1.and.k.eq.1) then
             z = z0 + h
             r = r0
          else if(i.eq.1.and.k.eq.2) then
             z = zts + h
             r = rt
          else
             z = z + step
          end if
          
          ! Given the Zenith distance (Z) find R:
          rg = r
          do j = 1,4
             if(k.eq.1) then
                call troposphere_model(r0,t0,a,rg,tg,n,dndr)
             else if(k.eq.2) then
                call stratosfeer_model(rt,tt,nt,a(2),rg,n,dndr)
             end if
             rg = rg - ( (rg*n - sk0/sin(z*dgr))/(n + rg*dndr) )
          end do
          r = rg
          
          ! Find Refractive index and Integrand at R:
          if(k.eq.1) then
             call troposphere_model(r0,t0,a,r,t,n,dndr)
          else if(k.eq.2) then
             call stratosfeer_model(rt,tt,nt,a(2),r,n,dndr)
          end if
          
          f = refi(r,n,dndr)
          
          if(istart.eq.0.and.mod(i,2).eq.0) then
             fe = fe + f
          else
             fo = fo + f
          end if
       end do
       
       ! Evaluate the integrand using Simpson's Rule:
       refp = h*(fb + 4.d0*fo + 2.d0*fe + ff)/3.d0
       
       if(abs(refp-refo).gt.0.5d0*eps/3600.d0) then
          is = 2*is
          in = is
          step = h
          h = h/2.d0
          fe = fe + fo
          fo = 0.d0
          refo = refp
          if(istart.eq.0) istart = 1
          goto 200
       end if
       if(k.eq.1) reft = refp
    end do
    
    aref = reft + refp
    
  end function aref
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  The refraction integrand
  !!
  !! \param r      The current distance from the centre of the Earth in metres
  !! \param n      The refractive index at R
  !! \param dndr   The rate the refractive index is changing at R
  !!
  !! \retval refi  The integrand of the refraction function
  
  function refi(r,n,dndr)
    use SUFR_kinds, only: double
    implicit none
    real(double), intent(in) :: r,n,dndr
    real(double) :: refi
    
    refi = r*dndr/(n + r*dndr)
  end function refi
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Atmospheric model for the troposphere
  !!
  !! \param r0     The height of the observer from the centre of the Earth
  !! \param t0     The temperature at the observer in Kelvin
  !! \param a      Constants defined at the observer
  !! \param r      The current distance from the centre of the Earth in metres
  !!
  !! \retval t     The temperature at R in Kelvin
  !! \retval n     The refractive index at R
  !! \retval dndr  The rate the refractive index is changing at R
  
  subroutine troposphere_model(r0,t0, a,r, t,n,dndr)
    use SUFR_kinds, only: double
    implicit none
    real(double), intent(in) :: r0,t0, a(10),r
    real(double), intent(out) :: t,n,dndr
    real(double) :: tt0, tt01,tt02
    
    t    = t0 - a(1)*(r-r0)
    tt0  = t/t0
    tt01 = tt0**(a(3)-2.d0)
    tt02 = tt0**(a(4)-2.d0)
    
    n    = 1.d0 + ( a(7)*tt01 - a(8)*tt02 )*tt0
    dndr = -a(9)*tt01 + a(10)*tt02
    
  end subroutine troposphere_model
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief Atmospheric model for the stratosphere
  !!
  !! \param  rt    The height of the tropopause from the centre of the Earth in metres
  !! \param  tt    The temperature at the tropopause in Kelvin
  !! \param  nt    The refractive index at the tropopause
  !! \param  a     Constant of the atmospheric model = G*MD/R
  !! \param  r     The current distance from the centre of the Earth in metres
  !!
  !! \retval n     The refractive index at R
  !! \retval dndr  The rate the refractive index is changing at R
  
  subroutine stratosfeer_model(rt,tt,nt, a,r, n,dndr)
    use SUFR_kinds, only: double
    implicit none
    real(double), intent(in) :: rt,tt,nt, a,r
    real(double), intent(out) :: n,dndr
    real(double) :: b
    
    b    = a/tt
    n    = 1.d0 + (nt - 1.d0)*exp(-b*(r-rt))
    dndr = -b*(nt-1.d0)*exp(-b*(r-rt))
    
  end subroutine stratosfeer_model
  !*********************************************************************************************************************************
  
end module TheSky_coordinates
!***********************************************************************************************************************************

