!> \file comets.f90  Contains modules and procedures to compute comet positions and properties for libTheSky


!  Copyright (c) 2002-2013  Marc van der Sluys - marc.vandersluys.nl
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
!> \brief Procedures for comets

module TheSky_comets
  implicit none
  save
  
contains
  
  !*********************************************************************************************************************************
  !> \brief Calculate heliocentric xyz coordinates for comet com at t1 (in J.mill DT).  Usually called by cometgc() below.
  !! 
  !! \param t1   Time in Julian millennia DT
  !! \param com  Comet ID
  !! 
  !! \retval x   Heliocentric x coordinate
  !! \retval y   Heliocentric y coordinate
  !! \retval z   Heliocentric z coordinate
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch. 30, 33-35.  Equation and page numbers refer to this book.
  !!
  !! \todo  Use "hyperbolic method" for 0.98 < e < 1 as well? - see CHECK
  
  subroutine cometxyz(t1,com,x,y,z)
    use SUFR_kinds, only: double
    use SUFR_constants, only: nlpname,enpname
    
    use TheSky_nutation, only: nutation
    use TheSky_cometdata, only: cometNames, cometElems, comepoche
    
    implicit none
    real(double), intent(in) :: t1
    integer, intent(in) :: com
    real(double), intent(out) :: x,y,z
    
    integer :: i,j,j1,j2
    real(double) :: k,jj,del,  tp,q,a,e,o1,o2,in,nu,r,t,jde,te,  m,n,ee,ee1,de,  qq,g,s,s0,s1,w,q2,q3,dq3
    real(double) :: eps,dpsi,eps0,deps,  ff,gg,hh,pp,qqq,rr,a1,a2,b1,b2,c1,c2
    
    jde = t1*365250.d0 + 2451545.d0       ! t1 in Julian Millennia after 2000.0 in dynamical time
    
    k   = 0.01720209895d0  ! Gaussian gravitational constant
    j1  = 0
    j2  = 0
    del = 1.d-10           ! Convergence criterion
    
    comepoche   = cometElems(com,1)      ! J2000.0
    q           = cometElems(com,2)      ! Perihelion distance (AU?)
    e           = cometElems(com,3)      ! Eccentricity
    in          = cometElems(com,4)      ! Inclination
    o1          = cometElems(com,5)      ! Argument of perihelion (lowercase omega)
    o2          = cometElems(com,6)      ! Longitude of ascending node (uppercase omega)
    tp          = cometElems(com,7)      ! Perihelion date (JD)
    nlpname(11) = trim(cometNames(com))  ! Warning: cometNames is much longer than nlpname
    enpname(11) = trim(cometNames(com))  ! Warning: cometNames is much longer than enpname
    
    !write(0,'(I9,9F15.5)') com, o1,o2,in,q,e,tp,comepoche
    !write(6,'(I9,2x,A50,9F15.5)') com, trim(cometNames(com)), q,e,in*r2d, o1*r2d,o2*r2d, tp,comepoche
    
    
    
    te = (comepoche - 2451545.d0)/365250.d0
    t = jde - tp  ! In days since perihelion
    nu = 0.d0     ! Avoid 'may be used uninitialized' warnings
    
    !call precess_orb(comepoche,jde,in,o1,o2)
    
    
    if(e.lt.1.0d0) then
       a = q/(1.d0 - e)   ! Semi-major axis,        Eq. 33.6a
       n = k*a**(-1.5d0)  ! Mean motion (rad/day),  Eq. 33.6b
       m = n*t            ! Mean anomaly - M = 0 at perihelion (t=0)
       ee = m             ! Excentric anomaly - initial value
       
       ! Solve Kepler's equation for mildly elliptic orbits ("first method", p.196):
       if(e.lt.0.5d0) then
          do i=1,nint(1e7)
             ee1 = m + e*sin(ee)  ! Eq. 30.6
             
             if(abs(ee-ee1).lt.del) exit
             ee = ee1
          end do
          ee = ee1
          
          
          !  Solve Kepler's equation for strongly elliptic orbits ("second method", p.199):
       else  ! CHECK or 0.5 - 0.98?
          do i=1,nint(1e7)
             de = (m + e*sin(ee) - ee) / (1.d0 - e*cos(ee))  ! Eq. 30.7
             ee = ee + de
             if(abs(de).lt.del) exit
          end do
       end if
       
       if(i.ge.nint(1e7)) write(0,'(A,I0,A,F8.5)') '  cometxyz():  WARNING:  Kepler solution did not converge for comet '// &
            trim(cometNames(com))//' (',com,'), with e =',e
       
       nu = 2.d0 * atan( sqrt( (1.d0+e)/(1.d0-e) ) * tan(ee/2.d0) )  ! True anomaly, elliptic orbits, Eq. 30.1
       
       
    else  ! e >= 1: Parabolic or hyperbolic orbits:
       
       ! Get a starting value for s, for parabolic or hyperbolic orbits:
       w = 3.d0*k/(q*sqrt(2.d0*q))*t  ! Eq. 34.1
       s = w/3.d0  !0.d0  Starting value for auxillary variable s = tan(nu/2)
       do i=1,1000
          s1 = (2.d0*s**3+w)/(3.d0*(s**2+1.d0))  ! Eq. 34.4
          if(abs(s-s1).lt.del) exit
          s = s1
       end do  ! i
       s = s1
       
       
       ! Parabolic orbit:
       if(e.eq.1.d0) then
          !write(6,*)'iter, delta: ',i,s-s0
          nu = 2.d0*atan(s)    ! True anomaly, parabolic orbits, Eq. 34.2
          
          
       ! Hyperbolic orbit:
       else  ! e > 1   CHECK - use this too for 0.98 < e < 1 ?
          qq = k/(2.d0*q)*sqrt((1.d0+e)/q)  ! Eq. "35.0a" - Q
          g  = (1.d0-e)/(1.d0+e)            ! Eq. "35.0b" - gamma
          q2 = qq*t                         ! Eq.  35.1a
          
          iloop: do i=1,nint(1e7)           ! Program p.246, loop lines 40-64
             s0 = s
             q3 = q2 + 2.d0*g*s**3/3.d0     ! Eq. 35.1b, but version from program on next page (246), line 42
             
             j1loop: do j=3,1000            ! Eq. 35.1c, d, ...  -  program p.246, loop lines 44-56
                jj = dble(j)
                dq3 = (-1)**(j-1) * g**(j-2) * ((jj-1.d0) - jj*g) * s**(2*j-1) / (2.d0*jj-1.d0)
                q3 = q3 + dq3               ! Program p.246, loop line 52
                if(abs(dq3).lt.del) exit  j1loop  ! j
                
                if(j.eq.1000) j1 = j1+1
                if(j1.eq.100) then
                   write(6,*) "COMETS: j1-iteration didn't converge"
                   return
                end if
             end do  j1loop  ! j
             
             j2loop: do j=1,1000                       ! Program p.246, loop lines 60-62
                s1 = s
                s = (2.d0/3.d0*s**3+q3) / (s**2+1.d0)  ! Program p.246, loop line 60
                if(abs(s-s1).lt.del) exit  j2loop  ! j
                
                if(j.eq.1000) j2 = j2+1
                if(j2.eq.100) then
                   write(6,*) "COMETS: j2-iteration didn't converge",abs(s-s1)
                   return
                end if
             end do  j2loop  ! j
             
             if(abs(s-s0).lt.del) then
                nu = 2.d0*atan(s)    ! True anomaly, hyperbolic orbits, Eq. "35.2" / 34.2
                exit  iloop  ! i
             end if
             
             if(i.ge.nint(1e7)) then
                write(6,*) "COMETS: i-Iteration didn't converge: ",abs(s-s1)
                return
             end if
             
          end do  iloop  ! i
       end if  ! Hyperbolic orbit
       
    end if  ! e >= 1  ! Parabolic or hyperbolic orbit
    
    !r = a * (1.d0 - e*cos(ee))  ! Eq. 30.2
    r = q * (1.d0+e) / (1.d0 + e*cos(nu))  ! Eq. 30.4
    
    call nutation(te,dpsi,eps0,deps)
    eps = eps0 !+ deps
    
    ! Eq. 33.7:
    ff  =  cos(o2)           ! F
    gg  =  sin(o2)*cos(eps)  ! G
    hh  =  sin(o2)*sin(eps)  ! H
    pp  = -sin(o2)*cos(in)                              ! P
    qqq =  cos(o2)*cos(in)*cos(eps) - sin(in)*sin(eps)  ! Q
    rr  =  cos(o2)*cos(in)*sin(eps) + sin(in)*cos(eps)  ! R
    
    ! Eq. 33.8:
    a1 = atan2(ff, pp)          ! A
    a2 = sqrt(ff*ff + pp*pp)    ! a
    b1 = atan2(gg, qqq)         ! B
    b2 = sqrt(gg*gg + qqq*qqq)  ! b
    c1 = atan2(hh, rr)          ! C
    c2 = sqrt(hh*hh + rr*rr)    ! c
    
    ! Eq. 33.9:
    x = r * a2 * sin(a1+o1+nu)
    y = r * b2 * sin(b1+o1+nu)
    z = r * c2 * sin(c1+o1+nu)
    
  end subroutine cometxyz
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calc geocentric lbr coordinates for comet com at t1 (in Julian millennia Dynamical Time)
  !! 
  !! \note
  !! This a computational detour to allow correction for aberration, 
  !!   fk5 and nutation in planet_position (and keep overview there)
  !!
  !! \param  t    Apparent time (taking light time into account) in Julian millennia Dynamical Time
  !! \param  t0   True time in Julian millennia DT
  !! \param  com  Comet ID
  !!
  !! \retval r    Apparent heliocentric distance of comet
  !! \retval l    Apparent geocentric ecliptic longitude of comet
  !! \retval b    Apparent geocentric ecliptic latitude of comet
  !! \retval d    Apparent geocentric distance of comet
  !! 
  
  subroutine cometgc(t,t0, com, r,l,b,d)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi
    use SUFR_angles, only: rev
    
    use TheSky_coordinates, only: ecl_spher_2_eq_rect, fk5, precess_xyz, eq_2_ecl
    use TheSky_nutation, only: nutation
    use TheSky_vsop, only: vsop_lbr
    use TheSky_cometdata, only: comepoche
    
    implicit none
    real(double), intent(in) :: t,t0
    integer, intent(in) :: com
    real(double),intent(out) :: l,b,r,d
    
    real(double) :: l0,b0,r0,x0,y0,z0,x,y,z
    real(double) :: ra,dec,dpsi,eps0,deps,eps,jd1,jd2
    
    
    call vsop_lbr(t0,3, l0,b0,r0)  ! Earth
    call fk5(t, l0,b0)
    call nutation(t,dpsi,eps0,deps)
    eps = eps0 + deps
    !call calcsunxyz(t,l0,b0,r0,x0,y0,z0)
    call cometxyz(t,com, x,y,z)                          ! Heliocentric equatorial rectangular coordinates
    call ecl_spher_2_eq_rect(rev(l0+pi),-b0,r0, eps0, x0,y0,z0)  ! l0+pi,-b0,r0 is geocentric SPHERICAL ECLIPTICAL pos. of Sun,
    !                                                              convert to geocentric RECTANGULAR EQUATORIAL position of Sun
    jd1 = t*365250.d0 + 2451545.d0                       ! From equinox of date...
    jd2 = comepoche                                      !  ... to equinox of elements
    call precess_xyz(jd1,jd2, x0,y0,z0)                  ! Precess geocentric position of the Sun
    
    
    r = sqrt(x*x + y*y + z*z)                            ! Heliocentric distance comet
    
    ! Heliocentric -> geocentric position of the comet:
    x = x + x0
    y = y + y0
    z = z + z0
    
    d   = sqrt(x*x+y*y+z*z)  ! 'Delta' - geocentric distance
    ra  = atan2(y,x)         ! Right ascension
    dec = asin(z/d)          ! Declination
    
    !call nutation(t,dpsi,eps0,deps)
    !eps = eps0 + deps
    call eq_2_ecl(ra,dec,eps, l,b)
    
  end subroutine cometgc
  !*********************************************************************************************************************************
  
  
  
end module TheSky_comets
!***********************************************************************************************************************************

