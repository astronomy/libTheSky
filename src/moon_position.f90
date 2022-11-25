!> \file moon_position.f90  Core procedures that calculate the position and magnitude of the Moon for libTheSky


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
!> \brief Procedures for the Moon

module TheSky_moon
  implicit none
  save
  
contains
  
  !*********************************************************************************************************************************
  !> \brief  Calculate the apparent geocentric ecliptical position of the Moon, using the Lunar Solution ELP 2000-82B
  !! 
  !! \param  tjj  Time for calculation in Julian millenia after 2000.0
  !!
  !! \param ll   Apparent geocentric ecliptical longitude (output)
  !! \param bb   Apparent geocentric ecliptical latitude (output)
  !! \param rr   Apparent geocentric distance (output)
  !!
  !! \note
  !! - This is supposed to be the ELP2000-85 version, with the corrections from the 1998 paper (mean arguments up to t^4,
  !!     etc,).  However, the accuracy seems to be less than promised for historical calculations (0.1° rather than 0.01°
  !!     for CE 0).  The subroutine elp_mpp02_lbr() below is supposed to give better results (but then again, so is this
  !!     one).
  !! - Differences compared to ELP-MPP02 (ELP82b minus ELP-MPP02):
  !!   - longitude: +0.00001° in CE 1975, ~+0.11° in CE 0/4000 and ~+0.65° in 3000 BCE (systematicly drifting to positive values)
  !!                - replacing w(1,2)  = -5.8883 with -6.8084 removes the systematic drift to +0.65° @3000 BCE; (data.f90 / readmoondata())
  !!                  the drift now oscillates around +-0.012°;
  !!   - latitude:  ~0.00000° in CE 1990, ~+/-0.01° in CE 0/4000 and +/- ~0.06° (now ~0.006°) in 3000 BCE
  !!   - distance:  +/-0.03 km in CE 2000, +/- ~30 km in CE 0/4000 and +/- ~300 km (now ~20km) in 3000 BCE
  !!
  !! \see
  !!  - Chapront-Touzé & Chapront, A&A, 124, 50 (1983)
  !!  - Chapront-Touzé & Chapront, A&A, 190, 342 (1988)
  !!  - ftp://cdsarc.u-strasbg.fr/pub/cats/VI/79/
  
  subroutine elp82b_lbr(tjj, ll,bb,rr)
    use SUFR_kinds, only: double, dbl
    use SUFR_constants, only: as2r, au, km  !,  r2d
    use SUFR_angles, only: rev
    
    use TheSky_moondata, only: t, nterm,nrang, pc1,pc2,pc3, per1,per2,per3, w, ath,a0
    
    implicit none
    real(double), intent(in) :: tjj
    real(double), intent(out) :: ll,bb,rr
    
    integer :: iv,itab,nt,k,j
    real(double) :: r(3),x,y, pa,t4,t8  !, lon,lat,dist
    
    
    r = 0.0_dbl
    x = 0.0_dbl
    y = 0.0_dbl
    
    t(1) = tjj*10.0_dbl  ! In Julian centuries since 2000.0
    t(2) = t(1)**2       ! t^2
    t(3) = t(2)*t(1)     ! t^3
    t(4) = t(2)**2       ! t^4
    
    t4 = t(4)   ! t^4
    t8 = t4**2  ! t^8
    
    
    do iv=1,3  ! 3 variables (lon, lat, dist)
       r(iv) = 0.0_dbl
       
       do itab = 1,12  ! 12 tables: ELP1,4,7,10,13,16,19,22,25,28,31,34 for lon
          do nt = 1,nterm(iv,itab)
             
             select case(iv)  ! variable
                
             case(1)  ! lon
                if(itab.eq.1) then
                   x = pc1(1,nt)
                   y = pc1(2,nt)
                   do k = 1,4
                      y = y + pc1(k+2,nt)*t(k)
                   end do
                else
                   j = nrang(1,itab-1) + nt
                   x = per1(1,j)
                   y = per1(2,j) + per1(3,j)*t(1)
                end if
                
             case(2)  ! lat
                if(itab.eq.1) then
                   x = pc2(1,nt)
                   y = pc2(2,nt)
                   do k = 1,4
                      y = y + pc2(k+2,nt)*t(k)
                   end do
                else
                   j = nrang(2,itab-1) + nt
                   x = per2(1,j)
                   y = per2(2,j) + per2(3,j)*t(1)
                end if
                
             case(3)  ! dist
                if(itab.eq.1) then
                   x = pc3(1,nt)
                   y = pc3(2,nt)
                   do k = 1,4
                      y = y + pc3(k+2,nt)*t(k)
                   end do
                else
                   j = nrang(3,itab-1) + nt
                   x = per3(1,j)
                   y = per3(2,j) + per3(3,j)*t(1)
                end if
             end select
             
             if(itab.eq.3.or.itab.eq.5 .or. itab.eq.7.or.itab.eq.9) x = x*t(1)
             if(itab.eq.12) x = x*t(2)
             r(iv) = r(iv) + x*sin(y)
          end do  ! nt
       end do  ! itab
    end do  ! iv
    
    
    ! Change of units:
    r(1) = r(1) * as2r + w(1,0) + w(1,1)*t(1) + w(1,2)*t(2) + w(1,3)*t(3) + w(1,4)*t(4)  ! Add mean longitude (see ELP PS-file, p.3)
    r(2) = r(2) * as2r
    r(3) = r(3) * a0 / ath
    
    !lon = r(1)
    !lat = r(2)
    !dist = r(3)
    !write(*,'(A, F15.9,2F14.7,9F14.5)') 'ELP82b:    ', t(1), rev(lon)*r2d, lat*r2d, dist, &
    !     rev(w(1,0) + w(1,1)*t(1) + w(1,2)*t(2) + w(1,3)*t(3) + w(1,4)*t(4))*r2d, &
    !     rev(w(1,0))*r2d, rev(w(1,1)*t(1))*r2d, rev(w(1,2)*t(2))*r2d, rev(w(1,3)*t(3))*r2d, rev(w(1,4)*t(4))*r2d
    
    ! Precess from J2000 to EoD, see ELP PS-file, p12;  Laskar 1986 - note: longitude only!:
    pa =   2.438174835301452d-2  * t(1)     + 5.391128133938040d-6  * t(2)     + 3.733065344543427d-10 * t(3)    &
         - 1.140766591650738d-10 * t4       - 8.753311012432670d-14 * t4*t(1)  + 8.460483549042511d-16 * t4*t(2) &
         + 6.348635154129373d-18 * t4*t(3)  + 1.175188363009515E-20 * t8       - 2.307228308400282d-22 * t8*t(1) &
         - 4.198486478408581d-25 * t8*t(2)
    
    r(1) = r(1) + pa
    
    ll = rev(r(1))
    bb = r(2)
    rr = r(3)/au*km
    
    ! The original code converts to rectangular coordinates and applies (more?) precession.
    ! Since I want spherical coordinates, I don't do that here.  This gives similar (though not identical) results.
    
    !write(98,*) ll,bb,rr
    
  end subroutine elp82b_lbr
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the spherical lunar coordinates using the ELP2000/MPP02 lunar theory in the dynamical mean ecliptic and
  !!           equinox of J2000.
  !!
  !! \param jd    Julian day to compute Moon position for
  !! \param mode  Index of the corrections to the constants: 0-Fit to LLR observations, 1-Fit to DE405 1950-2060 (historical)
  !!
  !! \param  lon  Ecliptic longitude (rad) (output)
  !! \param  lat  Ecliptic latitude (rad) (output)
  !! \param  rad  Distance (AU) (output)
  !!
  !! \note  See Lunar Solution ELP 2000/MPP02, Chapront & Francou (2002): ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/
  
  subroutine elp_mpp02_lbr(jd, mode, lon,lat,rad)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000, au, km  !, r2d
    use SUFR_angles, only: rev
    use TheSky_coordinates, only: precess_ecl
    
    implicit none
    real(double), intent(in) :: jd
    integer, intent(in) :: mode
    real(double), intent(out) :: lon,lat,rad
    integer :: ierr
    real(double) :: xyz(3),vxyz(3)
    
    call elp_mpp02_xyz(jd, mode, xyz,vxyz, ierr)
    
    ! Compute ecliptic l,b,r:
    rad = sqrt(sum(xyz**2))
    lon = atan2(xyz(2),xyz(1))
    lat = asin(xyz(3)/rad)
    
    call precess_ecl(jd2000,jd, lon,lat)
    
    rad = rad*km/au  ! km -> AU
    
  end subroutine elp_mpp02_lbr
  !*********************************************************************************************************************************
  
  
  !***************************************************************************************************
  !> \brief  Compute the rectangular lunar coordinates using the ELP/MPP02 lunar theory in the dynamical mean ecliptic and equinox of J2000.
  !!
  !! \param jd    Julian day to compute Moon position for
  !! \param mode  Index of the corrections to the constants: 0-Fit to LLR observations, 1-Fit to DE405 1950-2060 (historical)
  !! 
  !! \param xyz   Geocentric rectangular coordinates: (output)
  !!               - xyz(1) : Position X (km)
  !!               - xyz(2) : Position Y (km)
  !!               - xyz(3) : Position Z (km)
  !! \param vxyz  Geocentric rectangular velocities: (output)
  !!               - vxyz(1) : Velocity X' (km/day)
  !!               - vxyz(2) : Velocity Y' (km/day)
  !!               - vxyz(3) : Velocity Z' (km/day)
  !! \param ierr  File error index - ierr=0: no error, ierr=1: file error (output)
  !!
  !! \note
  !!  - The subroutine elp_mpp02() uses two modules:
  !!    - elp_mpp02_constants:  Constants of the solution ELP/MPP02 (input),
  !!    - elp_mpp02_series:     Series of the solution ELP/MPP02 (input).
  !!
  !!  - The nominal values of some constants have to be corrected.  There are two sets of corrections, which can be selected
  !!    using the parameter 'mode' (used in elp_mpp02_initialise()).
  !!    - mode=0, the constants are fitted to LLR observations provided from 1970 to 2001; it is the default value;
  !!    - mode=1, the constants are fitted to DE405 ephemeris over one century (1950-2060); the lunar angles W1, W2, W3
  !!              receive also additive corrections to the secular coefficients ('historical mode').
  !!    When the mode is changed, the constants will be reinitialised and the data file reread.
  !!
  !!  - Solutions (discussed) in the paper:
  !!    - ELP (original):
  !!      - ELP2000-82: using VSOP82 (1983)
  !!      - ELP2000-85: new mean lunar arguments, higher truncation level, longer time range (1988)
  !!      - ELP2000-82B, here called "ELP": ELP2000-82, using mean lunar arguments from ELP2000-85 (19??)
  !!    - ELP/MPP01:  using latest planetary perturbations from MPP01 and VSOP2000, but simpler than MPP01
  !!    - ELP/MPP02:  ELP/MPP01, but for some arguments back to ELP + different selection of perturbations + lower truncation.  Good fit with DE 405 in [1950,2060]
  !!    - ELP/MPP02*: improved secular arguments, better long-term comparison to DE 405/406 [-3000,2500]
  !!    - ELP/MPP02(LLR): ELP/MPP02(*?), optimised for lunar ranging since 1970
  !!    - ELPa: ELP + few Poisson terms (tested in the current study only?)
  !!    - ELPa*: ELPa + better secular arguments (as in ELP/MPP02*)
  !!  - It is not entirely clear which version is given below, but we can hope it is ELP/MPP02*.  However, the subroutine
  !!      elp82b_lbr() above is known to underperform (by a factor of 10) in accuracy.
  
  subroutine elp_mpp02_xyz(jd, mode, xyz,vxyz, ierr)
    use SUFR_kinds, only: double
    use SUFR_constants, only: r2as, jd2000  !, r2d!,as2r
    use SUFR_system, only: quit_program_error
    use SUFR_angles, only: rev
    use TheSky_data, only: elp_mpp02_initialise_and_read_files
    
    use TheSky_elp_mpp02_series, only: cmpb,fmpb,nmpb,   cper,fper,nper
    use TheSky_elp_mpp02_constants, only: w, p1,p2,p3,p4,p5, q1,q2,q3,q4,q5
    
    implicit none
    integer, intent(in) :: mode
    real(double), intent(in) :: jd
    integer, intent(out) :: ierr
    real(double), intent(out) :: xyz(3),vxyz(3)
    
    
    real(double), parameter :: a405=384747.9613701725d0, aelp=384747.980674318d0, sc=36525.d0  ! Moon mean distance for DE405 und ELP; Julian century in days
    
    integer :: it,iLine,iVar, k
    real(double) :: rjd, t(-1:4),v(6)  !, lon,lat,dist
    real(double) :: cbeta,clamb,cw, ppw,ppw2,ppwqpw,ppwra,pw,pw2,pwqw,pwra, qpw,qpw2,qpwra,qw,qw2,qwra
    real(double) :: ra,rap,sbeta,slamb,sw, x,x1,x2,x3, xp,xp1,xp2,xp3, y,yp
    
    ! Initialise data and read files if needed:
    call elp_mpp02_initialise_and_read_files(mode, ierr)
    if(ierr.ne.0) call quit_program_error('Could not read ELP-MPP02 files',0)
    
    
    ! Initialization of time powers:
    rjd  = jd - jd2000  ! Reduced JD - JD since 2000
    t(0) = 1.d0
    t(1) = rjd/sc       ! t: time since 2000 in Julian centuries
    t(2) = t(1)**2      ! t^2
    t(3) = t(2)*t(1)    ! t^3
    t(4) = t(2)**2      ! t^4
    
    ! Evaluation of the series: substitution of time in the series
    v = 0.d0
    do iVar=1,3  ! iVar=1,2,3: Longitude, Latitude, Distance
       
       ! Main Problem series:
       do iLine=nmpb(iVar,2),nmpb(iVar,3)
          x = cmpb(iLine)
          y = fmpb(0,iLine)
          yp = 0.d0
          
          do k=1,4
             y  = y  +   fmpb(k,iLine)*t(k)
             yp = yp + k*fmpb(k,iLine)*t(k-1)
          end do  ! k
          
          v(iVar)   = v(iVar)   + x*sin(y)
          v(iVar+3) = v(iVar+3) + x*yp*cos(y)
       end do  ! iLine
       
       ! Perturbations series:
       do it=0,3
          do iLine=nper(iVar,it,2),nper(iVar,it,3)
             x = cper(iLine)
             y = fper(0,iLine)
             xp = 0.d0
             yp = 0.d0
             if(it.ne.0) xp = it * x * t(it-1)
             
             do k=1,4
                y = y   +   fper(k,iLine)*t(k)
                yp = yp + k*fper(k,iLine)*t(k-1)
             end do  ! k
             
             v(iVar)   = v(iVar)   + x * t(it) * sin(y)
             v(iVar+3) = v(iVar+3) + xp*sin(y) + x*t(it)*yp*cos(y)
          end do  ! iLine
       end do  ! it
       
    end do  ! iVar
    
    
    ! Compute the spherical coordinates for the mean inertial ecliptic and equinox of date:
    v(1)   = v(1)/r2as + w(1,0) + w(1,1)*t(1) + w(1,2)*t(2) + w(1,3)*t(3) + w(1,4)*t(4)  ! Longitude + mean longitude (rad)
    !v(1)   = rev(v(1)/r2as) + w(1,0) + rev(w(1,1)*t(1)) + rev(w(1,2)*t(2)) + rev(w(1,3)*t(3)) + rev(w(1,4)*t(4))  ! Longitude + mean longitude (rad)
    v(2)   = v(2)/r2as                                                                   ! Latitude (rad)
    v(3)   = v(3) * a405 / aelp                                                          ! Distance (km)
    
    !lon = v(1)
    !lat = v(2)
    !dist = v(3)
    !lon = lon + (5029.0966d0*t(1) + 1.1120d0*t(2) + 0.000077d0*t(3) - 0.00002353d0*t(4)  -  0.29965d0*t(1)) * as2r  ! Precession from J2000 to EoD(?), but only in longitude!
    !write(*,'(A, F15.9,2F14.7,9F14.5)') 'ELP-MPP02: ', t(1), rev(lon)*r2d, lat*r2d, dist, &
    !     rev(w(1,0) + w(1,1)*t(1) + w(1,2)*t(2) + w(1,3)*t(3) + w(1,4)*t(4))*r2d, &
    !     rev(w(1,0))*r2d, rev(w(1,1)*t(1))*r2d, rev(w(1,2)*t(2))*r2d, rev(w(1,3)*t(3))*r2d, rev(w(1,4)*t(4))*r2d
    
    v(1) = rev(v(1))  ! This adds a bit of CPU time, but also alters the outcome of the cos/sin, similarly to
    !                   taking the 5 rev()s before adding up the terms when computing v(1), which might
    !                   indicate that this gives a better result, especially for dates far from 2000
    
    
    ! Compute the rectangular coordinates (for the EoD?):
    clamb  = cos(v(1))
    slamb  = sin(v(1))
    cbeta  = cos(v(2))
    sbeta  = sin(v(2))
    cw     = v(3)*cbeta
    sw     = v(3)*sbeta
    
    x1     = cw*clamb
    x2     = cw*slamb
    x3     = sw
    
    ! Is this simply precession in rectangular coordinates from EoD to J2000?
    pw     = (p1 + p2*t(1) + p3*t(2) + p4*t(3) + p5*t(4)) * t(1)
    qw     = (q1 + q2*t(1) + q3*t(2) + q4*t(3) + q5*t(4)) * t(1)
    
    ra     = 2*sqrt(1.d0 - pw**2 - qw**2)
    pwqw   = 2*pw*qw
    pw2    = 1.d0 - 2*pw**2
    qw2    = 1.d0 - 2*qw**2
    pwra   = pw*ra
    qwra   = qw*ra
    
    xyz(1) =  pw2*x1  + pwqw*x2 + pwra*x3
    xyz(2) =  pwqw*x1 + qw2*x2  - qwra*x3
    xyz(3) = -pwra*x1 + qwra*x2 + (pw2+qw2-1.d0)*x3
    
    !xyz(1) = x1
    !xyz(2) = x2
    !xyz(3) = x3
    
    
    ! Compute the rectangular velocities for the equinox J2000:
    v(4)   = v(4)/r2as + w(1,1) + 2*w(1,2)*t(1) + 3*w(1,3)*t(2) + 4*w(1,4)*t(3)
    v(5)   = v(5)/r2as
    
    xp1    = (v(6)*cbeta - v(5)*sw)*clamb - v(4)*x2
    xp2    = (v(6)*cbeta - v(5)*sw)*slamb + v(4)*x1
    xp3    = v(6)*sbeta  + v(5)*cw
    
    ppw    = p1 + (2*p2 + 3*p3*t(1) + 4*p4*t(2) + 5*p5*t(3)) * t(1)
    qpw    = q1 + (2*q2 + 3*q3*t(1) + 4*q4*t(2) + 5*q5*t(3)) * t(1)
    ppw2   = -4*pw*ppw
    qpw2   = -4*qw*qpw
    ppwqpw = 2*(ppw*qw + pw*qpw)
    rap    = (ppw2+qpw2)/ra
    ppwra  = ppw*ra + pw*rap
    qpwra  = qpw*ra + qw*rap
    
    vxyz(1) = (pw2*xp1 + pwqw*xp2 + pwra*xp3  +  ppw2*x1 + ppwqpw*x2 + ppwra*x3) / sc
    vxyz(2) = (pwqw*xp1 + qw2*xp2 - qwra*xp3  +  ppwqpw*x1 + qpw2*x2 - qpwra*x3) / sc
    vxyz(3) = (-pwra*xp1 + qwra*xp2 + (pw2+qw2-1.d0)*xp3  -  ppwra*x1 + qpwra*x2 + (ppw2+qpw2)*x3) / sc
    
  end subroutine elp_mpp02_xyz
  !***************************************************************************************************
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Quick, lower-accuracy lunar coordinates; ~600x faster than ELP
  !!
  !! \param jd    Julian day for computation
  !! \param calc  Calculate: 1: l,b,r, 2: & ra,dec, 3: & gmst,agst, 4: & az,alt, nt: number of terms <=60
  !! \param nt    Number of terms to use
  !!
  !! \see
  !! - Meeus, Astronomical Algorithms, 1998, Ch. 22 and 47
  !! - Simon et al, A&A, 282, p.663 (1994)
  
  subroutine moonpos_la(jd, calc,nt)
    use SUFR_kinds, only: double
    use SUFR_constants, only: au, pi,d2r, jd2000, pland,earthr
    use SUFR_angles, only: rev, rev2
    use SUFR_system, only: quit_program_warning
    
    use TheSky_sun, only: sunpos_la
    use TheSky_planetdata, only: planpos, nPlanpos, moonla_arg,moonla_lrb
    use TheSky_coordinates, only: ecl_2_eq, eq2horiz
    use TheSky_datetime, only: calc_deltat, calc_gmst
    use TheSky_local, only: lon0,lat0
    
    implicit none
    real(double), intent(in) :: jd
    integer, intent(in) :: calc,nt
    
    integer :: iLine,iArg, mal(4,60),mab(4,60)
    real(double) :: jde,deltat, tjc,tjc2,tjc3,tjc4, l,b,r,  lm,d,ms,mm,f,e,esl(60),esb(60),a1,a2,a3
    real(double) :: ls,omg,ra,dec,eps,eps0,deps,gmst,agst,lst,dpsi,az,alt,hh, args(4),argl,argb
    real(double) :: moondat(nPlanpos), hcl0,hcb0,hcr0, gcl,gcb,delta, elon,pa,illfr
    
    if(sum(moonla_lrb(1:2,1)) .ne. -14616581) call quit_program_warning('moonpos_la(): moon_la.dat not read', 0)
    
    deltat = calc_deltat(jd)
    jde = jd + deltat/86400.d0
    
    tjc   = (jde-jd2000)/36525.d0  ! Julian Centuries after 2000.0 in dynamical time, the T in Meeus, p.163
    tjc2  = tjc**2
    tjc3  = tjc*tjc2
    tjc4  = tjc2**2
    
    ! Moon's mean longitude, Meeus p.338:
    lm = rev(3.8103408236d0 + 8399.7091116339958d0*tjc - 2.755176757d-5*tjc2 + 3.239043d-8*tjc3 - 2.6771d-10*tjc4)
    
    ! Delauney arguments [d, ms, mm, f] (Meeus p.144) = [D, l', l, F] in Simon et al. 1994, Sect. 3.5:
    !d  = rev(5.19846946025d0 + 7771.37714617d0*tjc - 3.340909d-5*tjc2    + 9.2114446d-8*tjc3)     ! Moon's mean elongation
    !ms = rev(6.24003588115d0 + 628.301956024d0*tjc - 2.79776d-6*tjc2     - 5.8177641733d-8*tjc3)  ! Sun's mean anomaly
    !mm = rev(2.3555483693d0  + 8328.69142288d0*tjc + 1.517947757d-4*tjc2 + 3.102807559d-7*tjc3)   ! Moon's mean anomaly
    !f  = rev(1.62790192912d0 + 8433.46615832d0*tjc - 6.42717497d-5*tjc2  + 5.3329949d-8*tjc3)     ! Moon's argument of latitute
    
    ! Meeus p.338, compares better to ELP-MPP02:
    d  = rev(5.1984665298d0  + 7771.377144834d0*tjc - 3.2845d-5*tjc2  + 3.197347d-8*tjc3    - 1.5436512d-10*tjc4)  ! Moon's mean elongation
    ms = rev(6.240060127d0   + 628.301955167d0*tjc  - 2.681d-6*tjc2   + 7.1267017d-10*tjc3)                        ! Sun's mean anomaly
    mm = rev(2.355555637d0   + 8328.691424759d0*tjc + 1.52566d-4*tjc2 + 2.5041d-7*tjc3      - 1.18633d-9*tjc4)     ! Moon's mean anomaly
    f  = rev(1.627905158d0   + 8433.466158061d0*tjc - 6.3773d-5*tjc2  - 4.94988d-9*tjc3     + 2.02167d-11*tjc4)    ! Moon's argument of latitute
    args = [d,ms,mm,f]  ! Delauney arguments
    
    e  = 1.d0 - 0.002516d0*tjc - 0.0000074d0*tjc2
    
    l  = 0.d0
    r  = 0.d0
    b  = 0.d0
    
    mal = moonla_arg(1:4,:)
    mab = moonla_arg(5:8,:)
    esl = e**abs(moonla_arg(2,:))
    esb = e**abs(moonla_arg(6,:))
    
    do iLine=1,min(nt,60)
       argl = 0.d0
       argb = 0.d0
       do iArg=1,4
          argl = argl + mal(iArg,iLine)*args(iArg)
          argb = argb + mab(iArg,iLine)*args(iArg)
       end do
       l = l + sin(argl) * moonla_lrb(1,iLine) * esl(iLine)
       r = r + cos(argl) * moonla_lrb(2,iLine) * esl(iLine)
       b = b + sin(argb) * moonla_lrb(3,iLine) * esb(iLine)
    end do
    
    
    ! Perturbations by other planets, and flattening of the Earth:
    ! Meeus, p.338:
    a1 = rev(2.090032d0  +     2.301199d0 * tjc)  ! Influence from Venus
    a2 = rev(0.926595d0  + 8364.7398477d0 * tjc)  ! Influence from Jupiter
    a3 = rev(5.4707345d0 + 8399.6847253d0 * tjc)
    
    ! Meeus, p.342:
    l = l + 3958*sin(a1) + 1962*sin(lm-f) + 318*sin(a2)
    b = b - 2235*sin(lm) + 382*sin(a3) + 175*sin(a1-f) + 175*sin(a1+f) + 127*sin(lm-mm) - 115*sin(lm+mm)
    
    
    ! Compute nutation:
    omg  = 2.18243858558d0 - 33.7570459367d0*tjc + 3.6142278d-5*tjc2 + 3.87850944888d-8*tjc3   ! Moon's mean lon. of asc.node, Meeus p.144
    ls   = 4.89506386655d0 + 62.84528862d0*tjc                                                 ! Mean long. Sun, Meeus p.144
    dpsi = -8.338795d-5*sin(omg) - 6.39954d-6*sin(2*ls) - 1.115d-6*sin(2*lm) + 1.018d-6*sin(2*omg)
    
    l = rev( dble(l) * 1.d-6*d2r + lm + dpsi)
    b = rev2(dble(b) * 1.d-6*d2r)
    r = (dble(r*100) + 3.8500056d10)/au
    
    
    ! Store results:
    planpos(1)   = l  ! Geocentric ecliptic longitude
    planpos(2)   = b  ! Geocentric ecliptic latitude
    planpos(4)   = r  ! Geocentric distance
    
    planpos(7)   = 5.77551830441d-3 * r  ! Light time in days
    planpos(12)  = pland(0)/(r*au)       ! Apparent diameter of the Moon
       
    planpos(40)  = jde   ! JDE
    planpos(46)  = tjc   ! App. dyn. time in Julian Centuries since 2000.0
    planpos(47)  = dpsi  ! Nutation in longitude
    
    if(calc.eq.1) return
    
    
    ! Obliquity of the ecliptic:
    eps0 = 0.409092804222d0 - 2.26965525d-4*tjc - 2.86d-9*tjc2 + 8.78967d-9*tjc3            ! Mean obliquity of the ecliptic (Lieske et al, 1977)
    deps = 4.468d-5*cos(omg) + 2.76d-6*cos(2*ls) + 4.848d-7*cos(2*lm) - 4.36d-7*cos(2*omg)  ! Nutation in obliquity
    eps  = eps0 + deps                                                                      ! True obliquity of the ecliptic
    
    call ecl_2_eq(l,b,eps, ra,dec)  ! Ecliptical -> equatorial coordinates
    
    ! Store results:
    planpos(5)  = ra    ! Geocentric right ascension
    planpos(6)  = dec   ! Geocentric declination
    
    planpos(48) = eps   ! True obliquity of the ecliptic; corrected for nutation
    planpos(50) = eps0  ! Mean obliquity of the ecliptic; without nutation
    
    if(calc.eq.2) return
    
    
    ! Sidereal time:
    gmst = calc_gmst(jd)                 ! Greenwich mean sidereal time
    agst = rev(gmst + dpsi*cos(eps))     ! Correction for equation of the equinoxes -> Gr. apparent sid. time
    lst  = rev(agst + lon0)              ! Local apparent sidereal time, lon0 > 0 for E
    
    planpos(44) = lst                    ! Local APPARENT sidereal time
    planpos(45) = agst                   ! Greenwich APPARENT sidereal time (in radians)
    planpos(49) = gmst                   ! Greenwich MEAN sidereal time (in radians)
    
    if(calc.eq.3) return
    
    
    call eq2horiz(ra,dec,agst, hh,az,alt)  ! Equatorial -> horizontal coordinates
    
    planpos(8)  = hh   ! Geocentric hour angle
    planpos(9)  = az   ! Geocentric azimuth
    planpos(10) = alt  ! Geocentric altitude
    
    if(calc.eq.4) return
    
    
    moondat = planpos  ! Store Moon data
    
    ! Get some Sun data:
    call sunpos_la(jd,1)
    hcl0 =  rev(planpos(1)+pi) ! Heliocentric longitude of the Earth - want geocentric lon of Sun?
    hcb0 = -planpos(2)         ! Heliocentric latitude of the Earth  - want geocentric lat of Sun?
    hcr0 =  planpos(3)
    
    planpos = moondat   ! Restore Moon data:
    gcl   = planpos(1)  ! l
    gcb   = planpos(2)  ! b
    delta = planpos(4)  ! r
    
    
    elon = acos(cos(gcb)*cos(hcb0)*cos(gcl-rev(hcl0+pi)))
    pa = atan2( hcr0*sin(elon) , delta - hcr0*cos(elon) )            ! Phase angle
    illfr = 0.5d0*(1.d0 + cos(pa))                                   ! Illuminated fraction of the disc
    
    planpos(11) = rev(elon)           ! Elongation
    planpos(13) = moonmagn(pa,delta)  ! Magnitude
    planpos(14) = illfr               ! Illuminated fraction
    planpos(15) = rev(pa)             ! Phase angle
    planpos(16) = atan2(sin(planpos(8)),tan(lat0)*cos(planpos(6)) - sin(planpos(6))*cos(planpos(8)))  ! Parallactic angle
    planpos(17) = asin(earthr/(planpos(4)*au))                                                        ! Horizontal parallax
    
    planpos(41) = rev(hcl0+pi)        ! Geocentric, true L,B,R for the Sun
    planpos(42) = rev2(-hcb0)
    planpos(43) = hcr0
    
  end subroutine moonpos_la
  !*********************************************************************************************************************************
    
  
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate the magnitude of the Moon
  !!
  !! \param pa     Phase angle (rad)
  !! \param delta  Geocentric(?) distance (AU)
  !!
  !! \see
  !! - Allen, 1976, par.66
  !! - http://cleardarksky.com/others/BenSugerman/star.htm
  
  function moonmagn(pa, delta)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: pa,delta
    real(double) :: moonmagn
    
    moonmagn = -12.73d0 + 1.489690267d0*abs(pa) + 0.04310727d0*pa**4  ! Allen, 1976, par.66 (d2r)
    
    ! Correct for variable distance and the 'opposition effect' or 'opposition surge' which occurs when pa < 7d:
    moonmagn = moonmagn - 2.5*log10( (2.5696d-3/delta)**2  *  max(1.d0, 1.35d0 - 2.864789d0*abs(pa) ) )  ! (1-1.35)/2.865 ~ 7d
    
  end function moonmagn
  !*********************************************************************************************************************************
  
  
  
end module TheSky_moon
!***********************************************************************************************************************************
