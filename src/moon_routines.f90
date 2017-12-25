!> \file moon_routines.f90  Procedures that calculate the physical data, phases and age of the Moon for libTheSky


!  Copyright (c) 2002-2017  Marc van der Sluys - marc.vandersluys.nl
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

module TheSky_moonroutines
  implicit none
  save
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Get physical data for the Moon: librations, position angles, selenographic position of the Sun
  !!
  !! \param  jd    Julian day for computation
  !!
  !! \retval libl  Libration (physical+optical) in longitude
  !! \retval libb  Libration (physical+optical) in latitude
  !!
  !! \retval pa    Position angle of the Moon's axis/north pole
  !! \retval blpa  Position angle of the Moon's bright limb
  !!
  !! \retval sunl  Selenographic longitude of the Sun
  !! \retval sunb  Selenographic latitude of the Sun
  !!
  !!
  !! \see  Meeus, Astronomical Algorithms, 1998, Ch. 53
  !!
  !! \note  This routine does NOT save Moon data in planpos !!!
  
  subroutine moonphys(jd, libl,libb, pa,blpa, sunl,sunb)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi, d2r
    use SUFR_angles, only: rev, rev2
    
    use TheSky_planets, only: planet_position
    use TheSky_planetdata, only: planpos, nplanpos
    
    implicit none
    real(double), intent(in) :: jd
    real(double), intent(out) :: libl,libb, pa,blpa, sunl,sunb
    real(double) :: moonpos(nplanpos),storepos(nplanpos), tjc,tjc2,tjc3,tjc4, dpsi,eps,lm,bm,  omg,mal,mem,mas,mam,ee,k1,k2,in,ww,aa
    real(double) :: rho,sig,tau,libl2,libb2,vv,xx,yy,om,  a0,d0,l0,r0,lh,bh,sunl2,sunb2
    
    storepos = planpos   ! Store current data
    
    call planet_position(jd,0)  ! Moon
    moonpos = planpos
    
    lm = moonpos(1)      ! Geocentric longitude
    bm = moonpos(2)      ! Geocentric latitude
    
    tjc    = moonpos(46)   ! Apparent dynamical time in Julian Centuries since 2000.0
    tjc2   = tjc**2        ! t^2
    tjc3   = tjc2*tjc      ! t^3
    tjc4   = tjc2**2       ! t^4
    
    dpsi = moonpos(47)   ! Nutation in longitude
    eps  = moonpos(48)   ! True obliquity of the ecliptic, corrected for nutation
    
    
    omg = 2.1824390725d0  - 33.7570464271d0    *tjc  + 3.622256d-5   *tjc2   + 3.7337958d-8 *tjc3  - 2.879321d-10  *tjc4  ! Moon's longitude of mean ascending node
    mal = 1.62790515798d0 + 8433.46615806092d0 *tjc  - 6.37725855d-5 *tjc2   - 4.9498844d-9 *tjc3  + 2.0216715d-11 *tjc4  ! Moon's argument of latitude (F)
    
    ! Lower accuracy(?) from Meeus Ch. 22:
    mem = 5.19846946025d0  + 7771.37714617d0 *tjc    - 3.340909d-5    *tjc2  + 9.2114446d-8    *tjc3  ! Mean elongation of the Moon
    mas = 6.24003588115d0  + 628.301956024d0 *tjc    - 2.79776d-6     *tjc2  - 5.8177641733d-8 *tjc3  ! Mean anomaly of the Sun (or Earth)
    mam = 2.3555483693d0   + 8328.69142288d0 *tjc    + 1.517947757d-4 *tjc2  + 3.102807559d-7  *tjc3  ! Mean anomaly of the Moon
    ee  = 1.d0  - 0.002516d0 *tjc  - 0.0000074 *tjc2
    k1  = 2.09003d0  + 2.301199d0 *tjc
    k2  = 1.2664d0   + 0.352312d0 *tjc
    
    
    ! Optical librations:
    in = 2.69203d-2         ! Moon equator inclination to ecliptic, Meeus p.372
    
    ! Meeus, Eq. 53.1:
    ww  = lm - dpsi - omg
    aa  = atan2( sin(ww)*cos(bm)*cos(in) - sin(bm)*sin(in),  cos(ww)*cos(bm) )
    
    libl = aa - mal   ! Along with AA: - = W on celestial sphere, E in selenographic coordinates - opposed to Sterrengids
    libb = asin( -sin(ww)*cos(bm)*sin(in) - sin(bm)*cos(in) )
    
    
    ! Physical librations:
    ! Meeus, p.373, in degrees:
    rho = &
         - 2.752d-2*cos(mam)    - 2.245d-2*sin(mal)        + 6.84d-3*cos(mam-2*mal) &
         - 2.93d-3*cos(2*mal)   - 8.5d-4*cos(2*(mal-mem))  - 5.4d-4*cos(mam-2*mem)  &
         - 2.d-4*sin(mam+mal)   - 2.d-4*cos(mam+2*mal)     - 2.d-4*cos(mam-mal)     &
         + 1.4d-4*cos(mam+2*(mal-mem))
    
    sig = &
         - 2.816d-2*sin(mam)    + 2.244d-2*cos(mal)        - 6.82d-3*sin(mam-2*mal) &
         - 2.79d-3*sin(2*mal)   - 8.3d-4*sin(2*(mal-mem))  + 6.9d-4*sin(mam-2*mem)  &
         + 4.d-4*cos(mam+mal)   - 2.5d-4*sin(2*mam)        - 2.3d-4*sin(mam+2*mal)  &
         + 2.d-4*cos(mam-mal)   + 1.9d-4*sin(mam-mal)                               &
         + 1.3d-4*sin(mam+2*(mal-mem))                     - 1.d-4*cos(mam-3*mal)
    
    tau =  2.520d-2*ee*sin(mas)         + 4.74d-3*sin(2*(mam-mal))  - 4.67d-3*sin(mam)        &
         + 3.96d-3*sin(k1)              + 2.76d-3*sin(2*(mam-mem))  + 1.96d-3*sin(omg)        &
         - 1.83d-3*cos(mam-mal)         + 1.15d-3*sin(mam-2*mem)    - 9.6d-4*sin(mam-mem)     &
         + 4.6d-4*sin(2*(mal-mem))      - 3.9d-4*sin(mam-mal)       - 3.2d-4*sin(mam-mas-mem) &
         + 2.7d-4*sin(2*(mam-mem)-mas)  + 2.3d-4*sin(k2)            - 1.4d-4*sin(2*mem)       &
         + 1.4d-4*cos(2*(mam-mal))      - 1.2d-4*sin(mam-2*mal)     - 1.2d-4*sin(2*mam)       &
         + 1.1d-4*sin(2*(mam-mas-mem))
    
    rho = rho*d2r
    sig = sig*d2r
    tau = tau*d2r
    
    ! Meeus, Eq. 53.2:
    libl2 = -tau + (rho*cos(aa) + sig*sin(aa)) * tan(libb)
    libb2 = sig*cos(aa) - rho*sin(aa)
    
    ! Total librations:
    libl = libl + libl2  ! Along with AA ( - = W on celestial sphere, E in selenographic coordinates)
    libb = libb + libb2
    
    
    ! Position Angle of the axis:
    ! Meeus, p.374:
    vv = omg + dpsi + sig/sin(in)
    xx = sin(in+rho) * sin(vv)
    yy = sin(in+rho) * cos(vv) * cos(eps)  -  cos(in+rho) * sin(eps)
    om = atan2(xx,yy)
    pa = asin( sqrt(xx*xx+yy*yy) * cos(moonpos(5)-om) / cos(libb) )
    
    
    call planet_position(jd,3)  ! Sun - CHECK - need full accuracy?
    l0 = planpos(1)  ! Geocentric ecliptic longitude of the Sun
    r0 = planpos(3)  ! Geocentric distance of the Sun
    a0 = planpos(5)  ! Geocentric right ascension of the Sun
    d0 = planpos(6)  ! Geocentric declination of the Sun
    
    ! Position angle of the bright limb:
    blpa = atan2( cos(d0)*sin(a0-moonpos(5)),  sin(d0)*cos(moonpos(6)) - cos(d0)*sin(moonpos(6))*cos(a0-moonpos(5)) )
    
    
    ! Selenographic position of the Sun:
    ! Meeus, p.376:
    
    ! Heliocentric l,b:
    lh = l0 + pi + moonpos(4)/r0*cos(bm)*sin(l0-lm)
    bh = moonpos(4)/r0*bm
    ww  = lh - dpsi - omg
    aa  = atan2( sin(ww)*cos(bh)*cos(in) - sin(bh)*sin(in),  cos(ww)*cos(bh) )
    
    ! Selenographic coordinates of the Sun:
    sunl  = aa - mal
    sunb  = asin( -sin(ww) * cos(bh) * sin(in)  -  sin(bh) * cos(in))
    sunl2 = -tau + (rho*cos(aa) + sig*sin(aa)) * tan(sunb)
    sunb2 = sig*cos(aa) - rho*sin(aa)
    sunl  = rev(sunl + sunl2)
    sunb  = rev2(sunb + sunb2)
    
    
    ! Restore current data:
    planpos = storepos
    
  end subroutine moonphys
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculates Julian Day of phase of the Moon for the desired phase k0
  !!
  !! \param k0   Desired phase:  x.00: New Moon - k = x.75: Last Quarter.  k=0 ~ 2000.0
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch. 49
  !!
  !! \note
  !! - First version: June 2003
  
  function moonphase(k0)
    use SUFR_kinds, only: double, dbl
    use SUFR_constants, only: d2r
    use SUFR_angles, only: rev
    
    use TheSky_local, only: deltat
    
    implicit none
    real(double), intent(in) :: k0
    real(double) :: moonphase,k
    real(double) :: t,t2,jde, ee,mas,mam,ff,omg, aa(14),ww
    integer :: phase,i
    
    k = k0
    phase = nint((k - int(k))*4.0_dbl)
    if(phase.lt.0) phase = phase + 4
    if(phase.gt.3) phase = phase - 4
    
    t = k/1236.85d0  ! Time in Julian Centuries
    t2 = t*t
    
    ! First approximation:
    jde = 2451550.09766d0 + 29.530588861d0*k + t2*(0.00015437 + t*(-0.000000150 + 0.00000000073*t))  ! Meeus, Eq. 49.1
    
    ! Return very rough estimate (~day):
    !moonphase = jde - deltat/86400.d0
    !return
    
    ! Correction parameters:
    ee  = 1.0d0      + t * (-0.002516d0 - 0.0000074d0 * t)  ! Meeus, Eq. 47.6
    
    ! Mean anomaly of the Sun, Meeus, Eq. 49.4:
    mas = 2.5534d0   + 29.10535669d0  * k  +  t2 * (-0.0000218d0 - 0.00000011d0 * t)
    
    ! Mean anomaly of the Moon, Meeus, Eq. 49.5:
    mam = 201.5643d0 + 385.81693528d0 * k  +  t2 * (0.0107438d0  + t*(0.00001239d0 - 0.000000058d0 * t))
    
    ! Moon's argument of latitude, Meeus, Eq. 49.6:
    ff  = 160.7108d0 + 390.67050274d0 * k  +  t2 * (-0.0016341d0*t * (-0.00000227d0 + 0.000000011d0 * t))
    
    ! Longitude of ascending node of lunar orbit, Meeus, Eq. 49.7:
    omg = 124.7746d0 - 1.56375580d0   * k  +  t2 * (0.0020691d0 + 0.00000215d0 * t)
    
    
    ! Planetary arguments:
    ! Meeus, p.351:
    aa(1)  = 299.77d0  +   0.107408d0 * k  -  0.009173d0 * t2
    aa(2)  = 251.88d0  +   0.016321d0 * k
    aa(3)  = 251.83d0  +  26.651886d0 * k
    aa(4)  = 349.42d0  +  36.412478d0 * k
    aa(5)  =  84.66d0  +  18.206239d0 * k
    aa(6)  = 141.74d0  +  53.303771d0 * k
    aa(7)  = 207.14d0  +   2.453732d0 * k
    aa(8)  = 154.84d0  +   7.306860d0 * k
    aa(9)  =  34.52d0  +  27.261239d0 * k
    aa(10) = 207.19d0  +   0.121824d0 * k
    aa(11) = 291.34d0  +   1.844379d0 * k
    aa(12) = 161.72d0  +  24.198154d0 * k
    aa(13) = 239.56d0  +  25.513099d0 * k
    aa(14) = 331.55d0  +   3.592518d0 * k
    
    ! Convert degrees -> radians:
    mas = rev(mas*d2r)
    mam = rev(mam*d2r)
    ff  = rev(ff*d2r)
    omg = rev(omg*d2r)
    
    do i=1,14
       aa(i)  = rev(aa(i)*d2r)
    end do
    
    
    
    !*******************************************************************************************************************************
    ! New Moon (phase 0), Meeus p.351:
    if(phase.eq.0) then
       jde = jde &
            - 0.40720d0 * sin(mam)            +  0.17241d0 * ee * sin(mas)            &
            + 0.01608d0 * sin(2*mam)          +  0.01039d0 * sin(2*ff)                &
            + 0.00739d0 * ee * sin(mam-mas)   -  0.00514d0 * ee * sin(mam+mas)        &
            + 0.00208d0 * ee**2 * sin(2*mas)  -  0.00111d0 * sin(mam-2*ff)            &
            - 0.00057d0 * sin(mam+2*ff)       +  0.00056d0 * ee * sin(2*mam+mas)      &
            - 0.00042d0 * sin(3*mam)          +  0.00042d0 * ee * sin(mas+2*ff)       &
            + 0.00038d0 * ee * sin(mas-2*ff)  -  0.00024d0 * ee * sin(2*mam-mas)      &
            - 0.00017d0 * sin(omg)            -  0.00007d0 * sin(mam+2*mas)           &
            + 0.00004d0 * sin(2*mam-2*ff)     +  0.00004d0 * sin(3*mas)               &
            + 0.00003d0 * sin(mam+mas-2*ff)   +  0.00003d0 * sin(2*mam+2*ff)          &
            - 0.00003d0 * sin(mam+mas+2*ff)   +  0.00003d0 * sin(mam-mas+2*ff)        &
            - 0.00002d0 * sin(mam-mas-2*ff)   -  0.00002d0 * sin(3*mam+mas)           &
            + 0.00002d0 * sin(4*mam)
    end if
    
    
    !*******************************************************************************************************************************
    ! Full Moon  (phase 2), Meeus p.351:
    if(phase.eq.2) then      
       jde = jde &
            - 0.40614d0 * sin(mam)            +  0.17302d0 * ee * sin(mas)            &
            + 0.01614d0 * sin(2*mam)          +  0.01043d0 * sin(2*ff)                &
            + 0.00734d0 * ee * sin(mam-mas)   -  0.00515d0 * ee * sin(mam+mas)        &
            + 0.00209d0 * ee**2 * sin(2*mas)  -  0.00111d0 * sin(mam-2*ff)            &
            - 0.00057d0 * sin(mam+2*ff)       +  0.00056d0 * ee * sin(2*mam+mas)      &
            - 0.00042d0 * sin(3*mam)          +  0.00042d0 * ee * sin(mas+2*ff)       &
            + 0.00038d0 * ee * sin(mas-2*ff)  -  0.00024d0 * ee * sin(2*mam-mas)      &
            - 0.00017d0 * sin(omg)            -  0.00007d0 * sin(mam+2*mas)           &
            + 0.00004d0 * sin(2*mam-2*ff)     +  0.00004d0 * sin(3*mas)               &
            + 0.00003d0 * sin(mam+mas-2*ff)   +  0.00003d0 * sin(2*mam+2*ff)          &
            - 0.00003d0 * sin(mam+mas+2*ff)   +  0.00003d0 * sin(mam-mas+2*ff)        &
            - 0.00002d0 * sin(mam-mas-2*ff)   -  0.00002d0 * sin(3*mam+mas)           &
            + 0.00002d0 * sin(4*mam)
    end if
    
    
    !*******************************************************************************************************************************
    ! First, last quarter  (phase 1,3), Meeus p.352:
    if(phase.eq.1.or.phase.eq.3) then      
       jde = jde &
            - 0.62801d0 * sin(mam)             +  0.17172d0 * ee * sin(mas)           &
            - 0.01183d0 * ee * sin(mam+mas)    +  0.00862d0 * sin(2*mam)              &
            + 0.00804d0 * sin(2*ff)            +  0.00454d0 * ee * sin(mam-mas)       &
            + 0.00204d0 * ee**2 * sin(2*mas)   -  0.00180d0 * sin(mam-2*ff)           &
            - 0.00070d0 * sin(mam+2*ff)        -  0.00040d0 * sin(3*mam)              &
            - 0.00034d0 * ee * sin(2*mam-mas)  +  0.00032d0 * ee * sin(mas+2*ff)      &
            + 0.00032d0 * ee * sin(mas-2*ff)   -  0.00028d0 * ee**2 * sin(mam+2*mas)  &
            + 0.00027d0 * ee * sin(2*mam+mas)  -  0.00017d0 * sin(omg)                &
            - 0.00005d0 * sin(mam-mas-2*ff)                                           &
            + 0.00004d0 * sin(2*mam+2*ff)      -  0.00004d0 * sin(mam+mas+2*ff)       &
            + 0.00004d0 * sin(mam-2*mas)       +  0.00003d0 * sin(mam+mas-2*ff)       &
            + 0.00003d0 * sin(3*mas)           +  0.00002d0 * sin(2*mam-2*ff)         &
            + 0.00002d0 * sin(mam-mas+2*ff)    -  0.00002d0 * sin(3*mam+mas)
       
       ww = 0.00306 &
            - 0.00038d0 * ee * cos(mas)        +  0.00026d0 * cos(mam)                &
            - 0.00002d0 * cos(mam-mas)         +  0.00002d0 * cos(mam+mas)            &
            + 0.00002d0 * cos(2*ff)
       
       if(phase.eq.3) ww = -ww
       jde = jde + ww
    end if
    
    
    ! Final corrections, Meeus p.352:
    jde = jde &
         + 0.000325d0 * sin(aa(1))   +  0.000165d0 * sin(aa(2))   &
         + 0.000164d0 * sin(aa(3))   +  0.000126d0 * sin(aa(4))   &
         + 0.000110d0 * sin(aa(5))   +  0.000062d0 * sin(aa(6))   &
         + 0.000060d0 * sin(aa(7))   +  0.000056d0 * sin(aa(8))   &
         + 0.000047d0 * sin(aa(9))   +  0.000042d0 * sin(aa(10))  &
         + 0.000040d0 * sin(aa(11))  +  0.000037d0 * sin(aa(12))  &
         + 0.000035d0 * sin(aa(13))  +  0.000023d0 * sin(aa(14))
    
    
    moonphase = jde - deltat/86400.d0
    
  end function moonphase
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the age of the Moon for a given JD
  !!
  !! \param  jd   Julian day for computation
  !! \retval age  Age of the Moon (since last New Moon) in days
  
  subroutine moon_age(jd, age)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000
    
    implicit none
    real(double), intent(in) :: jd
    real(double), intent(out) :: age
    real(double) :: djd, jd0, k0
    
    
    ! Moon age:
    djd = jd - jd2000 - 0.5d0
    
    ! Overshoot by a month, to make sure you don't miss one.  Age.gt.29.530589d0 doesn't work, since this is the MEAN month:
    k0 = floor(djd/365.25*12.3685d0) + 1.d0
    age = -1.d0
    
    do while(age.lt.0.d0)
       jd0 = moonphase(k0)
       age = jd-jd0
       k0 = k0 - 1.d0  ! Previous New Moon
    end do
    
  end subroutine moon_age
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the next lunar phase after a given JD
  !!
  !! \param  JD      Julian day for computation
  !! \retval phase   Next lunar phase: 0-NM ... 3-LQ
  !! \retval JDnext  Julian day of next lunar phase
  
  subroutine moon_phase_next(JD, phase, JDnext)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000
  use TheSky_datetime, only: printdate1
    
    implicit none
    real(double), intent(in) :: JD
    integer, intent(out) :: phase
    real(double), intent(out) :: JDnext
    real(double) :: dJD, k0
    
    
    ! Undershoot by half a month, to make sure you don't miss anything:
    dJD = JD - jd2000 - 0.5d0
    k0 = floor(dJD/365.25*12.3685d0) - 0.5d0
    JDnext = JD - 30.d0
    
    do while(JDnext.lt.JD)
       k0 = k0 + 0.25d0  ! Next phase
       JDnext = moonphase(k0)
    end do
    
    phase = nint( (k0 - floor(k0))*4.d0 )  ! 0-3 for NM - LQ
    
  end subroutine moon_phase_next
  !*********************************************************************************************************************************
  
  
end module TheSky_moonroutines
!***********************************************************************************************************************************
