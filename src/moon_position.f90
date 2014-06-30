!> \file moon_position.f90  Core procedures that calculate the position and magnitude of the Moon for libTheSky


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
  !! \retval ll   Apparent geocentric ecliptical longitude
  !! \retval bb   Apparent geocentric ecliptical latitude
  !! \retval rr   Apparent geocentric distance
  !!
  !! \see ftp://cdsarc.u-strasbg.fr/pub/cats/VI/79/
  
  subroutine moon_lbr(tjj, ll,bb,rr)
    use SUFR_kinds, only: double, dbl
    use SUFR_constants, only: as2r, au, km
    use SUFR_angles, only: rev
    use TheSky_moondata, only: t, nterm,nrang, pc1,pc2,pc3, per1,per2,per3, w, ath,a0
    
    implicit none
    real(double), intent(in) :: tjj
    real(double), intent(out) :: ll,bb,rr
    
    integer :: iv,itab,nt,k,j
    real(double) :: r(3),x,y, pa,t4,t8
    
    
    r(1) = 0.0_dbl
    r(2) = 0.0_dbl
    r(3) = 0.0_dbl
    
    x = 0.0_dbl
    y = 0.0_dbl
    
    t(1) = tjj*10.0_dbl  ! In centuries
    t(2) = t(1)*t(1)     ! t^2
    t(3) = t(2)*t(1)     ! t^3
    t(4) = t(3)*t(1)     ! t^4
    
    t4 = t(4)   ! t^4
    t8 = t4*t4  ! t^8
    
    
    do iv=1,3
       r(iv) = 0.0_dbl
       
       do itab = 1,12
          do nt = 1,nterm(iv,itab)
             
             select case(iv)
                
             case(1)
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
                
             case(2)
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
                
             case(3)
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
    
    pa =   2.438174835301452d-2  * t(1)     + 5.391128133938040d-6  * t(2)     + 3.733065344543427d-10 * t(3)    &
         - 1.140766591650738d-10 * t4       - 8.753311012432670d-14 * t4*t(1)  + 8.460483549042511d-16 * t4*t(2) &
         + 6.348635154129373d-18 * t4*t(3)  + 1.175188363009515E-20 * t8       - 2.307228308400282d-22 * t8*t(1) &
         - 4.198486478408581d-25 * t8*t(2)
    
    r(1) = r(1) + pa   ! Precess to equinox of date, see ELP PS-file, p12
    
    !!lbr2xyz: spherical to rectangular coordinates
    !x1 = r(3)*cos(r(2))
    !x2 = x1*sin(r(1))
    !x1 = x1*cos(r(1))
    !x3 = r(3)*sin(r(2))
    !xx = rev(x1)
    !yy = x2
    !zz = x3
    
    
    ll = rev(r(1))
    bb = r(2)
    rr = r(3)/au*km
    
    !write(98,*)ll,bb,rr
    
  end subroutine moon_lbr
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Quick, lower-accuracy lunar coordinates; ~600x faster than ELP
  !!
  !! \param jd    Julian day for computation
  !! \param calc  Calculate: 1: l,b,r, 2: & ra,dec, 3: & gmst,agst, 4: & az,alt, nt: number of terms <=60
  !! \param nt    Number of terms to use
  !!
  !! \see  Meeus, Astronomical Algorithms, 1998, Ch. 22 and 47
  !!
  !! \todo
  !! - Improve some parts (search for 'improve' below): use Meeus Ch.47/p.338 rather than Ch.22/p.144
  
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
    
    integer :: i,j,mal(4,60),mab(4,60)
    real(double) :: jde,deltat,t,t2,t3,l,b,r,  lm,d,ms,mm,f,e,esl(60),esb(60),a1,a2,a3
    real(double) :: ls,omg,ra,dec,eps,eps0,deps,gmst,agst,lst,dpsi,az,alt,hh, args(4),argl,argb
    real(double) :: moondat(nPlanpos), hcl0,hcb0,hcr0, gcl,gcb,delta, elon,pa,illfr
    
    if(sum(moonla_lrb(1:2,1)) .ne. -14616581) call quit_program_warning('moonpos_la(): moon_la.dat not read', 0)
    
    deltat = calc_deltat(jd)
    jde = jd + deltat/86400.d0
    t   = (jde-jd2000)/36525.d0  ! Julian Centuries after 2000.0 in dynamical time, the T in Meeus, p.163
    t2  = t*t
    t3  = t*t2
    
    ! These can be improved somewhat (by using Meeus Ch.47/p.338 rather than Ch.22/p.144):
    ! In fact, these values may be 'incompatible' with the periodic terms used
    lm = rev(3.8103417d0 + 8399.709113d0*t)                                                 ! Moon's mean longitude, Meeus p.144
    d  = rev(5.19846946025d0 + 7771.37714617d0*t - 3.340909d-5*t2    + 9.2114446d-8*t3)     ! Moon's mean elongation, Meeus p.144
    ms = rev(6.24003588115d0 + 628.301956024d0*t - 2.79776d-6*t2     - 5.8177641733d-8*t3)  ! Sun's mean anomaly, Meeus p.144
    mm = rev(2.3555483693d0  + 8328.69142288d0*t + 1.517947757d-4*t2 + 3.102807559d-7*t3)   ! Moon's mean anomaly, Meeus p.144
    f  = rev(1.62790192912d0 + 8433.46615832d0*t - 6.42717497d-5*t2  + 5.3329949d-8*t3)     ! Moon's argument of latit., Meeus p.144
    e  = 1.d0 - 0.002516d0*t - 0.0000074*t2
    args = (/d,ms,mm,f/)
    
    ! Meeus, p.338:
    a1 = rev(2.090032d0  + 2.301199d0 * t)
    a2 = rev(0.926595d0  + 8364.7398477d0 * t)
    a3 = rev(5.4707345d0 + 8399.6847253d0 * t)
    
    l  = 0.d0
    r  = 0.d0
    b  = 0.d0
    
    mal = moonla_arg(1:4,:)
    mab = moonla_arg(5:8,:)
    esl = e**abs(moonla_arg(2,:))
    esb = e**abs(moonla_arg(6,:))
    
    do i=1,min(nt,60)
       argl = 0.d0
       argb = 0.d0
       do j=1,4
          argl = argl + mal(j,i)*args(j)
          argb = argb + mab(j,i)*args(j)
       end do
       l = l + sin(argl) * moonla_lrb(1,i) * esl(i)
       r = r + cos(argl) * moonla_lrb(2,i) * esl(i)
       b = b + sin(argb) * moonla_lrb(3,i) * esb(i)
    end do
    
    ! Meeus, p.342:
    l = l + 3958*sin(a1) + 1962*sin(lm-f) + 318*sin(a2)
    b = b - 2235*sin(lm) + 382*sin(a3) + 175*sin(a1-f) + 175*sin(a1+f) + 127*sin(lm-mm) - 115*sin(lm+mm)
    
    ! Compute nutation:
    omg  = 2.18243858558d0 - 33.7570459367d0*t + 3.6142278d-5*t2 + 3.87850944888d-8*t3   ! Moon's mean lon. of asc.node, Meeus p.144
    ls   = 4.89506386655d0 + 62.84528862d0*t                                             ! Mean long. Sun, Meeus p.144
    dpsi = -8.338795d-5*sin(omg) - 6.39954d-6*sin(2*ls) - 1.115d-6*sin(2*lm) + 1.018d-6*sin(2*omg)
    
    l = rev(dble(l)*1.d-6*d2r + lm + dpsi)
    b = rev2(dble(b)*1.d-6*d2r)
    r = (dble(r*100) + 3.8500056d10)/au
    
    ! Store results:
    planpos(1)   = l  ! Geocentric ecliptic longitude
    planpos(2)   = b  ! Geocentric ecliptic latitude
    planpos(4)   = r  ! Geocentric distance
    
    planpos(7)   = 5.77551830441d-3 * r  ! Light time in days
    planpos(12)  = pland(0)/(r*au)       ! Apparent diameter of the Moon
       
    planpos(40)  = jde   ! JDE
    planpos(46)  = t     ! App. dyn. time in Julian Centuries since 2000.0
    planpos(47)  = dpsi  ! Nutation in longitude
    
    if(calc.eq.1) return
    
    
    ! Obliquity of the ecliptic:
    eps0 = 0.409092804222d0 - 2.26965525d-4*t - 2.86d-9*t2 + 8.78967d-9*t2*t                ! Mean obliquity of the ecliptic
    deps = 4.468d-5*cos(omg) + 2.76d-6*cos(2*ls) + 4.848d-7*cos(2*lm) - 4.36d-7*cos(2*omg)  ! Nutation in obliquity
    eps  = eps0 + deps                                                                      ! True obliquity of the ecliptic
    
    call ecl_2_eq(l,b,eps, ra,dec)  ! Ecliptical -> equatorial coordinates
    
    ! Store results:
    planpos(5)  = ra    ! Geocentric right ascension
    planpos(6)  = dec   ! Geocentric declination
    
    planpos(48) = eps   ! True obliquity of the ecliptic; corrected for nutation
    planpos(50) = eps0  ! Mean obliquity of the ecliptic; without nutation
    
    if(calc.eq.2) return
    
    
    ! Siderial time:
    gmst = calc_gmst(jd)                 ! Greenwich mean siderial time
    agst = rev(gmst + dpsi*cos(eps))     ! Correction for equation of the equinoxes -> Gr. apparent sid. time
    lst  = rev(agst + lon0)              ! Local apparent siderial time, lon0 > 0 for E
    
    planpos(44) = lst                    ! Local APPARENT siderial time
    planpos(45) = agst                   ! Greenwich APPARENT siderial time (in radians)
    planpos(49) = gmst                   ! Greenwich MEAN siderial time (in radians)
    
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
  !! \param pa     Phase angle
  !! \param delta  Geocentric(?) distance
  !!
  !! \see  http://cleardarksky.com/others/BenSugerman/star.htm, but converted to radians
  
  function moonmagn(pa, delta)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: pa,delta
    real(double) :: moonmagn,magn,dist,ill
    
    magn = -12.73d0 + 1.48969d0*abs(pa) + 0.043107d0*pa**4  ! Allen, 1976
    dist = 2.5696d-3/delta                                  ! Average distance over current distance, in AU
    
    ! Correct for variable distance and 'opposition effect' if pa < 7d:
    ill = 10.d0**(-0.4d0*(magn+16.57d0))  *  dist**2  *  max(1.d0,(1.35d0 - 2.864789d0*abs(pa)))
    moonmagn = -2.5d0*log10(ill) - 16.57d0
    
  end function moonmagn
  !*********************************************************************************************************************************
  
  
  
end module TheSky_moon
!***********************************************************************************************************************************
