!> \file riset.f90  Compute rise, transit and set times, or beginning and end of twilight for libTheSky


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
!> \brief Procedures to determine rise, transit and set times

module TheSky_riset
  implicit none
  save
  
contains
  
  !*********************************************************************************************************************************
  !> \brief  Rise, transit and set times routine for Sun, Moon, planets and asteroids.
  !!         Based on riset_ipol(), but recomputes positions (low accuracy first, then full accuracy) rather than interpolating.
  !!
  !! \param jd     Julian day number
  !! \param pl     Planet/object number  (-1 - 9: ES, Moon, Mer-Plu; >10000 asteroids)
  !! \param sa0    Altitude to return rise/set data for (degrees; 0. is actual rise/set).  sa0>90: compute transit only
  !!
  !! \retval rt    Rise time (hours)
  !! \retval tt    Transit time (hours)
  !! \retval st    Set time (hours)
  !! \retval rh    Rising wind direction (rad)
  !! \retval ta    Transit altitude (rad)
  !! \retval sh    Setting wind direction (rad)
  !! 
  !! \param ltime  Passed to planet_position(). If .true., include light time, doubling the CPU time while gaining a bit of accur.
  !!
  !! \note
  !! - for sa0 = 0.d0, rise and set times are computed
  !! - for (sa0.ne.0), the routine calculates when alt=sa0 is reached
  !! - use sa0=-6,-12,-18 for the sun for twilight calculations (sa0 is expressed in degrees).
  !! - Returns times, rise/set azimuth and transit altitude
  !! - Moon transit error ~0.15s (?)
  !!
  !! \todo
  !! - This version sometimes finds the answer m(i)>1, which is the answer for the next day...
  !!   - (only?) solution: use a solver on JD for az,alt
  !!
  !! \note Speed wrt riset_ipol(), transit only, when using low-accuracy approximations:
  !! - Moon: 270x faster (used to be slowest by factor 3-5 wrt planets, factor 6.5 wrt Sun)
  !! - Sun:  340x faster - full accuracy (VSOP): ~2-3x SLOWER!
  !! - Mer:   20% slower
  !! - Ven:    4% slower
  !! - Mar:   10% slower
  !! - Jup:   25% slower
  !! - Nep:   30% slower
  !!
  !! \see
  !! - riset_ipol()
  !! - Meeus, Astronomical algorithms, Ch.15, but with geographic longitude east of Greenwich defined as > 0
  
  
  subroutine riset(jd,pl,  rt,tt,st, rh,ta,sh,  sa0, ltime)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi,pi2, d2r,am2r, enpname
    use SUFR_angles, only: rev
    use SUFR_date_and_time, only: cal2jd,jd2cal
    use SUFR_numerics, only: deq
    
    use TheSky_planets, only: planet_position, planet_position_la
    use TheSky_planetdata, only: planpos
    use TheSky_local, only: tz, lat0,lon0,deltat
    
    implicit none
    integer, intent(in) :: pl
    real(double), intent(in) :: jd, sa0
    real(double), intent(out) :: rt,tt,st,rh,ta,sh
    logical, intent(in), optional :: ltime
    
    integer :: mi,mj,yr,mnt,tc,mmax
    real(double) :: dy,day0,  jd0,jd1,m(3),  ra,dec
    real(double) :: sa,ch0,h0,agst,th0,dm,corr(3),accur,  ha,alt,azalt(3)
    character :: event(3)*(13)
    logical :: use_vsop, lltime
    
    lltime = .false.                    ! Call planet_position() IGNORING light time by default (faster, lower accuracy)
    if(present(ltime)) lltime = ltime
    
    ! Use the old interpolation routine for all but Moon and Sun:
    if(pl.ne.0.and.pl.ne.3) then
       call riset_ipol(jd,pl, rt,tt,st, rh,ta,sh, sa0, lltime)
       return
    end if
    
    
    alt = 0.d0;  ha = 0.d0;  h0 = 0.d0;  azalt = 0.d0;  m = 0.d0;  corr = 0.d0
    tc = 0        ! 0: geocentric, 1: topocentric, seems to give wrong results (for the Moon), see also different sa
    event = (/'Transit time ','Rise time    ','Set time     '/)
    
    sa = -0.5667d0*d2r                      ! Standard altitude for planets
    if(pl.eq.3) sa = sa - 16.d0*am2r        ! Compensate for radius of Sun, -16 arcminutes  ! sa = -0.8333d0*d2r - default for Sun
    if(pl.eq.0) sa = 0.125d0*d2r            ! Approximate standard altitude for Moon
    
    if(abs(sa0).gt.1.d-9) sa = sa0*d2r
    
    mmax = 3                   ! Maximum m: transit only: mmax=1, +rise/set: mmax=3
    if(sa0.gt.90.d0) mmax = 1  ! Compute only transit
    
    call jd2cal(jd, yr,mnt,dy)
    
    day0 = dble(int(dy))-tz/24.d0
    jd0  = cal2jd(yr,mnt,day0)
    
    call planet_position_la(jd0, pl, 3, 60)  ! Compute low-accuracy positions - calc=2 computes ra,dec, calc=3 computes agst
    
    ra   = planpos(5+tc*20)
    dec  = planpos(6+tc*20)
    agst = planpos(45)
    
    if(pl.eq.0) then
       if(tc.eq.0) then
          sa = asin(4.26345d-5/planpos(4))*0.7275d0-0.5667d0*d2r  ! Exact altitude for Moon, ~0.1-0.2deg, for geocentric coordinates
       else  !tc.eq.1:
          sa = -0.8333d0*d2r  ! For Moon, in combination with topocentric coordinates
       end if
    end if
    
    
    if(mmax.eq.3) then  ! Computing transit, rise and set
       ch0 = (sin(sa)-sin(lat0)*sin(dec)) / (cos(lat0)*cos(dec))
       if(abs(ch0).gt.1.d0) then  ! Body never rises/sets
          mmax = 1                ! Compute transit time and altitude only
       else
          h0 = rev(2*acos(ch0))/2.d0
       end if
    end if
    
    
    m(1) = rev(ra - lon0 - agst)/pi2 + corr(1)  ! Transit time in days; m(1)=m0, m(2)=m1, m(3)=m2 in Meeus, but lon0 > 0 for E
    if(mmax.eq.3) then
       m(2) = rev(m(1)*pi2 - h0)/pi2 + corr(2)  ! Rise time in days
       m(3) = rev(m(1)*pi2 + h0)/pi2 + corr(3)  ! Set time in days
    end if
    
    
    do mi=1,mmax
       mj = 0
       accur = 1.d-4       ! Accuracy.  Initially 1d-4, later 1d-6 ~ 0.1s. Don't make this smaller than 1d-16
       use_vsop = .false.  ! Initially
       
       dm = huge(dm)
       do while(abs(dm).ge.accur .or. .not.use_vsop)
          th0 = agst + 6.300388092591991d0*m(mi)  ! Meeus, p.103
          jd1 = jd0 + m(mi) + deltat/86400.d0
          
          if(abs(dm).le.accur) then
             use_vsop = .true.
             accur = 1.d-6  ! Changing this to 1.d-5 (~1s) speeds the code yearly_moon_table code up by ~30%
          end if
          
          if(use_vsop) then
             ! Accuracy of 1 min = 0.25 deg = 4e-3 rad.  1.d-6,1.d-2 saves ~43% CPU time for the Sun, <1 round-off errors/year
             call planet_position(jd1,pl, LBaccur=1.d-6,Raccur=1.d-2, ltime=lltime)  ! Uses truncated VSOP87
          else
             call planet_position_la(jd1,pl, 2,60)  ! Computes low-accuracy positions - calc=2 computes ra,dec - use 60 terms
          end if
          
          ra  = planpos(5+tc*20)  ! Right ascension
          dec = planpos(6+tc*20)  ! Declination
          
          ha  = rev(th0 + lon0 - ra + pi) - pi                         ! Hour angle
          alt = asin(sin(lat0)*sin(dec) + cos(lat0)*cos(dec)*cos(ha))  ! Altitude;  Meeus, Eq.13.6
          
          ! Correction to transit/rise/set times:
          if(mi.eq.1) then  ! Transit
             dm = -ha/pi2
          else              ! Rise/set
             dm = (alt-sa)/(pi2*cos(dec)*cos(lat0)*sin(ha))
          end if
          m(mi) = m(mi) + dm
          
          mj = mj+1
          if(mj.gt.30) exit  ! do-while loop
       end do  ! do while(abs(dm).ge.accur .or. .not.use_vsop)
       
       
       if(mj.gt.30) then  ! Convergence failed
          if(pl.ne.3.or.nint(sa0).ne.-18) write(0,'(A,F10.3,A)') '  * WARNING:  riset():  Riset failed to converge: '// &
               trim(enpname(pl))//'  '//trim(event(mi)),sa0,'d *'
          
          m(mi) = 0.d0
          azalt(mi) = 0.d0
       else               ! Result converged, store it
          if(mi.eq.1) then
             azalt(mi) = alt  ! Transit altitude
          else
             azalt(mi) = atan2(sin(ha),(cos(ha)*sin(lat0)-tan(dec)*cos(lat0)))  ! Rise,set hour angle -> azimuth
          end if
       end if
       
       if(pl.eq.0 .and. m(mi).gt.1.d0) then  ! Moon;  in case after m_i = m_i+1, m_f > 1
          m(mi) = 0.d0
          azalt(mi) = 0.d0
       end if
       
       if(m(mi).lt.0.d0 .and. deq(sa0,0.d0)) then
          m(mi) = 0.d0
          azalt(mi) = 0.d0
       end if
       
    end do  ! mi
    
    
    ! Store results:
    rt = m(2)*24  ! Rise time - days -> hours
    tt = m(1)*24  ! Transit time - days -> hours
    st = m(3)*24  ! Set time - days -> hours
    
    rh = azalt(2)  ! Rise azimuth
    ta = azalt(1)  ! Transit altitude
    sh = azalt(3)  ! Set azimuth
    
  end subroutine riset
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Old routine for rise, transit and set times routine for planets and asteroids - don't use this for Sun and Moon
  !!         Computes three sets of planet positions, and interpolates
  !!
  !! \param jd    Julian day number
  !! \param pl    Planet/object number  (-1 - 9: ES, Moon, Mer-Plu; >10000 asteroids)
  !! \param sa0   Altitude to return rise/set data for (degrees; 0. is actual rise/set).  sa0>90: compute transit only
  !!
  !! \retval rt   Rise time (hours)
  !! \retval tt   Transit time (hours)
  !! \retval st   Set time (hours)
  !! \retval rh   Rising wind direction (rad)
  !! \retval ta   Transit altitude (rad)
  !! \retval sh   Setting wind direction (rad)
  !! 
  !! \param ltime  Passed to planet_position(). If .true., include light time, doubling the CPU time while gaining a bit of accur.
  !! 
  !! \note
  !! - for sa0 = 0.d0, rise and set times are computed
  !! - for (sa0.ne.0), the routine calculates when alt=sa0 is reached
  !! - use sa0=-6,-12,-18 for the sun for twilight calculations (sa0 is expressed in degrees).
  !! - transit error for Moon <= ~10s (due to interpolation); typically ~4s?
  !! - code gets confused around midnight - it sometimes returns a small negative m(i), which is on the day before(!)
  !!
  !! \see
  !! - Meeus, Astronomical algorithms, Ch.15, but with geographic longitude east of Greenwich defined as > 0
  
  subroutine riset_ipol(jd,pl, rt,tt,st, rh,ta,sh, sa0, ltime)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi,pi2, d2r,am2r, enpname
    use SUFR_system, only: warn
    use SUFR_angles, only: rev
    use SUFR_date_and_time, only: cal2jd,jd2cal
    use SUFR_numerics, only: deq
    
    use TheSky_planets, only: planet_position
    use TheSky_planetdata, only: planpos
    use TheSky_local, only: tz, lat0,lon0,deltat
    
    implicit none
    integer, intent(in) :: pl
    real(double), intent(in) :: jd, sa0
    real(double), intent(out) :: rt,tt,st, rh,ta,sh
    logical, intent(in) :: ltime
    
    integer :: mi,mj,yr,mnt,tc,mmax, indic
    real(double) :: dy,day0,  jd0,jd1,jd2,m(3)
    real(double) :: ra0,dec0,ra1,dec1,ra2,dec2,ra,dec, sa,ch0,h0,agst,th0,n,dm,corr(3),accur,  ha,alt,azalt(3)
    character :: event(3)*(13)
    save :: indic
    
    if((pl.eq.0.or.pl.eq.3) .and. indic.ne.12345) then
       call warn("riset_ipol():  Don't use this routine for Sun or Moon - use riset() instead.", 0)
       indic = 12345
    end if
    
    rt=0.d0; tt=0.d0; st=0.d0;  rh=0.d0; ta=0.d0; sh=0.d0
    alt = 0.d0;  ha = 0.d0;  corr = 0.d0;  h0 = 0.d0;  dec = 0.d0
    tc = 0        ! 0: geocentric, 1: topocentric, seems to give wrong results (for the Moon), see also different sa
    event = (/'Transit time ','Rise time    ','Set time     '/)
    accur = 1.d-6          ! Accuracy.  1d-6 ~ 0.1s. Don't make this smaller than 1d-16
    
    
    sa = -0.5667d0*d2r                      ! Standard altitude for planets
    if(pl.eq.3) sa = sa - 16.0d0*am2r       ! Compensate for radius of Sun, -16 arcminutes  ! sa = -0.8333d0*d2r - default for Sun
    if(pl.eq.0) sa = 0.125d0*d2r            ! Approximate standard altitude for Moon
    
    if(abs(sa0).gt.1.d-9) sa = sa0*d2r
    
    call jd2cal(jd,yr,mnt,dy)
    
    day0 = dble(int(dy))-tz/24.d0
    jd0  = cal2jd(yr,mnt,day0-1.d0)
    jd1  = cal2jd(yr,mnt,day0)
    jd2  = cal2jd(yr,mnt,day0+1.d0)
    
    call planet_position(jd0,pl, ltime=ltime)
    ra0  = planpos(5+tc*20)
    dec0 = planpos(6+tc*20)
    
    call planet_position(jd1,pl, ltime=ltime)
    ! Exact altitude for moon, ~0.1-0.2deg, for geocentric coordinates:
    if(tc.eq.0.and.pl.eq.0) sa = asin(4.26345d-5/planpos(4))*0.7275d0-0.5667d0*d2r
    
    ! For Moon, in combination with topocentric coordinates:
    if(tc.eq.1.and.pl.eq.0) sa = -0.8333d0*d2r
    ra1  = planpos(5+tc*20)
    dec1 = planpos(6+tc*20)
    agst = planpos(45)
    
    call planet_position(jd2,pl, ltime=ltime)
    ra2  = planpos(5+tc*20)
    dec2 = planpos(6+tc*20)
    
    if(abs(ra1-ra0).gt.2..or.abs(ra2-ra1).gt.2.) then  ! ra flips 2pi <-> 0
       if(ra0.ge.4.5.and.ra1.ge.4.5.and.ra2.le.2.) ra2 = ra2+2*pi
       if(ra0.ge.4.5.and.ra1.le.2.and.ra2.le.2.)   ra0 = ra0-2*pi
       if(ra0.le.2.and.ra1.ge.4.5.and.ra2.ge.4.5)  ra0 = ra0+2*pi
       if(ra0.le.4.5.and.ra1.le.2.and.ra2.ge.4.5)  ra2 = ra2-2*pi
       !write(6,'(A)') 'RA flipped'
    end if
    
    ch0 = (sin(sa)-sin(lat0)*sin(dec2)) / (cos(lat0)*cos(dec2))  ! Meeus, Eq.15.1
    
    m = 0.
    azalt = 0.
    mmax = 3
    if(abs(ch0).gt.1.d0) then  ! Body never rises/sets
       mmax = 1                ! Compute transit time and altitude only
    else
       h0 = rev(2*acos(ch0))/2.d0  ! cos(H0) -> H0
    end if
    
    
    m(1) = rev(ra2 - lon0 - agst)/pi2 + corr(1)  ! Transit time in days;  Meeus Eq. 15.2a, but lon0 > 0 for E
    if(mmax.eq.3) then
       m(2) = rev(m(1)*pi2 - h0)/pi2 + corr(2)   ! Rise time in days;  Meeus Eq.15.2b
       m(3) = rev(m(1)*pi2 + h0)/pi2 + corr(3)   ! Set time in days;  Meeus Eq.15.2c
    end if
    
    mi = 1  ! i=1,mmax
    do while(mi.le.mmax)
       mj = 0
       
       dm = huge(dm)
       do while(abs(dm).ge.accur)
          th0 = agst + 6.300388092591991d0*m(mi)  ! Meeus, p.103
          n   = m(mi) + deltat/86400.d0
          ra  = mipol(ra0,ra1,ra2,n)     ! Interpolate right ascension
          dec = mipol(dec0,dec1,dec2,n)  ! Interpolate declination
          
          ha = rev(th0 + lon0 - ra + pi) - pi                          ! Hour angle;  Meeus p.103
          alt = asin(sin(lat0)*sin(dec) + cos(lat0)*cos(dec)*cos(ha))  ! Meeus, Eq.13.6
          
          ! Correction to transit/rise/set times:
          if(mi.eq.1) then  ! Transit
             dm = -ha/pi2
          else              ! Rise/set
             dm = (alt-sa)/(pi2*cos(dec)*cos(lat0)*sin(ha))
          end if
          m(mi) = m(mi) + dm
          
          mj = mj+1
          if(mj.gt.30) exit  ! do-while loop
       end do  ! do while(abs(dm).ge.accur)
       
       
       if(mj.gt.30) then  ! Convergence failed
          if(pl.ne.3.or.nint(sa0).ne.-18) write(0,'(A,F10.3,A)') '  * WARNING:  riset_ipol():  Riset failed to converge: '// &
               trim(enpname(min(pl,19)))//'  '//trim(event(mi)),sa0,'d'
          
          m(mi) = 0.d0
          azalt(mi) = 0.d0
       else               ! Result converged, store it
          if(mi.eq.1) then
             azalt(mi) = alt  ! Transit altitude
          else
             azalt(mi) = atan2(sin(ha),(cos(ha)*sin(lat0)-tan(dec)*cos(lat0)))  ! Rise,set hour angle -> azimuth
          end if
       end if
       
       
       if(pl.eq.0) then
          if(m(mi).lt.0.) then    ! Sometimes, (at least) the Moon doesn't r,t,s twice, because m_i slightly > 1, and rev(m_i) << 1,
             corr(mi) = corr(mi) + 1.d0     ! then initially m_f < 0, though if m_i = m_i + 1,  then 0 < m_f < 1
             
             ! Restart the do loop:
             mi = 1
             m(1) = rev(ra2 - lon0 - agst)/pi2 + corr(1)  ! Transit; m(1)=m0, m(2)=m1, m(3)=m2 in Meeus, but lon0 > 0 for E
             if(mmax.eq.3) then
                m(2) = rev(m(1)*pi2 - h0)/pi2 + corr(2)   ! Rise
                m(3) = rev(m(1)*pi2 + h0)/pi2 + corr(3)   ! Set
             end if
             cycle  ! i
          end if  ! if(m(mi).lt.0.)
          
          if(m(mi).gt.1.d0) then  ! in case after m_i = m_i+1, m_f > 1
             m(mi) = 0.d0
             azalt(mi) = 0.d0
          end if  ! if(m(mi).gt.1.d0)
       end if  ! if(pl.eq.0)
       
       
       if(m(mi).lt.0. .and. deq(sa0,0.d0)) then
          m(mi) = 0.d0
          azalt(mi) = 0.d0
       end if
       
       mi = mi+1  ! mi=1,mmax
    end do  ! do while(mi.le.mmax)
    
    
    ! Store results:
    rt = m(2)*24  ! Rise time - days -> hours
    tt = m(1)*24  ! Transit time - days -> hours
    st = m(3)*24  ! Set time - days -> hours
    
    rh = azalt(2)  ! Rise azimuth
    ta = azalt(1)  ! Transit altitude
    sh = azalt(3)  ! Set azimuth
    
  end subroutine riset_ipol
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief Rise, transit and set times routine for an object with fixed ra & dec
  !!
  !! \param jd    Julian day number
  !! \param ra    Right ascension (rad)
  !! \param dec   Declination (rad)
  !! \param sa0   Altitude to return rise/set data for (degrees; 0. is actual rise/set).  sa0>90: compute transit only
  !!
  !! \retval rt   Rise time (hours)
  !! \retval tt   Transit time (hours)
  !! \retval st   Set time (hours)
  !! \retval rh   Rising wind direction (rad)
  !! \retval ta   Transit altitude (rad)
  !! \retval sh   Setting wind direction (rad)
  !!
  !! for sa0 = 0., rise and set times are computed
  !! for sa0.ne.0, the routine calculates when alt=sa0 is reached
  !! use sa0=-6,-12,-18 for the sun for twilight calculations
  !! (sa0 is expressed in degrees).
  !!
  !! \see
  !! - Meeus, Astronomical algorithms, Ch.15, but with geographic longitude east of Greenwich defined as > 0
  
  subroutine riset_ad(jd, ra,dec, rt,tt,st,rh,ta,sh, sa0)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi2, d2r
    use SUFR_angles, only: rev, rev2
    use SUFR_date_and_time, only: cal2jd, jd2cal
    use SUFR_numerics, only: deq
    
    use TheSky_planets, only: planet_position
    use TheSky_planetdata, only: planpos
    use TheSky_local, only: tz, lat0,lon0 !,deltat
    
    implicit none
    real(double), intent(in) :: jd, ra,dec, sa0
    real(double), intent(out) :: rt,tt,st,rh,ta,sh
    integer :: mi,mj,yr,mnt,mmax
    real(double) :: dy,day0,jd1,m(3),  sa,ch0,h0,agst,th0,dm,accur,  ha,alt,azalt(3) !,n
    character :: event(3)*(13)
    
    alt = 0.d0;  ha = 0.d0;  h0 = 0.d0
    
    accur = 1.d-6  ! Accuracy.  1d-6 ~ 0.1s. Don't make this smaller than 1d-16
    
    event = (/'Transit time ','Rise time    ','Set time     '/)
    
    sa = sa0*d2r
    mmax = 3
    if(sa0.gt.90.d0) mmax = 1
    
    call jd2cal(jd,yr,mnt,dy)
    day0 = dble(int(dy))-tz/24.d0
    jd1  = cal2jd(yr,mnt,day0)
    call planet_position(jd1,3)
    agst = planpos(45)
    
    m = 0.
    azalt = 0.
    if(mmax.eq.3) then
       ch0 = (sin(sa)-sin(lat0)*sin(dec))/(cos(lat0)*cos(dec))
       if(abs(ch0).gt.1.d0) then  ! Body never rises/sets
          mmax = 1                ! Compute transit time and altitude only
       else
          h0 = rev(2*acos(ch0))/2.d0
       end if
    end if
    
    m(1) = rev(ra - lon0 - agst)/pi2  ! Transit;  m(1)=m0, m(2)=m1, m(3)=m2 in Meeus, but lon0 > 0 for E
    if(mmax.eq.3) then
       m(2) = rev(m(1)*pi2 - h0)/pi2  ! Rise
       m(3) = rev(m(1)*pi2 + h0)/pi2  ! Set
    end if
    
    
    do mi=1,mmax
       mj = 0
       dm = huge(dm)
       do while(abs(dm).ge.accur)
          th0 = agst + 6.300388092591991d0*m(mi)  ! Meeus, p.103
          
          ha = rev2(th0 + lon0 - ra)                                   ! Hour angle;  Meeus p.103
          alt = asin(sin(lat0)*sin(dec) + cos(lat0)*cos(dec)*cos(ha))  ! Altitude;  Meeus Eq. 13.6
          
          ! Correction on transit/rise/set times:
          if(mi.eq.1) then  ! Transit
             dm = -ha/pi2
          else              ! Rise/set
             dm = (alt-sa)/(pi2*cos(dec)*cos(lat0)*sin(ha))
          end if
          m(mi) = m(mi) + dm
          
          mj = mj+1
          if(mj.gt.30) exit  ! do-while loop
       end do  ! do while(abs(dm).ge.accur)
       
       
       if(mj.gt.30) then  ! Convergence failed
          write(0,'(A,F10.3,A)') '  * WARNING:  riset_ad():  Riset failed to converge: '//trim(event(mi)),sa0,'d *'
          
          m(mi) = 0.d0
          azalt(mi) = 0.d0
       else               ! Result converged, store it
          if(mi.eq.1) then
             azalt(mi) = alt  ! Transit altitude
          else
             azalt(mi) = atan2(sin(ha),(cos(ha)*sin(lat0)-tan(dec)*cos(lat0)))  ! Rise,set hour angle -> azimuth
          end if
       end if
       
       if(m(mi).lt.0. .and. deq(sa0,0.d0)) then
          m(mi) = 0.d0
          azalt(mi) = 0.d0
       end if
       
    end do  ! mi
    
    
    ! Store results:
    rt = m(2)*24  ! Rise time - days -> hours
    tt = m(1)*24  ! Transit time - days -> hours
    st = m(3)*24  ! Set time - days -> hours
    
    rh = azalt(2)  ! Rise azimuth
    ta = azalt(1)  ! Transit altitude
    sh = azalt(3)  ! Set azimuth
    
  end subroutine riset_ad
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief Compute the date (m,d) at which an object with ra can be observed best, i.e., it transits at midnight  
  !!
  !! \param  jd00  Julian day?
  !! \param  ra0   Right ascension
  !!
  !! \retval m     Month of best observation
  !! \retval d     Day of month of best observation
  !!
  !! \note 
  !! - VERY SLOW!!! -> use best_obs_date_ra instead!!!  (but not (yet) in libTheSky...)
  !! - riset() should be faster now (04/2010)
  
  subroutine best_obs_date(jd00,ra0, m,d)
    use SUFR_kinds, only: double
    use SUFR_constants, only: r2h
    use SUFR_angles, only: rev,rv12
    use SUFR_date_and_time, only: jd2cal
    use SUFR_numerics, only: deq
    
    use TheSky_planets, only: planet_position
    use TheSky_planetdata, only: planpos, nplanpos
    use TheSky_local, only: tz
    use TheSky_datetime, only: gettz
    
    implicit none
    real(double), intent(in) :: jd00,ra0
    integer, intent(out) :: m,d
    real(double) :: jd0,oldjd,z(nplanpos)
    real(double) :: rt,stt,ott,st,rh,ta,sh,midnight,hh,dt,dd
    integer :: y,iter
    
    oldjd = -99.d0
    jd0 = jd00
    
    dd = 9999.d0
    iter = 0
    do while(abs(dd).gt.0.5d0)
       call riset(jd0,3, rt,stt,st,rh,ta,sh, 0.d0)  ! Sun transit time
       midnight = mod(stt+36.d0,24.d0)
       
       call planet_position(jd0,3)                  ! Position of the Sun
       z = planpos
       
       hh  = rev(z(44) - ra0)                       ! Hour angle is LST - RA
       ott = -hh*r2h                                ! Time of object transit = -ha
       dt  = rv12(ott-midnight)                     ! Difference between transit time and midnight on this day
       dd  = dt/24.d0*365.25d0                      ! Estimate for difference in days between today and day we're looking for
       jd0 = dble(nint(jd0 + dd + 0.5d0)) - 0.5d0   ! Force midnight UT
       tz  = gettz(jd0)
       jd0 = jd0 - tz/24.d0                         ! Local midnight
       
       ! Consider result 'converged':
       if(deq(jd0,oldjd)) exit
       if(iter.gt.5 .and. abs(dd).lt.0.7d0) exit
       if(iter.gt.10 .and. abs(dd).lt.1.d0) exit
       if(iter.gt.100) exit
       
       oldjd = jd0
       iter = iter + 1
    end do
    
    call jd2cal(jd0, y,m,dd)
    d = floor(dd)
    
  end subroutine best_obs_date
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Interpolate m for riset
  
  function mipol(y1,y2,y3, n)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: y1,y2,y3,n
    real(double) :: mipol, a,b,c
    
    a = y2-y1
    b = y3-y2
    c = b-a
    mipol = y2 + n/2.d0 * (a + b + c*n)
    
  end function mipol
  !*********************************************************************************************************************************
  
  
end module TheSky_riset
!***********************************************************************************************************************************
