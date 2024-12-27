!> \file riset.f90  Compute rise, transit and set times, or beginning and end of twilight for libTheSky


!  Copyright (c) 2002-2023  Marc van der Sluys - marc.vandersluys.nl
!   
!  This file is part of the libTheSky package, 
!  see: https://libthesky.sf.net/
!  
!  This is free software: you can redistribute it and/or modify it under the terms of the
!  European Union Public Licence 1.2 (EUPL 1.2).
!  
!  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!  See the EU Public Licence for more details.
!  
!  You should have received a copy of the European Union Public Licence along with this code.
!  If not, see <https://www.eupl.eu/1.2/en/>.



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
  !! \param  jd        Julian day number
  !! \param  pl        Planet/object number  (-1 - 9: ES, Moon, Mer-Plu; >10000 asteroids)
  !! \param  rsAlt     Altitude to return rise/set data for (degrees; 0. is actual rise/set).  rsAlt>90: compute transit only
  !! 
  !! \param rt        Rise time (hours) (output)
  !! \param tt        Transit time (hours) (output)
  !! \param st        Set time (hours) (output)
  !! 
  !! \param rh        Rising wind direction (rad) (output)
  !! \param ta        Transit altitude (rad) (output)
  !! \param sh        Setting wind direction (rad) (output)
  !! 
  !! \param  ltime     Passed to planet_position(). If .true., include light time, doubling the CPU time while gaining a bit of 
  !!                     accuracy (optional; default: false)
  !! \param  cWarn     Warn upon convergence failure (optional; default: true)
  !! 
  !! \param converge  Number of iterations needed to converge (optional) (output)
  !!
  !! \note
  !! - for rsAlt = 0.d0, rise and set times are computed
  !! - for (rsAlt.ne.0), the routine calculates when alt=rsAlt is reached
  !! - use rsAlt=-6,-12,-18 for the sun for twilight calculations (rsAlt is expressed in degrees).
  !! - Returns times, rise/set azimuth and transit altitude
  !! - Moon transit error ~0.15s (?)
  !! 
  !! \todo
  !! - This version sometimes finds the answer tmRad(i)>1, which is the answer for the next day...
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
  
  
  subroutine riset(jd,pl,  rt,tt,st, rh,ta,sh,  rsAlt, ltime, cWarn, converge)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi2, d2r,am2r, r2h, enpname, earthr,AU
    use SUFR_angles, only: rev, rev2
    use SUFR_date_and_time, only: cal2jd,jd2cal
    use SUFR_numerics, only: deq0
    
    use TheSky_local, only: tz, lat0,lon0,deltat
    use TheSky_planets, only: planet_position, planet_position_la
    use TheSky_planetdata, only: planpos
    use TheSky_coordinates, only: refract
    
    implicit none
    integer, intent(in) :: pl
    real(double), intent(in) :: jd, rsAlt
    real(double), intent(out) :: rt,tt,st,rh,ta,sh
    logical, intent(in), optional :: ltime, cWarn
    integer, intent(out), optional :: converge(3)
    
    integer :: evi,iter,yr,mnt,tc,evMax, lconverge(3)
    real(double) :: dy,day0,  jd0,jd1,tmRad(3),  ra,dec
    real(double) :: rsa,cosH0,h0,agst0,th0,dTmRad,accur,  ha,alt,azAlt(3)
    character :: event(3)*(13)
    logical :: use_vsop, lltime, lcWarn
    
    lltime = .false.                    ! Call planet_position() IGNORING light time by default (faster, lower accuracy)
    if(present(ltime)) lltime = ltime
    
    lcWarn = .true.                     ! Warn upon convergence failure
    if(present(cWarn)) lcWarn = cWarn
    
    ! Use the old interpolation routine for all but Moon and Sun:
    if(pl.ne.0.and.pl.ne.3) then  ! A planet
       call riset_ipol(jd,pl, rt,tt,st, rh,ta,sh, rsAlt, lltime, lcWarn)
       return
    end if
    
    
    alt=0.d0; ha=0.d0; h0=0.d0; azAlt=0.d0; tmRad=0.d0
    tc = 1        ! 0: geocentric, 1: topocentric;  see also different rsa for the Moon below
    event = ['Transit time ','Rise time    ','Set time     ']
    
    if(abs(rsAlt).gt.1.d-9) then
       rsa = rsAlt*d2r  ! Use a user-specified altitude
    else
       rsa = -0.5667d0*d2r                      ! Standard rise/set altitude for planets
       if(pl.eq.3) rsa = rsa - 16.d0*am2r       ! Correct for radius of Sun, -16 arcminutes  ! rsa = -0.8333d0*d2r - default for Sun
       if(pl.eq.0) rsa = 0.125d0*d2r            ! Approximate standard rise/set altitude for Moon, including parallax
    end if
    
    evMax = 3                     ! 'Maximum' event to compute: transit only: evMax=1, +rise/set: evMax=3
    if(rsAlt.gt.90.d0) evMax = 1  ! Compute transit time and altitude only
    
    call jd2cal(jd, yr,mnt,dy)
    
    day0 = dble(int(dy))-tz/24.d0  ! Midnight local time, needed for agst0
    jd0  = cal2jd(yr,mnt,day0)
    
    call planet_position_la(jd0, pl, 3, 60)  ! Compute low-accuracy positions - calc=2 computes ra,dec, calc=3 computes AGST
    
    ra    = planpos(5+tc*20)
    dec   = planpos(6+tc*20)
    agst0 = planpos(45)       ! AGST for midnight
    
    if(pl.eq.0) then  ! Moon
       if(tc.eq.0) then  ! Geocentric
          rsa = asin(earthr/(planpos(4)*AU))*0.7275d0-0.5667d0*d2r  ! Exact altitude for Moon, ~0.1-0.2deg, for geocentric coord.
       else  ! tc.eq.1 - topocentric:
          rsa = -0.8333d0*d2r  ! For Moon, in combination with topocentric coordinates
       end if
    end if
    
    
    if(evMax.eq.3) then  ! Computing transit, rise and set
       cosH0 = (sin(rsa)-sin(lat0)*sin(dec)) / (cos(lat0)*cos(dec))  ! Cosine of the hour angle of rise/set; Meeus, Eq.15.1
       if(abs(cosH0).gt.1.d0) then  ! Body never rises/sets
          evMax = 1                 ! Compute transit time and altitude only
       else
          h0 = rev(2*acos(cosH0))/2.d0  ! Hour angle of rise/set
       end if
    end if
    
    
    tmRad(1) = rev(ra - lon0 - agst0)  ! Transit time in radians; tmRad(1)=m0, tmRad(2)=m1, tmRad(3)=m2 in Meeus, but lon0 > 0 for E
    if(evMax.eq.3) then
       tmRad(2) = rev(tmRad(1) - h0)   ! Rise time in radians
       tmRad(3) = rev(tmRad(1) + h0)   ! Set time in radians
    end if
    
    
    do evi=1,evMax         ! Transit, rise, set
       iter = 0
       accur = 1.d-3       ! Accuracy.  Initially 1d-3, later 1d-5 ~ 0.1s. Don't make this smaller than 1d-16
       use_vsop = .false.  ! Initially
       
       dTmRad = huge(dTmRad)
       do while(abs(dTmRad).ge.accur .or. .not.use_vsop)
          th0 = agst0 + 1.002737909350795d0*tmRad(evi)   ! Solar day in sidereal days in 2000; Expl.Suppl.tt.Astr.Almanac 3rdEd 
          jd1 = jd0 + tmRad(evi)/pi2 + deltat/86400.d0   !                                    Eq.3.17 (removed '...37...' typo)
          
          if(abs(dTmRad).le.accur) then
             use_vsop = .true.
             accur = 1.d-5  ! 1d-5~0.1s.  Changing this to 1.d-4 (~1s) speeds the code yearly_moon_table code up by ~30%.  Don't make this smaller than 1d-16.
          end if
          
          if(use_vsop) then
             ! Accuracy of 1 min = 0.25 deg = 4e-3 rad.  1.d-6,1.d-2 saves ~43% CPU time for the Sun, <1 round-off errors/year
             call planet_position(jd1,pl, LBaccur=1.d-6,Raccur=1.d-2, ltime=lltime)  ! Uses truncated VSOP87
          else
            call planet_position_la(jd1,pl, 2,60)  ! Computes low-accuracy positions - calc=2 computes ra,dec - use 60 terms
          end if
          
          ra  = planpos(5+tc*20)  ! Right ascension
          dec = planpos(6+tc*20)  ! Declination
          
          ha  = rev2(th0 + lon0 - ra)                                  ! Hour angle
          ! ha  = planpos(8+tc*20)                                       ! Hour angle -/- DeltaT
          alt = asin(sin(lat0)*sin(dec) + cos(lat0)*cos(dec)*cos(ha))  ! Altitude;  Meeus, Eq.13.6
          
          ! Correction to transit/rise/set times:
          if(evi.eq.1) then  ! Transit
             dTmRad = -ha
          else               ! Rise/set
             dTmRad = (alt-rsa)/(cos(dec)*cos(lat0)*sin(ha))
          end if
          tmRad(evi) = tmRad(evi) + dTmRad
          
          iter = iter+1
          if(iter.gt.30) exit  ! do-while loop
       end do  ! do while(abs(dTmRad).ge.accur .or. .not.use_vsop)
       
       
       if(iter.gt.30) then  ! Convergence failed
          if(lcWarn .and. (pl.ne.3.or.nint(rsAlt).ne.-18)) write(0,'(A,F10.3,A)') '  * WARNING:  riset():  Riset failed to '// &
               'converge: '//trim(enpname(pl))//'  '//trim(event(evi)),rsAlt,'d *'
          
          tmRad(evi) = 0.d0
          azAlt(evi) = 0.d0
       else                    ! Result converged, store it
          if(evi.eq.1) then
             azAlt(evi) = alt  ! Transit altitude
          else
             azAlt(evi) = atan2(sin(ha),(cos(ha)*sin(lat0)-tan(dec)*cos(lat0)))  ! Rise,set hour angle -> azimuth
          end if
       end if
       
       if(pl.eq.0 .and. tmRad(evi).gt.pi2) then  ! Moon;  in case after m_i = m_i+1, m_f > 1
          tmRad(evi) = 0.d0
          azAlt(evi) = 0.d0
       end if
       
       if(tmRad(evi).lt.0.d0 .and. deq0(rsAlt)) then
          tmRad(evi) = 0.d0
          azAlt(evi) = 0.d0
       end if
       
       lconverge(evi) = iter  ! Number of iterations needed to converge
       
    end do  ! evi
    
    
    ! Store results:
    tmRad = tmRad * r2h                ! Times radians -> hours
    
    tt = tmRad(1)                      ! Transit time
    rt = tmRad(2)                      ! Rise time
    st = tmRad(3)                      ! Set time
    
    ta = azAlt(1) + refract(azAlt(1))  ! Transit altitude + refraction
    rh = azAlt(2)                      ! Rise azimuth
    sh = azAlt(3)                      ! Set azimuth
    
    if(present(converge)) converge = lconverge  ! Number of iterations needed to converge
    
  end subroutine riset
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Old routine for rise, transit and set times routine for planets and asteroids - don't use this for Sun and Moon
  !!         Computes three sets of planet positions, and interpolates
  !!
  !! \param jd     Julian day number
  !! \param pl     Planet/object number  (-1 - 9: ES, Moon, Mer-Plu; >10000 asteroids)
  !! \param rsAlt  Altitude to return rise/set data for (degrees; 0. is actual rise/set).  rsAlt>90: compute transit only
  !!
  !! \param rt   Rise time (hours) (output)
  !! \param tt   Transit time (hours) (output)
  !! \param st   Set time (hours) (output)
  !! \param rh   Rising wind direction (rad) (output)
  !! \param ta   Transit altitude (rad) (output)
  !! \param sh   Setting wind direction (rad) (output)
  !! 
  !! \param ltime  Passed to planet_position(). If .true., include light time, doubling the CPU time while gaining a bit of accur.
  !! \param cWarn  Warn upon convergence failure (optional; default: true)
  !! 
  !! \note
  !! - for rsAlt = 0.d0, rise and set times are computed
  !! - for (rsAlt.ne.0), the routine calculates when alt=rsAlt is reached
  !! - use rsAlt=-6,-12,-18 for the sun for twilight calculations (rsAlt is expressed in degrees).
  !! - transit error for Moon <= ~10s (due to interpolation); typically ~4s?
  !! - code gets confused around midnight - it sometimes returns a small negative tmdy(i), which is on the day before(!)
  !!
  !! \see
  !! - Meeus, Astronomical algorithms, Ch.15, but with geographic longitude east of Greenwich defined as > 0
  
  subroutine riset_ipol(jd,pl, rt,tt,st, rh,ta,sh, rsAlt, ltime, cWarn)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi,pi2, d2r,am2r, enpname, earthr,AU
    use SUFR_system, only: warn
    use SUFR_angles, only: rev, rev2
    use SUFR_date_and_time, only: cal2jd,jd2cal
    use SUFR_numerics, only: deq0
    
    use TheSky_local, only: tz, lat0,lon0,deltat
    use TheSky_planets, only: planet_position
    use TheSky_planetdata, only: planpos
    use TheSky_coordinates, only: refract
    
    implicit none
    integer, intent(in) :: pl
    real(double), intent(in) :: jd, rsAlt
    real(double), intent(out) :: rt,tt,st, rh,ta,sh
    logical, intent(in) :: ltime
    logical, intent(in), optional :: cWarn
    
    integer :: evi,iter,yr,mnt,tc,evMax, indic
    real(double) :: dy,day0,  jd0,jd1,jd2,tmdy(3),dTmdy(3)
    real(double) :: ra0,dec0,ra1,dec1,ra2,dec2,ra,dec, rsa,cosH0,h0,agst0,th0,n,dtm,accur,  ha,alt,azAlt(3)
    character :: event(3)*(13)
    logical :: lcWarn
    save :: indic
    
    lcWarn = .true.
    if(present(cWarn)) lcWarn = cWarn
    
    if((pl.eq.0.or.pl.eq.3) .and. indic.ne.12345) then
       call warn("riset_ipol():  Don't use this routine for Sun or Moon - use riset() instead.", 0)
       indic = 12345
    end if
    
    rt=0.d0; tt=0.d0; st=0.d0;  rh=0.d0; ta=0.d0; sh=0.d0
    alt = 0.d0;  ha = 0.d0;  dTmdy = 0.d0;  h0 = 0.d0;  dec = 0.d0
    tc = 1        ! 0: geocentric, 1: topocentricl;  see also different rsa for the Moon below
    event = ['Transit time ','Rise time    ','Set time     ']
    accur = 1.d-6          ! Accuracy.  1d-6 ~ 0.1s. Don't make this smaller than 1d-16
    
    
    rsa = -0.5667d0*d2r                     ! Standard rise/set altitude for planets
    if(pl.eq.3) rsa = rsa - 16.0d0*am2r     ! Compensate for radius of Sun, -16 arcminutes  ! rsa = -0.8333d0*d2r - default for Sun
    if(pl.eq.0) rsa = 0.125d0*d2r           ! Approximate standard altitude for Moon
    
    if(abs(rsAlt).gt.1.d-9) rsa = rsAlt*d2r
    
    call jd2cal(jd,yr,mnt,dy)
    
    day0 = dble(int(dy))-tz/24.d0  ! Midnight local time, needed for agst0
    jd0  = cal2jd(yr,mnt,day0-1.d0)
    jd1  = cal2jd(yr,mnt,day0)
    jd2  = cal2jd(yr,mnt,day0+1.d0)
    
    call planet_position(jd0,pl, LBaccur=1.d-6,Raccur=1.d-2, ltime=ltime)
    ra0  = planpos(5+tc*20)
    dec0 = planpos(6+tc*20)
    
    call planet_position(jd1,pl, LBaccur=1.d-6,Raccur=1.d-2, ltime=ltime)
    ra1   = planpos(5+tc*20)
    dec1  = planpos(6+tc*20)
    agst0 = planpos(45)       ! AGST for midnight
    
    if(pl.eq.0) then
       if(tc.eq.0) then  ! Geocentric
          rsa = asin(earthr/(planpos(4)*AU))*0.7275d0-0.5667d0*d2r  ! Exact altitude for Moon, ~0.1-0.2deg, for geocentric coord.
       else  ! tc.eq.1 - topocentric:
          rsa = -0.8333d0*d2r  ! For Moon, in combination with topocentric coordinates
       end if
    end if
    
    call planet_position(jd2,pl, LBaccur=1.d-6,Raccur=1.d-2, ltime=ltime)
    ra2  = planpos(5+tc*20)
    dec2 = planpos(6+tc*20)
    
    if(abs(ra1-ra0).gt.2..or.abs(ra2-ra1).gt.2.) then  ! ra flips 2pi <-> 0
       if(ra0.ge.4.5.and.ra1.ge.4.5.and.ra2.le.2.) ra2 = ra2+2*pi
       if(ra0.ge.4.5.and.ra1.le.2.and.ra2.le.2.)   ra0 = ra0-2*pi
       if(ra0.le.2.and.ra1.ge.4.5.and.ra2.ge.4.5)  ra0 = ra0+2*pi
       if(ra0.le.4.5.and.ra1.le.2.and.ra2.ge.4.5)  ra2 = ra2-2*pi
       ! write(6,'(A)') 'RA flipped'
    end if
    
    cosH0 = (sin(rsa)-sin(lat0)*sin(dec2)) / (cos(lat0)*cos(dec2))  ! Cosine of the hour angle of rise/set; Meeus, Eq.15.1
    
    tmdy = 0.d0
    azAlt = 0.d0
    evMax = 3                    ! 'Maximum' event to compute: transit only: evMax=1, +rise/set: evMax=3
    if(abs(cosH0).gt.1.d0) then  ! Body never rises/sets
       evMax = 1                 ! Compute transit time and altitude only
    else
       h0 = rev(2*acos(cosH0))/2.d0  ! Hour angle of rise/set
    end if
    
    
    tmdy(1) = rev(ra2 - lon0 - agst0)/pi2 + dTmdy(1)    ! Transit time in days;  Meeus Eq. 15.2a, but lon0 > 0 for E
    if(evMax.eq.3) then
       tmdy(2) = rev(tmdy(1)*pi2 - h0)/pi2 + dTmdy(2)   ! Rise time in days;  Meeus Eq.15.2b
       tmdy(3) = rev(tmdy(1)*pi2 + h0)/pi2 + dTmdy(3)   ! Set time in days;  Meeus Eq.15.2c
    end if
    
    evi = 1  ! i=1,evMax
    do while(evi.le.evMax)
       iter = 0
       
       dtm = huge(dtm)
       do while(abs(dtm).ge.accur)
          th0 = agst0 + 6.300388092591991d0*tmdy(evi)  ! Meeus, p.103
          n   = tmdy(evi) + deltat/86400.d0
          ra  = rsIpol(ra0,ra1,ra2,n)     ! Interpolate right ascension
          dec = rsIpol(dec0,dec1,dec2,n)  ! Interpolate declination
          
          ha = rev2(th0 + lon0 - ra)                                   ! Hour angle;  Meeus p.103
          alt = asin(sin(lat0)*sin(dec) + cos(lat0)*cos(dec)*cos(ha))  ! Meeus, Eq.13.6
          
          ! Correction to transit/rise/set times:
          if(evi.eq.1) then  ! Transit
             dtm = -ha/pi2
          else              ! Rise/set
             dtm = (alt-rsa)/(pi2*cos(dec)*cos(lat0)*sin(ha))
          end if
          tmdy(evi) = tmdy(evi) + dtm
          
          iter = iter+1
          if(iter.gt.30) exit  ! do-while loop
       end do  ! do while(abs(dtm).ge.accur)
       
       
       if(iter.gt.30) then  ! Convergence failed
          if(lcWarn .and. (pl.ne.3.or.nint(rsAlt).ne.-18)) write(0,'(A,F10.3,A)') '  * WARNING:  riset_ipol():  Riset failed '// &
               'to converge: '//trim(enpname(min(pl,19)))//'  '//trim(event(evi)),rsAlt,'d'
          
          tmdy(evi) = 0.d0
          azAlt(evi) = 0.d0
       else               ! Result converged, store it
          if(evi.eq.1) then
             azAlt(evi) = alt  ! Transit altitude
          else
             azAlt(evi) = atan2(sin(ha),(cos(ha)*sin(lat0)-tan(dec)*cos(lat0)))  ! Rise,set hour angle -> azimuth
          end if
       end if
       
       
       if(pl.eq.0) then
          if(tmdy(evi).lt.0.) then    ! Sometimes, (at least) the Moon doesn't r,t,s twice, because m_i slightly > 1, and rev(m_i) << 1,
             dTmdy(evi) = dTmdy(evi) + 1.d0     ! then initially m_f < 0, though if m_i = m_i + 1,  then 0 < m_f < 1
             
             ! Restart the do loop:
             evi = 1
             tmdy(1) = rev(ra2 - lon0 - agst0)/pi2 + dTmdy(1)  ! Transit; tmdy(1)=m0, tmdy(2)=m1, tmdy(3)=m2 in Meeus, but lon0 > 0 for E
             if(evMax.eq.3) then
                tmdy(2) = rev(tmdy(1)*pi2 - h0)/pi2 + dTmdy(2)   ! Rise
                tmdy(3) = rev(tmdy(1)*pi2 + h0)/pi2 + dTmdy(3)   ! Set
             end if
             cycle  ! i
          end if  ! if(tmdy(evi).lt.0.)
          
          if(tmdy(evi).gt.1.d0) then  ! in case after m_i = m_i+1, m_f > 1
             tmdy(evi) = 0.d0
             azAlt(evi) = 0.d0
          end if  ! if(tmdy(evi).gt.1.d0)
       end if  ! if(pl.eq.0)
       
       
       if(tmdy(evi).lt.0. .and. deq0(rsAlt)) then
          tmdy(evi) = 0.d0
          azAlt(evi) = 0.d0
       end if
       
       evi = evi+1  ! evi=1,evMax
    end do  ! do while(evi.le.evMax)
    
    
    ! Store results:
    rt = tmdy(2)*24                    ! Rise time - days -> hours
    tt = tmdy(1)*24                    ! Transit time - days -> hours
    st = tmdy(3)*24                    ! Set time - days -> hours
    
    rh = azAlt(2)                      ! Rise azimuth
    ta = azAlt(1) + refract(azAlt(1))  ! Transit altitude + refraction
    sh = azAlt(3)                      ! Set azimuth
    
  end subroutine riset_ipol
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief Rise, transit and set times routine for an object with fixed ra & dec
  !!
  !! \param jd     Julian day number
  !! \param ra     Right ascension (rad)
  !! \param dec    Declination (rad)
  !! \param rsAlt  Altitude to return rise/set data for (degrees; 0. is actual rise/set).  rsAlt>90: compute transit only
  !!
  !! \param rt   Rise time (hours) (output)
  !! \param tt   Transit time (hours) (output)
  !! \param st   Set time (hours) (output)
  !! \param rh   Rising wind direction (rad) (output)
  !! \param ta   Transit altitude (rad) (output)
  !! \param sh   Setting wind direction (rad) (output)
  !!
  !! \param cWarn  Warn upon convergence failure (optional; default: true)
  !!
  !! for rsAlt = 0., rise and set times are computed
  !! for rsAlt.ne.0, the routine calculates when alt=rsAlt is reached
  !! use rsAlt=-6,-12,-18 for the sun for twilight calculations
  !! (rsAlt is expressed in degrees).
  !!
  !! \see
  !! - Meeus, Astronomical algorithms, Ch.15, but with geographic longitude east of Greenwich defined as > 0
  
  subroutine riset_ad(jd, ra,dec, rt,tt,st,rh,ta,sh, rsAlt, cWarn)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi2, d2r
    use SUFR_angles, only: rev, rev2
    use SUFR_date_and_time, only: cal2jd, jd2cal
    use SUFR_numerics, only: deq0
    
    use TheSky_planets, only: planet_position
    use TheSky_planetdata, only: planpos
    use TheSky_local, only: tz, lat0,lon0
    
    implicit none
    real(double), intent(in) :: jd, ra,dec, rsAlt
    real(double), intent(out) :: rt,tt,st,rh,ta,sh
    logical, intent(in), optional :: cWarn
    
    integer :: evi,iter,yr,mnt,evMax
    real(double) :: dy,day0,jd1,tmdy(3),  rsa,cosH0,h0,agst0,th0,dtm,accur,  ha,alt,azAlt(3) !,n
    character :: event(3)*(13)
    logical :: lcWarn
    
    lcWarn = .true.
    if(present(cWarn)) lcWarn = cWarn
    
    alt = 0.d0;  ha = 0.d0;  h0 = 0.d0
    
    accur = 1.d-6  ! Accuracy.  1d-6 ~ 0.1s. Don't make this smaller than 1d-16
    
    event = ['Transit time ','Rise time    ','Set time     ']
    
    rsa = rsAlt*d2r
    evMax = 3                     ! 'Maximum' event to compute: transit only: evMax=1, +rise/set: evMax=3
    if(rsAlt.gt.90.d0) evMax = 1  ! Compute transit time and altitude only
    
    call jd2cal(jd,yr,mnt,dy)
    day0 = dble(int(dy))-tz/24.d0  ! Midnight local time, needed for agst0
    jd1  = cal2jd(yr,mnt,day0)
    call planet_position(jd1,3)
    agst0 = planpos(45)            ! AGST for midnight
    
    tmdy = 0.d0
    azAlt = 0.d0
    if(evMax.eq.3) then
       cosH0 = (sin(rsa)-sin(lat0)*sin(dec))/(cos(lat0)*cos(dec))  ! Cosine of the hour angle of rise/set; Meeus, Eq.15.1
       if(abs(cosH0).gt.1.d0) then  ! Body never rises/sets
          evMax = 1                 ! Compute transit time and altitude only
       else
          h0 = rev(2*acos(cosH0))/2.d0  ! Hour angle of rise/set
       end if
    end if
    
    tmdy(1) = rev(ra - lon0 - agst0)/pi2    ! Transit;  tmdy(1)=m0, tmdy(2)=m1, tmdy(3)=m2 in Meeus, but lon0 > 0 for E
    if(evMax.eq.3) then
       tmdy(2) = rev(tmdy(1)*pi2 - h0)/pi2  ! Rise
       tmdy(3) = rev(tmdy(1)*pi2 + h0)/pi2  ! Set
    end if
    
    
    do evi=1,evMax
       iter = 0
       dtm = huge(dtm)
       do while(abs(dtm).ge.accur)
          th0 = agst0 + 6.300388092591991d0*tmdy(evi)  ! Meeus, p.103
          
          ha = rev2(th0 + lon0 - ra)                                   ! Hour angle;  Meeus p.103
          alt = asin(sin(lat0)*sin(dec) + cos(lat0)*cos(dec)*cos(ha))  ! Altitude;  Meeus Eq. 13.6
          
          ! Correction on transit/rise/set times:
          if(evi.eq.1) then  ! Transit
             dtm = -ha/pi2
          else              ! Rise/set
             dtm = (alt-rsa)/(pi2*cos(dec)*cos(lat0)*sin(ha))
          end if
          tmdy(evi) = tmdy(evi) + dtm
          
          iter = iter+1
          if(iter.gt.30) exit  ! do-while loop
       end do  ! do while(abs(dtm).ge.accur)
       
       
       if(iter.gt.30) then  ! Convergence failed
          if(lcWarn) write(0,'(A,F10.3,A)') '  * WARNING:  riset_ad():  Riset failed to converge: '//trim(event(evi)),rsAlt,'d *'
          
          tmdy(evi) = 0.d0
          azAlt(evi) = 0.d0
       else               ! Result converged, store it
          if(evi.eq.1) then
             azAlt(evi) = alt  ! Transit altitude
          else
             azAlt(evi) = atan2(sin(ha),(cos(ha)*sin(lat0)-tan(dec)*cos(lat0)))  ! Rise,set hour angle -> azimuth
          end if
       end if
       
       if(tmdy(evi).lt.0. .and. deq0(rsAlt)) then
          tmdy(evi) = 0.d0
          azAlt(evi) = 0.d0
       end if
       
    end do  ! evi
    
    
    ! Store results:
    rt = tmdy(2)*24  ! Rise time - days -> hours
    tt = tmdy(1)*24  ! Transit time - days -> hours
    st = tmdy(3)*24  ! Set time - days -> hours
    
    rh = azAlt(2)  ! Rise azimuth
    ta = azAlt(1)  ! Transit altitude
    sh = azAlt(3)  ! Set azimuth
    
  end subroutine riset_ad
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief Compute the date (m,d) at which an object with ra can be observed best, i.e., it transits at midnight  
  !!
  !! \param  jd00  Julian day?
  !! \param  ra0   Right ascension
  !!
  !! \param m     Month of best observation (output)
  !! \param d     Day of month of best observation (output)
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
  !> \brief  Interpolate rise/set times for riset_ipol()
  
  function rsIpol(y1,y2,y3, n)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: y1,y2,y3,n
    real(double) :: rsIpol, a,b,c
    
    a = y2-y1
    b = y3-y2
    c = b-a
    rsIpol = y2 + n/2.d0 * (a + b + c*n)
    
  end function rsIpol
  !*********************************************************************************************************************************
  
  
end module TheSky_riset
!***********************************************************************************************************************************
