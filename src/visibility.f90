!> \file visibility.f90  Contains procedures to determine the visibility of objects for libTheSky


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


!  best_planet_visibility:              Find the best moment (JD) to observe a planet on a given day (JD)
!  planet_visibility_tonight:           Compute when a given planet is visible in a given night
!  comet_invisible:                     Determine whether a comet is invisible, given a magnitude and altitude limit
!  transitalt:                          Compute the transit altitude for a given geographic latitude and declination
!  get_dra_obj:                         Compute the difference between a given right ascension and the RA of the Sun
!  best_obs_date_ra:                    Compute the best date to observe an object with a given right ascension



!***********************************************************************************************************************************
!> \brief  Procedures to determine the visibility of objects

module TheSky_visibility
  implicit none
  save

contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Find a rough indication for the best moment (JD) to observe a planet on a given day (jd).  Start out by demanding 
  !!           that Sun alt < -6deg, but relax this if problematic.  This often results in the moment where Sun alt is indeed -6deg.
  !!
  !! \param  jd                      Julian day to compute best moment for
  !! \param  pl                      Planet to compute visibility for (0-Moon, 1-Mercury, etc.)
  !! \retval best_planet_visibility  Best moment to observe the planet (JD)
  
  function best_planet_visibility(jd,pl) 
    use SUFR_kinds, only: double
    use SUFR_constants, only: r2d
    use SUFR_date_and_time, only: jd2cal,cal2jd
    use SUFR_angles, only: rv12
    
    use TheSky_local, only: tz
    use TheSky_riset, only: riset
    use TheSky_planets, only: planet_position
    use TheSky_planetdata, only: planpos, nplanpos
    
    implicit none
    real(double), intent(in) :: jd
    integer, intent(in) :: pl
    integer :: yr,mn
    real(double) :: best_planet_visibility,dy,jd1, sun_maxalt,smaxalt
    real(double) :: srt,stt,sst,srh,sta,ssh, mrt,mtt,mst,mrh,mta,msh, sundat(nplanpos)
    
    ! Maximum altitude of the Sun (starting value, may be relaxed):
    sun_maxalt = -6.d0  ! ~ -6, -9, -12, the first gives nicer 'twilight' maps?
    
    call riset(jd,pl, mrt,mtt,mst, mrh,mta,msh, 0.d0)  ! Planet
    
    ! First, try moment of planet transit:
    call jd2cal(jd,yr,mn,dy)
    dy = floor(dy) + (mtt-tz)/24.d0
    jd1 = cal2jd(yr,mn,dy)
    
    call planet_position(jd1,3)  ! Sun
    sundat = planpos
    
    if(sundat(30)*r2d.lt.sun_maxalt) then    ! It is "night" at planet transit
       call planet_position(jd1,pl)          ! Planet
       best_planet_visibility = jd1
       return
    end if
    
    
    ! Then try relaxed conditions for the Sun, and demand that planet is at least half as far above the horizon as the Sun is below
    ! (can presumably still see Venus with Sun at -4d and Venus at +2 ...?)
    smaxalt = 2*sun_maxalt
    
    do while(smaxalt.lt.-0.25d0)
       smaxalt = smaxalt/2.d0                               ! Relax conditions
       call riset(jd,3, srt,stt,sst, srh,sta,ssh, smaxalt)  ! Sun
       
       if(abs(rv12(mtt-sst)).lt.abs(rv12(mtt-srt))) then    ! Evening
          dy = floor(dy) + (sst-tz)/24.d0
       else                                                 ! Morning
          dy = floor(dy) + (srt-tz)/24.d0
       end if
       
       jd1 = cal2jd(yr,mn,dy)
       call planet_position(jd1,pl)                         ! Planet
       if(planpos(30)*r2d.gt.-smaxalt/2.d0) exit
    end do
    
    best_planet_visibility = jd1
    
  end function best_planet_visibility
  !*********************************************************************************************************************************
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute when a given planet is visible in a given night
  !!
  !! \param  jd      Julian day for calculation
  !! \param  pl      Planet ID (1-2, 4-8 for Mer-Ven, Mar-Nep)
  !! \param  sunalt  Sun altitude below which planet is visible (deg; e.g. -6.d0)
  !! \param  plalt   Planet altitude above which planet is visible (deg; e.g. 5.d0)
  !! \param  comp    Compute: 1-compute twilight/planet @plalt events only (1,4),  2-include actual rise/set times (2,5)
  !!                 11, 12 = 1, 2 + compute only for today, not tomorrow
  !!
  !! \retval rts     "Rise times":       1-twilight (Sun @sunalt), 2-Sun rise/set, 4-planet @plalt, 5-planet rise/set
  !! \retval tts     Transit times:      1-2, of Sun,  4-5 of planet
  !! \retval sts     "Set times":        1-twilight (Sun @sunalt), 2-Sun rise/set, 4-planet @plalt, 5-planet rise/set
  !! \retval tas     Transit altitudes:  1-2, of Sun,  4-5 of planet
  !! \retval plvis   Planet visibility times (hours):  1-begin, 2-end;  plvis(1)*plvis(2) = 0 when invisible
  !!
  
  subroutine planet_visibility_tonight(jd, pl, sunalt, plalt, comp,  rts, tts, sts, tas, plvis)
    use SUFR_kinds, only: double
    use SUFR_constants, only: r2d
    use SUFR_angles, only: rv12
    use SUFR_sorting, only: sorted_index_list
    
    use TheSky_riset, only: riset
    
    implicit none
    real(double), intent(in) :: jd, sunalt, plalt
    integer, intent(in) :: pl, comp
    real(double), intent(out) :: rts(5), tts(5), sts(5), tas(5), plvis(2)
    
    integer :: irs, ix, maxi, ind(4)
    real(double) :: alt, times(4), st,rh,sh
    logical :: sup,pup,pvis,pvisold
    
    rts=0.d0; tts=0.d0; sts=0.d0; tas=0.d0; plvis=0.d0
    
    ! On day JD, we want set times for JD and rise times for JD+1:
    
    ! Sun:
    maxi = 1                              ! exclude rise/set times
    if(comp.eq.2.or.comp.eq.12) maxi = 2  ! include rise/set times
    alt = sunalt
    do irs=1,maxi   ! 1-twilight (sun @sunalt), 2-Sun rise/set
       if(irs.eq.2) alt = 0.d0
       
       ! Today:    s.set sts:
       call riset(jd, 3,  rts(irs), tts(irs), sts(irs), rh, tas(irs), sh,  alt)
       
       ! Tomorrow: s.rise rts(2) / st. twl rts(1):
       if(comp.lt.10) call riset(jd+1.d0, 3,  rts(irs), tts(irs), st, rh, tas(irs), sh,  alt)
    end do  ! irs=1,2
    
    
    maxi = 4                              ! exclude rise/set times
    if(comp.eq.2.or.comp.eq.12) maxi = 5  ! include rise/set times
    alt = plalt
    do irs=4,maxi   ! 4-planet @plalt, 5-planet rise/set
       if(irs.eq.5) alt = 0.d0
       
       ! Today: planet set:
       call riset(jd, pl,  rts(irs), tts(irs), sts(irs), rh, tas(irs), sh,  alt)
       
       ! Tomorrow: planet rise, transit:
       if(comp.lt.10) call riset(jd+1.d0, pl,  rts(irs), tts(irs), st, rh, tas(irs), sh,  alt)
    end do  ! irs=4,5
    
    
    
    ! When is the planet visible?
    !   plvis(1): time planet becomes visible, plvis(2): time planet becomes invisible
    
    times = (/rv12(rts(1)), rv12(sts(1)),  rv12(rts(4)), rv12(sts(4))/)  ! rv12: -12 - 12h; noon - noon
    call sorted_index_list(times,ind)
    
    
    ! Starting at noon, Sun is up, planet may be:
    sup = .true.   ! Sun up at noon
    pup = .false.  ! Planet not up at noon?
    
    if(pl.gt.9.and.abs(rts(4)*sts(4)).lt.1.d-20) then  ! Comet or asteroid doesn't "rise" or "set" (cross plalt)
       if(tas(4)*r2d.gt.plalt) then  ! Transit alt > plalt: object is always "up"
          plvis(1) = rv12(sts(1))    ! Planet becomes visible at evening twilight
          plvis(2) = rv12(rts(1))    ! Planet becomes invisible at morning twilight
       else
          plvis = 0.d0               ! Planet is never visible
       end if
       return
    end if
    
    if(rv12(rts(4)).ge.0.d0 .and. rv12(rts(4)).gt.rv12(sts(4))) pup = .true.  ! Planet up at noon
    
    pvis = pup.and..not.sup  ! planet visible - should be false
    pvisold = pvis           ! idem
    
    do ix=1,4
       select case(ind(ix))
       case(1)         ! Sunrise/dawn
          sup = .true.
       case(2)         ! Sunset/dusk
          sup = .false.
       case(3)         ! Planet rise
          pup = .true.
       case(4)         ! Planet set
          pup = .false.
       end select
       pvis = pup.and..not.sup  ! Planet visible?
       
       if(pvis.and..not.pvisold) plvis(1) = times(ind(ix))  ! Planet becomes visible
       if(.not.pvis.and.pvisold) plvis(2) = times(ind(ix))  ! Planet becomes invisible
       
       pvisold = pvis
    end do
    
  end subroutine planet_visibility_tonight
  !*********************************************************************************************************************************










  !*********************************************************************************************************************************
  !> \brief  Cheap function to determine whether a comet is invisible at a given time, based on a magnitude and altitude limit
  !!
  !! \param jd        Julian day for calculation
  !! \param cometID   Comet ID
  !! \param mlim      Maximum magnitude, below which the comet is defined as visible (deg)
  !! \param minalt    Minimum transit altitude, above which the comet is defined as visible (deg)
  !!
  !! \note  The comet may still be invisible if comet_invisible=.false., e.g. when it's too close to the Sun
  
  function comet_invisible(jd, cometID, mlim, minalt)
    use SUFR_kinds, only: double
    use SUFR_constants, only: r2d
    use TheSky_planetdata, only: planpos
    use TheSky_comets, only: cometgc
    use TheSky_cometdata, only: cometElems
    use TheSky_local, only: lat0
    use TheSky_coordinates, only: ecl_2_eq
    
    implicit none
    real(double), intent(in) :: jd, mlim, minalt
    integer, intent(in) :: cometID
    real(double) :: tjm, hcr,gcl,gcb,delta, magn, eps, ra,dec, maxalt
    logical :: comet_invisible
    
    comet_invisible = .false.
    
    tjm = (jd-2451545.d0)/365250.d0                    ! Julian Millennia after 2000.0 (not dyn. time - approx)
    call cometgc(tjm,tjm, cometID, hcr,gcl,gcb,delta)
    
    ! Typical difference with full method: ~10^-4 mag:
    magn = cometElems(cometID,8) + 5*log10(delta) + 2.5d0*cometElems(cometID,9)*log10(hcr)
    if(magn.gt.mlim) then
       comet_invisible = .true.
       return
    end if
    
    eps = planpos(48)
    call ecl_2_eq(gcl,gcb,eps, ra,dec)    ! RA, Dec
    
    maxalt = transitalt(lat0, dec) * r2d
    if(maxalt.lt.minalt) comet_invisible = .true.
    
  end function comet_invisible
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the transit altitude of an object with given declination for an observer with a given geographic latitude
  !!
  !! \param  lat  Geographic latitude of the observer (-pi/2 - pi/2; rad)
  !! \param  dec  Declination of the object (-pi/2 - pi/2; rad)
  
  function transitalt(lat, dec)
    use SUFR_kinds, only: double
    implicit none
    real(double), intent(in) :: lat, dec
    real(double) :: transitalt
    
    transitalt = asin(sin(lat)*sin(dec) + cos(lat)*cos(dec))
    
  end function transitalt
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the best date in the year to observe an object with a given right ascension
  !!
  !! \param  year      Year (CE)
  !! \param  RA        Right ascension of the object
  !! \param  accuracy  Get an accurate (but more expensive) result: 0-no, 1-yes
  !! \retval mon       Month of year
  !! \retval dy        Day of month
  
  subroutine best_obs_date_ra(year, RA, accuracy,  mon, dy)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi
    use SUFR_date_and_time, only: doy2md, cal2jd, jd2cal
    use SUFR_solvers, only: root_solver
    
    implicit none
    integer, intent(in) :: year,accuracy
    real(double), intent(in) :: RA
    integer, intent(out) :: mon,dy
    
    integer :: yr,best_obs_doy0,doy1
    real(double) :: jd0,lRA,dRA,djd,jd1,dy1
    
    common /best_obs_ra/ dRA,lRA
    
    lRA = RA
    dRA = pi     ! Find rev2(lRA-sunRA-dRA) = 0 -> dRA=pi = 'opposition'
    yr = year
    
    best_obs_doy0 = 266                                                  ! Approximate doy of opposition of RA=0
    doy1 = mod(nint(lRA/0.017202791805d0) + best_obs_doy0 + 10*366,366)  ! Cheap day of year, +- 5-7? days
    call doy2md(doy1, yr,mon,dy)
    
    if(accuracy.ge.1) then
       jd0 = cal2jd(yr,mon,dble(dy))
       djd = 10.d0
       
       jd1 = root_solver(get_dRA_obj, jd0-djd,jd0+djd, 1.d-1)
       call jd2cal(jd1,yr,mon,dy1)
       dy = nint(dy1)
    end if
    
  end subroutine best_obs_date_ra
  !*********************************************************************************************************************************
  
  !*********************************************************************************************************************************
  !> \brief  Compute the difference between a given right ascension and the RA of the Sun,  used privately by best_obs_date_ra()
  !!
  !! \param jd  Julian day
  !! - Get (RA_obj - RA_sun) for object with RA ra, passed via the common block best_obs_RA
  
  function get_dRA_obj(jd)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev2
    
    use TheSky_planets, only: planet_position
    use TheSky_planetdata, only: planpos
    
    implicit none
    real(double), intent(in) :: jd
    real(double) :: get_dRA_obj,dRA,RA,sunRA
    common /best_obs_RA/ dRA,RA
    
    call planet_position(jd,3)
    sunRA = planpos(5)
    get_dRA_obj = rev2(RA - sunRA - dRA)
    
  end function get_dRA_obj
  !*********************************************************************************************************************************
  
  
  
end module TheSky_visibility
!***********************************************************************************************************************************
