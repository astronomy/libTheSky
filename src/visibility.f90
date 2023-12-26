!> \file visibility.f90  Contains procedures to determine the visibility of objects for libTheSky


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


!  best_planet_xsmag:             Find the moment (JD) of best excess magnitude (mag-lim.mag) for a planet
!  best_planet_visibility:        Find the best moment (JD) to observe a planet on a given day (JD)
!  planet_visibility_tonight:     Compute when a given planet is visible in a given night
!  comet_invisible:               Determine whether a comet is invisible, given a magnitude and altitude limit
!
!  transitalt:                    Compute the transit altitude for a given geographic latitude and declination
!  best_obs_date_ra:              Compute the best date to observe an object with a given right ascension
!  get_dra_obj:                   Compute the difference between a given right ascension and the RA of the Sun
!
!  airmass:                       Compute the airmass for a celestial object with a given TRUE altitude
!  airmass_la:                    Compute the airmass for a celestial object with given APP alt; low-accuracy alt. for airmass()
!  extinction_magPam:             Compute the extinction in magnitudes per unit airmass for an observer with given elevation
!  extinction_mag:                Compute the extinction in magnitudes for a given object altitude and observer elevation
!  extinction_fac:                Compute the extinction factor for an observer with given elevation and object with given altitude
!
!  limmag_full:                   Calculate limiting magnitude, full function
!  limmag_extinction:             Calculate limiting magnitude based on Sun altitude and Moon phase
!  limmag_skyBrightness:          Calculate sky brightness based on Sun altitude and Moon phase
!  limmag_jd:                     Calculate limiting magnitude based on JD and object altitude, wrapper for limmag_full()
!  limmag_zenith_jd:              Calculate limiting magnitude for the local zenith, based on JD and observer's location
!  limmag_jd_pl:                  Calculate limiting magnitude based on JD and planet ID, wrapper for limmag_jd()
!  limmag_sun_airmass:            Calculate limiting magnitude based on Sun altitude and object altitude (airmass); assume New Moon
!  limmag_sun:                    Calculate limiting magnitude, based on the altitude of the Sun only
!
!  pl_xsmag:                      Compute the excess magnitude (mag-lim.mag) for planet pl at JD, considering Sun, Moon and airmass
!  pl_xsmag_pl:                   Compute the excess magnitude at JD, wrapper for pl_xsmag() for solvers
!  pl_xsmag_la:                   Compute the excess magnitude, considering airmass and Sun altitude only
!  pl_xsmag_la_pl:                Compute the excess magnitude, wrapper for pl_xsmag_la() for solvers
!
!  aperture:                      Aperture needed to observe an object with given excess magnitude
!  skybrightness_mas2_from_mlim:  Convert naked-eye limiting magnitude to sky surface brightness in magnitudes per square arcsecond
!  skybrightness_mlim_from_mas2:  Convert sky surface brightness in magnitudes per square arcsecond to naked-eye limiting magnitude
!  skybrightness_mlim_from_cdm2:  Convert sky brightness in candela per square meter to naked-eye visual limiting magnitude
!  skybrightness_cdm2_from_mlim:  Convert limiting visual magnitude to sky brightness in candela per square meter
!  skybrightness_cdm2_from_nL:    Convert sky brightness in nanolambert to candela per square meter
!  skybrightness_nL_from_cdm2:    Convert sky brightness in candela per square meter to nanolambert

!***********************************************************************************************************************************
!> \brief  Procedures to determine the visibility of objects

module TheSky_visibility
  implicit none
  save

contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Find the moment (JD) of optimal excess magnitude for a planet, i.e. the lowest mag - lim.mag
  !!
  !! \param  jdin     Julian day to compute best moment for
  !! \param  plID     Planet to compute visibility for (0-Moon, 1-Mercury, etc.)
  !!
  !! \param jdout    Moment of best excess magnitude (output)
  !! \param xsmag    Excess magnitude at that moment (magnitude - limiting magnitude; <0: ~visible to the naked eye) (output)
  !!
  !! \param  rsCWarn  Print convergence warnings in riset()  (optional; default = true)
  
  subroutine best_planet_xsmag(jdin,plID, jdout, xsmag, rsCWarn)
    use SUFR_kinds, only: double
    use SUFR_date_and_time, only: cal2jd
    use SUFR_solvers, only: minimum_solver
    
    use TheSky_planetdata, only: pl0
    use TheSky_datetime, only: gettz
    
    implicit none
    real(double), intent(in) :: jdin
    integer, intent(in) :: plID
    real(double), intent(out) :: jdout
    real(double), intent(out), optional :: xsmag
    logical, intent(in), optional :: rsCWarn
    
    integer :: status
    real(double) :: tz, jd1,jd2, plvis(2), lxsmag
    logical :: lrsCWarn
    
    lrsCWarn = .true.  ! riset() convergence warnings on by default
    if(present(rsCWarn)) lrsCWarn = rsCWarn
    
    call planet_visibility_tonight(jdin, plID, 0.d0, 0.d0, 11,  plvis, rsCWarn=lrsCWarn)  ! Sun<0d; Planet>0d, comp=11 twl/today
    
    tz = gettz(jdin)
    jd1 = floor(jdin+0.5d0) - 0.5d0 + (plvis(1)-tz)/24.d0
    jd2 = floor(jdin+0.5d0) - 0.5d0 + (plvis(2)-tz)/24.d0
    if(plvis(1).gt.12.d0) jd1 = jd1 - 1.d0  ! After noon, use the previous day
    if(plvis(2).gt.12.d0) jd2 = jd2 - 1.d0
    
    pl0 = plID  ! Needed by pl_xsmag_pl()
    lxsmag = minimum_solver(pl_xsmag_pl, jd1, (jd1+jd2)/2.d0, jd2, 1.d-3,  jdout, status=status)  ! Use full routine
    
    if(present(xsmag)) xsmag = lxsmag
    
  end subroutine best_planet_xsmag
  !*********************************************************************************************************************************
   
   
   
  !*********************************************************************************************************************************
  !> \brief  Find a rough indication for the best moment (JD) to observe a planet on a given day (jd).  Start out by demanding 
  !!           that Sun alt < -6deg, but relax this if problematic.  This often results in the moment where Sun alt is indeed -6deg.
  !!
  !! \param  jd                      Julian day to compute best moment for
  !! \param  pl                      Planet to compute visibility for (0-Moon, 1-Mercury, etc.)
  !!
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
  !! \param  sunAlt  Sun altitude below which planet is visible (deg; e.g. -6.d0)
  !! \param  plalt   Planet altitude above which planet is visible (deg; e.g. 5.d0)
  !! \param  comp    Compute: 1-compute twilight/planet at plalt events only (1,4),  2-include actual rise/set times (2,5)
  !!                 11, 12 = 1, 2 + compute only for today, not tomorrow
  !!
  !! \param plvis   Planet visibility times (hours):    1-begin, 2-end;  plvis(1)*plvis(2) = 0 when invisible (output)
  !! \param plazs   Planet visibility azimuths (rad):   1-begin of visibility,  2-end of visibility (output)
  !!
  !! \param rts     "Rise times":       1-twilight (Sun at sunAlt), 2-Sun rise/set, 4-planet at plalt, 5-planet rise/set (output)
  !! \param tts     Transit times:      1-2, of Sun,  4-5 of planet (output)
  !! \param sts     "Set times":        1-twilight (Sun at sunAlt), 2-Sun rise/set, 4-planet at plalt, 5-planet rise/set (output)
  !!
  !! \param ras     Rise azimuths:      1-2, of Sun,  4-5 of planet (output)
  !! \param tas     Transit altitudes:  1-2, of Sun,  4-5 of planet (output)
  !! \param sas     Set azimuths:       1-2, of Sun,  4-5 of planet (output)
  !!
  !! \param  rsCWarn  Print convergence warnings in riset()  (optional; default = true)
  
  subroutine planet_visibility_tonight(jd, pl, sunAlt, plalt, comp,   plvis, plazs,   rts, tts, sts,  ras, tas, sas, rsCWarn)
    use SUFR_kinds, only: double
    use SUFR_constants, only: r2d
    use SUFR_angles, only: rv12
    use SUFR_sorting, only: sorted_index_list
    
    use TheSky_riset, only: riset
    
    implicit none
    real(double), intent(in) :: jd, sunAlt, plalt
    integer, intent(in) :: pl, comp
    real(double), intent(out) :: plvis(2)
    real(double), intent(out), optional :: plazs(2),  rts(5), tts(5), sts(5),  ras(5), tas(5), sas(5)
    logical, intent(in), optional :: rsCWarn
    
    integer :: irs, ix, maxi, ind(4)
    real(double) :: alt, times(4),azs(4), st,sh,  lrts(5), ltts(5), lsts(5), lras(5),ltas(5),lsas(5)
    logical :: sup,pup,pvis,pvisold, lrsCWarn
    
    lrts=0.d0; ltts=0.d0; lsts=0.d0; lras=0.d0;ltas=0.d0;lsas=0.d0; plvis=0.d0
    if(present(plazs)) plazs=0.d0
    
    lrsCWarn = .true.  ! riset() convergence warnings on by default
    if(present(rsCWarn)) lrsCWarn = rsCWarn
    
    ! On day JD, we want set times for JD and rise times for JD+1:
    
    ! Sun:
    maxi = 1                              ! exclude rise/set times
    if(comp.eq.2.or.comp.eq.12) maxi = 2  ! include rise/set times
    alt = sunAlt
    do irs=1,maxi   ! 1-twilight (sun @sunAlt), 2-Sun rise/set
       if(irs.eq.2) alt = 0.d0
       
       ! Today:    s.set lsts:
       call riset(jd, 3,   lrts(irs), ltts(irs), lsts(irs),   lras(irs), ltas(irs), lsas(irs),   alt,  cWarn=lrsCWarn)
       
       ! Tomorrow: s.rise lrts(2) / st. twl lrts(1):
       if(comp.lt.10) call riset(jd+1.d0, 3,   lrts(irs), ltts(irs), st,   lras(irs), ltas(irs), sh,   alt,  cWarn=lrsCWarn)
    end do  ! irs=1,2
    
    
    maxi = 4                              ! exclude rise/set times
    if(comp.eq.2.or.comp.eq.12) maxi = 5  ! include rise/set times
    alt = plalt
    do irs=4,maxi   ! 4-planet @plalt, 5-planet rise/set
       if(irs.eq.5) alt = 0.d0
       
       ! Today: planet set:
       call riset(jd, pl,   lrts(irs), ltts(irs), lsts(irs),   lras(irs), ltas(irs), lsas(irs),   alt,  cWarn=lrsCWarn)
       
       ! Tomorrow: planet rise, transit:
       if(comp.lt.10) call riset(jd+1.d0, pl,   lrts(irs), ltts(irs), st,   lras(irs), ltas(irs), sh,   alt,  cWarn=lrsCWarn)
    end do  ! irs=4,5
    
    
    
    ! When is the planet visible?
    !   plvis(1): time planet becomes visible, plvis(2): time planet becomes invisible
    
    times = (/rv12(lrts(1)), rv12(lsts(1)),  rv12(lrts(4)), rv12(lsts(4))/)  ! rv12: -12 - 12h; noon - noon: Sun rise,set, pl r,s
    azs   = (/lras(1), lsas(1),  lras(4), lsas(4)/)
    call sorted_index_list(times,ind)
    
    ! Special case: ind = [2 4 3 1] - Ss Ps Pr Sr: Planet is visible twice.  Program takes latter, make sure it's the longest:
    if(ind(1).eq.2.and.ind(2).eq.4.and.ind(3).eq.3.and.ind(4).eq.1 .and. times(4)-times(2).gt.times(1)-times(3)) then
       ind = [3,1,2,4]
       
       ! If the planet is invisible for less than an hour, ignore it:
       if(times(2)-times(1) .lt. 1.d0) times = [times(2),times(2), times(2),times(1)]  ! Only the last two array elements count
    end if
    
    
    
    ! Starting at noon, Sun is up, planet may be:
    sup = .true.   ! Sun up at noon
    pup = .false.  ! Planet not up at noon?
    
    if(pl.gt.9.and.abs(lrts(4)*lsts(4)).lt.1.d-20) then  ! Comet or asteroid doesn't "rise" or "set" (cross plalt)
       if(ltas(4)*r2d.gt.plalt) then  ! Transit alt > plalt: object is always "up"
          plvis(1) = rv12(lsts(1))    ! Planet becomes visible at evening twilight
          plvis(2) = rv12(lrts(1))    ! Planet becomes invisible at morning twilight
          
          if(present(plazs)) then
             plazs(1) = lsas(1)       ! Azimuth where planet becomes visible at evening twilight
             plazs(2) = lras(1)       ! Azimuth where planet becomes invisible at morning twilight
          end if
       else
          plvis = 0.d0               ! Planet is never visible
          if(present(plazs)) plazs = 0.d0
       end if
       return
    end if
    
    if(rv12(lrts(4)).ge.0.d0 .and. rv12(lrts(4)).gt.rv12(lsts(4))) pup = .true.  ! Planet up at noon
    
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
       
       if(pvis.and..not.pvisold) then  ! Planet becomes visible
          plvis(1) = times(ind(ix))
          if(present(plazs)) plazs(1) = azs(ind(ix))
       end if
       if(.not.pvis.and.pvisold) then  ! Planet becomes invisible
          plvis(2) = times(ind(ix))
          if(present(plazs)) plazs(2) = azs(ind(ix))
       end if
       
       pvisold = pvis
    end do
    
    ! Assign values to optional dummy arguments:
    if(present(rts)) rts = lrts
    if(present(tts)) tts = ltts
    if(present(sts)) sts = lsts
    
    if(present(ras)) ras = lsas  ! Rise azimuths
    if(present(tas)) tas = ltas  ! Tranist altitudes
    if(present(sas)) sas = lras  ! Set azimuths
    
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
    use SUFR_constants, only: r2d, jd2000
    use TheSky_planetdata, only: planpos
    use TheSky_comets, only: cometgc
    use TheSky_cometdata, only: cometElems, cometDiedAtP
    use TheSky_local, only: lat0
    use TheSky_coordinates, only: ecl_2_eq
    
    implicit none
    real(double), intent(in) :: jd, mlim, minalt
    integer, intent(in) :: cometID
    real(double) :: tjm, hcr,gcl,gcb,delta, magn, eps, ra,dec, maxalt
    logical :: comet_invisible
    
    comet_invisible = .false.
    
    tjm = (jd-jd2000)/365250.d0   ! Julian Millennia after 2000.0 (not dyn. time - approx)
    call cometgc(tjm,tjm, cometID, hcr,gcl,gcb,delta)
    
    ! Typical difference with full method: ~10^-4 mag:
    if(cometDiedAtP(cometID).ne.0 .and. jd.gt.cometElems(cometID,7)) then
       magn = 99.9d0                                                                           ! Comet died at perihelion
    else
       magn = cometElems(cometID,8) + 5*log10(delta) + 2.5d0*cometElems(cometID,9)*log10(hcr)  ! m = H + 5log(d) + 2.5*G*log(r)
    end if
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
  !! \param mon       Month of year (output)
  !! \param dy        Day of month (output)
  
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
  !!
  !! \note
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
    
    call planet_position(jd,3)            ! Compute Sun position
    sunRA = planpos(5)                    ! Sun's right ascension
    get_dRA_obj = rev2(RA - sunRA - dRA)
    
  end function get_dRA_obj
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the airmass for a celestial object with a given TRUE altitude
  !!
  !! \param alt  TRUE altitude of object (radians)
  !!
  !! \note
  !! - Results are 1 <= airmass <~ 38.35; return ~1000 for h<-34'
  !! - Maximum supposed error (at the horizon) of 0.0037 air mass
  !!
  !! \see Young (1994); https://en.wikipedia.org/wiki/Air_mass_%28astronomy%29#Interpolative_formulas
  
  function airmass(alt)
    use SUFR_kinds, only: double
    use SUFR_constants, only: am2r
    
    implicit none
    real(double), intent(in) :: alt
    real(double) :: airmass,sinh,sinh2
    
    if(alt.lt.-34*am2r) then  ! h<-34' - true altitude can be <0
       airmass = 1000.d0 * (0.15d0 + abs(alt))  ! Very bad (adds at least ~30 magnitudes due to extinction), 
       !                                          but still worse when farther below the horizon - for solvers
    else
       sinh = sin(alt)
       sinh2 = sinh**2
       airmass = (1.002432d0*sinh2 + 0.148386d0*sinh + 0.0096467d0) / &
            (sinh2*sinh + 0.149864d0*sinh2 + 0.0102963d0*sinh + 0.000303978d0)
       airmass = max( airmass, 1.d0 )
    end if
    
  end function airmass
  !*********************************************************************************************************************************
  
   
  !*********************************************************************************************************************************
  !> \brief  Compute the airmass for a celestial object with a given apparent altitude; low-accuracy alternative for airmass()
  !!
  !! \param alt  Apparent altitude of object (radians)
  !!
  !! \note
  !! - Results are 1 <= airmass <~ 37.92; return ~1000 for h<0
  !!
  !! \see Kasten and Young (1989); http://en.wikipedia.org/wiki/Airmass_%28astronomy%29#Interpolative_formulas
  
  function airmass_la(alt)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pio2, r2d
    
    implicit none
    real(double), intent(in) :: alt
    real(double) :: airmass_la,z,zdeg
    
    if(alt.lt.0.d0) then
       airmass_la = 1000.d0 * (0.15d0 + abs(alt))  ! Very bad (adds at least ~30 magnitudes due to extinction), 
       !                                             but still worse when farther below the horizon - for solvers
    else
       z = min(pio2 - alt, pio2)  ! Zenith angle
       zdeg = z*r2d
       airmass_la = max( 1.d0 / ( cos(z) + 0.50572d0*(96.07995d0-zdeg)**(-1.6364d0) ) ,  1.d0 )
       
    end if
    
  end function airmass_la
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the extinction in magnitdes per unit airmass for an observer with given elevation
  !!
  !! \param ele  Evelation of the observer above sea level (metres)
  !!
  !! \note  
  !!   - The magnitude of an object corrected for airmass should be  m' = m + extinction_magPam(ele) * airmass(alt) 
  !!     (see function extinction_mag())
  !!   - This assumes night vision, using the rods in the human eye's retina.  Hence, this will not hold for the Sun!
  !!     (There is some correspondence when only considering ~420-580 nm, which is roughly the sensitivity band of the rods.)
  !!
  !! \see  Green, ICQ 14, 55 (1992),  http://www.icq.eps.harvard.edu/ICQExtinct.html, (for lambda=510 nm!), based on
  !!       Hayes & Latham, ApJ 197, 593 (1975): http://esoads.eso.org/abs/1975ApJ...197..593H (wavelength dependent, for Wega)
  
  function extinction_magPam(ele)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: ele
    real(double) :: extinction_magPam, Aoz,Aray,Aaer
    
    Aoz  = 0.016d0                         ! Ozone  (Schaefer 1992)
    Aray = 0.1451d0 * exp(-ele/7996.d0)    ! Rayleigh scattering (for lambda=510 nm), Eq.2
    Aaer = 0.120d0  * exp(-ele/1500.d0)    ! Aerosol scattering (for lambda = 0.51 micron), Eq.4
    
    extinction_magPam = Aoz + Aray + Aaer  ! Total extinction in magnitudes per unit air mass (~0.2811 at sea level)
    
  end function extinction_magPam
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the extinction in magnitdes for an observer with given elevation and an object with given TRUE altitude
  !!
  !! \param alt  Altitude of object (radians)
  !! \param ele  Evelation of the observer above sea level (metres; optional)
  !!
  !! \note - The magnitude of an object corrected for airmass should be  m' = m + extinction_mag(alt,ele)
  !!       - Assumed is the response of the rods in the human retina, hence night vision
  !!
  !! \see   Functions extinction_magPam() and airmass()
  
  function extinction_mag(alt, ele)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: alt
    real(double), intent(in), optional :: ele
    real(double) :: extinction_mag, elel
    
    elel = 0.d0                  ! Observer is at sea level by default
    if(present(ele)) elel = ele  ! User-specified observer elevation
    
    extinction_mag = extinction_magPam(elel) * airmass(alt)  ! Extinction in magnitudes per unit air mass x airmass
    
  end function extinction_mag
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the extinction factor for an observer with given elevation and an object with given altitude
  !!
  !! \param alt  Altitude of object (radians)
  !! \param ele  Evelation of the observer above sea level (metres; optional)
  !!
  !! \note - extinction_fac = 1: no extinction, extinction_fac > 1 extinction.
  !!       - Hence, the flux, corrected for extinction, should be  f' = f / extinction_fac(alt,ele)
  !!       - Assumed is the response of the rods in the human retina, hence night vision
  !!
  !! \see  function extinction_mag()
  
  function extinction_fac(alt, ele)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: alt
    real(double), intent(in), optional :: ele
    real(double) :: extinction_fac, elel
    
    elel = 0.d0                  ! Observer is at sea level by default
    if(present(ele)) elel = ele  ! User-specified observer elevation
    
    extinction_fac = 10.d0**( extinction_mag(alt,elel) / 2.5d0)  ! Extinction in magnitudes -> factor
    
  end function extinction_fac
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate limiting magnitude based on Sun and Moon positions and phase, and on the altitude of the observed object
  !!
  !! \param  year        Year of observation (for solar cycle)
  !! \param  month       Month of observation (for approximate Sun RA)
  !! \param  obsElev     Elevation of the observer (m)
  !! \param  obsLat      Latitude of the observer (rad)
  !!
  !! \param  sunAlt      Altitude of the Sun (rad)
  !! \param  sunElon     Elongation object-Sun (rad)
  !!
  !! \param  moonPhase   Phase of the Moon (fraction)
  !! \param  moonAlt     Altitude of the Moon (rad)
  !! \param  moonElon    Elongation object-Moon (rad)
  !!
  !! \param  objAlt      Altitude of the observed object (rad)
  !!
  !! \param  humid       Relative humidity  (%; optional, default: 70)
  !! \param  airTemp     Air temperature  (degrees Celsius; optional, default: 10)
  !! \param  snrat       Snellen ratio for vision  (optional, default: 1)
  !!
  !! \param  verbosity   Verbosity  (optional, default: 0=quiet)
  !!
  !! \retval limmag_full  Limiting magnitude - return >= 99 if the object is below the horizon
  !!
  !! \see BASIC program VISLIMIT.BAS by Bradley E. Schaefer, Sky and Telescope, May 1998, p.57,
  !!      in turn based on Schaefer Sky&Tel 78, 522 (1989).
  
  function limmag_full(year,month, obsElev,obsLat, sunAlt,sunElon, moonPhase,moonAlt,moonElon, objAlt, humid,airTemp,snrat, verbosity)
    use SUFR_kinds, only: double
    
    implicit none
    integer, intent(in) :: year, month
    real(double), intent(in) :: obselev,obslat, sunAlt,sunElon, moonPhase,moonAlt,moonElon, objalt
    integer, intent(in), optional :: verbosity
    real(double), intent(in), optional :: humid,airTemp,snrat
    
    integer :: verb
    real(double) :: limmag_full,  skyBr(5),extCoef(5),extMag(5),  relHum,temp,sn, c1,c2,th
    
    if(objAlt.lt.0.d0) then
       limmag_full = 99.d0 - objAlt  ! Object below the horizon; limiting magnitude meaningless
       return
    end if
    
    
    ! Optional arguments:
    relHum = 70.d0                       ! Relative humidity (%), is higher at night than during the day, annual average at 7:00
    if(present(humid)) relHum = humid    ! in S-E USA ~83%, see http://www.sercc.com/climateinfo/historical/avgrh.html
    
    temp = 10.d0                         ! Temperature (deg.C)
    if(present(airTemp)) temp = airTemp
    
    sn = 1.d0                            ! Snellen ratio quantifies vision and observing experience
    if(present(snrat)) sn = snrat
    
    verb = 0
    if(present(verbosity)) verb = verbosity  ! Verbosity
    
    ! Compute the extinction and sky brightness:
    call limmag_extinction(month, obsLat,obsElev, objAlt, relHum,temp, 3,3,  extCoef, extMag)  ! Compute only V band (3)
    call limmag_skyBrightness(year, sunAlt,sunElon, moonPhase,moonAlt,moonElon, objAlt, 3,3,  extCoef,  skyBr)  ! Only V band (3)
    
    
    ! Visual limiting magnitude (V=3 in UBVRI):
    !> \see Knoll, Tousey and Hulburt, J.Opt.Soc.Am. 36, 480 (1946); Hecht J.Opt.Soc.Am. 37, 59 (1947), Eq.3; 
    !!      Schaefer PASP 102, 212 (1990), Eq.2
    if(skyBr(3).lt.1500.d0) then  ! 1500 nanoLamberts; log(B) ~ 3.17
       c1 = 10.d0**(-9.8d0)
       c2 = 10.d0**(-1.9d0)
    else
       c1 = 10.d0**(-8.35d0)
       c2 = 10.d0**(-5.9d0)
    end if
    th = c1*(1.d0 + sqrt(c2*skyBr(3)))**2  ! Visual sky brightness in foot-candles
    
    limmag_full = -16.57d0 - 2.5d0*log10(th) - extMag(3) + 5*log10(sn)  ! Limiting visual magnitude (V)
    
    if(verb>0) then
       write(*,'(A,F12.5)')  'Extinction coefficient (V band):                       ', extCoef(3)
       write(*,'(A,F12.5)')  'Extinction in V magnitudes:                            ', extMag(3)
       write(*,'(A,ES12.5)') 'Sky surface brightness (V) in erg/s/cm^2/μ/arcsec^2:   ', skyBr(3)*1.02d-15  !   Should this be 1.02d-15 iso 1.11d-15?  Or 4.54d-16? - http://www.archaeocosmology.org/eng/vislimitbas.htm
       write(*,'(A,F12.5)')  'Sky surface brightness (V) in nL:                      ', skyBr(3)
       write(*,'(A,ES12.5)') 'Sky brightness (V) in foot candles:                    ', th
       write(*,'(A,F12.5)')  'V magnitude limit:                                     ', limmag_full
    end if
    
  end function limmag_full
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate limiting magnitude based on Sun altitude and Moon phase
  !!
  !! \param  month    Month of observation (for approximate Sun RA)
  !!
  !! \param  obsLat   Latitude of the observed (rad)
  !! \param  obsElev  Elevation of the observer (m)
  !!
  !! \param  objAlt   Altitude of the observed object (rad)
  !!
  !! \param  relHum   Relative humidity
  !! \param  temp     Air temperature
  !!
  !! \param  band1    First of UBVRI bands to compute (1-5); compute band1-band2
  !! \param  band2    Last of UBVRI bands to compute (1-5); compute band1-band2
  !!
  !! \param extCoef  Atmospheric extinction for UBVRI (output)
  !! \param extMag   Delta magnitude due to atmospheric extinction for UBVRI (output)
  !!
  !! \see BASIC program VISLIMIT.BAS by Bradley E. Schaefer, Sky and Telescope, May 1998, p.57,
  !!      in turn based on Schaefer Sky&Tel 78, 522 (1989).
  
  subroutine limmag_extinction(month, obsLat,obsElev, objAlt, relHum,temp, band1,band2,  extCoef,extMag)
    use SUFR_kinds, only: double
    use SUFR_constants, only: d2r, earthr
    
    implicit none
    integer, intent(in) :: month, band1,band2
    real(double), intent(in) :: obsLat, objAlt, obsElev, relHum, temp
    real(double), intent(out) :: extCoef(5),extMag(5)
    
    integer :: iBand, band1l,band2l, sl
    real(double) :: wa(5), oz(5),wt(5),  ra,xg,xa,xo,kr,ka,ko,kw, Re_km, sinObjAlt
    
    ! UBVRI bands are 1-5:
    band1l = max(1,band1)
    band2l = min(5,band2)
    
    extCoef = 0.d0;  extMag = 0.d0
    wa = [ 0.365d0,   0.44d0,   0.55d0,    0.7d0,    0.9d0 ]  ! Gas/aerosols
    oz = [ 0.000d0,  0.000d0,  0.031d0,  0.008d0,  0.000d0 ]  ! Ozone
    wt = [ 0.074d0,  0.045d0,  0.031d0,  0.020d0,  0.015d0 ]
    
    
    ra = (month-3)*30*d2r                       ! Derive a rough Sun RA from the month number
    sl = nint(obsLat/abs(obsLat+tiny(obsLat)))  ! N: +1;  S: -1
    
    ! Airmass for each component:
    Re_km = earthr*1.d-5  ! Earth's radius cm -> km
    sinObjAlt = sin(objAlt)
    xg = 1.d0 / ( sinObjAlt + 0.0286d0 * exp(-10.5d0*sinObjAlt) )       ! Gas
    xa = 1.d0 / ( sinObjAlt + 0.0123d0 * exp(-24.5d0*sinObjAlt) )       ! Aerosol
    xo = 1.d0 / sqrt( 1.d0 - (cos(objAlt)/(1.d0 + (20.d0/Re_km)))**2 )  ! Ozone
    
    ! UBVRI (1-5) extinction for each component:
    do iBand=band1l,band2l
       kr = 0.1066d0 * exp(-obsElev/8200.d0) * (wa(iBand)/0.55d0)**(-4)
       ka = 0.1d0 * (wa(iBand)/0.55d0)**(-1.3d0) * exp(-obsElev/1500.d0)
       ka = ka * (1.d0 - 0.32d0/log(relHum/100.d0))**1.33d0 * (1.d0 + 0.33d0 * sl * sin(ra))  ! Correction from http://www.archaeocosmology.org/eng/vislimitbas.htm
       ko = oz(iBand) * (3.d0 + 0.4d0 * (obsLat*cos(ra) - cos(3*obsLat)) )/3.d0
       kw = wt(iBand) * 0.94d0 * (relHum/100.d0) * exp(temp/15.d0 - obsElev/8200.d0)
       
       extCoef(iBand)  = kr + kw + ka + ko           ! Extinction in magnitudes per unit air mass (?)
       extMag(iBand)   = (kr+kw)*xg + ka*xa + ko*xo  ! Total extinction in magnitudes (?)
    end do  ! iBand
    
  end subroutine limmag_extinction
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate sky brightness based on Sun altitude and Moon phase
  !!
  !! \param year       Year CE (for the solar cycle)
  !!
  !! \param objAlt     Altitude of the observed object (rad)
  !! \param moonAlt    Moon altitude (rad)
  !! \param sunAlt     Sun altitude (rad)
  !!
  !! \param moonPhase  Moon illuminated fraction (0-1)
  !! \param moonElon   Moon elongation from observed object (rad)
  !! \param sunElon    Sun elongation from observed object (rad)
  !!
  !! \param band1      First of UBVRI bands to compute (1-5); compute band1-band2
  !! \param band2      Last of UBVRI bands to compute (1-5); compute band1-band2
  !!
  !! \param extCoef    Extinction for UBVRI
  !!
  !! \param skyBr      Sky brightness for UBVRI in nanoLambert (output)
  !!
  !! \see BASIC program VISLIMIT.BAS by Bradley E. Schaefer, Sky and Telescope, May 1998, p.57,
  !!      in turn based on Schaefer Sky&Tel 78, 522 (1989).
  
  subroutine limmag_skyBrightness(year, sunAlt,sunElon, moonPhase,moonAlt,moonElon, objAlt, band1,band2, extCoef,  skyBr)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi,pi2,pio2, r2d
    use TheSky_moon, only: moonmagn
    
    implicit none
    integer, intent(in) :: year, band1,band2
    real(double), intent(in), value :: objAlt, moonAlt,sunAlt, moonPhase,moonElon, sunElon  ! May be called as planpos()
    real(double), intent(in) :: extCoef(5)
    real(double), intent(out) :: skyBr(5)
    
    integer :: iBand, band1l,band2l
    real(double) :: mo(5),bo(5),cm(5),ms(5),  objAM,moonAM,sunAM,bn, moonPA,mm,c3, moonElonM, fm,bm,bt,c4,fs,bd
    
    ! UBVRI bands are 1-5:
    band1l = max(1, band1)
    band2l = min(5, band2)
    
    
    ! Initialise data:
    skyBr = 0.d0
    mo = [ -10.93d0, -10.45d0, -11.05d0, -11.90d0, -12.70d0 ]  ! Moon magnitudes(?)
    bo = [  8.0d-14,  7.0d-14,  1.0d-13,  1.0d-13,  3.0d-13 ]  ! Dark sky (erg/s/cm^2/μ/arcsec^2)
    cm = [   1.36d0,   0.91d0,   0.00d0,  -0.76d0,  -1.17d0 ]  ! UBVRI colour correction for the Moon
    ms = [ -25.96d0, -26.09d0, -26.74d0, -27.26d0, -27.55d0 ]  ! Sun magnitude
    
    ! Air mass for the observed object:
    objAM  = 1.d0/(sin(objAlt) + 0.025d0*exp(-11.d0*sin(objAlt)))
    
    ! Air mass for the Moon:
    if(moonAlt.lt.0.d0) then  ! Moon below the horizon
       moonAM = 40.d0
    else
       moonAM = 1.d0/(sin(moonAlt) + 0.025d0*exp(-11.d0*sin(moonAlt)))
    end if
    
    ! Air mass for the Sun:
    if(sunAlt.lt.0.d0) then  ! Sun below the horizon
       sunAM = 40.d0
    else
       sunAM = 1.d0/(sin(sunAlt) + 0.025d0*exp(-11.d0*sin(sunAlt)))
    end if
    
    
    ! UBVRI (1-5) sky brightness for each band (erg/s/cm^2/μ/arcsec^2):
    do iBand=band1l,band2l
       
       ! Dark night-sky brightness:
       bn = bo(iBand) * (1.d0 + 0.3d0*cos(pi2*dble(year-1992)/11.d0))  ! Solar cycle
       bn = bn * (0.4d0 + 0.6d0/sqrt(1.d0-0.96d0*cos(objAlt)**2))
       bn = bn * 10.d0**(-0.4d0*extCoef(iBand)*objAM)
       
       ! Moonlight brightness:
       if(moonAlt.ge.0.d0) then                           ! Moon is up
          moonPA = (1.d0 - moonPhase)*pi                           ! Approximate phase angle of the Moon (rad)
          mm = moonmagn(moonPA,2.5696d-3) + cm(iBand)              ! Moon magnitude for mean distance (2.57e-3 AU) + colour correction
          
          c3 = 10.d0**(-0.4d0*extCoef(iBand)*moonAM)
          
          moonElonM = max(moonElon,4.36d-3)                        ! Elongation Moon-object (rad) - can't be <~0.25d ~ 4.36e-3 rad
          fm = 1.888628d4/moonElonM**2 + 10.d0**(6.15d0 - moonElonM/0.6981317d0)
          fm = fm + 10.d0**5.36d0 * (1.06d0 + cos(moonElonM)**2)
          
          bm = 10.d0**( -0.4d0 * (mm - mo(iBand) + 43.27d0) )      ! Brightness of the Moon (erg/s/cm^2/μ/arcsec^2)
          bm = bm * (1.d0 - 10.d0**(-0.4d0*extCoef(iBand)*objAM))
          bm = bm * (fm*c3 + 440000.d0*(1.d0-c3))
       else
          bm = 0.d0  ! No contribution from the Moon if below the horizon
       end if
       
       ! Twilight brightness (erg/s/cm^2/μ/arcsec^2):
       bt = 10.d0**( -0.4d0*(ms(iBand) - mo(iBand) + 32.5d0 - sunAlt*r2d - ((pio2-objAlt)/(pi2*extCoef(iBand))) ) )
       bt = bt * (1.745329252d0/sunElon) * (1.d0 - 10.d0**((-0.4d0*extCoef(iBand))*objAM))    ! Correction from http://www.archaeocosmology.org/eng/vislimitbas.htm
       
       ! Daylight brightness (erg/s/cm^2/μ/arcsec^2):
       c4 = 10.d0**(-0.4d0*extCoef(iBand)*sunAM)
       fs = 1.888628d4/sunElon**2 + 10.d0**(6.15d0 - sunElon/0.6981317d0)
       fs = fs + 10.d0**5.36d0 * (1.06d0 + cos(sunElon)**2)
       bd = 10.d0**(-0.4d0*(ms(iBand) - mo(iBand) + 43.27d0))
       bd = bd * (1.d0 - 10.d0**(-0.4d0*extCoef(iBand)*objAM))
       bd = bd * (fs*c4 + 440000.d0*(1.d0-c4))
       
       ! Total sky brightness (erg/s/cm^2/μ/arcsec^2):
       if(bd.lt.bt) then
          skyBr(iBand) = bn + bd + bm  ! Daylight
       else
          skyBr(iBand) = bn + bt + bm  ! Twilight
       end if
       
       skyBr(iBand) = skyBr(iBand) / 1.02d-15  ! Convert sky brightness from erg/s/cm2/μ/arcsec^2 to nanolamberts.  Should this be 1.02d-15 iso 1.11d-15?  Or 4.54d-16? - http://www.archaeocosmology.org/eng/vislimitbas.htm
       
    end do  ! iBand
    
  end subroutine limmag_skyBrightness
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate limiting magnitude based on JD and object altitude, wrapper for limmag_full()
  !!
  !! \param  jd         Julian day
  !! \param  objRA      Right ascension of the observed object (rad)
  !! \param  objDec     Declination of the observed object (rad)
  !! \param  objAlt     Altitude of the observed object (rad)
  !!
  !! \param  lat        Latitude of the observer (optional; rad)
  !! \param  height     Height/altitude of the observer above sea level (optional; metres)
  !!
  !! \retval limmag_jd  Limiting magnitude
  !!
  !! \note  Using observer's location from module TheSky_local; overruled by optional dummy variables lat and height
  
  function limmag_jd(jd, objRA,objDec,objAlt, lat,height)
    use SUFR_kinds, only: double
    use SUFR_angles, only: asep
    use SUFR_date_and_time, only: jd2cal
    
    use TheSky_planets, only: planet_position_la
    use TheSky_planetdata, only: planpos, nplanpos
    use TheSky_local, only: lat0, height0=>height
    
    implicit none
    real(double), intent(in), value :: jd, objRA,objDec,objAlt  ! Often called as planpos(i) -> pass by value
    real(double), intent(in), optional :: lat,height
    integer :: year, month
    real(double) :: limmag_jd,  llat,lheight,  day, sunAlt,sunElon, moonPhase,moonAlt,moonElon, planpos0(nplanpos)
    
    ! Use values from TheSky_local by default:
    llat = lat0
    lheight = height0
    
    ! Overrule with values from optional arguments if present:
    if(present(lat)) llat = lat
    if(present(height)) lheight = height
    
    
    call jd2cal(jd, year,month,day)
    
    planpos0 = planpos                                        ! Save current contents
    
    call planet_position_la(jd, 3, 4, 0)                      ! Sun position - 4: need alitude
    sunAlt = planpos(31)                                      ! Sun altitude
    sunElon = asep(objRA, planpos(5),  objDec, planpos(6))    ! Sun elongation
    
    call planet_position_la(jd, 0, 5, 60)                     ! Moon position - 4: need k;  60 terms
    moonPhase = planpos(14)                                   ! Moon phase (illuminated fraction, k)
    moonAlt   = planpos(31)                                   ! Moon altitude
    moonElon  = asep(objRA, planpos(5),  objDec, planpos(6))  ! Moon elongation
    
    limmag_jd = limmag_full(year,month, lheight,llat, sunAlt,sunElon, moonPhase,moonAlt,moonElon, objAlt)
    
    planpos = planpos0                                        ! Restore original contents
    
  end function limmag_jd
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate limiting magnitude for the local zenith, based on JD and observer's location
  !!
  !! \param  jd         Julian day
  !! \param  lat        Latitude of the observer (rad)
  !! \param  lon        Longitude of the observer (rad)
  !! \param  height     Height/altitude of the observer above sea level (metres)
  !!
  !! \retval limmag_zenith_jd  Limiting magnitude
  
  function limmag_zenith_jd(jd, lat,lon,height)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pio2
    
    use TheSky_datetime, only: calc_gmst
    use TheSky_coordinates, only: horiz2eq
    
    implicit none
    real(double), intent(in) :: jd, lat,lon,height
    real(double) :: limmag_zenith_jd,  objAlt, gmst,  hh,objRA,objDec
    
    ! Compute the RA and dec corresponding to the local zenith at JD, lat and lon:
    objAlt = pio2  ! Altitude of the zenith = pi/2
    gmst = calc_gmst(jd)  ! Greenwich mean sidereal time
    call horiz2eq(0.d0,objAlt, gmst,  hh,objRA,objDec,  lat,lon)  ! Use azimuth = 0, gmst iso agst
    
    ! Compute the limiting magnitude for the desired JD and the RA, dec and alt corresponding to the zenith:
    limmag_zenith_jd = limmag_jd(jd, objRA,objDec, objAlt, lat,height)
    
  end function limmag_zenith_jd
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate limiting magnitude based on JD and planet ID, wrapper for limmag_jd()
  !!
  !! \param  jd            Julian day
  !! \param  pl            Planet ID
  !!
  !! \retval limmag_jd_pl  Limiting magnitude
  !!
  !! \note  Using observer's location from module TheSky_local
  
  function limmag_jd_pl(jd, pl)
    use SUFR_kinds, only: double
    
    use TheSky_planets, only: planet_position_la
    use TheSky_planetdata, only: planpos, nplanpos
    
    implicit none
    real(double), intent(in) :: jd
    integer, intent(in) :: pl
    real(double) :: limmag_jd_pl,  objRA,objDec,objAlt, planpos0(nplanpos)
    
    planpos0 = planpos                     ! Save
    
    call planet_position_la(jd, pl, 6,10)  ! Planet position
    objRA    = planpos(5)                  ! Object right ascension
    objDec   = planpos(6)                  ! Object declination
    objAlt   = planpos(31)                 ! Object altitude - refracted
    
    limmag_jd_pl = limmag_jd(jd, objRA,objDec,objAlt)
    
    planpos = planpos0                     ! Restore
    
  end function limmag_jd_pl
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate limiting magnitude based on Sun altitude and object altitude (airmass).  Assumes New Moon.
  !!
  !! \param  month      Month of the year
  !! \param  sunAlt     Altitude of the Sun (rad)
  !! \param  sunAz      Azimuth of the Sun (rad)
  !! \param  objAlt     Altitude of the observed object (rad)
  !! \param  objAz      Azimuth of the observed object (rad)
  !!
  !! \param  lat        Latitude of the observer (optional; rad)
  !! \param  height     Height/altitude of the observer above sea level (optional; metres)
  !!
  !! \retval limmag_sun_airmass  Limiting magnitude
  !!
  !! \note  Using observer's location from module TheSky_local; overruled by optional dummy variables lat and height
  
  function limmag_sun_airmass(month, sunAlt,sunAz, objAlt,objAz, lat,height)
    use SUFR_kinds, only: double
    use SUFR_constants, only: d2r
    use SUFR_angles, only: asep
    use SUFR_date_and_time, only: jd2cal
    
    use TheSky_local, only: lat0, height0=>height
    
    implicit none
    integer, intent(in) :: month
    real(double), intent(in), value :: sunAlt,sunAz, objAlt,objAz  ! Often called as planpos(i) -> pass by value
    real(double), intent(in), optional :: lat,height
    integer :: year
    real(double) :: limmag_sun_airmass,  llat,lheight,  sunElon, moonPhase,moonAlt,moonElon
    
    ! Use values from TheSky_local by default:
    llat = lat0
    lheight = height0
    
    ! Overrule with values from optional arguments if present:
    if(present(lat)) llat = lat
    if(present(height)) lheight = height
    
    year = 1998  ! Average year for solar activity (1992+5.5)
    
    sunElon = asep(objAz, sunAz,  objAlt, sunAlt)        ! Sun elongation
    
    ! No Moon:
    moonPhase = 0.d0
    moonAlt   = -89*d2r
    moonElon  = 179*d2r
    
    limmag_sun_airmass = limmag_full(year,month, lheight,llat, sunAlt,sunElon, moonPhase,moonAlt,moonElon, objAlt)
    
  end function limmag_sun_airmass
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief Calculate limiting magnitude, based on the altitude of the Sun.  Simplified version of limmag_full()
  !!
  !! \param  sunAlt  Altitude of the Sun (rad)
  !!
  !! \retval limmag_sun  Limiting magnitude
  !!
  !! \note 
  !! - Mag depends only on sunAlt in rad
  !! - assume object in zenith, no moon (am=180), humidity 50%, T=10degC, lat=45deg, sn=1
  !! 
  !! \see http://www.go.ednet.ns.ca/~larry/astro/vislimit.html, 
  !!      a JavaScript version of a BASIC program by Bradley E. Schaefer, Sky and Telescope, May 1998, p.57
  
  function limmag_sun(sunAlt)
    use SUFR_kinds, only: double
    use SUFR_constants, only: r2d
    
    implicit none
    real(double), intent(in) :: sunAlt
    real(double) :: limmag_sun,  dm,bn,bt,bd,bl,c1,c2,th
    
    dm = 0.285667465769191d0         ! Extinction
    bn = 7.686577723230466d-14       ! Dark night sky brightness: no solar cycle, star in zenith
    bt = 10.d0**(-0.4d0*(-26.74d0 + 11.05d0 + 32.5d0 - sunAlt*r2d )) * 0.12852345982053d0  ! Twilight brightness (erg/s/cm^2/μ/arcsec^2)
    bd = 9.456022312552874d-7        ! Daylight brightness
    bl = (bn + min(bd,bt))/1.02d-15  ! Total sky brightness in V, convert from erg/s/cm2/μ/arcsec^2 to nanolamberts.  Should this be 1.02d-15 iso 1.11d-15?  Or 4.54d-16? - http://www.archaeocosmology.org/eng/vislimitbas.htm
    
    ! Visual limiting magnitude
    if(bl.lt.1500.d0) then
       c1 = 10.d0**(-9.8d0)
       c2 = 10.d0**(-1.9d0)
    else
       c1 = 10.d0**(-8.350001d0)
       c2 = 10.d0**(-5.9d0)
    end if
    th = c1*(1.d0 + sqrt(c2*bl))**2
    
    limmag_sun = -16.57d0 - 2.5d0*log10(th) - dm   ! Limiting magnitude
    
  end function limmag_sun
  !*********************************************************************************************************************************
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the excess magnitude for planet pl at JD, considering Sun, Moon and airmass
  !!
  !! \param jd  Julian day for moment of interest
  !! \param pl  Planet ID (0-Moon, 1-Mer, 8-Nep, >10-comet
  !!
  !! \note
  !! - The excess magnitude is defined as  the magnitude of an object  MINUS  the limiting magnitude
  !!   - xsmag < 0 - object is visible (in theory!) with the naked eye
  !! - The magnitude is corrected for extinction due to airmass, but for sea level
  !! - The limiting magnitude here solely depends on the Sun's altitude, simplified expression (especially, no Moon!)
  !! - This routine should be used as an indication only - no hard facts...
  
  function pl_xsmag(jd, pl)
    use SUFR_kinds, only: double
    
    use TheSky_planets, only: planet_position_la
    use TheSky_planetdata, only: planpos
    
    implicit none
    real(double), intent(in) :: jd
    integer, intent(in) :: pl
    real(double) :: pl_xsmag,  objRA,objDec,objAlt
    
    call planet_position_la(jd, pl, 4,60)  ! Compute planet position
    objRA  = planpos(5)                    ! Object right ascension
    objDec = planpos(6)                    ! Object declination
    objAlt = planpos(31)                   ! Object altitude - refracted
    
    if(objAlt.lt.0.d0) then  ! Object is below the horizon
       pl_xsmag = 99.d0 - objAlt  ! Huge contrast, and worse if further below the horizon, for solvers
    else
       pl_xsmag = planpos(13) - limmag_jd(jd, objRA,objDec,objAlt)  ! Magnitude - limiting magnitude
    end if
    
  end function pl_xsmag
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the excess magnitude at JD, wrapper for pl_xsmag() for solvers.
  !!         The planet ID pl0 is passed through module planetdata.
  !!
  !! \param jd  Julian day for moment of interest
  
  function pl_xsmag_pl(jd)
    use SUFR_kinds, only: double
    use TheSky_planetdata, only: pl0
    
    implicit none
    real(double), intent(in) :: jd
    real(double) :: pl_xsmag_pl
    
    pl_xsmag_pl = pl_xsmag(jd,pl0)
    
  end function pl_xsmag_pl
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the excess magnitude for planet pl at JD, considering airmass and Sun alt - low-accuracy version of pl_xsmag()
  !!
  !! \param jd  Julian day for moment of interest
  !! \param pl  Planet ID (0-Moon, 1-Mer, 8-Nep, >10-comet
  !!
  !! \note
  !! - The excess magnitude is defined as  the magnitude of an object  MINUS  the limiting magnitude
  !!   - xsmag < 0 - object is visible (in theory!) with the naked eye
  !! - The magnitude is corrected for extinction due to airmass, but for sea level
  !! - The limiting magnitude here solely depends on the Sun's altitude, simplified expression (especially, no Moon!)
  
  function pl_xsmag_la(jd, pl)
    use SUFR_kinds, only: double
    use TheSky_planets, only: planet_position_la
    use TheSky_planetdata, only: planpos
    
    implicit none
    real(double), intent(in) :: jd
    integer, intent(in) :: pl
    real(double) :: pl_xsmag_la, extinction_magPam, sunAlt
    
    call planet_position_la(jd,  3, 4,0)  ! Sun position
    sunAlt = planpos(31)                  ! Refracted
    
    call planet_position_la(jd, pl, 0,0)  ! Planet position
    
    if(planpos(31).lt.0.d0) then  ! Object is below the horizon
       pl_xsmag_la = 99.d0 - planpos(31) ! Huge contrast, and worst if further below the horizon, for solvers
    else
       extinction_magPam = 0.2811d0  ! Extinction in magnitudes per unit airmass, at sea level
       pl_xsmag_la = (planpos(13) + extinction_magPam * airmass(planpos(30))) - limmag_sun(sunAlt)  ! (m + ext) - limmag; TRUE alt!
    end if
    
  end function pl_xsmag_la
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the excess magnitude at JD, wrapper for pl_xsmag_la() for solvers.
  !!         The planet ID pl0 is passed through module planetdata.
  !!
  !! \param jd  Julian day for moment of interest
  
  function pl_xsmag_la_pl(jd)
    use SUFR_kinds, only: double
    use TheSky_planetdata, only: pl0
    
    implicit none
    real(double), intent(in) :: jd
    real(double) :: pl_xsmag_la_pl
    
    pl_xsmag_la_pl = pl_xsmag_la(jd,pl0)
    
  end function pl_xsmag_la_pl
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Aperture diameter in centimetres needed to observe an object with given excess magnitude (0: no instrument needed)
  !!
  !! \param xsmag  Excess magnitude:  magnitude of an object  MINUS  limiting magnitude  (>0: too weak for naked eye)
  !! \param pupil  Pupil diameter (mm, default: 6)
  !! \param tc     Transmission coefficient of the instrument (default: 0.8 = 80%)
  !!
  !! \note  This routine should be used as an indication only - no hard facts...
  !!
  !! \see
  !! - http://iovs.arvojournals.org/article.aspx?articleid=2161149 - Fig.3a
  !!   - average pupil diameter at low light: 7 -> 5mm for 20- -> 80-year olds
  !!   - use 6mm as default value
  !! - http://www.telescope-optics.net/functions.htm
  !!   - default tc for amateur telescopes: ~80%
  
  function aperture(xsmag, pupil, tc)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: xsmag
    real(double), intent(in), optional :: pupil, tc
    real(double) :: aperture, lpupil, ltc
    
    ! Pupil size:
    lpupil = 6.d0  ! Default pupil diameter: 6mm
    if(present(pupil)) lpupil = pupil
    
    ! Light-transmission coefficient:
    ltc = 0.8d0  ! Default tc for amateur telescopes ~80% (http://www.telescope-optics.net/functions.htm)
    if(present(tc)) ltc = tc
    
    if(xsmag.le.0.d0) then
       aperture = 0.d0
    else
       aperture = lpupil * sqrt(10.d0**(xsmag/2.5d0) / ltc) / 10.d0  ! Aperture diameter in cm
    end if
    
  end function aperture
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert naked-eye visual limiting magnitude (V) to sky surface brightness in (B) magnitudes per square arcsecond
  !!
  !! \param  Mlim                          Naked-eye visual limiting magnitude (V)
  !! \retval skybrightness_mas2_from_mlim  Sky surface brightness in (B) magnitudes per square arcsecond
  !!
  !! \see   http://adsabs.harvard.edu/abs/1990PASP..102..212S
  !! \note  Sky surface brightness is sometimes referred to as sqm (which is in fact a device to measure it)
  
  elemental function skybrightness_mas2_from_mlim(Mlim)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: Mlim
    real(double) :: skybrightness_mas2_from_mlim
    
    if(Mlim.lt.7.93d0) then
       skybrightness_mas2_from_mlim = 21.58d0 - 5*log10(10.d0**(1.586d0-Mlim/5.d0) - 1.d0)
    else
       skybrightness_mas2_from_mlim = 99.99d0
    end if
        
  end function skybrightness_mas2_from_mlim
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert sky surface brightness in (B) magnitudes per square arcsecond to naked-eye visual limiting magnitude (V)
  !!
  !! \param  skyBright                     Sky surface brightness in (B) magnitudes per square arcsecond
  !! \retval skybrightness_mlim_from_mas2  Naked-eye visual limiting magnitude (V)
  !!
  !! \see   http://adsabs.harvard.edu/abs/1990PASP..102..212S
  !! \note  Sky surface brightness is sometimes referred to as sqm (which is in fact a device to measure it)
  
  elemental function skybrightness_mlim_from_mas2(skyBright)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: skyBright
    real(double) :: skybrightness_mlim_from_mas2
    
    skybrightness_mlim_from_mas2 = 7.93d0 - 5*log10(10.d0**(4.316d0-skyBright/5.d0) + 1.d0)
    
  end function skybrightness_mlim_from_mas2
  !*********************************************************************************************************************************
  
  !*********************************************************************************************************************************
  !> \brief  Convert sky brightness in candela per square meter to naked-eye visual limiting magnitude (V)
  !!
  !! \param  cdm2    Sky brightness in cd/m^2
  !! \retval Naked-eye visual limiting magnitude (V)
  !!
  !! \note
  !!    Adapted from Weaver (1947) and Garsteng (1986):
  !!    - shift to match data by Weaver (1947): 7.930 -> 7.706;  4.305 -> 4.115;
  !!    - use candela per square meter instead of nanolambert.
  
  elemental function skybrightness_mlim_from_cdm2(cdm2)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: cdm2
    real(double) :: skybrightness_mlim_from_cdm2, nL, Mlim
    
    nL = skybrightness_nL_from_cdm2(cdm2)  ! Convert from candela per square meter to nanolambert
    
    Mlim = 7.706d0 - 5*log10(1.d0 + 0.1122d0   * sqrt(nL))                      ! b <= 1479 nL = Mlim>4.0
    if(nL.gt.1479) Mlim = 4.115d0 - 5*log10(1.d0 + 0.001122d0 * sqrt(nL))    ! b >  1479 nL = Mlim<4.0
    
    skybrightness_mlim_from_cdm2 = Mlim
  end function skybrightness_mlim_from_cdm2
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert limiting visual magnitude to sky brightness in candela per square meter
  !!
  !! \param mlim  Naked-eye visual limiting magnitude (V)
  !! \retval      Sky brightness in cd/m^2
  !!
  !! \note  Inverse of skybrightness_mlim_from_cdm2()
  
  elemental function skybrightness_cdm2_from_mlim(Mlim)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: Mlim
    real(double) :: skybrightness_cdm2_from_mlim, nL, cdm2
    
    nL = ((10.d0**(-(mlim - 7.706d0)/5.d0) - 1.d0)/0.1122d0)**2                      ! Mlim >= 4.0; b <= 1479 nL
    if(mlim.lt.4) nL = ((10.d0**(-(mlim  - 4.115d0)/5.d0) - 1.d0)/0.001122d0)**2  ! Mlim < 4.0;  b >  1479 nL
    
    cdm2 = skybrightness_cdm2_from_nL(nL)  ! Convert nanolambert to candela per square meter
    
    skybrightness_cdm2_from_mlim = cdm2
  end function skybrightness_cdm2_from_mlim
  !*********************************************************************************************************************************
  
    
  
  !*********************************************************************************************************************************
  !> \brief  Convert sky brightness in nanolambert to candela per square meter
  !!
  !! \param  nL  Sky brightness in nanolambert (nL)
  !! \retval     Sky brightness in candela per square meter (cd/m^2)
  
  elemental function skybrightness_cdm2_from_nL(nL)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi
    
    implicit none
    real(double), intent(in) :: nL
    real(double) :: skybrightness_cdm2_from_nL
    
    skybrightness_cdm2_from_nL = nL * 1d-5/pi
  end function skybrightness_cdm2_from_nL
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert sky brightness in candela per square meter to nanolambert
  !!
  !! \param  cdm2  Sky brightness in candela per square meter (cd/m^2)
  !! \retval       Sky brightness in nanolambert (nL)
  
  elemental function skybrightness_nL_from_cdm2(cdm2)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi
    
    implicit none
    real(double), intent(in) :: cdm2
    real(double) :: skybrightness_nL_from_cdm2
    
    skybrightness_nL_from_cdm2 = cdm2 * 1d5 * pi
  end function skybrightness_nL_from_cdm2
  !*********************************************************************************************************************************
  
  
end module TheSky_visibility
!***********************************************************************************************************************************
