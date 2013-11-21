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
!  limmag_sunmoon                       Calculate limiting magnitude based on Sun and Moon
!  limmag                               Calculate limiting magnitude, based on the altitude of the Sun.  Simplif'd limmag_sunmoon()



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
  
  
  !*********************************************************************************************************************************
  !> \brief Compute the airmass for a celestial object with a given altitude
  !!
  !! \param alt  Altitude of object (radians)
  !!
  !! - Results are 1 <= airmass <~ 38; return 1.d99 for h<0
  !!
  !! \see Kasten and Young (1989); http://en.wikipedia.org/wiki/Airmass#Interpolative_formulas
  
  function airmass(alt)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pio2, r2d
    
    implicit none
    real(double), intent(in) :: alt
    real(double) :: airmass,z,zdeg
    
    if(alt.lt.0.d0) then
       airmass = 1000.d0 * (0.15d0 + abs(alt))  ! Very bad (adds at least ~30 magnitudes due to extinction), 
       !                                          but still worse when farther below the horizon
    else
       z = min(pio2 - alt, pio2)  ! Zenith angle
       zdeg = z*r2d
       airmass = max( 1.d0 / ( cos(z) + 0.50572d0*(96.07995d0-zdeg)**(-1.6364d0) ) ,  1.d0 )
    end if
    
  end function airmass
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the extinction in magnitdes per unit airmass for an observer with given elevation
  !!
  !! \param ele  Evelation of the observer above sea level (metres)
  !!
  !! \note  The magnitude of an object corrected for airmass should be  m' = m + airmass_ext(ele) * airmass(alt)
  !!
  !! \see  Green, ICQ 14, 55 (1992),  http://www.icq.eps.harvard.edu/ICQExtinct.html
  !!
  
  function airmass_ext(ele)
    use SUFR_kinds, only: double
    
    implicit none
    real(double), intent(in) :: ele
    real(double):: airmass_ext, Aoz,Aray,Aaer
    
    Aoz  = 0.016d0                       ! Ozone  (Schaefer 1992)
    Aray = 0.1451d0 * exp(-ele/7996.d0)  ! Rayleigh scattering, Eq.2
    Aaer = 0.120d0  * exp(-ele/1500.d0)  ! Aerosol scattering, Eq.4
    
    airmass_ext = Aoz + Aray + Aaer      ! Total extinction in magnitudes per unit air mass
    
  end function airmass_ext
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate limiting magnitude based on Sun altitude and Moon phase
  !!
  !! \param  year        Year of observation
  !! \param  month       Month of observation
  !! \param  obselev     Elevation of the observer (m)
  !! \param  obslat      Latitude of the observer (rad)
  !!
  !! \param  sunalt      Altitude of the Sun (rad)
  !! \param  sunelon     Elongation object-Sun (rad)
  !!
  !! \param  moonphase   Phase of the Moon (fraction?)
  !! \param  moonalt     Altitude of the Moon (rad)
  !! \param  moonelon    Elongation object-Moon (rad)
  !!
  !! \param  objalt      Altitude of the observed object (rad)
  !!
  !! \retval limmag_sunmoon  Limiting magnitude
  !!
  !! \see http://www.go.ednet.ns.ca/~larry/astro/vislimit.html, 
  !!      a JavaScript version of a BASIC program by Bradley E. Schaefer, Sky and Telescope, May 1998, p.57
  
  function limmag_sunmoon(year,month, obselev,obslat, sunalt,sunelon, moonphase,moonalt,moonelon, objalt)
    use SUFR_kinds, only: double
    use SUFR_constants, only: d2r,r2d
    
    implicit none
    integer, intent(in) :: year, month
    real(double), intent(in) :: obselev,obslat, sunalt,sunelon, moonphase,moonalt,moonelon, objalt
    
    integer :: i,m,y
    real(double) :: limmag_sunmoon,  b(5),k(5),dm(5),wa(5),mo(5),oz(5),wt(5),bo(5),cm(5),ms(5),  am,zm,rm,zs,rs,rh,te,la,al,sn,z
    real(double) :: lt,ra,sl,zz,xg,xa,xo,kr,ka,ko,kw,  x,xm,xs,bn,mm,c3,fm,bm,hs,bt,c4,fs,bd,bl,c1,c2,th,mn
    
    
    b = 0.d0;  k = 0.d0;  dm = 0.d0
    wa = (/0.365d0, 0.44d0, 0.55d0, 0.7d0, 0.9d0/)
    mo = (/-10.93d0, -10.45d0, -11.05d0, -11.90d0, -12.70d0/)
    oz = (/0.000d0, 0.000d0, 0.031d0, 0.008d0, 0.000d0/)
    wt = (/0.074d0, 0.045d0, 0.031d0, 0.020d0, 0.015d0/)
    bo = (/8.0d-14, 7.0d-14, 1.0d-13, 1.0d-13, 3.0d-13/)
    cm = (/1.36d0, 0.91d0, 0.00d0, -0.76d0, -1.17d0/)
    ms = (/-25.96d0, -26.09d0, -26.74d0, -27.26d0, -27.55d0/)
    
    am = (1.d0 - moonphase)*180  ! Phase angle moon (deg)
    zm = 90.d0 - moonalt*r2d     ! Zenith angle moon (deg)
    rm = moonelon*r2d            ! Elongation moon-star (deg)
    zs = 90.d0 - sunalt*r2d      ! Zenith angle sun (deg)
    rs = sunelon*r2d             ! Elongation sun-star (deg)
    rh = 50.d0                   ! Relative humidity (%)
    te = 10.d0                   ! Temperature (deg.C)
    la = obslat*r2d              ! Observer's latitude (deg)
    al = obselev                 ! Observer's elevation (m)
    m  = month                   ! Month (for approximate sun ra)
    y  = year                    ! Year (for solar cyle...)
    sn = 1.d0                    ! Snellen ratio quantifies vision
    z  = 90.d0 - objalt*r2d      ! Zenith angle object (deg)
    
    
    !*** Extinction ************************************************************
    lt = la*d2r
    ra = (m-3)*30*d2r  ! Derive rough Sun RA from month
    sl = la/abs(la)
    
    ! Airmass for each component:
    zz = z*d2r
    xg = 1.d0/(cos(zz) + 0.0286d0*exp(-10.5d0*cos(zz)))           ! Gas
    xa = 1.d0/(cos(zz) + 0.0123d0*exp(-24.5d0*cos(zz)))           ! Aerosol
    xo = 1.d0/sqrt(1.d0 - (sin(zz)/(1.d0 + (20.d0/6378.d0)))**2)  ! Ozone
    
    ! UBVRI extinction for each component:
    do i=1,5
       kr = 0.1066d0*exp(-al/8200.d0)*(wa(i)/0.55d0)**(-4)
       ka = 0.1d0*(wa(i)/0.55d0)**(-1.3) * exp(-al/1500.d0)
       ka = ka * (1.d0-0.32d0/log(rh/100.d0))**1.33d0 * (1.d0 + 0.33d0*sl*sin(ra))
       ko = oz(i)*(3.d0 + 0.4d0*(lt*cos(ra)-cos(3*lt)))/3.d0
       kw = wt(i)*0.94d0*(rh/100.d0)*exp(te/15.d0)*exp(-al/8200.d0)
       k(i) = kr + ka + ko + kw
       dm(i) = kr*xg + ka*xa + ko*xo + kw*xg
    end do
    
    
    !*** Sky *******************************************************************
    x  = 1.d0/(cos(zz) + 0.025d0*exp(-11*cos(zz)))
    xm = 1.d0/(cos(zm*d2r) + 0.025d0*exp(-11*cos(zm*d2r)))
    if(zm.gt.90.d0) xm = 40.d0
    xs = 1.d0/(cos(zs*d2r) + 0.025d0*exp(-11*cos(zs*d2r)))
    if(zs.gt.90.d0) xs = 40.d0
    
    do i=1,5
       ! Dark night sky brightness:
       bn = bo(i) * (1.d0 + 0.3d0*cos(6.283d0*dble(y-1992)/11.d0))
       bn = bn * (0.4d0 + 0.6d0/sqrt(1.d0-0.96d0*sin(zz)**2))
       bn = bn * 10.d0**(-0.4d0*k(i)*x)
       
       ! Moonlight brightness:
       mm = -12.73d0 + 0.026d0*abs(am) + 4.d-9*am**4  ! Moon magnitude
       mm = mm + cm(i)
       c3 = 10.d0**(-0.4d0*k(i)*xm)
       fm = 6.2d7/(rm*rm) + 10.d0**(6.15d0 - rm/40.d0)
       fm = fm + 10.d0**5.36d0 * (1.06d0 + cos(rm*d2r)**2)
       bm = 10.d0**(-0.4d0*(mm-mo(i)+43.27d0))
       bm = bm * (1.d0 - 10.d0**(-0.4d0*k(i)*x))
       bm = bm * (fm*c3 + 440000.d0*(1.d0-c3))
       
       ! Twilight brightness:
       hs = 90.d0 - zs
       bt = 10.d0**(-0.4d0*(ms(i) - mo(i) + 32.5d0 - hs - (z/(360.d0*k(i))) ))
       bt = bt* (100.d0/rs) * (1.d0 - 10.d0**((-0.4d0*k(i))*x))
       
       ! Daylight brightness:
       c4 =  10.d0**(-0.4d0*k(i)*xs)
       fs = 6.2d7/(rs*rs) + 10.d0**(6.15d0 - rs/40.d0)
       fs = fs + 10.d0**5.36d0 * (1.06d0 + cos(rs*d2r)**2)
       bd = 10.d0**(-0.4d0*(ms(i) -mo(i) + 43.27d0))
       bd = bd * (1.d0 - 10.d0**(-0.4d0*k(i)*x))
       bd = bd * (fs*c4 + 440000.d0*(1.d0-c4))
       
       ! Total sky brightness:
       if(bd.lt.bt) then
          b(i) = bn + bd
       else
          b(i) = bn + bt
       end if
       
       if(zm.lt.90.d0) b(i) = b(i) + bm
       b(i) = b(i)*1.d12  ! Convert from ergs to picoergs
    end do
    
    ! Visual limiting magnitude:
    bl = b(3)/1.11d-3           ! Sky brightness in V, convert to nanolamberts
    if(bl.lt.1500.d0) then
       c1 = 10.d0**(-9.8d0)
       c2 = 10.d0**(-1.9d0)
    else
       c1 = 10.d0**(-8.350001d0)
       c2 = 10.d0**(-5.9d0)
    end if
    th = c1*(1.d0 + sqrt(c2*bl))**2
    mn = -16.57d0 - 2.5d0*log10(th) - dm(3) + 5*log10(sn)  ! Magnitude
    
    limmag_sunmoon = mn
    
  end function limmag_sunmoon
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief Calculate limiting magnitude, based on the altitude of the Sun.  Simplified version of limmag_sunmoon()
  !!
  !! \param  sunalt  Altitude of the Sun (rad)
  !!
  !! \retval limmag  Limiting magnitude
  !!
  !! \note 
  !! - Mag depends only on sunalt in rad
  !! - assume object in zenith, no moon (am=180), humidity 50%, T=10degC, lat=45deg, sn=1
  !! 
  !! \see http://www.go.ednet.ns.ca/~larry/astro/vislimit.html, 
  !!      a JavaScript version of a BASIC program by Bradley E. Schaefer, Sky and Telescope, May 1998, p.57
  
  
  function limmag(sunalt)
    use SUFR_kinds, only: double
    use SUFR_constants, only: r2d
    
    implicit none
    real(double), intent(in) :: sunalt
    real(double) :: limmag
    real(double) :: dm,bn,bt,bd,bl,c1,c2,th,mn
    
    dm = 0.285667465769191d0     ! Extinction
    bn = 7.686577723230466d-14   ! Dark night sky brightness: no solar cycle, star in zenith
    bt = 10.d0**(-0.4d0*(-26.74d0 + 11.05d0 + 32.5d0 - sunalt*r2d )) * 0.12852345982053d0  ! Twilight brightness
    bd = 9.456022312552874d-7    ! Daylight brightness
    bl = 1.d12*(bn + min(bd,bt))/1.11d-3   ! Total sky brightness in V, convert to nanolamberts
    
    ! Visual limiting magnitude
    if(bl.lt.1500.d0) then
       c1 = 10.d0**(-9.8)
       c2 = 10.d0**(-1.9)
    else
       c1 = 10.d0**(-8.350001)
       c2 = 10.d0**(-5.9)
    end if
    th = c1*(1.d0 + sqrt(c2*bl))**2
    mn = -16.57d0 - 2.5d0*log10(th) - dm   ! Magnitude
    
    limmag = mn
    
  end function limmag
  !*********************************************************************************************************************************
  
  
  
  
  
end module TheSky_visibility
!***********************************************************************************************************************************
