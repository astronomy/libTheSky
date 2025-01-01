!> \file planets.f90  Procedures to compute planet positions and more for libTheSky
 

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
!> \brief Procedures for planets

module TheSky_planets
  implicit none
  save
  
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute the position, distance, etc of a planet
  !!
  !! \param jd          Julian date of computation
  !! \param pl          Planet number: Moon=0, Merc=1,Nep=8,Plu=9 3=Sun, 0=Moon, >10 for other objects
  !!
  !! \param lat         Latitude of the observer (rad, optional)
  !! \param lon         Longitude of the observer (rad, optional)
  !! \param hgt         Altitude/elevation of the observer above sea level (metres, optional)
  !! 
  !! \param LBaccur     Desired accuracy of the heliocentric L,B in VSOP87 (rad, optional; defaults to full accuracy)
  !! \param Raccur      Desired accuracy of the heliocentric R in VSOP87 (AU, optional; defaults to full accuracy)
  !!
  !! \param ltime       Correct for light time (optional; defaults to true). .false. saves ~50% in CPU time at the cost of some accuracy.
  !! \param aber        Correct for aberration (optional; defaults to true).
  !! \param to_fk5      Convert coordinates to the FK5 system (optional; defaults to true).
  !! \param assume_jde  Assume JD provided is actually JDE (optional; defaults to false).
  !! 
  !! \param lunar_theory  Choose Lunar theory:  1: ELP82b,  2: ELP-MPP02/LLR,  3: ELP-MPP02/DE405 ('historical' - default)
  !! \param nutat         IAU nutation model to use: 0 (no nutation!), 1980 or 2000 (default).
  !! \param magmdl        Planet-magnitude model: 1: Müller (1893), 2: Meeus p.286, 3: AA 1992, 4: AA 2012 (optional, defaults to 4).
  !! \param verbosity     Verbosity for debug output (0-3).  Defaults to 0: silent.
  !! 
  !!
  !! \note
  !! - lat0 and lon0 can be provided through the module TheSky_local (rad, rad and m), or through the optional arguments.
  !!   Note that using the latter will update the former!
  !! - results are returned in the array planpos() in the module TheSky_planetdata
  !!
  
  subroutine planet_position(jd,pl, lat,lon,hgt, LBaccur,Raccur, ltime,aber,to_fk5,assume_jde, &
    lunar_theory,nutat,magmdl, verbosity)
    
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi,pi2, r2d,r2h, au,earthr,pland, enpname, jd2000
    use SUFR_system, only: warn, quit_program_error
    use SUFR_angles, only: rev, rev2
    use SUFR_dummy, only: dumdbl1,dumdbl2
    use SUFR_text, only: d2s
    use SUFR_numerics, only: deq0
    
    use TheSky_vsop, only: vsop87d_lbr
    use TheSky_local, only: lon0,lat0,height, deltat
    use TheSky_coordinates, only: ecl_2_eq, eq2horiz, hc_spher_2_gc_rect, geoc2topoc_ecl, rect_2_spher
    use TheSky_coordinates, only: precess_ecl, fk5, aberration_ecl, refract
    use TheSky_nutation, only: nutation, nutation2000
    use TheSky_sun, only: sunmagn
    use TheSky_moon, only: elp82b_lbr, elp_mpp02_lbr, moonmagn
    use TheSky_comets, only: cometgc
    use TheSky_planetdata, only: planpos, pl0
    use TheSky_cometdata, only: cometElems, cometDiedAtP
    use TheSky_asteroids, only: asteroid_magn, asteroid_lbr
    use TheSky_datetime, only: calc_deltat, calc_gmst
    
    implicit none
    real(double), intent(in) :: jd
    integer, intent(in) :: pl
    real(double), intent(in), optional :: lat,lon,hgt, LBaccur,Raccur
    logical, intent(in), optional :: ltime,aber,to_fk5, assume_jde
    integer, intent(in), optional :: lunar_theory, nutat,magmdl, verbosity
    
    integer :: j, llunar_theory, lnutat,lmagmdl, lverb
    real(double) :: jdl,jde,jde_lt, tjm,tjm0,  llat,llon,lhgt, lLBaccur,lRaccur,  dpsi,eps0,deps,eps,tau,tau1
    real(double) :: hcl0,hcb0,hcr0, hcl,hcb,hcr, hcl00,hcb00,hcr00, sun_gcl,sun_gcb, gcx,gcy,gcz, gcx0,gcy0,gcz0, dhcr
    real(double) :: gcl,gcb,delta,gcl0,gcb0,delta0
    real(double) :: ra,dec,gmst,agst,lst,hh,az,alt,elon,  topra,topdec,topl,topb,topdiam,topdelta,tophh,topaz,topalt
    real(double) :: diam,illfr,pa,magn,  parang,parang_ecl,hp, rES1,rES2
    logical :: lltime,laber,lto_fk5
    
    ! Make sure these are always defined:
    gcl0 = 0.d0;  gcb0 = 0.d0;  hcr00 = 0.d0
    gcx  = 0.d0;  gcy  = 0.d0;  gcz  = 0.d0
    gcx0 = 0.d0;  gcy0 = 0.d0;  gcz0 = 0.d0
    topdelta = 0.d0;  delta0 = 0.d0
    diam = 0.d0
    magn = 0.d0  ! Make sure magn is always defined
    
    pl0 = pl                           ! Remember which planet was computed last
    
    
    ! Handle optional variables:
    ! Geographic location:
    llat = lat0
    llon = lon0
    lhgt = height
    if(present(lat)) llat = lat
    if(present(lon)) llon = lon
    if(present(hgt)) lhgt = hgt
    
    ! Coordinate accuracy:
    lLBaccur = 0.d0  ! Seting L,B,R accuracy equal to VSOP87 truncation still skips some terms!
    lRaccur = 0.d0
    if(present(LBaccur)) lLBaccur = LBaccur
    if(present(Raccur))  lRaccur  = Raccur
    
    ! Light time, aberration and FK5 switches:
    lltime = .true.                       ! Correct for light time by default
    if(present(ltime)) lltime = ltime
    
    laber = .true.                        ! Correct for aberration by default
    if(present(aber)) laber = aber
    
    lto_fk5 = .true.                      ! Convert to FK5 by default
    if(present(to_fk5)) lto_fk5 = to_fk5
    
    ! Select lunar theory:
    ! llunar_theory = 1  ! Use (corrected) ELP82b as default
    llunar_theory = 3  ! Use ELP-MPP02/DE405 ('historical') as default.  This should be the "default default" :-)
    if(present(lunar_theory)) llunar_theory = lunar_theory
    if(llunar_theory.lt.1 .or. llunar_theory.gt.3) &
         call quit_program_error('planet_position(): lunar_theory must be 1, 2 or 3', 1)
    
    ! Nutation model:
    lnutat = 2000  ! 2000 IAU nutation model
    if(present(nutat)) then
       lnutat = nutat
       if(lnutat.ne.1980 .and. lnutat.ne.2000 .and. lnutat.ne.0) &
            call quit_program_error('planet_position(): nutat must be 2000, 1980 or 0', 1)
    end if
    
    ! Planet-magnitude model:
    lmagmdl = 4  ! AA 2012
    if(present(magmdl)) then
       lmagmdl = magmdl
       if(lmagmdl.lt.1 .or. lmagmdl.gt.4) &
            call quit_program_error('planet_position(): magmdl must be 1-4', 1)
    end if
    
    ! Verbosity:
    lverb = 0  ! Silent
    if(present(verbosity)) lverb = verbosity
    
    
    ! Calc JDE and tjm:
    jdl = jd  ! Local version of JD
    deltat = calc_deltat(jdl)
    ! deltat = 0.d0
    ! write(*,'(/,A,/)') '*** WARNING: libTheSky/planets.f90: setting DeltaT to 0 ***'
    
    if(present(assume_jde) .and. assume_jde) then
       jde = jdl  ! JD provided is actually JDE
       jdl = jde - deltat/86400.d0
    else
       jde = jdl + deltat/86400.d0
    end if
    
    tjm    = (jde-jd2000)/365250.d0                                    ! Julian Millennia after 2000.0 in dyn. time, tau in Meeus
    tjm0   = tjm
    if(lverb.gt.0) then
       print*,'JD:     ', d2s(jdl,8)      ! ms accuracy
       print*,'DeltaT: ', d2s(deltat,3)  ! ms accuracy
       print*,'JDE:    ', d2s(jde,8)     ! ms accuracy
       print*,'t_jm:   ', d2s(tjm,10)    ! ~s accuracy
    end if
    
    call vsop87d_lbr(tjm,3, hcl0,hcb0,hcr0, lLBaccur,lRaccur)             ! Calculate the Earth's true heliocentric l,b,r
    
    if(lverb.gt.3) then
       print*
       ! print*, 'Heliocentric longitude Earth, FK4: ', d2s(modulo(hcl0,pi2)*r2d, 9), ' deg'
       ! print*, 'Heliocentric latitude Earth,  FK4: ', d2s(hcb0*r2d, 9), ' deg'
       ! print*, 'Heliocentric distance Earth:       ', d2s(hcr0, 9), ' au'
       print*
       print*, 'Earth:'
       print*, 't_jm: ', d2s(tjm, 10)
       print*, 'Heliocentric longitude Earth: ', d2s(modulo(hcl0,pi2)*r2d, 9), ' deg'
       print*, 'Heliocentric latitude  Earth: ', d2s(hcb0*r2d, 9), ' deg'
       print*, 'Heliocentric distance  Earth: ', d2s(hcr0, 9), ' au'
    end if
    
    
    ! Iterate to get the light time tau, hence apparent positions:
    tau = 0.d0;  tau1 = 0.d0                                           ! Initial values for light time in days
    j = 0                                                              ! Count iterations; allows escape in case of infinite loop
    do while(deq0(tau1) .or. abs((tau1-tau)/tau1).gt.1.d-10)           ! Moon's tau ~10^-5; 1.d-10~10^-5 sec  (1.d-7~10^-2 sec)
       
       tau = tau1                                                      ! On first iteration, tau=0 to get true positions
       tjm = tjm0 - tau/365250.d0                                      ! Iterate to calculate light time
       jde_lt = jd2000 + tjm*365250                                    ! JDE, corrected for light time
       hcl = 0.d0; hcb = 0.d0; hcr = 0.d0                              ! Heliocentric coordinates
       
       
       ! Compute l,b,r for planet, Pluto, asteroid, Earth's shadow:
       if(pl.gt.0.and.pl.lt.9) call vsop87d_lbr(tjm,pl, hcl,hcb,hcr, lLBaccur,lRaccur)  ! Heliocentric l,b,r
       if(pl.eq.9) call plutolbr(tjm*10.d0, hcl,hcb,hcr)               ! This is for 2000.0, precess 10 lines below
       if(pl.gt.10000) call asteroid_lbr(tjm,pl-10000, hcl,hcb,hcr)    ! Heliocentric lbr for asteroids
       if(pl.eq.-1) then                                               ! Centre of Earth's shadow
          call vsop87d_lbr(tjm,3, hcl,hcb,hcr, lLBaccur,lRaccur)       ! = heliocentric coordinates of the Earth...
          
          select case(llunar_theory)
          case(1)
             call elp82b_lbr(tjm, dumdbl1,dumdbl2,dhcr)                   ! only want dhcr
          case(2)
             call elp_mpp02_lbr(jde_lt, 0, dumdbl1,dumdbl2,dhcr)          ! only want dhcr - 0: LLR mode
          case(3)
             call elp_mpp02_lbr(jde_lt, 1, dumdbl1,dumdbl2,dhcr)          ! only want dhcr - 1: DE405 ('historical') mode
          end select
          
          hcr = hcr + dhcr                                             ! + distance Earth-Moon
       end if
       
       
       ! Compute geocentric ecliptical position for the Moon:
       if(pl.eq.0) then  ! Moon
          select case(llunar_theory)
          case(1)
             call elp82b_lbr(tjm, gcl,gcb,delta)                 ! Get apparent geocentric coordinates of the Moona
          case(2)
             call elp_mpp02_lbr(jde_lt, 0, gcl,gcb,delta)        ! Get apparent geocentric coordinates of the Moon.  Mode=0: LLR mode
          case(3)
             call elp_mpp02_lbr(jde_lt, 1, gcl,gcb,delta)        ! Get apparent geocentric coordinates of the Moon.  Mode=1: DE405 ('historical') mode
          end select
       end if
       
       if(pl.ne.0.and.pl.lt.10 .or. pl.gt.10000) then                  ! Planet, asteroid, Earth's shadow
          ! For Neptune's birthday:
          ! print*,'TESTING!!!'
          ! call precess_ecl(jde,jd2000,hcl,hcb)                        ! from JoD to J2000.0
          
          call hc_spher_2_gc_rect(hcl,hcb,hcr, hcl0,hcb0,hcr0, gcx,gcy,gcz)  ! Convert heliocentric l,b,r to geocentric x,y,z
          call rect_2_spher(gcx,gcy,gcz, gcl,gcb,delta)                      ! Convert geocentric x,y,z to geocentric l,b,r
          if(pl.eq.3) delta = hcr
          if(pl.eq.9) call precess_ecl(jd2000,jde, gcl,gcb)            ! Pluto:  from J2000.0 to JoD
       end if
       
       if(pl.gt.10.and.pl.lt.10000) &
            call cometgc(tjm,tjm0, pl, hcr,gcl,gcb,delta)              ! Calc geocentric l,b,r coordinates for a comet
       
       
       ! Store 'true' position, for tau=0:
       if(j.eq.0) then
          if(abs(tau).gt.1.d-10) call warn('planet_position():  tau != 0 on first light-time iteration', 0)
          
          hcl00  = hcl                                                 ! Heliocentric l,b,r
          hcb00  = hcb
          hcr00  = hcr

          gcx0   = gcx                                                 ! Geocentric x,y,z
          gcy0   = gcy
          gcz0   = gcz
          
          delta0 = delta                                               ! Geocentric l,b,r
          gcl0   = gcl
          gcb0   = gcb
       end if
       
       tau1 = 5.77551830441d-3 * delta                                 ! Light time in days
       
       if(lverb.gt.3) then
          print*
          print*
          print*, 'Iteration:  ', j
          ! print*,'t_jm: ', d2s(tjm, 10)
          ! print*
          ! print*, 'Heliocentric longitude planet, FK4: ', d2s(modulo(hcl00,pi2)*r2d, 9), ' deg'
          ! print*, 'Heliocentric latitude planet,  FK4: ', d2s(hcb00*r2d, 9), ' deg'
          ! print*, 'Heliocentric distance planet:       ', d2s(hcr00, 9), ' au'
          print*
          print*, 'Planet:'
          print*, 't_jm:                   ', d2s(tjm, 10)
          print*, 'Heliocentric longitude: ', d2s(modulo(hcl,pi2)*r2d, 9), ' deg'
          print*, 'Heliocentric latitude:  ', d2s(hcb*r2d, 9), ' deg'
          print*, 'Heliocentric distance:  ', d2s(hcr, 9), ' au'
          print*
          print*, 'Geocentric x: ', d2s(gcx, 9), ' au'
          print*, 'Geocentric y: ', d2s(gcy, 9), ' au'
          print*, 'Geocentric z: ', d2s(gcz, 9), ' au'
          print*
          print*, 'Geocentric longitude: ', d2s(modulo(gcl,pi2)*r2d, 9), ' deg'
          print*, 'Geocentric latitude:  ', d2s(gcb*r2d, 9), ' deg'
          print*, 'True geocentric distance:      ', d2s(delta0, 9), ' au'
          print*, 'Apparent geocentric distance:  ', d2s(delta, 9), ' au'
          print*
          print*, 'JDE_i:       ', d2s(jde_lt, 9)
          print*, 't_jm:        ', d2s(tjm, 10)
          print*, 'Light time:  ', d2s(tau,9),  ' day'
          print*, 'Light time1: ', d2s(tau1,9), ' day'
          print*
       end if
       
       if(.not.lltime) exit                                            ! Do not take into account light time
       if(pl.eq.3) exit                                                ! Want true position for Sun(!)
       
       j = j+1
       if(j.ge.30) then
          call warn('planet_position():  Light time failed to converge for '//trim(enpname(pl)), 0)
          exit
       end if
    end do
    
    tau = tau1
    tjm = tjm0                                                  ! Still Julian Millennia since 2000.0 in dynamical time
    delta = delta0                                              ! Want the TRUE distance, not APPARENT!
    
    if(lverb.gt.3) then
       print*
       print*, 'Number of iterations: ', j
       print*, 'Final light time: ', d2s(tau, 9)
       print*, 'Final t_Jm:       ', d2s(tjm, 10)
       print*
    end if
    
    if(pl.eq.-1) then  ! Earth's shadow
       select case(llunar_theory)
       case(1)
          call elp82b_lbr(tjm, dumdbl1,dumdbl2,delta)           ! Geocentric distance of the Moon for Earth shadow; only need delta
       case(2) 
          call elp_mpp02_lbr(jde_lt, 0, dumdbl1,dumdbl2,delta)  ! Geocentric distance of the Moon for Earth shadow; only need delta - LLR mode
       case(3)
          call elp_mpp02_lbr(jde_lt, 1, dumdbl1,dumdbl2,delta)  ! Geocentric distance of the Moon for Earth shadow; only need delta - DE405 ('historical') mode
       end select
    end if
    
    ! Earth -> Sun.  Note that we want the TRUE position!
    
    ! (As for other planets, we'd use the true position for the Earth (when the light arrives) and the
    ! APPARENT position of the Sun (when the light left).  However, since the positions are heliocentric, the
    ! latter is always (0,0,0).
    if(pl.eq.3) then
       gcl   =  rev(hcl0+pi)
       gcb   = -hcb0
       delta =  hcr0
    end if
    
    
    if(lverb.gt.2) then
       print*
       print*
       ! print*,'t_jm: ', d2s(tjm, 10)
       ! print*
       ! print*, 'Heliocentric longitude planet, FK4: ', d2s(modulo(hcl00,pi2)*r2d, 9), ' deg'
       ! print*, 'Heliocentric latitude planet,  FK4: ', d2s(hcb00*r2d, 9), ' deg'
       ! print*, 'Heliocentric distance planet:       ', d2s(hcr00, 9), ' au'
       print*
       print*, 'Planet after light-time iterations, before correction for nutation, aberration, FK5:'
       print*, 'Final t_jm:                   ', d2s(tjm, 10)
       print*, 'Final heliocentric longitude: ', d2s(modulo(hcl,pi2)*r2d, 9), ' deg'
       print*, 'Final heliocentric latitude:  ', d2s(hcb*r2d, 9), ' deg'
       print*, 'Final heliocentric distance:  ', d2s(hcr, 9), ' au'
       print*
       print*, 'Final geocentric x: ', d2s(gcx, 9), ' au'
       print*, 'Final geocentric y: ', d2s(gcy, 9), ' au'
       print*, 'Final geocentric z: ', d2s(gcz, 9), ' au'
       print*
       print*, 'Final geocentric longitude:        ', d2s(modulo(gcl,pi2)*r2d, 9), ' deg'
       print*, 'Final geocentric latitude:         ', d2s(gcb*r2d, 9), ' deg'
       print*, 'Final (true) geocentric distance:  ', d2s(delta, 9), ' au'
       print*
       print*, 'Final jde_i:       ', d2s(jde_lt, 9)
       print*, 'Final t_jm:        ', d2s(tjm, 10)
       print*, 'Final light time:  ', d2s(tau,9),  ' day'
       print*, 'Final light time1: ', d2s(tau1,9), ' day'
       print*
    end if
    
    ! Compute and correct for nutation:
    if(lnutat.eq.1980) then
       call nutation(tjm, dpsi,eps0,deps)       ! IAU 1980 nutation model: dpsi: nutation in longitude, deps: in obliquity
    else !  if(lnutat.eq.2000) then
       call nutation2000(jdl, dpsi,deps, eps0)  ! IAU 2000 nutation model (originally doesn't provide eps0)
       
       if(lnutat.eq.0) then                     ! No nutation(!)
          dpsi = 0.d0
          deps = 0.d0
       end if
    end if
    eps = eps0 + deps                           ! Correct for nutation in obliquity: mean -> true obliquity of the ecliptic
    
    
    ! Correct for aberration and nutation, and convert to FK5:
    if(pl.lt.10) then  ! I suppose this should also happen for comets, but this compares better...
       if(laber .and. (pl.ne.0)) call aberration_ecl(tjm,hcl0, gcl,gcb)  ! Aberration - not for the Moon
       if(lto_fk5) call fk5(tjm, gcl,gcb)
       gcl = gcl + dpsi  ! Nutation in longitude
    end if
    
    ! sun_gcl,sun_gcb give gcl,gcb as if pl.eq.3
    sun_gcl = rev(hcl0+pi)
    sun_gcb = -hcb0
    
    ! Aberration, FK5 and nutation for sun_gcl,sun_gcb:
    if(laber) call aberration_ecl(tjm,hcl0, sun_gcl,sun_gcb)
    if(lto_fk5)  call fk5(tjm, sun_gcl,sun_gcb)
    sun_gcl = sun_gcl + dpsi  ! Nutation in longitude
    
    ! FK5 conversion for l,b, hcl0,hcb0, hcl00,hcb00:
    if(lto_fk5) then
       call fk5(tjm, hcl,hcb)
       call fk5(tjm, hcl0,hcb0)
       call fk5(tjm, hcl00,hcb00)
    end if
    
    if(lverb.gt.3) then
       print*
       print*, 'Convert to the FK5 system:        ', lto_fk5
       print*, 'True heliocentric longitude, FK5: ', d2s(modulo(hcl00,pi2)*r2d, 9), ' deg'
       print*, 'True heliocentric latitude,  FK5: ', d2s(hcb00*r2d, 9), ' deg'
       print*, 'True heliocentric distance:       ', d2s(hcr00, 9), ' au'
       print*
    end if
    
    
    ! Convert geocentric ecliptical to equatorial coordinates:
    call ecl_2_eq(gcl,gcb,eps, ra,dec)    ! RA, Dec
    
    
    ! Sidereal time:
    gmst = calc_gmst(jdl)                 ! Greenwich mean sidereal time
    agst = rev(gmst + dpsi*cos(eps))      ! Correction for equation of the equinoxes -> Greenwich apparent sidereal time
    lst  = rev(agst + llon)               ! Local apparent sidereal time, llon > 0 for E
    
    
    ! Apparent diameter:
    if(pl.ge.0.and.pl.lt.10) diam = atan(pland(pl)/(delta*au))
    !if(pl.lt.10.and.pl.ge.0) diam = atan(planr(pl)/(2*delta*au))*2
    if(pl.eq.-1) then 
       call earthshadow(delta,hcr0, rES1,rES2)  ! Calc Earth shadow radii at Moon distance:  (pen)umbra radius - geocentric
       diam = 2*rES1                            ! Umbra diameter to planpos(12)
    end if
    
    ! Convert heliocentric to topocentric coordinates:
    call geoc2topoc_ecl(gcl,gcb, delta,diam/2.d0, eps,lst, topl,topb,topdiam, lat=llat,hgt=lhgt)  ! Geocentric to topoc: l, b, diam
      if(pl.eq.-1) topdiam = rES2                                                            ! Earth penumbra radius at Moon distance
      topdiam = 2*topdiam                                                                    ! Was radius, now diameter
    if(pl.ge.0.and.pl.lt.10) topdelta = pland(pl)/tan(topdiam)/au
    
    
    ! Convert equatorial to horizontal coordinates:
    call eq2horiz(ra,dec,agst, hh,az,alt, lat=llat,lon=llon)
    
    
    ! Elongation:
    elon = acos(cos(gcb)*cos(hcb0)*cos(gcl-rev(hcl0+pi)))
    if(pl.gt.10) elon = acos((hcr0**2 + delta**2 - hcr**2)/(2*hcr0*delta))        ! For comets (only?) - CHECK: looks like phase angle...
    
    ! Convert topocentric coordinates: ecliptical -> equatorial -> horizontal:
    ! call geoc2topoc_eq(ra,dec,delta,hh,topra,topdec)                            ! Geocentric to topocentric 
    call ecl_2_eq(topl,topb,eps, topra,topdec)                                   ! Topocentric l,b -> RA,dec - identical result (AA 1999-01-02/22); probably cheaper
    call eq2horiz(topra,topdec,agst, tophh,topaz,topalt, lat=llat,lon=llon)
    
    ! Phase angle:
    pa = 0.d0  ! Make sure pa is always defined
    if(pl.eq.0) then
       pa = atan2( hcr0*sin(elon) , delta - hcr0*cos(elon) )         ! Moon
    else if(pl.gt.0) then
       pa = acos( (hcr**2 + delta**2 - hcr0**2) / (2*hcr*delta) )    ! Planets
    end if
    illfr = 0.5d0*(1.d0 + cos(pa))                                   ! Illuminated fraction of the disc
    
    
    ! Compute magnitude:
    if(pl.gt.0.and.pl.lt.10) then
       magn = planet_magnitude(pl,hcr,delta,pa, lmagmdl)
       if(pl.eq.6) magn = magn + dsatmagn(tjm*10., gcl,gcb)  ! Correct Saturn's magnitude for rings
    end if
    
    if(pl.eq.0) magn = moonmagn(pa,delta)                             ! Moon
    if(pl.gt.10000) magn = asteroid_magn(pl-10000,delta,hcr,pa)       ! Asteroid magnitude (valid for |pa|<120°)
    
    ! Comet magnitude:
    if(pl.gt.10 .and. pl.lt.10000) then
       diam = 0.d0
       if(cometDiedAtP(pl) .and. jdl.gt.cometElems(pl,7)) then        ! Comet died at perihelion and is now dead
          magn = 99.9d0
       else
          magn = cometElems(pl,8) + 5*log10(delta) + 2.5d0*cometElems(pl,9)*log10(hcr)  ! m = H + 5log(d) + 2.5*G*log(r)
          magn = magn + comet_scatter_magnitude_correction(pa)
       end if
       topdelta = delta
    end if
    
    
    ! Parallactic angle:
    !parang     = atan2(sin(hh),tan(llat)*cos(dec) - sin(dec)*cos(hh))                      ! Parallactic angle (between zenith and celestial north pole): geocentric - Meeus Eq.14.1
    parang     = atan2(sin(tophh),tan(llat)*cos(topdec) - sin(topdec)*cos(tophh))           ! Parallactic angle (between zenith and celestial north pole): topocentric
    parang_ecl = atan( (cos(topl)*tan(eps)) / (sin(topl)*sin(topb)*tan(eps) - cos(topb)) )  ! "Parallactic angle" (between ECLIPTICAL and celestial north pole): topocentric - Meeus p.100
    
    ! Horizontal parallax:
    hp = asin(earthr/(delta*au))                                                ! Horizontal parallax
    
    
    ! Overrule some values: 
    if(pl.eq.0) hcr    = 0.d0                                                   ! No heliocentric distance  for the Moon
    if(pl.eq.3) then  ! Sun:
       elon  = 0.d0                                                             ! Elongation
       pa    = 0.d0                                                             ! Phase angle
       magn  = sunmagn(delta)                                                   ! Magnitude; changes ~0.05 over half a year...
       illfr = 1.d0                                                             ! Illuminated fraction
    end if
    
    
    
    
    ! Save variables:
    ! Geocentric:
    planpos     = 0.d0              ! 1-12 are geocentric, repeated as topocentric in 21-32
    
    planpos(1)  = rev(gcl)          ! Ecliptical longitude
    planpos(2)  = rev2(gcb)         ! Ecliptical latitude
    planpos(3)  = hcr               ! Distance to the Sun (heliocentric!)
    
    planpos(4)  = delta             ! True(!) geocentric distance
    planpos(5)  = rev(ra)           ! R.A.
    planpos(6)  = dec               ! Declination
    
    planpos(7)  = tau               ! Light time in days
    
    planpos(8)  = rev(hh)           ! Hour angle
    planpos(9)  = rev(az)           ! Azimuth
    planpos(10) = alt               ! Altitude
    
    planpos(11) = rev(elon)         ! Elongation
    planpos(12) = diam              ! Apparent diameter
    
    if(lverb.gt.1) then
       print*
       print*, 'Final apparent geocentric longitude:  ', d2s(rev(gcl)*r2d, 9), ' deg'
       print*, 'Final apparent geocentric latitude:   ', d2s(rev2(gcb)*r2d, 9), ' deg'
       print*, 'Final apparent geocentric distance:   ', d2s(delta, 9), ' au'
       print*, 'Final apparent heliocentric distance: ', d2s(hcr, 9), ' au'
       print*
       print*, 'Final apparent geocentric right ascension:  ', d2s(rev(ra)*r2h, 9), ' hr'
       print*, 'Final apparent geocentric declination:      ', d2s(dec*r2d, 9), ' deg'
       print*
    end if

    
    planpos(13) = magn              ! Apparent visual magnitude
    planpos(14) = illfr             ! Illuminated fraction
    planpos(15) = rev(pa)           ! Phase angle
    
    planpos(16) = parang            ! Topocentric parallactic angle (between celestial pole and zenith)
    planpos(17) = hp                ! Horizontal parallax
    planpos(18) = parang_ecl        ! Topocentric ecliptical "parallactic angle" (between celestial and ecliptical pole)
    
    ! Topocentric:
    planpos(21) = rev(topl)   
    planpos(22) = rev2(topb)
    planpos(23) = delta0            ! True geocentric distance (== delta?)
    
    planpos(24) = topdelta
    planpos(25) = rev(topra)
    planpos(26) = topdec
    
    planpos(28) = rev(tophh)        ! Hour angle
    planpos(29) = rev(topaz)
    planpos(30) = topalt
    planpos(31) = topalt + refract(topalt)  ! Topocentric altitude, corrected for (cheap) atmospheric refraction
    
    planpos(32) = topdiam
    
    ! Heliocentric:
    planpos(33) = rev(hcl)          ! Heliocentric apparent l,b,r
    planpos(34) = rev2(hcb) 
    planpos(35) = hcr
    
    planpos(36) = rev(hcl00)        ! Heliocentric true l,b,r, converted to FK5
    planpos(37) = rev2(hcb00)
    planpos(38) = hcr00
    
    ! Other variables:
    planpos(39) = dble(pl)          ! Remember for which planet data were computed
    planpos(40) = jde               ! JDE for these data
    
    planpos(41) = rev(hcl0+pi+dpsi) ! Geocentric, true L,B,R for the Sun, in FK5, corrected for nutation
    planpos(42) = rev2(-hcb0)
    planpos(43) = hcr0
    
    planpos(44) = lst               ! Local APPARENT sidereal time
    planpos(45) = rev(agst)         ! Greenwich APPARENT sidereal time (in radians)
    planpos(46) = tjm0 * 10         ! Apparent dynamical time in Julian Centuries since 2000.0
    
    planpos(47) = dpsi              ! Nutation in longitude
    planpos(48) = eps               ! True obliquity of the ecliptic; corrected for nutation
    planpos(49) = rev(gmst)         ! Greenwich MEAN sidereal time (in radians)
    planpos(50) = eps0              ! Mean obliquity of the ecliptic; without nutation
    
    ! Geocentric:
    planpos(51) = rev(sun_gcl)      ! Apparent geocentric longitude of the Sun (variable was treated as if pl.eq.3)
    planpos(52) = rev2(sun_gcb)     ! Apparent geocentric latitude of the Sun (variable was treated as if pl.eq.3)

    planpos(58) = jdl               ! JD used for calculation
    planpos(59) = deltat            ! DeltaT used for calculation
    
    planpos(61) = gcx               ! Apparent geocentric x,y,z
    planpos(62) = gcy
    planpos(63) = gcz
    
    planpos(64) = gcx0              ! True geocentric x,y,z
    planpos(65) = gcy0
    planpos(66) = gcz0
    
    planpos(67) = gcl0              ! True geocentric l,b,r
    planpos(68) = gcb0
    planpos(69) = delta0            ! (==delta?)
    
  end subroutine planet_position
  !*********************************************************************************************************************************
  
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate Pluto's position l,b,r at time t
  !! 
  !! \param  t   Dynamical time in Julian Centuries after 2000.0
  !! \param l   Heliocentric longitude (rad) (output)
  !! \param b   Heliocentric latitude (rad) (output)
  !! \param r   Heliocentric distance (AU?) (output)
  
  subroutine plutolbr(t,l,b,r)
    use SUFR_kinds, only: double
    use SUFR_constants, only: d2r
    use TheSky_planetdata, only: pluc,plul,plub,plur
    
    implicit none
    real(double), intent(in) :: t
    real(double), intent(out) :: l,b,r
    integer :: i
    real(double) :: j,s,p,a, cosa,sina
    
    j =  34.35d0 + 3034.9057d0 * t
    s =  50.08d0 + 1222.1138   * t
    p = 238.96d0 + 144.96d0    * t
    
    l = 0.d0
    b = 0.d0
    r = 0.d0
    
    do i=1,43
       a = (pluc(i,1)*j + pluc(i,2)*s + pluc(i,3)*p) * d2r
       cosa = cos(a)
       sina = sin(a)
       l = l + plul(i,1)*sina + plul(i,2)*cosa
       b = b + plub(i,1)*sina + plub(i,2)*cosa
       r = r + plur(i,1)*sina + plur(i,2)*cosa
    end do
    
    l = ( l*1.d-6 + 238.958116d0  + 144.96d0*t ) * d2r
    b = ( b*1.d-6 -   3.908239d0               ) * d2r
    r =   r*1.d-7 +  40.7241346d0
    
  end subroutine plutolbr
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate orbital elements for the planets
  !!
  !! \param jd  Julian day of calculation
  !!
  !! \note
  !! Results are returned through the module planetdata, for planet pl and element el:
  !! - plelems(pl,el):     EoD
  !! - plelems2000(pl,el): J2000.0
  !!
  !! \note
  !! Elements (el):
  !! - 1: L     - mean longitude
  !! - 2: a     - semi-major axis
  !! - 3: e     - eccentricity
  !! - 4: i     - inclination
  !! - 5: Omega - longitude of ascending node
  !! - 6: pi    - longitude of perihelion
  
  subroutine planetelements(jd)
    use SUFR_kinds, only: double
    use SUFR_constants, only: d2r, jd2000
    use SUFR_system, only: quit_program_error
    use SUFR_angles, only: rev
    use SUFR_numerics, only: deq0
    
    use TheSky_local, only: deltat
    use TheSky_planetdata, only: plelems, plelems2000, plelemdata
    use TheSky_datetime, only: calc_deltat
    
    implicit none
    real(double), intent(in) :: jd
    integer :: el,pl,te
    real(double) :: jde,tm,tms(0:3)
    
    if(deq0(plelemdata(1,1,1,1)))  call quit_program_error('planetelements():  did you forget to call readplanetelements()?', 0)
    
    deltat = calc_deltat(jd)
    jde = jd + deltat/86400.d0
    tm = (jde-jd2000)/36525.d0  ! Julian Centuries after 2000.0 in dynamical time
    
    
    ! Powers of tm:
    do te=0,3
       tms(te) = tm**te
    end do
    
    
    do pl=1,8
       do el=1,6
          
          plelems(pl,el) = 0.d0
          plelems2000(pl,el) = 0.d0
          
          do te=0,3
             plelems(pl,el)     = plelems(pl,el)     + plelemdata(1,pl,el,te) * tms(te)
             plelems2000(pl,el) = plelems2000(pl,el) + plelemdata(2,pl,el,te) * tms(te)
          end do
          
          if(el.ne.2.and.el.ne.3) then
             plelems(pl,el)     = rev(plelems(pl,el)     * d2r)
             plelems2000(pl,el) = rev(plelems2000(pl,el) * d2r)
          end if
          
       end do
    end do
    
  end subroutine planetelements
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate planet magnitude
  !! 
  !! \param pl          Planet ID
  !! \param hc_dist     Distance from the Sun (AU)
  !! \param gc_dist     Distance from the Earth (AU)
  !! \param phang       Phase angle (rad)
  !! \param model       Model to use: 1: Müller (1893), 2: Meeus p.286, 3: AA 1992, 4: AA 2012 (optional, defaults to 4).
  !! 
  !! \retval planet_magnitude  Apparent visual planet magnitiude
  !! 
  !! \see
  !! - Müller, POPot 8, 193 (1893) - https://ui.adsabs.harvard.edu/abs/1893POPot...8..193M/ - p.341/149 - Model 1.
  !! - Meeus, Astronomical Algorithms, 1998, Ch.41 - Model 2 (origin unsure).
  !! - Expl.Supl.tt.Astr.Almanac, 2nd Ed. (1992) - Model 3.
  !! - Expl.Supl.tt.Astr.Almanac, 3rd Ed. (2012) - Model 4 (default).
  
  function planet_magnitude(pl, hc_dist,gc_dist, phang, model)
    use SUFR_kinds, only: double
    use SUFR_constants, only: r2d
    use SUFR_system, only: quit_program_error
    
    implicit none
    integer, intent(in) :: pl
    integer, intent(in), optional :: model
    integer :: lpl, lmodel
    real(double), intent(in) :: hc_dist,gc_dist,phang
    real(double) :: planet_magnitude,pa,pa2,pa3,a0(9),a1(9),a2(9),a3(9)
    
    ! Optional argument:
    lmodel = 4  ! Model 4: Expl.Supl.tt.Astr.Almanac, 3rd Ed. (2012)
    if(present(model)) lmodel = model
    
    ! PA must be in degrees:
    pa  = phang*r2d
    if(lmodel.eq.1 .and. pl.eq.1) pa = pa - 50.d0  ! For Mercury (deg) -- ONLY FOR MODEL 1 !!!
    pa2 = pa*pa
    pa3 = pa*pa2
    
    lpl = pl  ! Local pl
    
    if(lmodel.eq.1) then       ! Model 1: Müller (1893)
       !       Mer      Ven       Ear   Mars     Jup     Sat     Ur      Nep     Pl
       a0 = [-1.16d0,  4.d0,     0.d0, 1.3d0,   8.93d0, 8.68d0, 6.85d0, 7.05d0, 1.d0] * (-1)
       a1 = [ 2.838d0, 1.322d0,  0.d0, 1.486d0, 0.d0,   0.d0,   0.d0,   0.d0,   0.d0] * 1.d-2
       a2 = [ 1.023d0, 0.d0,     0.d0, 0.d0,    0.d0,   0.d0,   0.d0,   0.d0,   0.d0] * 1.d-4
       a3 = [ 0.d0,    0.4247d0, 0.d0, 0.d0,    0.d0,   0.d0,   0.d0,   0.d0,   0.d0] * 1.d-6
       
    else if(lmodel.eq.2) then  ! Model 2: Meeus p.286
       !        Mer     Ven     Ear   Mars    Jup    Sat     Ur      Nep     Pl
       a0 = [ 0.42d0, 4.40d0, 0.d0, 1.52d0, 9.40d0, 8.88d0, 7.19d0, 6.87d0, 1.00d0] * (-1)
       a1 = [ 3.80d0, 0.09d0, 0.d0, 1.6d0,  0.5d0,  0.d0,   0.d0,   0.d0,   0.d0]   * 1.d-2
       a2 = [-2.73d0, 2.39d0, 0.d0, 0.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0]   * 1.d-4
       a3 = [ 2.d0,  -0.65d0, 0.d0, 0.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0]   * 1.d-6
       
    else if(lmodel.eq.3) then  ! Model 3: Expl.Supl.tt.Astr.Almanac, 2nd Ed. (1992)
       !        Mer     Ven     Ear   Mars    Jup    Sat     Ur      Nep     Pl
       a0 = [ 0.36d0, 4.29d0, 0.d0, 1.52d0, 9.25d0, 8.88d0, 7.19d0, 6.87d0, 1.00d0] * (-1)
       a1 = [ 3.80d0, 0.09d0, 0.d0, 1.6d0,  0.5d0,  0.d0,   0.d0,   0.d0,   0.d0]   * 1.d-2
       a2 = [-2.73d0, 2.39d0, 0.d0, 0.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0]   * 1.d-4
       a3 = [ 2.d0,  -0.65d0, 0.d0, 0.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0]   * 1.d-6
       
    else if(lmodel.eq.4) then  ! Model 4: Expl.Supl.tt.Astr.Almanac, 3rd Ed. (2012)
       !       Mer      Ven       Ven2    Mars     Jup     Sat     Ur      Nep     Pl
       a0 = [ 0.60d0,  4.47d0,  -0.98d0, 1.52d0,  9.40d0, 8.88d0, 7.19d0, 6.87d0, 1.01d0] * (-1)
       a1 = [ 4.98d0,  1.03d0,  -1.02d0, 1.6d0,   0.5d0,  4.4d0,  0.2d0,  0.d0,   0.d0] * 1.d-2
       a2 = [-4.88d0,  0.57d0,   0.d0,   0.d0,    0.d0,   0.d0,   0.d0,   0.d0,   0.d0] * 1.d-4
       a3 = [ 3.02d0,  0.13d0,   0.d0,   0.d0,    0.d0,   0.d0,   0.d0,   0.d0,   0.d0] * 1.d-6
       
       if(pl.eq.2 .and. pa.gt.163.6d0) lpl = 3  ! Venus 2
    else
       call quit_program_error('planet_magnitude():  model must be 1-4', 0)
    end if
    
    planet_magnitude = 5*log10(hc_dist*gc_dist) + a0(lpl) + a1(lpl)*pa + a2(lpl)*pa2 + a3(lpl)*pa3
    
  end function planet_magnitude
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate Saturn's magnitude
  !!
  !! \param t   Dynamical time in Julian centuries after J2000.0
  !!
  !! \param gl  Geocentric longitude (rad)
  !! \param gb  Geocentric latitude (rad)
  !! \param d   Geocentric distance (AU)
  !!
  !! \param l   Heliocentric longitude (rad)
  !! \param b   Heliocentric latitude (rad)
  !! \param r   Heliocentric distance (AU)
  !! 
  !! \param magmdl  Model to use:  1: Müller (1893),  2: Meeus p.286;  optional, default: 2.
  !!
  !! \see
  !! - Müller, POPot 8, 193 (1893) - https://ui.adsabs.harvard.edu/abs/1893POPot...8..193M/ - p.341/149
  !! - Meeus, Astronomical Algorithms, 1998, Ch.41
  
  function satmagn(t, gl,gb,d, l,b,r, magmdl)
    use SUFR_kinds, only: double
    use SUFR_constants, only: d2r,r2d
    use SUFR_angles, only: rev2
    
    implicit none
    real(double), intent(in) :: t, gl,gb,d, l,b,r
    integer, intent(in), optional :: magmdl
    real(double) :: satmagn, in,om,bbb,n,ll,bb,u1,u2,du
    integer :: lmagmdl
    
    lmagmdl = 2
    if(present(magmdl)) lmagmdl = magmdl
    
    in = 0.490005d0  - 2.2686d-4*t + 7.d-8*t**2      ! Radians
    om = 2.9584809d0 + 0.0243418d0*t + 7.19d-6*t**2  ! Radians
    
    bbb = asin( sin(in)*cos(gb)*sin(gl-om) - cos(in)*sin(gb) )
    n  = (113.6655d0 + 0.8771d0*t)*d2r
    ll = l - 0.01759d0*d2r/r
    bb = b - 7.64d-4*d2r * cos(l-n)/r
    
    u1 = atan2( sin(in)*sin(bb) + cos(in)*cos(bb)*sin(ll-om) , cos(bb)*cos(ll-om) )
    u2 = atan2( sin(in)*sin(gb) + cos(in)*cos(gb)*sin(gl-om) , cos(gb)*cos(gl-om) )
    du = abs(rev2(u1-u2))  ! rev2() is needed because u1 and u2 jump from pi to -pi at different moments
    
    if(lmagmdl.eq.1) then
       ! Model 1 - Müller (1983):
       satmagn = -8.68d0 + 5.d0*log10(d*r) + 0.044d0*abs(du)*r2d - 2.60d0*sin(abs(bbb)) + 1.25d0*(sin(bbb))**2
    else
       ! Method 2 - Meeus p.286:  - constant difference of 0.2m...(?)
       satmagn = -8.88d0 + 5.d0*log10(d*r) + 0.044d0*abs(du)*r2d - 2.60d0*sin(abs(bbb)) + 1.25d0*(sin(bbb))**2
    end if
    
  end function satmagn
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate a correction to Saturn's magnitude due to its rings.
  !! 
  !! \param tjc    Dynamical time in Julian centuries after J2000.0.
  !! 
  !! \param glon   Geocentric longitude (rad).
  !! \param glat   Geocentric latitude (rad).
  !! 
  !! \see Müller, POPot 8, 193 (1893) - https://ui.adsabs.harvard.edu/abs/1893POPot...8..193M/ - p.341/149
  
  function dsatmagn(tjc, glon,glat)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev2
    
    implicit none
    real(double), intent(in) :: tjc, glon,glat
    real(double) :: dsatmagn, inc,omg, bbb  ! , n,ll,bb, u1,u2,du, bbb
    
    inc = 0.490005d0  - 2.2686d-4*tjc   + 7.d-8*tjc**2    ! Radians
    omg = 2.9584809d0 + 0.0243418d0*tjc + 7.19d-6*tjc**2  ! Radians
    
    bbb = asin( sin(inc)*cos(glat)*sin(glon-omg) - cos(inc)*sin(glat) )
    
    dsatmagn = -2.60d0*sin(abs(bbb)) + 1.25d0*(sin(bbb))**2
    
  end function dsatmagn
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate the umbra and penumbra geocentric radii of the Earth's shadow
  !! 
  !! \param dm0  Distance of the Moon (AU)
  !! \param ds0  Distance of the Sun  (AU)
  !!
  !! \param r1  Umbra radius at distance of the Moon (rad) (output)
  !! \param r2  Penumbra radius at distance of the Moon (rad) (output)
  !!
  !! \see  Expl.sup. to the Astronomical Almanac, p.428
  
  subroutine earthshadow(dm0,ds0, r1,r2)
    use SUFR_kinds, only: double
    use SUFR_constants, only: au, earthr,planr
    
    implicit none
    real(double), intent(in) :: dm0,ds0
    real(double), intent(out) :: r1,r2
    real(double) :: rs,re,ds,dm,p1,ps,ss
    
    dm = dm0*au  ! AU -> cm
    ds = ds0*au  ! AU -> cm
    
    re = earthr
    rs = planr(3)
    
    p1 = asin(re/dm) * 0.998340d0        ! Moon parallax, compensated for Earth's atmosphere
    ps = asin(re/ds)                     ! Solar parallax
    ss = asin(rs/ds)                     ! Solar semi-diameter
    
    r1 = 1.02d0 * (p1 + ps - ss)
    r2 = 1.02d0 * (p1 + ps + ss)
    
  end subroutine earthshadow
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute physical data for Jupiter
  !!
  !! \param  jd      Julian day for computation
  !!
  !! \param de      Jovocentric latitude of the Earth = inclination of planet axis to plane of the sky as seen from Earth (output)
  !! \param ds      Jovocentric latitude of the Sun   = inclination of planet axis to plane of the sky as seen from Earth (output)
  !! \param omg1    Longitude of central meridian of System I (equator+-10deg) (output)
  !! \param omg2    Longitude of central meridian of System II (output)
  !!
  !! \param dphase  Phase correction (output)
  !! \param pa      Position angle of Jupiter's north pole (from N to E) (output)
  !! \param in      Inclination of Jupiter's rotation axis to the orbital plane (output)
  !! \param om      Longitude of node of Jupiter's equator on ecliptic (output)
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.43, p.293-295
  
  subroutine jupiterphys(jd,  de,ds, omg1,omg2, dphase, pa, in,om)
    use SUFR_kinds, only: double
    use SUFR_constants, only: d2r, pi, jd1900
    use SUFR_angles, only: rev
    
    use TheSky_local, only: deltat
    use TheSky_planetdata, only: planpos
    use TheSky_coordinates, only: ecl_2_eq, eq_2_ecl
    use TheSky_datetime, only: calc_deltat
    
    implicit none
    real(double), intent(in) :: jd
    real(double), intent(out) :: de,ds, omg1,omg2, dphase, pa, in,om
    real(double) :: jde,d,t1,t2,t3
    real(double) :: alp0,del0,w1,w2,x,y,z,eps0,eps,dpsi,l,b,alps,dels, r,r0,l0
    real(double) :: u,v,alp,del,dze,delta,l1,b1
    
    ! To get data for the exact (rounded off) moment of Meeus' example:
    !d = 15690.00068d0
    !jde = d + 2433282.5d0
    !jd = jde - deltat/86400.d0
    
    call planet_position(jd,5)
    call planetelements(jd) 
    
    deltat = calc_deltat(jd)
    jde    = jd + deltat/86400.d0
    
    ! Meeus, step 1, p.293:
    d      = jde - 2433282.5d0  ! Since 1950 (!)
    t1     = d/36525.d0  ! Time since 1950 (!) in Julian centuries
    
    t2     = (jde - jd1900)/36525.d0           ! Time since J1900.0 (!), in Julian centuries
    in = (3.120262d0 + 0.0006d0 * t2) * d2r    ! Inclination of Jupiter's rotation axis to the orbital plane (Meeus, p.311)
    
    t3 = jde - 2443000.5d0 - planpos(7)        ! Days since 1976-08-10 (Meeus, p.305), minus light-travel time
    om = (316.518203d0 - 2.08362d-6*t3) * d2r  ! Long. of node of Jupiter's equator on ecliptic, psi in Meeus, p.306
    om = om + (1.3966626d0*t1 + 3.088d-4*t1**2)*d2r  ! Correct for precession (since B1950, but J1950 will do), Meeus p.311
    
    alp0 = rev((268.d0 + 0.1061d0*t1)*d2r)  ! Right ascension of Jupiter's north pole
    del0 = rev((64.5d0 - 0.0164d0*t1)*d2r)  ! Declination of Jupiter's north pole
    
    ! Meeus, step 2, p.294:
    w1 = rev((17.710d0 + 877.90003539d0*d)*d2r)  ! Longitude system I
    w2 = rev((16.838d0 + 870.27003539d0*d)*d2r)  ! Longitude system II
    
    x = planpos(61)  ! Apparent geocentric x
    y = planpos(62)  ! Apparent geocentric y
    z = planpos(63)  ! Apparent geocentric z
    
    l = planpos(33)  ! Heliocentric apparent l
    b = planpos(34)  ! Heliocentric apparent b
    r = planpos(35)  ! Heliocentric apparent r
    delta = planpos(4)  ! Apparent geocentric distance
    
    l0 = rev(planpos(41)+pi)  ! Heliocentric, true coordinates of the Earth
    r0 = planpos(43)
    
    ! Meeus, step 8, p.294:
    eps0 = planpos(50)  ! Mean obliquity of the ecliptic; without nutation
    eps  = planpos(48)  ! True obliquity of the ecliptic; corrected for nutation
    dpsi = planpos(47)  ! Nutation in longitude
    
    ! Meeus, step 9, p.294:
    call ecl_2_eq(l,b,eps, alps,dels)
    
    ! Meeus, step 10, p.294:
    ds = -sin(del0)*sin(dels) - cos(del0)*cos(dels)*cos(alp0-alps)  ! Planetocentric declination of the Sun
    
    ! Meeus, step 11, p.294:
    u = y*cos(eps0) - z*sin(eps0)
    v = y*sin(eps0) + z*cos(eps0)
    
    alp = rev(atan2(u,x))
    del = atan2(v,sqrt(x**2+u**2))
    dze = rev( atan2( sin(del0)*cos(del)*cos(alp0-alp) - sin(del)*cos(del0) , cos(del)*sin(alp0-alp) ) )
    
    ! Meeus, step 12, p.294:
    de = -sin(del0)*sin(del) - cos(del0)*cos(del)*cos(alp0-alp)  ! Planetocentric declination of the Earth
    
    ! Meeus, step 14, p.295 - phase correction - same sign as sin(l-l0):
    dphase = abs( (2*r*delta + r0**2 - r**2 - delta**2)/(4*r*delta) ) * sin(l-l0)/abs(sin(l-l0))
    
    ! Meeus, step 13, p.295:
    omg1 = rev(w1 - dze - 5.07033d0*d2r * delta + dphase)  ! Longitude of central meridian of System I (equator +- 10deg)
    omg2 = rev(w2 - dze - 5.02626d0*d2r * delta + dphase)  ! Longitude of central meridian of System II (higher lat)
    
    ! Meeus, steps 16, 17, p.295:
    call eq_2_ecl(alp0,del0,eps0, l1,b1)      ! Has to be eps0
    call ecl_2_eq(l1+dpsi,b1,eps, alp0,del0)  ! Has to be eps
    
    ! Meeus, step 18, p.295:
    pa = atan2( cos(del0)*sin(alp0-alp) , sin(del0)*cos(del) - cos(del0)*sin(del)*cos(alp0-alp) )  ! PA of Jupiter's north pole
    
  end subroutine jupiterphys
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute physical data for Saturn
  !!
  !! \param  jd    Julian day of computation
  !!
  !! \param be    Saturnicentric latitude of the Earth = inclination of planet axis to plane of the sky as seen from Earth (output)
  !! \param bs    Saturnicentric latitude of the Sun   = inclination of planet axis to plane of the sky as seen from the Sun (output)
  !!
  !! \param pa    Position angle of Saturn's north pole (from N to E) (output)
  !! \param in    Inclination of Saturn's rotation axis to the orbital plane (output)
  !! \param om    Longitude of node of Saturn's equator on ecliptic (output)
  !!
  !! \param ar    Projected major axis of ring (NOT semi-!!!) (output)
  !! \param br    Projected minor axes of ring (NOT semi-!!!), = ar * sin(be) (output)
  !! \param du    Difference between saturnicentric longitudes of Sun and Earth (Sun - Earth, NOT abs!!!) (output)
  !! 
  !! \param pa_s  Position angle of Saturn's north pole (from N to E), as seen from the SUN (output)
  !! \param ar_s  Projected major axis of ring (NOT semi-!!!), as seen from the SUN (output)
  !! \param br_s  Projected minor axes of ring (NOT semi-!!!), as seen from the SUN, = ar_s * sin(bs) (output)
  !!
  !!
  !! \see Meeus, Astronomical Algorithms, 1998, Ch.45  (The ring of Saturn)
  !!
  !! \todo
  !! - check true/apparent coordinates (see CHECK)
  
  subroutine saturnphys(jd, be,bs, pa, in,om, ar,br, du,  pa_s, ar_s,br_s)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi, as2r,d2r
    use SUFR_angles, only: rev
    
    use TheSky_planetdata, only: planpos, plelems
    use TheSky_coordinates, only: ecl_2_eq, hc_spher_2_gc_rect, rect_2_spher
    
    implicit none
    real(double), intent(in) :: jd
    real(double), intent(out) :: be,bs, pa, in,om, ar,br, du,  pa_s, ar_s,br_s
    integer :: loop
    real(double) :: t, lam,bet, d,l,b,r, l0,b0,r0, x,y,z, dpsi,eps, ascn
    real(double) :: l2,b2, u1,u2, lam0,bet0, dl,db, ra,dec, ra0,dec0,  sinbe
    
    call planet_position(jd,6)  ! Position of Saturn
    
    ! Meeus, step 3:
    l   =  planpos(33)  ! Heliocentric longitude of Saturn
    b   =  planpos(34)  ! Heliocentric latitude of Saturn
    r   =  planpos(35)  ! Heliocentric distance of Saturn
    
    t    = planpos(46)  ! App. dyn. time in Julian Centuries since 2000.0
    
    ! Meeus, step 10:
    dpsi = planpos(47)  ! Nutation in longitude
    eps  = planpos(48)  ! True obliquity of the ecliptic; corrected for nutation
    
    do loop=1,2  ! From Sun, Earth:
       
       ! Meeus, step 2:
       if(loop.eq.1) then  ! Seen from the Sun:
          lam = l
          bet = b
          d = r
       else  ! Seen from the Earth:
          l0  =  rev(planpos(41)+pi)  ! Geocentric, true L,B,R for the Earth, in FK5  -  CHECK - need apparent?
          b0  = -planpos(42)
          r0  =  planpos(43)
          
          ! Meeus, Ch. 45, step 5:
          call hc_spher_2_gc_rect(l,b,r, l0,b0,r0, x,y,z)
          call rect_2_spher(x,y,z, lam,bet,d)
       end if
       
       ! Meeus, step 1:
       in = (28.075216d0  - 0.012998d0*t + 4.d-6 * t**2)   * d2r  ! Inclination of rotation axis,               Eq. 45.1a
       om = (169.508470d0 + 1.394681d0*t + 4.12d-4 * t**2) * d2r  ! Longitude of asc. node of equator - Omega,  Eq. 45.1b
       
       ! Meeus, step 6:
       sinbe = sin(in)*cos(bet)*sin(lam-om) - cos(in)*sin(bet)
       ar    = 375.35d0*as2r/d  ! Major axis of outer edge of outer ring
       br    = ar*sinbe         ! Minor axis of outer edge of outer ring
       be    = asin(sinbe)      ! B
       
       if(loop.eq.2) then  ! Seen from Earth:
          ! Meeus, step 7:
          call planetelements(jd)
          ascn = plelems(6,5)                     ! Longitude of Saturn's ascending node
          l2 = l - 0.01759d0*d2r / r              ! Correct for the Sun's aberration on Saturn - lon
          b2 = b - 7.64d-4*d2r * cos(l-ascn) / r  ! Correct for the Sun's aberration on Saturn - lat
          
          ! Meeus, step 8:
          bs = asin(sin(in)*cos(b2)*sin(l2-om) - cos(in)*sin(b2))  ! B'
          
          ! Meeus, step 9:
          u1 = atan2(sin(in)*sin(b2)  + cos(in)*cos(b2)*sin(l2-om),   cos(b2)*cos(l2-om))    ! Saturnicentric longitude of the Sun
          u2 = atan2(sin(in)*sin(bet) + cos(in)*cos(bet)*sin(lam-om), cos(bet)*cos(lam-om))  ! Saturnicentric longitude of the Earth
          du = u1-u2
       end if
       
       ! Meeus, step 11:
       lam0 = om - pi/2.d0  ! Ecliptical longitude of Saturn's north pole
       bet0 = pi/2.d0 - in  ! Ecliptical latitude of Saturn's north pole
       
       if(loop.eq.2) then  ! Seen from the Earth:
          ! Meeus, step 12 - correct for the aberration of Saturn:
          dl  = 0.005693d0*d2r * cos(l0-lam)/cos(bet)
          db  = 0.005693d0*d2r * sin(l0-lam)*sin(bet)
          lam = lam + dl
          bet = bet + db
          
          ! Meeus, step 13 - correct for nutation in longitude:
          lam0 = lam0 + dpsi
          lam  = lam  + dpsi
       end if
       
       ! Meeus, step 14:
       call ecl_2_eq(lam,bet,eps,   ra,dec)
       call ecl_2_eq(lam0,bet0,eps, ra0,dec0)
       
       ! Meeus, step 15:
       pa = atan2(cos(dec0)*sin(ra0-ra),sin(dec0)*cos(dec) - cos(dec0)*sin(dec)*cos(ra0-ra))  ! Position angle
       
       if(loop.eq.1) then  ! Save data for the Sun's pov:
          pa_s = pa
          ar_s = ar
          br_s = br
       end if
    end do  ! loop=1,2 - from Sun, Earth
    
  end subroutine saturnphys
  !*********************************************************************************************************************************
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute low-accuracy planet positions.  Sun and Moon have dedicated routines, for planets abridged VSOP87 is used.
  !! 
  !! \param jd    Julian Day of computation
  !! \param pl    Planet number:  0: Moon,  1-2: Mer-Ven,  3: Sun,  4-9: Mar-Plu
  !!
  !! \param calc  Calculate  1: l,b,r,diam, 2: + ra,dec, 3: + gmst,agst, 4: + az,alt, 5: + elon, mag, k, pa, parang, hp, GC LBR,
  !!                         6: + topocentric positions (Sun and Moon only, optional - default = 5)
  !! \param nt    Number of terms to use for the calculation (has an effect for Moon only; nt<=60); 
  !!              a smaller nt gives faster, but less accurate results (optional, default=60)
  !!
  !! \param lat   Latitude of the observer (rad, optional)
  !! \param lon   Longitude of the observer (rad, optional)
  
  subroutine planet_position_la(jd, pl, calc,nt, lat,lon)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev,rev2
    
    use TheSky_planetdata, only: planpos
    use TheSky_moon, only: moonpos_la
    use TheSky_sun, only: sunpos_la
    use TheSky_coordinates, only: geoc2topoc_ecl, geoc2topoc_eq, eq2horiz, refract
    use TheSky_local, only: lat0,lon0
    
    implicit none
    real(double), intent(in) :: jd
    integer, intent(in) :: pl
    integer, intent(in), optional :: calc,nt
    real(double), intent(in), optional :: lat,lon
    
    integer :: lcalc,lnt
    real(double) :: dh_ref, parAng, dRA_ref
    
    lcalc = 5
    if(present(calc)) lcalc = calc
    
    if(present(lat)) lat0 = lat
    if(present(lon)) lon0 = lon
    
    
    planpos = 0.d0
    planpos(39) = dble(pl)
    
    select case(pl)
    case(0)  ! Moon
       
       lnt = 60
       if(present(nt)) lnt = nt
       call moonpos_la(jd, lcalc,lnt)
              
    case(3)  ! Sun
       
       call sunpos_la(jd, lcalc)
       
    case default  ! Planets
       
       call planet_position(jd,pl, LBaccur=1.d-6,Raccur=1.d-4, ltime=.false.)  ! Truncated VSOP87, no light-time correction
       
    end select
    
    
    ! Fill topocentric planpos elements (21-32) with geocentric values (1-12) for Moon and Sun:
    if(pl.eq.0 .or. pl.eq.3) then
       planpos(21:32) = planpos(1:12)
       planpos(31) = planpos(30)  ! Altitude, "corrected" for refraction - except that it is not corrected in the cheap version
       
       
       ! Simple correction for horizontal parallax, Meeus p.281, assuming rho=1:
       if(lcalc.eq.5) planpos(30) = planpos(10) - asin(sin(planpos(17))*cos(planpos(10)))  ! alt' = alt - Hp*cos(alt)
       
       ! More precise conversion to topocentric coordinates:
       if(lcalc.ge.6) then  ! This takes 15% more CPU time for the Moon (full accuracy; nt=60) and 110% more for the Sun!
          call geoc2topoc_ecl(planpos(1),planpos(2), planpos(4),planpos(12)/2.d0, planpos(48),planpos(44), & ! l,b, del,r, eps,lst
               planpos(21),planpos(22),planpos(32))  ! L,B, r
          planpos(32) = planpos(32)*2.d0                    ! Apparent diameter
          planpos(24) = planpos(4)*planpos(12)/planpos(32)  ! Distance
          
          call geoc2topoc_eq(planpos(5),planpos(6), planpos(4),planpos(8), planpos(25),planpos(26))  ! ra,dec, del,hh, topo ra,dec
          
          call eq2horiz(planpos(25),planpos(26),planpos(45), planpos(28),planpos(29),planpos(30))  ! topo ra,dec,agst, -> hh,az,alt
       end if
       
       if(lcalc.ge.5) then
          dh_ref = refract(planpos(30))       ! Corrected for atmospheric refraction in altitude
          planpos(31) = planpos(30) + dh_ref  ! Topocentric altitude, corrected for atmospheric refraction
          
          parAng = atan2(sin(planpos(28)),tan(lat0)*cos(planpos(26)) - sin(planpos(26))*cos(planpos(28)))  ! Topoc. paral. angle
          dRA_ref = dh_ref / cos(planpos(26)) * sin(parAng)
          planpos(25) = rev(planpos(25) + dRA_ref)  ! Correct topocentric RA for refraction
          planpos(28) = rev(planpos(28) - dRA_ref)  ! Correct topocentric HA for refraction
          planpos(26) = planpos(26) + dh_ref*cos(parAng)
       end if
    end if
    
  end subroutine planet_position_la
  !*********************************************************************************************************************************
  
  !*********************************************************************************************************************************
  !> \brief  Correction for comet magnitude due to forward and backscattering
  !! 
  !! \param phaseAng    Phase angle of the comet (rad)
  !!
  !! \retval            Magnitude correction, to be added to the magnitude (usually negative)
  !!
  !! \see https://ui.adsabs.harvard.edu/abs/2007ICQ....29..119M
  !!
    
  function comet_scatter_magnitude_correction(phaseAng)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi  ! , r2d
    use SUFR_angles, only: rev2
    
    implicit none
    real(double), intent(in) :: phaseAng
    real(double) :: comet_scatter_magnitude_correction, theta,d90,gf,gb,k,phi
    
    theta = rev2(pi - phaseAng)  ! Scattering angle = 180° - phase angle
    
    d90   =  0.3d0   ! Dust-to-gas light ratio in the coma at theta=90°.  1 for a "normal" comet, ~10 for a "dusty" one.  3.3 ~ half way.  0.3 taken from http://astro.vanbuitenen.nl/comet/2020F8
    gf    =  0.9d0   ! Forward-scattering asymmetry factor
    gb    = -0.6d0   ! Backscattering asymmetry factor
    k     =  0.95d0  ! Partition coefficient between forward and backscattering
    
    phi   = d90/(1+d90) * (k*((1+gf**2)/(1+gf**2 - 2*gf*cos(theta)))**1.5d0 + &
         k*((1+gb**2)/(1+gb**2 - 2*gb*cos(theta)))**1.5d0 + 1.d0/d90)            ! Eq.1
    
    comet_scatter_magnitude_correction = -2.5d0*log10(phi)
    
    ! write(*,'(99F9.3)') phaseAng*r2d,theta*r2d,d90,gf,gb,k,phi,comet_scatter_magnitude_correction
  end function comet_scatter_magnitude_correction
  !*********************************************************************************************************************************
  
  
end module TheSky_planets
!***********************************************************************************************************************************
