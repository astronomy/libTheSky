!> \file modules.f90  Modules used by libTheSky


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
!> \brief  Constants used in libTheSky

module TheSky_constants
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  integer, parameter :: deltat_nmax = 1000   ! Need ~430 until 2000
  real(double) :: deltat_values(deltat_nmax), deltat_years(deltat_nmax)
  real(double) :: deltat_accel, deltat_change, deltat_0, deltat_forced, jd1820
  integer :: deltat_n, deltat_minyr, deltat_maxyr, TheSky_verbosity
  
  real(double) :: nutationdat(9,63)
  
  character :: TheSkydir*(99),library_name*(99)
  
end module TheSky_constants
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Local parameters for libTheSky: location, date, time

module TheSky_local
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  integer, parameter :: maxlocs = 100
  real(double) :: lat0,lon0,height,deltat,tz,tz0,second,day
  integer :: year,month,hour,minute,dsttp
  
end module TheSky_local
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Planet data, needed to compute planet positions

module TheSky_planetdata
  use SUFR_kinds, only: double, long
  implicit none
  save
  private :: double, long
  
  integer, parameter :: nplanpos=100
  integer, parameter :: nasteroids=1000    ! Nasteroids is actually much larger; look at the first Nasteroids asteroids only
  
  integer(long) :: moonla_lrb(3,60)
  integer :: pluc(43,3),plul(43,2),plub(43,2),plur(43,2), pl0
  integer :: moonla_arg(8,60), VSOPnls(3,8), vsopNblk(0:5,3,8)
  
  real(double) :: VSOPdat(4,6827,10), VSOPtruncs(3,8)
  real(double) :: plelems(8,6),plelems2000(8,6),plelemdata(2,8,6,0:3)
  real(double) :: asterelems(nasteroids,9)
  
  !> \brief  Planpos[] is an array with many different types of coordinates and related variables:
  !!
  !! Geocentric (mostly) coordinates:
  !! - 1   =  gcl             =  Apparent geocentric ecliptic longitude
  !! - 2   =  gcb             =  Apparent geocentric ecliptic latitude
  !! - 3   =  r               =  True heliocentric distance
  !! - 4   =  delta           =  Apparent geocentric distance
  !! - 5   =  ra              =  Apparent geocentric R.A.
  !! - 6   =  dec             =  Apparent geocentric declination
  !! - 7   =  tau             =  Geocentric light time in days
  !! - 8   =  hh              =  Apparent geocentric hour angle
  !! - 9   =  az              =  Apparent geocentric azimuth
  !! - 10  =  alt             =  Apparent geocentric altitude
  !! - 11  =  elon            =  Apparent TOPOCENTRIC elongation, using topocentric altitude and correcting for refraction
  !! - 12  =  diam            =  Apparent geocentric apparent diameter
  !! 
  !! Magnitude, phase, parallax:
  !! - 13  =  magn            =  Apparent visual magnitude
  !! - 14  =  k               =  Illuminated fraction
  !! - 15  =  pa              =  Phase angle
  !! - 16  =  parang          =  Topocentric parallactic angle
  !! - 17  =  hp              =  Horizontal parallax
  !! 
  !! Topocentric (mostly) coordinates:
  !! - 21  =  topl            =  Apparent topocentric ecliptic longitude
  !! - 22  =  topb            =  Apparent topocentric ecliptic latitude
  !! - 23  =  delta0          =  True geocentric distance
  !! - 24  =  topdelta        =  Apparent topocentric distance
  !! - 25  =  topra           =  Apparent topocentric R.A.
  !! - 26  =  topdec          =  Apparent topocentric declination
  !!
  !! - 28 = tophh             =  Apparent topocentric hour angle
  !! - 29 = topaz             =  Apparent topocentric azimuth
  !! - 30 = topalt            =  Apparent topocentric altitude
  !! - 31 = topalt + refract  =  Apparent topocentric altitude, corrected for refraction
  !! - 32 = topdiam           =  Apparent topocentric apparent diameter
  !! 
  !! Heliocentric coordinates:
  !! - 33  =  l               =  Apparent heliocentric ecliptic longitude
  !! - 34  =  b               =  Apparent heliocentric ecliptic latitude
  !! - 35  =  r               =  Apparent heliocentric distance
  !! - 36  =  hcl00           =  True heliocentric ecliptic longitude, FK5
  !! - 37  =  hcb00           =  True heliocentric ecliptic latitude, FK5
  !! - 38  =  hcr00           =  True heliocentric distance
  !! 
  !! Other variables:
  !! - 39  =  pl              =  Planet for which data were computed
  !! - 40  =  jde             =  JDE for the instance of calculation
  !! - 41  =  l0+pi           =  True geocentric ecliptic longitude for the Sun, in FK5
  !! - 42  =  -b0             =  True geocentric ecliptic latitude for the Sun, in FK5
  !! - 43  =  r0              =  True geocentric distance for the Sun
  !! - 44  =  lst             =  Local APPARENT siderial time (radians)
  !! - 45  =  agst            =  Greenwich APPARENT siderial time (radians)
  !! - 46  =  t0 * 10.d0      =  Apparent dynamical time in Julian Centuries since 2000.0
  !! - 47  =  dpsi            =  Nutation in longitude
  !! - 48  =  eps             =  True obliquity of the ecliptic, corrected for nutation
  !! - 49  =  gmst            =  Greenwich MEAN siderial time (radians)
  !! - 50  =  eps0            =  Mean obliquity of the ecliptic, without nutation
  !! 
  !! More geocentric coordinates:
  !! - 51 = sun_gcl,b         =  Apparent geocentric longitude of the Sun (variable was treated as if pl.eq.3)
  !! - 52 = sun_gcl,b         =  Apparent geocentric latitude of the Sun (variable was treated as if pl.eq.3)
  !! 
  !! - 61 = gcx               =  Apparent geocentric x
  !! - 62 = gcy               =  Apparent geocentric y
  !! - 63 = gcz               =  Apparent geocentric z
  !! 
  !! - 64 = gcx0              =  True geocentric x
  !! - 65 = gcy0              =  True geocentric y
  !! - 66 = gcz0              =  True geocentric z
  !!
  !! - 67 = gcl0              =  True geocentric ecliptic longitude
  !! - 68 = gcb0              =  True geocentric ecliptic latitude
  !! - 69 = delta0            =  True geocentric distance
  
  real(double) :: planpos(nplanpos)
  
  character :: plcon(0:19)*(3),asternames(nasteroids)*(18)
  
end module TheSky_planetdata
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  ELP-82B Moon data, needed to compute Moon positions

module TheSky_moondata
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  real(double), parameter :: c1=1.d0/60.d0, c2=1.d0/3600.d0
  real(double), parameter :: ath=384747.9806743165d0, a0=384747.9806448954d0
  
  real(double) :: p1,p2,p3,p4,p5, q1,q2,q3,q4,q5
  real(double) :: w(3,0:4),eart(0:4),peri(0:4),p(8,0:1), del(4,0:4),zeta(0:1),t(0:4)
  real(double) :: pre(3),coef(7),zone(6), pc1(6,1023),pc2(6,918),pc3(6,704), per1(3,19537),per2(3,6766),per3(3,8924)
  
  integer :: ilu(4),ipla(11),nterm(3,12),nrang(3,0:12)
  
  real(double) :: prec0
  integer :: ideb
  
  data ideb/0/,prec0/-1.d0/,t/1.d0,4*0.d0/
  
end module TheSky_moondata
!***********************************************************************************************************************************


!***************************************************************************************************
! \brief  Constants for the ELP-MPP02 lunar theory

module TheSky_elp_mpp02_constants
  use SUFR_kinds, only: double
  
  implicit none
  real(double) :: w(3,0:4),eart(0:4),peri(0:4),  zeta(0:4),del(4,0:4)
  real(double) :: p(8,0:4),delnu,dele,delg,delnp,delep,dtasm,am,   p1,p2,p3,p4,p5,q1,q2,q3,q4,q5
  
end module TheSky_elp_mpp02_constants
!***************************************************************************************************


!***************************************************************************************************
! \brief  Data series for the ELP-MPP02 lunar theory

module TheSky_elp_mpp02_series
  use SUFR_kinds, only: double
  
  implicit none
  integer, parameter, private :: max1=2645, max2=33256
  integer :: nmpb(3,3), nper(3,0:3,3)
  real(double) :: cmpb(max1),fmpb(0:4,max1),   cper(max2),fper(0:4,max2)
  
end module TheSky_elp_mpp02_series
!***************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Data to compute comet positions

module TheSky_cometdata
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  integer, parameter :: nCometsMax=10000   ! Max. number of comets
  integer :: nComets, cometDatFile, cometDiedAtP(nCometsMax)
  real(double) :: cometelems(nCometsMax,9), comepoche
  character :: cometnames(nCometsMax)*(60)
  
end module TheSky_cometdata
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Star and basic constellation data

module TheSky_stardata
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  ! Selected bright stars close to the ecliptic:
  integer, parameter :: nstars=17,nconstel=88,nconid=357
  real(double) :: starra(nstars),stardec(nstars),starl(nstars),starb(nstars)
  real :: starmags(nstars),starrads(nstars)
  character :: starnames(nstars)*(10),starnamesnl(nstars)*(11),starcons(nstars)*(10),starconsnl(nstars)*(10),starconsabr(nstars)*(3)
  
  ! Constellations:
  integer :: conid(nconid)
  real(double) :: conidral(nconid),conidrau(nconid),coniddecl(nconid)
  character :: conabr(nconstel)*(3), conidabr(nconid)*(3)
  character :: latconnames(nconstel)*(19),genconnames(nconstel)*(19),nlconnames(nconstel)*(17),enconnames(nconstel)*(18)
  
end module TheSky_stardata
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Data from the Bright Star Catalogue (BSC)

module TheSky_BSCdata
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  integer, parameter :: n_bsc=9110
  integer :: bsc_sao(n_bsc),bsc_vm_indx(n_bsc)
  real(double) :: bsc_ra(n_bsc),bsc_dec(n_bsc),bsc_pma(n_bsc),bsc_pmd(n_bsc),bsc_rv(n_bsc),bsc_vm(n_bsc)
  real(double) :: bsc_par(n_bsc),bsc_bv(n_bsc),bsc_ub(n_bsc),bsc_ri(n_bsc)
  character :: bsc_name(n_bsc)*(10),bsc_abbr(n_bsc)*(10),bsc_mult(n_bsc),bsc_var(n_bsc),bsc_sptype(n_bsc)*(20)
  
end module TheSky_BSCdata
!***********************************************************************************************************************************


