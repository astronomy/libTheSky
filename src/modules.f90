!> \file modules.f90  Modules used by libTheSky


!  Copyright (c) 2002-2025  Marc van der Sluys - Nikhef/Utrecht University - marc.vandersluys.nl
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
  
  integer, parameter :: deltat_nmax = 1000    !< Maximum number of Delta-T measurements.  Need ~430 until 2000
  real(double) :: deltat_values(deltat_nmax)  !< Values of DeltaT
  real(double) :: deltat_years(deltat_nmax)   !< Years for DeltaT values
  real(double) :: deltat_accel                !< Acceleration for DeltaT parabola
  real(double) :: deltat_change               !< Change for DeltaT parabola
  real(double) :: deltat_0                    !< Zero point for DeltaT parabola
  real(double) :: deltat_forced               !< Forced value for DeltaT, overriding computation
  real(double) :: jd1820                      !< JD of 1820.0, for DeltaT
  
  integer :: deltat_n                         !< Actual number of DeltaT measurements
  integer :: deltat_minyr                     !< Start year of DeltaT measurements
  integer :: deltat_maxyr                     !< End year of DeltaT measurements
  
  integer :: TheSky_verbosity                 !< Verbosity of libTheSky output
  
  real(double) :: nutationdat(9,63)           !< Data for simple nutation function
  
  character :: TheSkyDatadir*(99)             !< Directory containing data files for libTheSky
  character :: library_name*(99)              !< Name of this library
  
end module TheSky_constants
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Local parameters for libTheSky: location, date, time

module TheSky_local
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  integer, parameter :: maxlocs = 100  !< Maximum number of observation locations
  real(double) :: lat0    !< Latitude of the observer (rad)
  real(double) :: lon0    !< Longitude of the observer (rad)
  real(double) :: height  !< Altitude of the observer above sea level (m)
  real(double) :: deltat  !< Current value of DeltaT (s)
  real(double) :: tz      !< Current value of time zone, taking into account DST (hours; >0 is east of Greenwich)
  real(double) :: tz0     !< Standard value of time zone, without DST ("winter time"; hours; >0 is east of Greenwich)
  real(double) :: second  !< Seconds of time of current instant
  real(double) :: day     !< Day of month of current instant, with decimals if desired
  
  integer :: year    !< Year CE of current instant
  integer :: month   !< Month of year of current instant
  integer :: hour    !< Hour of time of current instant
  integer :: minute  !< Minute of time of current instant
  integer :: dsttp   !< DST type for current location (0: none, 1: EU, 2: USA+Canada (after 2007)
  
end module TheSky_local
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Planet data, needed to compute planet positions

module TheSky_planetdata
  use SUFR_kinds, only: double, long
  implicit none
  save
  private :: double, long
  
  integer, parameter :: nasteroids=1000    !< Number of entries in the asteroids array.  Nasteroids is actually much larger; look at the first Nasteroids asteroids only
  
  integer :: moonla_arg(8,60)              !< Arguments for the low-accuracy (la) Moon position
  integer(long) :: moonla_lrb(3,60)        !< L,B,R data for the low-accuracy (la) Moon position
  
  integer :: pluc(43,3)  !< Constants for the periodic terms for the position of Pluto
  integer :: plul(43,2)  !< Constants for the longitude of Pluto
  integer :: plub(43,2)  !< Constants for the latitude of Pluto
  integer :: plur(43,2)  !< Constants for the distance of Pluto
  
  integer :: pl0         !< Remember a special planet
  
  integer :: VSOPnls(3,8)       !< Numbers of lines in the VSOP input files (l,b,r x 8 pl)
  integer :: vsopNblk(0:5,3,8)  !< Line number in the VSOP data where the next block of (Planet, Variable (LBR), Power) starts
  
  real(double) :: VSOPdat(4,6827,10)  !< Periodic terms for VSOP87
  real(double) :: VSOPtruncs(3,8)     !< Truncuate VSOP87 terms at these accuracies
  
  real(double) :: plelems(8,6)           !< Planet orbital elements for Equation of Data
  real(double) :: plelems2000(8,6)       !< Planet orbital elements for J2000
  real(double) :: plelemdata(2,8,6,0:3)  !< Data to compute planet orbital elements
  
  real(double) :: asterelems(nasteroids,9)  !< Asteroid orbital elements
  character :: asternames(nasteroids)*(18)  !< Names of the asteroids
  
  integer, parameter :: nplanpos=100       !< Number of entries in the planpos array
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
  character :: plcon(0:19)*(3)  !< Constellation abbreviation for the planets (0-9) and bright stars
  
end module TheSky_planetdata
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  ELP 2000-82B Moon data, needed to compute Moon positions

module TheSky_moondata
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  real(double), parameter :: c1=1.d0/60.d0            !< Constant for ELP 2000-82B theory (arcminutes to degrees)
  real(double), parameter :: c2=1.d0/3600.d0          !< Constant for ELP 2000-82B theory (arcseconds to degrees)
  real(double), parameter :: ath=384747.9806743165d0  !< Constant for ELP 2000-82B theory (orbital separation?)
  real(double), parameter :: a0=384747.9806448954d0   !< Constant for ELP 2000-82B theory (orbital separation?)
  
  real(double) :: p1  !< Precession sine coefficient for ELP 2000-82B theory
  real(double) :: p2  !< Precession sine coefficient for ELP 2000-82B theory
  real(double) :: p3  !< Precession sine coefficient for ELP 2000-82B theory
  real(double) :: p4  !< Precession sine coefficient for ELP 2000-82B theory
  real(double) :: p5  !< Precession sine coefficient for ELP 2000-82B theory
  
  real(double) :: q1  !< Precession cosine coefficient for ELP 2000-82B theory
  real(double) :: q2  !< Precession cosine coefficient for ELP 2000-82B theory
  real(double) :: q3  !< Precession cosine coefficient for ELP 2000-82B theory
  real(double) :: q4  !< Precession cosine coefficient for ELP 2000-82B theory
  real(double) :: q5  !< Precession cosine coefficient for ELP 2000-82B theory
  
  real(double) :: w(3,0:4)   !< Constants for mean longitude
  real(double) :: eart(0:4)  !< Earth-Moon barycentre (EMB) elements
  real(double) :: peri(0:4)  !< Mean longitude of the perihelion of the Earth-Moon barycentre (EMB)
  real(double) :: p(8,0:1)   !< Planetary arguments: mean longitudes and mean motions
  
  real(double) :: del(4,0:4)  !< Delaunay's variables (https://en.wikipedia.org/wiki/Orbital_elements#Delaunay_variables)
  real(double) :: zeta(0:1)   !< Mean longitude (w) + rate precession (?)
  real(double) :: t(0:4)      !< Array for time^0, ..., time^4
  
  real(double) :: pre(3)         !< CHECK: does this actually do anything?
  real(double) :: coef(7)        !< CHECK: Coefficients in ELP 2000-82B data file (float)
  real(double) :: zone(6)        !< CHECK: Something in ELP 2000-82B theory
  real(double) :: pc1(6,1023)    !< CHECK: Something in ELP 2000-82B theory
  real(double) :: pc2(6,918)     !< CHECK: Something in ELP 2000-82B theory
  real(double) :: pc3(6,704)     !< CHECK: Something in ELP 2000-82B theory
  real(double) :: per1(3,19537)  !< CHECK: Something in ELP 2000-82B theory
  real(double) :: per2(3,6766)   !< CHECK: Something in ELP 2000-82B theory
  real(double) :: per3(3,8924)   !< CHECK: Something in ELP 2000-82B theory
  
  integer :: ilu(4)              !< CHECK: Coefficients in ELP 2000-82B data file (integer)
  integer :: ipla(11)            !< CHECK: Coefficients in ELP 2000-82B data file (integer)
  integer :: nterm(3,12)         !< CHECK: Number of terms? in ELP 2000-82B data file
  integer :: nrang(3,0:12)       !< CHECK: Number of terms? in ELP 2000-82B data file
  
  real(double) :: prec0          !< CHECK: Something in ELP 2000-82B theory
  integer :: ideb                !< Memorise whether this routine has been run before
  
  data ideb/0/, prec0/-1.d0/, t/1.d0,4*0.d0/
  
end module TheSky_moondata
!***********************************************************************************************************************************


!***************************************************************************************************
! \brief  Constants for the ELP-MPP02 lunar theory

module TheSky_elp_mpp02_constants
  use SUFR_kinds, only: double
  
  implicit none
  real(double) :: w(3,0:4)    !< Constants for mean longitude
  real(double) :: eart(0:4)   !< Earth-Moon barycentre (EMB) elements
  real(double) :: peri(0:4)   !< Mean longitude of the perihelion of the Earth-Moon barycentre (EMB)
  real(double) :: zeta(0:4)   !< Mean longitude (w) + rate precession (?)
  real(double) :: del(4,0:4)  !< Delaunay's variables (https://en.wikipedia.org/wiki/Orbital_elements#Delaunay_variables)
  
  real(double) :: p(8,0:4)   !< Planetary arguments: mean longitudes, mean motions, ...?
  
  real(double) :: delnu  !< Corrections of the constants (fit to DE200/LE200) (?)
  real(double) :: dele   !< Corrections of the constants (fit to DE200/LE200) (?)
  real(double) :: delg   !< Corrections of the constants (fit to DE200/LE200) (?)
  real(double) :: delnp  !< Corrections of the constants (fit to DE200/LE200) (?)
  real(double) :: delep  !< Corrections of the constants (fit to DE200/LE200) (?)
  
  real(double) :: am      !< Ratio of the mean motions (EMB / Moon)
  real(double) :: dtasm   !< 2/3 of product of (ratio of the semi-major axes) x (ratio of the mean motions) (Moon / EMB) ???
  
  real(double) :: p1  !< Precession sine coefficient for ELP MPP02 theory
  real(double) :: p2  !< Precession sine coefficient for ELP MPP02 theory
  real(double) :: p3  !< Precession sine coefficient for ELP MPP02 theory
  real(double) :: p4  !< Precession sine coefficient for ELP MPP02 theory
  real(double) :: p5  !< Precession sine coefficient for ELP MPP02 theory
  
  real(double) :: q1  !< Precession cosine coefficient for ELP MPP02 theory
  real(double) :: q2  !< Precession cosine coefficient for ELP MPP02 theory
  real(double) :: q3  !< Precession cosine coefficient for ELP MPP02 theory
  real(double) :: q4  !< Precession cosine coefficient for ELP MPP02 theory
  real(double) :: q5  !< Precession cosine coefficient for ELP MPP02 theory
  
  
end module TheSky_elp_mpp02_constants
!***************************************************************************************************


!***************************************************************************************************
! \brief  Data series for the ELP-MPP02 lunar theory

module TheSky_elp_mpp02_series
  use SUFR_kinds, only: double
  
  implicit none
  integer, parameter, private :: max1=2645, max2=33256
  integer :: nmpb(3,3)      !< Number of lines in the file for the main problem
  integer :: nper(3,0:3,3)  !< Number of lines in the file for the perturbation
  
  real(double) :: cmpb(max1)      !< Coefficients? for the main problem
  real(double) :: fmpb(0:4,max1)  !< Factors? for the main problem
  real(double) :: cper(max2)      !< Coefficients? for the perturbation series
  real(double) :: fper(0:4,max2)  !< Factors? for the perturbation series
  
end module TheSky_elp_mpp02_series
!***************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Data to compute comet positions

module TheSky_cometdata
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  integer, parameter :: nCometsMax=10000    !< Size of comet database
  integer :: nComets                        !< Actual number of comets in database
  integer :: cometDatFile                   !< Data file to use 1: comets.dat (MANY comets, no magnitude info), 2: comets_mpc.dat (currently visible comets + magn. info)
  logical :: cometDiedAtP(nCometsMax)       !< This comet died at perihelion (true/false)
  real(double) :: cometelems(nCometsMax,9)  !< Orbital elements of the comets: 1: JD of epoch (often J2000), 2: Perihelion distance (AU?), 3: Eccentricity, 4: Inclination, 5: Argument of perihelion (omega), 6: Longitude of ascending node (OMEGA) i, 7: JD of perihelion
  real(double) :: comepoche                 !< JD of epoch (often J2000) == cometelems(i,1)
  character :: cometnames(nCometsMax)*(60)  !< Names of the comets
  
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
  integer, parameter :: nstars=17    !< Number of bright stars close to the ecliptic (including Pleiades)
  
  real(double) :: starra(nstars)   !< Right ascensions of the bright, ecliptical stars
  real(double) :: stardec(nstars)  !< Declinations of the bright, ecliptical stars
  real(double) :: starl(nstars)    !< Ecliptic longitudes of the bright, ecliptical stars
  real(double) :: starb(nstars)    !< Ecliptic latitudes of the bright, ecliptical stars
  
  real :: starmags(nstars)  !< Magnitudes of the bright, ecliptical stars
  real :: starrads(nstars)  !< Radii of the bright, ecliptical stars (only non-zero for Pleiades)
  
  character :: starnames(nstars)*(10)    !< English names of the bright, ecliptical stars
  character :: starnamesnl(nstars)*(11)  !< Dutch names of the bright, ecliptical stars
  character :: starcons(nstars)*(10)     !< Latin/English constellation names for the bright, ecliptical stars
  character :: starconsnl(nstars)*(10)   !< Dutch constellation names for the bright, ecliptical stars
  character :: starconsabr(nstars)*(3)   !< Constellation abbreviations for the bright, ecliptical stars
  
  ! Constellations:
  integer, parameter :: nconstel=88  !< Number of constellations
  integer, parameter :: nconid=357   !< Number of data points for constellation ID 
  integer :: conid(nconid)           !< Constellation ID
  
  real(double) :: conidral(nconid)   !< Constellation lower RA boundary for ID
  real(double) :: conidrau(nconid)   !< Constellation uppwer RA boundary for ID
  real(double) :: coniddecl(nconid)  !< Constellation lower declination boundary for ID
  
  character :: conabr(nconstel)*(3)  !< Abbreviations of the constellations
  character :: conidabr(nconid)*(3)  !< Abbreviations of the constellation IDs
  
  character :: latconnames(nconstel)*(19)  !< Latin constellation names
  character :: genconnames(nconstel)*(19)  !< Genitives of the Latin constellation names
  character :: nlconnames(nconstel)*(17)   !< Dutch constellation names
  character :: enconnames(nconstel)*(18)   !< English constellation names
  
end module TheSky_stardata
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Data from the Bright Star Catalogue (BSC)

module TheSky_BSCdata
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  integer, parameter :: n_bsc=9110  !< Size of the Bright Star Catalogue (BSC)
  integer :: bsc_sao(n_bsc)         !< SAO numbers of the BSC stars
  integer :: bsc_vm_indx(n_bsc)     !< Index array for sorting to visual magnitude
  
  real(double) :: bsc_ra(n_bsc)     !< Right ascensions of the BSC stars
  real(double) :: bsc_dec(n_bsc)    !< Declinations of the BSC stars
  real(double) :: bsc_pma(n_bsc)    !< Proper motions in RA of the BSC stars
  real(double) :: bsc_pmd(n_bsc)    !< Proper motions in declination of the BSC stars
  real(double) :: bsc_rv(n_bsc)     !< Radial velocities of the BSC stars
  real(double) :: bsc_vm(n_bsc)     !< Visual magnitudes of the BSC stars
  
  real(double) :: bsc_par(n_bsc)    !< Parallaxes of the BSC stars
  real(double) :: bsc_bv(n_bsc)     !< B-V colours of the BSC stars
  real(double) :: bsc_ub(n_bsc)     !< U-B colours of the BSC stars
  real(double) :: bsc_ri(n_bsc)     !< R-I colours of the BSC stars
  
  character :: bsc_name(n_bsc)*(10)    !< Proper names of the BSC stars
  character :: bsc_abbr(n_bsc)*(10)    !< Abbreviated names/codes for the BSC stars
  character :: bsc_mult(n_bsc)         !< Multiplicity codes for the BSC stars
  character :: bsc_var(n_bsc)          !< Variability codes for the BSC stars
  character :: bsc_sptype(n_bsc)*(20)  !< Spectral types of the BSC stars
  
end module TheSky_BSCdata
!***********************************************************************************************************************************


