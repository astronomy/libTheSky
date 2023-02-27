!> \file nutation.f90  Routines to compute nutation for libTheSky
!!
!! \todo  Check difference in outcome between rev.198 and 199


!  Copyright (c) 2002-2023  Marc van der Sluys - marc.vandersluys.nl
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
!> \brief Procedures for nutation

module TheSky_nutation
  implicit none
  save
  
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate nutation - cheap routine from Meeus
  !! 
  !! \param t     Time in Julian Millennia since 2000.0 in dynamical time
  !! \param dpsi  Nutation in longitude (output)
  !! \param eps0  Obliquity of the ecliptic (output)
  !! \param deps  Nutation in obliquity (output)
  !!
  !! \see  Meeus, Astronomical Algorithms, 1998, Ch. 22
  
  subroutine nutation(t, dpsi,eps0,deps)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi
    use TheSky_constants, only: nutationdat
    
    implicit none
    real(double), intent(in) :: t
    real(double), intent(out) :: dpsi,eps0,deps
    real(double) :: tt,tt2,tt3,u, d,ms,mm,f,omg,tmpvar, nu(9,63), conv
    integer :: i
    
    tt  = t*10.d0  ! Julian Centuries since 2000.0 in dynamical time
    tt2 = tt**2
    tt3 = tt*tt2
    
    d   = 5.19846946025d0  +  7771.37714617d0 *tt  -  3.340909d-5    *tt2  +  9.2114446d-8     *tt3  ! D in Meeus, p.144
    ms  = 6.24003588115d0  +  628.301956024d0 *tt  -  2.79776d-6     *tt2  -  5.8177641733d-8  *tt3  ! M in Meeus, p.144
    mm  = 2.3555483693d0   +  8328.69142288d0 *tt  +  1.517947757d-4 *tt2  +  3.102807559d-7   *tt3  ! M' in Meeus, p.144
    f   = 1.62790192912d0  +  8433.46615832d0 *tt  -  6.42717497d-5  *tt2  +  5.3329949d-8     *tt3  ! F in Meeus, p.144
    omg = 2.18243858558d0  -  33.7570459367d0 *tt  +  3.6142278d-5   *tt2  +  3.87850944888d-8 *tt3  ! Omega in Meeus, p.144
    
    dpsi=0.d0
    deps=0.d0
    nu = nutationdat
    do i=1,63  ! Use data from the IAU1980 model (Seidelmann 1981; see nutation.dat)
       tmpvar = nu(1,i)*d + nu(2,i)*ms + nu(3,i)*mm + nu(4,i)*f + nu(5,i)*omg
       dpsi = dpsi + (nu(6,i) + nu(7,i)*tt) * sin(tmpvar)
       deps = deps + (nu(8,i) + nu(9,i)*tt) * cos(tmpvar)
    end do
    
    conv = pi/(1.d4*3600.d0*180.d0)    ! Convert from 0.0001" to radians
    dpsi = dpsi*conv
    deps = deps*conv
    
    u   = t/10.d0  ! Dynamical time since 2000 in units of 10,000 years
    eps0 = 0.409092804222d0 - 0.022693789d0*u - 7.5146d-6*u*u + 0.0096926375d0*u**3 - 2.49097d-4*u**4 - 0.0012104343d0*u**5 -  &
         1.893197d-4*u**6 + 3.452d-5*u**7 + 1.3512d-4*u**8 + 2.8071d-5*u**9 + 1.1878d-5*u**10  ! Laskar, A&A 157 59 (1986), Tab.8.
    
  end subroutine nutation
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute nutation using the IAU 2000 model
  !! 
  !! \param jd        Julian day of computation
  !! \param dpsi_tot  Total nutation in longitude (output)
  !! \param deps_tot  Total nutation in obliquity (output)
  !!
  !! \see http://geoweb.mit.edu/~tah/mhb2000/
  
  subroutine nutation2000(jd, dpsi_tot, deps_tot)
    use SUFR_kinds, only: double
    implicit none
    
    !     Subroutine to compute the complete MHB_2000 nutation series
    !     with the associated corrections for planetary nutations,
    !     the freely excited free-core-nutation (valid for 1979-2001.42),
    !     the precession constant change and a rate of change of oblquity.)
    
    ! USAGE:
    !     call  MHB_2000( jd, dpsi_ls, deps_ls, dpsi_plan, deps_plan,
    !    .        dpsi_fcn , deps_fcn , dPsi_prec, deps_prec,
    !    .        dpsi_tot , deps_tot )
    
    ! where:
    !     <jd>    is the full julian date including fractional part of
    !             of the day (REAL(DOUBLE) input)
    !     <dpsi_ls> and <deps_ls> are the luni-solar nutation in
    !             longitude and oblquity (mas) (REAL(DOUBLE) OUTPUT)
    !     <dpsi_plan> and <deps_plan> are the contributions to the
    !             nutations in longitude and obliquity due direct
    !             planetary nutations and the perturbations of the
    !             lunar and terrestrial orbits (mas). (REAL(DOUBLE) OUTPUT)
    !     <dpsi_fcn> and <deps_fcn> are the contributions to the
    !             nutations in longitude and obliquity due the free-
    !             excitation of the Free-core-nutation (mas).  These
    !             values are valid for 1979-2000. (REAL(DOUBLE) OUTPUT)
    !     <dpsi_prec> and <deps_prec> are the contributions to the
    !             nutations in longitude and obliquity due changes in
    !             the precession constant and rate of change of
    !             obliquity (mas) (REAL(DOUBLE) OUTPUT).
    !     <dpsi_tot> and <deps_tot> are the total nutations in longitude
    !             and obliquity including the correction for the precession
    !             constant (when precession is computed using the IAU 1976
    !             precession constant), and are obtained by summing all
    !             of the above corrections (mas) (REAL(DOUBLE) OUTPUT).
    
    ! RESTRICTIONS: if <jd> is less than 2000000.0 this routine
    !               assumes an MJD has been passed and the time
    !               used will be converted to JD.  A warning
    !               message will be printed.  See individual modules
    !               for further restrictions.
    
    ! PASSED VARIABLES
    !
    ! INPUT Values
    ! jd     - Time at which value needed. (jd + fraction of day)
    
    ! OUTPUT Values
    !     dpsi_ls and deps_ls      - luni-solar nutation in
    !             longitude and oblquity (mas) (REAL* OUTPUT)
    !     dpsi_plan and deps_plan  - contributions to the
    !             nutations in longitude and obliquity due direct
    !             planetary nutations and the perturbations of the
    !             lunar and terrestrial orbits (mas). (REAL* OUTPUT)
    !     dpsi_fcn and deps_fcn    - contributions to the
    !             nutations in longitude and obliquity due the free-
    !             excitation of the Free-core-nutation (mas).  These
    !             values are valid for 1988-1994. (REAL* OUTPUT)
    !     dpsi_prec and deps_prec  - contributions to the
    !             nutations in longitude and obliquity due changes in
    !             the precession constant and rate of change of
    !             obliquity (mas) (REAL* OUTPUT).
    !     dpsi_tot and deps_tot    - total nutations in longitude
    !             and obliquity including the correction for the precession
    !             constant (when precession is computed using the IAU 1976
    !             precession constant), and are obtained by summing all
    !             of the above corrections (mas) (REAL* OUTPUT).
    
    
    real(double), intent(in) :: jd
    real(double), intent(out) :: dpsi_tot, deps_tot
    real(double) :: dpsi_ls, deps_ls, dpsi_plan, deps_plan, dpsi_fcn ,  deps_fcn, dpsi_prec, deps_prec, pi, mas2r
    
    !---------------------------------------------------------------
    
    !     Call each of the routines needed for each contribution.
    
    !     Luni-solar nutation
    call ls_nut( jd, dpsi_ls, deps_ls )
    
    !     Planetary nutation
    call plan_nut ( jd, dpsi_plan, deps_plan )
    
    !     Freely excited FCN (NOTE: No warning message is printed
    !     if the JD is out of the range of 1979-2000)
    call fcn_nut ( jd, dpsi_fcn , deps_fcn )
    
    !     Precession and obliquity rate contributions (NOTE: IAU-1976
    !     precession constant assumed to be used in the basic calculation
    !     of precession).
    
    call prec_nut( jd, dpsi_prec, deps_prec )
    
    !     Now add up all of the terms to get the total nutation angles
    
    dpsi_tot = dpsi_ls + dpsi_plan + dpsi_fcn + dpsi_prec
    deps_tot = deps_ls + deps_plan + deps_fcn + deps_prec
    
    ! Convert from mas to radians:
    pi = 4*atan(1.d0)
    mas2r = pi/(180.d0*3.6d6)
    dpsi_tot = dpsi_tot*mas2r
    deps_tot = deps_tot*mas2r
    
  end subroutine nutation2000
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  
  !TITLE LS_NUT
  
  !*********************************************************************************************************************************
  subroutine ls_nut( jd, dpsi_ls, deps_ls )
    use SUFR_kinds, only: double
    implicit none
    
    !     Routine to compute the MHB_2000 luni-solar contributions
    !     the nutations in longitude and obliquity.  The MHB_2000 is
    !     based on:
    
    !     (1) The Souchay and Kinoshita Rigid Earth nutation series
    !     SKRE1997.   There are many duplicate argument terms in this
    !     series and all of these have been compacted in signle argument
    !     terms.
    
    !     (2) Value of the Retrograde FCN resonance factors from
    !     the Mathews et al., 2000 , nutation formulation (full complex
    !     estimates, and scaling parameter R from the same
    !     theory.  
    
    !     (3) The effects of annual modulation of geodetic precession.
    !     The correction applied is
    !      0  1  0  0  0 -0.150 (correction to in-phase nutation in
    !                         longitude).
    
    !     (4) A prograde annual nutation has been estimated along with the
    !      resonance coefficients.  This probably reflects the influence of
    !      S1 atmospheric tide.
    
    !     (5) The free RFCN mode was estimated once every two years for the
    !      data after 1984.  (See values commented in eval_ls_nut.  For the
    !      last 6 years the values seem to be resonably stable.  
    
    !     (6) The new Simons et al., fundamental arguments are used in this
    !      version.  (The largest change from KSV_1995_1 was 0.007 mas for
    !      semiannual nutation.  All other changes, including the 18.6 year
    !      nutation were 0.001-0.002 mas.)
    
    !     REFERENCES:
    ! NEW Version based on: Corrections and new developments in rigid Earth
    !     nutation theory: Lunisolar influence including indirect planetary
    !     effects, J. Souchay and H. Kinioshita, Astron. and Astrophys., 1995.
    ! (Version here based on data files: SKRE1997.DPSI and SKRE1997.DEPS
    !  and generated with ks_plan.f)
    !     Souchay, J., and H. Kinoshita, Corrections and new developments in
    !         rigid Earth nutation theory: I. Lunisolar influence including indirect
    !         planetary effects, Astron. Astrophys.,  312, 1017--1030, 1996.
    !     Souchay, J., and H. Kinoshita, Corrections and new developments in
    !         rigid Earth nutation theory: II. Influence of second-order geopotential
    !         and direct planetray effect, Astron. Astrophys., 318, 639--652, 1997.
    !     Souchay, J., B. Loysel, H, Kinoshita, and M. Folgueira, Corrections
    !         and new developments in rigid Earth nutation theory: III. final tables
    !         REN-2000 including crossed-nutation and spin-orbit coupling effects,
    !         Astron. Astrophys. Suppl., 135, 111-131, 1999.    
    !     Kinoshita, H., and J. Souchay, The theory of the nutations for
    !         the rigid Earth at the second order, Celes. Mech. and Dynam.
    !         Astron., 48, 187--266, 1990.
    !     Mathews, P. M., B. A. Buffett, and T. A. Herring, Modeling of Nutation-Precession:
    !         Insights into the Earth's Interior and New Nonrigid Earth Nutation Series
    !         to be submitted, J. Geophys. Res, 2000.
    !     Simon, J. L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    !         Francou, G., Laskar, J., 1994, "Numerical Expressions for 
    !         Precession Formulae and Mean Elements for the Moon and
    !         Planets," Astron. Astrophys., 282, pp. 663-683.
    
    
    
    ! USAGE:
    !     call ls_nut( jd, dpsi_ls, deps_ls )
    !     where <jd>    is a full julian date with fractional part
    !                   of the day added (REAL(DOUBLE) INPUT)
    !     and <dpsi_ls> and <deps_ls> are the nutations
    !                   in longitude and obliquity in milliarcsec.
    !                   (REAL(DOUBLE) OUTPUT)
    
    ! RESTRICTIONS: if <jd> is less than 2000000.0 this routine
    !               assumes an MJD has been passed and the time
    !               used will be converted to JD.  A warning
    !               message will be printed.
    
    ! PASSED VARIABLES
    ! 
    ! INPUT Values
    ! jd     - Time at which value needed. (jd + fraction of day)
    
    ! OUTPUT Values
    ! dpsi_ls  - The nutation in longitude (mas).
    ! deps_ls  - The nutation in obliquity (mas).
    
    
    
    real(double), intent(in) :: jd
    real(double), intent(out) :: dpsi_ls, deps_ls
    
    ! LOCAL VARIABLES
    
    !   epoch       - Julian date (jd passed in unless the JD
    !                 appears to be an MJD in which case it is
    !                 converted to JD (2 400 000.5d0 added)
    !   ls_arg(5)   - The arguments for the Luni-solar nutations.
    !                 (l, l', F, D and Omega).  All in Radians.
    
    
    real(double) :: epoch, ls_arg(5)
    
    !***** Check to make sure user passed JD and not MJD.  Correct
    !     problem and warn the user.
    ! MvdS: remove this 'solution'
    !if( jd .lt.2000000.0d0  ) then
    !             write(*,100) jd
    !    100      format('**WARNING** MJD apparently passed to SD_COMP',
    !        .          ' Value (',F10.2,') converted to JD')
    !   epoch = jd + 2 400 000.5d0
    !else
    !   epoch = jd
    !end if
    
    epoch = jd
    
    !***** Get the fundamental arguments at this epoch
    
    call ls_angles( epoch, ls_arg)
    
    !     Now compute the luni-solare nutations by summing over all
    !     terms in the series.
    
    call eval_ls_nut( epoch, ls_arg, dpsi_ls, deps_ls )
    
  end subroutine ls_nut
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  
  !TITLE LS_ANGLES
  
  !*********************************************************************************************************************************
  subroutine ls_angles( epoch, ls_arg )
    use SUFR_kinds, only: double
    implicit none
    
    !     Routine to compute the value of the fundamental argument
    !     for Brown's arguments.  Arguments based on the IERS
    !     standards.
    
    ! MOD TAH 960206: Changed arguments to use Simons et al., 1994 
    ! values:
    !     Simon, J. L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    !          Francou, G., Laskar, J., 1994, "Numerical Expressions for 
    !          Precession Formulae and Mean Elements for the Moon and
    !          Planets," Astron. Astrophys., 282, pp. 663-683.
    
    
    ! PHYSICAL CONSTANTS
    
    !   pi          - Define here to full precision
    !   rad_to_deg  - Conversion from radians to degs.
    !   DJ2000      - Julian date of J2000
    !   sec360      - number of seconds in 360 degreees.
    
    
    real(double) :: pi, rad_to_deg, DJ2000, sec360
    
    parameter ( pi            = 3.1415926535897932D0 )
    parameter ( DJ2000        = 2451545.d0           )
    parameter ( sec360        = 1296000.d0           )
    
    !     Computed quanities
    parameter ( rad_to_deg    = 180.d0   /pi         )
    
    !-------------------------------------------------------------------
    
    ! PASSED VARIABLES
    
    ! INPUT
    ! epoch  - Julian date for arguments (jd + fraction of day, REAL(DOUBLE))
    
    ! OUTPUT
    ! ls_arg(5) -  Brown's arguments (radians, REAL(DOUBLE))
    
    
    real(double), intent(in) :: epoch
    real(double), intent(out) :: ls_arg(5)
    
    ! LOCAL VARIABLES
    !      cent             - Julian centuries to DJ2000.
    !      el,eld           - Mean longitude of moon minus mean
    !                       - longitude of moon's perigee (arcsec)
    !      elc(5)           - Coefficients for computing el
    !      elp,elpd         - Mean longitude of the sun minus mean
    !                       - longitude of sun perigee (arcsec)
    !      elpc(5)          - Coeffiecents for computing elp
    !      f,fd             - Moon's mean longitude minus omega (sec)
    !      fc(5)            - Coefficients for computing f
    !      d,dd             - Mean elongation of the moon from the
    !                       - sun (arcsec)
    !      dc(5)            - coefficients for computing d
    !      om,omd           - longitude of the ascending node of the
    !                       - moon's mean orbit on the elliptic
    !                       - measured from the mean equinox of date
    !      omc(5)           - Coefficients for computing om.
    
    
    real(double) :: cent, el, elc(5), elp, elpc(5), f, fc(5), d, dc(5), om, omc(5)  ! ,eld , elpd ,fd ,dd ,omd
    
    !****  DATA statements for the fundamental arguments.
    !     Simons et al., 1994 values
    
    data elc    /    -0.00024470d0,    0.051635d0,   31.8792d0,  1717915923.2178d0,   485868.249036d0/
    data elpc   /    -0.00001149d0,    +0.000136d0,  -0.5532d0,  129596581.0481d0,   1287104.79305d0/
    data fc     /     0.00000417d0,    -0.001037d0,  -12.7512d0, 1739527262.8478d0,    335779.526232d0/
    data dc     /    -0.00003169d0,     0.006593d0,   -6.3706d0, 1602961601.2090d0,   1072260.70369d0/
    ! MOD TAH MHB_2000: 960606: Replaced <Om> with expression from b.3 of 
    !     Simon et al., 1994 since b.3 is computed with new precession constant
    !     (Only the rate changes).   
    data omc    /    -0.00005939d0,       0.007702d0,    7.4722d0,  -6962890.5431d0,     450160.398036d0/
    
    
    !****  Get the number of centuries to current time
    
    cent = (epoch-dj2000) / 36525.d0
    
    !****  Compute angular arguments and their time derivatives
    ! New formulas adding in the higher order term.
    
    el = elc(1) * cent**4 + elc(2) * cent**3 + elc(3) * cent**2 + elc(4) * cent + elc(5)
    el = mod( el, sec360 )
    !eld = 4.d0 * elc(1) * cent**3 + 3.d0 * elc(2) * cent**2 +  2.d0 * elc(3) * cent    +        elc(4) 
    
    elp = elpc(1) * cent**4 + elpc(2) * cent**3 + elpc(3) * cent**2 + elpc(4) * cent + elpc(5)
    elp = mod( elp, sec360 )
    !elpd = 4.d0 * elpc(1) * cent**3 + 3.d0 * elpc(2) * cent**2 +  2.d0 * elpc(3) * cent    +        elpc(4) 
    
    f = fc(1) * cent**4 + fc(2) * cent**3 + fc(3) * cent**2 + fc(4) * cent + fc(5)
    f = mod( f, sec360 )
    !fd = 4.d0 * fc(1) * cent**3 + 3.d0 * fc(2) * cent**2 +  2.d0 * fc(3) * cent    +        fc(4) 
    
    d = dc(1) * cent**4 + dc(2) * cent**3 + dc(3) * cent**2 + dc(4) * cent + dc(5)
    d = mod( d, sec360 )
    !dd = 4.d0 * dc(1) * cent**3 + 3.d0 * dc(2) * cent**2 +  2.d0 * dc(3) * cent    +        dc(4) 
    
    om = omc(1) * cent**4 + omc(2) * cent**3 + omc(3) * cent**2 + omc(4) * cent + omc(5)
    om = mod( om, sec360 )
    !omd = 4.d0 * omc(1) * cent**3 + 3.d0 * omc(2) * cent**2 +  2.d0 * omc(3) * cent    +        omc(4) 
    
    
    
    !****  Now save the values.  Convert values from arcseconds to radians
    
    ls_arg(1) = el / (3600.d0*rad_to_deg)
    ls_arg(2) = elp/ (3600.d0*rad_to_deg)
    ls_arg(3) = f  / (3600.d0*rad_to_deg)
    ls_arg(4) = d  / (3600.d0*rad_to_deg)
    ls_arg(5) = om / (3600.d0*rad_to_deg)
    
  end subroutine ls_angles
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  !TITLE EVAL_LS_NUT
  
  !*********************************************************************************************************************************
  subroutine eval_ls_nut( epoch, ls_arg, dpsi_ls, deps_ls )
    use SUFR_kinds, only: double
    implicit none
    
    !     Routine to compute the nutations in longitude and obliquity
    !     by summing over all terms in the nutations series.
    
    ! NOTE: ls_angles must be called before routine.
    
    ! PARAMETERS:
    
    ! num_ls  - Number of terms in the nutations series
    
    integer :: num_ls
    
    parameter ( num_ls      =  678)
    
    !   DJ2000      - Julian date of J2000
    !   pi          - Pi.
    
    
    real(double) :: pi, DJ2000
    
    parameter ( pi            = 3.1415926535897932D0 )
    parameter ( DJ2000        = 2451545.d0           )
    
    ! PASSED PARAMETERS:
    
    ! INPUT:
    ! epoch      - Julian date at which nutation angles are needed.
    ! ls_arg(5)  - Five arguments for the nutatoions (l,l',F,D and Om)
    !              computed at the epoch that the nutations need to be
    !              evaluated (rad) (REAL(DOUBLE))
    
    ! OUTPUT:
    ! dpsi_ls, deps_ls   - nutations in longitude and obliquity (mas)
    !              (REAL(DOUBLE))
    
    
    real(double), intent(in) :: epoch, ls_arg(5)
    real(double), intent(out) :: dpsi_ls, deps_ls
    
    ! LOCAL VARIABLES:
    
    !  i and j   - Counters to loop over the coeffients and the argumemts
    
    
    integer :: i,j
    
    !  arg       - Final summed argumemt for the nutations
    !              contributions (rads)
    !  cent      - Number of centuries since J2000.
    !  dpsi_lsu and deps_lsu - Nutations in longitude and oblquity
    !              in micro-arc-sec (units that the data statements
    !              are in)
    !  carg, sarg -- Cosine and sine of arguments.
    
    
    real(double) :: arg, cent, dpsi_lsu, deps_lsu, carg, sarg
    
    !     RFCN Freq.  -1.00231810920 cyc per sidreal day, Period   430.2082 solar days
    !  Units now in 0.1 microarcsec
    !  IX01-IX68(11,10)  -- Invidual declarations of the coefficents
    !                           of the nutation series so no data statement has more than 10 lines
    !                           The first 5 values are the arguments for l lp F D Om
    !                           The remaining elements are:
    !                            6 - Nutation in longitude psi (sin, uas)
    !                            7 - dpsi/dt (uasec/cent)
    !                            8 - Nutation in oblquity eps (cos, uas)
    !                            9 - deps/dt (uas/cent)
    !                           10 - Out-of-phase longitude (cos, uas)
    !                           11 - Out-of-phase obliquity (sin, uas)
    
    integer :: IX01(11,10), IX02(11,10), IX03(11,10), IX04(11,10),  &
         IX05(11,10), IX06(11,10), IX07(11,10), IX08(11,10),  IX09(11,10), IX10(11,10), IX11(11,10), IX12(11,10),  &
         IX13(11,10), IX14(11,10), IX15(11,10), IX16(11,10),  IX17(11,10), IX18(11,10), IX19(11,10), IX20(11,10),  &
         IX21(11,10), IX22(11,10), IX23(11,10), IX24(11,10),  IX25(11,10), IX26(11,10), IX27(11,10), IX28(11,10),  &
         IX29(11,10), IX30(11,10), IX31(11,10), IX32(11,10),  IX33(11,10), IX34(11,10), IX35(11,10), IX36(11,10),  &
         IX37(11,10), IX38(11,10), IX39(11,10), IX40(11,10),  IX41(11,10), IX42(11,10), IX43(11,10), IX44(11,10),  &
         IX45(11,10), IX46(11,10), IX47(11,10), IX48(11,10),  IX49(11,10), IX50(11,10), IX51(11,10), IX52(11,10),  &
         IX53(11,10), IX54(11,10), IX55(11,10), IX56(11,10),  IX57(11,10), IX58(11,10), IX59(11,10), IX60(11,10),  &
         IX61(11,10), IX62(11,10), IX63(11,10), IX64(11,10),  IX65(11,10), IX66(11,10), IX67(11,10), IX68(11, 8)
    
    integer :: nutc_int(11,678)
    equivalence (nutc_int(1,  1),IX01(1,1))
    equivalence (nutc_int(1, 11),IX02(1,1))
    equivalence (nutc_int(1, 21),IX03(1,1))
    equivalence (nutc_int(1, 31),IX04(1,1))
    equivalence (nutc_int(1, 41),IX05(1,1))
    equivalence (nutc_int(1, 51),IX06(1,1))
    equivalence (nutc_int(1, 61),IX07(1,1))
    equivalence (nutc_int(1, 71),IX08(1,1))
    equivalence (nutc_int(1, 81),IX09(1,1))
    equivalence (nutc_int(1, 91),IX10(1,1))
    equivalence (nutc_int(1,101),IX11(1,1))
    equivalence (nutc_int(1,111),IX12(1,1))
    equivalence (nutc_int(1,121),IX13(1,1))
    equivalence (nutc_int(1,131),IX14(1,1))
    equivalence (nutc_int(1,141),IX15(1,1))
    equivalence (nutc_int(1,151),IX16(1,1))
    equivalence (nutc_int(1,161),IX17(1,1))
    equivalence (nutc_int(1,171),IX18(1,1))
    equivalence (nutc_int(1,181),IX19(1,1))
    equivalence (nutc_int(1,191),IX20(1,1))
    equivalence (nutc_int(1,201),IX21(1,1))
    equivalence (nutc_int(1,211),IX22(1,1))
    equivalence (nutc_int(1,221),IX23(1,1))
    equivalence (nutc_int(1,231),IX24(1,1))
    equivalence (nutc_int(1,241),IX25(1,1))
    equivalence (nutc_int(1,251),IX26(1,1))
    equivalence (nutc_int(1,261),IX27(1,1))
    equivalence (nutc_int(1,271),IX28(1,1))
    equivalence (nutc_int(1,281),IX29(1,1))
    equivalence (nutc_int(1,291),IX30(1,1))
    equivalence (nutc_int(1,301),IX31(1,1))
    equivalence (nutc_int(1,311),IX32(1,1))
    equivalence (nutc_int(1,321),IX33(1,1))
    equivalence (nutc_int(1,331),IX34(1,1))
    equivalence (nutc_int(1,341),IX35(1,1))
    equivalence (nutc_int(1,351),IX36(1,1))
    equivalence (nutc_int(1,361),IX37(1,1))
    equivalence (nutc_int(1,371),IX38(1,1))
    equivalence (nutc_int(1,381),IX39(1,1))
    equivalence (nutc_int(1,391),IX40(1,1))
    equivalence (nutc_int(1,401),IX41(1,1))
    equivalence (nutc_int(1,411),IX42(1,1))
    equivalence (nutc_int(1,421),IX43(1,1))
    equivalence (nutc_int(1,431),IX44(1,1))
    equivalence (nutc_int(1,441),IX45(1,1))
    equivalence (nutc_int(1,451),IX46(1,1))
    equivalence (nutc_int(1,461),IX47(1,1))
    equivalence (nutc_int(1,471),IX48(1,1))
    equivalence (nutc_int(1,481),IX49(1,1))
    equivalence (nutc_int(1,491),IX50(1,1))
    equivalence (nutc_int(1,501),IX51(1,1))
    equivalence (nutc_int(1,511),IX52(1,1))
    equivalence (nutc_int(1,521),IX53(1,1))
    equivalence (nutc_int(1,531),IX54(1,1))
    equivalence (nutc_int(1,541),IX55(1,1))
    equivalence (nutc_int(1,551),IX56(1,1))
    equivalence (nutc_int(1,561),IX57(1,1))
    equivalence (nutc_int(1,571),IX58(1,1))
    equivalence (nutc_int(1,581),IX59(1,1))
    equivalence (nutc_int(1,591),IX60(1,1))
    equivalence (nutc_int(1,601),IX61(1,1))
    equivalence (nutc_int(1,611),IX62(1,1))
    equivalence (nutc_int(1,621),IX63(1,1))
    equivalence (nutc_int(1,631),IX64(1,1))
    equivalence (nutc_int(1,641),IX65(1,1))
    equivalence (nutc_int(1,651),IX66(1,1))
    equivalence (nutc_int(1,661),IX67(1,1))
    equivalence (nutc_int(1,671),IX68(1,1))
    data IX01/   0,   0,   0,   0,   1,-172064161,-174666,   92052331,  9086, 33386, 15377, &
         0,   0,   2,  -2,   2, -13170906,   -1675,   5730336, -3015,-13696, -4587, &
         0,   0,   2,   0,   2,  -2276413,    -234,   978459,  -485,  2796,  1374, &
         0,   0,   0,   0,   2,   2074554,     207,   -897492,   470,  -698,  -291, &
         0,   1,   0,   0,   0,   1475877,   -3633,   73871,  -184, 11817, -1924, &
         0,   1,   2,  -2,   2,   -516821,    1226,   224386,  -677,  -524,  -174, &
         1,   0,   0,   0,   0,    711159,      73,   -6750,     0,  -872,   358, &
         0,   0,   2,   0,   1,   -387298,    -367,   200728,    18,   380,   318, &
         1,   0,   2,   0,   2,   -301461,     -36,   129025,   -63,   816,   367, &
         0,  -1,   2,  -2,   2,    215829,    -494,   -95929,   299,   111,   132 /
    data IX02/   0,   0,   2,  -2,   1,    128227,    137,   -68982,    -9,   181,    39, &
         -1,   0,   2,   0,   2,    123457,      11,   -53311,    32,    19,    -4, &
         -1,   0,   0,   2,   0,    156994,      10,   -1235,     0,  -168,    82, &
         1,   0,   0,   0,   1,     63110,      63,   -33228,     0,    27,    -9, &
         -1,   0,   0,   0,   1,    -57976,     -63,   31429,     0,  -189,   -75, &
         -1,   0,   2,   2,   2,    -59641,     -11,   25543,   -11,   149,    66, &
         1,   0,   2,   0,   1,    -51613,     -42,   26366,     0,   129,    78, &
         -2,   0,   2,   0,   1,     45893,      50,   -24236,   -10,    31,    20, &
         0,   0,   0,   2,   0,     63384,      11,   -1220,     0,  -150,    29, &
         0,   0,   2,   2,   2,    -38571,      -1,   16452,   -11,   158,    68 /
    data IX03/   0,  -2,   2,  -2,   2,     32481,      0,   -13870,     0,     0,     0, &
         -2,   0,   0,   2,   0,    -47722,       0,   477,     0,   -18,   -25, &
         2,   0,   2,   0,   2,    -31046,      -1,   13238,   -11,   131,    59, &
         1,   0,   2,  -2,   2,     28593,       0,   -12338,    10,    -1,    -3, &
         -1,   0,   2,   0,   1,     20441,      21,   -10758,     0,    10,    -3, &
         2,   0,   0,   0,   0,     29243,       0,   -609,     0,   -74,    13, &
         0,   0,   2,   0,   0,     25887,       0,   -550,     0,   -66,    11, &
         0,   1,   0,   0,   1,    -14053,     -25,   8551,    -2,    79,   -45, &
         -1,   0,   0,   2,   1,     15164,      10,   -8001,     0,    11,    -1, &
         0,   2,   2,  -2,   2,    -15794,      72,   6850,   -42,   -16,    -5 /
    data IX04/   0,   0,  -2,   2,   0,     21783,      0,   -167,     0,    13,    13, &
         1,   0,   0,  -2,   1,    -12873,     -10,   6953,     0,   -37,   -14, &
         0,  -1,   0,   0,   1,    -12654,      11,   6415,     0,    63,    26, &
         -1,   0,   2,   2,   1,    -10204,       0,   5222,     0,    25,    15, &
         0,   2,   0,   0,   0,     16707,     -85,   168,    -1,   -10,    10, &
         1,   0,   2,   2,   2,     -7691,       0,   3268,     0,    44,    19, &
         -2,   0,   2,   0,   0,    -11024,       0,   104,     0,   -14,     2, &
         0,   1,   2,   0,   2,      7566,     -21,   -3250,     0,   -11,    -5, &
         0,   0,   2,   2,   1,     -6637,     -11,   3353,     0,    25,    14, &
         0,  -1,   2,   0,   2,     -7141,      21,   3070,     0,     8,     4 /
    data IX05/   0,   0,   0,   2,   1,     -6302,    -11,   3272,     0,     2,     4, &
         1,   0,   2,  -2,   1,      5800,      10,   -3045,     0,     2,    -1, &
         2,   0,   2,  -2,   2,      6443,       0,   -2768,     0,    -7,    -4, &
         -2,   0,   0,   2,   1,     -5774,     -11,   3041,     0,   -15,    -5, &
         2,   0,   2,   0,   1,     -5350,       0,   2695,     0,    21,    12, &
         0,  -1,   2,  -2,   1,     -4752,     -11,   2719,     0,    -3,    -3, &
         0,   0,   0,  -2,   1,     -4940,     -11,   2720,     0,   -21,    -9, &
         -1,  -1,   0,   2,   0,      7350,       0,   -51,     0,    -8,     4, &
         2,   0,   0,  -2,   1,      4065,       0,   -2206,     0,     6,     1, &
         1,   0,   0,   2,   0,      6579,       0,   -199,     0,   -24,     2 /
    data IX06/   0,   1,   2,  -2,   1,      3579,      0,   -1900,     0,     5,     1, &
         1,  -1,   0,   0,   0,      4725,       0,   -41,     0,    -6,     3, &
         -2,   0,   2,   0,   2,     -3075,       0,   1313,     0,    -2,    -1, &
         3,   0,   2,   0,   2,     -2904,       0,   1233,     0,    15,     7, &
         0,  -1,   0,   2,   0,      4348,       0,   -81,     0,   -10,     2, &
         1,  -1,   2,   0,   2,     -2878,       0,   1232,     0,     8,     4, &
         0,   0,   0,   1,   0,     -4230,       0,   -20,     0,     5,    -2, &
         -1,  -1,   2,   2,   2,     -2819,       0,   1207,     0,     7,     3, &
         -1,   0,   2,   0,   0,     -4056,       0,   40,     0,     5,    -2, &
         0,  -1,   2,   2,   2,     -2647,       0,   1129,     0,    11,     5 /
    data IX07/  -2,   0,   0,   0,   1,     -2294,      0,   1266,     0,   -10,    -4, &
         1,   1,   2,   0,   2,      2481,       0,   -1062,     0,    -7,    -3, &
         2,   0,   0,   0,   1,      2179,       0,   -1129,     0,    -2,    -2, &
         -1,   1,   0,   1,   0,      3276,       0,   -9,     0,     1,     0, &
         1,   1,   0,   0,   0,     -3389,       0,   35,     0,     5,    -2, &
         1,   0,   2,   0,   0,      3339,       0,   -107,     0,   -13,     1, &
         -1,   0,   2,  -2,   1,     -1987,       0,   1073,     0,    -6,    -2, &
         1,   0,   0,   0,   2,     -1981,       0,   854,     0,     0,     0, &
         -1,   0,   0,   1,   0,      4026,       0,   -553,     0,  -353,  -139, &
         0,   0,   2,   1,   2,      1660,       0,   -710,     0,    -5,    -2 /
    data IX08/  -1,   0,   2,   4,   2,     -1521,      0,   647,     0,     9,     4, &
         -1,   1,   0,   1,   1,      1314,       0,   -700,     0,     0,     0, &
         0,  -2,   2,  -2,   1,     -1283,       0,   672,     0,     0,     0, &
         1,   0,   2,   2,   1,     -1331,       0,   663,     0,     8,     4, &
         -2,   0,   2,   2,   2,      1383,       0,   -594,     0,    -2,    -2, &
         -1,   0,   0,   0,   2,      1405,       0,   -610,     0,     4,     2, &
         1,   1,   2,  -2,   2,      1290,       0,   -556,     0,     0,     0, &
         -2,   0,   2,   4,   2,     -1214,       0,   518,     0,     5,     2, &
         -1,   0,   4,   0,   2,      1146,       0,   -490,     0,    -3,    -1, &
         2,   0,   2,  -2,   1,      1019,       0,   -527,     0,    -1,    -1 /
    data IX09/   2,   0,   2,   2,   2,     -1100,      0,   465,     0,     9,     4, &
         1,   0,   0,   2,   1,      -970,       0,   496,     0,     2,     1, &
         3,   0,   0,   0,   0,      1575,       0,   -50,     0,    -6,     0, &
         3,   0,   2,  -2,   2,       934,       0,   -399,     0,    -3,    -1, &
         0,   0,   4,  -2,   2,       922,       0,   -395,     0,    -1,    -1, &
         0,   1,   2,   0,   1,       815,       0,   -422,     0,    -1,    -1, &
         0,   0,  -2,   2,   1,       834,       0,   -440,     0,     2,     1, &
         0,   0,   2,  -2,   3,      1248,       0,   -170,     0,     0,     1, &
         -1,   0,   0,   4,   0,      1338,       0,   -39,     0,    -5,     0, &
         2,   0,  -2,   0,   1,       716,       0,   -389,     0,    -2,    -1 /
    data IX10/  -2,   0,   0,   4,   0,      1282,      0,   -23,     0,    -3,     1, &
         -1,  -1,   0,   2,   1,       742,       0,   -391,     0,     1,     0, &
         -1,   0,   0,   1,   1,      1020,       0,   -495,     0,   -25,   -10, &
         0,   1,   0,   0,   2,       715,       0,   -326,     0,    -4,     2, &
         0,   0,  -2,   0,   1,      -666,       0,   369,     0,    -3,    -1, &
         0,  -1,   2,   0,   1,      -667,       0,   346,     0,     1,     1, &
         0,   0,   2,  -1,   2,      -704,       0,   304,     0,     0,     0, &
         0,   0,   2,   4,   2,      -694,       0,   294,     0,     5,     2, &
         -2,  -1,   0,   2,   0,     -1014,       0,   4,     0,    -1,    -1, &
         1,   1,   0,  -2,   1,      -585,       0,   316,     0,    -2,    -1 /
    data IX11/  -1,   1,   0,   2,   0,      -949,      0,   8,     0,     1,    -1, &
         -1,   1,   0,   1,   2,      -595,       0,   258,     0,     0,     0, &
         1,  -1,   0,   0,   1,       528,       0,   -279,     0,     0,     0, &
         1,  -1,   2,   2,   2,      -590,       0,   252,     0,     4,     2, &
         -1,   1,   2,   2,   2,       570,       0,   -244,     0,    -2,    -1, &
         3,   0,   2,   0,   1,      -502,       0,   250,     0,     3,     2, &
         0,   1,  -2,   2,   0,      -875,       0,   29,     0,     1,     0, &
         -1,   0,   0,  -2,   1,      -492,       0,   275,     0,    -3,    -1, &
         0,   1,   2,   2,   2,       535,       0,   -228,     0,    -2,    -1, &
         -1,  -1,   2,   2,   1,      -467,       0,   240,     0,     1,     1 /
    data IX12/   0,  -1,   0,   0,   2,       591,      0,   -253,     0,     0,     0, &
         1,   0,   2,  -4,   1,      -453,       0,   244,     0,    -1,    -1, &
         -1,   0,  -2,   2,   0,       766,       0,   9,     0,     1,     0, &
         0,  -1,   2,   2,   1,      -446,       0,   225,     0,     2,     1, &
         2,  -1,   2,   0,   2,      -488,       0,   207,     0,     2,     1, &
         0,   0,   0,   2,   2,      -468,       0,   201,     0,     0,     0, &
         1,  -1,   2,   0,   1,      -421,       0,   216,     0,     1,     1, &
         -1,   1,   2,   0,   2,       463,       0,   -200,     0,     0,     0, &
         0,   1,   0,   2,   0,      -673,       0,   14,     0,     2,     0, &
         0,  -1,  -2,   2,   0,       658,       0,   -2,     0,     0,     0 /
    data IX13/   0,   3,   2,  -2,   2,      -438,      0,   188,     0,     0,     0, &
         0,   0,   0,   1,   1,      -390,       0,   205,     0,     0,     0, &
         -1,   0,   2,   2,   0,       639,     -11,   -19,     0,    -2,     0, &
         2,   1,   2,   0,   2,       412,       0,   -176,     0,    -2,    -1, &
         1,   1,   0,   0,   1,      -361,       0,   189,     0,     0,     0, &
         1,   1,   2,   0,   1,       360,       0,   -185,     0,    -1,    -1, &
         2,   0,   0,   2,   0,       588,       0,   -24,     0,    -3,     0, &
         1,   0,  -2,   2,   0,      -578,       0,   5,     0,     1,     0, &
         -1,   0,   0,   2,   2,      -396,       0,   171,     0,     0,     0, &
         0,   1,   0,   1,   0,       565,       0,   -6,     0,    -1,     0 /
    data IX14/   0,   1,   0,  -2,   1,      -335,      0,   184,     0,    -1,    -1, &
         -1,   0,   2,  -2,   2,       357,       0,   -154,     0,     1,     0, &
         0,   0,   0,  -1,   1,       321,       0,   -174,     0,     1,     0, &
         -1,   1,   0,   0,   1,      -301,       0,   162,     0,    -1,     0, &
         1,   0,   2,  -1,   2,      -334,       0,   144,     0,     0,     0, &
         1,  -1,   0,   2,   0,       493,       0,   -15,     0,    -2,     0, &
         0,   0,   0,   4,   0,       494,       0,   -19,     0,    -2,     0, &
         1,   0,   2,   1,   2,       337,       0,   -143,     0,    -1,    -1, &
         0,   0,   2,   1,   1,       280,       0,   -144,     0,    -1,     0, &
         1,   0,   0,  -2,   2,       309,       0,   -134,     0,     1,     0 /
    data IX15/  -1,   0,   2,   4,   1,      -263,      0,   131,     0,     2,     1, &
         1,   0,  -2,   0,   1,       253,       0,   -138,     0,     1,     0, &
         1,   1,   2,  -2,   1,       245,       0,   -128,     0,     0,     0, &
         0,   0,   2,   2,   0,       416,       0,   -17,     0,    -2,     0, &
         -1,   0,   2,  -1,   1,      -229,       0,   128,     0,     0,     0, &
         -2,   0,   2,   2,   1,       231,       0,   -120,     0,     0,     0, &
         4,   0,   2,   0,   2,      -259,       0,   109,     0,     2,     1, &
         2,  -1,   0,   0,   0,       375,       0,   -8,     0,    -1,     0, &
         2,   1,   2,  -2,   2,       252,       0,   -108,     0,     0,     0, &
         0,   1,   2,   1,   2,      -245,       0,   104,     0,     1,     0 /
    data IX16/   1,   0,   4,  -2,   2,       243,      0,   -104,     0,    -1,     0, &
         -1,  -1,   0,   0,   1,       208,       0,   -112,     0,     1,     0, &
         0,   1,   0,   2,   1,       199,       0,   -102,     0,     0,     0, &
         -2,   0,   2,   4,   1,      -208,       0,   105,     0,     1,     0, &
         2,   0,   2,   0,   0,       335,       0,   -14,     0,    -2,     0, &
         1,   0,   0,   1,   0,      -325,       0,   7,     0,     1,     0, &
         -1,   0,   0,   4,   1,      -187,       0,   96,     0,     0,     0, &
         -1,   0,   4,   0,   1,       197,       0,   -100,     0,    -1,     0, &
         2,   0,   2,   2,   1,      -192,       0,   94,     0,     2,     1, &
         0,   0,   2,  -3,   2,      -188,       0,   83,     0,     0,     0 /
    data IX17/  -1,  -2,   0,   2,   0,       276,      0,   -2,     0,     0,     0, &
         2,   1,   0,   0,   0,      -286,       0,   6,     0,     1,     0, &
         0,   0,   4,   0,   2,       186,       0,   -79,     0,    -1,     0, &
         0,   0,   0,   0,   3,      -219,       0,   43,     0,     0,     0, &
         0,   3,   0,   0,   0,       276,       0,   2,     0,     0,     0, &
         0,   0,   2,  -4,   1,      -153,       0,   84,     0,    -1,     0, &
         0,  -1,   0,   2,   1,      -156,       0,   81,     0,     0,     0, &
         0,   0,   0,   4,   1,      -154,       0,   78,     0,     1,     0, &
         -1,  -1,   2,   4,   2,      -174,       0,   75,     0,     1,     0, &
         1,   0,   2,   4,   2,      -163,       0,   69,     0,     2,     1 /
    data IX18/  -2,   2,   0,   2,   0,      -228,      0,   1,     0,     0,     0, &
         -2,  -1,   2,   0,   1,        91,       0,   -54,     0,    -4,    -2, &
         -2,   0,   0,   2,   2,       175,       0,   -75,     0,     0,     0, &
         -1,  -1,   2,   0,   2,      -159,       0,   69,     0,     0,     0, &
         0,   0,   4,  -2,   1,       141,       0,   -72,     0,     0,     0, &
         3,   0,   2,  -2,   1,       147,       0,   -75,     0,     0,     0, &
         -2,  -1,   0,   2,   1,      -132,       0,   69,     0,     0,     0, &
         1,   0,   0,  -1,   1,       159,       0,   -54,     0,   -28,    11, &
         0,  -2,   0,   2,   0,       213,       0,   -4,     0,     0,     0, &
         -2,   0,   0,   4,   1,       123,       0,   -64,     0,     0,     0 /
    data IX19/  -3,   0,   0,   0,   1,      -118,      0,   66,     0,    -1,     0, &
         1,   1,   2,   2,   2,       144,       0,   -61,     0,    -1,     0, &
         0,   0,   2,   4,   1,      -121,       0,   60,     0,     1,     0, &
         3,   0,   2,   2,   2,      -134,       0,   56,     0,     1,     1, &
         -1,   1,   2,  -2,   1,      -105,       0,   57,     0,     0,     0, &
         2,   0,   0,  -4,   1,      -102,       0,   56,     0,     0,     0, &
         0,   0,   0,  -2,   2,       120,       0,   -52,     0,     0,     0, &
         2,   0,   2,  -4,   1,       101,       0,   -54,     0,     0,     0, &
         -1,   1,   0,   2,   1,      -113,       0,   59,     0,     0,     0, &
         0,   0,   2,  -1,   1,      -106,       0,   61,     0,     0,     0 /
    data IX20/   0,  -2,   2,   2,   2,      -129,      0,   55,     0,     1,     0, &
         2,   0,   0,   2,   1,      -114,       0,   57,     0,     0,     0, &
         4,   0,   2,  -2,   2,       113,       0,   -49,     0,    -1,     0, &
         2,   0,   0,  -2,   2,      -102,       0,   44,     0,     0,     0, &
         0,   2,   0,   0,   1,       -94,       0,   51,     0,     0,     0, &
         1,   0,   0,  -4,   1,      -100,       0,   56,     0,    -1,     0, &
         0,   2,   2,  -2,   1,        87,       0,   -47,     0,     0,     0, &
         -3,   0,   0,   4,   0,       161,       0,   -1,     0,     0,     0, &
         -1,   1,   2,   0,   1,        96,       0,   -50,     0,     0,     0, &
         -1,  -1,   0,   4,   0,       151,       0,   -5,     0,    -1,     0 /
    data IX21/  -1,  -2,   2,   2,   2,      -104,      0,   44,     0,     0,     0, &
         -2,  -1,   2,   4,   2,      -110,       0,   48,     0,     0,     0, &
         1,  -1,   2,   2,   1,      -100,       0,   50,     0,     1,     0, &
         -2,   1,   0,   2,   0,        92,       0,   12,     0,    -5,    -2, &
         -2,   1,   2,   0,   1,        82,       0,   -45,     0,     0,     0, &
         2,   1,   0,  -2,   1,        82,       0,   -45,     0,     0,     0, &
         -3,   0,   2,   0,   1,       -78,       0,   41,     0,     0,     0, &
         -2,   0,   2,  -2,   1,       -77,       0,   43,     0,     0,     0, &
         -1,   1,   0,   2,   2,         2,       0,   54,     0,     0,     0, &
         0,  -1,   2,  -1,   2,        94,       0,   -40,     0,     0,     0 /
    data IX22/  -1,   0,   4,  -2,   2,       -93,      0,   40,     0,     0,     0, &
         0,  -2,   2,   0,   2,       -83,       0,   40,     0,    10,    -2, &
         -1,   0,   2,   1,   2,        83,       0,   -36,     0,     0,     0, &
         2,   0,   0,   0,   2,       -91,       0,   39,     0,     0,     0, &
         0,   0,   2,   0,   3,       128,       0,   -1,     0,     0,     0, &
         -2,   0,   4,   0,   2,       -79,       0,   34,     0,     0,     0, &
         -1,   0,  -2,   0,   1,       -83,       0,   47,     0,     0,     0, &
         -1,   1,   2,   2,   1,        84,       0,   -44,     0,     0,     0, &
         3,   0,   0,   0,   1,        83,       0,   -43,     0,     0,     0, &
         -1,   0,   2,   3,   2,        91,       0,   -39,     0,     0,     0 /
    data IX23/   2,  -1,   2,   0,   1,       -77,      0,   39,     0,     0,     0, &
         0,   1,   2,   2,   1,        84,       0,   -43,     0,     0,     0, &
         0,  -1,   2,   4,   2,       -92,       0,   39,     0,     1,     0, &
         2,  -1,   2,   2,   2,       -92,       0,   39,     0,     1,     0, &
         0,   2,  -2,   2,   0,       -94,       0,   0,     0,     0,     0, &
         -1,  -1,   2,  -1,   1,        68,       0,   -36,     0,     0,     0, &
         0,  -2,   0,   0,   1,       -61,       0,   32,     0,     0,     0, &
         1,   0,   2,  -4,   2,        71,       0,   -31,     0,     0,     0, &
         1,  -1,   0,  -2,   1,        62,       0,   -34,     0,     0,     0, &
         -1,  -1,   2,   0,   1,       -63,       0,   33,     0,     0,     0 /
    data IX24/   1,  -1,   2,  -2,   2,       -73,      0,   32,     0,     0,     0, &
         -2,  -1,   0,   4,   0,       115,       0,   -2,     0,     0,     0, &
         -1,   0,   0,   3,   0,      -103,       0,   2,     0,     0,     0, &
         -2,  -1,   2,   2,   2,        63,       0,   -28,     0,     0,     0, &
         0,   2,   2,   0,   2,        74,       0,   -32,     0,     0,     0, &
         1,   1,   0,   2,   0,      -103,       0,   3,     0,    -3,    -1, &
         2,   0,   2,  -1,   2,       -69,       0,   30,     0,     0,     0, &
         1,   0,   2,   1,   1,        57,       0,   -29,     0,     0,     0, &
         4,   0,   0,   0,   0,        94,       0,   -4,     0,     0,     0, &
         2,   1,   2,   0,   1,        64,       0,   -33,     0,     0,     0 /
    data IX25/   3,  -1,   2,   0,   2,       -63,      0,   26,     0,     0,     0, &
         -2,   2,   0,   2,   1,       -38,       0,   20,     0,     0,     0, &
         1,   0,   2,  -3,   1,       -43,       0,   24,     0,     0,     0, &
         1,   1,   2,  -4,   1,       -45,       0,   23,     0,     0,     0, &
         -1,  -1,   2,  -2,   1,        47,       0,   -24,     0,     0,     0, &
         0,  -1,   0,  -1,   1,       -48,       0,   25,     0,     0,     0, &
         0,  -1,   0,  -2,   1,        45,       0,   -26,     0,     0,     0, &
         -2,   0,   0,   0,   2,        56,       0,   -25,     0,     0,     0, &
         -2,   0,  -2,   2,   0,        88,       0,   2,     0,     0,     0, &
         -1,   0,  -2,   4,   0,       -75,       0,   0,     0,     0,     0 /
    data IX26/   1,  -2,   0,   0,   0,        85,      0,   0,     0,     0,     0, &
         0,   1,   0,   1,   1,        49,       0,   -26,     0,     0,     0, &
         -1,   2,   0,   2,   0,       -74,       0,   -1,     0,    -3,    -1, &
         1,  -1,   2,  -2,   1,       -39,       0,   21,     0,     0,     0, &
         1,   2,   2,  -2,   2,        45,       0,   -20,     0,     0,     0, &
         2,  -1,   2,  -2,   2,        51,       0,   -22,     0,     0,     0, &
         1,   0,   2,  -1,   1,       -40,       0,   21,     0,     0,     0, &
         2,   1,   2,  -2,   1,        41,       0,   -21,     0,     0,     0, &
         -2,   0,   0,  -2,   1,       -42,       0,   24,     0,     0,     0, &
         1,  -2,   2,   0,   2,       -51,       0,   22,     0,     0,     0 /
    data IX27/   0,   1,   2,   1,   1,       -42,      0,   22,     0,     0,     0, &
         1,   0,   4,  -2,   1,        39,       0,   -21,     0,     0,     0, &
         -2,   0,   4,   2,   2,        46,       0,   -18,     0,     0,     0, &
         1,   1,   2,   1,   2,       -53,       0,   22,     0,     0,     0, &
         1,   0,   0,   4,   0,        82,       0,   -4,     0,     0,     0, &
         1,   0,   2,   2,   0,        81,       0,   -4,     0,    -1,     0, &
         2,   0,   2,   1,   2,        47,       0,   -19,     0,     0,     0, &
         3,   1,   2,   0,   2,        53,       0,   -23,     0,     0,     0, &
         4,   0,   2,   0,   1,       -45,       0,   22,     0,     0,     0, &
         -2,  -1,   2,   0,   0,       -44,       0,   -2,     0,     0,     0 /
    data IX28/   0,   1,  -2,   2,   1,       -33,      0,   16,     0,     0,     0, &
         1,   0,  -2,   1,   0,       -61,       0,   1,     0,     0,     0, &
         0,  -1,  -2,   2,   1,        28,       0,   -15,     0,     0,     0, &
         2,  -1,   0,  -2,   1,       -38,       0,   19,     0,     0,     0, &
         -1,   0,   2,  -1,   2,       -33,       0,   21,     0,     0,     0, &
         1,   0,   2,  -3,   2,       -60,       0,   0,     0,     0,     0, &
         0,   1,   2,  -2,   3,        48,       0,   -10,     0,     0,     0, &
         0,   0,   2,  -3,   1,        27,       0,   -14,     0,     0,     0, &
         -1,   0,  -2,   2,   1,        38,       0,   -20,     0,     0,     0, &
         0,   0,   2,  -4,   2,        31,       0,   -13,     0,     0,     0 /
    data IX29/  -2,   1,   0,   0,   1,       -29,      0,   15,     0,     0,     0, &
         -1,   0,   0,  -1,   1,        28,       0,   -15,     0,     0,     0, &
         2,   0,   2,  -4,   2,       -32,       0,   15,     0,     0,     0, &
         0,   0,   4,  -4,   4,        45,       0,   -8,     0,     0,     0, &
         0,   0,   4,  -4,   2,       -44,       0,   19,     0,     0,     0, &
         -1,  -2,   0,   2,   1,        28,       0,   -15,     0,     0,     0, &
         -2,   0,   0,   3,   0,       -51,       0,   0,     0,     0,     0, &
         1,   0,  -2,   2,   1,       -36,       0,   20,     0,     0,     0, &
         -3,   0,   2,   2,   2,        44,       0,   -19,     0,     0,     0, &
         -3,   0,   2,   2,   1,        26,       0,   -14,     0,     0,     0 /
    data IX30/  -2,   0,   2,   2,   0,       -60,      0,   2,     0,     0,     0, &
         2,  -1,   0,   0,   1,        35,       0,   -18,     0,     0,     0, &
         -2,   1,   2,   2,   2,       -27,       0,   11,     0,     0,     0, &
         1,   1,   0,   1,   0,        47,       0,   -1,     0,     0,     0, &
         0,   1,   4,  -2,   2,        36,       0,   -15,     0,     0,     0, &
         -1,   1,   0,  -2,   1,       -36,       0,   20,     0,     0,     0, &
         0,   0,   0,  -4,   1,       -35,       0,   19,     0,     0,     0, &
         1,  -1,   0,   2,   1,       -37,       0,   19,     0,     0,     0, &
         1,   1,   0,   2,   1,        32,       0,   -16,     0,     0,     0, &
         -1,   2,   2,   2,   2,        35,       0,   -14,     0,     0,     0 /
    data IX31/   3,   1,   2,  -2,   2,        32,      0,   -13,     0,     0,     0, &
         0,  -1,   0,   4,   0,        65,       0,   -2,     0,     0,     0, &
         2,  -1,   0,   2,   0,        47,       0,   -1,     0,     0,     0, &
         0,   0,   4,   0,   1,        32,       0,   -16,     0,     0,     0, &
         2,   0,   4,  -2,   2,        37,       0,   -16,     0,     0,     0, &
         -1,  -1,   2,   4,   1,       -30,       0,   15,     0,     0,     0, &
         1,   0,   0,   4,   1,       -32,       0,   16,     0,     0,     0, &
         1,  -2,   2,   2,   2,       -31,       0,   13,     0,     0,     0, &
         0,   0,   2,   3,   2,        37,       0,   -16,     0,     0,     0, &
         -1,   1,   2,   4,   2,        31,       0,   -13,     0,     0,     0 /
    data IX32/   3,   0,   0,   2,   0,        49,      0,   -2,     0,     0,     0, &
         -1,   0,   4,   2,   2,        32,       0,   -13,     0,     0,     0, &
         1,   1,   2,   2,   1,        23,       0,   -12,     0,     0,     0, &
         -2,   0,   2,   6,   2,       -43,       0,   18,     0,     0,     0, &
         2,   1,   2,   2,   2,        26,       0,   -11,     0,     0,     0, &
         -1,   0,   2,   6,   2,       -32,       0,   14,     0,     0,     0, &
         1,   0,   2,   4,   1,       -29,       0,   14,     0,     0,     0, &
         2,   0,   2,   4,   2,       -27,       0,   12,     0,     0,     0, &
         1,   1,  -2,   1,   0,        30,       0,   0,     0,     0,     0, &
         -3,   1,   2,   1,   2,       -11,       0,   5,     0,     0,     0 /
    data IX33/   2,   0,  -2,   0,   2,       -21,      0,   10,     0,     0,     0, &
         -1,   0,   0,   1,   2,       -34,       0,   15,     0,     0,     0, &
         -4,   0,   2,   2,   1,       -10,       0,   6,     0,     0,     0, &
         -1,  -1,   0,   1,   0,       -36,       0,   0,     0,     0,     0, &
         0,   0,  -2,   2,   2,        -9,       0,   4,     0,     0,     0, &
         1,   0,   0,  -1,   2,       -12,       0,   5,     0,     0,     0, &
         0,  -1,   2,  -2,   3,       -21,       0,   5,     0,     0,     0, &
         -2,   1,   2,   0,   0,       -29,       0,   -1,     0,     0,     0, &
         0,   0,   2,  -2,   4,       -15,       0,   3,     0,     0,     0, &
         -2,  -2,   0,   2,   0,       -20,       0,   0,     0,     0,     0 /
    data IX34/  -2,   0,  -2,   4,   0,        28,      0,   0,     0,     0,    -2, &
         0,  -2,  -2,   2,   0,        17,       0,   0,     0,     0,     0, &
         1,   2,   0,  -2,   1,       -22,       0,   12,     0,     0,     0, &
         3,   0,   0,  -4,   1,       -14,       0,   7,     0,     0,     0, &
         -1,   1,   2,  -2,   2,        24,       0,   -11,     0,     0,     0, &
         1,  -1,   2,  -4,   1,        11,       0,   -6,     0,     0,     0, &
         1,   1,   0,  -2,   2,        14,       0,   -6,     0,     0,     0, &
         -3,   0,   2,   0,   0,        24,       0,   0,     0,     0,     0, &
         -3,   0,   2,   0,   2,        18,       0,   -8,     0,     0,     0, &
         -2,   0,   0,   1,   0,       -38,       0,   0,     0,     0,     0 /
    data IX35/   0,   0,  -2,   1,   0,       -31,      0,   0,     0,     0,     0, &
         -3,   0,   0,   2,   1,       -16,       0,   8,     0,     0,     0, &
         -1,  -1,  -2,   2,   0,        29,       0,   0,     0,     0,     0, &
         0,   1,   2,  -4,   1,       -18,       0,   10,     0,     0,     0, &
         2,   1,   0,  -4,   1,       -10,       0,   5,     0,     0,     0, &
         0,   2,   0,  -2,   1,       -17,       0,   10,     0,     0,     0, &
         1,   0,   0,  -3,   1,         9,       0,   -4,     0,     0,     0, &
         -2,   0,   2,  -2,   2,        16,       0,   -6,     0,     0,     0, &
         -2,  -1,   0,   0,   1,        22,       0,   -12,     0,     0,     0, &
         -4,   0,   0,   2,   0,        20,       0,   0,     0,     0,     0 /
    data IX36/   1,   1,   0,  -4,   1,       -13,      0,   6,     0,     0,     0, &
         -1,   0,   2,  -4,   1,       -17,       0,   9,     0,     0,     0, &
         0,   0,   4,  -4,   1,       -14,       0,   8,     0,     0,     0, &
         0,   3,   2,  -2,   2,         0,       0,   -7,     0,     0,     0, &
         -3,  -1,   0,   4,   0,        14,       0,   0,     0,     0,     0, &
         -3,   0,   0,   4,   1,        19,       0,   -10,     0,     0,     0, &
         1,  -1,  -2,   2,   0,       -34,       0,   0,     0,     0,     0, &
         -1,  -1,   0,   2,   2,       -20,       0,   8,     0,     0,     0, &
         1,  -2,   0,   0,   1,         9,       0,   -5,     0,     0,     0, &
         1,  -1,   0,   0,   2,       -18,       0,   7,     0,     0,     0 /
    data IX37/   0,   0,   0,   1,   2,        13,      0,   -6,     0,     0,     0, &
         -1,  -1,   2,   0,   0,        17,       0,   0,     0,     0,     0, &
         1,  -2,   2,  -2,   2,       -12,       0,   5,     0,     0,     0, &
         0,  -1,   2,  -1,   1,        15,       0,   -8,     0,     0,     0, &
         -1,   0,   2,   0,   3,       -11,       0,   3,     0,     0,     0, &
         1,   1,   0,   0,   2,        13,       0,   -5,     0,     0,     0, &
         -1,   1,   2,   0,   0,       -18,       0,   0,     0,     0,     0, &
         1,   2,   0,   0,   0,       -35,       0,   0,     0,     0,     0, &
         -1,   2,   2,   0,   2,         9,       0,   -4,     0,     0,     0, &
         -1,   0,   4,  -2,   1,       -19,       0,   10,     0,     0,     0 /
    data IX38/   3,   0,   2,  -4,   2,       -26,      0,   11,     0,     0,     0, &
         1,   2,   2,  -2,   1,         8,       0,   -4,     0,     0,     0, &
         1,   0,   4,  -4,   2,       -10,       0,   4,     0,     0,     0, &
         -2,  -1,   0,   4,   1,        10,       0,   -6,     0,     0,     0, &
         0,  -1,   0,   2,   2,       -21,       0,   9,     0,     0,     0, &
         -2,   1,   0,   4,   0,       -15,       0,   0,     0,     0,     0, &
         -2,  -1,   2,   2,   1,         9,       0,   -5,     0,     0,     0, &
         2,   0,  -2,   2,   0,       -29,       0,   0,     0,     0,     0, &
         1,   0,   0,   1,   1,       -19,       0,   10,     0,     0,     0, &
         0,   1,   0,   2,   2,        12,       0,   -5,     0,     0,     0 /
    data IX39/   1,  -1,   2,  -1,   2,        22,      0,   -9,     0,     0,     0, &
         -2,   0,   4,   0,   1,       -10,       0,   5,     0,     0,     0, &
         2,   1,   0,   0,   1,       -20,       0,   11,     0,     0,     0, &
         0,   1,   2,   0,   0,       -20,       0,   0,     0,     0,     0, &
         0,  -1,   4,  -2,   2,       -17,       0,   7,     0,     0,     0, &
         0,   0,   4,  -2,   4,        15,       0,   -3,     0,     0,     0, &
         0,   2,   2,   0,   1,         8,       0,   -4,     0,     0,     0, &
         -3,   0,   0,   6,   0,        14,       0,   0,     0,     0,     0, &
         -1,  -1,   0,   4,   1,       -12,       0,   6,     0,     0,     0, &
         1,  -2,   0,   2,   0,        25,       0,   0,     0,     0,     0 /
    data IX40/  -1,   0,   0,   4,   2,       -13,      0,   6,     0,     0,     0, &
         -1,  -2,   2,   2,   1,       -14,       0,   8,     0,     0,     0, &
         -1,   0,   0,  -2,   2,        13,       0,   -5,     0,     0,     0, &
         1,   0,  -2,  -2,   1,       -17,       0,   9,     0,     0,     0, &
         0,   0,  -2,  -2,   1,       -12,       0,   6,     0,     0,     0, &
         -2,   0,  -2,   0,   1,       -10,       0,   5,     0,     0,     0, &
         0,   0,   0,   3,   1,        10,       0,   -6,     0,     0,     0, &
         0,   0,   0,   3,   0,       -15,       0,   0,     0,     0,     0, &
         -1,   1,   0,   4,   0,       -22,       0,   0,     0,     0,     0, &
         -1,  -1,   2,   2,   0,        28,       0,   -1,     0,     0,     0 /
    data IX41/  -2,   0,   2,   3,   2,        15,      0,   -7,     0,     0,     0, &
         1,   0,   0,   2,   2,        23,       0,   -10,     0,     0,     0, &
         0,  -1,   2,   1,   2,        12,       0,   -5,     0,     0,     0, &
         3,  -1,   0,   0,   0,        29,       0,   -1,     0,     0,     0, &
         2,   0,   0,   1,   0,       -25,       0,   1,     0,     0,     0, &
         1,  -1,   2,   0,   0,        22,       0,   0,     0,     0,     0, &
         0,   0,   2,   1,   0,       -18,       0,   0,     0,     0,     0, &
         1,   0,   2,   0,   3,        15,       0,   3,     0,     0,     0, &
         3,   1,   0,   0,   0,       -23,       0,   0,     0,     0,     0, &
         3,  -1,   2,  -2,   2,        12,       0,   -5,     0,     0,     0 /
    data IX42/   2,   0,   2,  -1,   1,        -8,      0,   4,     0,     0,     0, &
         1,   1,   2,   0,   0,       -19,       0,   0,     0,     0,     0, &
         0,   0,   4,  -1,   2,       -10,       0,   4,     0,     0,     0, &
         1,   2,   2,   0,   2,        21,       0,   -9,     0,     0,     0, &
         -2,   0,   0,   6,   0,        23,       0,   -1,     0,     0,     0, &
         0,  -1,   0,   4,   1,       -16,       0,   8,     0,     0,     0, &
         -2,  -1,   2,   4,   1,       -19,       0,   9,     0,     0,     0, &
         0,  -2,   2,   2,   1,       -22,       0,   10,     0,     0,     0, &
         0,  -1,   2,   2,   0,        27,       0,   -1,     0,     0,     0, &
         -1,   0,   2,   3,   1,        16,       0,   -8,     0,     0,     0 /
    data IX43/  -2,   1,   2,   4,   2,        19,      0,   -8,     0,     0,     0, &
         2,   0,   0,   2,   2,         9,       0,   -4,     0,     0,     0, &
         2,  -2,   2,   0,   2,        -9,       0,   4,     0,     0,     0, &
         -1,   1,   2,   3,   2,        -9,       0,   4,     0,     0,     0, &
         3,   0,   2,  -1,   2,        -8,       0,   4,     0,     0,     0, &
         4,   0,   2,  -2,   1,        18,       0,   -9,     0,     0,     0, &
         -1,   0,   0,   6,   0,        16,       0,   -1,     0,     0,     0, &
         -1,  -2,   2,   4,   2,       -10,       0,   4,     0,     0,     0, &
         -3,   0,   2,   6,   2,       -23,       0,   9,     0,     0,     0, &
         -1,   0,   2,   4,   0,        16,       0,   -1,     0,     0,     0 /
    data IX44/   3,   0,   0,   2,   1,       -12,      0,   6,     0,     0,     0, &
         3,  -1,   2,   0,   1,        -8,       0,   4,     0,     0,     0, &
         3,   0,   2,   0,   0,        30,       0,   -2,     0,     0,     0, &
         1,   0,   4,   0,   2,        24,       0,   -10,     0,     0,     0, &
         5,   0,   2,  -2,   2,        10,       0,   -4,     0,     0,     0, &
         0,  -1,   2,   4,   1,       -16,       0,   7,     0,     0,     0, &
         2,  -1,   2,   2,   1,       -16,       0,   7,     0,     0,     0, &
         0,   1,   2,   4,   2,        17,       0,   -7,     0,     0,     0, &
         1,  -1,   2,   4,   2,       -24,       0,   10,     0,     0,     0, &
         3,  -1,   2,   2,   2,       -12,       0,   5,     0,     0,     0 /
    data IX45/   3,   0,   2,   2,   1,       -24,      0,   11,     0,     0,     0, &
         5,   0,   2,   0,   2,       -23,       0,   9,     0,     0,     0, &
         0,   0,   2,   6,   2,       -13,       0,   5,     0,     0,     0, &
         4,   0,   2,   2,   2,       -15,       0,   7,     0,     0,     0, &
         0,  -1,   1,  -1,   1,         0,       0,   0,     0, -1988, -1679, &
         -1,   0,   1,   0,   3,         0,       0,   0,     0,   -63,   -27, &
         0,  -2,   2,  -2,   3,        -4,       0,   0,     0,     0,     0, &
         1,   0,  -1,   0,   1,         0,       0,   0,     0,     5,     4, &
         2,  -2,   0,  -2,   1,         5,       0,   -3,     0,     0,     0, &
         -1,   0,   1,   0,   2,         0,       0,   0,     0,   364,   176 /
    data IX46/  -1,   0,   1,   0,   1,         0,      0,   0,     0, -1044,  -891, &
         -1,  -1,   2,  -1,   2,        -3,       0,   1,     0,     0,     0, &
         -2,   2,   0,   2,   2,         4,       0,   -2,     0,     0,     0, &
         -1,   0,   1,   0,   0,         0,       0,   0,     0,   330,     0, &
         -4,   1,   2,   2,   2,         5,       0,   -2,     0,     0,     0, &
         -3,   0,   2,   1,   1,         3,       0,   -2,     0,     0,     0, &
         -2,  -1,   2,   0,   2,        -3,       0,   1,     0,     0,     0, &
         1,   0,  -2,   1,   1,        -5,       0,   2,     0,     0,     0, &
         2,  -1,  -2,   0,   1,         3,       0,   -1,     0,     0,     0, &
         -4,   0,   2,   2,   0,         3,       0,   0,     0,     0,     0 /
    data IX47/  -3,   1,   0,   3,   0,         3,      0,   0,     0,     0,     0, &
         -1,   0,  -1,   2,   0,         0,       0,   0,     0,     5,     0, &
         0,  -2,   0,   0,   2,         0,       0,   1,     0,     0,     0, &
         0,  -2,   0,   0,   2,         4,       0,   -2,     0,     0,     0, &
         -3,   0,   0,   3,   0,         6,       0,   0,     0,     0,     0, &
         -2,  -1,   0,   2,   2,         5,       0,   -2,     0,     0,     0, &
         -1,   0,  -2,   3,   0,        -7,       0,   0,     0,     0,     0, &
         -4,   0,   0,   4,   0,       -12,       0,   0,     0,     0,     0, &
         2,   1,  -2,   0,   1,         5,       0,   -3,     0,     0,     0, &
         2,  -1,   0,  -2,   2,         3,       0,   -1,     0,     0,     0 /
    data IX48/   0,   0,   1,  -1,   0,        -5,      0,   0,     0,     0,     0, &
         -1,   2,   0,   1,   0,         3,       0,   0,     0,     0,     0, &
         -2,   1,   2,   0,   2,        -7,       0,   3,     0,     0,     0, &
         1,   1,   0,  -1,   1,         7,       0,   -4,     0,     0,     0, &
         1,   0,   1,  -2,   1,         0,       0,   0,     0,   -12,   -10, &
         0,   2,   0,   0,   2,         4,       0,   -2,     0,     0,     0, &
         1,  -1,   2,  -3,   1,         3,       0,   -2,     0,     0,     0, &
         -1,   1,   2,  -1,   1,        -3,       0,   2,     0,     0,     0, &
         -2,   0,   4,  -2,   2,        -7,       0,   3,     0,     0,     0, &
         -2,   0,   4,  -2,   1,        -4,       0,   2,     0,     0,     0 /
    data IX49/  -2,  -2,   0,   2,   1,        -3,      0,   1,     0,     0,     0, &
         -2,   0,  -2,   4,   0,         0,       0,   0,     0,     0,     0, &
         1,   2,   2,  -4,   1,        -3,       0,   1,     0,     0,     0, &
         1,   1,   2,  -4,   2,         7,       0,   -3,     0,     0,     0, &
         -1,   2,   2,  -2,   1,        -4,       0,   2,     0,     0,     0, &
         2,   0,   0,  -3,   1,         4,       0,   -2,     0,     0,     0, &
         -1,   2,   0,   0,   1,        -5,       0,   3,     0,     0,     0, &
         0,   0,   0,  -2,   0,         5,       0,   0,     0,     0,     0, &
         -1,  -1,   2,  -2,   2,        -5,       0,   2,     0,     0,     0, &
         -1,   1,   0,   0,   2,         5,       0,   -2,     0,     0,     0 /
    data IX50/   0,   0,   0,  -1,   2,        -8,      0,   3,     0,     0,     0, &
         -2,   1,   0,   1,   0,         9,       0,   0,     0,     0,     0, &
         1,  -2,   0,  -2,   1,         6,       0,   -3,     0,     0,     0, &
         1,   0,  -2,   0,   2,        -5,       0,   2,     0,     0,     0, &
         -3,   1,   0,   2,   0,         3,       0,   0,     0,     0,     0, &
         -1,   1,  -2,   2,   0,        -7,       0,   0,     0,     0,     0, &
         -1,  -1,   0,   0,   2,        -3,       0,   1,     0,     0,     0, &
         -3,   0,   0,   2,   0,         5,       0,   0,     0,     0,     0, &
         -3,  -1,   0,   2,   0,         3,       0,   0,     0,     0,     0, &
         2,   0,   2,  -6,   1,        -3,       0,   2,     0,     0,     0 /
    data IX51/   0,   1,   2,  -4,   2,         4,      0,   -2,     0,     0,     0, &
         2,   0,   0,  -4,   2,         3,       0,   -1,     0,     0,     0, &
         -2,   1,   2,  -2,   1,        -5,       0,   2,     0,     0,     0, &
         0,  -1,   2,  -4,   1,         4,       0,   -2,     0,     0,     0, &
         0,   1,   0,  -2,   2,         9,       0,   -3,     0,     0,     0, &
         -1,   0,   0,  -2,   0,         4,       0,   0,     0,     0,     0, &
         2,   0,  -2,  -2,   1,         4,       0,   -2,     0,     0,     0, &
         -4,   0,   2,   0,   1,        -3,       0,   2,     0,     0,     0, &
         -1,  -1,   0,  -1,   1,        -4,       0,   2,     0,     0,     0, &
         0,   0,  -2,   0,   2,         9,       0,   -3,     0,     0,     0 /
    data IX52/  -3,   0,   0,   1,   0,        -4,      0,   0,     0,     0,     0, &
         -1,   0,  -2,   1,   0,        -4,       0,   0,     0,     0,     0, &
         -2,   0,  -2,   2,   1,         3,       0,   -2,     0,     0,     0, &
         0,   0,  -4,   2,   0,         8,       0,   0,     0,     0,     0, &
         -2,  -1,  -2,   2,   0,         3,       0,   0,     0,     0,     0, &
         1,   0,   2,  -6,   1,        -3,       0,   2,     0,     0,     0, &
         -1,   0,   2,  -4,   2,         3,       0,   -1,     0,     0,     0, &
         1,   0,   0,  -4,   2,         3,       0,   -1,     0,     0,     0, &
         2,   1,   2,  -4,   2,        -3,       0,   1,     0,     0,     0, &
         2,   1,   2,  -4,   1,         6,       0,   -3,     0,     0,     0 /
    data IX53/   0,   1,   4,  -4,   4,         3,      0,   0,     0,     0,     0, &
         0,   1,   4,  -4,   2,        -3,       0,   1,     0,     0,     0, &
         -1,  -1,  -2,   4,   0,        -7,       0,   0,     0,     0,     0, &
         -1,  -3,   0,   2,   0,         9,       0,   0,     0,     0,     0, &
         -1,   0,  -2,   4,   1,        -3,       0,   2,     0,     0,     0, &
         -2,  -1,   0,   3,   0,        -3,       0,   0,     0,     0,     0, &
         0,   0,  -2,   3,   0,        -4,       0,   0,     0,     0,     0, &
         -2,   0,   0,   3,   1,        -5,       0,   3,     0,     0,     0, &
         0,  -1,   0,   1,   0,       -13,       0,   0,     0,     0,     0, &
         -3,   0,   2,   2,   0,        -7,       0,   0,     0,     0,     0 /
    data IX54/   1,   1,  -2,   2,   0,        10,      0,   0,     0,     0,     0, &
         -1,   1,   0,   2,   2,         3,       0,   -1,     0,     0,     0, &
         1,  -2,   2,  -2,   1,        10,       0,   6,     0,    13,    -5, &
         0,   0,   1,   0,   2,         0,       0,   0,     0,    30,    14, &
         0,   0,   1,   0,   1,         0,       0,   0,     0,  -162,  -138, &
         0,   0,   1,   0,   0,         0,       0,   0,     0,    75,     0, &
         -1,   2,   0,   2,   1,        -7,       0,   4,     0,     0,     0, &
         0,   0,   2,   0,   2,        -4,       0,   2,     0,     0,     0, &
         -2,   0,   2,   0,   2,         4,       0,   -2,     0,     0,     0, &
         2,   0,   0,  -1,   1,         5,       0,   -2,     0,     0,     0 /
    data IX55/   3,   0,   0,  -2,   1,         5,      0,   -3,     0,     0,     0, &
         1,   0,   2,  -2,   3,        -3,       0,   0,     0,     0,     0, &
         1,   2,   0,   0,   1,        -3,       0,   2,     0,     0,     0, &
         2,   0,   2,  -3,   2,        -4,       0,   2,     0,     0,     0, &
         -1,   1,   4,  -2,   2,        -5,       0,   2,     0,     0,     0, &
         -2,  -2,   0,   4,   0,         6,       0,   0,     0,     0,     0, &
         0,  -3,   0,   2,   0,         9,       0,   0,     0,     0,     0, &
         0,   0,  -2,   4,   0,         5,       0,   0,     0,     0,     0, &
         -1,  -1,   0,   3,   0,        -7,       0,   0,     0,     0,     0, &
         -2,   0,   0,   4,   2,        -3,       0,   1,     0,     0,     0 /
    data IX56/  -1,   0,   0,   3,   1,        -4,      0,   2,     0,     0,     0, &
         2,  -2,   0,   0,   0,         7,       0,   0,     0,     0,     0, &
         1,  -1,   0,   1,   0,        -4,       0,   0,     0,     0,     0, &
         -1,   0,   0,   2,   0,         4,       0,   0,     0,     0,     0, &
         0,  -2,   2,   0,   1,        -6,       0,   3,     0,    -3,     1, &
         -1,   0,   1,   2,   1,         0,       0,   0,     0,    -3,    -2, &
         -1,   1,   0,   3,   0,        11,       0,   0,     0,     0,     0, &
         -1,  -1,   2,   1,   2,         3,       0,   -1,     0,     0,     0, &
         0,  -1,   2,   0,   0,        11,       0,   0,     0,     0,     0, &
         -2,   1,   2,   2,   1,        -3,       0,   2,     0,     0,     0 /
    data IX57/   2,  -2,   2,  -2,   2,        -1,      0,   3,     0,     3,    -1, &
         1,   1,   0,   1,   1,         4,       0,   -2,     0,     0,     0, &
         1,   0,   1,   0,   1,         0,       0,   0,     0,   -13,   -11, &
         1,   0,   1,   0,   0,         3,       0,   0,     0,     6,     0, &
         0,   2,   0,   2,   0,        -7,       0,   0,     0,     0,     0, &
         2,  -1,   2,  -2,   1,         5,       0,   -3,     0,     0,     0, &
         0,  -1,   4,  -2,   1,        -3,       0,   1,     0,     0,     0, &
         0,   0,   4,  -2,   3,         3,       0,   0,     0,     0,     0, &
         0,   1,   4,  -2,   1,         5,       0,   -3,     0,     0,     0, &
         4,   0,   2,  -4,   2,        -7,       0,   3,     0,     0,     0 /
    data IX58/   2,   2,   2,  -2,   2,         8,      0,   -3,     0,     0,     0, &
         2,   0,   4,  -4,   2,        -4,       0,   2,     0,     0,     0, &
         -1,  -2,   0,   4,   0,        11,       0,   0,     0,     0,     0, &
         -1,  -3,   2,   2,   2,        -3,       0,   1,     0,     0,     0, &
         -3,   0,   2,   4,   2,         3,       0,   -1,     0,     0,     0, &
         -3,   0,   2,  -2,   1,        -4,       0,   2,     0,     0,     0, &
         -1,  -1,   0,  -2,   1,         8,       0,   -4,     0,     0,     0, &
         -3,   0,   0,   0,   2,         3,       0,   -1,     0,     0,     0, &
         -3,   0,  -2,   2,   0,        11,       0,   0,     0,     0,     0, &
         0,   1,   0,  -4,   1,        -6,       0,   3,     0,     0,     0 /
    data IX59/  -2,   1,   0,  -2,   1,        -4,      0,   2,     0,     0,     0, &
         -4,   0,   0,   0,   1,        -8,       0,   4,     0,     0,     0, &
         -1,   0,   0,  -4,   1,        -7,       0,   3,     0,     0,     0, &
         -3,   0,   0,  -2,   1,        -4,       0,   2,     0,     0,     0, &
         0,   0,   0,   3,   2,         3,       0,   -1,     0,     0,     0, &
         -1,   1,   0,   4,   1,         6,       0,   -3,     0,     0,     0, &
         1,  -2,   2,   0,   1,        -6,       0,   3,     0,     0,     0, &
         0,   1,   0,   3,   0,         6,       0,   0,     0,     0,     0, &
         -1,   0,   2,   2,   3,         6,       0,   -1,     0,     0,     0, &
         0,   0,   2,   2,   2,         5,       0,   -2,     0,     0,     0 /
    data IX60/  -2,   0,   2,   2,   2,        -5,      0,   2,     0,     0,     0, &
         -1,   1,   2,   2,   0,        -4,       0,   0,     0,     0,     0, &
         3,   0,   0,   0,   2,        -4,       0,   2,     0,     0,     0, &
         2,   1,   0,   1,   0,         4,       0,   0,     0,     0,     0, &
         2,  -1,   2,  -1,   2,         6,       0,   -3,     0,     0,     0, &
         0,   0,   2,   0,   1,        -4,       0,   2,     0,     0,     0, &
         0,   0,   3,   0,   3,         0,       0,   0,     0,   -26,   -11, &
         0,   0,   3,   0,   2,         0,       0,   0,     0,   -10,    -5, &
         -1,   2,   2,   2,   1,         5,       0,   -3,     0,     0,     0, &
         -1,   0,   4,   0,   0,       -13,       0,   0,     0,     0,     0 /
    data IX61/   1,   2,   2,   0,   1,         3,      0,   -2,     0,     0,     0, &
         3,   1,   2,  -2,   1,         4,       0,   -2,     0,     0,     0, &
         1,   1,   4,  -2,   2,         7,       0,   -3,     0,     0,     0, &
         -2,  -1,   0,   6,   0,         4,       0,   0,     0,     0,     0, &
         0,  -2,   0,   4,   0,         5,       0,   0,     0,     0,     0, &
         -2,   0,   0,   6,   1,        -3,       0,   2,     0,     0,     0, &
         -2,  -2,   2,   4,   2,        -6,       0,   2,     0,     0,     0, &
         0,  -3,   2,   2,   2,        -5,       0,   2,     0,     0,     0, &
         0,   0,   0,   4,   2,        -7,       0,   3,     0,     0,     0, &
         -1,  -1,   2,   3,   2,         5,       0,   -2,     0,     0,     0 /
    data IX62/  -2,   0,   2,   4,   0,        13,      0,   0,     0,     0,     0, &
         2,  -1,   0,   2,   1,        -4,       0,   2,     0,     0,     0, &
         1,   0,   0,   3,   0,        -3,       0,   0,     0,     0,     0, &
         0,   1,   0,   4,   1,         5,       0,   -2,     0,     0,     0, &
         0,   1,   0,   4,   0,       -11,       0,   0,     0,     0,     0, &
         1,  -1,   2,   1,   2,         5,       0,   -2,     0,     0,     0, &
         0,   0,   2,   2,   3,         4,       0,   0,     0,     0,     0, &
         1,   0,   2,   2,   2,         4,       0,   -2,     0,     0,     0, &
         -1,   0,   2,   2,   2,        -4,       0,   2,     0,     0,     0, &
         -2,   0,   4,   2,   1,         6,       0,   -3,     0,     0,     0 /
    data IX63/   2,   1,   0,   2,   1,         3,      0,   -2,     0,     0,     0, &
         2,   1,   0,   2,   0,       -12,       0,   0,     0,     0,     0, &
         2,  -1,   2,   0,   0,         4,       0,   0,     0,     0,     0, &
         1,   0,   2,   1,   0,        -3,       0,   0,     0,     0,     0, &
         0,   1,   2,   2,   0,        -4,       0,   0,     0,     0,     0, &
         2,   0,   2,   0,   3,         3,       0,   0,     0,     0,     0, &
         3,   0,   2,   0,   2,         3,       0,   -1,     0,     0,     0, &
         1,   0,   2,   0,   2,        -3,       0,   1,     0,     0,     0, &
         1,   0,   3,   0,   3,         0,       0,   0,     0,    -5,    -2, &
         1,   1,   2,   1,   1,        -7,       0,   4,     0,     0,     0 /
    data IX64/   0,   2,   2,   2,   2,         6,      0,   -3,     0,     0,     0, &
         2,   1,   2,   0,   0,        -3,       0,   0,     0,     0,     0, &
         2,   0,   4,  -2,   1,         5,       0,   -3,     0,     0,     0, &
         4,   1,   2,  -2,   2,         3,       0,   -1,     0,     0,     0, &
         -1,  -1,   0,   6,   0,         3,       0,   0,     0,     0,     0, &
         -3,  -1,   2,   6,   2,        -3,       0,   1,     0,     0,     0, &
         -1,   0,   0,   6,   1,        -5,       0,   3,     0,     0,     0, &
         -3,   0,   2,   6,   1,        -3,       0,   2,     0,     0,     0, &
         1,  -1,   0,   4,   1,        -3,       0,   2,     0,     0,     0, &
         1,  -1,   0,   4,   0,        12,       0,   0,     0,     0,     0 /
    data IX65/  -2,   0,   2,   5,   2,         3,      0,   -1,     0,     0,     0, &
         1,  -2,   2,   2,   1,        -4,       0,   2,     0,     0,     0, &
         3,  -1,   0,   2,   0,         4,       0,   0,     0,     0,     0, &
         1,  -1,   2,   2,   0,         6,       0,   0,     0,     0,     0, &
         0,   0,   2,   3,   1,         5,       0,   -3,     0,     0,     0, &
         -1,   1,   2,   4,   1,         4,       0,   -2,     0,     0,     0, &
         0,   1,   2,   3,   2,        -6,       0,   3,     0,     0,     0, &
         -1,   0,   4,   2,   1,         4,       0,   -2,     0,     0,     0, &
         2,   0,   2,   1,   1,         6,       0,   -3,     0,     0,     0, &
         5,   0,   0,   0,   0,         6,       0,   0,     0,     0,     0 /
    data IX66/   2,   1,   2,   1,   2,        -6,      0,   3,     0,     0,     0, &
         1,   0,   4,   0,   1,         3,       0,   -2,     0,     0,     0, &
         3,   1,   2,   0,   1,         7,       0,   -4,     0,     0,     0, &
         3,   0,   4,  -2,   2,         4,       0,   -2,     0,     0,     0, &
         -2,  -1,   2,   6,   2,        -5,       0,   2,     0,     0,     0, &
         0,   0,   0,   6,   0,         5,       0,   0,     0,     0,     0, &
         0,  -2,   2,   4,   2,        -6,       0,   3,     0,     0,     0, &
         -2,   0,   2,   6,   1,        -6,       0,   3,     0,     0,     0, &
         2,   0,   0,   4,   1,        -4,       0,   2,     0,     0,     0, &
         2,   0,   0,   4,   0,        10,       0,   0,     0,     0,     0 /
    data IX67/   2,  -2,   2,   2,   2,        -4,      0,   2,     0,     0,     0, &
         0,   0,   2,   4,   0,         7,       0,   0,     0,     0,     0, &
         1,   0,   2,   3,   2,         7,       0,   -3,     0,     0,     0, &
         4,   0,   0,   2,   0,         4,       0,   0,     0,     0,     0, &
         2,   0,   2,   2,   0,        11,       0,   0,     0,     0,     0, &
         0,   0,   4,   2,   2,         5,       0,   -2,     0,     0,     0, &
         4,  -1,   2,   0,   2,        -6,       0,   2,     0,     0,     0, &
         3,   0,   2,   1,   2,         4,       0,   -2,     0,     0,     0, &
         2,   1,   2,   2,   1,         3,       0,   -2,     0,     0,     0, &
         4,   1,   2,   0,   2,         5,       0,   -2,     0,     0,     0 /
    data IX68/  -1,  -1,   2,   6,   2,        -4,      0,   2,     0,     0,     0, &
         -1,   0,   2,   6,   1,        -4,       0,   2,     0,     0,     0, &
         1,  -1,   2,   4,   1,        -3,       0,   2,     0,     0,     0, &
         1,   1,   2,   4,   2,         4,       0,   -2,     0,     0,     0, &
         3,   1,   2,   2,   2,         3,       0,   -1,     0,     0,     0, &
         5,   0,   2,   0,   1,        -3,       0,   1,     0,     0,     0, &
         2,  -1,   2,   4,   2,        -3,       0,   1,     0,     0,     0, &
         2,   0,   2,   4,   1,        -3,       0,     2,     0,     0,     0      /
    
    !****  Initialize the values and sum over the series
    
    dpsi_lsu = 0.0d0
    deps_lsu = 0.0d0
    
    cent = (epoch-DJ2000) / 36525.d0
    
    do i = num_ls, 1, -1
       
       !         Sum the mulitpliers by the arguments to the argument of
       !         nutation
       arg = 0.d0
       do j = 1,5
          
          !            Sum into the argument for nutation.
          arg = arg + nutc_int(j,i)*ls_arg(j)
       end do
       
       arg = mod(arg, 2.d0*pi)
       carg = cos(arg)
       sarg = sin(arg)
       
       !****      Now add contributions to dpsi and deps
       dpsi_lsu = dpsi_lsu + (nutc_int( 6,i)+ nutc_int(7,i)*cent)*sarg + nutc_int(10,i)*carg
       deps_lsu = deps_lsu + (nutc_int( 8,i)+ nutc_int(9,i)*cent)*carg + nutc_int(11,i)*sarg
       
    end do
    
    !     Convert values from 0.1 micro-arc-sec to mill-arc-second
    dpsi_ls = dpsi_lsu * 1.d-4
    deps_ls = deps_lsu * 1.d-4
    
  end subroutine eval_ls_nut
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  !TITLE OUT_PLAN_NUT
  
  !*********************************************************************************************************************************
  subroutine out_plan_nut   
    use SUFR_kinds, only: double
    implicit none
    
    !     Routine to write the planetary contribution to the nutations
    !     to stdout.
    !     
    
    ! USAGE:
    !     call out_plan_nut
    
    ! APPOXIMATIONS: The Oppolzer terms have not been added (should be
    !                < 0.005 mas), and
    !                Contributions from a non-rigid Earth have not been
    !                computed.  For many of these terms the contribution
    !                arises from the perturbation of the Earth's orbit and
    !                therefore there will be not deformation effects.
    
    
    ! PASSED VARIABLES
    !   NONE
    
    ! LOCAL VARIABLES
    
    !   epoch       - Dummy Julian date (used in call to plan_angles).
    !   plan_arg(14) - Values of the planetary arguments 
    !           L, L', F, D, Om, Mercury, Venus, Earth, Mars, 
    !           Jupiter, Saturn, Uranus, Uranus(?), pa.  (rads)
    !   plan_rat(14) - Rates of changes of the planetary arguments
    !                  (rad/year).  Used to get periods for the
    !                 terms.
    !   dpsi, deps   - Dummy returns for nutations in long and oblquity.
    
    real(double) :: epoch, plan_arg(14), plan_rat(14), dpsi, deps
    
    !****  Set the epoch to a dummy value
    
    epoch = 2400000.5d0 
    
    !***** Get the fundamental arguments at this epoch
    
    call plan_angles( epoch, plan_arg, plan_rat)
    
    !     Now compute the contributions of the planetery ntations by 
    !     summing over the series.
    
    call eval_plan_nut( plan_arg, plan_rat, dpsi, deps, 'NO ' )
    
  end subroutine out_plan_nut
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  !TITLE PLAN_NUT                                
  
  !*********************************************************************************************************************************
  subroutine plan_nut( jd, dpsi, deps )   
    use SUFR_kinds, only: double
    implicit none
    
    !     Routine to compute the planetary contribution to the nutations.
    !     Coefficents from Tables XIV to XIX of Kinoshita, H. and J. Souchay,
    !     Nutations for the rigid Earth, Celes. Mech. and Dynam. Astron,
    !
    ! NEW Version based on: Corrections and new developments in rigid Earth
    !     nutation theory: Lunisolar influence including indirect planetary
    !     effects, J. Souchay and H. Kinioshita, Astron. and Astrophys., 1995.
    ! (Version here based on data files: SKRE1997.DPSI and SKRE1997.DEPS
    !  and generated with ks_plan_SKRE97.f)
    
    ! Arguments based on Souchay, Loysel, Kinoshita, Folgueira, Corrections
    !    and new developments in rigid earth nutation theory, Aston. Astrophys.
    !    Suppl. Ser, 135, 111-131, (1999)
    !     
    
    ! USAGE:
    !     call plan_nut( jd, dpsi, deps )   
    !     where <jd>    is a full julian date with fractional part
    !                   of the day added (REAL(DOUBLE) INPUT)
    !     and <dpsi> and <deps> are the contributions to the nutations
    !                   in longitude and obliquity in milliarcsec.
    !                   (REAL(DOUBLE) OUTPUT)
    
    ! RESTRICTIONS: if <jd> is less than 2000000.0 this routine
    !               assumes an MJD has been passed and the time
    !               used will be converted to JD.  A warning 
    !               message will be printed.
    ! APPOXIMATIONS: The Oppolzer terms have not been added (should be
    !                < 0.005 mas), and
    !                Contributions from a non-rigid Earth have not been
    !                computed.  For many of these terms the contribution
    !                arises from the perturbation of the Earth's orbit and
    !                therefore there will be not deformation effects.
    
    
    ! PASSED VARIABLES
    !
    ! INPUT Values
    ! jd     - Time at which value needed. (jd + fraction of day) 
    
    ! OUTPUT Values
    ! dpsi   - Contribution to the nutation in longitude (mas).  Should
    !          be added to standard nutation in longitude values.
    ! deps   - Contribution to the nutation in obliquity (mas).  Should
    !          be added to standard nutation in obliquity values.
    
    real(double), intent(in) :: jd
    real(double), intent(out) :: dpsi, deps     
    
    ! LOCAL VARIABLES
    
    !   epoch       - Julian date (jd passed in unless the JD 
    !                 appears to be an MJD in which case it is 
    !                 converted to JD (2 400 000.5d0 added) 
    !   plan_arg(10) - Values of the planetary arguments  (rads) 
    !   plan_rat(14) - Rates of planetary arguments (rads/yr)
    
    real(double) :: epoch, plan_arg(14), plan_rat(14)
    
    !***** Check to make sure user passed JD and not MJD.  Correct
    !     problem and warn the user.
    ! MvdS: remove this 'solution'
    !if( jd .lt.2000000.0d0  ) then
    !             write(*,100) jd    
    !    100      format('**WARNING** MJD apparently passed to SD_COMP',
    !        .          ' Value (',F10.2,') converted to JD')
    !   epoch = jd + 2 400 000.5d0 
    !else
    !   epoch = jd
    !end if
    epoch = jd
    
    !***** Get the fundamental arguments at this epoch
    
    call plan_angles( epoch, plan_arg, plan_rat)
    
    !     Now compute the contributions of the planetery ntations by 
    !     summing over the series.
    
    call eval_plan_nut( plan_arg, plan_rat, dpsi, deps, 'NO ' )
    
  end subroutine plan_nut
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  
  !TITLE EVAL_PLAN_NUT
  
  !*********************************************************************************************************************************
  subroutine eval_plan_nut( plan_arg, plan_rat, dpsi, deps, out ) 
    use SUFR_kinds, only: double
    implicit none
    
    !     Routine to compute the planetary nutations by summing over the
    !     KS1990 coefficients.  The coefficients and their arguments are
    !     saved here in are integers in micro-arc-seconds.  
    
    ! NOTE: plan_angles must be called before routine.
    
    ! PARAMETERS:
    
    ! num_plan  - Number of contributions to the planetary nutations
    
    integer :: num_plan 
    parameter ( num_plan      = 687 )
    
    real(double) :: pi
    parameter ( pi            = 3.1415926535897932D0 )
    
    ! PASSED PARAMETERS:
    
    ! INPUT:
    ! plan_arg(14)  - 14 planetary arguments including pa as given
    !                 (KS1990.) (rad)
    ! plan_rat(14)  - 14 Planetary argument rates (rad/yr)
    
    ! OUTPUT:
    ! dpsi, deps   - Contributions to nutations in longitude and 
    !                obliquity (mas).
    
    real(double), intent(in) :: plan_arg(14), plan_rat(14)
    real(double), intent(out) :: dpsi, deps
    
    character, intent(in) :: out*(*)
    
    ! LOCAL VARIABLES:
    
    
    !  i and j   - Counters to loop over the coeffients and the argumemts
    
    integer :: i,j
    
    !  arg       - Final summed argumemt for the nutations contributions (rads)
    !  dargdt    - Rate of change of the argument (rads/yr)
    !  period    - Period of the nutation in days.
    !  amp       - Total Amplitude of the planetary nutation.  (To be used for
    !              sorting size)
    !   carg, sarg   - cosine and sin of arguments.
    
    real(double) :: arg, dargdt, period, amp, carg, sarg
    
    
    ! ks_plan: Series based on:skre97_rigid.plan
    ! IX01-IX69(18,10) - Integer values of the planetary arguments and 
    !                    values (0.1 micro-arc-seconds for values)
    integer :: IX01(18,10), IX02(18,10), IX03(18,10), IX04(18,10), IX05(18,10), IX06(18,10)
    integer :: IX07(18,10), IX08(18,10), IX09(18,10), IX10(18,10), IX11(18,10), IX12(18,10)
    integer :: IX13(18,10), IX14(18,10), IX15(18,10), IX16(18,10), IX17(18,10), IX18(18,10)
    integer :: IX19(18,10), IX20(18,10), IX21(18,10), IX22(18,10)
    integer :: IX23(18,10), IX24(18,10), IX25(18,10), IX26(18,10), IX27(18,10), IX28(18,10)
    integer :: IX29(18,10), IX30(18,10), IX31(18,10), IX32(18,10)
    integer :: IX33(18,10), IX34(18,10), IX35(18,10), IX36(18,10), IX37(18,10), IX38(18,10)
    integer :: IX39(18,10), IX40(18,10), IX41(18,10), IX42(18,10)
    integer :: IX43(18,10), IX44(18,10), IX45(18,10), IX46(18,10), IX47(18,10), IX48(18,10)
    integer :: IX49(18,10), IX50(18,10), IX51(18,10), IX52(18,10)
    integer :: IX53(18,10), IX54(18,10), IX55(18,10), IX56(18,10), IX57(18,10), IX58(18,10)
    integer :: IX59(18,10), IX60(18,10), IX61(18,10), IX62(18,10)
    integer :: IX63(18,10), IX64(18,10), IX65(18,10), IX66(18,10), IX67(18,10), IX68(18,10)
    integer :: IX69(18, 7)
    
    integer :: Plan_int(18,687)
    
    
    equivalence (Plan_int(1,  1),IX01)
    equivalence (Plan_int(1, 11),IX02)
    equivalence (Plan_int(1, 21),IX03)
    equivalence (Plan_int(1, 31),IX04)
    equivalence (Plan_int(1, 41),IX05)
    equivalence (Plan_int(1, 51),IX06)
    equivalence (Plan_int(1, 61),IX07)
    equivalence (Plan_int(1, 71),IX08)
    equivalence (Plan_int(1, 81),IX09)
    equivalence (Plan_int(1, 91),IX10)
    equivalence (Plan_int(1,101),IX11)
    equivalence (Plan_int(1,111),IX12)
    equivalence (Plan_int(1,121),IX13)
    equivalence (Plan_int(1,131),IX14)
    equivalence (Plan_int(1,141),IX15)
    equivalence (Plan_int(1,151),IX16)
    equivalence (Plan_int(1,161),IX17)
    equivalence (Plan_int(1,171),IX18)
    equivalence (Plan_int(1,181),IX19)
    equivalence (Plan_int(1,191),IX20)
    equivalence (Plan_int(1,201),IX21)
    equivalence (Plan_int(1,211),IX22)
    equivalence (Plan_int(1,221),IX23)
    equivalence (Plan_int(1,231),IX24)
    equivalence (Plan_int(1,241),IX25)
    equivalence (Plan_int(1,251),IX26)
    equivalence (Plan_int(1,261),IX27)
    equivalence (Plan_int(1,271),IX28)
    equivalence (Plan_int(1,281),IX29)
    equivalence (Plan_int(1,291),IX30)
    equivalence (Plan_int(1,301),IX31)
    equivalence (Plan_int(1,311),IX32)
    equivalence (Plan_int(1,321),IX33)
    equivalence (Plan_int(1,331),IX34)
    equivalence (Plan_int(1,341),IX35)
    equivalence (Plan_int(1,351),IX36)
    equivalence (Plan_int(1,361),IX37)
    equivalence (Plan_int(1,371),IX38)
    equivalence (Plan_int(1,381),IX39)
    equivalence (Plan_int(1,391),IX40)
    equivalence (Plan_int(1,401),IX41)
    equivalence (Plan_int(1,411),IX42)
    equivalence (Plan_int(1,421),IX43)
    equivalence (Plan_int(1,431),IX44)
    equivalence (Plan_int(1,441),IX45)
    equivalence (Plan_int(1,451),IX46)
    equivalence (Plan_int(1,461),IX47)
    equivalence (Plan_int(1,471),IX48)
    equivalence (Plan_int(1,481),IX49)
    equivalence (Plan_int(1,491),IX50)
    equivalence (Plan_int(1,501),IX51)
    equivalence (Plan_int(1,511),IX52)
    equivalence (Plan_int(1,521),IX53)
    equivalence (Plan_int(1,531),IX54)
    equivalence (Plan_int(1,541),IX55)
    equivalence (Plan_int(1,551),IX56)
    equivalence (Plan_int(1,561),IX57)
    equivalence (Plan_int(1,571),IX58)
    equivalence (Plan_int(1,581),IX59)
    equivalence (Plan_int(1,591),IX60)
    equivalence (Plan_int(1,601),IX61)
    equivalence (Plan_int(1,611),IX62)
    equivalence (Plan_int(1,621),IX63)
    equivalence (Plan_int(1,631),IX64)
    equivalence (Plan_int(1,641),IX65)
    equivalence (Plan_int(1,651),IX66)
    equivalence (Plan_int(1,661),IX67)
    equivalence (Plan_int(1,671),IX68)
    equivalence (Plan_int(1,681),IX69)
    
    data IX01/   0,   0,   0,   0,   0,   0,   0,   8, -16,   4,   5,   0,   0,   0,  1440,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -8,  16,  -4,   -5,   0,   0,   2,    56,  -117,   -42,   -40, &
         0,   0,   0,   0,   0,   0,   0,   8, -16,   4,   5,   0,   0,   2,   125,   -43,     0,   -54, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   2,   2,     0,     5,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -4,   8,  -1,   -5,   0,   0,   2,     3,    -7,    -3,     0, &
         0,   0,   0,   0,   0,   0,   0,   4,  -8,   3,   0,   0,   0,   1,     3,     0,     0,    -2, &
         0,   0,   1,  -1,   1,   0,   0,   3,  -8,   3,   0,   0,   0,   0,  -114,     0,     0,    61, &
         -1,   0,   0,   0,   0,   0,  10,  -3,   0,   0,   0,   0,   0,   0,  -219,    89,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,   6,  -3,   0,   2,    -3,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   4,  -8,   3,   0,   0,   0,   0,  -462,  1604,     0,     0 /
    
    data IX02/   0,   0,   1,  -1,   1,   0,   0,  -5,   8,  -3,   0,   0,   0,   0,    99,     0,     0,   -53, &
         0,   0,   0,   0,   0,   0,   0,  -4,   8,  -3,   0,   0,   0,   1,    -3,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,   0,   4,  -8,   1,   5,   0,   0,   2,     0,     6,     2,     0, &
         0,   0,   0,   0,   0,   0,  -5,   6,   4,   0,   0,   0,   0,   2,     3,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   -5,   0,   0,   2,   -12,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   -5,   0,   0,   1,    14,  -218,   117,     8, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   2,   -5,   0,   0,   0,    31,  -481,  -257,   -17, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   -5,   0,   0,   0,  -491,   128,     0,     0, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,  -2,   5,   0,   0,   0, -3084,  5123,  2735,  1647, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,   5,   0,   0,   1, -1444,  2409, -1286,  -771 /
    
    data IX03/   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,   5,   0,   0,   2,    11,   -24,   -11,    -9, &
         2,   0,  -1,  -1,   0,   0,   0,   3,  -7,   0,   0,   0,   0,   0,    26,    -9,     0,     0, &
         1,   0,   0,  -2,   0,   0,  19, -21,   3,   0,   0,   0,   0,   0,   103,   -60,     0,     0, &
         0,   0,   1,  -1,   1,   0,   2,  -4,   0,  -3,   0,   0,   0,   0,     0,   -13,    -7,     0, &
         1,   0,   0,  -1,   1,   0,   0,  -1,   0,   2,   0,   0,   0,   0,   -26,   -29,   -16,    14, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,  -4,   10,   0,   0,   0,     9,   -27,   -14,    -5, &
         -2,   0,   0,   2,   1,   0,   0,   2,   0,   0,   -5,   0,   0,   0,    12,     0,     0,    -6, &
         0,   0,   0,   0,   0,   0,   3,  -7,   4,   0,   0,   0,   0,   0,    -7,     0,     0,     0, &
         0,   0,  -1,   1,   0,   0,   0,   1,   0,   1,   -1,   0,   0,   0,     0,    24,     0,     0, &
         -2,   0,   0,   2,   1,   0,   0,   2,   0,  -2,   0,   0,   0,   0,   284,     0,     0,  -151 /
    
    data IX04/  -1,   0,   0,   0,   0,   0,  18, -16,   0,   0,   0,   0,   0,   0,   226,   101,     0,     0, &
         -2,   0,   1,   1,   2,   0,   0,   1,   0,  -2,   0,   0,   0,   0,     0,    -8,    -2,     0, &
         -1,   0,   1,  -1,   1,   0,  18, -17,   0,   0,   0,   0,   0,   0,     0,    -6,    -3,     0, &
         -1,   0,   0,   1,   1,   0,   0,   2,  -2,   0,   0,   0,   0,   0,     5,     0,     0,    -3, &
         0,   0,   0,   0,   0,   0,  -8,  13,   0,   0,   0,   0,   0,   2,   -41,   175,    76,    17, &
         0,   0,   2,  -2,   2,   0,  -8,  11,   0,   0,   0,   0,   0,   0,     0,    15,     6,     0, &
         0,   0,   0,   0,   0,   0,  -8,  13,   0,   0,   0,   0,   0,   1,   425,   212,  -133,   269, &
         0,   0,   1,  -1,   1,   0,  -8,  12,   0,   0,   0,   0,   0,   0,  1200,   598,   319,  -641, &
         0,   0,   0,   0,   0,   0,   8, -13,   0,   0,   0,   0,   0,   0,   235,   334,     0,     0, &
         0,   0,   1,  -1,   1,   0,   8, -14,   0,   0,   0,   0,   0,   0,    11,   -12,    -7,    -6 /
    
    data IX05/   0,   0,   0,   0,   0,   0,   8, -13,   0,   0,   0,   0,   0,   1,     5,    -6,     3,     3, &
         -2,   0,   0,   2,   1,   0,   0,   2,   0,  -4,   5,   0,   0,   0,    -5,     0,     0,     3, &
         -2,   0,   0,   2,   2,   0,   3,  -3,   0,   0,   0,   0,   0,   0,     6,     0,     0,    -3, &
         -2,   0,   0,   2,   0,   0,   0,   2,   0,  -3,   1,   0,   0,   0,    15,     0,     0,     0, &
         0,   0,   0,   0,   1,   0,   3,  -5,   0,   2,   0,   0,   0,   0,    13,     0,     0,    -7, &
         -2,   0,   0,   2,   0,   0,   0,   2,   0,  -4,   3,   0,   0,   0,    -6,    -9,     0,     0, &
         0,   0,  -1,   1,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   266,   -78,     0,     0, &
         0,   0,   0,   0,   1,   0,   0,  -1,   2,   0,   0,   0,   0,   0,  -460,  -435,  -232,   246, &
         0,   0,   1,  -1,   2,   0,   0,  -2,   2,   0,   0,   0,   0,   0,     0,    15,     7,     0, &
         -1,   0,   1,   0,   1,   0,   3,  -5,   0,   0,   0,   0,   0,   0,    -3,     0,     0,     2 /
    
    data IX06/  -1,   0,   0,   1,   0,   0,   3,  -4,   0,   0,   0,   0,   0,   0,     0,   131,     0,     0, &
         -2,   0,   0,   2,   0,   0,   0,   2,   0,  -2,   -2,   0,   0,   0,     4,     0,     0,     0, &
         -2,   0,   2,   0,   2,   0,   0,  -5,   9,   0,   0,   0,   0,   0,     0,     3,     0,     0, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   0,   0,  -1,   0,   0,     0,     4,     2,     0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,     0,     3,     0,     0, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   0,   0,   0,   2,   0,   -17,   -19,   -10,     9, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   1,    -9,   -11,     6,    -5, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   2,    -6,     0,     0,     3, &
         -1,   0,   0,   1,   0,   0,   0,   3,  -4,   0,   0,   0,   0,   0,   -16,     8,     0,     0, &
         0,   0,  -1,   1,   0,   0,   0,   1,   0,   0,   2,   0,   0,   0,     0,     3,     0,     0 /
    
    data IX07/   0,   0,   1,  -1,   2,   0,   0,  -1,   0,   0,   2,   0,   0,   0,    11,    24,    11,    -5, &
         0,   0,   0,   0,   1,   0,   0,  -9,  17,   0,   0,   0,   0,   0,    -3,    -4,    -2,     1, &
         0,   0,   0,   0,   2,   0,  -3,   5,   0,   0,   0,   0,   0,   0,     3,     0,     0,    -1, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,  -1,   2,   0,   0,   0,     0,    -8,    -4,     0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   -2,   0,   0,   0,     0,     3,     0,     0, &
         1,   0,   0,  -2,   0,   0,  17, -16,   0,  -2,   0,   0,   0,   0,     0,     5,     0,     0, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   1,   -3,   0,   0,   0,     0,     3,     2,     0, &
         -2,   0,   0,   2,   1,   0,   0,   5,  -6,   0,   0,   0,   0,   0,    -6,     4,     2,     3, &
         0,   0,  -2,   2,   0,   0,   0,   9, -13,   0,   0,   0,   0,   0,    -3,    -5,     0,     0, &
         0,   0,   1,  -1,   2,   0,   0,  -1,   0,   0,   1,   0,   0,   0,    -5,     0,     0,     2 /
    
    data IX08/   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,     4,    24,    13,    -2, &
         0,   0,  -1,   1,   0,   0,   0,   1,   0,   0,   1,   0,   0,   0,   -42,    20,     0,     0, &
         0,   0,  -2,   2,   0,   0,   5,  -6,   0,   0,   0,   0,   0,   0,   -10,   233,     0,     0, &
         0,   0,  -1,   1,   1,   0,   5,  -7,   0,   0,   0,   0,   0,   0,    -3,     0,     0,     1, &
         -2,   0,   0,   2,   0,   0,   6,  -8,   0,   0,   0,   0,   0,   0,    78,   -18,     0,     0, &
         2,   0,   1,  -3,   1,   0,  -6,   7,   0,   0,   0,   0,   0,   0,     0,     3,     1,     0, &
         0,   0,   0,   0,   2,   0,   0,   0,   0,   1,   0,   0,   0,   0,     0,    -3,    -1,     0, &
         0,   0,  -1,   1,   1,   0,   0,   1,   0,   1,   0,   0,   0,   0,     0,    -4,    -2,     1, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   0,   0,   2,   0,   0,     0,    -8,    -4,    -1, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   1,     0,    -5,     3,     0 /
    
    data IX09/   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   2,    -7,     0,     0,     3, &
         0,   0,   0,   0,   0,   0,   0,  -8,  15,   0,   0,   0,   0,   2,   -14,     8,     3,     6, &
         0,   0,   0,   0,   0,   0,   0,  -8,  15,   0,   0,   0,   0,   1,     0,     8,    -4,     0, &
         0,   0,   1,  -1,   1,   0,   0,  -9,  15,   0,   0,   0,   0,   0,     0,    19,    10,     0, &
         0,   0,   0,   0,   0,   0,   0,   8, -15,   0,   0,   0,   0,   0,    45,   -22,     0,     0, &
         1,   0,  -1,  -1,   0,   0,   0,   8, -15,   0,   0,   0,   0,   0,    -3,     0,     0,     0, &
         2,   0,   0,  -2,   0,   0,   2,  -5,   0,   0,   0,   0,   0,   0,     0,    -3,     0,     0, &
         -2,   0,   0,   2,   0,   0,   0,   2,   0,  -5,   5,   0,   0,   0,     0,     3,     0,     0, &
         2,   0,   0,  -2,   1,   0,   0,  -6,   8,   0,   0,   0,   0,   0,     3,     5,     3,    -2, &
         2,   0,   0,  -2,   1,   0,   0,  -2,   0,   3,   0,   0,   0,   0,    89,   -16,    -9,   -48 /
    
    data IX10/  -2,   0,   1,   1,   0,   0,   0,   1,   0,  -3,   0,   0,   0,   0,     0,     3,     0,     0, &
         -2,   0,   1,   1,   1,   0,   0,   1,   0,  -3,   0,   0,   0,   0,    -3,     7,     4,     2, &
         -2,   0,   0,   2,   0,   0,   0,   2,   0,  -3,   0,   0,   0,   0,  -349,   -62,     0,     0, &
         -2,   0,   0,   2,   0,   0,   0,   6,  -8,   0,   0,   0,   0,   0,   -15,    22,     0,     0, &
         -2,   0,   0,   2,   0,   0,   0,   2,   0,  -1,   -5,   0,   0,   0,    -3,     0,     0,     0, &
         -1,   0,   0,   1,   0,   0,   0,   1,   0,  -1,   0,   0,   0,   0,   -53,     0,     0,     0, &
         -1,   0,   1,   1,   1,   0, -20,  20,   0,   0,   0,   0,   0,   0,     5,     0,     0,    -3, &
         1,   0,   0,  -2,   0,   0,  20, -21,   0,   0,   0,   0,   0,   0,     0,    -8,     0,     0, &
         0,   0,   0,   0,   1,   0,   0,   8, -15,   0,   0,   0,   0,   0,    15,    -7,    -4,    -8, &
         0,   0,   2,  -2,   1,   0,   0, -10,  15,   0,   0,   0,   0,   0,    -3,     0,     0,     1 /
    
    data IX11/   0,   0,  -1,   1,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0,   -21,   -78,     0,     0, &
         0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0,   0,    20,   -70,   -37,   -11, &
         0,   0,   1,  -1,   2,   0,   0,  -1,   0,   1,   0,   0,   0,   0,     0,     6,     3,     0, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,  -2,   4,   0,   0,   0,     5,     3,     2,    -2, &
         2,   0,   0,  -2,   1,   0,  -6,   8,   0,   0,   0,   0,   0,   0,   -17,    -4,    -2,     9, &
         0,   0,  -2,   2,   1,   0,   5,  -6,   0,   0,   0,   0,   0,   0,     0,     6,     3,     0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   -1,   0,   0,   1,    32,    15,    -8,    17, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   0,   -1,   0,   0,   0,   174,    84,    45,   -93, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,    11,    56,     0,     0, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   0,   1,   0,   0,   0,   -66,   -12,    -6,    35 /
    
    data IX12/   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   1,    47,     8,     4,   -25, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   2,     0,     8,     4,     0, &
         0,   0,   2,  -2,   1,   0,   0,  -9,  13,   0,   0,   0,   0,   0,    10,   -22,   -12,    -5, &
         0,   0,   0,   0,   1,   0,   0,   7, -13,   0,   0,   0,   0,   0,    -3,     0,     0,     2, &
         -2,   0,   0,   2,   0,   0,   0,   5,  -6,   0,   0,   0,   0,   0,   -24,    12,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   9, -17,   0,   0,   0,   0,   0,     5,    -6,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -9,  17,   0,   0,   0,   0,   2,     3,     0,     0,    -2, &
         1,   0,   0,  -1,   1,   0,   0,  -3,   4,   0,   0,   0,   0,   0,     4,     3,     1,    -2, &
         1,   0,   0,  -1,   1,   0,  -3,   4,   0,   0,   0,   0,   0,   0,     0,    29,    15,     0, &
         0,   0,   0,   0,   2,   0,   0,  -1,   2,   0,   0,   0,   0,   0,    -5,    -4,    -2,     2 /
    
    data IX13/   0,   0,  -1,   1,   1,   0,   0,   0,   2,   0,   0,   0,   0,   0,     8,    -3,    -1,    -5, &
         0,   0,  -2,   2,   0,   1,   0,  -2,   0,   0,   0,   0,   0,   0,     0,    -3,     0,     0, &
         0,   0,   0,   0,   0,   0,   3,  -5,   0,   2,   0,   0,   0,   0,    10,     0,     0,     0, &
         -2,   0,   0,   2,   1,   0,   0,   2,   0,  -3,   1,   0,   0,   0,     3,     0,     0,    -2, &
         -2,   0,   0,   2,   1,   0,   3,  -3,   0,   0,   0,   0,   0,   0,    -5,     0,     0,     3, &
         0,   0,   0,   0,   1,   0,   8, -13,   0,   0,   0,   0,   0,   0,    46,    66,    35,   -25, &
         0,   0,  -1,   1,   0,   0,   8, -12,   0,   0,   0,   0,   0,   0,   -14,     7,     0,     0, &
         0,   0,   2,  -2,   1,   0,  -8,  11,   0,   0,   0,   0,   0,   0,     0,     3,     2,     0, &
         -1,   0,   0,   1,   0,   0,   0,   2,  -2,   0,   0,   0,   0,   0,    -5,     0,     0,     0, &
         -1,   0,   0,   0,   1,   0,  18, -16,   0,   0,   0,   0,   0,   0,   -68,   -34,   -18,    36 /
    
    data IX14/   0,   0,   1,  -1,   1,   0,   0,  -1,   0,  -1,   1,   0,   0,   0,     0,    14,     7,     0, &
         0,   0,   0,   0,   1,   0,   3,  -7,   4,   0,   0,   0,   0,   0,    10,    -6,    -3,    -5, &
         -2,   0,   1,   1,   1,   0,   0,  -3,   7,   0,   0,   0,   0,   0,    -5,    -4,    -2,     3, &
         0,   0,   1,  -1,   2,   0,   0,  -1,   0,  -2,   5,   0,   0,   0,    -3,     5,     2,     1, &
         0,   0,   0,   0,   1,   0,   0,   0,   0,  -2,   5,   0,   0,   0,    76,    17,     9,   -41, &
         0,   0,   0,   0,   1,   0,   0,  -4,   8,  -3,   0,   0,   0,   0,    84,   298,   159,   -45, &
         1,   0,   0,   0,   1,   0, -10,   3,   0,   0,   0,   0,   0,   0,     3,     0,     0,    -1, &
         0,   0,   2,  -2,   1,   0,   0,  -2,   0,   0,   0,   0,   0,   0,    -3,     0,     0,     2, &
         -1,   0,   0,   0,   1,   0,  10,  -3,   0,   0,   0,   0,   0,   0,    -3,     0,     0,     1, &
         0,   0,   0,   0,   1,   0,   0,   4,  -8,   3,   0,   0,   0,   0,   -82,   292,   156,    44 /
    
    data IX15/   0,   0,   0,   0,   1,   0,   0,   0,   0,   2,   -5,   0,   0,   0,   -73,    17,     9,    39, &
         0,   0,  -1,   1,   0,   0,   0,   1,   0,   2,   -5,   0,   0,   0,    -9,   -16,     0,     0, &
         2,   0,  -1,  -1,   1,   0,   0,   3,  -7,   0,   0,   0,   0,   0,     3,     0,    -1,    -2, &
         -2,   0,   0,   2,   0,   0,   0,   2,   0,   0,   -5,   0,   0,   0,    -3,     0,     0,     0, &
         0,   0,   0,   0,   1,   0,  -3,   7,  -4,   0,   0,   0,   0,   0,    -9,    -5,    -3,     5, &
         -2,   0,   0,   2,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   0,  -439,     0,     0,     0, &
         1,   0,   0,   0,   1,   0, -18,  16,   0,   0,   0,   0,   0,   0,    57,   -28,   -15,   -30, &
         -2,   0,   1,   1,   1,   0,   0,   1,   0,  -2,   0,   0,   0,   0,     0,    -6,    -3,     0, &
         0,   0,   1,  -1,   2,   0,  -8,  12,   0,   0,   0,   0,   0,   0,    -4,     0,     0,     2, &
         0,   0,   0,   0,   1,   0,  -8,  13,   0,   0,   0,   0,   0,   0,   -40,    57,    30,    21 /
    
    data IX16/   0,   0,   0,   0,   0,   0,   0,   1,  -2,   0,   0,   0,   0,   1,    23,     7,     3,   -13, &
         0,   0,   1,  -1,   1,   0,   0,   0,  -2,   0,   0,   0,   0,   0,   273,    80,    43,  -146, &
         0,   0,   0,   0,   0,   0,   0,   1,  -2,   0,   0,   0,   0,   0,  -449,   430,     0,     0, &
         0,   0,   1,  -1,   1,   0,   0,  -2,   2,   0,   0,   0,   0,   0,    -8,   -47,   -25,     4, &
         0,   0,   0,   0,   0,   0,   0,  -1,   2,   0,   0,   0,   0,   1,     6,    47,    25,    -3, &
         -1,   0,   0,   1,   1,   0,   3,  -4,   0,   0,   0,   0,   0,   0,     0,    23,    13,     0, &
         -1,   0,   0,   1,   1,   0,   0,   3,  -4,   0,   0,   0,   0,   0,    -3,     0,     0,     2, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   0,   -2,   0,   0,   0,     3,    -4,    -2,    -2, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   0,   2,   0,   0,   0,   -48,  -110,   -59,    26, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   1,    51,   114,    61,   -27 /
    
    data IX17/   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   2,  -133,     0,     0,    57, &
         0,   0,   1,  -1,   0,   0,   3,  -6,   0,   0,   0,   0,   0,   0,     0,     4,     0,     0, &
         0,   0,   0,   0,   1,   0,  -3,   5,   0,   0,   0,   0,   0,   0,   -21,    -6,    -3,    11, &
         0,   0,   1,  -1,   2,   0,  -3,   4,   0,   0,   0,   0,   0,   0,     0,    -3,    -1,     0, &
         0,   0,   0,   0,   1,   0,   0,  -2,   4,   0,   0,   0,   0,   0,   -11,   -21,   -11,     6, &
         0,   0,   2,  -2,   1,   0,  -5,   6,   0,   0,   0,   0,   0,   0,   -18,  -436,  -233,     9, &
         0,   0,  -1,   1,   0,   0,   5,  -7,   0,   0,   0,   0,   0,   0,    35,    -7,     0,     0, &
         0,   0,   0,   0,   1,   0,   5,  -8,   0,   0,   0,   0,   0,   0,     0,     5,     3,     0, &
         -2,   0,   0,   2,   1,   0,   6,  -8,   0,   0,   0,   0,   0,   0,    11,    -3,    -1,    -6, &
         0,   0,   0,   0,   1,   0,   0,  -8,  15,   0,   0,   0,   0,   0,    -5,    -3,    -1,     3 /
    
    data IX18/  -2,   0,   0,   2,   1,   0,   0,   2,   0,  -3,   0,   0,   0,   0,   -53,    -9,    -5,    28, &
         -2,   0,   0,   2,   1,   0,   0,   6,  -8,   0,   0,   0,   0,   0,     0,     3,     2,     1, &
         1,   0,   0,  -1,   1,   0,   0,  -1,   0,   1,   0,   0,   0,   0,     4,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   -5,   0,   0,   0,     0,    -4,     0,     0, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,  -1,   0,   0,   0,   0,   -50,   194,   103,    27, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   0,   0,   0,   1,   -13,    52,    28,     7, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   -91,   248,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   1,     6,    49,    26,    -3, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   1,   0,   0,   0,   0,    -6,   -47,   -25,     3, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   1,     0,     5,     3,     0 /
    
    data IX19/   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   2,    52,    23,    10,   -23, &
         0,   0,   1,  -1,   2,   0,   0,  -1,   0,   0,   -1,   0,   0,   0,    -3,     0,     0,     1, &
         0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   -1,   0,   0,   0,     0,     5,     3,     0, &
         0,   0,  -1,   1,   0,   0,   0,   1,   0,   0,   -1,   0,   0,   0,    -4,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -7,  13,   0,   0,   0,   0,   2,    -4,     8,     3,     2, &
         0,   0,   0,   0,   0,   0,   0,   7, -13,   0,   0,   0,   0,   0,    10,     0,     0,     0, &
         2,   0,   0,  -2,   1,   0,   0,  -5,   6,   0,   0,   0,   0,   0,     3,     0,     0,    -2, &
         0,   0,   2,  -2,   1,   0,   0,  -8,  11,   0,   0,   0,   0,   0,     0,     8,     4,     0, &
         0,   0,   2,  -2,   1,  -1,   0,   2,   0,   0,   0,   0,   0,   0,     0,     8,     4,     1, &
         -2,   0,   0,   2,   0,   0,   0,   4,  -4,   0,   0,   0,   0,   0,    -4,     0,     0,     0 /
    
    data IX20/   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   -2,   0,   0,   0,    -4,     0,     0,     0, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   0,   3,   0,   0,   0,    -8,     4,     2,     4, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   1,     8,    -4,    -2,    -4, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   2,     0,    15,     7,     0, &
         -2,   0,   0,   2,   0,   0,   3,  -3,   0,   0,   0,   0,   0,   0,  -138,     0,     0,     0, &
         0,   0,   0,   0,   2,   0,   0,  -4,   8,  -3,   0,   0,   0,   0,     0,    -7,    -3,     0, &
         0,   0,   0,   0,   2,   0,   0,   4,  -8,   3,   0,   0,   0,   0,     0,    -7,    -3,     0, &
         2,   0,   0,  -2,   1,   0,   0,  -2,   0,   2,   0,   0,   0,   0,    54,     0,     0,   -29, &
         0,   0,   1,  -1,   2,   0,   0,  -1,   0,   2,   0,   0,   0,   0,     0,    10,     4,     0, &
         0,   0,   1,  -1,   2,   0,   0,   0,  -2,   0,   0,   0,   0,   0,    -7,     0,     0,     3 /
    
    data IX21/   0,   0,   0,   0,   1,   0,   0,   1,  -2,   0,   0,   0,   0,   0,   -37,    35,    19,    20, &
         0,   0,  -1,   1,   0,   0,   0,   2,  -2,   0,   0,   0,   0,   0,     0,     4,     0,     0, &
         0,   0,  -1,   1,   0,   0,   0,   1,   0,   0,   -2,   0,   0,   0,    -4,     9,     0,     0, &
         0,   0,   2,  -2,   1,   0,   0,  -2,   0,   0,   2,   0,   0,   0,     8,     0,     0,    -4, &
         0,   0,   1,  -1,   1,   0,   3,  -6,   0,   0,   0,   0,   0,   0,    -9,   -14,    -8,     5, &
         0,   0,   0,   0,   0,   0,   3,  -5,   0,   0,   0,   0,   0,   1,    -3,    -9,    -5,     3, &
         0,   0,   0,   0,   0,   0,   3,  -5,   0,   0,   0,   0,   0,   0,  -145,    47,     0,     0, &
         0,   0,   1,  -1,   1,   0,  -3,   4,   0,   0,   0,   0,   0,   0,   -10,    40,    21,     5, &
         0,   0,   0,   0,   0,   0,  -3,   5,   0,   0,   0,   0,   0,   1,    11,   -49,   -26,    -7, &
         0,   0,   0,   0,   0,   0,  -3,   5,   0,   0,   0,   0,   0,   2, -2150,     0,     0,   932 /
    
    data IX22/   0,   0,   2,  -2,   2,   0,  -3,   3,   0,   0,   0,   0,   0,   0,   -12,     0,     0,     5, &
         0,   0,   0,   0,   0,   0,  -3,   5,   0,   0,   0,   0,   0,   2,    85,     0,     0,   -37, &
         0,   0,   0,   0,   0,   0,   0,   2,  -4,   0,   0,   0,   0,   1,     4,     0,     0,    -2, &
         0,   0,   1,  -1,   1,   0,   0,   1,  -4,   0,   0,   0,   0,   0,     3,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   0,   2,  -4,   0,   0,   0,   0,   0,   -86,   153,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -2,   4,   0,   0,   0,   0,   1,    -6,     9,     5,     3, &
         0,   0,   1,  -1,   1,   0,   0,  -3,   4,   0,   0,   0,   0,   0,     9,   -13,    -7,    -5, &
         0,   0,   0,   0,   0,   0,   0,  -2,   4,   0,   0,   0,   0,   1,    -8,    12,     6,     4, &
         0,   0,   0,   0,   0,   0,   0,  -2,   4,   0,   0,   0,   0,   2,   -51,     0,     0,    22, &
         0,   0,   0,   0,   0,   0,  -5,   8,   0,   0,   0,   0,   0,   2,   -11,  -268,  -116,     5 /
    
    data IX23/   0,   0,   2,  -2,   2,   0,  -5,   6,   0,   0,   0,   0,   0,   0,     0,    12,     5,     0, &
         0,   0,   0,   0,   0,   0,  -5,   8,   0,   0,   0,   0,   0,   2,     0,     7,     3,     0, &
         0,   0,   0,   0,   0,   0,  -5,   8,   0,   0,   0,   0,   0,   1,    31,     6,     3,   -17, &
         0,   0,   1,  -1,   1,   0,  -5,   7,   0,   0,   0,   0,   0,   0,   140,    27,    14,   -75, &
         0,   0,   0,   0,   0,   0,  -5,   8,   0,   0,   0,   0,   0,   1,    57,    11,     6,   -30, &
         0,   0,   0,   0,   0,   0,   5,  -8,   0,   0,   0,   0,   0,   0,   -14,   -39,     0,     0, &
         0,   0,   1,  -1,   2,   0,   0,  -1,   0,  -1,   0,   0,   0,   0,     0,    -6,    -2,     0, &
         0,   0,   0,   0,   1,   0,   0,   0,   0,  -1,   0,   0,   0,   0,     4,    15,     8,    -2, &
         0,   0,  -1,   1,   0,   0,   0,   1,   0,  -1,   0,   0,   0,   0,     0,     4,     0,     0, &
         0,   0,   2,  -2,   1,   0,   0,  -2,   0,   1,   0,   0,   0,   0,    -3,     0,     0,     1 /
    
    data IX24/   0,   0,   0,   0,   0,   0,   0,  -6,  11,   0,   0,   0,   0,   2,     0,    11,     5,     0, &
         0,   0,   0,   0,   0,   0,   0,   6, -11,   0,   0,   0,   0,   0,     9,     6,     0,     0, &
         0,   0,   0,   0,   0,  -1,   0,   4,   0,   0,   0,   0,   0,   2,    -4,    10,     4,     2, &
         0,   0,   0,   0,   0,   1,   0,  -4,   0,   0,   0,   0,   0,   0,     5,     3,     0,     0, &
         2,   0,   0,  -2,   1,   0,  -3,   3,   0,   0,   0,   0,   0,   0,    16,     0,     0,    -9, &
         -2,   0,   0,   2,   0,   0,   0,   2,   0,   0,   -2,   0,   0,   0,    -3,     0,     0,     0, &
         0,   0,   2,  -2,   1,   0,   0,  -7,   9,   0,   0,   0,   0,   0,     0,     3,     2,    -1, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   4,   -5,   0,   0,   2,     7,     0,     0,    -3, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   -25,    22,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   1,    42,   223,   119,   -22 /
    
    data IX25/   0,   0,   1,  -1,   1,   0,   0,  -1,   0,   2,   0,   0,   0,   0,   -27,  -143,   -77,    14, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   1,     9,    49,    26,    -5, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   2, -1166,     0,     0,   505, &
         0,   0,   2,  -2,   2,   0,   0,  -2,   0,   2,   0,   0,   0,   0,    -5,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   5,   0,   0,   2,    -6,     0,     0,     3, &
         0,   0,   0,   0,   1,   0,   3,  -5,   0,   0,   0,   0,   0,   0,    -8,     0,     1,     4, &
         0,   0,  -1,   1,   0,   0,   3,  -4,   0,   0,   0,   0,   0,   0,     0,    -4,     0,     0, &
         0,   0,   2,  -2,   1,   0,  -3,   3,   0,   0,   0,   0,   0,   0,   117,     0,     0,   -63, &
         0,   0,   0,   0,   1,   0,   0,   2,  -4,   0,   0,   0,   0,   0,    -4,     8,     4,     2, &
         0,   0,   2,  -2,   1,   0,   0,  -4,   4,   0,   0,   0,   0,   0,     3,     0,     0,    -2 /
    
    data IX26/   0,   0,   1,  -1,   2,   0,  -5,   7,   0,   0,   0,   0,   0,   0,    -5,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,   0,   3,  -6,   0,   0,   0,   0,   0,     0,    31,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -3,   6,   0,   0,   0,   0,   1,    -5,     0,     1,     3, &
         0,   0,   1,  -1,   1,   0,   0,  -4,   6,   0,   0,   0,   0,   0,     4,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   0,  -3,   6,   0,   0,   0,   0,   1,    -4,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,   0,  -3,   6,   0,   0,   0,   0,   2,   -24,   -13,    -6,    10, &
         0,   0,  -1,   1,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,     3,     0,     0,     0, &
         0,   0,   0,   0,   1,   0,   2,  -3,   0,   0,   0,   0,   0,   0,     0,   -32,   -17,     0, &
         0,   0,   0,   0,   0,   0,   0,  -5,   9,   0,   0,   0,   0,   2,     8,    12,     5,    -3, &
         0,   0,   0,   0,   0,   0,   0,  -5,   9,   0,   0,   0,   0,   1,     3,     0,     0,    -1 /
    
    data IX27/   0,   0,   0,   0,   0,   0,   0,   5,  -9,   0,   0,   0,   0,   0,     7,    13,     0,     0, &
         0,   0,  -1,   1,   0,   0,   0,   1,   0,  -2,   0,   0,   0,   0,    -3,    16,     0,     0, &
         0,   0,   2,  -2,   1,   0,   0,  -2,   0,   2,   0,   0,   0,   0,    50,     0,     0,   -27, &
         -2,   0,   1,   1,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,     0,    -5,    -3,     0, &
         0,   0,  -2,   2,   0,   0,   3,  -3,   0,   0,   0,   0,   0,   0,    13,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,  -6,  10,   0,   0,   0,   0,   0,   1,     0,     5,     3,     1, &
         0,   0,   0,   0,   0,   0,  -6,  10,   0,   0,   0,   0,   0,   2,    24,     5,     2,   -11, &
         0,   0,   0,   0,   0,   0,  -2,   3,   0,   0,   0,   0,   0,   2,     5,   -11,    -5,    -2, &
         0,   0,   0,   0,   0,   0,  -2,   3,   0,   0,   0,   0,   0,   1,    30,    -3,    -2,   -16, &
         0,   0,   1,  -1,   1,   0,  -2,   2,   0,   0,   0,   0,   0,   0,    18,     0,     0,    -9 /
    
    data IX28/   0,   0,   0,   0,   0,   0,   2,  -3,   0,   0,   0,   0,   0,   0,     8,   614,     0,     0, &
         0,   0,   0,   0,   0,   0,   2,  -3,   0,   0,   0,   0,   0,   1,     3,    -3,    -1,    -2, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   0,   1,     6,    17,     9,    -3, &
         0,   0,   1,  -1,   1,   0,   0,  -1,   0,   3,   0,   0,   0,   0,    -3,    -9,    -5,     2, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   0,   1,     0,     6,     3,    -1, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   0,   2,  -127,    21,     9,    55, &
         0,   0,   0,   0,   0,   0,   0,   4,  -8,   0,   0,   0,   0,   0,     3,     5,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -4,   8,   0,   0,   0,   0,   2,    -6,   -10,    -4,     3, &
         0,   0,  -2,   2,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   0,     5,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -4,   7,   0,   0,   0,   0,   2,    16,     9,     4,    -7 /
    
    data IX29/   0,   0,   0,   0,   0,   0,   0,  -4,   7,   0,   0,   0,   0,   1,     3,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   0,   4,  -7,   0,   0,   0,   0,   0,     0,    22,     0,     0, &
         0,   0,   0,   0,   1,   0,  -2,   3,   0,   0,   0,   0,   0,   0,     0,    19,    10,     0, &
         0,   0,   2,  -2,   1,   0,   0,  -2,   0,   3,   0,   0,   0,   0,     7,     0,     0,    -4, &
         0,   0,   0,   0,   0,   0,   0,  -5,  10,   0,   0,   0,   0,   2,     0,    -5,    -2,     0, &
         0,   0,   0,   0,   1,   0,  -1,   2,   0,   0,   0,   0,   0,   0,     0,     3,     1,     0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   4,   0,   0,   0,   2,    -9,     3,     1,     4, &
         0,   0,   0,   0,   0,   0,   0,  -3,   5,   0,   0,   0,   0,   2,    17,     0,     0,    -7, &
         0,   0,   0,   0,   0,   0,   0,  -3,   5,   0,   0,   0,   0,   1,     0,    -3,    -2,    -1, &
         0,   0,   0,   0,   0,   0,   0,   3,  -5,   0,   0,   0,   0,   0,   -20,    34,     0,     0 /
    
    data IX30/   0,   0,   0,   0,   0,   0,   1,  -2,   0,   0,   0,   0,   0,   1,   -10,     0,     1,     5, &
         0,   0,   1,  -1,   1,   0,   1,  -3,   0,   0,   0,   0,   0,   0,    -4,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,   1,  -2,   0,   0,   0,   0,   0,   0,    22,   -87,     0,     0, &
         0,   0,   0,   0,   0,   0,  -1,   2,   0,   0,   0,   0,   0,   1,    -4,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,  -1,   2,   0,   0,   0,   0,   0,   2,    -3,    -6,    -2,     1, &
         0,   0,   0,   0,   0,   0,  -7,  11,   0,   0,   0,   0,   0,   2,   -16,    -3,    -1,     7, &
         0,   0,   0,   0,   0,   0,  -7,  11,   0,   0,   0,   0,   0,   1,     0,    -3,    -2,     0, &
         0,   0,  -2,   2,   0,   0,   4,  -4,   0,   0,   0,   0,   0,   0,     4,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   2,  -3,   0,   0,   0,   0,   0,   -68,    39,     0,     0, &
         0,   0,   2,  -2,   1,   0,  -4,   4,   0,   0,   0,   0,   0,   0,    27,     0,     0,   -14 /
    
    data IX31/   0,   0,  -1,   1,   0,   0,   4,  -5,   0,   0,   0,   0,   0,   0,     0,    -4,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,  -1,   0,   0,   0,   0,   0,   -25,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,  -4,   7,   0,   0,   0,   0,   0,   1,   -12,    -3,    -2,     6, &
         0,   0,   1,  -1,   1,   0,  -4,   6,   0,   0,   0,   0,   0,   0,     3,     0,     0,    -1, &
         0,   0,   0,   0,   0,   0,  -4,   7,   0,   0,   0,   0,   0,   2,     3,    66,    29,    -1, &
         0,   0,   0,   0,   0,   0,  -4,   6,   0,   0,   0,   0,   0,   2,   490,     0,     0,  -213, &
         0,   0,   0,   0,   0,   0,  -4,   6,   0,   0,   0,   0,   0,   1,   -22,    93,    49,    12, &
         0,   0,   1,  -1,   1,   0,  -4,   5,   0,   0,   0,   0,   0,   0,    -7,    28,    15,     4, &
         0,   0,   0,   0,   0,   0,  -4,   6,   0,   0,   0,   0,   0,   1,    -3,    13,     7,     2, &
         0,   0,   0,   0,   0,   0,   4,  -6,   0,   0,   0,   0,   0,   0,   -46,    14,     0,     0 /
    
    data IX32/  -2,   0,   0,   2,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,    -5,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,     2,     1,     0,     0, &
         0,   0,  -1,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,     0,    -3,     0,     0, &
         0,   0,   0,   0,   1,   0,   1,  -1,   0,   0,   0,   0,   0,   0,   -28,     0,     0,    15, &
         0,   0,   0,   0,   0,   0,   0,  -1,   0,   5,   0,   0,   0,   2,     5,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   0,   1,  -3,   0,   0,   0,   0,   0,     0,     3,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -1,   3,   0,   0,   0,   0,   2,   -11,     0,     0,     5, &
         0,   0,   0,   0,   0,   0,   0,  -7,  12,   0,   0,   0,   0,   2,     0,     3,     1,     0, &
         0,   0,   0,   0,   0,   0,  -1,   1,   0,   0,   0,   0,   0,   2,    -3,     0,     0,     1, &
         0,   0,   0,   0,   0,   0,  -1,   1,   0,   0,   0,   0,   0,   1,    25,   106,    57,   -13 /
    
    data IX33/   0,   0,   1,  -1,   1,   0,  -1,   0,   0,   0,   0,   0,   0,   0,     5,    21,    11,    -3, &
         0,   0,   0,   0,   0,   0,   1,  -1,   0,   0,   0,   0,   0,   0,  1485,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   1,  -1,   0,   0,   0,   0,   0,   1,    -7,   -32,   -17,     4, &
         0,   0,   1,  -1,   1,   0,   1,  -2,   0,   0,   0,   0,   0,   0,     0,     5,     3,     0, &
         0,   0,   0,   0,   0,   0,   0,  -2,   5,   0,   0,   0,   0,   2,    -6,    -3,    -2,     3, &
         0,   0,   0,   0,   0,   0,   0,  -1,   0,   4,   0,   0,   0,   2,    30,    -6,    -2,   -13, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,  -4,   0,   0,   0,   0,    -4,     4,     0,     0, &
         0,   0,   0,   0,   1,   0,  -1,   1,   0,   0,   0,   0,   0,   0,   -19,     0,     0,    10, &
         0,   0,   0,   0,   0,   0,   0,  -6,  10,   0,   0,   0,   0,   2,     0,     4,     2,    -1, &
         0,   0,   0,   0,   0,   0,   0,  -6,  10,   0,   0,   0,   0,   0,     0,     3,     0,     0 /
    
    data IX34/   0,   0,   2,  -2,   1,   0,   0,  -3,   0,   3,   0,   0,   0,   0,     4,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   0,  -3,   7,   0,   0,   0,   0,   2,     0,    -3,    -1,     0, &
         -2,   0,   0,   2,   0,   0,   4,  -4,   0,   0,   0,   0,   0,   0,    -3,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -5,   8,   0,   0,   0,   0,   2,     5,     3,     1,    -2, &
         0,   0,   0,   0,   0,   0,   0,   5,  -8,   0,   0,   0,   0,   0,     0,    11,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -1,   0,   3,   0,   0,   0,   2,   118,     0,     0,   -52, &
         0,   0,   0,   0,   0,   0,   0,  -1,   0,   3,   0,   0,   0,   1,     0,    -5,    -3,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,  -3,   0,   0,   0,   0,   -28,    36,     0,     0, &
         0,   0,   0,   0,   0,   0,   2,  -4,   0,   0,   0,   0,   0,   0,     5,    -5,     0,     0, &
         0,   0,   0,   0,   0,   0,  -2,   4,   0,   0,   0,   0,   0,   1,    14,   -59,   -31,    -8 /
    
    data IX35/   0,   0,   1,  -1,   1,   0,  -2,   3,   0,   0,   0,   0,   0,   0,     0,     9,     5,     1, &
         0,   0,   0,   0,   0,   0,  -2,   4,   0,   0,   0,   0,   0,   2,  -458,     0,     0,   198, &
         0,   0,   0,   0,   0,   0,  -6,   9,   0,   0,   0,   0,   0,   2,     0,   -45,   -20,     0, &
         0,   0,   0,   0,   0,   0,  -6,   9,   0,   0,   0,   0,   0,   1,     9,     0,     0,    -5, &
         0,   0,   0,   0,   0,   0,   6,  -9,   0,   0,   0,   0,   0,   0,     0,    -3,     0,     0, &
         0,   0,   0,   0,   1,   0,   0,   1,   0,  -2,   0,   0,   0,   0,     0,    -4,    -2,    -1, &
         0,   0,   2,  -2,   1,   0,  -2,   2,   0,   0,   0,   0,   0,   0,    11,     0,     0,    -6, &
         0,   0,   0,   0,   0,   0,   0,  -4,   6,   0,   0,   0,   0,   2,     6,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   0,   4,  -6,   0,   0,   0,   0,   0,   -16,    23,     0,     0, &
         0,   0,   0,   0,   1,   0,   3,  -4,   0,   0,   0,   0,   0,   0,     0,    -4,    -2,     0 /
    
    data IX36/   0,   0,   0,   0,   0,   0,   0,  -1,   0,   2,   0,   0,   0,   2,    -5,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,  -2,   0,   0,   0,   0,  -166,   269,     0,     0, &
         0,   0,   0,   0,   1,   0,   0,   1,   0,  -1,   0,   0,   0,   0,    15,     0,     0,    -8, &
         0,   0,   0,   0,   0,   0,  -5,   9,   0,   0,   0,   0,   0,   2,    10,     0,     0,    -4, &
         0,   0,   0,   0,   0,   0,   0,   3,  -4,   0,   0,   0,   0,   0,   -78,    45,     0,     0, &
         0,   0,   0,   0,   0,   0,  -3,   4,   0,   0,   0,   0,   0,   2,     0,    -5,    -2,     0, &
         0,   0,   0,   0,   0,   0,  -3,   4,   0,   0,   0,   0,   0,   1,     7,     0,     0,    -4, &
         0,   0,   0,   0,   0,   0,   3,  -4,   0,   0,   0,   0,   0,   0,    -5,   328,     0,     0, &
         0,   0,   0,   0,   0,   0,   3,  -4,   0,   0,   0,   0,   0,   1,     3,     0,     0,    -2, &
         0,   0,   0,   0,   1,   0,   0,   2,  -2,   0,   0,   0,   0,   0,     5,     0,     0,    -2 /
    
    data IX37/   0,   0,   0,   0,   1,   0,   0,  -1,   0,   2,   0,   0,   0,   0,     0,     3,     1,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   -3,   0,   0,   0,    -3,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   -5,   0,   0,   0,    -3,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -1,   0,   1,   0,   0,   0,   1,     0,    -4,    -2,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,  -1,   0,   0,   0,   0, -1223,   -26,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,  -1,   0,   0,   0,   1,     0,     7,     3,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,  -3,   5,   0,   0,   0,     3,     0,     0,     0, &
         0,   0,   0,   0,   1,   0,  -3,   4,   0,   0,   0,   0,   0,   0,     0,     3,     2,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   -2,   0,   0,   0,    -6,    20,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   2,  -2,   0,   0,   0,   0,   0,  -368,     0,     0,     0 /
    
    data IX38/   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   -1,   0,   0,   0,   -75,     0,     0,     0, &
         0,   0,   0,   0,   1,   0,   0,  -1,   0,   1,   0,   0,   0,   0,    11,     0,     0,    -6, &
         0,   0,   0,   0,   1,   0,   0,  -2,   2,   0,   0,   0,   0,   0,     3,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,  -8,  14,   0,   0,   0,   0,   0,   2,    -3,     0,     0,     1, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   2,   -5,   0,   0,   0,   -13,   -30,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   5,  -8,   3,   0,   0,   0,   0,    21,     3,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   5,  -8,   3,   0,   0,   0,   2,    -3,     0,     0,     1, &
         0,   0,   0,   0,   0,   0,   0,  -1,   0,   0,   0,   0,   0,   1,    -4,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,     8,   -27,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   3,  -8,   3,   0,   0,   0,   0,   -19,   -11,     0,     0 /
    
    data IX39/   0,   0,   0,   0,   0,   0,   0,  -3,   8,  -3,   0,   0,   0,   2,    -4,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,  -2,   5,   0,   0,   2,     0,     5,     2,     0, &
         0,   0,   0,   0,   0,   0,  -8,  12,   0,   0,   0,   0,   0,   2,    -6,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,  -8,  12,   0,   0,   0,   0,   0,   0,    -8,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   -2,   0,   0,   0,    -1,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   1,   0,   0,   2,   -14,     0,     0,     6, &
         0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,     6,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   2,   -74,     0,     0,    32, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   2,   0,   0,   2,     0,    -3,    -1,     0, &
         0,   0,   2,  -2,   1,   0,  -5,   5,   0,   0,   0,   0,   0,   0,     4,     0,     0,    -2 /
    
    data IX40/   0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0,     8,    11,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   1,     0,     3,     2,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   2,  -262,     0,     0,   114, &
         0,   0,   0,   0,   0,   0,   3,  -6,   0,   0,   0,   0,   0,   0,     0,    -4,     0,     0, &
         0,   0,   0,   0,   0,   0,  -3,   6,   0,   0,   0,   0,   0,   1,    -7,     0,     0,     4, &
         0,   0,   0,   0,   0,   0,  -3,   6,   0,   0,   0,   0,   0,   2,     0,   -27,   -12,     0, &
         0,   0,   0,   0,   0,   0,   0,  -1,   4,   0,   0,   0,   0,   2,   -19,    -8,    -4,     8, &
         0,   0,   0,   0,   0,   0,  -5,   7,   0,   0,   0,   0,   0,   2,   202,     0,     0,   -87, &
         0,   0,   0,   0,   0,   0,  -5,   7,   0,   0,   0,   0,   0,   1,    -8,    35,    19,     5, &
         0,   0,   1,  -1,   1,   0,  -5,   6,   0,   0,   0,   0,   0,   0,     0,     4,     2,     0 /
    
    data IX41/   0,   0,   0,   0,   0,   0,   5,  -7,   0,   0,   0,   0,   0,   0,    16,    -5,     0,     0, &
         0,   0,   2,  -2,   1,   0,   0,  -1,   0,   1,   0,   0,   0,   0,     5,     0,     0,    -3, &
         0,   0,   0,   0,   0,   0,   0,  -1,   0,   1,   0,   0,   0,   0,     0,    -3,     0,     0, &
         0,   0,   0,   0,   0,  -1,   0,   3,   0,   0,   0,   0,   0,   2,     1,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   2,   0,   0,   0,   2,   -35,   -48,   -21,    15, &
         0,   0,   0,   0,   0,   0,   0,  -2,   6,   0,   0,   0,   0,   2,    -3,    -5,    -2,     1, &
         0,   0,   0,   0,   1,   0,   2,  -2,   0,   0,   0,   0,   0,   0,     6,     0,     0,    -3, &
         0,   0,   0,   0,   0,   0,   0,  -6,   9,   0,   0,   0,   0,   2,     3,     0,     0,    -1, &
         0,   0,   0,   0,   0,   0,   0,   6,  -9,   0,   0,   0,   0,   0,     0,    -5,     0,     0, &
         0,   0,   0,   0,   0,   0,  -2,   2,   0,   0,   0,   0,   0,   1,    12,    55,    29,    -6 /
    
    data IX42/   0,   0,   1,  -1,   1,   0,  -2,   1,   0,   0,   0,   0,   0,   0,     0,     5,     3,     0, &
         0,   0,   0,   0,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,  -598,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   1,    -3,   -13,    -7,     1, &
         0,   0,   0,   0,   0,   0,   0,   1,   0,   3,   0,   0,   0,   2,    -5,    -7,    -3,     2, &
         0,   0,   0,   0,   0,   0,   0,  -5,   7,   0,   0,   0,   0,   2,     3,     0,     0,    -1, &
         0,   0,   0,   0,   0,   0,   0,   5,  -7,   0,   0,   0,   0,   0,     5,    -7,     0,     0, &
         0,   0,   0,   0,   1,   0,  -2,   2,   0,   0,   0,   0,   0,   0,     4,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   0,   4,  -5,   0,   0,   0,   0,   0,    16,    -6,     0,     0, &
         0,   0,   0,   0,   0,   0,   1,  -3,   0,   0,   0,   0,   0,   0,     8,    -3,     0,     0, &
         0,   0,   0,   0,   0,   0,  -1,   3,   0,   0,   0,   0,   0,   1,     8,   -31,   -16,    -4 /
    
    data IX43/   0,   0,   1,  -1,   1,   0,  -1,   2,   0,   0,   0,   0,   0,   0,     0,     3,     1,     0, &
         0,   0,   0,   0,   0,   0,  -1,   3,   0,   0,   0,   0,   0,   2,   113,     0,     0,   -49, &
         0,   0,   0,   0,   0,   0,  -7,  10,   0,   0,   0,   0,   0,   2,     0,   -24,   -10,     0, &
         0,   0,   0,   0,   0,   0,  -7,  10,   0,   0,   0,   0,   0,   1,     4,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,   0,   0,   0,    27,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,  -4,   8,   0,   0,   0,   0,   0,   2,    -3,     0,     0,     1, &
         0,   0,   0,   0,   0,   0,  -4,   5,   0,   0,   0,   0,   0,   2,     0,    -4,    -2,     0, &
         0,   0,   0,   0,   0,   0,  -4,   5,   0,   0,   0,   0,   0,   1,     5,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   4,  -5,   0,   0,   0,   0,   0,   0,     0,    -3,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   2,   -13,     0,     0,     6 /
    
    data IX44/   0,   0,   0,   0,   0,   0,   0,  -2,   0,   5,   0,   0,   0,   2,     5,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   0,   0,   2,   -18,   -10,    -4,     8, &
         0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,    -4,   -28,     0,     0, &
         0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   2,    -5,     6,     3,     2, &
         0,   0,   0,   0,   0,   0,  -9,  13,   0,   0,   0,   0,   0,   2,    -3,     0,     0,     1, &
         0,   0,   0,   0,   0,   0,   0,  -1,   5,   0,   0,   0,   0,   2,    -5,    -9,    -4,     2, &
         0,   0,   0,   0,   0,   0,   0,  -2,   0,   4,   0,   0,   0,   2,    17,     0,     0,    -7, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,  -4,   0,   0,   0,   0,    11,     4,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,  -2,   7,   0,   0,   0,   0,   2,     0,    -6,    -2,     0, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,  -3,   0,   0,   0,   0,    83,    15,     0,     0 /
    
    data IX45/   0,   0,   0,   0,   0,   0,  -2,   5,   0,   0,   0,   0,   0,   1,    -4,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,  -2,   5,   0,   0,   0,   0,   0,   2,     0,  -114,   -49,     0, &
         0,   0,   0,   0,   0,   0,  -6,   8,   0,   0,   0,   0,   0,   2,   117,     0,     0,   -51, &
         0,   0,   0,   0,   0,   0,  -6,   8,   0,   0,   0,   0,   0,   1,    -5,    19,    10,     2, &
         0,   0,   0,   0,   0,   0,   6,  -8,   0,   0,   0,   0,   0,   0,    -3,     0,     0,     0, &
         0,   0,   0,   0,   1,   0,   0,   2,   0,  -2,   0,   0,   0,   0,    -3,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,   0,  -3,   9,   0,   0,   0,   0,   2,     0,    -3,    -1,     0, &
         0,   0,   0,   0,   0,   0,   0,   5,  -6,   0,   0,   0,   0,   0,     3,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   5,  -6,   0,   0,   0,   0,   2,     0,    -6,    -2,     0, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   0,   393,     3,     0,     0 /
    
    data IX46/   0,   0,   0,   0,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   1,    -4,    21,    11,     2, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   2,    -6,     0,    -1,     3, &
         0,   0,   0,   0,   0,   0,  -5,  10,   0,   0,   0,   0,   0,   2,    -3,     8,     4,     1, &
         0,   0,   0,   0,   0,   0,   0,   4,  -4,   0,   0,   0,   0,   0,     8,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   4,  -4,   0,   0,   0,   0,   2,    18,   -29,   -13,    -8, &
         0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   0,   0,   0,   1,     8,    34,    18,    -4, &
         0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,   0,   0,   0,   0,    89,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,   0,   0,   0,   1,     3,    12,     6,    -1, &
         0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,   0,   0,   0,   2,    54,   -15,    -7,   -24, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   -3,   0,   0,   0,     0,     3,     0,     0 /
    
    data IX47/   0,   0,   0,   0,   0,   0,   0,  -5,  13,   0,   0,   0,   0,   2,     3,     0,     0,    -1, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,  -1,   0,   0,   0,   0,     0,    35,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,  -1,   0,   0,   0,   2,  -154,   -30,   -13,    67, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   -2,   0,   0,   0,    15,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   -2,   0,   0,   1,     0,     4,     2,     0, &
         0,   0,   0,   0,   0,   0,   0,   3,  -2,   0,   0,   0,   0,   0,     0,     9,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   3,  -2,   0,   0,   0,   0,   2,    80,   -71,   -31,   -35, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   -1,   0,   0,   2,     0,   -20,    -9,     0, &
         0,   0,   0,   0,   0,   0,   0,  -6,  15,   0,   0,   0,   0,   2,    11,     5,     2,    -5, &
         0,   0,   0,   0,   0,   0,  -8,  15,   0,   0,   0,   0,   0,   2,    61,   -96,   -42,   -27 /
    
    data IX48/   0,   0,   0,   0,   0,   0,  -3,   9,  -4,   0,   0,   0,   0,   2,    14,     9,     4,    -6, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   2,   -5,   0,   0,   2,   -11,    -6,    -3,     5, &
         0,   0,   0,   0,   0,   0,   0,  -2,   8,  -1,   -5,   0,   0,   2,     0,    -3,    -1,     0, &
         0,   0,   0,   0,   0,   0,   0,   6,  -8,   3,   0,   0,   0,   2,   123,  -415,  -180,   -53, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,     0,     0,     0,   -35, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,    -5,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   1,     7,   -32,   -17,    -4, &
         0,   0,   1,  -1,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,     0,    -9,    -5,     0, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   1,     0,    -4,     2,     0, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   2,   -89,     0,     0,    38 /
    
    data IX49/   0,   0,   0,   0,   0,   0,   0,  -6,  16,  -4,   -5,   0,   0,   2,     0,   -86,   -19,    -6, &
         0,   0,   0,   0,   0,   0,   0,  -2,   8,  -3,   0,   0,   0,   2,     0,     0,   -19,     6, &
         0,   0,   0,   0,   0,   0,   0,  -2,   8,  -3,   0,   0,   0,   2,  -123,  -416,  -180,    53, &
         0,   0,   0,   0,   0,   0,   0,   6,  -8,   1,   5,   0,   0,   2,     0,    -3,    -1,     0, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,  -2,   5,   0,   0,   2,    12,    -6,    -3,    -5, &
         0,   0,   0,   0,   0,   0,   3,  -5,   4,   0,   0,   0,   0,   2,   -13,     9,     4,     6, &
         0,   0,   0,   0,   0,   0,  -8,  11,   0,   0,   0,   0,   0,   2,     0,   -15,    -7,     0, &
         0,   0,   0,   0,   0,   0,  -8,  11,   0,   0,   0,   0,   0,   1,     3,     0,     0,    -1, &
         0,   0,   0,   0,   0,   0,  -8,  11,   0,   0,   0,   0,   0,   2,   -62,   -97,   -42,    27, &
         0,   0,   0,   0,   0,   0,   0,  11,   0,   0,   0,   0,   0,   2,   -11,     5,     2,     5 /
    
    data IX50/   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   1,   0,   0,   2,     0,   -19,    -8,     0, &
         0,   0,   0,   0,   0,   0,   3,  -3,   0,   2,   0,   0,   0,   2,    -3,     0,     0,     1, &
         0,   0,   2,  -2,   1,   0,   0,   4,  -8,   3,   0,   0,   0,   0,     0,     4,     2,     0, &
         0,   0,   1,  -1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,     0,     3,     0,     0, &
         0,   0,   2,  -2,   1,   0,   0,  -4,   8,  -3,   0,   0,   0,   0,     0,     4,     2,     0, &
         0,   0,   0,   0,   0,   0,   0,   1,   2,   0,   0,   0,   0,   2,   -85,   -70,   -31,    37, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   1,   0,   0,   0,   2,   163,   -12,    -5,   -72, &
         0,   0,   0,   0,   0,   0,  -3,   7,   0,   0,   0,   0,   0,   2,   -63,   -16,    -7,    28, &
         0,   0,   0,   0,   0,   0,   0,   0,   4,   0,   0,   0,   0,   2,   -21,   -32,   -14,     9, &
         0,   0,   0,   0,   0,   0,  -5,   6,   0,   0,   0,   0,   0,   2,     0,    -3,    -1,     0 /
    
    data IX51/   0,   0,   0,   0,   0,   0,  -5,   6,   0,   0,   0,   0,   0,   1,     3,     0,     0,    -2, &
         0,   0,   0,   0,   0,   0,   5,  -6,   0,   0,   0,   0,   0,   0,     0,     8,     0,     0, &
         0,   0,   0,   0,   0,   0,   5,  -6,   0,   0,   0,   0,   0,   2,     3,    10,     4,    -1, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,   2,   0,   0,   0,   2,     3,     0,     0,    -1, &
         0,   0,   0,   0,   0,   0,   0,  -1,   6,   0,   0,   0,   0,   2,     0,    -7,    -3,     0, &
         0,   0,   0,   0,   0,   0,   0,   7,  -9,   0,   0,   0,   0,   2,     0,    -4,    -2,     0, &
         0,   0,   0,   0,   0,   0,   2,  -1,   0,   0,   0,   0,   0,   0,     6,    19,     0,     0, &
         0,   0,   0,   0,   0,   0,   2,  -1,   0,   0,   0,   0,   0,   2,     5,  -173,   -75,    -2, &
         0,   0,   0,   0,   0,   0,   0,   6,  -7,   0,   0,   0,   0,   2,     0,    -7,    -3,     0, &
         0,   0,   0,   0,   0,   0,   0,   5,  -5,   0,   0,   0,   0,   2,     7,   -12,    -5,    -3 /
    
    data IX52/   0,   0,   0,   0,   0,   0,  -1,   4,   0,   0,   0,   0,   0,   1,    -3,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,  -1,   4,   0,   0,   0,   0,   0,   2,     3,    -4,    -2,    -1, &
         0,   0,   0,   0,   0,   0,  -7,   9,   0,   0,   0,   0,   0,   2,    74,     0,     0,   -32, &
         0,   0,   0,   0,   0,   0,  -7,   9,   0,   0,   0,   0,   0,   1,    -3,    12,     6,     2, &
         0,   0,   0,   0,   0,   0,   0,   4,  -3,   0,   0,   0,   0,   2,    26,   -14,    -6,   -11, &
         0,   0,   0,   0,   0,   0,   0,   3,  -1,   0,   0,   0,   0,   2,    19,     0,     0,    -8, &
         0,   0,   0,   0,   0,   0,  -4,   4,   0,   0,   0,   0,   0,   1,     6,    24,    13,    -3, &
         0,   0,   0,   0,   0,   0,   4,  -4,   0,   0,   0,   0,   0,   0,    83,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   4,  -4,   0,   0,   0,   0,   0,   1,     0,   -10,    -5,     0, &
         0,   0,   0,   0,   0,   0,   4,  -4,   0,   0,   0,   0,   0,   2,    11,    -3,    -1,    -5 /
    
    data IX53/   0,   0,   0,   0,   0,   0,   0,   2,   1,   0,   0,   0,   0,   2,     3,     0,     1,    -1, &
         0,   0,   0,   0,   0,   0,   0,  -3,   0,   5,   0,   0,   0,   2,     3,     0,     0,    -1, &
         0,   0,   0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,    -4,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   1,     5,   -23,   -12,    -3, &
         0,   0,   0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   2,  -339,     0,     0,   147, &
         0,   0,   0,   0,   0,   0,  -9,  12,   0,   0,   0,   0,   0,   2,     0,   -10,    -5,     0, &
         0,   0,   0,   0,   0,   0,   0,   3,   0,  -4,   0,   0,   0,   0,     5,     0,     0,     0, &
         0,   0,   2,  -2,   1,   0,   1,  -1,   0,   0,   0,   0,   0,   0,     3,     0,     0,    -1, &
         0,   0,   0,   0,   0,   0,   0,   7,  -8,   0,   0,   0,   0,   2,     0,    -4,    -2,     0, &
         0,   0,   0,   0,   0,   0,   0,   3,   0,  -3,   0,   0,   0,   0,    18,    -3,     0,     0 /
    
    data IX54/   0,   0,   0,   0,   0,   0,   0,   3,   0,  -3,   0,   0,   0,   2,     9,   -11,    -5,    -4, &
         0,   0,   0,   0,   0,   0,  -2,   6,   0,   0,   0,   0,   0,   2,    -8,     0,     0,     4, &
         0,   0,   0,   0,   0,   0,  -6,   7,   0,   0,   0,   0,   0,   1,     3,     0,     0,    -1, &
         0,   0,   0,   0,   0,   0,   6,  -7,   0,   0,   0,   0,   0,   0,     0,     9,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   6,  -6,   0,   0,   0,   0,   2,     6,    -9,    -4,    -2, &
         0,   0,   0,   0,   0,   0,   0,   3,   0,  -2,   0,   0,   0,   0,    -4,   -12,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   3,   0,  -2,   0,   0,   0,   2,    67,   -91,   -39,   -29, &
         0,   0,   0,   0,   0,   0,   0,   5,  -4,   0,   0,   0,   0,   2,    30,   -18,    -8,   -13, &
         0,   0,   0,   0,   0,   0,   3,  -2,   0,   0,   0,   0,   0,   0,     0,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   3,  -2,   0,   0,   0,   0,   0,   2,     0,  -114,   -50,     0 /
    
    data IX55/   0,   0,   0,   0,   0,   0,   0,   3,   0,  -1,   0,   0,   0,   2,     0,     0,     0,    23, &
         0,   0,   0,   0,   0,   0,   0,   3,   0,  -1,   0,   0,   0,   2,   517,    16,     7,  -224, &
         0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   -2,   0,   0,   2,     0,    -7,    -3,     0, &
         0,   0,   0,   0,   0,   0,   0,   4,  -2,   0,   0,   0,   0,   2,   143,    -3,    -1,   -62, &
         0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   -1,   0,   0,   2,    29,     0,     0,   -13, &
         0,   0,   2,  -2,   1,   0,   0,   1,   0,  -1,   0,   0,   0,   0,    -4,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,  -8,  16,   0,   0,   0,   0,   0,   2,    -6,     0,     0,     3, &
         0,   0,   0,   0,   0,   0,   0,   3,   0,   2,   -5,   0,   0,   2,     5,    12,     5,    -2, &
         0,   0,   0,   0,   0,   0,   0,   7,  -8,   3,   0,   0,   0,   2,   -25,     0,     0,    11, &
         0,   0,   0,   0,   0,   0,   0,  -5,  16,  -4,   -5,   0,   0,   2,    -3,     0,     0,     1 /
    
    data IX56/   0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   0,   0,   0,   2,     0,     4,     2,     0, &
         0,   0,   0,   0,   0,   0,   0,  -1,   8,  -3,   0,   0,   0,   2,   -22,    12,     5,    10, &
         0,   0,   0,   0,   0,   0,  -8,  10,   0,   0,   0,   0,   0,   2,    50,     0,     0,   -22, &
         0,   0,   0,   0,   0,   0,  -8,  10,   0,   0,   0,   0,   0,   1,     0,     7,     4,     0, &
         0,   0,   0,   0,   0,   0,  -8,  10,   0,   0,   0,   0,   0,   2,     0,     3,     1,     0, &
         0,   0,   0,   0,   0,   0,   0,   2,   2,   0,   0,   0,   0,   2,    -4,     4,     2,     2, &
         0,   0,   0,   0,   0,   0,   0,   3,   0,   1,   0,   0,   0,   2,    -5,   -11,    -5,     2, &
         0,   0,   0,   0,   0,   0,  -3,   8,   0,   0,   0,   0,   0,   2,     0,     4,     2,     0, &
         0,   0,   0,   0,   0,   0,  -5,   5,   0,   0,   0,   0,   0,   1,     4,    17,     9,    -2, &
         0,   0,   0,   0,   0,   0,   5,  -5,   0,   0,   0,   0,   0,   0,    59,     0,     0,     0 /
    
    data IX57/   0,   0,   0,   0,   0,   0,   5,  -5,   0,   0,   0,   0,   0,   1,     0,    -4,    -2,     0, &
         0,   0,   0,   0,   0,   0,   5,  -5,   0,   0,   0,   0,   0,   2,    -8,     0,     0,     4, &
         0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,    -3,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   1,     4,   -15,    -8,    -2, &
         0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   2,   370,    -8,     0,  -160, &
         0,   0,   0,   0,   0,   0,   0,   7,  -7,   0,   0,   0,   0,   2,     0,     0,    -3,     0, &
         0,   0,   0,   0,   0,   0,   0,   7,  -7,   0,   0,   0,   0,   2,     0,     3,     1,     0, &
         0,   0,   0,   0,   0,   0,   0,   6,  -5,   0,   0,   0,   0,   2,    -6,     3,     1,     3, &
         0,   0,   0,   0,   0,   0,   7,  -8,   0,   0,   0,   0,   0,   0,     0,     6,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   5,  -3,   0,   0,   0,   0,   2,   -10,     0,     0,     4 /
    
    data IX58/   0,   0,   0,   0,   0,   0,   4,  -3,   0,   0,   0,   0,   0,   2,     0,     9,     4,     0, &
         0,   0,   0,   0,   0,   0,   1,   2,   0,   0,   0,   0,   0,   2,     4,    17,     7,    -2, &
         0,   0,   0,   0,   0,   0,  -9,  11,   0,   0,   0,   0,   0,   2,    34,     0,     0,   -15, &
         0,   0,   0,   0,   0,   0,  -9,  11,   0,   0,   0,   0,   0,   1,     0,     5,     3,     0, &
         0,   0,   0,   0,   0,   0,   0,   4,   0,  -4,   0,   0,   0,   2,    -5,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,   0,   4,   0,  -3,   0,   0,   0,   2,   -37,    -7,    -3,    16, &
         0,   0,   0,   0,   0,   0,  -6,   6,   0,   0,   0,   0,   0,   1,     3,    13,     7,    -2, &
         0,   0,   0,   0,   0,   0,   6,  -6,   0,   0,   0,   0,   0,   0,    40,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   6,  -6,   0,   0,   0,   0,   0,   1,     0,    -3,    -2,     0, &
         0,   0,   0,   0,   0,   0,   0,   4,   0,  -2,   0,   0,   0,   2,  -184,    -3,    -1,    80 /
    
    data IX59/   0,   0,   0,   0,   0,   0,   0,   6,  -4,   0,   0,   0,   0,   2,    -3,     0,     0,     1, &
         0,   0,   0,   0,   0,   0,   3,  -1,   0,   0,   0,   0,   0,   0,    -3,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   3,  -1,   0,   0,   0,   0,   0,   1,     0,   -10,    -6,    -1, &
         0,   0,   0,   0,   0,   0,   3,  -1,   0,   0,   0,   0,   0,   2,    31,    -6,     0,   -13, &
         0,   0,   0,   0,   0,   0,   0,   4,   0,  -1,   0,   0,   0,   2,    -3,   -32,   -14,     1, &
         0,   0,   0,   0,   0,   0,   0,   4,   0,   0,   -2,   0,   0,   2,    -7,     0,     0,     3, &
         0,   0,   0,   0,   0,   0,   0,   5,  -2,   0,   0,   0,   0,   2,     0,    -8,    -4,     0, &
         0,   0,   0,   0,   0,   0,   0,   4,   0,   0,   0,   0,   0,   0,     3,    -4,     0,     0, &
         0,   0,   0,   0,   0,   0,   8,  -9,   0,   0,   0,   0,   0,   0,     0,     4,     0,     0, &
         0,   0,   0,   0,   0,   0,   5,  -4,   0,   0,   0,   0,   0,   2,     0,     3,     1,     0 /
    
    data IX60/   0,   0,   0,   0,   0,   0,   2,   1,   0,   0,   0,   0,   0,   2,    19,   -23,   -10,     2, &
         0,   0,   0,   0,   0,   0,   2,   1,   0,   0,   0,   0,   0,   1,     0,     0,     0,   -10, &
         0,   0,   0,   0,   0,   0,   2,   1,   0,   0,   0,   0,   0,   1,     0,     3,     2,     0, &
         0,   0,   0,   0,   0,   0,  -7,   7,   0,   0,   0,   0,   0,   1,     0,     9,     5,    -1, &
         0,   0,   0,   0,   0,   0,   7,  -7,   0,   0,   0,   0,   0,   0,    28,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   4,  -2,   0,   0,   0,   0,   0,   1,     0,    -7,    -4,     0, &
         0,   0,   0,   0,   0,   0,   4,  -2,   0,   0,   0,   0,   0,   2,     8,    -4,     0,    -4, &
         0,   0,   0,   0,   0,   0,   4,  -2,   0,   0,   0,   0,   0,   0,     0,     0,    -2,     0, &
         0,   0,   0,   0,   0,   0,   4,  -2,   0,   0,   0,   0,   0,   0,     0,     3,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   5,   0,  -4,   0,   0,   0,   2,    -3,     0,     0,     1 /
    
    data IX61/   0,   0,   0,   0,   0,   0,   0,   5,   0,  -3,   0,   0,   0,   2,    -9,     0,     1,     4, &
         0,   0,   0,   0,   0,   0,   0,   5,   0,  -2,   0,   0,   0,   2,     3,    12,     5,    -1, &
         0,   0,   0,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   2,    17,    -3,    -1,     0, &
         0,   0,   0,   0,   0,   0,  -8,   8,   0,   0,   0,   0,   0,   1,     0,     7,     4,     0, &
         0,   0,   0,   0,   0,   0,   8,  -8,   0,   0,   0,   0,   0,   0,    19,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   5,  -3,   0,   0,   0,   0,   0,   1,     0,    -5,    -3,     0, &
         0,   0,   0,   0,   0,   0,   5,  -3,   0,   0,   0,   0,   0,   2,    14,    -3,     0,    -1, &
         0,   0,   0,   0,   0,   0,  -9,   9,   0,   0,   0,   0,   0,   1,     0,     0,    -1,     0, &
         0,   0,   0,   0,   0,   0,  -9,   9,   0,   0,   0,   0,   0,   1,     0,     0,     0,    -5, &
         0,   0,   0,   0,   0,   0,  -9,   9,   0,   0,   0,   0,   0,   1,     0,     5,     3,     0 /
    
    data IX62/   0,   0,   0,   0,   0,   0,   9,  -9,   0,   0,   0,   0,   0,   0,    13,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   6,  -4,   0,   0,   0,   0,   0,   1,     0,    -3,    -2,     0, &
         0,   0,   0,   0,   0,   0,   0,   6,   0,   0,   0,   0,   0,   2,     2,     9,     4,     3, &
         0,   0,   0,   0,   0,   0,   0,   6,   0,   0,   0,   0,   0,   0,     0,     0,     0,    -4, &
         0,   0,   0,   0,   0,   0,   0,   6,   0,   0,   0,   0,   0,   0,     8,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   6,   0,   0,   0,   0,   0,   1,     0,     4,     2,     0, &
         0,   0,   0,   0,   0,   0,   0,   6,   0,   0,   0,   0,   0,   2,     6,     0,     0,    -3, &
         0,   0,   0,   0,   0,   0,   0,   6,   0,   0,   0,   0,   0,   0,     6,     0,     0,     0, &
         0,   0,   0,   0,   0,   0,   0,   6,   0,   0,   0,   0,   0,   1,     0,     3,     1,     0, &
         0,   0,   0,   0,   0,   0,   0,   6,   0,   0,   0,   0,   0,   2,     5,     0,     0,    -2 /
    
    data IX63/   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,     3,     0,     0,    -1, &
         1,   0,   0,  -2,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   0,    -3,     0,     0,     0, &
         1,   0,   0,  -2,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,     6,     0,     0,     0, &
         1,   0,   0,  -2,   0,   0,   0,   1,   0,  -1,   0,   0,   0,   0,     7,     0,     0,     0, &
         1,   0,   0,  -2,   0,   0,   1,  -1,   0,   0,   0,   0,   0,   0,    -4,     0,     0,     0, &
         -1,   0,   0,   0,   0,   0,   3,  -3,   0,   0,   0,   0,   0,   0,     4,     0,     0,     0, &
         -1,   0,   0,   0,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   0,     6,     0,     0,     0, &
         -1,   0,   0,   2,   0,   0,   0,   4,  -8,   3,   0,   0,   0,   0,     0,    -4,     0,     0, &
         1,   0,   0,  -2,   0,   0,   0,   4,  -8,   3,   0,   0,   0,   0,     0,    -4,     0,     0, &
         -2,   0,   0,   2,   0,   0,   0,   4,  -8,   3,   0,   0,   0,   0,     5,     0,     0,     0 /
    
    data IX64/  -1,   0,   0,   0,   0,   0,   0,   2,   0,  -3,   0,   0,   0,   0,    -3,     0,     0,     0, &
         -1,   0,   0,   0,   0,   0,   0,   1,   0,  -1,   0,   0,   0,   0,     4,     0,     0,     0, &
         -1,   0,   0,   0,   0,   0,   1,  -1,   0,   0,   0,   0,   0,   0,    -5,     0,     0,     0, &
         -1,   0,   0,   2,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,     4,     0,     0,     0, &
         1,   0,  -1,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,     0,     3,     0,     0, &
         -1,   0,   0,   2,   0,   0,   0,   2,   0,  -3,   0,   0,   0,   0,    13,     0,     0,     0, &
         -2,   0,   0,   0,   0,   0,   0,   2,   0,  -3,   0,   0,   0,   0,    21,    11,     0,     0, &
         1,   0,   0,   0,   0,   0,   0,   4,  -8,   3,   0,   0,   0,   0,     0,    -5,     0,     0, &
         -1,   0,   1,  -1,   1,   0,   0,  -1,   0,   0,   0,   0,   0,   0,     0,    -5,    -2,     0, &
         1,   0,   1,  -1,   1,   0,   0,  -1,   0,   0,   0,   0,   0,   0,     0,     5,     3,     0 /
    
    data IX65/  -1,   0,   0,   0,   0,   0,   0,   4,  -8,   3,   0,   0,   0,   0,     0,    -5,     0,     0, &
         -1,   0,   0,   2,   1,   0,   0,   2,   0,  -2,   0,   0,   0,   0,    -3,     0,     0,     2, &
         0,   0,   0,   0,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   0,    20,    10,     0,     0, &
         -1,   0,   0,   2,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   0,   -34,     0,     0,     0, &
         -1,   0,   0,   2,   0,   0,   3,  -3,   0,   0,   0,   0,   0,   0,   -19,     0,     0,     0, &
         1,   0,   0,  -2,   1,   0,   0,  -2,   0,   2,   0,   0,   0,   0,     3,     0,     0,    -2, &
         1,   0,   2,  -2,   2,   0,  -3,   3,   0,   0,   0,   0,   0,   0,    -3,     0,     0,     1, &
         1,   0,   2,  -2,   2,   0,   0,  -2,   0,   2,   0,   0,   0,   0,    -6,     0,     0,     3, &
         1,   0,   0,   0,   0,   0,   1,  -1,   0,   0,   0,   0,   0,   0,    -4,     0,     0,     0, &
         1,   0,   0,   0,   0,   0,   0,   1,   0,  -1,   0,   0,   0,   0,     3,     0,     0,     0 /
    
    data IX66/   0,   0,   0,  -2,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,     3,     0,     0,     0, &
         0,   0,   0,  -2,   0,   0,   0,   1,   0,  -1,   0,   0,   0,   0,     4,     0,     0,     0, &
         0,   0,   2,   0,   2,   0,  -2,   2,   0,   0,   0,   0,   0,   0,     3,     0,     0,    -1, &
         0,   0,   2,   0,   2,   0,   0,  -1,   0,   1,   0,   0,   0,   0,     6,     0,     0,    -3, &
         0,   0,   2,   0,   2,   0,  -1,   1,   0,   0,   0,   0,   0,   0,    -8,     0,     0,     3, &
         0,   0,   2,   0,   2,   0,  -2,   3,   0,   0,   0,   0,   0,   0,     0,     3,     1,     0, &
         0,   0,   0,   2,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   0,    -3,     0,     0,     0, &
         0,   0,   1,   1,   2,   0,   0,   1,   0,   0,   0,   0,   0,   0,     0,    -3,    -2,     0, &
         1,   0,   2,   0,   2,   0,   0,   1,   0,   0,   0,   0,   0,   0,   126,   -63,   -27,   -55, &
         -1,   0,   2,   0,   2,   0,  10,  -3,   0,   0,   0,   0,   0,   0,    -5,     0,     1,     2 /
    
    data IX67/   0,   0,   1,   1,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,    -3,    28,    15,     2, &
         1,   0,   2,   0,   2,   0,   0,   1,   0,   0,   0,   0,   0,   0,     5,     0,     1,    -2, &
         0,   0,   2,   0,   2,   0,   0,   4,  -8,   3,   0,   0,   0,   0,     0,     9,     4,     1, &
         0,   0,   2,   0,   2,   0,   0,  -4,   8,  -3,   0,   0,   0,   0,     0,     9,     4,    -1, &
         -1,   0,   2,   0,   2,   0,   0,  -4,   8,  -3,   0,   0,   0,   0,  -126,   -63,   -27,    55, &
         2,   0,   2,  -2,   2,   0,   0,  -2,   0,   3,   0,   0,   0,   0,     3,     0,     0,    -1, &
         1,   0,   2,   0,   1,   0,   0,  -2,   0,   3,   0,   0,   0,   0,    21,   -11,    -6,   -11, &
         0,   0,   1,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,     0,    -4,     0,     0, &
         -1,   0,   2,   0,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   -21,   -11,    -6,    11, &
         -2,   0,   2,   2,   2,   0,   0,   2,   0,  -2,   0,   0,   0,   0,    -3,     0,     0,     1 /
    
    data IX68/   0,   0,   2,   0,   2,   0,   2,  -3,   0,   0,   0,   0,   0,   0,     0,     3,     1,     0, &
         0,   0,   2,   0,   2,   0,   1,  -1,   0,   0,   0,   0,   0,   0,     8,     0,     0,    -4, &
         0,   0,   2,   0,   2,   0,   0,   1,   0,  -1,   0,   0,   0,   0,    -6,     0,     0,     3, &
         0,   0,   2,   0,   2,   0,   2,  -2,   0,   0,   0,   0,   0,   0,    -3,     0,     0,     1, &
         -1,   0,   2,   2,   2,   0,   0,  -1,   0,   1,   0,   0,   0,   0,     3,     0,     0,    -1, &
         1,   0,   2,   0,   2,   0,  -1,   1,   0,   0,   0,   0,   0,   0,    -3,     0,     0,     1, &
         -1,   0,   2,   2,   2,   0,   0,   2,   0,  -3,   0,   0,   0,   0,    -5,     0,     0,     2, &
         2,   0,   2,   0,   2,   0,   0,   2,   0,  -3,   0,   0,   0,   0,    24,   -12,    -5,   -11, &
         1,   0,   2,   0,   2,   0,   0,  -4,   8,  -3,   0,   0,   0,   0,     0,     3,     1,     0, &
         1,   0,   2,   0,   2,   0,   0,   4,  -8,   3,   0,   0,   0,   0,     0,     3,     1,     0 /
    
    data IX69/   1,   0,   1,   1,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,     0,     3,     2,     0, &
         0,   0,   2,   0,   2,   0,   0,   1,   0,   0,   0,   0,   0,   0,   -24,   -12,    -5,    10, &
         2,   0,   2,   0,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,     4,     0,    -1,    -2, &
         -1,   0,   2,   2,   2,   0,   0,   2,   0,  -2,   0,   0,   0,   0,    13,     0,     0,    -6, &
         -1,   0,   2,   2,   2,   0,   3,  -3,   0,   0,   0,   0,   0,   0,     7,     0,     0,    -3, &
         1,   0,   2,   0,   2,   0,   1,  -1,   0,   0,   0,   0,   0,   0,     3,     0,     0,    -1, &
         0,   0,   2,   2,   2,   0,   0,   2,   0,  -2,   0,   0,   0,   0,     3,     0,     0,    -1    /
    
    
    !****  Initialize the values and sum over the series
    
    dpsi = 0.0d0
    deps = 0.0d0
    
    do i = num_plan, 1, -1
       
       !         Sum the mulitpliers by the arguments to the argument of
       !         nutation
       arg = 0.d0
       dargdt = 0.d0
       do j = 1,14
          
          !            Planetary values
          arg = arg + Plan_int(j,i)*plan_arg(j)
          dargdt = dargdt + Plan_int(j,i)*plan_rat(j)
       end do
       
       arg = mod(arg, 2.d0*pi)
       carg = cos(arg)
       sarg = sin(arg)
       ! MOD TAH 010819: Removed the d4 from the period since rates
       !         are now in rads/yr rather than rads/1.d4 yrs. 
       !         Output period is in Julian days.
       period = (2*pi/dargdt)*365.25d0
       
       !****      Now add contributions to dpsi and deps
       dpsi = dpsi + (Plan_int(15,i)*sarg + Plan_int(16,i)*carg) * 1.d-4
       deps = deps + (Plan_int(17,i)*sarg + Plan_int(18,i)*carg) * 1.d-4
       
       if( out(1:3).eq.'YES' ) then
          amp = sqrt( ( dble(Plan_int(15,i)**2 + Plan_int(16,i)**2)*0.4d0**2  +  Plan_int(17,i)**2 + Plan_int(18,i)**2) * 1.d-8 )
          write(*,175) i, (Plan_int(j,i), j=1,14), period, (Plan_int(j,i)*1.d-4, j=15,18), amp
175       format(i4,14(1x,I3),1x,F11.2,1x,2(F8.4,1x,F8.4,2x), 1x,F8.4)
       end if
       
    end do
    
  end subroutine eval_plan_nut
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  !TITLE 'PLAN_ANGLES'
  
  !*********************************************************************************************************************************
  subroutine plan_angles( epoch, plan_arg, plan_rat )
    use SUFR_kinds, only: double
    implicit none
    
    
    !     Routine to compute of planetary arguments for planetary
    !     nutation.  The longitudes of the major planets is computed
    !     according to:
    ! 
    ! Arguments based on Souchay, Loysel, Kinoshita, Folgueira, Corrections
    !    and new developments in rigid earth nutation theory, Aston. Astrophys.
    !    Suppl. Ser, 135, 111-131, (1999)
    !     
    
    
    
    ! PHYSICAL CONSTANTS NEEDED FOR SD_COMP
    
    !   pi          - Define here to full precision
    !   DJ2000      - Julian date of J2000
    
    !real(double), parameter :: pi = 3.1415926535897932d0
    real(double), parameter :: DJ2000 = 2451545.d0
    
    
    ! PASSED VARIABLES
    
    ! INPUT
    ! epoch  - Julian date plus fraction of a day
    !
    ! OUTPUT
    ! plan_arg(14)  - Planetary arguments for longtitudes of
    !           L, L', F, D, Om, Mercury, Venus, Earth, Mars, 
    !           Jupiter, Saturn, Uranus, Uranus(?), pa.
    ! plan_rat(14)  - Planetary argument rates (rads/year)
    
    real(double), intent(in) :: epoch
    real(double), intent(out) :: plan_arg(14), plan_rat(14)
    
    ! LOCAL VARIABLES 
    !      cent         - Centuries since J2000.
    !      nl           - Mecurcy longitude (rads)
    !      nlc(2)       - coefficients for computing nl
    !      vl           - Venus longitude (rads)
    !      vlc(2)       - coefficients for computing vl
    !      tl           - Earth longitude (rads)
    !      tlc(2)       - coefficients for computing tl
    !      ml           - Mars  longitude (rads)
    !      mlc(2)       - coefficients for computing ml
    !      jl           - Jupliter longitude (rads)
    !      jlc(2)       - coefficients for computing jl
    !      sl           - Saturn longitude (rads)
    !      slc(2)       - coefficients for computing sl
    !      ul           - Uranus longitude (rads)
    !      ulc(2)       - coefficients for computing ul
    !      xl           - Neptune longitude (rads)
    !      xlc(2)       - coefficients for compulting xl
    !      lm           - Mean longitude of moon minus mean longitude
    !                     of perigee (rads)
    !      lmc(2)       - Coefficients for computing lm
    !      sl           - Sun longitude (rads)
    !      slc(2)       - coefficients for computing sl
    !      Fr           - Moon's mean longitude minus Om (rad)
    !      Frc(2)       - Coefficients for computing Fr
    !      Dr           - Mean elongation of the Moon from the Sun (rads)
    !      drc(2)       - Coefficients for computing dc
    !      Om           - Longitude of the ascending node of the moon's
    !                     mean orbit on the elliptic (rad)
    !      Omc(2)       - Coefficients for computing Om
    !      pa           - pa (rads)
    !      pac(2)       - Coefficients for computing pa from KS1990.
    !                     (Values converted from rates and acceleration
    !                      in 1000's year to centuries).
    
    real(double) :: cent, nl, nlc(2), vl, vlc(2), tl, tlc(2), ml, mlc(2), &
         jl, jlc(2), sl, slc(2), ul, ulc(2), xl, xlc(2), &
         pa, pac(2), Dr, drc(2), Fr, frc(2),  &
         lm, lmc(2), ls, lsc(2), Om, Omc(2)
    
    data nlc / 4.402608842d0, 2608.7903141574d0 /
    data vlc / 3.176146697d0, 1021.3285546211d0 /
    data tlc / 1.753470314d0,  628.3075849991d0 /
    data mlc / 6.203480913d0,  334.0612426700d0 /
    data jlc / 0.599546497d0,   52.9690962641d0 /
    data slc / 0.874016757d0,   21.3299104960d0 /
    data ulc / 5.481293871d0,    7.4781598567d0 /
    data xlc / 5.321159000d0,    3.8127774000d0 /
    data lmc / 2.355555980d0, 8328.6914269554d0 /
    data lsc / 6.240060130d0,  628.301955d0 /
    data frc / 1.627905234d0, 8433.466158131d0 /
    data drc / 5.198466741d0, 7771.3771468121d0 /
    data omc / 2.182439200d0,  -33.757045d0 /
    data pac / 0.2438175d-1,       0.000538691d-2 /
    
    !***** Get number of Centuries since J2000
    
    cent = (epoch-DJ2000) / 36525.d0
    
    !     Compute arguments 
    nl  = nlc(1) + nlc(2)*cent
    vl  = vlc(1) + vlc(2)*cent
    tl  = tlc(1) + tlc(2)*cent
    ml  = mlc(1) + mlc(2)*cent
    jl  = jlc(1) + jlc(2)*cent
    sl  = slc(1) + slc(2)*cent
    ul  = ulc(1) + ulc(2)*cent
    xl  = xlc(1) + xlc(2)*cent
    lm  = lmc(1) + lmc(2)*cent
    ls  = lsc(1) + lsc(2)*cent
    fr  = frc(1) + frc(2)*cent
    dr  = drc(1) + drc(2)*cent
    om  = omc(1) + omc(2)*cent
    pa  = pac(1)*cent + pac(2)*cent**2
    
    !****  Now save the values
    plan_arg( 1) = lm
    plan_arg( 2) = ls
    plan_arg( 3) = fr 
    plan_arg( 4) = dr 
    plan_arg( 5) = om 
    plan_arg( 6) = nl
    plan_arg( 7) = vl
    plan_arg( 8) = tl
    plan_arg( 9) = ml
    plan_arg(10) = jl 
    plan_arg(11) = sl 
    plan_arg(12) = ul 
    plan_arg(13) = xl 
    plan_arg(14) = pa
    
    !****  Now save the rates (useful if periods are wanted)
    ! MOD TAH 010819: Changed *100.d0 to /100.d0 to make rates
    !     in radians/year rather than radians/10000 year.  Also
    !     need to correct period calculation in eval_plan_nut
    plan_rat( 1) = lmc(2)/100.d0
    plan_rat( 2) = lsc(2)/100.d0
    plan_rat( 3) = frc(2)/100.d0
    plan_rat( 4) = drc(2)/100.d0
    plan_rat( 5) = omc(2)/100.d0
    plan_rat( 6) = nlc(2)/100.d0
    plan_rat( 7) = vlc(2)/100.d0
    plan_rat( 8) = tlc(2)/100.d0
    plan_rat( 9) = mlc(2)/100.d0
    plan_rat(10) = jlc(2)/100.d0
    plan_rat(11) = slc(2)/100.d0
    plan_rat(12) = ulc(2)/100.d0
    plan_rat(13) = xlc(2)/100.d0
    plan_rat(14) = pac(1)/100.d0 + 2*pac(2)*cent/100.d0
    
    
  end subroutine plan_angles
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  
  
  !TITLE FCN_NUT
  
  !*********************************************************************************************************************************
  subroutine fcn_nut ( jd, dpsi_fcn, deps_fcn )
    use SUFR_kinds, only: double
    implicit none
    
    !     Routine to compute the consttributions of the freely excited
    !     FCN mode to the nutations in longitude and obliquity.
    
    ! USAGE:
    !     call fcn_nut( jd, dpsi_fcn, deps_fcn )
    !     where <jd>    is a full julian date with fractional part
    !                   of the day added (REAL(DOUBLE) INPUT)
    !     and <dpsi_fcn> and <deps_fcn> are the contributions to the nutations
    !                   in longitude and obliquity in milliarcsec.
    !                   (REAL(DOUBLE) OUTPUT)
    
    ! RESTRICTIONS: if <jd> is less than 2000000.0 this routine
    !               assumes an MJD has been passed and the time
    !               used will be converted to JD.  A warning
    !               message will be printed.
    
    ! RESTRICTIONS: This term represents as free excitation mode and
    !               therefore will change with time (in much the same
    !               way that the Chandler Wobble changes).  The
    !               frequency of the FCN used here is accurate, but
    !               coefficients used will depend on time.  The values
    !               are interpolated over the 1979-2000 interval with
    !               the first or last values being used outside of these
    !               times.
    
    ! PARAMETERS:
    
    !   DJ2000      - Julian date of J2000
    !   solar_to_sidereal   - Conversion from solar days to sidereal
    !                 days.
    !   num_fcn     - Number of FCN amplitudes (linear interpolation between
    !                 values).
    
    
    real(double) :: DJ2000, pi, solar_to_sidereal
    
    integer :: num_fcn
    
    parameter ( pi            = 3.1415926535897932D0 )
    parameter ( DJ2000        = 2451545.d0           )
    parameter ( solar_to_sidereal = 1.00273790935d0 )
    parameter ( num_fcn       = 11 )
    
    ! PASSED VARIABLES
    !
    ! INPUT Values
    ! jd     - Time at which value needed. (jd + fraction of day)
    
    ! OUTPUT Values
    ! dpsi_fcn   - Contribution to the nutation in longitude (mas).  Should
    !          be added to standard nutation in longitude values.
    ! deps_fcn   - Contribution to the nutation in obliquity (mas).  Should
    !          be added to standard nutation in obliquity values.
    
    
    real(double), intent(in) :: jd
    real(double), intent(out) :: dpsi_fcn, deps_fcn
    
    ! LOCAL VARIABLES
    
    !   epoch       - Julian date (jd passed in unless the JD
    !                 appears to be an MJD in which case it is
    !                 converted to JD (2 400 000.5d0 added)
    
    !   fcn_freq    - Freqency of the FCN mode (cycles per sidreal day).
    !   fcn_arg     - Argument for the fcn mode computed from J2000 (rad).
    !   fcn_tabl(2, num_fcn)  - Amplitude for the fcn free exciation.  These
    !                 are converted to nutation in longitude and obliquity.  (mas)
    !   fcn_ampl(2) - Interpolated FCN amplitude (mas)
    !   fcn_jd(num_fcn) - Starting epochs for the fcn amplitudes (JD).  These
    !                 epochs are the midpoints of the intervals overwhich the
    !                 FCN free terms have been estimated.
    !   sine        - Sine of the mean obliquity of the ecliptic.  (A constant
    !                 value can be used here since the changes are small for
    !                 this constribtion i.e., between 1980 and 2000 the error
    !                 in the nutation in longitude is only 0.05 micro-arc-sec.
    !                 Values based on 8438.14059" at J2000 (IERS 2000 Conventions)
    !   dt          - Time difference between epoch and tabular interval (days)
    !   dt_tab      - Time difference in tables values.
    !   dfcn_amp(2) - Change in FCN amplitude between tabular points (mas)
    
    
    real(double) :: epoch, fcn_freq, fcn_arg, fcn_ampl(2), fcn_tabl(2,num_fcn), sine, fcn_jd(num_fcn), dt, dt_tab, dfcn_amp(2)
    
    !   i           - A counter used in do loop to find the correct pair of
    !                 amplitudes to use
    
    integer :: i
    
    data  fcn_freq  /   -1.00231810920d0 /
    
    ! FCN estimated values (by time range, peicewise linear function)
    ! Node date     Cos       +-       Sin        +- 
    ! 1979/ 1/ 1  -0.0620   0.1256    -0.1346   0.1293   mas 
    ! 1984/ 1/ 1   0.0447   0.0302    -0.1679   0.0309   mas 
    ! 1986/ 1/ 1   0.2406   0.0163    -0.2759   0.0159   mas 
    ! 1988/ 1/ 1   0.1183   0.0127    -0.2163   0.0128   mas 
    ! 1990/ 1/ 1   0.0479   0.0084    -0.1965   0.0083   mas 
    ! 1992/ 1/ 1  -0.0796   0.0071    -0.1321   0.0071   mas 
    ! 1994/ 1/ 1  -0.0075   0.0057    -0.1150   0.0057   mas 
    ! 1996/ 1/ 1  -0.0128   0.0058    -0.0998   0.0058   mas 
    ! 1998/ 1/ 1  -0.0263   0.0059    -0.1122   0.0059   mas 
    ! 2000/ 1/ 1   0.0519   0.0071     0.0081   0.0070   mas 
    ! 2001/ 6/ 1   0.2100   0.0162     0.1401   0.0163   mas 
    
    ! These are the estimated values with their standard deviations.
    ! Due to the large standard deviation of thr 1979 node, we replace
    ! its value with the 1984 value, thus maintaining a constant amplitude
    ! between 1979 and 1984. 
    
    !     Time dependent values 
    data  fcn_jd   / 2443874.5d0,  2445700.5d0, 2446431.5d0, 2447161.5d0,  2447892.5d0, 2448622.5d0, &
         2449353.5d0,  2450083.5d0, 2450814.5d0, 2451544.5d0,  2452061.5d0 /
    
    
    data  fcn_tabl /     -0.062d0,    -0.135d0,   0.045d0,    -0.168d0,   0.241d0,    -0.276d0,   0.118d0,    -0.216d0,   &
         0.048d0,    -0.197d0,   -0.080d0,    -0.132d0,  -0.007d0,    -0.115d0,   -0.013d0,    -0.100d0,   &
         -0.026d0,    -0.112d0,   0.052d0,     0.008d0,   0.210d0,     0.140d0 /  
    
    data  sine      /   0.3977769687d0 /
    
    fcn_ampl(1:2) = 0.d0
    
    !***** Check to make sure user passed JD and not MJD.  Correct
    !     problem and warn the user.
    ! MvdS: remove this 'solution'
    !if( jd .lt.2000000.0d0  ) then
    !          write(*,100) jd
    ! 100      format('**WARNING** MJD apparently passed to FCN_NUT',
    !     .          ' Value (',F10.2,') converted to JD')
    !   epoch = jd + 2 400 000.5d0
    !else
    !   epoch = jd
    !end if
    epoch = jd
    
    !****  Find out which table values we should use.
    if( epoch.le.fcn_jd(1) ) then
       fcn_ampl(1) = fcn_tabl(1,1)
       fcn_ampl(2) = fcn_tabl(2,1)
    else if( epoch.ge. fcn_jd(num_fcn) ) then
       fcn_ampl(1) = fcn_tabl(1,num_fcn)
       fcn_ampl(2) = fcn_tabl(2,num_fcn)
    else
       do i = 1, num_fcn-1
          if( epoch.ge.fcn_jd(i) .and. epoch.lt.fcn_jd(i+1) ) then
             dt = epoch - fcn_jd(i)
             dt_tab = fcn_jd(i+1) - fcn_jd(i)
             dfcn_amp(1) = fcn_tabl(1,i+1) - fcn_tabl(1,i)
             dfcn_amp(2) = fcn_tabl(2,i+1) - fcn_tabl(2,i)
             fcn_ampl(1) = fcn_tabl(1,i) + (dfcn_amp(1)/dt_tab)*dt
             fcn_ampl(2) = fcn_tabl(2,i) + (dfcn_amp(2)/dt_tab)*dt
          end if
       end do
    end if
    
    
    !***** Get the argument for the FCN mode at this times
    
    fcn_arg = -2*pi*(1.d0+fcn_freq)*solar_to_sidereal*(epoch-DJ2000)
    
    dpsi_fcn = (-fcn_ampl(1)*sin(fcn_arg) +  fcn_ampl(2)*cos(fcn_arg))/sine
    
    deps_fcn = (-fcn_ampl(1)*cos(fcn_arg) -  fcn_ampl(2)*sin(fcn_arg))
    
  end subroutine fcn_nut
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  
  
  !TITLE PREC_NUT
  
  !*********************************************************************************************************************************
  subroutine prec_nut( jd, dpsi_prec, deps_prec )
    use SUFR_kinds, only: double
    implicit none
    
    !     Routine to evaluate the corrections to the nutations in longitude
    !     and obliquity due to the corrections to the IAU-1976 Luni-solar
    !     precession constant and the secular rate of change of the obliquity
    !     of the ecliptic.
    
    ! PARAMETERS:
    
    !   DJ2000      - Julian date of J2000
    
    
    real(double) :: DJ2000
    
    parameter ( DJ2000        = 2451545.d0           )
    
    ! PASSED VARIABLES
    !
    ! INPUT Values
    ! jd     - Time at which value needed. (jd + fraction of day)
    
    ! OUTPUT Values
    ! dpsi_prec   - Contribution to the nutation in longitude (mas).  Should
    !          be added to standard nutation in longitude values. Value
    !          valid only when the IAU-1976 precession constant used to
    !          compute the transformation to mean system.
    ! deps_prec   - Contribution to the nutation in obliquity (mas).  Should
    !          be added to standard nutation in obliquity values.
    
    
    real(double), intent(in) :: jd
    real(double), intent(out) :: dpsi_prec, deps_prec
    
    ! LOCAL VARIABLES
    
    !   epoch       - Julian date (jd passed in unless the JD
    !                 appears to be an MJD in which case it is
    !                 converted to JD (2 400 000.5d0 added)
    !   cent        - Number of Julian centuries since J2000.0
    !   DpsiDt      - Correction to precession constant as a
    !                 linear rate of change of nutation in
    !                 longitude. (arc-second/century)
    !   DepsDt      - Correction to rate of change of oblquity
    !                 (arc-second/century)
    
    
    
    real(double) :: epoch, cent, DpsiDt,  DepsDt
    !
    !     Theoretical estimaate of adjustment to precession constant is
    !     -0.29965 "/cent; estimate is -0.29965 +- 0.0004 "/cent.
    
    data  DpsiDt  /  -0.29965d0  /
    data  DepsDt  /  -0.02524d0  /
    
    !***** Check to make sure user passed JD and not MJD.  Correct
    !     problem and warn the user.
    ! MvdS: remove this 'solution'
    !if( jd .lt.2000000.0d0  ) then
    !             write(*,100) jd
    !    100      format('**WARNING** MJD apparently passed to SD_COMP',
    !        .          ' Value (',F10.2,') converted to JD')
    !   epoch = jd + 2 400 000.5d0
    !else
    !   epoch = jd
    !end if
    epoch = jd
    
    !****  Compute the number of centuries
    
    cent = (epoch - DJ2000)/36525.d0
    
    dpsi_prec = DpsiDt*cent*1000.d0
    deps_prec = DepsDt*cent*1000.d0
    
  end subroutine prec_nut
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  
  
  
  !TITLE ls_iau76
  
  !*********************************************************************************************************************************
  subroutine ls_iau76( epoch, ls_arg )
    use SUFR_kinds, only: double
    implicit none
    
    !     Routine to compute the value of the fundamental argument
    !     for Brown's arguments.  Arguments based on the IERS
    !     standards.
    
    ! PHYSICAL CONSTANTS
    
    !   pi          - Define here to full precision
    !   rad_to_deg  - Conversion from radians to degs.
    !   DJ2000      - Julian date of J2000
    !   sec360      - number of seconds in 360 degreees.
    
    
    real(double) :: pi, rad_to_deg, DJ2000, sec360
    
    parameter ( pi            = 3.1415926535897932D0 )
    parameter ( DJ2000        = 2451545.d0           )
    parameter ( sec360        = 1296000.d0           )
    
    !     Computed quanities
    parameter ( rad_to_deg    = 180.d0   /pi         )
    
    !-------------------------------------------------------------------
    
    ! PASSED VARIABLES
    
    ! INPUT
    ! epoch  - Julian date for arguments (jd + fraction of day, REAL(DOUBLE))
    
    ! OUTPUT
    ! ls_arg(5) -  Brown's arguments (radians, REAL(DOUBLE))
    
    
    real(double), intent(in) :: epoch
    real(double), intent(out) :: ls_arg(5)
    
    ! LOCAL VARIABLES
    !      cent             - Julian centuries to DJ2000.
    !      el,eld           - Mean longitude of moon minus mean
    !                       - longitude of moon's perigee (arcsec)
    !      elc(5)           - Coefficients for computing el
    !      elp,elpd         - Mean longitude of the sun minus mean
    !                       - longitude of sun perigee (arcsec)
    !      elpc(5)          - Coeffiecents for computing elp
    !      f,fd             - Moon's mean longitude minus omega (sec)
    !      fc(5)            - Coefficients for computing f
    !      d,dd             - Mean elongation of the moon from the
    !                       - sun (arcsec)
    !      dc(5)            - coefficients for computing d
    !      om,omd           - longitude of the ascending node of the
    !                       - moon's mean orbit on the elliptic
    !                       - measured from the mean equinox of date
    !      omc(5)           - Coefficients for computing om.
    
    
    real(double) :: cent, el, elc(5), elp, elpc(5),  f, fc(5), d, dc(5), om, omc(5)  ! ,eld , elpd ,fd ,dd ,omd
    
    !****  DATA statements for the fundamental arguments.
    
    data elc    /     0.064d0,    31.310d0,    715922.633d0,  485866.733d0,     1325.0d0 /
    data elpc   /    -0.012d0,    -0.577d0,   1292581.224d0,  1287099.804d0,      99.0d0 /
    data fc     /     0.011d0,   -13.257d0,    295263.137d0,  335778.877d0,     1342.0d0 /
    data dc     /     0.019d0,    -6.891d0,   1105601.328d0,  1072261.307d0,    1236.0d0 /
    data omc    /     0.008d0,     7.455d0,   -482890.539d0,  450160.280d0,       -5.0d0 /
    
    !****  Get the number of centuries to current time
    
    cent = (epoch-dj2000) / 36525.d0
    
    !****  Compute angular arguments
    el = elc(1) * cent**3 + elc(2) * cent**2 + elc(3) * cent   + elc(4) + mod( elc(5) * cent, 1.d0 ) * sec360
    el = mod( el, sec360 )
    !eld = 3.d0 * elc(1) * cent**2 + 2.d0 * elc(2) * cent + elc(3)   + elc(5) * sec360
    
    elp = elpc(1) * cent**3 + elpc(2) * cent**2 + elpc(3) * cent   + elpc(4) + mod( elpc(5) * cent, 1.d0 ) * sec360
    elp = mod( elp, sec360 )
    !elpd = 3.d0 * elpc(1) * cent**2 + 2.d0 * elpc(2) * cent + elpc(3)   + elpc(5) * sec360
    
    f = fc(1) * cent**3 + fc(2) * cent**2 + fc(3) * cent   + fc(4) + mod( fc(5) * cent, 1.d0 ) * sec360
    f = mod( f, sec360 )
    !fd = 3.d0 * fc(1) * cent**2 + 2.d0 * fc(2) * cent + fc(3)   + fc(5) * sec360
    
    d = dc(1) * cent**3 + dc(2) * cent**2 + dc(3) * cent   + dc(4) + mod( dc(5) * cent, 1.d0 ) * sec360
    d = mod( d, sec360 )
    !dd = 3.d0 * dc(1) * cent**2 + 2.d0 * dc(2) * cent + dc(3)   + dc(5) * sec360
    
    om = omc(1) * cent**3 + omc(2) * cent**2 + omc(3) * cent   + omc(4) + mod( omc(5) * cent, 1.d0 ) * sec360
    om = mod( om, sec360 )
    !omd = 3.d0 * omc(1) * cent**2 + 2.d0 * omc(2) * cent + omc(3)   + omc(5) * sec360
    
    
    !****  Now save the values.  Convert values from arcseconds to radians
    
    ls_arg(1) = el / (3600.d0*rad_to_deg)
    ls_arg(2) = elp/ (3600.d0*rad_to_deg)
    ls_arg(3) = f  / (3600.d0*rad_to_deg)
    ls_arg(4) = d  / (3600.d0*rad_to_deg)
    ls_arg(5) = om / (3600.d0*rad_to_deg)
    
  end subroutine ls_iau76
  !*********************************************************************************************************************************
  
  
end module TheSky_nutation
!***********************************************************************************************************************************

