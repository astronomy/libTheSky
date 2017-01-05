!> \file modules.f90  Modules used by libTheSky


!  Copyright (c) 2002-2017  Marc van der Sluys - marc.vandersluys.nl
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
!> \brief  Constants used in libTheSky

module TheSky_constants
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  integer, parameter :: deltat_nmax = 1000   ! Need ~430 until 2000
  real(double) :: deltat_values(deltat_nmax), deltat_years(deltat_nmax)
  real(double) :: deltat_accel, deltat_change, deltat_0
  integer :: deltat_n, deltat_minyr, deltat_maxyr
  
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
  
  real(double) :: VSOPdat(4,6827,10), planpos(nplanpos), VSOPtruncs(3,8)
  real(double) :: plelems(8,6),plelems2000(8,6),plelemdata(2,8,6,0:3)
  real(double) :: asterelems(nasteroids,9)
  
  character :: plcon(0:19)*(3),asternames(nasteroids)*(18)
  
end module TheSky_planetdata
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Moon data, needed to compute Moon positions

module TheSky_moondata
  use SUFR_kinds, only: double
  implicit none
  save
  private :: double
  
  real(double), parameter :: c1=1.d0/60.d0, c2=1.d0/3600.d0
  real(double), parameter :: ath=384747.9806743165d0, a0=384747.9806448954d0
  
  real(double) :: p1,p2,p3,p4,p5,q1,q2,q3,q4,q5
  real(double) :: w(3,0:4),eart(0:4),peri(0:4),p(8,0:1), del(4,0:4),zeta(0:1),t(0:4)
  real(double) :: pre(3),coef(7),zone(6), pc1(6,1023),pc2(6,918),pc3(6,704), per1(3,19537),per2(3,6766),per3(3,8924)
  
  integer :: ilu(4),ipla(11),nterm(3,12),nrang(3,12)
  
  real(double) :: prec0
  integer :: ideb
  
  data ideb/0/,prec0/-1.d0/,t/1.d0,4*0.d0/
  
end module TheSky_moondata
!***********************************************************************************************************************************


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


