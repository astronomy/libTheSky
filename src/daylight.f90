!> \file daylight.f90  Procedures that deal with daylight for libTheSky


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
!> \brief Daylight procedures

module TheSky_daylight
  implicit none
  save
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute diffuse radiation on an inclined surface using the Perez 1987 model
  !!
  !! \param  DoY          Day of year (Nday)
  !! \param  alt          Altitude of the Sun (radians)
  !!
  !! \param  surfIncl     Surface inclination wrt horizontal (radians) - 0 = horizontal, pi/2 = vertical
  !! \param  theta        Angle between surface normal vector and Sun position vector (radians)
  !!
  !! \param  Gbeam_n      Beam (direct) normal radiation (W/m2; in the direction of the Sun)
  !! \param  Gdif_hor     Diffuse radiation on a horizontal surface (W/m2)
  !!
  !! \retval Gdif_inc     Diffuse irradiation on the inclined surface (W/m2)
  !!
  !! \retval Gdif_inc_cs  Diffuse irradiation on the inclined surface - circumsolar part (optional; W/m2)
  !! \retval Gdif_inc_hz  Diffuse irradiation on the inclined surface - horizon-band part (optional; W/m2)
  !!
  !! \see Perez et al. Solar Energy Vol. 39, Nr. 3, p. 221 (1987) - references to equations and tables are to this paper.
  !! Most equations can be found in the Nomenclature section at the end of the paper (p.230).  We use a and c here, not b and d.
  
  subroutine diffuse_radiation_Perez87(DoY, alt, surfIncl, theta, Gbeam_n,Gdif_hor,  Gdif_inc,   Gdif_inc_cs, Gdif_inc_hz)
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi2,pio2, d2r,r2d
    
    implicit none
    integer, intent(in) :: DoY
    real(double), intent(in) :: alt, surfIncl, theta, Gbeam_n,Gdif_hor
    real(double), intent(out) :: Gdif_inc
    real(double), intent(out), optional :: Gdif_inc_cs, Gdif_inc_hz
    integer :: f11,f12,f13, f21,f22,f23
    real(double) :: zeta, AM0rad, Mair,Delta,epsilon, alpha, psiC,psiH, chiC,chiH, F(6),F1,F2, A,C,  Gdif_inc_csl, Gdif_inc_hzl
    
    
    ! *** Compute the brightness coefficients for the circumsolar (F1) and horizon (F2) regions ***
    
    ! 'External' (AM0) radiation:
    AM0rad = 1370.d0 * (1.d0 + 0.00333d0 * cos(pi2/365.d0 * DoY))
    
    ! Air mass:
    if(alt .lt. -3.885d0*d2r) then
       Mair = 99.d0
    else if(alt .lt. 10.d0*d2r) then
       Mair = 1.d0 / ( sin(alt) + 0.15d0 * (alt*r2d + 3.885d0)**(-1.253d0) )
    else
       Mair = 1.d0 / sin(alt)
    end if
    Delta = Gdif_hor * Mair / AM0rad  ! Brightness of overcast sky - par. 2.2.4 (a)
    
    
    ! Cloudliness: epsilon;  epsilon ~ 1: overcast, epsilon -> infinity: clear  (epsilon ~ 1/fraction of covered sky)
    !   Needed for correct row in Table 1
    if(Gdif_hor.le.0.d0) then  ! Division by zero
       if(Gbeam_n.le.0.d0) then  ! No direct light: 0/0
          epsilon = 0.d0      !   -> completely overcast - first row of Table 1
       else                   ! Some direct daylight: x/0 = large
          epsilon = 99.d0     !   -> completely clear, should be >11 for last row of Table 1
       end if
    else
       epsilon = (Gdif_hor + Gbeam_n) / Gdif_hor  ! Overcast: epsilon ~ 1,  clear: epsilon -> infinity
    end if
    
    
    ! Table 1
    f11=1;  f12=2;  f13=3;  f21=4; f22=5; f23=6
    if(epsilon .le. 1.056d0) then
       F = [-0.011d0,  0.748d0, -0.080d0, -0.048d0,  0.073d0, -0.024d0]
    else if(epsilon .le. 1.253d0) then
       F = [-0.038d0,  1.115d0, -0.109d0, -0.023d0,  0.106d0, -0.037d0]
    else if(epsilon .le. 1.586d0) then
       F = [ 0.166d0,  0.909d0, -0.179d0,  0.062d0, -0.021d0, -0.050d0]
    else if(epsilon .le. 2.134d0) then
       F = [ 0.419d0,  0.646d0, -0.262d0,  0.140d0, -0.167d0, -0.042d0]
    else if(epsilon .le. 3.230d0) then
       F = [ 0.710d0,  0.025d0, -0.290d0,  0.243d0, -0.511d0, -0.004d0]
    else if(epsilon .le. 5.980d0) then
       F = [ 0.857d0, -0.370d0, -0.279d0,  0.267d0, -0.792d0,  0.076d0]
    else if(epsilon .le. 10.080d0) then
       F = [ 0.734d0, -0.073d0, -0.228d0,  0.231d0, -1.180d0,  0.199d0]
    else
       F = [ 0.421d0, -0.661d0,  0.097d0,  0.119d0, -2.125d0,  0.446d0]
    end if
    
    zeta = pio2 - alt  ! Zenith angle = pi/2 - alt
    F1 = F(f11)  +  F(f12) * Delta  +  F(f13) * zeta  ! Circumsolar brightness coefficient
    F2 = F(f21)  +  F(f22) * Delta  +  F(f23) * zeta  ! Horizon brightness coefficient
    
    
    
    
    ! *** Compute the mean solid angles occupied by the circumsolar region (C and A, needed for Eq.8) ***
    
    alpha = 25.d0*d2r  ! Half angle of the circumsolar region (degrees -> radians; below Eq.9)
    
    
    ! Solid angle of the circumsolar region weighted by incidence on the HORIZONTAL (variable C, subscript H; 
    !   see Nomenclature, under c):
    ! psiH:
    if(zeta .gt. pio2 - alpha) then
       psiH = 0.5d0 * (pio2 - zeta + alpha) / alpha  ! Dimensionless ratio
    else
       psiH = 1.d0
    end if
    
    ! chiH:
    if(zeta .lt. pio2 - alpha) then
       chiH = cos(zeta)  ! = sin(alt)
    else
       chiH = psiH * sin(psiH*alpha)
    end if
    
    C = 2 * (1.d0 - cos(alpha)) * chiH  ! Solid angle of the circumsolar region, weighted by HORIZONTAL incidence
    
    
    ! Solid angle of the circumsolar region weighted by incidence on the SLOPE (variable A, subscript C; 
    !   see Nomenclature, under c):
    ! psiC:
    psiC = 0.5d0 * (pio2 - theta + alpha) / alpha
    
    ! chiC:
    if(theta .lt. pio2 - alpha) then
       chiC = psiH * cos(theta)
    else if(theta .lt. pio2 + alpha) then
       chiC = psiH * psiC * sin(psiC*alpha)
    else
       chiC = 0.d0
    end if
    
    A = 2 * (1.d0 - cos(alpha)) * chiC  ! Solid angle of the circumsolar region, weighted by SLOPE incidence
    
    
    
    ! Diffuse radiation from circumsolar (F1) and horizon (F2) regions on the inclined surface (Eq.8):
    Gdif_inc_csl = Gdif_hor * ( 0.5d0 * (1.d0 + cos(surfIncl)) * (1.d0 - F1)  +  F1 * A/C )  ! Circumsolar
    Gdif_inc_hzl = Gdif_hor * ( F2 * sin(surfIncl) )                                         ! Horizon band
    
    Gdif_inc_csl = max(Gdif_inc_csl, 0.d0)  ! Components are sometimes negative
    Gdif_inc_hzl = max(Gdif_inc_hzl, 0.d0)
    Gdif_inc = Gdif_inc_csl + Gdif_inc_hzl
    
    ! Assign optional return values:
    if(present(Gdif_inc_cs)) Gdif_inc_cs = Gdif_inc_csl
    if(present(Gdif_inc_hz)) Gdif_inc_hz = Gdif_inc_hzl
    
  end subroutine diffuse_radiation_Perez87
  !*********************************************************************************************************************************
  
  
  
  
end module TheSky_daylight
!***********************************************************************************************************************************
