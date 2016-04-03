!> \file vsop.f90  Core VSOP87 routines for libTheSky


!  Copyright (c) 2002-2015  Marc van der Sluys - marc.vandersluys.nl
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
!> \brief Procedures for VSOP87
!!
!! \note
!! Accuracies 1900-2100 in milliarcseconds (VSOP87 paper):
!! - Me: 1,  Ve: 6,  Ea: 5,  Ma: 23,  Ju: 20,  Sa: 100,  Ur: 16,  Ne: 30
!!
!! \note
!! Accuracy < 1" in JD2000 +- range in kyr (Wikipedia):
!! - Me-Ma: 4, Ju,Sa: 2, Ur,Ne: 6  (i.e., Earth position more accurate than 1" between -2000 and +6000)
!!
!! \see Bretagnon, P.; Francou, G. Planetary theories in rectangular and spherical variables - VSOP 87 solutions 
!!      http://esoads.eso.org/abs/1988A%26A...202..309B

module TheSky_vsop
  implicit none
  save
  
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate true heliocentric ecliptic coordinates l,b,r for planet pl and the mean ecliptic and equinox of date 
  !!         using VSOP87D
  !!
  !! \param  tm   Dynamical time in Julian Millennia after J2000.0, (tau in Meeus)
  !! \param  pl   Planet to compute position for
  !!
  !! \retval lon  Heliocentric longitude (rad)
  !! \retval lat  Heliocentric latitude (rad)
  !! \retval rad  Heliocentric distance (AU)
  !!
  !! \param  LBaccur  Desired accuracy of L,B (rad, optional)
  !! \param  Raccur   Desired accuracy of R (AU, optional)
  !!
  !! \see http://esoads.eso.org/abs/1988A%26A...202..309B
  
  subroutine vsop_lbr(tm,pl, lon,lat,rad, LBaccur,Raccur)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev
    use TheSky_planetdata, only: VSOPnls, VSOPdat, vsopNblk, VSOPtruncs
    
    implicit none
    real(double), intent(in) :: tm
    integer, intent(in) :: pl
    real(double), intent(out) :: lon,lat,rad
    real(double), intent(in), optional :: LBaccur,Raccur
    
    integer :: li, pow, Nli,Nle, var, nTerm, skip
    real(double) :: fac, lbr(3), accur, desired_accuracy(3)
    
    desired_accuracy = VSOPtruncs(1:3, pl)  ! Set accuracy equal to VSOP87 truncation
    if(present(LBaccur)) desired_accuracy(1:2) = (/LBaccur,LBaccur/)
    if(present(Raccur))  desired_accuracy(3)   = Raccur
    !desired_accuracy = 1.d-9  ! 5e-9 rad = 1 mas = VSOP87 accuracy for Mercury in 1900-2100
    !desired_accuracy = 0.d0  ! Use all available terms
    
    lbr = 0.d0
    !$omp parallel do private(Nli,Nle,skip,li,pow,fac,nTerm)
    do var=1,3  ! L,B,R
       Nli = 0
       if(var.ge.2) Nli = VSOPnls(1,pl)
       if(var.eq.3) Nli = Nli + VSOPnls(2,pl)
       Nle = Nli + VSOPnls(var,pl)
       
       skip = -1  ! Switch to skip terms if sufficient accuracy is reached
       do li=Nli+1,Nle
          pow = nint(VSOPdat(1,li,pl))
          
          if(pow.ne.skip) then
             fac = tm**pow * VSOPdat(2,li,pl)
             lbr(var) = lbr(var) + fac * cos( VSOPdat(3,li,pl) + VSOPdat(4,li,pl)*tm )
             
             ! Determine current accuracy (A*T^n):
             if(mod(li,3).eq.0) then  ! Reduce overhead, experimentally optimised to 3
                nTerm = li - vsopNblk(pow,1,pl)+1
                accur = 2*sqrt(dble(nTerm)) * abs(fac)
                if(accur.lt.desired_accuracy(var)) skip = pow  ! Skip the remaining terms for this power
             end if
          end if
          
       end do  ! li
    end do  ! var
    !$omp end parallel do
    
    lon = lbr(1)
    lat = lbr(2)
    rad = lbr(3)
    
  end subroutine vsop_lbr
  !*********************************************************************************************************************************
  
  
  
end module TheSky_vsop
!***********************************************************************************************************************************
