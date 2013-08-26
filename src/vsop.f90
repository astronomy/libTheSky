!> \file vsop.f90  Core VSOP87 routines for libTheSky


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


!***********************************************************************************************************************************
!> \brief Procedures for VSOP87

module TheSky_vsop
  implicit none
  save
  
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate l,b,r for planet pl using VSOP87  (true/apparent???)
  !!
  !! \param  tm   Dynamical time in Julian Millennia after 2000.0, the tau in Meeus
  !! \param  pl   Planet to compute position for
  !!
  !! \retval lon  Heliocentric longitude (rad)
  !! \retval lat  Heliocentric latitude (rad)
  !! \retval rad  Heliocentric distance (AU?)
  
  subroutine vsop_lbr(tm,pl, lon,lat,rad)
    use SUFR_kinds, only: double
    use SUFR_angles, only: rev
    use TheSky_planetdata, only: vsopdat
    
    implicit none
    real(double), intent(in) :: tm
    integer, intent(in) :: pl
    real(double), intent(out) :: lon,lat,rad
    
    integer :: nls(3,8), nl(3), li  !, omp_get_thread_num
    real(double) :: dat(4,6827)
    
    nls = reshape( (/ 2808,1620,2399, 671,426,585, 1080,349,997, 2393,915, &    ! Number of lines in VSOP input files (l,b,r x 8 pl)
         2175,1484,530,1469,2358,966,2435,1578,516,1897,681,290,959/), &
         (/3,8/))
    
    dat = vsopdat(:,:,pl)
    nl  = nls(:,pl)
    
    lon = 0.d0
    lat = 0.d0
    rad = 0.d0
    
    !call omp_set_num_threads(3)
    ! $omp parallel shared(pl,dat,nl,tm,vsopdat,nls,lon,lat,rad) private(i)
    ! $omp sections
    ! $omp section
    do li=1,nl(1)
       !write(91,'(I6,I4)') i!,omp_get_thread_num()
       lon = lon + tm**nint(dat(1,li)) * dat(2,li) * cos( dat(3,li) + dat(4,li)*tm )
    end do
    lon = rev(lon)
    
    ! $omp section
    do li=nl(1),nl(1)+nl(2)
       !write(92,'(I6,I4)') i!,omp_get_thread_num()
       lat = lat + tm**nint(dat(1,li)) * dat(2,li) * cos( dat(3,li) + dat(4,li)*tm )
    end do
    
    ! $omp section
    do li=nl(1)+nl(2),nl(1)+nl(2)+nl(3)
       !write(93,'(I6,I4)') i!,omp_get_thread_num()
       rad = rad + tm**nint(dat(1,li)) * dat(2,li) * cos( dat(3,li) + dat(4,li)*tm )
    end do
    ! $omp end sections nowait
    ! $omp end parallel
    
    !print*,n
  end subroutine vsop_lbr
  !*********************************************************************************************************************************
  
  
  
end module TheSky_vsop
!***********************************************************************************************************************************
