!> \file stars.f90  Star procedures for libTheSky


!  Copyright (c) 2002-2026  Marc van der Sluys - Nikhef/Utrecht University - marc.vandersluys.nl
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
!> \brief Procedures for stars

module TheSky_stars
  implicit none
  save
  
contains
  
  !*********************************************************************************************************************************
  !> \brief  Calculate the apparent position of selected bright stars close to the ecliptic, correcting for proper motion, 
  !!         precession, nutation and aberration
  !!
  !! \param jd  Julian day
  
  subroutine calcstars(jd)  
    use SUFR_kinds, only: double
    use SUFR_constants, only: pi, ras2r, jd2000
    
    use TheSky_nutation, only: nutation
    use TheSky_coordinates, only: precess_eq, aberration_eq, eq_2_ecl
    use TheSky_stardata, only: starnames,starnamesnl, starmags,starrads,starcons,starconsnl,starconsabr
    use TheSky_stardata, only: nstars, starra,stardec,starl,starb
    
    implicit none
    real(double), intent(in) :: jd
    
    integer :: i
    real(double) :: t,dpsi,eps0,deps,eps,  ra(nstars),dec(nstars),pma(nstars),pmd(nstars)
    real(double) :: da1(nstars),dd1(nstars),da2(nstars),dd2(nstars),  l(nstars),b(nstars)
    
    
    ! Data shared in module stardata  -  Add 11 Pleiades members: nstars 6->17:
    starnames = (/'Pleiades  ','Aldebaran ','Pollux    ','Regulus   ','Spica     ','Antares   ','16 Tau    ','Electra   ', &
         '18 Tau    ','Taygete   ','Maia      ','21 Tau    ','22 Tau    ','Merope    ','Alcyone   ','Atlas     ','28 Tau    '/)
    starnamesnl = (/'de Pleiaden','Aldebaran  ','Pollux     ','Regulus    ','Spica      ','Antares    ','16 Tau     ', &
         'Electra    ','18 Tau     ','Taygete    ','Maia       ','21 Tau     ','22 Tau     ','Merope     ','Alcyone    ', &
         'Atlas      ','28 Tau     '/)
    starmags = (/1.6,0.85,1.14,1.35,0.98,0.96,5.46,3.70,5.64,4.30,3.87,5.76,6.43,4.18,2.87,3.63,5.09/)
    starrads = (/110.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)*60*ras2r
    starcons = (/'Taurus ','Taurus ','Gemini ','Leo    ','Virgo  ','Scorpio','Taurus ','Taurus ','Taurus ','Taurus ','Taurus ', &
         'Taurus ','Taurus ','Taurus ','Taurus ','Taurus ','Taurus '/)
    starconsnl = (/'Stier     ','Stier     ','Tweelingen','Leeuw     ','Maagd     ','Schorpioen','Stier     ','Stier     ', &
         'Stier     ','Stier     ','Stier     ','Stier     ','Stier     ','Stier     ','Stier     ','Stier     ','Stier     '/)
    starconsabr = (/'Tau','Tau','Gem','Leo','Vir','Sco','Tau','Tau','Tau','Tau','Tau','Tau','Tau','Tau','Tau','Tau','Tau'/)
    
    
    ! RA & Dec in Rad, FK5 2000.0/2000.0:
    ra  = (/0.99047d0,1.2039309324d0,2.0303233601d0,2.6545229429d0,3.5133171901d0,4.3171054224d0,0.980890d0,0.981202d0,0.982453d0, &
         0.982657d0,0.985355d0,0.985704d0,0.986322d0,0.987536d0,0.992591d0,0.999906d0,1.000015d0/)
    dec = (/0.42092d0,0.2881416664d0,0.4891494426d0,0.2088671634d0,-0.1948018168d0,-0.4613254715d0,0.423931d0,0.420857d0, &
         0.433525d0,0.427034d0,0.425298d0,0.428561d0,0.428095d0,0.417977d0,0.420712d0,0.419810d0,0.421264d0/)
    
    
    ! Proper motion in mas/yr->rad/yr:
    pma = (/0.,62.78,-625.69,-249.40,-42.50,-10.16,0.011,0.019,0.020,0.018,0.020,0.011,0.014,0.021,0.019,0.018,0.013/)/6.48d8*pi
    pmd = (/0.,-189.35,-45.96,4.91,-31.73,-23.21,-0.046,-0.046,-0.047,-0.045,-0.046,-0.042,-0.044,-0.045,-0.046,-0.047,-0.050/) &
         /6.48d8*pi
    
    
    ! PM = change of epoch:
    t = (jd-jd2000)/365.250d0  ! Julian years since 2000.0
    ra  = ra  + pma*t
    dec = dec + pmd*t
    
    
    ! Precession = change of equinox (from J2000.0 to EoD):
    do i=1,nstars
       call precess_eq(jd2000,jd,ra(i),dec(i))
    end do
    
    
    ! Nutation, 1st order, Meeus Eq.23.1:
    t = (jd-jd2000)/365250.d0  ! Julian millennia since 2000.0
    call nutation(t, dpsi,eps0,deps)
    eps = eps0 + deps
    
    do i=1,nstars
       da1(i) = (cos(eps0)+sin(eps0)*sin(ra(i))*tan(dec(i)))*dpsi - (cos(ra(i))*tan(dec(i)))*deps
       dd1(i) = (sin(eps0)*cos(ra(i)))*dpsi + sin(ra(i))*deps
    end do
    
    
    ! Aberration in ra,dec:
    do i=1,nstars
       call aberration_eq(jd, ra(i),dec(i), da2(i), dd2(i), eps0=eps0)
    end do
    
    ra  = ra  + da1 + da2
    dec = dec + dd1 + dd2
    
    
    ! Calculate ecliptic coordinates l,b:
    do i=1,nstars
       call eq_2_ecl(ra(i),dec(i),eps,l(i),b(i))
    end do
    
    
    ! Transport to module stardata:
    starra  = ra
    stardec = dec
    starl   = l
    starb   = b 
    
  end subroutine calcstars
  !*********************************************************************************************************************************
  
  
  
  
  
end module TheSky_stars
!***********************************************************************************************************************************
