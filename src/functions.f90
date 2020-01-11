!> \file functions.f90  Contains general functions for libTheSky


!  Copyright (c) 2002-2020  Marc van der Sluys - marc.vandersluys.nl
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


!  const_id:                            Identify the constellation a,d lies in

!  plsep:                               Calculates angular separation between two planets
!  plpa:                                Calculates position angle between two planets

!  conabr2conid:                        Convert a three-letter constellation abbreviation to a constellation ID number

!  fillfloat                            Print a floating-point number with leading zeroes in format (0...0)Ftot.dgt to a string

!  ds:                                  Returns distance as a nice string in AU
!  ds1:                                 Returns distance as a nice, smaller string in AU



!***********************************************************************************************************************************
!> \brief  Assorted procedures

module TheSky_functions
  implicit none
  save
  
contains
  
  !*********************************************************************************************************************************
  !> \brief  Set a location, using the variables in TheSky_local
  !!
  !! \param lat0    Latitude (radians, >0 is northern hemisphere)
  !! \param lon0    Longitude (radians, >0 is east of Greenwich)
  !! \param height  Altitude/elevation above sealevel, in metres
  !! \param tz0     Default timezone (standard time, no DST, >0 is east of Greenwich)
  !! \param dsttp   DST type: 0-none, 1-EU, 2-USA
  
  subroutine set_TheSky_location(lat0,lon0,height, tz0,dsttp)
    use SUFR_kinds, only: double
    use TheSky_local, only: llon0=>lon0,llat0=>lat0,lheight=>height,ltz=>tz,ltz0=>tz0,ldsttp=>dsttp
    
    implicit none
    real(double), intent(in) :: lat0,lon0,height, tz0
    integer, intent(in) :: dsttp
    
    llat0 = lat0
    llon0 = lon0
    lheight = height
    
    ltz  = tz0     ! Current timezone (standard time or DST)
    ltz0 = tz0     ! Default timezone (standard time, no DST)
    ldsttp = dsttp
    
  end subroutine set_TheSky_location
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Identify the constellation a,d lies in
  !!
  !! \param jd   Equinox in JD
  !! \param ra   RA in radians
  !! \param dec  Dec in radians
  
  function const_id(jd, ra,dec)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd1875, r2d,r2h
    use SUFR_angles, only: rev
    use SUFR_system, only: warn
    
    use TheSky_coordinates, only: precess_eq
    use TheSky_stardata, only: coniddecl, conidrau, conidral, conid
    
    implicit none
    real(double), intent(in) :: jd,ra,dec
    integer :: line,const_id
    real(double) :: lra,ldec
    
    lra = ra
    ldec = dec
    
    ! Precess position to 1875.0 equinox:
    call precess_eq(jd,jd1875,lra,ldec)
    lra  = rev(lra)*r2h
    ldec = ldec*r2d
    
    ! Find constellation such that the declination entered is higher than the lower boundary of the constellation 
    !  when the upper and lower right ascensions for the constellation bound the entered right ascension:
    line = 1
    
    do while(coniddecl(line).gt.ldec)
       line = line + 1
    end do
    
    
    do
       do while(conidrau(line).le.lra)
          line = line + 1
       end do
       
       do while(conidral(line).gt.lra)
          line = line + 1
       end do
       
       ! If constellation has been found, write result and exit/return.  Otherwise, continue the search by returning to rau.
       if(lra.ge.conidral(line).and.lra.lt.conidrau(line).and.coniddecl(line).le.ldec) then
          const_id = conid(line)
          exit
       else if(conidrau(line).gt.lra) then  ! Nothing found
          const_id = 0
          call warn('const_id():  constellation ID was not found', 0)
          exit
       end if  ! else cycle
    end do
    
  end function const_id
  !*********************************************************************************************************************************
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculates the angular separation between two planets
  !!
  !! \param jd0  Julian day for calculation
  !! \param p1   ID of planet 1
  !! \param p2   ID of planet 2
  !!
  !! \note
  !!  Uses asep()
  
  function plsep(jd0, p1,p2)
    use SUFR_kinds, only: double
    use SUFR_angles, only: asep
    
    use TheSky_planets, only: planet_position
    use TheSky_planetdata, only: planpos
    
    implicit none
    real(double), intent(in) :: jd0
    integer, intent(in) :: p1,p2
    real(double) :: plsep,jd,l1,l2,b1,b2
    
    jd = jd0
    call planet_position(jd,p1)
    l1  = planpos(25)
    b1 = planpos(26)
    
    call planet_position(jd,p2)
    l2  = planpos(25)
    b2 = planpos(26)
    
    plsep = asep(l1,l2, b1,b2)
    
  end function plsep
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculates the position angle of planet 2 with respect to planet 1, COUNTERCLOCKWISE from the north
  !!
  !! \param jd0  Julian day for calculation
  !! \param p1   ID of planet 1
  !! \param p2   ID of planet 2
  !!
  !! \note
  !!   Uses calpa()
  
  function plpa(jd0,p1,p2)
    use SUFR_kinds, only: double
    use SUFR_angles, only: calpa
    
    use TheSky_planets, only: planet_position
    use TheSky_planetdata, only: planpos
    
    implicit none
    real(double), intent(in) :: jd0
    integer, intent(in) :: p1,p2
    real(double) :: plpa,jd,l1,l2,b1,b2
    
    jd = jd0
    call planet_position(jd,p1)
    l1  = planpos(5)     ! RA and Dec
    b1 = planpos(6)
    
    call planet_position(jd,p2)
    l2  = planpos(5)
    b2 = planpos(6)
    
    plpa = calpa(l1,l2, b1,b2)
    
  end function plpa
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief Convert a three-letter constellation abbreviation to a constellation ID number
  !!
  !! \param  myconabr  Constellation abbreviation (e.g. And)
  !! \param myconid   Constellation ID number (output)
  
  subroutine conabr2conid(myconabr, myconid)
    use TheSky_stardata, only: conabr, nconstel
    
    implicit none
    character, intent(in) :: myconabr*(*)
    integer, intent(out) :: myconid 
    integer :: con
    
    myconid = 0
    do con=1,nconstel
       if(trim(conabr(con)).eq.trim(myconabr)) then
          myconid = con
          exit
       end if
    end do
    
  end subroutine conabr2conid
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Print a floating-point number with leading zeroes in format (0...0)Ftot.dgt to a string
  !!
  !! \param x    Number to be printed
  !! \param tot  Total length of printed number - Ftot.dgt
  !! \param dgt  Number of decimal digits - Ftot.dgt
  
  function fillfloat(x,tot,dgt)
    use SUFR_kinds, only: double
    implicit none
    real(double), intent(in) :: x
    integer, intent(in) :: tot,dgt
    character :: fillfloat*(99), fmt*(99)
    integer :: dec
    
    dec = tot-dgt-1    ! Number of decimal places before decimal sign
    if(x.lt.0.d0) then
       dec = dec-1     ! Need extra room for minus sign
       write(fmt,'(A,2(I3.3,A,I3.3,A))') '(A1,I',dec,'.',dec,',F',dgt+1,'.',dgt,')'
       write(fillfloat,trim(fmt)) '-', floor(abs(x)), abs(x)-floor(abs(x))
    else
       write(fmt,'(A,2(I3.3,A,I3.3,A))') '(I',dec,'.',dec,',F',dgt+1,'.',dgt,')'
       write(fillfloat,trim(fmt)) floor(x), abs(x)-floor(abs(x))
    end if
    
  end function fillfloat
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Print planet distance in AU as a nice string, but use km for the Moon
  !!
  !! \param d  Distance (AU)
  
  function ds(d)
    use SUFR_kinds, only: double
    use SUFR_constants, only: au
    
    implicit none
    real(double), intent(in) :: d
    character :: ds*(11)
    
    write(ds,'(F11.8)') d
    if(d.lt.0.01d0) write(ds,'(F11.4)') d*au/1.d5  ! For the Moon
    
  end function ds
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Print planet distance in AU as a nice string, but use km for the Moon - smaller string
  !!
  !! \param d  Distance (AU)
  
  function ds1(d)
    use SUFR_kinds, only: double
    use SUFR_constants, only: au
    
    implicit none
    real(double), intent(in) :: d
    character :: ds1*(9)
    
    write(ds1,'(F9.6)') d
    if(d.lt.0.01d0) write(ds1,'(F9.2)') d*au/1.d5  ! For the Moon
    
  end function ds1
  !*********************************************************************************************************************************
  
end module TheSky_functions
!***********************************************************************************************************************************
