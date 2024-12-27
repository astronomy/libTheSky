!> \file date_time.f90  Contains date and time procedures for libTheSky


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


!  set_date_and_time:                   Set global date/time variables (year, month, ..., minute, second) to specified values
!  get_date_and_time:                   Retrieve the current global date/time variables (year, ..., second, stored in TheSky_local)
!  set_date_and_time_to_system_clock:   Set global date/time variables (year, month, ..., minute, second) to system clock
!  set_date_and_time_to_jd2000:         Set global date/time variables (year, month, ..., minute, second) to JD2000.0

!  calctime:                            Calculates ut, jd and jde
!  calc_gmst:                           Calculate Greenwich Mean Sidereal Time in RAD!
!  calc_deltat:                         Calculates deltat from jd: SLOW!
!  calc_deltat_ymd:                     Calculates deltat from y,m,d, faster
!  find_deltat_in_range:                Find a precise value for DeltaT through linear interpolation in tabulated values

!  jd2dtm:                              Converts JD to y,m,d,h,m,s, input in UT, output in LT
!  jd2dthm:                             Converts JD to y,m,d,h,m, input in UT, output in LT (no seconds)
!  jd2ltime:                            Converts JD to time;  input JD in UT, output time in hours LT

!  printdate:                           Converts JD to y,m,d,h,m,s and prints it out
!  printdate1:                          prints date/time of JD (UT) without hard return
!  print_date_time_and_location         Display a 'program header' with date, time and location

!  dls:                                 Beginning and end of daylight savings for given year
!  gettz:                               Get the time zone (1 or 2) for JD
!  dow:                                 Calculates day of week
!  woy:                                 Calculates the week-of-year number



!***********************************************************************************************************************************
!> \brief Date and time procedures

module TheSky_datetime
  implicit none
  save
  
contains
  
  !*********************************************************************************************************************************
  !> \brief Set global date/time variables (year, month, ..., minute, second in TheSky_local) to specified values
  !!
  !! \param year    Year
  !! \param month   Month
  !! \param day     Day
  !!
  !! \param hour    Hour
  !! \param minute  Minute
  !! \param second  Second
  
  subroutine set_date_and_time(year,month,day, hour,minute, second)
    use SUFR_kinds, only: double
    use TheSky_local, only: lyear=>year,lmonth=>month,lday=>day, lhour=>hour,lminute=>minute,lsecond=>second
    
    implicit none
    integer, intent(in) :: year,month,day, hour,minute
    real(double), intent(in) :: second
    
    lyear   = year
    lmonth  = month
    lday    = dble(day)
    lhour   = hour
    lminute = minute
    lsecond = second
    
  end subroutine set_date_and_time
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Retrieve the current global date/time variables (year, month, ..., minute, second, stored in TheSky_local)
  !!
  !! \param year    Year (output)
  !! \param month   Month (output)
  !! \param day     Day (output)
  !!
  !! \param hour    Hour (output)
  !! \param minute  Minute (output)
  !! \param second  Second (output)
  
  subroutine get_date_and_time(year,month,day, hour,minute, second)
    use SUFR_kinds, only: double
    use TheSky_local, only: lyear=>year,lmonth=>month,lday=>day, lhour=>hour,lminute=>minute,lsecond=>second
    
    implicit none
    integer, intent(out) :: year,month,day, hour,minute
    real(double), intent(out) :: second
    
    year   = lyear
    month  = lmonth
    day    = nint(lday)
    
    hour   = lhour
    minute = lminute
    second = lsecond
    
  end subroutine get_date_and_time
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief Set global date/time variables (year, month, ..., minute, second in TheSky_local) to system clock
  !!
  !! \param ut  Return UT instead of local system time (optional).
  
  subroutine set_date_and_time_to_system_clock(ut)  
    use SUFR_date_and_time, only: system_clock_2_ymdhms, correct_time
    use TheSky_local, only: year,month,day, hour,minute,second, tz
    
    implicit none
    logical, optional, intent(in) :: ut
    integer :: dy
    
    call system_clock_2_ymdhms(year,month,dy, hour,minute,second, tz=tz)
    
    if(present(ut)) then
       if(ut) then
          hour = hour - tz
          tz = 0.d0
          call correct_time(year,month,dy, hour,minute,second)
       end if
    end if
    
    day = dble(dy)
    
  end subroutine set_date_and_time_to_system_clock
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Set global date/time variables (year, month, ..., minute, second in TheSky_local) to JD2000.0
  
  subroutine set_date_and_time_to_jd2000()
    implicit none
    
    call set_date_and_time(2000,1,1, 0,0,0.d0)
    
  end subroutine set_date_and_time_to_jd2000
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute UT, JD, JDE, DeltaT and TZ, using the date and (local) time and TZ stored in the module TheSky_local
  !!
  !! \param ut   Universal time (output)
  !! \param jd   Julian date (output)
  !! \param jde  Julian date (output)
  !!
  !! \note
  !! - Date and time are obtained from year, month, day, hour, minute, second, through the module TheSky_local, assuming LOCAL time
  !! - DeltaT and TZ are returned through the module TheSky_local
  !! - Computed JD is in UNIVERSAL time
  !!
  !! \todo  Remove recomputation of 'UT JD'!!!
  
  subroutine calctime(ut,jd,jde)
    use SUFR_kinds, only: double
    use SUFR_date_and_time, only: cal2jd
    use SUFR_numerics, only: dne
    use TheSky_local, only: year,month,day, hour,minute,second, tz,deltat
    
    implicit none
    real(double), intent(out) :: ut,jd,jde
    real(double) :: lt, oldtz
    
    ! Compute local time, UT and JD:
    lt = dble(hour) + (dble(minute) + second/60.d0)/60.d0
    ut = lt - tz
    jd = cal2jd(year,month,day+ut/24.d0)              ! This is the 'UT JD'
    
    ! Determine the correct timezone from the JD, and recompute UT and JD:
    oldtz = tz
    tz = gettz(jd)
    if(dne(tz, oldtz)) then                           ! CHECK - this may be a very stupid idea - provide UT hms and tz=0, and this
       ut = lt - tz                                   ! will currently return a non-UT ut, jd if dsttp = 1,2 or tz0 != 0
       jd = cal2jd(year,month,day+ut/24.d0)           ! This is the 'UT JD'
    end if
    
    deltat = calc_deltat_ymd(year,month,day)
    jde = jd + deltat/86400.d0
    
  end subroutine calctime
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate Greenwich Mean Sidereal Time for any instant, in radians
  !!
  !! \param  jd         Julian day of computation
  !! \param  DeltaT     ΔT (s) (optional - default: use DelTaT from TheSky_local)
  !! 
  !! \retval calc_gmst  Greenwich Mean Sidereal Time (radians)
  !! 
  !! \see Explanatory Supplement to the Astronomical Almanac, 3rd edition, Eq. 6.66 (2012)
  
  function calc_gmst(jd, DeltaT)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000
    use SUFR_angles, only: rev
    use TheSky_local, only: DeltaTmod=>DeltaT
    
    implicit none
    real(double), intent(in) :: jd
    real(double), intent(in), optional :: DeltaT
    real(double) :: calc_gmst, djd,djd2,djd4, gmst, DeltaTl
    
    DeltaTl = DeltaTmod                   ! Use Delta T from module by default
    if(present(DeltaT)) DeltaTl = DeltaT  ! Use Delta T from optional variable if present
    
    djd  = jd-jd2000                      ! Julian Days after 2000.0 UT
    djd2 = djd**2
    djd4 = djd2**2
    
    gmst = 4.89496121042905d0 + 6.30038809894828323d0*djd + 5.05711849d-15*djd2 - 4.378d-28*djd2*djd - 8.1601415d-29*djd4 &
         - 2.7445d-36*djd4*djd + 7.0855723730d-12*DeltaTl  ! Eq. 6.66
    
    calc_gmst = rev(gmst)          ! If corrected for equation of the equinoxes: agst = rev(gmst + dpsi*cos(eps))
    
  end function calc_gmst
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculate Greenwich Mean Sidereal Time for any instant, in radians, using Meeus
  !!
  !! \param  jd          Julian day of computation
  !! \retval gmst_meeus  Greenwich Mean Sidereal Time in radians
  !!
  !! \see Meeus, Sect. 12: Siderial time at Greenwich; Eq.12.4
  
  function gmst_meeus(jd)
    use SUFR_kinds, only: double
    use SUFR_constants, only: jd2000
    use SUFR_angles, only: rev
    
    implicit none
    real(double), intent(in) :: jd
    real(double) :: gmst_meeus, djd,tjc,tjc2, gmst
    
    djd  = jd-jd2000                ! Julian Days after 2000.0 UT
    tjc  = djd/36525.d0             ! Julian Centuries after 2000.0 UT
    tjc2 = tjc**2
    
    gmst = 4.89496121273579229d0 + 6.3003880989849575d0*dJD + 6.77070812713916d-6*tjc2 - 4.5087296615715d-10*tjc2*tjc  ! Eq. 12.4
    
    gmst_meeus = rev(gmst)          ! If corrected for equation of the equinoxes: agst = rev(gmst + dpsi*cos(eps))
    
  end function gmst_meeus
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Computes JD, DeltaT and TZ, from date/time variables in module TheSky_local
  !!
  !! \param JD   Julian date (output)
  !!
  !! \note
  !! - Date and time are obtained from year, month, day, hour, minute, second, through the module TheSky_local
  !! - DeltaT and TZ are returned through the module TheSky_local
  
  subroutine localtime2jd(jd)
    use SUFR_kinds, only: double
    use SUFR_date_and_time, only: cal2jd
    use TheSky_local, only: year,month,day, hour,minute,second, tz,deltat
    
    implicit none
    real(double), intent(out) :: jd
    real(double) :: lt,ut
    
    lt = hour + (minute + second/60.d0)/60.d0         ! LT
    ut = lt - tz                                      ! UT
    jd = cal2jd(year, month, day+ut/24.d0)            ! UT
    
    tz = gettz(jd)                                    ! Compute actual timezone
    ut = lt - tz                                      ! UT
    jd = cal2jd(year, month, day+ut/24.d0)            ! UT
    
    deltat = calc_deltat_ymd(year,month,day)
    ! jde = jd + deltat/86400.d0
    
  end subroutine localtime2jd
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute DeltaT for a given JD
  !!
  !! \param jd               Julian day
  !! \param force_recompute  Force recomputation of DeltaT, even if the year is the same as in the last call (done in calc_deltat_ymd())
  !! 
  !! \note VERY SLOW, use calc_deltat_ymd() if y,m,d are known
  
  function calc_deltat(jd, force_recompute)
    use SUFR_kinds, only: double
    use SUFR_date_and_time, only: jd2cal
    use SUFR_numerics, only: deq0
    use TheSky_constants, only: deltat_forced
    use TheSky_local, only: year,month,day
    
    implicit none
    real(double), intent(in) :: jd
    logical, intent(in), optional :: force_recompute
    
    real(double) :: calc_deltat,d
    integer :: y,m
    logical :: force_recompute_l
    
    
    ! If deltat_forced > -9999, use that value:
    if(deltat_forced .gt. -9999.d0) then
       calc_deltat = deltat_forced
       return
    end if
    
    
    ! Optional variable:
    force_recompute_l = .false.
    if(present(force_recompute)) force_recompute_l = force_recompute
    
    if(deq0(jd)) then  ! Use variables from module local
       y = year
       m = month
       d = day
    else
       call jd2cal(jd, y,m,d)  ! SLOW!
    end if
    
    calc_deltat = calc_deltat_ymd(y,m,d, force_recompute_l)
    
  end function calc_deltat
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute DeltaT for given y,m,d
  !!
  !! \param y                Year
  !! \param m                Month
  !! \param d                Day
  !! 
  !! \param force_recompute  Force recomputation of DeltaT, even if the year is the same as in the last call
  !! 
  !! \note
  !! - Faster than calc_deltat. Use this routine rather than calc_deltat() if y,m,d are known
  
  function calc_deltat_ymd(y,m,d, force_recompute)
    use SUFR_kinds, only: double
    use SUFR_numerics, only: deq
    use TheSky_constants, only: deltat_0, deltat_accel, deltat_change, deltat_minyr, deltat_maxyr, deltat_years, deltat_values, deltat_forced
    
    implicit none
    integer, intent(in) :: y,m
    real(double), intent(in) :: d
    logical, intent(in), optional :: force_recompute
    
    integer :: yr1,yr2
    real(double) :: calc_deltat_ymd, y0,dy, ddt
    real(double), save :: calc_deltat_old, y0_old
    logical :: force_recompute_l
    
    ! If deltat_forced > -9999, use that value:
    if(deltat_forced .gt. -9999.d0) then
       calc_deltat_ymd = deltat_forced
       return
    end if
    
    
    ! Optional variable:
    force_recompute_l = .false.
    if(present(force_recompute)) force_recompute_l = force_recompute
    
    calc_deltat_ymd = 60.d0
    
    y0 = (dble(m-1)+((d-1)/31.d0))/12.d0 + y  ! ~decimal year
    
    ! Return old value if same instance:
    if(deq(y0,y0_old) .and. .not. force_recompute_l) then
       calc_deltat_ymd = calc_deltat_old
       return
    end if
    
    
    yr1 = deltat_minyr
    yr2 = deltat_maxyr
    
    if(y.ge.yr1.and.y.le.yr2) then  ! Historical value is known; look it up and interpolate
       
       calc_deltat_ymd = find_deltat_in_range(y,y0)
       
    else                            ! Historical value is not known, extrapolate
       
       if(y.lt.yr1) then  ! Earlier than historical record
          dy = y0 - dble(deltat_minyr)
          ddt = (deltat_values(2)-deltat_values(1)) / (deltat_years(2)-deltat_years(1))  ! 1 and 2 are spaced by 100 yr
          calc_deltat_ymd = deltat_values(1)  +  ddt * dy  +  deltat_accel * dy**2
       end if
       
       if(y.gt.yr2) then  ! Future extrapolation
          dy = y0 - dble(deltat_maxyr)
          calc_deltat_ymd = deltat_0 + deltat_change * dy  +  deltat_accel * dy**2
       end if
       
    end if  ! if(y.gt.yr0.and.y.lt.yr1)
    
    
    calc_deltat_old = calc_deltat_ymd
    y0_old = y0
    
    ! write(0,'(A)')'  WARNING:  fixed DeltaT!!!'
    ! calc_deltat_ymd = 66.4d0
    
  end function calc_deltat_ymd
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Find a precise value for DeltaT through linear interpolation of the two adjacent tabulated values
  !!
  !! \param y   Current year
  !! \param y0  Current date, as fractional year
  
  function find_deltat_in_range(y,y0)
    use SUFR_kinds, only: double
    use TheSky_constants, only: deltat_n, deltat_years, deltat_values
    
    implicit none
    real(double) :: find_deltat_in_range
    integer, intent(in) :: y
    real(double),intent(in) :: y0
    integer :: i
    real(double) :: dt0,dt,yr0,yr,a
    
    do i=deltat_n-1,1,-1              ! i is usually near deltat_n, so start from the back
       if(deltat_years(i).le.y) exit
    end do
    i = max(i,1)  ! i can be 0 if whole do loop is run without match
    
    yr0 = deltat_years(i)
    dt0 = deltat_values(i)
    yr  = deltat_years(i+1)
    dt  = deltat_values(i+1)
    
    a = (dt-dt0)/dble(yr-yr0)
    find_deltat_in_range = dt0 + a*(y0-yr0)
    
  end function find_deltat_in_range
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Compute DeltaT for given JD, using a simple parabolic approximation
  !!
  !! \param jd       Julian day for computation
  !!
  !! \retval deltat  The value for DeltaT (= TDT-UT) in seconds.
  
  function calc_deltat_approx(jd)
    use SUFR_kinds, only: double
    use TheSky_constants, only: jd1820
    
    implicit none
    real(double), intent(in) :: jd
    real(double) :: calc_deltat_approx
    
    ! 12 + 0.5 * 1.8e-3/86400/(36525*86400) * ((jd-jd1820)*86400)**2  ! Comprehensible notation
    calc_deltat_approx = 12.d0 + 0.5d0 * 1.8d-3 / 36525.d0 * (jd-jd1820)**2  ! Simplified notation
    
  end function calc_deltat_approx
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert a Julian day (UT) to LOCAL date and time (h,m,s)
  !!
  !! \param  jd  Julian day (UT)
  !!
  !! \param yy  Year (CE, LT) (output)
  !! \param mm  Month (LT) (output)
  !! \param d   Day (LT) (output)
  !! \param h   Hour (LT) (output)
  !! \param m   Minute (LT) (output)
  !! \param s   Second (+ fraction, LT) (output)
  
  subroutine jd2dtm(jd,  yy,mm,d, h,m,s)
    use SUFR_kinds, only: double, dbl
    use SUFR_system, only: quit_program_error
    use SUFR_constants, only: mlen
    use SUFR_date_and_time, only: jd2cal, leapyr
    use TheSky_local, only: tz
    
    implicit none
    real(double), intent(in) :: jd
    integer, intent(out) :: yy,mm,d,h,m
    real(double), intent(out) :: s
    real(double) :: dd,tm
    
    call jd2cal(jd + tz/24.d0,  yy,mm,dd)  ! in LT
    if(mm.lt.1 .or. mm.gt.12) call quit_program_error('libTheSky: jd2dtm(): something went wrong in the jd2cal() call: month must be between 1 and 12.', 1)  ! Needed because of mlen()
    
    ! jd2cal returns zeroes if JD not defined (i.e., JD=-huge), and mlen(mm) is not defined - catch this:
    if(yy.eq.0.and.mm.eq.0) then
       d = 0
       h = 0
       m = 0
       s = 0.0_dbl
       return
    end if
    
    
    mlen(2) = 28 + leapyr(yy)
    
    d  = int(dd)
    tm = (dd - dble(d))*24.d0
    h  = int(tm)
    m  = int((tm-h)*60.d0)
    s  = (tm-h-m/60.d0)*3600.d0
    
    if(s.gt.59.999) then
       s = 0.d0
       m = m + 1
    end if
    if(m.eq.60) then
       m = 0
       h = h+1
    end if
    if(h.eq.24) then
       h = 0
       d = d+1
    end if
    if(d.gt.mlen(mm)) then
       d  = d - mlen(mm)
       mm = mm + 1
    end if
    if(mm.gt.12) then
       mm = mm - 12
       yy = yy + 1
    end if
    
  end subroutine jd2dtm
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert a Julian day (UT) to LOCAL date and time (h,m - no seconds)
  !!
  !! \param  jd  Julian day (UT)
  !! \param yy   Year (CE, LT) (output)
  !! \param mm   Month (LT) (output)
  !! \param d    Day (LT) (output)
  !! \param h    Hour (LT) (output)
  !! \param m    Minute (LT) (output)
  
  subroutine jd2dthm(jd, yy,mm, d,h,m)
    use SUFR_kinds, only: double
    use SUFR_system, only: quit_program_error
    use SUFR_constants, only: mlen
    use SUFR_date_and_time, only: jd2cal, leapyr
    use TheSky_local, only: tz
    
    implicit none
    real(double), intent(in) :: jd
    integer, intent(out) :: yy,mm,d,h,m
    real(double) :: dd,tm
    
    call jd2cal(jd+tz/24.d0, yy,mm,dd)  ! in LT
    if(mm.lt.1 .or. mm.gt.12) call quit_program_error('libTheSky: jd2dthm(): something went wrong in the jd2cal() call: month must be between 1 and 12.', 1)  ! Needed because of mlen()
    
    ! jd2cal returns zeroes if JD not defined (i.e., JD=-huge), and mlen(mm) is not defined - catch this:
    if(yy.eq.0 .and. m.eq.0) then
       d = 0
       h = 0
       m = 0
       return
    end if
    
    mlen(2) = 28 + leapyr(yy)
    
    d  = int(dd)
    tm = (dd - dble(d))*24
    h  = int(tm)
    m  = nint((tm-h)*60)
    
    if(m.ge.60) then
       m = m-60
       h = h+1
    end if
    if(h.ge.24) then
       h = h-24
       d = d+1
    end if
    if(d.gt.mlen(mm)) then
       d  = d - mlen(mm)
       mm = mm + 1
    end if
    if(mm.gt.12) then
       mm = mm - 12
       yy = yy + 1
    end if
    
  end subroutine jd2dthm
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert a Julian day (UT) to a local time (LT, h)
  !!
  !! \param jd0  Julian day (UT)
  
  function jd2ltime(jd0)
    use SUFR_kinds, only: double
    use SUFR_date_and_time, only: jd2cal
    use TheSky_local, only: tz
    
    implicit none
    real(double), intent(in) :: jd0
    real(double) :: jd1,dd,jd2ltime
    integer :: d,mm,yy
    
    jd1 = jd0 + tz/24.d0       ! UT -> LT
    call jd2cal(jd1,yy,mm,dd)
    d  = int(dd)
    jd2ltime = (dd - dble(d))*24.d0
    
  end function jd2ltime
  !*********************************************************************************************************************************
  
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Prints date/time of a given Julian day (UT) to standard output
  !!
  !! \param jd   Julian day (UT)
  !! \param jde  Julian day ephemeris time (dynamical time)
  !!
  !! \note
  !! - calls printdate1()
  
  subroutine printdate(jd, jde)
    use SUFR_kinds, only: double
    implicit none
    real(double), intent(in) :: jd
    real(double), intent(in), optional :: jde

    if(present(jde)) then
       call printdate1(jd, jde)
    else
       call printdate1(jd)
    end if
    write(*,*)''
    
  end subroutine printdate
  !*********************************************************************************************************************************
  
  !*********************************************************************************************************************************
  !> \brief  Prints date/time of a given Julian day (UT) to standard output, but without a newline
  !!
  !! \param jd   Julian day (UT)
  !! \param jde  Julian day ephemeris time (dynamical time)
  
  subroutine printdate1(jd, jde)
    use SUFR_kinds, only: double
    use SUFR_date_and_time, only: jd2cal
    
    implicit none
    real(double), intent(in) :: jd
    real(double), intent(in), optional :: jde
    real(double) :: dd,tm,s
    integer :: d,mm,yy,h,m
    
    call jd2cal(jd+1.d-10, yy,mm,dd)
    
    d  = int(dd)
    
    tm = (dd - dble(d))*24.d0
    h  = int(tm)
    m  = int((tm-h)*60.d0)
    s  = (tm-h-m/60.d0)*3600.d0
    
    if(s.gt.59.999d0) then
       s = s - 60.d0
       m = m + 1
    end if
    if(m.ge.60) then
       m = m - 60
       h = h + 1
    end if
    if(h.ge.24) then
       h = h - 24
       d = d + 1
    end if
    
    write(*,'(A,F0.6)', advance='no') '  JD:  ', jd
    if(present(jde)) write(*,'(A,F0.6)', advance='no') '  JDE:  ', jde
    write(*,'(A,I0,2I3.2, A,2I3.2,F7.3,A)', advance='no') '   date:  ',yy,mm,d, '   time:  ',h,m,s,'  UT'
    
  end subroutine printdate1
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Find the two Julian days of the beginning and the end of daylight-savings time in the EU for a given year
  !!
  !! \param  yr   Year (CE)
  !!
  !! \param jdb  Julian day of beginning of DST (output)
  !! \param jde  Julian day of end of DST (output)
  
  subroutine dls(yr, jdb,jde)
    use SUFR_kinds, only: double
    use SUFR_date_and_time, only: cal2jd
    
    implicit none
    integer, intent(in) :: yr
    real(double), intent(out) :: jdb,jde
    real(double) :: d1,d2
    integer :: lyr,m1,m2
    
    lyr = yr
    
    d1 = 31.d0
    m1 = 3
    d2 = 31.d0
    m2 = 10
    
    d1  = d1 - dble(dow(cal2jd(lyr,m1,d1)))
    jdb = cal2jd(lyr,m1,d1)
    d2  = d2 - dble(dow(cal2jd(lyr,m2,d2)))
    jde = cal2jd(lyr,m2,d2)
    
  end subroutine dls
  !*********************************************************************************************************************************
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Returns time zone: tz0 or tz0+1
  !!
  !! \param jd      Julian day (UT?)
  !! \param ltz0    Default time zone for the current location ('winter time')
  !! \param ldsttp  Daylight-savings time rules to use:  0: no DST (e.g. UT),  1: EU,  2: USA/Canada
  !! 
  !! \note
  !!  - currently implemented for no DST (dsttp=0), EU (dsttp=1) and USA/Canada >2007 (dsttp=2) only
  !!  - tz0 and dsttp can be provided through the module TheSky_local, or using the optional arguments.  Note
  !!    that using the latter will update the former!
  
  function gettz(jd, ltz0,ldsttp)
    use SUFR_kinds, only: double
    use SUFR_system, only: warn
    use SUFR_date_and_time, only: jd2cal, cal2jd
    use TheSky_local, only: dsttp, tz,tz0
    
    implicit none
    real(double), intent(in) :: jd
    real(double), intent(in), optional :: ltz0
    integer, intent(in), optional :: ldsttp
    
    real(double) :: gettz,dd,jd0,d0
    integer :: m,m1,m2,y
    
    
    ! Handle optional variables:
    if(present(ltz0))   tz0 = ltz0
    if(present(ldsttp)) dsttp = ldsttp
    
    gettz = 0.d0
    call jd2cal(jd,y,m,dd)
    
    if(dsttp.lt.0.or.dsttp.gt.2) then
       dsttp = 1
       call warn('gettz(): (l)dsttp must be 0-2; resetting to 1 (Europe)')
    end if
    
    
    if(dsttp.eq.1) then  ! Europe (Netherlands)
       if(y.lt.1977) then
          gettz = 0.d0
       else
          m1 = 3
          m2 = 10
          if(y.lt.1996) m2 = 9
          
          if(m.gt.m1.and.m.le.m2) gettz = 1.d0
          
          if(m.eq.m1.or.m.eq.m2) then
             jd0 = cal2jd(y,m,31.d0) + (m2 - 10)  
             d0  = 31.d0 - dble(dow(jd0)) + (m2 - 10)  
             jd0 = cal2jd(y,m,d0+1.d0/24.d0)  ! 0 UT = 1 MET = 2 MEZT
             if(m.eq.m1.and.jd.gt.jd0) gettz = 1.d0
             if(m.eq.m2.and.jd.gt.jd0) gettz = 0.d0
          end if  ! if(m.eq.m1.or.m.eq.m2)
       end if ! if(y.lt.1977)
    end if  ! if(dsttp.eq.1)
    
    
    if(dsttp.eq.2) then  ! USA/Canada
       
       !if(y.lt.2007) then
       !   m1 = 4        ! April
       !   maxdom1 = 7   ! Maximum day of month: first Sunday
       !   m2 = 10       ! October
       !   maxdom2 = 31  ! Maximum day of month: last Sunday
       !   m = m-1
       !end if
       
       call warn('DST rules not implemented prior to 2007!')
       
       if(y.ge.2007) then
          m1 = 3        ! March
          m2 = 11       ! November
          
          if(m.gt.m1.and.m.lt.m2) gettz = 1.d0
          
          if(m.eq.m1) then
             jd0 = cal2jd(y,m,1.999999d0)          ! 1st of the month, end of the day UT
             d0 = dble(14 - dow(jd0-1))            ! jd0-1: switch from 7 to 1 iso 6 to 0; last possible day of month: 14
             jd0 = cal2jd(y,m,d0+(2.d0-tz)/24.d0)  ! 2h LT
             if(jd.gt.jd0) gettz = 1.d0
          end if
          
          if(m.eq.m2) then
             jd0 = cal2jd(y,m,1.999999d0)          ! 1st of the month, end of the day UT
             d0 = dble(7 - dow(jd0-1))             ! jd0-1: switch from 7 to 1 iso 6 to 0; last possible day of month: 7
             jd0 = cal2jd(y,m,d0+(2.d0-tz)/24.d0)  ! 2h LT
             if(jd.lt.jd0) gettz = 1.d0
          end if
          
          ! call printdate(jd0)
       end if
    end if
    
    gettz = tz0 + gettz
    ! gettz = tz0 + 1.d0  ! Force DST
    
  end function gettz
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Display a banner with date, time and location of calculation.  Computes and returns UT, JD and JDE.
  !! 
  !! \param  op     Output unit (optional, default 6)
  !! \param  nlbef  Number of newlines before output (optional, default 0)
  !! \param  nlaf   Number of newlines after output (optional, default 0)
  !! 
  !! \param  ut     UT: Universal Time (optional)
  !! \param  jd     JD: Julian day (optional)
  !! \param  jde    JDE: Apparent Julian day (optional)
  !! 
  !! \param  locname  Location name (optional)
  !! \param  tzname   Time-zone name (optional)
  !! 
  !! \note
  !! - uses the module local to obtain the variables year, month, etc.
  !! - calls calctime(), gettz()
  
  subroutine print_date_time_and_location(op,nlbef,nlaf, ut,jd,jde, locname,tzname)
    use SUFR_kinds, only: double
    use SUFR_constants, only: r2d,r2h, endays,enmonths
    use SUFR_time2string, only: hms
    use SUFR_text, only: d2s
    
    use TheSky_local, only: year,month,day, hour,minute,second, tz, lat0,lon0,height,deltat
    
    
    implicit none
    integer, intent(in), optional :: op, nlbef,nlaf
    real(double), intent(out), optional :: ut,jd,jde
    character, intent(in), optional :: locname*(*),tzname*(*)
    integer :: il,  opl, nlbefl,nlafl
    real(double) :: gmst, lt,  utl,jdl,jdel, second_fr
    character :: locnamel*(999),tznamel*(99)
    
    ! Optional input parameters:
    opl = 6                            ! Output unit - local variable
    if(present(op)) opl = op
    nlbefl = 0                         ! Number of newlines before output - local variable
    if(present(nlbef)) nlbefl = nlbef
    nlafl = 0                          ! Number of newlines after output - local variable
    if(present(nlaf)) nlafl = nlaf
    locnamel =  ' '                    ! Name of the location - local variable
    if(present(locname)) locnamel = trim(locname)
    tznamel =  ' '                     ! Name of the time zone - local variable
    if(present(tzname)) tznamel = trim(tzname)//','
    
    
    call calctime(utl,jdl,jdel)
    tz = gettz(jdl)
    
    call calctime(utl,jdl,jdel)
    lt = hour + (minute + second/60.d0)/60.d0 + 1.d-50  ! so that 0 doesn't give --:--
    utl = utl + 1.d-50                                  ! so that 0 doesn't give --:--
    second_fr = second - floor(second)                  ! fraction of the second - works for second<0
    
    gmst = calc_gmst(jdl, DeltaT)*r2h  ! Greenwich mean sidereal time in hours
    
    do il=1,nlbefl
       write(opl,'(A)') ''  ! Newline
    end do
    
    write(opl,'(A20,3A,I3,A1,I5,  A9,A9,F4.3,2x,A,  A13,2(A4,A), A4,A,A1)') &
         'LOCAL:      Date: ',trim(endays(dow(jdl))),' ', trim(enmonths(month)),nint(day),',',year,  &
         'Time:',hms(lt),second_fr,'tz: '//trim(tznamel)//' '//d2s(tz,1),   &
         'Location:','l: ',d2s(lon0*r2d,4),'b: ',d2s(lat0*r2d,4), 'h: ',d2s(height,1),'m'
    
    write(opl,'(A17,A9,F4.3, A5,F15.6, A10,F0.2,A1, A6,F15.6, A7,F8.4, 5x,A)') &
         'UNIVERSAL:  UT:',hms(utl),second_fr, 'JD:',jdl, 'DeltaT: ',DeltaT,'s', &
         'JDE:',jdel, 'GMST:',gmst, trim(locnamel)
    
    do il=1,nlafl
       write(opl,'(A)') ''  ! Newline
    end do
    
    ! Optional output parameters:
    if(present(ut))  ut  = utl
    if(present(jd))  jd  = jdl
    if(present(jde)) jde = jdel
    
  end subroutine print_date_time_and_location
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Calculates day of week (0 - Sunday, ... 6)
  !!
  !! \param jd0  Julian day (UT)
  !!
  !! \note Input in UT, output in local time
  !! \todo Switch using dow(jd) to dow_ut(jd+tz/24.d0) to make it general
  
  function dow(jd0)
    use SUFR_kinds, only: double
    use TheSky_local, only: tz
    
    implicit none
    real(double), intent(in) :: jd0
    integer :: dow
    real(double) :: jd,jw
    
    jd  = dble(nint(jd0+tz/24.d0))-0.5d0  ! ~ JD
    jw  = (jd + 1.5d0)/7.d0               ! ~ 'Julian week'
    dow = nint(jd + 1.5d0 - floor(jw)*7)  ! Day of week (0-6)
    
  end function dow
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief   Calculate the week-of-year number
  !!
  !! \param jd  Julian day
  !!
  !! \todo 
  !! - CHECK: int or floor?
  !! - Depends on dow(), which depends on module local -> switch to dow_ut()
  
  function woy(jd)
    use SUFR_kinds, only: double
    use SUFR_date_and_time, only: doy
    
    implicit none
    real(double), intent(in) :: jd
    integer :: woy
    
    woy = int(dble(doy(jd) + 7 - dow(jd))/7.d0)  ! Use int or floor ?
    
  end function woy
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief   Calculate the date of Easter using Gauss' method.
  !!
  !! \param year   Julian/Gregorian year to compute date of Easter for.
  !!
  !! \param month  Julian/Gregorian month of year of Easter Sunday.
  !! \param day    Julian/Gregorian day of month of Easter Sunday.
  !! 
  !! \note
  !!  - Based on Gauß, C.F., 1816. Berichtigung zu dem Aufsatze: Berechnung des Osterfestes. Zeitschrift für
  !!             Astronomie und verwandte Wissenschaften, Bd, 1, p.158.
  !!  - Taken from https://en.wikipedia.org/wiki/Date_of_Easter
  !!
  
  subroutine easter_gauss(year,  month, day)
    implicit none
    integer, intent(in) :: year
    integer, intent(out) :: month, day
    integer :: a,b,c, k,p,q, M,N, d,e
    
    a = modulo(year, 19)
    b = modulo(year,  4)
    c = modulo(year,  7)
    
    if(year>1582) then  ! Gregorian calendar
       k = floor(year / 100.d0)
       p = floor((13 + 8*k) / 25.d0)
       q = floor(k / 4.d0)
       M = modulo(15 + k - p - q, 30)
       N = modulo( 4 + k     - q,  7)
    else  ! Julian calendar
       M = 15
       N = 6
    end if
    
    d = modulo(19*a + M, 30)
    e = modulo(2*b + 4*c + 6*d + N, 7)
    
    month = 4
    day = 22 + d + e
    
    
    ! Sort out the details:
    if((d.eq.29) .and. (e.eq.6)) then  ! Replace 26 with 19 April in some cases
       day = 19
    else if( (d.eq.28) .and. (e.eq.6) .and. (a.gt.10) ) then  ! Replace 25 with 18 April in some cases
       day = 18
    else if(day.gt.31) then  ! If day > 31, we're in April
       day = day - 31
    else                     ! If not, we're in March
       month = 3
    end if
    
  end subroutine easter_gauss
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief   Calculate the date of Passover or Pesach (Jewish Easter) using Gauss' method.
  !! 
  !! \param year   Julian/Gregorian year to compute date of Passover for.
  !! 
  !! \param month  Julian/Gregorian month of year of the DAY of Passover.
  !! \param day    Julian/Gregorian day of month of the DAY of Passover.
  !! 
  !! \note
  !! - Based on Gauss, C.F., 1802. Berechnung des Judischen Osterfestes. Werke, pp.80-81, as described in Meeus, Ch.9.
  !! - Passover is always on 15 Nisan.  Note that this is the DAY of Passover, the day AFTER the evening of the (first) Seder.
  !! - Valid for the Julian and Gregorian calendar.
  !! - The Jewish year is the Julian/Gregorian year + 3760.
  !! - For the years -4000 - +4000, this corresponds EXACTLY to Har’El, Z., 2005. Gauss Formula for the Julian Date
  !!   of Passover. A∆(A), 1, p.0.  Har’El also computes the day of week, which for the same years corresponds exactly
  !!   to dow(cal2jd(year,month,day)) with cal2jd from libSUFR/date_and_time, dow() from libTheSky/datetimeSky, and
  !!   year,month,day as in the current subroutine.  Note that Passover day can only fall on Sun, Tue, Thu or Sat.
  
  subroutine passover_gauss(year,  month, day)
    use SUFR_kinds, only: double
    
    implicit none
    integer, intent(in) :: year
    integer, intent(out) :: month, day
    integer :: C,S, a,b, j
    real(double) :: Q, r
    
    if(year.lt.1583) then         ! Use the Julian calendar
       S = 0
    else                          ! Gregorian calendar
       C = floor(year / 100.d0)
       S = floor((3*C - 5)/4.d0)  !   difference between Julian and Gregorian calendars, in days
    end if
    
    a = modulo(12*year + 12, 19)
    b = modulo(year, 4)
    
    Q = -1.904412361576d0 + 1.554241796621d0*a + 0.25d0*b - 0.003177794022d0*year + S
    j = modulo(floor(Q) + 3*year + 5*b + 2 - S, 7)
    r = Q - floor(Q)  ! Decimal fraction of Q
    
    
    if(j.eq.2 .or. j.eq.4 .or. j.eq.6) then
       day = floor(Q) + 23
    else if(j.eq.1 .and. a.gt.6 .and. r.gt.0.632870370d0) then
       day = floor(Q) + 24
    else if(j.eq.0 .and. a.gt.11 .and. r.gt.0.897723765d0) then
       day = floor(Q) + 23
    else
       day = floor(Q) + 22
    end if
    
    if(day.gt.31) then  ! We're in April
       month = 4
       day = day - 31
    else                ! We're in March
       month = 3
    end if
    
  end subroutine passover_gauss
  !*********************************************************************************************************************************
  
  
end module TheSky_datetime
!***********************************************************************************************************************************

