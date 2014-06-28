!> \file data.f90  Procedures to define constants and read data files for libTheSky


!  Copyright (c) 2002-2014  Marc van der Sluys - marc.vandersluys.nl
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
!> \brief Procedures to set constants and read data files

module TheSky_data
  implicit none
  save
  
  
contains


  !*********************************************************************************************************************************
  !> \brief Set the values of the constants used in this package
  
  subroutine set_TheSky_constants()
    use SUFR_constants, only: set_SUFR_constants, homedir
    use SUFR_system, only: error,quit_program_error
    use TheSky_constants, only: library_name, TheSkydir
    
    implicit none
    integer, parameter :: ntrials = 11
    integer :: try
    character :: tryPaths(ntrials)*(99)
    logical :: ex
    
    call set_SUFR_constants()  ! Set libSUFR constants
    
    library_name = 'libTheSky (libthesky.sourceforge.net)'
    
    
    ! Find the libTheSky path:
    trypaths(1)  = trim(TheSkydir)  ! Predefined?
    if(len_trim(TheSkydir).eq.0 .or. len_trim(TheSkydir).eq.len(TheSkydir)) trypaths(1) = './'  ! Current working dir
    
    ! Try the user's home directory first:
    trypaths(2)  = trim(homedir)//'/share/libTheSky/'
    trypaths(3)  = trim(homedir)//'/local/share/libTheSky/'
    trypaths(4)  = trim(homedir)//'/usr/share/libTheSky/'
    trypaths(5)  = trim(homedir)//'/usr/local/share/libTheSky/'
    trypaths(6)  = trim(homedir)//'/opt/share/libTheSky/'
    trypaths(7)  = trim(homedir)//'/opt/local/share/libTheSky/'
     
    ! Try default system paths:
    trypaths(8)  = '/usr/share/libTheSky/'
    trypaths(9)  = '/usr/local/share/libTheSky/'
    trypaths(10) = '/opt/share/libTheSky/'
    trypaths(11) = '/opt/local/share/libTheSky/'
    
    
    do try = 1,ntrials
       TheSkydir = tryPaths(try)
       inquire(file=trim(TheSkydir)//'planets.dat', exist=ex)
       if(ex) exit
    end do
    
    if(.not.ex) then
       call error('I could not find the libTheSky data directory.  I tried the following possibilities:', 0)
       do try= 1,ntrials
          write(0,'(I4,2x,A)') try, trim(trypaths(try))
       end do
       call quit_program_error('libTheSky: Data directory not found.', 1)
    end if
    
  end subroutine set_TheSky_constants
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Read input data files
  
  subroutine TheSky_readdata()
    implicit none
    
    call read_deltat()
    
    call readPlanetData()
    call readPluto()
    
    call readmoonpos_la()
    
    call readplanetelements()
    
    call readAsteroidElements()
    call readCometElements()
    call readNutation()
    call readConstellations()
    
  end subroutine TheSky_readdata
  !*********************************************************************************************************************************
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief Read deltat.dat with historical values for DeltaT
  !! 
  !! - Is needed by read_locations and hence called from there
  !! \see http://hemel.waarnemen.com/Computing/deltat.html
  
  subroutine read_deltat()
    use SUFR_system, only: find_free_io_unit, file_open_error_quit, file_read_error_quit
    use SUFR_dummy, only: dumstr
    use TheSky_constants, only: TheSkydir, deltat_0, deltat_accel, deltat_change, deltat_minyr, deltat_maxyr, deltat_n, deltat_nmax
    use TheSky_constants, only: deltat_years, deltat_values
    
    implicit none
    integer :: i,status, dn, ip
    character :: infile*(199)
    
    call find_free_io_unit(ip)
    infile = trim(TheSkydir)//'deltat.dat'
    open(ip, form='formatted', status='old', action='read', file=trim(infile), iostat=status)
    if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
    
    do i=1,7
       read(ip,*) dumstr
    end do
    
    do i=1,deltat_nmax
       read(ip,*, iostat=status) deltat_years(i),deltat_values(i)
       if(status.lt.0) exit
       if(status.gt.0) call file_read_error_quit(trim(infile), i, 1)  ! Exit status 1
    end do
    
    close(ip)
    
    deltat_n = i-1
    deltat_minyr = nint(deltat_years(1))
    deltat_maxyr = nint(deltat_years(deltat_n))
    
    
    ! Set the average acceleration of DeltaT, from -1500 to 2010 CE:
    !   see http://hemel.waarnemen.com/Computing/deltat.html
    deltat_accel = 3.975539d-3  ! Seconds / yr^2:  39.755 +- 0.048 s/cty^2 (fit.f90)
    
    ! Compute the average change in DeltaT in the last dn years:
    dn = 100
    deltat_change = (deltat_values(deltat_n)-deltat_values(deltat_n-dn)) / (deltat_years(deltat_n)-deltat_years(deltat_n-dn))
    
    ! Take the most recent tabulated value of DeltaT:
    deltat_0 = deltat_values(deltat_n)
    
  end subroutine read_deltat
  !*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Read constellation names: abbreviation, Latin, genitive, Dutch and English
  
  subroutine readconstellations()
    use SUFR_system, only: find_free_io_unit, file_open_error_quit, file_read_end_error
    use SUFR_dummy, only: dumstr9
    
    use TheSky_constants, only: TheSkydir
    use TheSky_stardata, only: conabr, conidral, conid, nconid, conidrau, conidabr, coniddecl
    use TheSky_stardata, only: nconstel, latconnames, genconnames, nlconnames, enconnames
    
    implicit none
    integer :: i,j, ip,status
    character :: infile*(199)
    
    ! Read constellation names (abbr., Lat, Lat.Gen., Nl, En):
    call find_free_io_unit(ip)
    infile = trim(TheSkydir)//'const_names.dat'
    open(ip, status='old', action='read', file=trim(infile), iostat=status)
    if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
    
    read(ip,*) dumstr9
    read(ip,*) dumstr9
    
    do i=1,nconstel
       read(ip,'(A3,2x,A19,2x,A19,2x,A17,2x,A18)', iostat=status) conabr(i),latconnames(i),genconnames(i), &
            nlconnames(i),enconnames(i)
       
       if(status.ne.0) call file_read_end_error(trim(infile), i, status, 1, 1)  ! stopcode=1, exitstatus=1
    end do
    
    close(ip)
    
    
    ! Read constellation ID data:
    infile = trim(TheSkydir)//'const_id.dat'
    open(ip, status='old', action='read', file=trim(infile), iostat=status)
    if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
    
    do i=1,nconid
       read(ip,'(2F8.4,F9.4,1x,A3)', iostat=status) conidral(i),conidrau(i),coniddecl(i),conidabr(i)
       if(status.ne.0) call file_read_end_error(trim(infile), i, status, 1, 1)  ! stopcode=1, exitstatus=1
       
       do j=1,nconstel
          if(conabr(j).eq.conidabr(i)) conid(i) = j  ! Assign a constellation number ID rather than abbrev.
       end do
    end do
    
    close(ip)
  end subroutine readconstellations
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Read nutation input files
  
  subroutine readnutation()
    use SUFR_system, only: find_free_io_unit, file_open_error_quit, file_read_end_error
    use SUFR_kinds, only: double
    use TheSky_constants, only: TheSkydir, nutationdat
    
    implicit none
    integer :: i,nu1(6),nu3, ip,status
    real(double) :: nu2,nu4
    character :: infile*(199)
    
    call find_free_io_unit(ip)
    infile = trim(TheSkydir)//'nutation.dat'
    open(ip, form='formatted', status='old', action='read', file=trim(infile), iostat=status)
    if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
    
    do i=1,63
       read(ip,'(5(I2),2x,I7,2x,F6.1,2x,I6,2x,F4.1)', iostat=status) nu1,nu2,nu3,nu4
       if(status.ne.0) call file_read_end_error(trim(infile), i, status, 1, 1)  ! stopcode=1, exitstatus=1
       
       nutationdat(1:6,i) = nu1
       nutationdat(7,i)   = nu2
       nutationdat(8,i)   = nu3
       nutationdat(9,i)   = nu4
    end do
    close(ip)
    
  end subroutine readnutation
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Reads VSOP planet data from planets.dat
  
  subroutine readplanetdata()
    use SUFR_system, only: find_free_io_unit, file_open_error_quit, file_read_end_error
    use TheSky_constants, only: TheSkydir
    use TheSky_planetdata, only: VSOPnls, VSOPdat
    
    implicit none
    integer :: pl, li, ip, status,  powr, TOTnls(8)
    character :: infile*(199)
    
    VSOPnls = reshape( (/ 2808,1620,2399, 671,426,585, 1080,349,997,  &    ! Number of lines in VSOP input files (l,b,r x 8 pl)
         2393,915,2175, 1484,530,1469, 2358,966,2435, 1578,516,1897, 681,290,959/),  (/3,8/))
    
    TOTnls(1:8) = sum(VSOPnls(:,1:8), 1)  ! Total number of lines per planet; sum along dimension 1
    
    
    ! Read planets.dat file:
    call find_free_io_unit(ip)
    infile = trim(TheSkydir)//'planets.dat'
    open(ip, form='formatted', status='old', action='read', file=trim(infile), iostat=status)
    if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
    
    VSOPdat = 0.d0
    do pl=1,8  ! Planet
       do li=1,TOTnls(pl)
          read(ip,'(I1,F18.11,F14.11,F20.11)', iostat=status) powr, VSOPdat(2:4,li,pl)
          if(status.ne.0) call file_read_end_error(trim(infile), li, status, 1, 1)      ! stopcode=1, exitstatus=1
          VSOPdat(1,li,pl) = dble(powr)
       end do  ! li
    end do  ! pl
    
    close(ip)
    
  end subroutine readplanetdata
  !*********************************************************************************************************************************
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Read periodic terms for the position of Pluto
  
  subroutine readpluto()
    use SUFR_system, only: find_free_io_unit, file_open_error_quit, file_read_end_error
    use TheSky_planetdata, only: pluc,plul,plub,plur
    use TheSky_constants, only: TheSkydir
    
    implicit none
    integer :: i, ip, status
    character :: infile*(199)
    
    call find_free_io_unit(ip)
    infile = trim(TheSkydir)//'pluto.dat'
    open(ip, form='formatted', status='old', action='read', file=trim(infile), iostat=status)
    if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
    
    do i=1,43
       read(ip,'(3I3,6I10)', iostat=status) pluc(i,1), pluc(i,2), pluc(i,3), plul(i,1), plul(i,2), plub(i,1), plub(i,2), &
            plur(i,1), plur(i,2)
       
       if(status.ne.0) call file_read_end_error(trim(infile), i, status, 1, 1)  ! stopcode=1, exitstatus=1
    end do
    close(ip)
    
  end subroutine readpluto
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Read the data needed to compute the orbital elements of the planets
  !!
  !! -  use planetelements() to compute them
  !! 
  !! - L, a, e, i, Omega, pi for Mer-Nep, mean equinox of date and J2000.0 (=2x6x8)
  !! - Terms in one row are a_0, a_1, a_2, a_3 and should be used as:  Sum(i=0,3) a_i t^i, with t the time in centuries since 2000.0
  !! - Angles are still expressed in degrees
  
  subroutine readplanetelements()
    use SUFR_system, only: find_free_io_unit, file_open_error_quit, file_read_end_error
    use TheSky_constants, only: TheSkydir
    use TheSky_planetdata, only: plelemdata
    
    implicit none
    integer :: eq,pl,el, ip,status
    character :: infile*(199)
    
    call find_free_io_unit(ip)
    infile = trim(TheSkydir)//'planetelements.dat'
    open(ip, form='formatted', status='old', action='read', file=trim(infile), iostat=status)
    if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
    
    outerloop: do eq=1,2  ! Equinoctes
       do pl=1,8  ! Planets
          do el=1,6  ! Elements
             if(eq.eq.2.and.(el.eq.2.or.el.eq.3)) cycle  ! a, e are not in the second part (J2000.0?)
             
             read(ip,'(F13.9, F18.10, 2F15.11)', iostat=status) plelemdata(eq, pl, el, 0:3)
             
             if(status.ne.0) then
                call file_read_end_error(trim(infile), 0, status, 0, 0)  ! line=0, stopcode=0, exitstatus=0
                exit outerloop
             end if
             
          end do  ! el
       end do  ! pl
    end do outerloop  ! eq
    close(ip)
    
    plelemdata(2,1:8,2,0:3) = plelemdata(1,1:8,2,0:3)  ! a_J2000 = a
    plelemdata(2,1:8,3,0:3) = plelemdata(1,1:8,3,0:3)  ! e_J2000 = e
    
  end subroutine readplanetelements
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Read ELP periodic terms from moondata.dat
  !!
  !! \see ftp://cdsarc.u-strasbg.fr/pub/cats/VI/79/
  !!
  !! Constants are defined in the module moondata
  
  subroutine readmoondata()
    use SUFR_kinds, only: double, dbl
    use SUFR_constants, only: pi, pio2, d2r, as2r, r2as
    use SUFR_system, only: find_free_io_unit,file_open_error_quit, file_read_end_error
    use SUFR_numerics, only: dne
    
    use TheSky_constants, only: TheSkydir
    use TheSky_moondata, only: nterm,nrang, pc1,pc2,pc3, per1,per2,per3, w, ath
    use TheSky_moondata, only: eart,peri,p,del,zeta,pre, coef,zone,ilu,ipla,ideb,c1,c2
    use TheSky_moondata, only: p1,p2,p3,p4,p5,q1,q2,q3,q4,q5, prec0
    
    implicit none
    integer :: i,ific,ir,itab,iv,j,k,iz, ip,status,  filen(36),fl,il
    real(double) :: prec,am,alpha,dtasm,precess,  delnu,dele,delg,delnp,delep,  xx,tgv,y,pha
    character :: infile*(199)
    
    filen = (/1023,918,704,347,316,237,14,11,8,14328,5233,6631,4384, &
         833,1715,170,150,114,226,188,169,3,2,2,6,4,5,20,12,14,11,4,10,28,13,19/)
    prec = 0.0_dbl
    
    ! Make sure these are always defined:
    delep = 0.0_dbl
    delnp = 0.0_dbl
    delg  = 0.0_dbl
    dele  = 0.0_dbl
    delnu = 0.0_dbl
    dtasm = 0.0_dbl
    am    = 0.0_dbl
    
    !*** Parameters:
    if(ideb.eq.0) then  ! then it's the first time you run this routine
       
       ideb = 1
       am = 0.074801329518d0
       alpha = 0.002571881335d0
       dtasm = 2.d0*alpha/(3.d0*am)
       
       ! Lunar arguments:
       w(1,0)  = (218.d0+18.d0*c1 + 59.95571d0*c2) * d2r
       w(2,0)  = (83.d0+21.d0*c1  + 11.67475d0*c2) * d2r
       w(3,0)  = (125.d0+2.d0*c1  + 40.39816d0*c2) * d2r
       eart(0) = (100.d0+27.d0*c1 + 59.22059d0*c2) * d2r
       peri(0) = (102.d0+56.d0*c1 + 14.42753d0*c2) * d2r
       
       w(1,1)  =  1732559343.73604d0 * as2r
       w(2,1)  =  14643420.2632d0    * as2r
       w(3,1)  = -6967919.3622d0     * as2r
       eart(1) =  129597742.2758d0   * as2r
       peri(1) =  1161.2283d0        * as2r
       
       w(1,2)  = -5.8883d0           * as2r
       w(2,2)  = -38.2776d0          * as2r
       w(3,2)  =  6.3622d0           * as2r
       eart(2) = -0.0202d0           * as2r
       peri(2) =  0.5327d0           * as2r
       
       w(1,3)  =  0.6604d-2          * as2r
       w(2,3)  = -0.45047d-1         * as2r
       w(3,3)  =  0.7625d-2          * as2r
       eart(3) =  0.9d-5             * as2r
       peri(3) = -0.138d-3           * as2r
       
       w(1,4)  = -0.3169d-4          * as2r
       w(2,4)  =  0.21301d-3         * as2r
       w(3,4)  = -0.3586d-4          * as2r
       eart(4) =  0.15d-6            * as2r
       peri(4) =  0.d0
       
       
       ! Precession constant:
       precess = 5029.0966d0 * as2r
       
       ! Planetary arguments:
       p(1,0) = (252.d0+15.d0*c1+3.25986d0*c2)  * d2r
       p(2,0) = (181.d0+58.d0*c1+47.28305d0*c2) * d2r
       p(3,0) = eart(0)
       p(4,0) = (355.d0+25.d0*c1+59.78866d0*c2) * d2r
       p(5,0) = (34.d0+21.d0*c1+5.34212d0*c2)   * d2r
       p(6,0) = (50.d0+4.d0*c1+38.89694d0*c2)   * d2r
       p(7,0) = (314.d0+3.d0*c1+18.01841d0*c2)  * d2r
       p(8,0) = (304.d0+20.d0*c1+55.19575d0*c2) * d2r
       p(1,1) = 538101628.68898d0               * as2r
       p(2,1) = 210664136.43355d0               * as2r
       p(3,1) = eart(1)
       p(4,1) = 68905077.59284d0                * as2r
       p(5,1) = 10925660.42861d0                * as2r
       p(6,1) = 4399609.65932d0                 * as2r
       p(7,1) = 1542481.19393d0                 * as2r
       p(8,1) = 786550.32074d0                  * as2r
       
       ! Corrections of the constants (fit to DE200/LE200):
       delnu = +0.55604d0/w(1,1)                * as2r
       dele  = +0.01789d0                       * as2r
       delg  = -0.08066d0                       * as2r
       delnp = -0.06424d0/w(1,1)                * as2r
       delep = -0.12879d0                       * as2r
       
       ! Delaunay's arguments:
       do i=0,4
          del(1,i) = w(1,i)  - eart(i)
          del(4,i) = w(1,i)  - w(3,i)
          del(3,i) = w(1,i)  - w(2,i)
          del(2,i) = eart(i) - peri(i)
       end do
       
       del(1,0) = del(1,0) + pi
       zeta(0)  = w(1,0)
       zeta(1)  = w(1,1) + precess
       
       ! Precession matrix:
       p1 =  0.10180391d-4
       p2 =  0.47020439d-6
       p3 = -0.5417367d-9
       p4 = -0.2507948d-11
       p5 =  0.463486d-14
       q1 = -0.113469002d-3
       q2 =  0.12372674d-6
       q3 =  0.1265417d-8
       q4 = -0.1371808d-11
       q5 = -0.320334d-14
       
    end if  ! if(ideb.eq.0)
    
    
    
    ! Read moondata file:
    call find_free_io_unit(ip)
    infile = trim(TheSkydir)//'moondata.dat'
    open(ip, form='formatted', status='old', action='read', file=trim(infile), iostat=status)
    if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
    
    if(dne(prec,prec0)) then
       prec0 = prec
       pre(1) = prec * r2as  -  1.d-12
       pre(2) = prec * r2as  -  1.d-12
       pre(3) = prec * ath
       
       do ific=1,36  ! Loop over all original ELP files
          fl = filen(ific)
          ir = 0
          itab = (ific+2)/3
          iv = mod(ific-1,3) + 1
          
          select case(ific)
          case(1:3)  ! Files: Main problem
             do il=1,fl
                read(ip,'(4I3,2x,F13.5,6(2x,F10.2))', iostat=status) ilu,coef
                if(status.ne.0) call file_read_end_error(trim(infile), i, status, 1, 1)  ! stopcode=1, exitstatus=1
                
                ir = ir+1
                tgv = coef(2)+dtasm*coef(6)
                if(ific.eq.3) coef(1) = coef(1)-2.d0*coef(1)*delnu/3.d0
                xx = coef(1) + tgv*(delnp-am*delnu) + coef(3)*delg + coef(4)*dele + coef(5)*delep
                zone(1) = xx
                do k=0,4
                   y = 0.d0
                   do i=1,4
                      y = y + ilu(i)*del(i,k)
                   end do
                   zone(k+2) = y
                end do
                if(iv.eq.3) zone(2) = zone(2) + pio2
                do i=1,6
                   if(iv.eq.1) pc1(i,ir) = zone(i)
                   if(iv.eq.2) pc2(i,ir) = zone(i)
                   if(iv.eq.3) pc3(i,ir) = zone(i)
                end do
             end do  ! do il=1,fl
             
             
          case(4:9,22:36)  ! Files: Tides - Relativity - Solar eccentricity
             do il=1,fl
                read(ip,'(5I3,1x,F9.5,1x,F9.5)', iostat=status) iz,ilu,pha,xx
                if(status.ne.0) call file_read_end_error(trim(infile), i, status, 1, 1)  ! stopcode=1, exitstatus=1
                
                ir = ir + 1
                zone(1) = xx
                do k=0,1
                   if(k.eq.0) y = pha*d2r
                   if(k.ne.0) y = 0.d0
                   y = y+iz*zeta(k)
                   do i=1,4
                      y = y+ilu(i)*del(i,k)
                   end do
                   zone(k+2) = y
                end do
                j = nrang(iv,itab-1)+ir
                do i=1,3
                   if(iv.eq.1) per1(i,j) = zone(i)
                   if(iv.eq.2) per2(i,j) = zone(i)
                   if(iv.eq.3) per3(i,j) = zone(i)
                end do
             end do  ! do il=1,fl
             
             
          case(10:21)  ! Files: Planetary perturbations
             do il=1,fl
                read(ip,'(11I3,1x,F9.5,1x,F9.5)', iostat=status) ipla,pha,xx
                if(status.ne.0) call file_read_end_error(trim(infile), i, status, 1, 1)  ! stopcode=1, exitstatus=1
                
                ir = ir+1
                zone(1) = xx
                if(ific.lt.16) then
                   do k=0,1
                      if(k.eq.0) y = pha*d2r
                      if(k.ne.0) y = 0.d0
                      y = y + ipla(9)*del(1,k) + ipla(10)*del(3,k) + ipla(11)*del(4,k)
                      do i=1,8
                         y = y+ipla(i)*p(i,k)
                      end do
                      zone(k+2) = y
                   end do
                else
                   do k=0,1
                      if(k.eq.0) y = pha*d2r
                      if(k.ne.0) y = 0.d0
                      do i=1,4
                         y = y+ipla(i+7)*del(i,k)
                      end do
                      do i=1,7
                         y = y+ipla(i)*p(i,k)
                      end do
                      zone(k+2) = y
                   end do
                end if
                j = nrang(iv,itab-1)+ir
                do i=1,3
                   if(iv.eq.1) per1(i,j) = zone(i)
                   if(iv.eq.2) per2(i,j) = zone(i)
                   if(iv.eq.3) per3(i,j) = zone(i)
                end do
             end do  ! do il=1,fl
             
          end select
          
          nterm(iv,itab) = ir
          if(itab.eq.1) then
             nrang(iv,itab) = 0
          else
             nrang(iv,itab) = nrang(iv,itab-1)+nterm(iv,itab)
          end if
          
       end do   ! do ific=1,36 - original ELP files
       
    end if  ! if(prec.ne.prec0)
    close(ip)
    
  end subroutine readmoondata
  !*********************************************************************************************************************************
  
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Read low-accuracy moon-position data from Meeus (moon_la.dat)
  
  subroutine readmoonpos_la()
    use SUFR_system, only: find_free_io_unit, file_open_error_quit, file_read_end_error
    use SUFR_dummy, only: dumstr
    use TheSky_constants, only: TheSkydir
    use TheSky_planetdata, only: moonla_arg,moonla_lrb
    
    implicit none
    integer :: i, ip, status
    character :: infile*(199)
    
    call find_free_io_unit(ip)
    infile = trim(TheSkydir)//'moon_la.dat'
    open(ip, form='formatted', status='old', action='read', file=trim(infile), iostat=status)
    if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
    
    read(ip,*) dumstr
    read(ip,*) dumstr
    read(ip,*) dumstr
    
    do i=1,60
       read(ip,'(I1,3I3,2I10)', iostat=status) moonla_arg(1:4,i),moonla_lrb(1:2,i)  ! Read terms for L and R
       if(status.ne.0) call file_read_end_error(trim(infile), 3+i, status, 1, 1)  ! stopcode=1, exitstatus=1
    end do
    
    read(ip,*) dumstr
    do i=1,60
       read(ip,'(I1,3I3,2I10)', iostat=status) moonla_arg(5:8,i),moonla_lrb(3,i)    ! Read terms for B
       if(status.ne.0) call file_read_end_error(trim(infile), 64+i, status, 1, 1)  ! stopcode=1, exitstatus=1
    end do
    close(ip)
    
  end subroutine readmoonpos_la
  !*********************************************************************************************************************************
  
  
  

  !*********************************************************************************************************************************
  !> \brief  Read orbital-element data for the asteroids
  !!
  !! \note 
  !! - data are passed via the module planetdata, in asterElems(*,1:9):
  !! - Epoch (JD), a, e, i, omega, Omega, M, H, G, for J2000.0
  !!
  !! \see
  !! - asteroids.dat: http://sf.net/projects/libthesky/files/asteroids.dat.bz2 (selection: H<15, a<100 - 25% of all bodies)
  !! - original: http://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR.gz
  
  subroutine readasteroidelements()
    use SUFR_constants, only: d2r
    use SUFR_system, only: find_free_io_unit, file_read_end_error
    use SUFR_dummy, only: dumstr9
    use TheSky_constants, only: TheSkydir
    use TheSky_planetdata, only: asternames, nasteroids, asterelems
    
    implicit none
    integer :: i,epoch, ip,status
    character :: infile*(199)
    
    call find_free_io_unit(ip)
    infile = trim(TheSkydir)//'asteroids.dat'  ! (Reduced) copy of http://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR.gz
    open(ip, form='formatted', status='old', action='read', file=trim(infile), iostat=status)
    if(status.ne.0) then
       write(0,'(/,A)') '  I could not open the file '//trim(infile)//' - you will not be able to compute asteroid data.'
       write(0,'(A,/)') '  Download http://sf.net/projects/libthesky/files/asteroids.dat.bz2 - unzip it and move it to '// &
            trim(infile)//' to fix this.'
       return
    end if
    
    read(ip,*) dumstr9
    read(ip,*) dumstr9
    
    do i=1,nasteroids
       read(ip,'(7x,A18,I5,F11.7,F11.8,3F10.5,F12.7,F6.2,F5.2)', iostat=status) asternames(i),epoch,asterelems(i,2:9)
       if(status.ne.0) call file_read_end_error(trim(infile), 2+i, status, 1, 1)  ! stopcode=1, exitstatus=1
       
       if(status.ne.0) then
          call file_read_end_error(trim(infile), 0, status, 0, 0)  ! line=0, stopcode=0, exitstatus=0
          exit
       end if
       
       asterelems(i,1) = dble(epoch) + 2400000.5d0
    end do
    
    close(ip)
    
    asterelems(1:nasteroids,4:7) = asterelems(1:nasteroids,4:7)*d2r
    
  end subroutine readasteroidelements
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Read orbital-element data for the comets
  !!
  !! \note
  !! - data are passed via the module comet_data, in cometElems(*,1:9)
  !! - colums:  1: Epoch (JD),  2: q,  3: e,  4: i,  5: omega,  6: Omega,  7: Tp, for J2000.0, 8: H, 9: G
  !!
  !! \see
  !! - comets.dat: http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET
  !! - comets_mpc.dat: http://www.minorplanetcenter.net/iau/Ephemerides/Comets/Soft00Cmt.txt
  
  subroutine readCometElements()
    use SUFR_kinds, only: double
    use SUFR_constants, only: d2r, currentJD
    use SUFR_system, only: find_free_io_unit, quit_program_error, file_open_error_quit, file_read_error_quit
    use SUFR_date_and_time, only: cal2jd
    use SUFR_dummy, only: dumstr9
    
    use TheSky_constants, only: TheSkydir
    use TheSky_cometdata, only: cometNames, cometElems, nCometsMax,nComets, cometDatFile, cometDiedAtP
    
    implicit none
    integer :: ci,epoch, yr,mnt,idy, status, ip
    real(double) :: dy
    character :: tpstr*(14),epstr*(8), infile*(199)
    
    cometDatFile = 2  ! 1: comets.dat (MANY comets, no magnitude info), 2: comets_mpc.dat (currently visible comets + magn. info)
    
    cometElems = 0.d0
    cometDiedAtP = 0   ! By default, comets do not die at perihelion - set to 1 for a given comet to change this
    nComets = 0
    
    call find_free_io_unit(ip)
    select case(cometDatFile)
    case(1)  ! comets.dat
       infile = trim(TheSkydir)//'comets.dat'  ! Copy of http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET
       open(ip, form='formatted', status='old', action='read', file=trim(infile), iostat=status)
       if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
       
       read(ip,*) dumstr9
       read(ip,*) dumstr9
       
       do ci=11,nCometsMax
          read(ip,'(A43,I8,F12.8,F11.8,3F10.5,1x,A14)', iostat=status) cometNames(ci), epoch, cometElems(ci,2), cometElems(ci,3), &
               cometElems(ci,4), cometElems(ci,5), cometElems(ci,6), tpstr
          
          if(status.lt.0) exit
          if(status.gt.0) call file_read_error_quit(trim(infile), ci-10+2, 1)
          
          nComets = nComets + 1
          cometElems(ci,1) = dble(epoch) + 2400000.5d0  ! JD of epoch
          read(tpstr,'(I4,I2,F8.5)', iostat=status) yr,mnt,dy
          if(status.ne.0) call quit_program_error('readCometElements():  Error reading perihelion date from string', 1)
          cometElems(ci,7) = cal2jd(yr,mnt,dy)  ! JD of perihelion
       end do
       
       close(ip)
       
    case(2)  ! comets_mpc.dat
       infile = trim(TheSkydir)//'comets_mpc.dat'  ! Copy of http://www.minorplanetcenter.net/iau/Ephemerides/Comets/Soft00Cmt.txt
       open(ip, form='formatted', status='old', action='read', file=trim(infile), iostat=status)
       if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
       
       do ci=11,nCometsMax
          read(ip,'(A12, I6,I3,F8.4, 2F10.6,3F10.4, 2x,A8, F6.1,F5.1,2x,A56)', iostat=status) dumstr9, yr,mnt,dy, &
               cometElems(ci,2:3), cometElems(ci,5:6),cometElems(ci,4), epstr, cometElems(ci,8:9), cometNames(ci)
          
          if(status.lt.0) exit
          if(status.gt.0) call file_read_error_quit(trim(infile), ci-10, 1)
          
          nComets = nComets + 1
          cometElems(ci,7) = cal2jd(yr,mnt,dy)  ! JD of perihelion
          
          read(epstr,'(I4,I2,I2)', iostat=status) yr,mnt,idy
          if(status.ne.0) call quit_program_error('readCometElements():  Error reading epoch from string', 1)
          cometElems(ci,1) = cal2jd(yr,mnt,dble(idy))  ! JD of epoch
          if(len_trim(epstr).ne.8) cometElems(ci,1) = currentJD  ! epstr is sometimes empty
       end do
       
       close(ip)
       
    case default
       call quit_program_error('readCometElements():  non-exisiting value of cometDatFile',0)
       
    end select
    
    
    cometElems(1:nComets,4:6) = cometElems(1:nComets,4:6)*d2r  ! omega, Omega, i
    
  end subroutine readCometElements
  !*********************************************************************************************************************************
  
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Read the Bright Star Catalogue
  
  subroutine read_bsc
    use SUFR_system, only: find_free_io_unit, file_open_error_quit, file_read_error_quit, file_read_end_error
    use SUFR_sorting, only: sorted_index_list
    use TheSky_constants, only: TheSkydir
    use TheSky_BSCdata, only: n_bsc, bsc_abbr,bsc_ra,bsc_dec,bsc_pma,bsc_pmd,bsc_rv,bsc_vm,bsc_par,bsc_mult,bsc_var
    use TheSky_BSCdata, only: bsc_sptype,bsc_sao, bsc_bv,bsc_ub,bsc_ri, bsc_vm_indx,bsc_name
    
    
    implicit none
    integer :: i,rv,nsn,snr, ip,status
    character :: snam*(10), infile*(199)
    
    call find_free_io_unit(ip)
    infile = trim(TheSkydir)//'bsc_data.dat'
    open(ip, action='read', form='formatted', status='old', file=trim(infile), iostat=status)
    if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
    
    do i=1,n_bsc
       read(ip,'(A10,1x,2F10.6,1x,2F7.3,I5,F6.2,F6.3,A2,A10,2x,A20,I8,1x,3F6.2)', iostat=status) &
            bsc_abbr(i),bsc_ra(i),bsc_dec(i),bsc_pma(i),bsc_pmd(i),rv,bsc_vm(i),bsc_par(i),bsc_mult(i),bsc_var(i), &
            bsc_sptype(i),bsc_sao(i),bsc_bv(i),bsc_ub(i),bsc_ri(i)
       
       if(status.ne.0) call file_read_end_error(trim(infile), i, status, 1, 1)  ! stopcode=1, exitstatus=1
       
       if(abs(bsc_vm(i)).lt.1.d-5) bsc_vm(i) = 99.d0   ! There are some non-stellar objects in the catalogue, with vm=0.0
       bsc_rv(i) = dble(rv)
    end do
    close(ip)
    
    ! Create an index sorted to visual magnitude:
    call sorted_index_list(bsc_vm(1:n_bsc),bsc_vm_indx(1:n_bsc))  ! size n_bsc
    
    
    ! Star names:
    bsc_name = '          '
    nsn = 79  ! Number of names in bsc_names.dat: actually 80, but Polaris is listed twice
    infile = trim(TheSkydir)//'bsc_names.dat'
    open(ip, action='read', form='formatted', status='old', file=trim(infile), iostat=status)
    if(status.ne.0) call file_open_error_quit(trim(infile), 1, 1)  ! 1-input file, 1-exit status
    
    do i=1,nsn
       read(ip,'(I4,2x,A10)', iostat=status) snr, snam
       if(status.lt.0) exit
       if(status.gt.0) call file_read_error_quit(trim(infile), i, 1)
       
       bsc_name(snr) = snam
    end do
    close(ip)
    
  end subroutine read_bsc
  !*********************************************************************************************************************************
  
  
  
end module TheSky_data
!***********************************************************************************************************************************
