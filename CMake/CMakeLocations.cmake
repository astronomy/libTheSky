##  CMakeLocations.cmake:
##  Set search locations for libraries and header files on Linux/Unix/BSD/MacOSX
##  Adapted from CMakeSettings.cmake 2818 2009-07-15 19:02:31Z baehren
##    by Marc van der Sluys - marc.vandersluys.nl
##  
##  Copyright (C) 2007 Lars B"ahren (bahren@astron.nl)
##  Copyright 2010-2013 Marc van der Sluys - marc.vandersluys.nl
##   
##  This file is part of the CMakeFiles package, 
##  see: http://cmakefiles.sf.net/
##   
##  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
##  
##  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public License along with this code.  If not, see 
##  <http://www.gnu.org/licenses/>.
##
##
##  Variables used through the configuration environment:
##    search_locations    --   general files (?)
##    bin_locations       --   executables
##    include_locations   --   header and module files
##    lib_locations       --   libraries


# User's home directory:
set( HOME $ENV{HOME} )
  

set( search_locations
  /usr
  /usr/local
  /opt
  /opt/local
  /sw
  ${HOME}/usr
  CACHE
  PATH
  "Directories to look for include files"
  FORCE
  )


## ---------------------------------------------------------------------------
## locations in which to look for applications/binaries

set( bin_locations
  /usr/bin
  /usr/local/bin
  /sw/bin
  ${HOME}/usr/bin
  ${HOME}/bin
  CACHE
  PATH
  "Extra directories to look for executable files"
  FORCE
  )


## ----------------------------------------------------------------------------
## locations in which to look for header and module files

set( include_locations
  /usr/include
  /usr/local/include
  /opt/include
  /opt/local/include
  /sw/include
  ${HOME}/usr/include
  ${HOME}/include
  CACHE
  PATH
  "Directories to look for include files"
  FORCE
  )


## ----------------------------------------------------------------------------
## locations in which to look for libraries

set( lib_locations
  /usr/local/lib64
  /usr/local/lib
  /usr/lib64
  /usr/lib
  /opt/lib
  /opt/local/lib
  /sw/lib
  ${HOME}/usr/lib
  ${HOME}/lib
  CACHE
  PATH
  "Directories to look for libraries"
  FORCE
  )


#message( "${search_locations}" )
#message( "${bin_locations}" )
#message( "${include_locations}" )
#message( "${lib_locations}" )

