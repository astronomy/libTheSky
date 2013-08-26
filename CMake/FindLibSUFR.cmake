##  FindLibSUFR.cmake:
##  Check for the presence of the LibSUFR headers and libraries
##  AF, 2010-12-27
##  
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
##  This CMake module defines the following variables:
##    LibSUFR_FOUND        =  Libraries and headers found; TRUE/FALSE
##    LibSUFR_INCLUDES     =  Path to the LibSUFR header files
##    LibSUFR_LIBRARIES    =  Path to all parts of the LibSUFR libraries
##    LibSUFR_LIBRARY_DIR  =  Path to the directory containing the LibSUFR libraries



# Standard locations where to look for required components:
include( CMakeLocations )

if( NOT LibSUFR_FIND_QUIETLY )
  message( STATUS "" )
  message( STATUS "Looking for LibSUFR..." )
endif( NOT LibSUFR_FIND_QUIETLY )


# Check for COMPILER-SPECIFIC header files:
find_path( LibSUFR_INCLUDES 
  NAMES sufr_constants.mod
  PATHS ${include_locations} ${lib_locations}
  PATH_SUFFIXES libSUFR/${Fortran_COMPILER_NAME}
  )

# If not found, check for GENERAL header files:
if( NOT LibSUFR_INCLUDES )
  find_path( LibSUFR_INCLUDES 
    NAMES sufr_constants.mod
    PATHS ${include_locations} ${lib_locations}
    PATH_SUFFIXES libSUFR
    )
endif( NOT LibSUFR_INCLUDES )





# Check for the libraries:
set( LibSUFR_LIBRARIES "" )

# Check for COMPILER-SPECIFIC libraries:
find_library( LibSUFR_LIBRARY
  NAMES SUFR_${Fortran_COMPILER_NAME}
  PATHS ${lib_locations}
  PATH_SUFFIXES libSUFR_${Fortran_COMPILER_NAME} libSUFR
  NO_DEFAULT_PATH
  )

# If not found, check for GENERAL libraries:
if( NOT LibSUFR_LIBRARY )
  find_library( LibSUFR_LIBRARY
    NAMES SUFR
    PATHS ${lib_locations}
    PATH_SUFFIXES libSUFR
    NO_DEFAULT_PATH
    )  
endif( NOT LibSUFR_LIBRARY )

# Libraries found?
if( LibSUFR_LIBRARY )
  
  list( APPEND LibSUFR_LIBRARIES ${LibSUFR_LIBRARY} )
  get_filename_component( LibSUFR_LIBRARY_DIR ${LibSUFR_LIBRARY} PATH )
  
endif( LibSUFR_LIBRARY )





# Headers AND libraries found?
if( LibSUFR_INCLUDES AND LibSUFR_LIBRARIES )
  
  # yes!
  set( LibSUFR_FOUND TRUE )
  
else( LibSUFR_INCLUDES AND LibSUFR_LIBRARIES )
  
  # no!
  set( LibSUFR_FOUND FALSE )
  
  if( NOT LibSUFR_FIND_QUIETLY )
    if( NOT LibSUFR_INCLUDES )
      message( STATUS "!! Unable to find LibSUFR header files!" )
    endif( NOT LibSUFR_INCLUDES )
    if( NOT LibSUFR_LIBRARIES )
      message( STATUS "!! Unable to find LibSUFR library files!" )
    endif( NOT LibSUFR_LIBRARIES )
  endif( NOT LibSUFR_FIND_QUIETLY )
  
endif( LibSUFR_INCLUDES AND LibSUFR_LIBRARIES )




# Headers AND libraries found!
if( LibSUFR_FOUND )
  
  if( NOT LibSUFR_FIND_QUIETLY )
    message( STATUS "Found components for LibSUFR:" )
    message( STATUS "* LibSUFR_INCLUDES  = ${LibSUFR_INCLUDES}" )
    message( STATUS "* LibSUFR_LIBRARIES = ${LibSUFR_LIBRARIES}" )
  endif( NOT LibSUFR_FIND_QUIETLY )
  
else( LibSUFR_FOUND )
  
  if( LibSUFR_FIND_REQUIRED )
    message( FATAL_ERROR "!! Could not find LibSUFR headers or libraries!" )
  else( LibSUFR_FIND_REQUIRED )
    message( STATUS "!! Could not find LibSUFR headers or libraries!" )
  endif( LibSUFR_FIND_REQUIRED )
  
endif( LibSUFR_FOUND )



# Mark as advanced options in ccmake:
mark_as_advanced( 
  LibSUFR_INCLUDES
  LibSUFR_LIBRARIES
  LibSUFR_LIBRARY
  LibSUFR_LIBRARY_DIR
  )

