##
##  CMakeFile for libTheSky
##  AF, 2013-06-27
##  
##  Copyright 2010-2013 Marc van der Sluys - marc.vandersluys.nl
##   
##  This file is part of the libTheSky package, 
##  see: https://libthesky.sf.net/
##   
##  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
##  
##  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public License along with this code.  If not, see 
##  <https://www.gnu.org/licenses/>.
##
##
##  To compile, from the directory that contains this file, do:
##    $ mkdir build; cd build
##    $ cmake ..
##    $ make install


cmake_minimum_required( VERSION 2.8 FATAL_ERROR )

# Set build type. Do this *before* we set the project name:
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo Profile."
    FORCE )
endif( NOT CMAKE_BUILD_TYPE )
set( CMAKE_CONFIGURATION_TYPES "${CMAKE_BUILD_TYPE}" CACHE INTERNAL "internal" )


# Project name and language:
project( libTheSky Fortran )


# Increase verbosity for debugging:
option( CMAKE_VERBOSE_MAKEFILE "Verbose makefile" on )

# Search in the CMake/ directory for CMake modules:
list( APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake )

# Various compile/optimisation options that we may want to enable:
include( SetCompileOptions )


# Place the products in their directories:
get_filename_component( Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME )
if( COMPILER_SPECIFIC_LIBS )
  set( MODULE_DIRECTORY "${CMAKE_SOURCE_DIR}/usr/include/libTheSky/${Fortran_COMPILER_NAME}" )
else( COMPILER_SPECIFIC_LIBS )
  set( MODULE_DIRECTORY "${CMAKE_SOURCE_DIR}/usr/include/libTheSky" )
endif( COMPILER_SPECIFIC_LIBS )

set( LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/usr/lib${LIB_SUFFIX}" )


# Find dependencies:
find_package( LibSUFR REQUIRED )
set( INCLUDE_FLAGS "-I${LibSUFR_INCLUDES}" )  # will be transferred to CompilerFlags


# Set source files:
include( FileList )

# Set FORTRAN compiler flags:
include( CompilerFlags_Fortran )


# Find the package version from the VERSION file and stick it into PKG_VERSION:
# The file VERSION is either in the root directory of the package, or in doc/:
set( VERSION_FILE ${CMAKE_SOURCE_DIR}/VERSION )
if( NOT EXISTS ${VERSION_FILE} )
  set( VERSION_FILE ${CMAKE_SOURCE_DIR}/doc/VERSION )
endif( NOT EXISTS ${VERSION_FILE} )

set( PKG_VERSION "0.0.0" )   # Indicates that the version number was not found
file(STRINGS ${VERSION_FILE} Lines )
foreach( Line ${Lines} )
  string(REGEX MATCH "Release version:.+" Line ${Line} )  # Returns the matching string
  string(COMPARE NOTEQUAL "${Line}" "" Matches )          # If the string is not empty, we have a match
  if( Matches )
    string(REPLACE "Release version:" "" Line ${Line})    # Remove text
    string(STRIP ${Line} Line)                            # Strip leading and trailing spaces
    set( PKG_VERSION ${Line} )
  endif( Matches )
endforeach()
message( STATUS "libTheSky version in the current package: ${PKG_VERSION}" )
message( STATUS "" )



# Create the file code_version.f90, which contains the code version number/hash and date:
if( NOT EXISTS ${CMAKE_SOURCE_DIR}/src/code_version.f90 OR CREATE_VERSION )
  add_custom_command(
    OUTPUT ${CMAKE_SOURCE_DIR}/src/code_version.f90
    COMMAND cd $(CMAKE_SOURCE_DIR)
    COMMAND . ${CMAKE_SOURCE_DIR}/code_version.sh ${CMAKE_SOURCE_DIR} src/code_version.f90 ${Fortran_COMPILER_NAME} ${OPT_FLAGS}
    )
endif( NOT EXISTS ${CMAKE_SOURCE_DIR}/src/code_version.f90 OR CREATE_VERSION )


# Create libraries:
if( CREATE_SHAREDLIB )
  add_library( "libTheSky_shared"  SHARED  ${SRC_FILES} )
  if( COMPILER_SPECIFIC_LIBS )
    set_target_properties( libTheSky_shared PROPERTIES OUTPUT_NAME "TheSky_${Fortran_COMPILER_NAME}" )
  else( COMPILER_SPECIFIC_LIBS )
    set_target_properties( libTheSky_shared PROPERTIES OUTPUT_NAME "TheSky" )
  endif( COMPILER_SPECIFIC_LIBS )
  set_target_properties( libTheSky_shared PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIRECTORY} )
endif( CREATE_SHAREDLIB )

# Use add_dependencies below in order to pretend the static lib depends on the shared.  This
#   seems to be necessary to ensure that it is built *after* shared.  If this doesn't
#   happen, multiple threads will be pouring data into the code_version file simultaneously.

if( CREATE_STATICLIB )
  add_library( "libTheSky_static"  STATIC  ${SRC_FILES} )
  if( COMPILER_SPECIFIC_LIBS )
    set_target_properties( libTheSky_static PROPERTIES OUTPUT_NAME "TheSky_${Fortran_COMPILER_NAME}" )
  else( COMPILER_SPECIFIC_LIBS )
    set_target_properties( libTheSky_static PROPERTIES OUTPUT_NAME "TheSky" )
  endif( COMPILER_SPECIFIC_LIBS )
  set_target_properties( libTheSky_static PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIRECTORY} )
  add_dependencies( libTheSky_static libTheSky_shared PKG_VERSION )  # Ensure it is built after the shared lib
endif( CREATE_STATICLIB )


add_custom_target( "libTheSky_data"
  mkdir -p ${CMAKE_SOURCE_DIR}/usr/share/libTheSky/
  COMMAND cp -vf ${CMAKE_SOURCE_DIR}/data/* ${CMAKE_SOURCE_DIR}/usr/share/libTheSky/ || echo "*** WARNING: No data files found - did you forget to download them?  This will not affect the creation of the libraries ***"
  COMMAND echo
  COMMENT "Installing data files"
  )


# Install to destination:
#   use DIRECTORY rather than TARGETS to allow include/<fortrancompiler>/
install( DIRECTORY usr/ DESTINATION ${CMAKE_INSTALL_PREFIX} )

# Install man pages:
install( DIRECTORY man/ DESTINATION ${CMAKE_INSTALL_PREFIX}/share/man )

# Install documents:
if( NOT EXISTS VERSION )
  install( FILES doc/CHANGELOG doc/LICENCE doc/README doc/VERSION DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/libthesky-${PKG_VERSION} )
else( NOT EXISTS VERSION )
  install( FILES CHANGELOG LICENCE README VERSION DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/libthesky-${PKG_VERSION} )
endif( NOT EXISTS VERSION )

