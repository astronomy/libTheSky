##  CompilerFlags_Fortran.cmake
##  Compiler flags for Fortran compilers
##  Currently, specific flags for gfortran, g95 and ifort are provided
##  Make sure to choose the correct ~last line (printing of the compiler options)
##  
##  Copyright 2010-2024 Marc van der Sluys - marc.vandersluys.nl
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


# Get compiler name:
get_filename_component( Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME )

# Are we on Linux?
if( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
  set( LINUX TRUE )
endif( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )


######################################################################################################################################################
#  Specific options per compiler:
######################################################################################################################################################
if( Fortran_COMPILER_NAME MATCHES "gfortran" )
  
  # Get compiler version:
  exec_program(
    ${CMAKE_Fortran_COMPILER}
    ARGS --version
    OUTPUT_VARIABLE _compiler_output)
  string(REGEX REPLACE ".* ([0-9]\\.[0-9]\\.[0-9]).*" "\\1"
    COMPILER_VERSION ${_compiler_output})
  
  
  set( CMAKE_Fortran_FLAGS_ALL "-std=f2008 -fall-intrinsics -pedantic" )
  set( CMAKE_Fortran_FLAGS_ALL "${CMAKE_Fortran_FLAGS_ALL} -frecursive" )  # don't complain about Array larger than -fmax-stack-var-size
  if( COMPILER_VERSION VERSION_GREATER "4.4.99" )
    set( CMAKE_Fortran_FLAGS_ALL "${CMAKE_Fortran_FLAGS_ALL} -fwhole-file" )  # >= v.4.5
  endif( COMPILER_VERSION VERSION_GREATER "4.4.99" )
  
  set( CMAKE_Fortran_FLAGS "-pipe -funroll-loops" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-pipe -funroll-loops" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-g -ffpe-trap=zero,invalid -fsignaling-nans -fbacktrace" )
  set( CMAKE_Fortran_FLAGS_PROFILE "-g -gp" )
  
  
  if(WANT_32BIT )
    set( BIT_FLAGS "-m32" )
  endif(WANT_32BIT )
  if(WANT_64BIT )
    set( BIT_FLAGS "-m64" )
  endif(WANT_64BIT )
  
  if(WANT_SSE42 )
    set( SSE_FLAGS "-msse4.2" )
  endif(WANT_SSE42 )
  
  if( WANT_OPENMP )
    set( OPENMP_FLAGS "-fopenmp" )
    message( STATUS "Compiling with OpenMP support" )
  endif( WANT_OPENMP )
  
  if( WANT_STATIC )
    set( STATIC_FLAGS "-static" )
  endif( WANT_STATIC )
  
  if( WANT_CHECKS )
    set( CHECK_FLAGS "-ffpe-trap=invalid,zero,overflow,underflow,denormal -fsignaling-nans -fbacktrace" )
    if( COMPILER_VERSION VERSION_GREATER "4.4.99" )
      set( CHECK_FLAGS "-fcheck=all ${CHECK_FLAGS}" )    # >= v.4.5
    else( COMPILER_VERSION VERSION_GREATER "4.4.99" )
      set( CHECK_FLAGS "-fbounds-check ${CHECK_FLAGS}" ) # <= v.4.4
    endif( COMPILER_VERSION VERSION_GREATER "4.4.99" )

    set( OPT_FLAGS "-O0" )
  else( WANT_CHECKS )
    set( OPT_FLAGS "-O2" )
    #set( OPT_FLAGS "-Ofast" )
  endif( WANT_CHECKS )
  
  if( WANT_WARNINGS )
    set( WARN_FLAGS "-Wall -Wextra -Wno-conversion -ffree-line-length-0" )
  endif( WANT_WARNINGS )
  if( STOP_ON_WARNING )
    set( WARN_FLAGS "${WARN_FLAGS} -Werror" )
  endif( STOP_ON_WARNING )
  
  if( WANT_LIBRARY )
    set( LIB_FLAGS "-fPIC -g" )
  endif( WANT_LIBRARY )
  
  
  # Package-specific flags:
  set( PACKAGE_FLAGS "" )
  
  
  ####################################################################################################################################################
elseif( Fortran_COMPILER_NAME MATCHES "g95" )
  
  # Get compiler version:
  exec_program(
    ${CMAKE_Fortran_COMPILER}
    ARGS --version
    OUTPUT_VARIABLE _compiler_output)
  string(REGEX REPLACE ".*g95 ([0-9]*\\.[0-9]*).*" "\\1"
    COMPILER_VERSION ${_compiler_output})
  
  
  set( CMAKE_Fortran_FLAGS "" )
  set( CMAKE_Fortran_FLAGS_RELEASE "" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-g" )
  
  if( WANT_CHECKS )
    set( CHECK_FLAGS "-fbounds-check -ftrace=full" )
    set( OPT_FLAGS "-O0" )
  else( WANT_CHECKS )
    set( CHECK_FLAGS "-fshort-circuit" )
    set( OPT_FLAGS "-O2" )
  endif( WANT_CHECKS )
  
  if( WANT_WARNINGS )
    # 102: module procedure not referenced,  136: module variable not used,  165: implicit interface
    # set( WARN_FLAGS "-Wall -Wextra -Wno=102,136,165" )
    set( WARN_FLAGS "-Wall -Wextra" )
    set( WARN_FLAGS "-std=f2003 -ffree-line-length-huge ${WARN_FLAGS}" )
  endif( WANT_WARNINGS )
  if( STOP_ON_WARNING )
    set( WARN_FLAGS "${WARN_FLAGS} -Werror" )
  endif( STOP_ON_WARNING )
  
  if( WANT_LIBRARY )
    set( LIB_FLAGS "-fPIC -g" )
  endif( WANT_LIBRARY )
  
  
  # Package-specific flags:
  set( PACKAGE_FLAGS "" )
  
  
  ####################################################################################################################################################
elseif( Fortran_COMPILER_NAME MATCHES "ifort" )
  
  # Get compiler version:
  exec_program(
    ${CMAKE_Fortran_COMPILER}
    ARGS --version
    OUTPUT_VARIABLE _compiler_output)
  string(REGEX REPLACE ".* ([0-9]*\\.[0-9]*\\.[0-9]*) .*" "\\1"
    COMPILER_VERSION ${_compiler_output})
  
  
  set( CMAKE_Fortran_FLAGS_ALL "-nogen-interfaces -heap-arrays 1024" )
  set( CMAKE_Fortran_FLAGS "-vec-guard-write -fpconstant -funroll-loops -align all -ip" )
  if( LINUX )
    set( CMAKE_Fortran_FLAGS_ALL "${CMAKE_Fortran_FLAGS_ALL} -mcmodel=large" )  # -mcmodel exists for Linux only...
  endif( LINUX )
  set( CMAKE_Fortran_FLAGS_RELEASE "-vec-guard-write -fpconstant -funroll-loops -align all -ip" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-g -traceback" )
  set( CMAKE_Fortran_FLAGS_PROFILE "-g -gp" )
  
  # !!! HOST_OPT overrules SSE42, as ifort would !!!
  if(WANT_SSE42 )
    set( SSE_FLAGS "-axSSE4.2,SSSE3" )
  endif(WANT_SSE42 )
  if(WANT_HOST_OPT)
    set (SSE_FLAGS "-xHost")
  endif(WANT_HOST_OPT)
  
  if(WANT_IPO )
    set( IPO_FLAGS "-ipo" )
  endif(WANT_IPO )
  
  if( WANT_OPENMP )
    set( OPENMP_FLAGS "-openmp -openmp-report0" )
    message( STATUS "Compiling with OpenMP support" )
  endif( WANT_OPENMP )
  
  if( WANT_STATIC )
    set( STATIC_FLAGS "-static" )
  endif( WANT_STATIC )
  
  if( WANT_CHECKS )
    set( CHECK_FLAGS "-ftrapuv -check all -check noarg_temp_created -check nostack -traceback" )  #  -check stack doesn't link in ifort 2013.3.174 14.0.3
    set( OPT_FLAGS "-O0" )
  else( WANT_CHECKS )
    set( OPT_FLAGS "-O2" )
  endif( WANT_CHECKS )
  
  if( WANT_WARNINGS )
    set( WARN_FLAGS "-warn all -std08 -diag-disable 6894,8290,8291,5268" )   # 8290,8291: format for F,ES: too many decimal places (for negative numbers), 5268: line longer than 132 characters
  endif( WANT_WARNINGS )
  
  if( STOP_ON_WARNING )
    # set( WARN_FLAGS "${WARN_FLAGS}" )  # No option found yet - remove 'unused' warning from CMake this way
  endif( STOP_ON_WARNING )
  
  if( WANT_LIBRARY )
    set( LIB_FLAGS "-fPIC -g" )  # -fPIC is also provided by CMake?
  endif( WANT_LIBRARY )
  
  
  # Package-specific flags:
  set( PACKAGE_FLAGS "" )
  
  
  ####################################################################################################################################################
else( Fortran_COMPILER_NAME MATCHES "gfortran" )
  
  
  message( "CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER} )
  message( "Fortran compiler: " ${Fortran_COMPILER_NAME} )
  message( "No optimized Fortran compiler flags are known, we just try -O2..." )
  set( CMAKE_Fortran_FLAGS "" )
  set( CMAKE_Fortran_FLAGS_RELEASE "" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-g" )
  set( OPT_FLAGS "-O2" )
  
  
  # Package-specific flags:
  set( PACKAGE_FLAGS "" )
  
  
endif( Fortran_COMPILER_NAME MATCHES "gfortran" )
######################################################################################################################################################





######################################################################################################################################################
#  Put everything together:
######################################################################################################################################################

set( USER_FLAGS "${OPT_FLAGS} ${LIB_FLAGS} ${CHECK_FLAGS} ${WARN_FLAGS} ${BIT_FLAGS} ${SSE_FLAGS} ${IPO_FLAGS} ${OPENMP_FLAGS} ${STATIC_FLAGS} ${INCLUDE_FLAGS} ${PACKAGE_FLAGS}" )

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_ALL} ${CMAKE_Fortran_FLAGS} ${USER_FLAGS}" )
set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_ALL} ${CMAKE_Fortran_FLAGS_RELEASE} ${USER_FLAGS}" )
set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_DEBUG}" )  # These are added to the RELEASE flags when releasing with debug info




######################################################################################################################################################
#  Report what's going on:
######################################################################################################################################################

message( STATUS "" )
message( STATUS "Using Fortran compiler:  ${Fortran_COMPILER_NAME} ${COMPILER_VERSION}  (${CMAKE_Fortran_COMPILER})" )

if( WANT_CHECKS )
  message( STATUS "Compiling with run-time checks:  ${CHECK_FLAGS}" )
endif( WANT_CHECKS )
if( WANT_WARNINGS )
  message( STATUS "Compiling with warnings:  ${WARN_FLAGS}" )
endif( WANT_WARNINGS )
if( WANT_LIBRARY )
  message( STATUS "Compiling with library options:  ${LIB_FLAGS}" )
endif( WANT_LIBRARY )
if( WANT_STATIC )
  message( STATUS "Linking statically:  ${STATIC_FLAGS}" )
endif( WANT_STATIC )


## Choose the one which is actually used: see CMAKE_BUILD_TYPE in CMakeLists.txt:
##   (Typically, RELEASE for executable code and RELEASE + RELWITHDEBINFO for libraries)
# message( STATUS "Compiler flags used:  ${CMAKE_Fortran_FLAGS_RELEASE}" )
message( STATUS "Compiler flags used:  ${CMAKE_Fortran_FLAGS_RELEASE} ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" )
message( STATUS "" )



