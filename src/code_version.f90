!> \file code_version.f90  Source file generated by CMake to report the libTheSky version used

!***********************************************************************************************************************************
!> \brief  Report libTheSky version number

module TheSky_version
  implicit none
  save
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Subroutine generated by CMake to report the code version used
  
  subroutine print_libTheSky_version(unit)
    implicit none
    integer, intent(in) :: unit
    character :: libTheSky_version*(99) = 'rev.57, hash eff674d (2014-06-28)'
    character :: compile_date*(99) = 'Sat 28 Jun 18:39:18 CEST 2014'
    character :: compiler*(99) = 'gfortran'
    character :: compiler_flags*(99) = '-O2'
    
    write(unit,'(A)', advance='no') 'libTheSky '//trim(libTheSky_version)//', compiled on '//trim(compile_date)//' with '// &
         trim(compiler)//' '//trim(compiler_flags)
    
  end subroutine print_libTheSky_version
  !*********************************************************************************************************************************
  
end module TheSky_version
!***********************************************************************************************************************************

