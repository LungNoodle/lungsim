!> \file
!> \author Merryn Tawhai
!> \brief This module contains definition of constants
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>
!>\Description
!> This module contains definition of constants (note that in the future this should be merged into a 'types' module
module other_consts
  use arrays, only: dp
  implicit none

  integer, parameter :: MAX_FILENAME_LEN = 255, MAX_STRING_LEN = 100, MAX_SUBNAME_LEN = 60

  real(dp), parameter :: PI = 3.14159265358979_dp
  real(dp), parameter :: TOLERANCE =EPSILON (EPSILON(1.0_dp))

  private
  public MAX_SUBNAME_LEN, MAX_STRING_LEN, MAX_FILENAME_LEN, PI,TOLERANCE
end module other_consts
