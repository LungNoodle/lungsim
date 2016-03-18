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
  use arrays,only: dp
  implicit none

  real(dp), parameter :: PI = 3.14159265358979_dp
  logical :: diagnostics_on
  logical :: test_cpu_time

end module other_consts
