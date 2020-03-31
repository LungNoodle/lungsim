module precision
  !
  !*Brief Description:* Defines the precision of real parameters.
  !
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  ! Defines the precision of real, via a parameter 'dp'.
  implicit none
  
  !Module parameters
  integer, parameter :: dp=kind(0.d0) ! for double precision  

  real(dp),parameter :: zero_tol = 1.0e-12_dp
  real(dp),parameter :: loose_tol = 1.0e-6_dp

  !Module types

  !Module variables

  !Interfaces
  private
  public dp
  public zero_tol
  public loose_tol
  
end module precision
