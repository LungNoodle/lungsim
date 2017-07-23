module tree_wave_propagation

!*Brief Description:* Simulating wave propagation in a 1D tree structure
!
!*LICENSE:*
!
!
!
!*Full Description:*
!Simulating wave propagation in a 1D tree structure
!
  use arrays, only: dp
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public evaluate_wave_propagation


contains
!
!##############################################################################
!
subroutine evaluate_wave_propagation(a0,no_freq,a,b)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SOLVE_WAVE_PROPAGATION: SOLVE_WAVE_PROPAGATION
  use diagnostics, only: enter_exit

  integer, intent(in):: no_freq
  real(dp), intent(in) :: a0
  real(dp), intent(in) :: a(no_freq)
  real(dp), intent(in) :: b(no_freq)
  !integer, intent(out) :: you
  !real(dp), intent(inout)  :: pass
  !local variables
  integer :: other,stuff

  character(len=60) :: sub_name

  sub_name = 'evalulate_wave_propagation'
  call enter_exit(sub_name,1)

  write(*,*) 'No freqs',a



      call enter_exit(sub_name,2)
end subroutine evaluate_wave_propagation

end module tree_wave_propagation
