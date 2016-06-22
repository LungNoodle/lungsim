module tree_wave_propagation_c
implicit none
private

contains

!!!###################################################################################

subroutine evaluate_wave_propagation_c(a0,no_freq,a,b) bind(C, name="evaluate_wave_propagation_c")

use tree_wave_propagation, only: evaluate_wave_propagation
use arrays,only:dp
implicit none
integer,intent(in) :: no_freq
real(dp),intent(in) :: a0
real(dp),intent(in) :: a(no_freq),b(no_freq)

#if defined _WIN32 && defined __INTEL_COMPILER
call so_evaluate_wave_propagation(a0,no_freq,a,b)
#else
call evaluate_wave_propagation(a0,no_freq,a,b)
#endif

end subroutine evaluate_wave_propagation_c

!###################################################################################
end module tree_wave_propagation_c
