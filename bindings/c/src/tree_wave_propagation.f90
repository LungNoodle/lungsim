module tree_wave_propagation_c
implicit none
private

contains

!!!###################################################################################

subroutine evaluate_wave_propagation_c() bind(C, name="evaluate_wave_propagation_c")

use tree_wave_propagation, only: evaluate_wave_propagation
implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
call so_evaluate_wave_propagation
#else
call evaluate_wave_propagation
#endif

end subroutine evaluate_wave_propagation_c

!###################################################################################
end module tree_wave_propagation_c
