module pressure_resistance_flow_c
implicit none
private

contains

!!!###################################################################################

subroutine evaluate_prq_c() bind(C, name="evaluate_prq_c")

use pressure_resistance_flow, only: evaluate_prq
implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
call so_evaluate_prq
#else
call evaluate_prq
#endif

end subroutine evaluate_prq_c

!###################################################################################
end module pressure_resistance_flow_c