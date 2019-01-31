module wave_transmission_c
implicit none
private

contains

!!!###################################################################################

subroutine evaluate_wave_transmission_c(n_time,heartrate,a0,no_freq,&
  a,b,n_adparams,admittance_param,n_model,model_definition) bind(C, name="evaluate_wave_transmission_c")

use wave_transmission, only: evaluate_wave_transmission
use arrays,only:dp
implicit none
integer, intent(in) :: n_time
real(dp), intent(in) :: heartrate
integer,intent(in) :: no_freq
real(dp),intent(in) :: a0
real(dp),intent(in) :: a(no_freq),b(no_freq)
integer, intent(in) :: n_adparams
real(dp), intent(in) :: admittance_param(n_adparams)
integer, intent(in) :: n_model
real(dp), intent(in) :: model_definition(n_model)

#if defined _WIN32 && defined __INTEL_COMPILER
call so_evaluate_wave_transmission(n_time,heartrate,a0,no_freq,a,b,n_adparams,admittance_param,&
  n_model,model_definition)
#else
call evaluate_wave_transmission(n_time,heartrate,a0,no_freq,a,b,n_adparams,admittance_param,&
  n_model,model_definition)
#endif

end subroutine evaluate_wave_transmission_c

!###################################################################################
end module wave_transmission_c
