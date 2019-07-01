module wave_transmission_c
implicit none
private

contains

!!!###################################################################################

subroutine evaluate_wave_transmission_c(grav_dirn,grav_factor,n_time,heartrate,&
  a0,no_freq,a,b,n_adparams,admittance_param,n_model,model_definition,cap_model) bind(C, &
name="evaluate_wave_transmission_c")
use iso_c_binding, only: c_ptr
use utils_c, only: strncpy
use other_consts, only: MAX_STRING_LEN
use wave_transmission, only: evaluate_wave_transmission
use arrays,only:dp
implicit none

integer,intent(in) :: grav_dirn
real(dp),intent(in) :: grav_factor
integer, intent(in) :: n_time
real(dp), intent(in) :: heartrate
integer,intent(in) :: no_freq
real(dp),intent(in) :: a0
real(dp),intent(in) :: a(no_freq),b(no_freq)
integer, intent(in) :: n_adparams
real(dp), intent(in) :: admittance_param(n_adparams)
integer, intent(in) :: n_model
real(dp), intent(in) :: model_definition(n_model)
integer, intent(in) :: cap_model



#if defined _WIN32 && defined __INTEL_COMPILER
call so_evaluate_wave_transmission(grav_dirn,grav_factor,n_time,heartrate,a0,no_freq,a,&
  b,n_adparams,admittance_param,n_model,model_definition,cap_model)
#else
call evaluate_wave_transmission(grav_dirn,grav_factor,n_time,heartrate,a0,no_freq,a,&
  b,n_adparams,admittance_param,n_model,model_definition,cap_model)
#endif

end subroutine evaluate_wave_transmission_c

!###################################################################################
end module wave_transmission_c
