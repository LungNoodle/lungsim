module wave_transmission_c
implicit none
private

contains

!!!###################################################################################

subroutine evaluate_wave_transmission_c(grav_dirn,grav_factor,n_time,heartrate,&
  a0,no_freq,a,b,n_adparams,admittance_param,n_model,model_definition,cap_model,&
  remodeling_grade,bc_type,bc_type_len,lobe_imped,lobe_imped_len) bind(C, name="evaluate_wave_transmission_c")

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
  real(dp), intent(in) :: remodeling_grade
  type(c_ptr), value, intent(in) :: bc_type, lobe_imped
  integer,intent(in) :: bc_type_len, lobe_imped_len
  character(len=MAX_STRING_LEN) :: bc_type_f, lobe_imped_f

  call strncpy(bc_type_f, bc_type, bc_type_len)
  call strncpy(lobe_imped_f, lobe_imped, lobe_imped_len)



  call evaluate_wave_transmission(grav_dirn,grav_factor,n_time,heartrate,a0,no_freq,a,&
    b,n_adparams,admittance_param,n_model,model_definition,cap_model,remodeling_grade,&
    bc_type_f,lobe_imped_f)

end subroutine evaluate_wave_transmission_c

!###################################################################################
end module wave_transmission_c
