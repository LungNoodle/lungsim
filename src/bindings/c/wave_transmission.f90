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

subroutine export_time_dependent_perfusion_c(flow_filename,flow_filename_len, &
  pressure_filename,pressure_filename_len) bind(C, name="export_time_dependent_perfusion_c")

  use iso_c_binding, only: c_ptr
  use utils_c, only: strncpy
  use other_consts, only: MAX_FILENAME_LEN
  use wave_transmission, only: export_time_dependent_perfusion
  implicit none

  type(c_ptr),value,intent(in) :: flow_filename,pressure_filename
  integer,intent(in) :: flow_filename_len,pressure_filename_len
  character(len=MAX_FILENAME_LEN) :: flow_filename_f,pressure_filename_f

  call strncpy(flow_filename_f,flow_filename,flow_filename_len)
  call strncpy(pressure_filename_f,pressure_filename,pressure_filename_len)
  call export_time_dependent_perfusion(flow_filename_f,pressure_filename_f)

end subroutine export_time_dependent_perfusion_c

!###################################################################################
end module wave_transmission_c
