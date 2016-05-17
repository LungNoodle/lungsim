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
  use arrays, only: dp
  implicit none

  integer, parameter :: MAX_FILENAME_LEN = 255, MAX_STRING_LEN = 100, MAX_SUBNAME_LEN = 60
  integer :: num_symm_gen

  real(dp), parameter :: PI = 3.14159265358979_dp
  real(dp) :: Palv,sheet_h0,K_cap,F_cap,F_sheet,sigma_cap,mu_c,alpha_c,area_scale
  real(dp) :: F_rec,sigma_rec,L_c,Plb_c,Pub_c,Pub_a_v,total_cap_area,L_art_terminal
  real(dp) :: H0,L_scale,L_vein_terminal,Ptm_max,R_art_terminal,R_vein_terminal
  real(dp) :: INITIAL_LPM,mu_app(13)
  real(dp) :: L_a(13),L_v(13),rad_a(13),rad_v(13),alpha_a,alpha_v

  private
  public MAX_SUBNAME_LEN, MAX_STRING_LEN, MAX_FILENAME_LEN, PI
end module other_consts
