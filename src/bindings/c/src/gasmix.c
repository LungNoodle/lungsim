
#include "gasmix.h"

void initial_gasmix_c(double *initial_concentration, double *inlet_concentration);
void solve_gasmix_c(int *fd, int *inr_itr_max, int *out_itr_max, double *diffusion_coeff, double *dt,
                      double *initial_volume, double *inlet_concentration, double *inlet_flow,
                      double *solve_tolerance, double *time_end, double *time_start,int *inspiration);
void transfer_flow_vol_from_units_c();

void initial_gasmix(double initial_concentration, double inlet_concentration)
{
  initial_gasmix_c(&initial_concentration, &inlet_concentration);
}

void solve_gasmix(int fd, int inr_itr_max, int out_itr_max, double diffusion_coeff, double dt, double initial_volume,
                  double inlet_concentration, double inlet_flow, double solve_tolerance, double time_end,
                  double time_start,int inspiration)
{
  solve_gasmix_c(&fd, &inr_itr_max, &out_itr_max, &diffusion_coeff, &dt, &initial_volume, &inlet_concentration, &inlet_flow,
                   &solve_tolerance, &time_end, &time_start, &inspiration);
}

void transfer_flow_vol_from_units()
{
  transfer_flow_vol_from_units_c();
}
