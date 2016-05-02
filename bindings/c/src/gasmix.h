
#ifndef AETHER_GASMIX_H
#define AETHER_GASMIX_H

#include "symbol_export.h"

SHO_PUBLIC void initial_gasmix(double initial_concentration, double inlet_concentration);
SHO_PUBLIC void solve_gasmix(int fd, int inr_itr_max, int out_itr_max, double diffusion_coeff,
                             double dt, double initial_volume, double inlet_concentration, double inlet_flow, double solve_tolerance, double time_end,
                             double time_start,int inspiration);
SHO_PUBLIC void transfer_flow_vol_from_units();

#endif /* AETHER_GASMIX_H */
