
#include "gas_exchange.h"

void initial_gasexchange_c(double *initial_concentration, double *surface_area, double *V_cap);

void steadystate_gasexchange_c(double *deadspace, double *p_i_o2, double *shunt_fraction, double *target_p_art_co2, double *target_p_ven_o2, double *VCO2, double *VO2);


void initial_gasexchange(double initial_concentration, double surface_area, double V_cap)
{
  initial_gasexchange_c(&initial_concentration, &surface_area, &V_cap);
}

void steadystate_gasexchange(double deadspace, double p_i_o2, double shunt_fraction, double target_p_art_co2, double target_p_ven_o2, double VCO2, double VO2)
{
  steadystate_gasexchange_c(&deadspace, &p_i_o2, &shunt_fraction, &target_p_art_co2, &target_p_ven_o2, &VCO2, &VO2);
}

