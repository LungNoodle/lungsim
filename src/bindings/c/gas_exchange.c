
#include "gas_exchange.h"

void steadystate_gasexchange_c(double *c_art_o2, double *c_ven_o2,
	double *p_art_co2, double *p_art_o2, double *p_i_o2, double *p_ven_co2, double *p_ven_o2, double *shunt_fraction,
	double *VCO2, double *VO2);

void steadystate_gasexchange(double *c_art_o2, double *c_ven_o2,
	double *p_art_co2, double *p_art_o2, double *p_i_o2, double *p_ven_co2, double *p_ven_o2, double *shunt_fraction,
	double *VCO2, double *VO2)
{
	steadystate_gasexchange_c(c_art_o2, c_ven_o2, p_art_co2, p_art_o2, p_i_o2, p_ven_co2, p_ven_o2, shunt_fraction, VCO2, VO2);
}

