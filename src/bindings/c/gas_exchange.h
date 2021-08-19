#ifndef AETHER_GAS_EXCHANGE_H
#define AETHER_GAS_EXCHANGE_H

#include "symbol_export.h"

SHO_PUBLIC void initial_gasexchange(double initial_concentration, double surface_area, double V_cap);
SHO_PUBLIC void steadystate_gasexchange(double c_art_o2, double c_ven_o2, double p_art_co2, double p_art_o2, 
	double p_i_o2, double p_ven_co2, double p_ven_o2, double shunt_fraction, double VCO2, double VO2);

#endif /* AETHER_GAS_EXCHANGE_H */
