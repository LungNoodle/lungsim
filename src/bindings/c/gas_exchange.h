#ifndef AETHER_GAS_EXCHANGE_H
#define AETHER_GAS_EXCHANGE_H

#include "symbol_export.h"

SHO_PUBLIC void initial_gasexchange(double initial_concentration, double surface_area, double V_cap);
SHO_PUBLIC void steadystate_gasexchange(double deadspace, double p_i_o2, double shunt_fraction, double target_p_art_co2, double target_p_ven_o2, double VCO2, double VO2);

#endif /* AETHER_GAS_EXCHANGE_H */
