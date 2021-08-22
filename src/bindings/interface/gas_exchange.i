
%module(package="aether") gas_exchange
%include symbol_export.h
%include gas_exchange.h

%{
#include "gas_exchange.h"
%}

void initial_gasexchange(double initial_concentration, double surface_area, double V_cap);
void steadystate_gasexchange(double deadspace, double p_i_o2, double shunt_fraction, double target_p_art_co2, double target_p_ven_o2, double VCO2, double VO2);

%include gas_exchange.h

