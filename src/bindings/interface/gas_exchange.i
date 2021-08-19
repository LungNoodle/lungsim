
%module(package="aether") gas_exchange
%include symbol_export.h
%include gas_exchange.h

%{
#include "gas_exchange.h"
%}

void initial_gasexchange(double initial_concentration, double surface_area, double V_cap);
void steadystate_gasexchange(double c_art_o2, double c_ven_o2, double p_art_co2, double p_art_o2, 
	double p_i_o2, double p_ven_co2, double p_ven_o2, double shunt_fraction, double VCO2, double VO2);

%include gas_exchange.h

