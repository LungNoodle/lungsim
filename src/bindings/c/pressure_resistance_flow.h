#ifndef AETHER_PRESSURE_RESISTANCE_FLOW_H
#define AETHER_PRESSURE_RESISTANCE_FLOW_H

#include "symbol_export.h"


SHO_PUBLIC void evaluate_prq(const char *mesh_type, const char *vessel_type, int grav_dirn, double grav_factor, const char *bc_type, double inlet_bc, double outlet_bc, double remodeling_grade);


#endif /* AETHER_PRESSURE_RESISTANCE_FLOW_H */
