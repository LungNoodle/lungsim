#ifndef AETHER_PRESSURE_RESISTANCE_FLOW_H
#define AETHER_PRESSURE_RESISTANCE_FLOW_H

#include "symbol_export.h"

SHO_PUBLIC void occlusion_list(int elemlist_len, int elemlist[]);
SHO_PUBLIC void compliance_list(int elemlist2_len, int elemlist2[]);
SHO_PUBLIC void set_perfusion_flows(int element_ids_len, int element_ids[], int flow_values_len, double flow_values[]);
SHO_PUBLIC void clear_perfusion_flows(void);
SHO_PUBLIC void evaluate_prq(const char *mesh_type, const char *vessel_type, int grav_dirn, double grav_factor, const char *bc_type, double inlet_bc, double outlet_bc, double remodeling_grade);


#endif /* AETHER_PRESSURE_RESISTANCE_FLOW_H */
