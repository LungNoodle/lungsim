#include "pressure_resistance_flow.h"

void evaluate_prq_c(const char *mesh_type, int *mesh_type_len, int *grav_dirn, double *grav_factor, const char *bc_type, int *bc_type_len,  double *inlet_bc, double *outlet_bc);

void evaluate_prq(const char *mesh_type, int grav_dirn, double grav_factor, const char *bc_type, double inlet_bc, double outlet_bc)
{
    int mesh_type_len = strlen(mesh_type);
	int bc_type_len = strlen(bc_type);
    evaluate_prq_c(mesh_type, &mesh_type_len, &grav_dirn, &grav_factor, bc_type, &bc_type_len, &inlet_bc, &outlet_bc);
}
