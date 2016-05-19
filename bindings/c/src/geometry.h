
#ifndef AETHER_GEOMETRY_H
#define AETHER_GEOMETRY_H

#include "symbol_export.h"

SHO_PUBLIC void add_mesh(const char *AIRWAY_MESHFILE);
SHO_PUBLIC void append_units();
SHO_PUBLIC void define_1d_elements(const char *ELEMFILE);
SHO_PUBLIC void define_mesh_geometry_test();
SHO_PUBLIC void define_node_geometry(const char *NODEFILE);
SHO_PUBLIC void define_rad_from_file(const char *FIELDFILE, const char *radius_type);
SHO_PUBLIC void define_rad_from_geom(const char *ORDER_SYSTEM, double CONTROL_PARAM, const char *START_FROM,
                                     double START_RAD, const char *GROUP_TYPE, const char *GROUP_OPTIONS);
SHO_PUBLIC void element_connectivity_1d();
SHO_PUBLIC void evaluate_ordering();
SHO_PUBLIC void set_initial_volume(int Gdirn, double COV, double total_volume, double Rmax, double Rmin);
SHO_PUBLIC void volume_of_mesh(double *volume_model, double *volume_tree);

#endif /* AETHER_GEOMETRY_H */
