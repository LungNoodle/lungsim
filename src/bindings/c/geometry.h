
#ifndef AETHER_GEOMETRY_H
#define AETHER_GEOMETRY_H

#include "symbol_export.h"

SHO_PUBLIC void add_mesh(const char *AIRWAY_MESHFILE);
SHO_PUBLIC void add_matching_mesh();
SHO_PUBLIC void append_units();
SHO_PUBLIC void define_1d_elements(const char *ELEMFILE);
SHO_PUBLIC void define_elem_geometry_2d(const char *ELEMFILE, const char *SF_OPTION);
SHO_PUBLIC void define_mesh_geometry_test();
SHO_PUBLIC void define_node_geometry(const char *NODEFILE);
SHO_PUBLIC void define_node_geometry_2d(const char *NODEFILE);
SHO_PUBLIC void define_data_geometry(const char *DATAFILE);
SHO_PUBLIC void group_elem_parent_term(int ne_parent);
SHO_PUBLIC void make_data_grid(int surface_elems, double spacing, int to_export, const char *filename, const char *groupname);
SHO_PUBLIC void make_2d_vessel_from_1d(int elemlist_len, int elemlist[]);
SHO_PUBLIC void define_rad_from_file(const char *FIELDFILE, const char *radius_type);
SHO_PUBLIC int get_local_node_f(const char *ndimenstion, const char *np_global);
SHO_PUBLIC void define_rad_from_geom(const char *ORDER_SYSTEM, double CONTROL_PARAM, const char *START_FROM,
                                     double START_RAD, const char *GROUP_TYPE, const char *GROUP_OPTIONS);
SHO_PUBLIC void element_connectivity_1d();
SHO_PUBLIC void evaluate_ordering();
SHO_PUBLIC void set_initial_volume(int Gdirn, double COV, double total_volume, double Rmax, double Rmin);
SHO_PUBLIC void volume_of_mesh(double *volume_model, double *volume_tree);
SHO_PUBLIC void write_elem_geometry_2d(const char *ELEMFILE);
SHO_PUBLIC void write_geo_file(int ntype, const char *GEOFILE);
SHO_PUBLIC void write_node_geometry_2d(const char *NODEFILE);

#endif /* AETHER_GEOMETRY_H */
