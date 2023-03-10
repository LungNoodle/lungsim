
#include "geometry.h"

#include "string.h"

void add_mesh_c(const char *AIRWAY_MESHFILE, int *filename_len);
void add_matching_mesh_c(void);
void append_units_c(void);
void define_1d_elements_c(const char *ELEMFILE, int *filename_len);
void define_elem_geometry_2d_c(const char *ELEMFILE, int *filename_len, const char *SF_OPTION, int *sf_option_len);
void define_mesh_geometry_test_c(void);
void define_node_geometry_c(const char *NODEFILE, int *filename_len);
void define_node_geometry_2d_c(const char *NODEFILE, int *filename_len);
void define_data_geometry_c(const char *DATAFILE, int *filename_len);
void import_node_geometry_2d_c(const char *NODEFILE, int *filename_len);
void import_ply_triangles_c(const char *ply_file, int *filename_len);
extern void make_data_grid_c(int *elemlist_len, int elemlist[], double *offset, double *spacing, const char *filename, int *filename_len, const char *groupname, int *groupname_len);
extern void make_2d_vessel_from_1d_c(int *elemlist_len, int elemlist[]);
void define_rad_from_file_c(const char *FIELDFILE, int *filename_len, const char *radius_type, int *radius_type_len);
int get_local_node_f_c(const char *ndimension, int *dimension_len, const char *np_global, int *np_global_len);
void define_rad_from_geom_c(const char *order_system, int *order_system_len, double *control_param,
                            const char *start_from, int *start_from_len, double *start_rad,
                            const char *group_type, int *group_type_len, const char *group_options, int *group_options_len);
void element_connectivity_1d_c(void);
void evaluate_ordering_c(void);
void volume_of_mesh_c(double *volume_model, double *volume_tree);
void write_elem_geometry_2d_c(const char *ELEMFILE, int *filename_len);
void write_geo_file_c(int *ntype, const char *GEOFILE, int *filename_len);
void write_node_geometry_2d_c(const char *NODEFILE, int *filename_len);

void add_mesh(const char *AIRWAY_MESHFILE)
{
  int filename_len = (int)strlen(AIRWAY_MESHFILE);
  add_mesh_c(AIRWAY_MESHFILE, &filename_len);
}

void add_matching_mesh()
{
	add_matching_mesh_c();
}

void append_units()
{
  append_units_c();
}

void define_1d_elements(const char *ELEMFILE)
{
  int filename_len = (int)strlen(ELEMFILE);
  define_1d_elements_c(ELEMFILE, &filename_len);
}

void define_elem_geometry_2d(const char *ELEMFILE, const char *SF_OPTION)
{
  int filename_len = (int)strlen(ELEMFILE);
  int sf_option_len = (int)strlen(SF_OPTION);
  define_elem_geometry_2d_c(ELEMFILE, &filename_len, SF_OPTION, &sf_option_len);
}

void define_mesh_geometry_test()
{
  define_mesh_geometry_test_c();
}

void define_node_geometry(const char *NODEFILE)
{
  int filename_len = (int)strlen(NODEFILE);
  define_node_geometry_c(NODEFILE, &filename_len);
}

void define_node_geometry_2d(const char *NODEFILE)
{
  int filename_len = (int)strlen(NODEFILE);
  define_node_geometry_2d_c(NODEFILE, &filename_len);
}

void define_data_geometry(const char *DATAFILE)
{
  int filename_len = (int)strlen(DATAFILE);
  define_data_geometry_c(DATAFILE, &filename_len);
}

void import_node_geometry_2d(const char *NODEFILE)
{
  int filename_len = (int)strlen(NODEFILE);
  import_node_geometry_2d_c(NODEFILE, &filename_len);
}

void import_ply_triangles(const char *ply_file)
{
  int filename_len = (int)strlen(ply_file);
  import_ply_triangles_c(ply_file, &filename_len);
}

void make_data_grid(int elemlist_len, int elemlist[], double offset, double spacing, const char *filename, const char *groupname)
{
  int filename_len = (int)strlen(filename);
  int groupname_len = (int)strlen(groupname);
  make_data_grid_c(&elemlist_len, elemlist, &offset, &spacing, filename, &filename_len, groupname, &groupname_len);
}

void make_2d_vessel_from_1d(int elemlist_len, int elemlist[])
{
  make_2d_vessel_from_1d_c(elemlist, &elemlist_len);
}

void define_rad_from_file(const char *FIELDFILE, const char *radius_type)
{
  int filename_len = (int)strlen(FIELDFILE);
  int radius_type_len = (int)strlen(radius_type);
  define_rad_from_file_c(FIELDFILE, &filename_len, radius_type, &radius_type_len);
}

int get_local_node_f(const char *ndimension, const char *np_global)
{
  int dimension_len = (int)strlen(ndimension);
  int np_global_len = (int)strlen(np_global);
  return get_local_node_f_c(ndimension, &dimension_len, np_global, &np_global_len);
}

void define_rad_from_geom(const char *order_system, double control_param, const char *start_from,
                          double start_rad, const char*group_type, const char *group_options)
{
  int order_system_len = (int)strlen(order_system);
  int start_from_len = (int)strlen(start_from);
  int group_type_len = (int)strlen(group_type);
  int group_options_len = (int)strlen(group_options);
  define_rad_from_geom_c(order_system, &order_system_len, &control_param, start_from, &start_from_len, &start_rad,
                         group_type, &group_type_len, group_options, &group_options_len);

}

void element_connectivity_1d()
{
  element_connectivity_1d_c();
}

void evaluate_ordering()
{
  evaluate_ordering_c();
}

void volume_of_mesh(double *volume_model, double *volume_tree)
{
  volume_of_mesh_c(volume_model, volume_tree);
}

void write_elem_geometry_2d(const char *ELEMFILE)
{
  int filename_len = (int)strlen(ELEMFILE);
  write_elem_geometry_2d_c(ELEMFILE, &filename_len);
}

void write_geo_file(int ntype, const char *GEOFILE)
{
  int filename_len = (int)strlen(GEOFILE);
  write_geo_file_c(&ntype, GEOFILE, &filename_len);
}

void write_node_geometry_2d(const char *NODEFILE)
{
  int filename_len = (int)strlen(NODEFILE);
  write_node_geometry_2d_c(NODEFILE, &filename_len);
}

