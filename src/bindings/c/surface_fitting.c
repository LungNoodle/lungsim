#include "surface_fitting.h"
#include "string.h"

void fit_surface_geometry_c(int *niterations,const char *fitting_file,int *fitting_file_len);
void reset_fitting_c();
void define_data_fit_group_c(const char *datafile, int *datafile_len, const char *groupname, int *groupname_len);
void initialise_fit_mesh_c(void);

void fit_surface_geometry(int niterations,const char *fitting_file)
{
  int filename_len = strlen(fitting_file);
  fit_surface_geometry_c(&niterations, fitting_file, &filename_len);
}

void reset_fitting()
{
  reset_fitting_c();
}

void define_data_fit_group(const char *datafile, const char *groupname)
{
  int datafile_len = strlen(datafile);
  int groupname_len = strlen(groupname);
  define_data_fit_group_c(datafile, &datafile_len, groupname, &groupname_len);
}

void initialise_fit_mesh()
{
  initialise_fit_mesh_c();
}

