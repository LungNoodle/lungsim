#include "surface_fitting.h"
#include "string.h"

void fit_surface_geometry_c(int *niterations,const char *fitting_file,int *fitting_file_len);
void initialise_fit_mesh_c(void);

void fit_surface_geometry(int niterations,const char *fitting_file)
{
  int filename_len = strlen(fitting_file);
  fit_surface_geometry_c(&niterations, fitting_file, &filename_len);
}

void initialise_fit_mesh()
{
  initialise_fit_mesh_c();
}

