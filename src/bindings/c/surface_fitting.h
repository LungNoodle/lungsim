#ifndef AETHER_SURFACE_FITTING_H
#define AETHER_SURFACE_FITTING_H

#include "symbol_export.h"

SHO_PUBLIC void fit_surface_geometry(int niterations, const char *fitting_file);
SHO_PUBLIC void initialise_fit_mesh();

#endif /* AETHER_SURFACE_FITTING_H */
