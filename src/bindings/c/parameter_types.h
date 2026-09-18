
#ifndef AETHER_PARAMETER_TYPES_H
#define AETHER_PARAMETER_TYPES_H

#include "symbol_export.h"

SHO_PUBLIC void update_lung(const char *param_name, double param_value);
SHO_PUBLIC void update_mechs(const char *param_name, double param_value);
SHO_PUBLIC void update_gasexchange(const char *param_name, double param_value);
SHO_PUBLIC void update_ventilation(const char *param_name, double param_value);
SHO_PUBLIC void update_cardiac(const char *param_name, double param_value);
SHO_PUBLIC void update_solve_gx(const char *param_name, double param_value);
SHO_PUBLIC void update_solve_v(const char *param_name, double param_value);
SHO_PUBLIC void update_species(const char *param_name);
SHO_PUBLIC void update_lymphatics(const char *param_name, double param_value);

#endif /* AETHER_PARAMETER_TYPES_H */
