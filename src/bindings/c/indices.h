#ifndef AETHER_INDICES_H
#define AETHER_INDICES_H

#include "symbol_export.h"
SHO_PUBLIC void define_problem_type(const char *PROBLEMTYPE);
SHO_PUBLIC void ventilation_indices();
SHO_PUBLIC void perfusion_indices();
SHO_PUBLIC int get_ne_radius();
SHO_PUBLIC int get_nj_conc1();

#endif /* AETHER_INDICES_H */
