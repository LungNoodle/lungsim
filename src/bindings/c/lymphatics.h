#ifndef AETHER_LYMPHATICS_H
#define AETHER_LYMPHATICS_H

#include "symbol_export.h"

SHO_PUBLIC void alveolar_capillary_flux(int num_nodes, int write_out);
SHO_PUBLIC void lymphatic_transport(const char *filename);

#endif /* AETHER_LYMPHATICS_H */
