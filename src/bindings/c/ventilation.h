#ifndef AETHER_VENTILATION_H
#define AETHER_VENTILATION_H

#include "symbol_export.h"

SHO_PUBLIC void evaluate_vent(const char *flow_filename, const char *pressure_filename);
SHO_PUBLIC void evaluate_uniform_flow();

#endif /* AETHER_VENTILATION_H */
