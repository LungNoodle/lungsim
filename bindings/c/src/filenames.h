#ifndef AETHER_FILENAMES_H
#define AETHER_FILENAMES_H

#include "symbol_export.h"

SHO_PUBLIC void read_geometry_main();
SHO_PUBLIC void read_geometry_evaluate_flow();
SHO_PUBLIC char *get_filename(const char *label);

#endif /* AETHER_FILENAMES_H */
