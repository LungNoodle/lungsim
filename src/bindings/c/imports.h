#ifndef AETHER_IMPORTS_H
#define AETHER_IMPORTS_H

#include "symbol_export.h"

SHO_PUBLIC void import_capillary(const char *FLOWFILE);
SHO_PUBLIC void import_terminal(const char *FLOWFILE);
SHO_PUBLIC void import_ventilation(const char *FLOWFILE);
SHO_PUBLIC void import_perfusion(const char *FLOWFILE);

#endif /* AETHER_IMPORTS_H */
