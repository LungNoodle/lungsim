
#ifndef AETHER_DIAGNOSTICS_H
#define AETHER_DIAGNOSTICS_H

#include "symbol_export.h"

SHO_PUBLIC void enter_exit(const char *sub_name, int place);
SHO_PUBLIC void set_diagnostics_on(int state);
SHO_PUBLIC int get_diagnostics_on();


#endif /* AETHER_DIAGNOSTICS_H */
