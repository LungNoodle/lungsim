
#ifndef AETHER_GROWTREE_H
#define AETHER_GROWTREE_H

#include "symbol_export.h"

SHO_PUBLIC void grow_tree(int parent_ne, int surface_elems, double angle_max, double angle_min, double branch_fraction, double length_limit, double shortest_length, double rotation_limit);

#endif /* AETHER_GROWTREE_H */
