
#ifndef AETHER_GROWTREE_H
#define AETHER_GROWTREE_H

#include "symbol_export.h"

SHO_PUBLIC void grow_tree(int surface_elems, int parent_ne, double angle_max, double angle_min, double branch_fraction, double length_limit, double shortest_length, double rotation_limit);
SHO_PUBLIC void grow_tree_wrap(int elemlist_len, int elemlist[], int parent_ne, double angle_max, double angle_min, double branch_fraction, double length_limit, double shortest_length, double rotation_limit);

#endif /* AETHER_GROWTREE_H */
