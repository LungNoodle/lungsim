
#include "growtree.h"

#include "string.h"

extern void grow_tree_c(int *elemlist_len, int elemlist[], int *parent_ne, double *angle_max, double *angle_min, double *branch_fraction, double *length_limit, double *shortest_length, double *rotation_limit, int *to_export, const char *filename, int *filename_len);

extern void smooth_1d_tree_c(int *num_elem_start, double *length_limit);

void grow_tree(int elemlist_len, int elemlist[], int parent_ne, double angle_max, double angle_min, double branch_fraction, double length_limit, double shortest_length, double rotation_limit, int to_export, const char *filename)
{
  int filename_len = strlen(filename);
  grow_tree_c(&elemlist_len, elemlist, &parent_ne, &angle_max, &angle_min, &branch_fraction, &length_limit, &shortest_length, &rotation_limit, &to_export, filename, &filename_len);
}

void smooth_1d_tree(int num_elem_start, double length_limit)
{
  smooth_1d_tree_c(&num_elem_start, &length_limit);
}
