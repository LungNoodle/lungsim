
#include "arrays.h"

void set_node_field_value_c(int *row, int *col, double *value);

void set_node_field_value(int row, int col, double value)
{
  set_node_field_value_c(&row, &col, &value);
}
