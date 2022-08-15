
#include "arrays.h"
#include "string.h"

void set_node_field_value_c(int *row, int *col, double *value);
void update_parameter_c(const char *parameter_name, int *parameter_name_len, double *parameter_value);

void set_node_field_value(int row, int col, double value)
{
  set_node_field_value_c(&row, &col, &value);
}

void update_parameter(const char *parameter_name, double parameter_value)
{
  int parameter_name_len = (int)strlen(parameter_name);
  update_parameter_c(parameter_name, &parameter_name_len, &parameter_value);
}

