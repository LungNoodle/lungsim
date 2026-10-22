
#include "ventilation.h"
#include <string.h>

void evaluate_vent_c(const char *flow_filename, int *flow_filename_len,
                     const char *pressure_filename, int *pressure_filename_len);
void evaluate_uniform_flow_c();

void evaluate_vent(const char *flow_filename, const char *pressure_filename)
{
  int flow_filename_len = strlen(flow_filename);
  int pressure_filename_len = strlen(pressure_filename);
  evaluate_vent_c(flow_filename, &flow_filename_len,
                  pressure_filename, &pressure_filename_len);
}

void evaluate_uniform_flow()
{
  evaluate_uniform_flow_c();
}

