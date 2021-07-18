
#include "ventilation.h"

void evaluate_vent_c(int *num_breaths, double *dt);
void evaluate_uniform_flow_c();
void two_unit_test_c();

void evaluate_vent(int num_breaths, double dt)
{
  evaluate_vent_c(&num_breaths, &dt);
}

void evaluate_uniform_flow()
{
  evaluate_uniform_flow_c();
}

void two_unit_test()
{
  two_unit_test_c();
}
