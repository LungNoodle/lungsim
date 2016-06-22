#include "tree_wave_propagation.h"

void evaluate_wave_propagation_c(double *a0, int *no_freq, double *a, double *b);

void evaluate_wave_propagation(double a0, int no_freq, double *a, double *b)
{
    evaluate_wave_propagation_c(&a0, &no_freq, a, b);
}

