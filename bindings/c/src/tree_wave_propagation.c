#include "tree_wave_propagation.h"

void evaluate_wave_propagation_c(int *no_freq,double *a0,double *a,double *b);

void evaluate_wave_propagation(int no_freq,double a0,double a,double b)
{
    evaluate_wave_propagation_c(&no_freq,&a0,&a,&b);
}

