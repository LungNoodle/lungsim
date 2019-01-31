#include "wave_transmission.h"

void evaluate_wave_transmission_c(int *n_time, double *heartrate, double *a0, int *no_freq, double *a, double *b, int *n_bcparams, double *admittance_param,int *n_model, double *model_definition);

void evaluate_wave_transmission(int n_time, double heartrate, double a0, int no_freq, double *a, double *b, int n_bcparams, double *admittance_param, int n_model, double *model_definition)
{
    evaluate_wave_transmission_c(&n_time, &heartrate, &a0, &no_freq, a, b, &n_bcparams, admittance_param, &n_model, model_definition);
}
