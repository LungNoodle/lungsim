#ifndef AETHER_WAVE_TRANSMISSION_H
#define AETHER_WAVE_TRANSMISSION_H

#include "symbol_export.h"

SHO_PUBLIC void evaluate_wave_transmission(int grav_dirn, double grav_factor, int n_time, double heartrate, double a0, int no_freq, double *a, double *b, int n_adparams, double *admittance_param, int n_model, double *model_definition, int cap_model, double remodeling_grade, const char *bc_type, const char *lobe_imped);


#endif /* AETHER_WAVE_TRANSMISSION_H */
