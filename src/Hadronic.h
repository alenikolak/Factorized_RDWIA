#ifndef HADRONIC_H
#define HADRONIC_H

#include <complex>

#include "Formfactors.h"
#include "S_tables.h"

int Get_Hadron(double TN, double pm, double cos, S_level *S, FormFactors *FF, complex<double> Hadron[4][4]);

int cosdep_Responses(double w, double q, S_level *S, FormFactors *FF);

int Inclusive_Responses(double w, double q, S_level *S, FormFactors *FF, double Responses[5]);

int Inclusive_Responses_th(double w, double q, S_level *S, FormFactors *FF, double Responses[5]);

int Get_Hadron_q(double TN, double q, double cos, S_level *S, FormFactors *FF, complex<double> Hadron[4][4]); 

#endif
