#ifndef CROSSSECTIONS_H
#define CROSSSECTIONS_H

#include "Formfactors.h"
#include "S_tables.h"
#include "Leptonic.h"

//Evaluates the differential cross section ds/dE_d\Omega_l d\Omega_N
//uses as input missing momentum pm instead of nucleon scattering angle
double cross_section(double E_i, double E_f, double costheta_f, double pm, double phi_N, LeptonTensor *Lep, FormFactors *FF, S_level *S);

//Gives the functions ABCDE to factorize the angular dependence
int Get_ABCDE(double E_i, double E_f, double costheta_f, double pm, LeptonTensor *Lep, FormFactors *FF, S_level *S, double ABCDE[5]); 

//Evaluates the differential cross sections ds/dE_d\Omega_l d\Omega_N
//uses as input \cos\theta_N, the scattering angle of the nucleon with respect to momentum transfer q.
//For any scattering angle there may be two physical final-states (with different nucleon momentum)
//The function returns the number of solutions and stores the values of the CS and nucleon momenta for each in the arrays CS and k_N_sols
int cross_section_cos(double E_i, double E_f, double costheta_f, double costhetaN, double phi_N, LeptonTensor *Lep, FormFactors *FF, S_level *S, double CS[2], double k_N_sols[2]);


//returns the inclusive cross section ds/dE_ld\Omega_l
// The partial cross sections \sigma_i, with i = {CC,CL,LL,T,Tp} are stored in the array CS_i
double CS_inclusive(double E_i, double E_f, double costheta_f, LeptonTensor *Lep, FormFactors *FF, S_level *S, double CS_i[5]);

#endif
