#ifndef LEPTONIC_H
#define LEPTONIC_H

#include "Constants.h"

typedef double (*prefactor)(double Q2);

double FX2_EM(double Q2);
double FX2_CC(double Q2);

//Dealing with the lepton side of the interaction
//Compute lepton tensor and lepton scattering kinematic variables
struct LeptonTensor
{
	//Inputs
	double mi = electronmass;
	double mf = electronmass;
	double helicity=0.; //chirality
	prefactor FX2 = FX2_EM;

	//outputs
	double FF;
	double L_s[4][4]; //Symmetric part of lepton tensor
	double L_a[4][4]; //Antisymetric part
	double w, q, Q2, ki, kf;

	//Set full lepton tensor and all kinematics
	int eval(double Ei, double Ef, double costheta_f);

	//Only evaluate upper diagonal non-zero terms of lepton tensor
	int eval_minim(double Ei, double Ef, double costheta_f);

	//Only evaluate the kinematic variables
	int eval_kinematics(double Ei, double Ef, double costheta_f);
	void set_zero();
	void set_elec();
	void set_CC_numu();
};





#endif
