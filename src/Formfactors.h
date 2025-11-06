#ifndef FORMFACTORS_H
#define FORMFACTORS_H


typedef void (*formfactors)(double Q2, double &F1, double &F2);

struct FormFactors
{
	double F1=0;
	double F2=0;
	double GA=0;
	double GP=0;

	formfactors F_vector;	
	formfactors F_axial;	

	void eval(double Q2)
	{
		F_vector(Q2, F1, F2);
		F_axial(Q2, GA, GP);
	}

};

//Proton and neutron formfactors from Kelly
void Fp_kelly(double Q2, double &F1, double &F2);
void Fn_kelly(double Q2, double &F1, double &F2);

//Isovector form factor Kelly
void FIV_kelly(double Q2, double &F1, double &F2);

//Dipole axial FF with Pion pole for GP
void FA_dipole(double Q2, double &GA, double &GP);

//Formfactors = 0
void Fzero(double Q2, double &G1, double &G2);


//For testing purposes:
void F1_only(double Q2, double &F1, double &F2);
void F2_only(double Q2, double &F1, double &F2);

void F2_one(double Q2, double &F1, double &F2);
void FA_only(double Q2, double &GA, double &GP);
void GP_only(double Q2, double &GA, double &GP);


#endif
