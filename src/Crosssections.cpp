#include "Crosssections.h"

#include "Constants.h"
#include "Hadronic.h"
#include "Kinematics.h"


//Evaluates the differential cross section ds/dE_d\Omega_l d\Omega_N
//uses as input missing momentum pm instead of nucleon scattering angle
double cross_section(double E_i, double E_f, double costheta_f, double pm, double phi_N, LeptonTensor *Lep, FormFactors *FF, S_level *S)
{

	//Evaluate lepton tensor and leptonic prefactors along with w, q , Q2, and the momenta ki, kf
	int ERR = Lep->eval_minim(E_i, E_f, costheta_f);

	if (ERR != 0){ return 0.;} //Invalid lepton kinematics, cross section is zero

	//Momentum transfer:
	double q = Lep->q; 

	//The outgoing nucleon kinetic energy
	double TN = E_i - E_f  - S->E_kappa - pm*pm/(2.*S->MB);

	if (TN < 0){return 0.;} 

	double kN2 = TN*(TN+2.*MN);
	double kN = sqrt(kN2);
	//Angle between pm and k_N
	double cos_pk = (q*q - kN2 - pm*pm)/(-2.*kN*q);

	//Get Hadron tensor elements, in q ~ z system with phi = 0
	complex<double> H[4][4]; //Overkill, could just define symmetric and antisymmetric
	ERR = Get_Hadron(TN, pm, cos_pk, S, FF, H);

	if (ERR != 0){return 0.;} //Hadron kinematics outside table range or otherwise invalid

	//The cross section then factorizes into factors A,B,C,D,E:
	
	//Symmetric contributions first
	double A = Lep->L_s[0][0] * abs(H[0][0]) + 2.*Lep->L_s[0][3] * real(H[0][3]) + Lep->L_s[3][3] * abs(H[3][3]); //L
	A += (Lep->L_s[1][1] + Lep->L_s[2][2])/2. * abs( H[1][1] + H[2][2]); //T

	double B = 2.*( Lep->L_s[0][1] * real(H[0][1]) + Lep->L_s[1][3] * real(H[1][3]));

	double C = (Lep->L_s[1][1] - Lep->L_s[2][2])/2. * abs( H[1][1] - H[2][2]);

	double D = -2.*( Lep->L_s[0][1] * real(H[0][2]) + Lep->L_s[1][3] * real(H[2][3]));

	double E = (Lep->L_s[2][2] - Lep->L_s[1][1])*real(H[1][2]);

	//Antisymmetric terms
	double h = Lep->helicity;
	if (h != 0){
		A += -2*h*Lep->L_a[1][2] * imag(H[1][2]);

		B += -2.*h*( Lep->L_a[0][2] * imag(H[0][2])  + Lep->L_a[2][3] * imag(H[2][3]) );

		D += -2.*h*( Lep->L_a[0][2] * imag(H[0][1])  + Lep->L_a[2][3] * imag(H[1][3]) );
	}

	double LepHad = A + B*cos(phi_N) + C*cos(2.*phi_N) + D*sin(phi_N) + E*sin(2.*phi_N);
	
	//Phase space factors and coupling
	double cos_qn = (kN - pm*cos_pk)/q; //nucleon angle with respect to q
	//Recoil factor
	double frec = abs( 1. +  (TN+MN)/S->MB *(1. - q/kN*cos_qn));
	
	double PS = (kN*(TN+MN))/frec/pow(2.*Pi, 5); //Not including MB/E_B ~ 1
	PS = PS * Lep->kf/E_i;

	return Lep->FF*PS*LepHad; //dsigma/ dEl d\Omega_l d\Omega_N;
}


//Evaluates the differential cross sections ds/dE_d\Omega_l d\Omega_N
//uses as input \cos\theta_N, the scattering angle of the nucleon with respect to momentum transfer q.
//For any scattering angle there may be two physical final-states (with different nucleon momentum)
//The function returns the number of solutions and stores the values of the CS and nucleon momenta for each in the arrays CS and k_N_sols
int cross_section_cos(double E_i, double E_f, double costheta_f, double costhetaN, double phi_N, LeptonTensor *Lep, FormFactors *FF, S_level *S, double CS[2], double k_N_sols[2])
{
	//When specifying the nucleon scattering angle, there can be two solutions for the nucleon energy
	//The function returns the number of solutions and the cross section and nucleon momentum for each one
	//Most often the 'small solution' is completely negligible or unphysical
	
	int ERR = Lep->eval_kinematics(E_i, E_f, costheta_f);

	if (ERR != 0){return 0;} //Bad lepton kinematics nr solutions is 0
		
	int nr_solutions = Kproton(S->MA, S->MB, Lep->w, Lep->q, costhetaN, k_N_sols);

	CS[0]=0.;
	CS[1]=0.;
	for (int i = 0 ; i < nr_solutions; i++)
	{
		//Calculate missing momentum for the solution
		double kN= k_N_sols[i];
		double pm = sqrt(kN*kN + Lep->q*Lep->q - 2.*kN*Lep->q*costhetaN);
		
		//Calculate CS
		CS[i] = cross_section(E_i, E_f, costheta_f, pm, phi_N, Lep, FF, S);
	}

	return nr_solutions;

}

//returns the inclusive cross section ds/dE_ld\Omega_l
// The partial cross sections \sigma_i, with i = {CC,CL,LL,T,Tp} are stored in the array CS_i
double CS_inclusive(double E_i, double E_f, double costheta_f, LeptonTensor *Lep, FormFactors *FF, S_level *S, double CS_i[5])
{
	//Evaluate leptonic part
	int ERR = Lep->eval_minim(E_i,E_f,costheta_f);
	if (ERR != 0){
		//invalid lepton kinematics
		CS_i[0] = 0.;
		CS_i[1] = 0.;
		CS_i[2] = 0.;
		CS_i[3] = 0.;
		CS_i[4] = 0.;
		return 0.;
	}

	//Get inclusive Responses
	double Responses[5] = {0.,0.,0.,0.,0.};
	Inclusive_Responses(Lep->w, Lep->q, S, FF, Responses);

	double PF = 1./pow(2.*Pi,2) * Lep->kf/E_i * Lep->FF;
	CS_i[0] = Responses[0] * Lep->L_s[0][0] * PF;
	CS_i[1] = 2.*Responses[1] * Lep->L_s[0][3] * PF;
	CS_i[2] = Responses[2] * Lep->L_s[3][3] * PF;

	CS_i[3] = Responses[3] * (Lep->L_s[1][1] + Lep->L_s[2][2])/2. * PF;
	CS_i[4] = -2.*Lep->helicity* Responses[4] * Lep->L_a[1][2] * PF;

	return CS_i[0] + CS_i[1] + CS_i[2] + CS_i[3] + CS_i[4];
}
