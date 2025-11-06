#include "Hadronic.h"
#include "Constants.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cstdarg>

int Get_Hadron(double TN, double pm, double cos, S_level *S, FormFactors *FF, complex<double> Hadron[4][4])
{

	//Set Hadron to zero just to be sure
	for (int mu = 0 ; mu < 4 ; mu++)
	{
		for (int nu = 0 ; nu < 4 ; nu++)
		{
			Hadron[mu][nu] = 0.;
		}
	}

	//Some kinematics
	double w = TN + S->E_kappa + pm*pm/2./S->MB; //non relativistic approximation for the recoil energy of system M_B, used thoughout

	double kN2 = TN*(TN+2*MN);
	double kN = sqrt(kN2);
	double q2 = kN2 - 2.*kN*pm*cos + pm*pm;

	double Q2 = q2 - w*w;
	if (Q2 < 0){ 
		//Unphysical region for scattering
		return -1; 
	}


	int ERR = S->eval(TN,pm,cos);
	if (ERR < 0){
		//Out of bounds, Hadron is zero
		return ERR;
	}

	//q = k_N - pm , system where all vectors are in the same plane and k_N ~ 3 
	double q = sqrt(q2);
	double sinq = -pm*sqrt(1.-cos*cos)/q; // Allways negative

	double cosq = (kN - pm*cos)/q; 

	double Q_z[4] = {w, q*sinq, 0., q*cosq}; // Vector q in k_N along z system obtained after rotation over -thetaN


	//Nucleon angles in q ~ z system
	double sin_qn = -sinq;
	double cos_qn = cosq;

	//Should be fine by contruction, but just to be sure:
	if (abs(cos_qn) > 1){
		return -1;
	}

	FF->eval(Q2);

	int nMJ = S->two_J+1;

	// g_{nm} Q^m
	double Q_z_l[4] = {w,-q*sinq, 0., -q*cosq};

	double hbc3 = pow(hbc,3); //to convert overlap to [MeV]^-3

	//Sum over all s and m, can be reduced
	for (int im = 0 ; im < nMJ; im++)
	{
		for (int is = 0 ; is < 2 ; is++)
		{

			complex<double> J_z[4]; //Hadron current in k_N~z system
			complex<double> J_l[4]; //Hadron current in q ~ z system
			for (int mu = 0 ; mu < 4 ; mu++)
			{
				J_z[mu] += FF->F1 * S->V[mu][is][im];
				J_z[mu] += -1.*FF->GA * S->A[mu][is][im];
				J_z[mu] += -1.*(FF->GP/2./MN) * Q_z[mu] * S->ps[is][im];

				for (int nu = 0 ; nu < 4 ; nu++)
				{
					J_z[mu] += (-1.*FF->F2/4./MN)*Q_z_l[nu]*S->T[mu][nu][is][im]; //Of course this is antisymmetric so can cut the loop in half
				} 

			}

			//Rotation to the q ~ z system over angle \theta_N
			J_l[0] = J_z[0];

			J_l[1] = cos_qn*J_z[1] + sin_qn*J_z[3];
			J_l[2] = J_z[2];
			J_l[3] = -sin_qn*J_z[1] + cos_qn*J_z[3];

			
			//Add to the Hadron tensor
			for (int mu = 0 ; mu < 4 ; mu++)
			{
				complex<double> daggerJmu = conj(J_l[mu])/hbc3; //Overlap has unit fm^{3/2} This makes H ~ Mev^{-3}, all units throughout are MeV

				for (int nu = 0 ; nu < 4 ; nu++)
				{
					Hadron[mu][nu] += daggerJmu*J_l[nu];  //In principle only 10 complex numbers are needed to define the tensor, but we keep the whole thing here.

				}

			}


		} //spin sum
	} // m sum
	

	return 0;	
}

int cosdep_Responses(double w, double q, S_level *S, FormFactors *FF)
{

	//We ignore the nuclear recoil:
	double TN = w - S->E_kappa;

	//Set the mass large to ignore recoil throughout:
	double MB = S->MB;
	S->MB = MB*1e6;

	int ntheta=50.;
	double Responses[5];
	for (int icos = 0 ; icos < ntheta ; icos++)
	{
		double costheta = -1 + icos*2./(ntheta-1.);
		complex<double> Hadron[4][4];
		Get_Hadron_q(TN, q , costheta, S, FF, Hadron);	

		cout << w <<  "   " << q  << "   " << costheta << "   ";

		//Inclusive responses:
		Responses[0] = abs(Hadron[0][0]);
		Responses[1] = real(Hadron[0][3]);
		Responses[2] = abs(Hadron[3][3]);
		Responses[3] = (abs(Hadron[1][1]) + abs(Hadron[2][2]) );
		Responses[4] = imag(Hadron[1][2]);

		for (int iR = 0 ; iR < 5 ; iR++)
		{
			cout << Responses[iR] << "   ";
		}

		//'B' factors, i.e. cos(\phi) dependence
		cout << real(Hadron[0][1]) << "  " <<  real(Hadron[1][3]) << "   " << imag(Hadron[0][2]) << "   " << imag(Hadron[2][3]) << "   ";

		//C factors: cos(2\pi)
		cout << abs(Hadron[1][1]) - abs(Hadron[2][2]) << "   ";

		//D factors : sin(\phi)
		cout << real(Hadron[0][2]) << "   " << real(Hadron[2][3]) << "   " << imag(Hadron[0][1]) << "   " << imag(Hadron[1][3]) << "   ";

		//E factor, sin(2\phi)
		cout << real(Hadron[1][2]) << endl;
//		
	}

	S->MB=MB; //Put the mass back

	return 0;
}


int Inclusive_Responses(double w, double q, S_level *S, FormFactors *FF, double Responses[5])
{
	//Do the integral over pm instead of \cos\theta since it's better phase space coverage
	//Additionally theres no two solutions for TN in this case
	double pm_max = S->pm_grid[S->N_pm-1]; //No point in going to larger pm than those in the grid
	double pm_min = S->pm_grid[0]; //Or smaller ones

	double R_CC = 0.;
	double R_CL = 0.;
	double R_LL = 0.;
	double R_T  = 0.;
	double R_Tp = 0.;
	complex<double> Hadron[4][4];
	double step = 1.;

	for (double pm = pm_min ; pm < pm_max ; pm +=step)
	{
		double TN = w - S->E_kappa - pm*pm/2./S->MB; 
		if (TN < 0){continue;}

		double kN2 = TN*(TN+2*MN);
		double kN = sqrt(kN2);
		double cos_pk = (q*q - kN2 - pm*pm)/(-2.*kN*pm);
		
		//Should be ok but check anyway
		if (abs(cos_pk) > 1){continue;}
	
		//This cosine is between kN and pm
		Get_Hadron(TN,pm,cos_pk, S, FF, Hadron);

		double Jac = pm/kN/q; //dcos_kq/dpm
	
		double PF = (kN*(TN+MN))/pow(2.*Pi, 2); // response is 1/MeV since Hadron ~ MeV^{-3} 

		//For the Recoil factor we need the angle between q and k
		double cos_kq = (pm*pm - q*q - kN2)/(-2.*kN*q);
		if (abs(cos_kq) > 1){continue;}
		double frec = abs( 1. +  (TN+MN)/S->MB *(1. - q/kN*cos_kq));

		Jac = Jac*PF/frec;

		R_CC += abs(Hadron[0][0])*Jac;
		R_CL += real(Hadron[0][3])*Jac;
		R_LL += abs(Hadron[3][3])*Jac;
		R_T += abs(Hadron[1][1] + Hadron[2][2])*Jac;
		R_Tp += imag(Hadron[2][3])*Jac;
	}

	Responses[0] = R_CC*step;
	Responses[1] = R_CL*step;
	Responses[2] = R_LL*step;
	Responses[3] = R_T*step;
	Responses[4] = R_Tp*step;

	return 0;

}


int Inclusive_Responses_th(double w, double q, S_level *S, FormFactors *FF, double Responses[5])
{
	//Do the integral over \theta_N to compute the inclusive response, neglecting nuclear recoil
	double TN = w - S->E_kappa; // Igore the recoil


	double R_CC = 0.;
	double R_CL = 0.;
	double R_LL = 0.;
	double R_T  = 0.;
	double R_Tp = 0.;

	if (TN < 0){
		Responses[0] = 0.;
		Responses[1] = 0.;
		Responses[2] = 0.;
		Responses[3] = 0.;
		Responses[4] = 0.;
		return 0;
	}

	double kN2 = TN*(TN+2.*MN);
	double kN = sqrt(kN2);
	double EN = TN + MN;

	double PF = (kN*EN)/pow(2.*Pi, 2); //  



	complex<double> Hadron[4][4];
	double step = 0.01;


	//To ignore recoil also further on we set MB large:
	double MBB = S->MB;
	S->MB = MBB*1.e6;


	int ntheta = 200; // probably overkill
	double th_step = Pi/(ntheta-1.);
	for (double ith = 0. ; ith < ntheta ; ith++)
	{
		double theta = ith*th_step;
		double costheta = cos(theta);
		Get_Hadron_q(TN, q, costheta, S, FF, Hadron);

		double Jac = th_step*sin(theta); 

		//We add following recoil factor using the actual mass of system B
		double frec = abs( 1. +  EN/MBB *(1. - q/kN*costheta));
		Jac = Jac/frec;

		R_CC += abs(Hadron[0][0])*Jac;
		R_CL += real(Hadron[0][3])*Jac;
		R_LL += abs(Hadron[3][3])*Jac;
		R_T += abs(Hadron[1][1] + Hadron[2][2])*Jac;
		R_Tp += imag(Hadron[2][3])*Jac;
	}

	Responses[0] = PF*R_CC;
	Responses[1] = PF*R_CL;
	Responses[2] = PF*R_LL;
	Responses[3] = PF*R_T;
	Responses[4] = PF*R_Tp;

	//Set the mass back to the original
	S->MB = MBB;	

	return 0;
}


int Get_Hadron_q(double TN, double q, double cos, S_level *S, FormFactors *FF, complex<double> Hadron[4][4])
{

	//cos is angle between k_N and q where q~3
	double kN2 = TN*(TN+2*MN);
	double kN = sqrt(kN2);

	double pm2 = kN2 - 2*q*kN*cos + q*q;
	double pm = sqrt(pm2);
	// p\cdotk = (kN sin)*kNsin + (cos k - q)*cosk = kn2 sin2 + cos2 k2 - qcos * k = (kn2 - kq cos)
	double cos_pk = (kN - q*cos)/pm;
	return Get_Hadron(TN, pm, cos_pk, S, FF, Hadron);
	
}

