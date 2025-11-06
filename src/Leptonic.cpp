#include "Leptonic.h"



int LeptonTensor::eval_kinematics(double Ei, double Ef, double costheta_f)
{
	if (Ef > Ei || Ei < mi || Ef < mf){
		return -1;
	}

	this->ki = sqrt(Ei*Ei - mi*mi);
	this->kf = sqrt(Ef*Ef - mf*mf);

	if (abs(costheta_f) > 1 ){
		return -1;
	}

	this->q = sqrt(ki*ki + kf*kf - 2.*ki*kf*costheta_f);
	this->w = Ei - Ef;
	this->Q2 = q*q - w*w;

	return 0;
}

int LeptonTensor::eval_minim(double Ei, double Ef, double costheta_f)
{

	int ERR = this->eval_kinematics(Ei,Ef,costheta_f);
	if (ERR < 0){
		this->set_zero();
		return ERR;
	}

	if (abs(costheta_f) > 1 ){
		this->set_zero();
		return -1;
	}

	double cosq = (ki - kf*costheta_f)/q;
	double sinq = kf/q * sqrt(1-costheta_f*costheta_f);

	double K_i[4] = {Ei, ki*sinq, 0. , ki*cosq};
	double K_f[4] = {Ef, K_i[1] , 0. , K_i[3] - q};

	double KidotKf = mi*mi - Ei*(Ei-Ef) + ki*q*cosq;
	double Mfac = mi*mf;	// g_{munu} * (mi*mf) usually neglected in relativistic limit
	double gmunuterm = Mfac - KidotKf;

	//diagonal parts
	this->L_s[0][0] = 2.*K_i[0]*K_f[0] + gmunuterm;
	this->L_s[1][1] = 2.*K_i[1]*K_f[1] - gmunuterm;
	this->L_s[2][2] =                  - gmunuterm;
	this->L_s[3][3] = 2*K_i[3] *K_f[3] - gmunuterm;

	//One time-like index, sign change because we are definining L_{mu\nu}
	this->L_s[0][1] = -K_i[0]*K_f[1] - K_i[1]*K_f[0];
	this->L_s[0][3] = -K_i[0]*K_f[3] - K_i[3]*K_f[0];

	this->L_s[1][3] = K_i[1]*K_f[3] + K_i[3]*K_f[1]; 

	//Antisymmetric terms:
	this->L_a[0][2] = - K_i[1]*K_f[3] + K_i[3]*K_f[1]; //Sign change
	this->L_a[1][2] =   K_i[0]*K_f[3] - K_i[3]*K_f[0];
	this->L_a[2][3] =   K_i[0]*K_f[1] - K_f[0]*K_f[1];

	this->FF = FX2(Q2);

	return 0;
}


double FX2_EM(double Q2)
{
	return pow(4*Pi*Alpha/Q2,2)/2.; //Factor 1/2 accounts for averaging over helicity initial lepton
}

double FX2_CC(double Q2)
{
	return pow(G_Fermi*Cabibbo,2); 

}

int LeptonTensor::eval(double Ei, double Ef, double costheta_f)
{

	this->set_zero();
	int ERR = this->eval_minim(Ei,Ef,costheta_f);
	
	if (ERR != 0){return ERR;}
	
	this->L_s[1][0] = this->L_s[0][1];
	this->L_s[3][0] = this->L_s[0][3];
	this->L_s[3][1] = this->L_s[1][3];

	this->L_a[2][0] = -1.*this->L_a[0][2];
	this->L_a[2][1] = -1.*this->L_a[1][2];
	this->L_a[3][2] = -1.*this->L_a[2][3];

	return 0;	

}

void LeptonTensor::set_zero()
{
	for (int mu = 0 ; mu < 4 ; mu++)
	{
		for (int nu = mu ; nu < 4 ; nu++)
		{
			this->L_s[mu][nu]=0.;
			this->L_a[mu][nu]=0.;
		}
	}
}

void LeptonTensor::set_elec()
{
	this->mi=electronmass;
	this->mf=electronmass;
	this->helicity = 0.;
	this->FX2 = FX2_EM;
}

void LeptonTensor::set_CC_numu()
{
	this->mi=0.;
	this->mf=muonmass;
	this->helicity = -1.;
	this->FX2 = FX2_CC;
}

